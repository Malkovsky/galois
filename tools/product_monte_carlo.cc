#include <fcntl.h>
#include <openssl/evp.h>
#include <openssl/rand.h>
#include <sys/file.h>
#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <optional>
#include <sstream>
#include <thread>
#include <utility>

#include "product_monte_carlo_data.h"
#include "product_monte_carlo_trials.h"

namespace mc {
namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;

std::string Timestamp() {
  const auto now = std::time(nullptr);
  std::tm tm{};
  gmtime_r(&now, &tm);
  std::ostringstream out;
  out << std::put_time(&tm, "%Y-%m-%dT%H:%M:%S+00:00");
  return out.str();
}

class File {
 public:
  File(const fs::path& path, int flags) {
    fd_ = ::open(path.c_str(), flags | O_CLOEXEC, 0600);
    if (fd_ < 0) {
      throw std::system_error(errno, std::generic_category(), path.string());
    }
#ifdef GF256_MC_TEST_HOOKS
    path_ = path;
#endif
  }
  ~File() { ::close(fd_); }
  File(const File&) = delete;
  File& operator=(const File&) = delete;
  void Write(std::string_view bytes) {
#ifdef GF256_MC_TEST_HOOKS
    if (path_.filename() == "journal.jsonl" &&
        std::getenv("MC_TEST_SYNC_ORDER")) {
      const auto record = Parse(std::string(bytes));
      if (record.contains("flip end")) {
        Require(U64(record.at("flip end")) <= synced_flip_end_,
                "journal referenced unsynced flip data");
      }
    }
#endif
    while (!bytes.empty()) {
      auto n = ::write(fd_, bytes.data(), bytes.size());
      if (n < 0 && errno == EINTR) {
        continue;
      }
      if (n <= 0) {
        throw std::system_error(n < 0 ? errno : EIO, std::generic_category(),
                                "write");
      }
      bytes.remove_prefix(static_cast<size_t>(n));
    }
  }
  void Sync() {
    int status;
    do {
      status = ::fsync(fd_);
    } while (status < 0 && errno == EINTR);
    if (status < 0) {
      throw std::system_error(errno, std::generic_category(), "fsync");
    }
#ifdef GF256_MC_TEST_HOOKS
    if (path_.filename() == "flips.bin") {
      synced_flip_end_ = Offset();
    }
#endif
  }
  uint64_t Offset() const {
    auto pos = ::lseek(fd_, 0, SEEK_CUR);
    Require(pos >= 0, "file offset unavailable");
    return static_cast<uint64_t>(pos);
  }
  void Lock() {
    if (::flock(fd_, LOCK_EX | LOCK_NB) < 0) {
      throw std::system_error(errno, std::generic_category(), "run lock");
    }
  }

 private:
  int fd_;
#ifdef GF256_MC_TEST_HOOKS
  fs::path path_;
  inline static uint64_t synced_flip_end_ = 0;
#endif
};

void SyncDirectory(const fs::path& path) {
  File dir(path.empty() ? fs::path(".") : path, O_RDONLY | O_DIRECTORY);
  dir.Sync();
}
void AtomicJson(const fs::path& path, const Json& value) {
  auto temp = path;
  temp += ".tmp";
  {
    File out(temp, O_WRONLY | O_CREAT | O_TRUNC | O_NOFOLLOW);
    out.Write(Dump(value, true) + "\n");
    out.Sync();
  }
  fs::rename(temp, path);
  SyncDirectory(path.parent_path());
}
std::string Hash(std::string_view data) {
  std::string result(32, '\0');
  unsigned size = 0;
  Require(EVP_Digest(data.data(), data.size(),
                     reinterpret_cast<unsigned char*>(result.data()), &size,
                     EVP_sha256(), nullptr) == 1 &&
              size == 32,
          "OpenSSL SHA256 failed");
  return result;
}
std::string Hex(std::string_view data) {
  constexpr char digits[] = "0123456789abcdef";
  std::string out;
  for (unsigned char c : data) {
    out += digits[c >> 4];
    out += digits[c & 15];
  }
  return out;
}
std::string ReadExact(std::istream& stream, size_t n) {
  std::string text(n, '\0');
  stream.read(text.data(), n);
  Require(static_cast<size_t>(stream.gcount()) == n, "truncated input");
  return text;
}
// Bounded line reads distinguish an incomplete final fragment from a complete
// malformed record. A corrupt interior/full line is never silently skipped.
std::optional<std::string> ReadLine(std::istream& stream, bool& incomplete) {
  std::string text;
  char c;
  while (stream.get(c)) {
    if (c == '\n') {
      return text;
    }
    Require(text.size() < 131072, "journal line exceeds schema size limit");
    text += c;
  }
  Require(stream.eof() && !stream.bad(), "journal read failed");
  incomplete = !text.empty();
  return std::nullopt;
}
void PutLE(std::string& bytes, size_t offset, uint64_t value, size_t width) {
  for (size_t i = 0; i < width; ++i) {
    bytes[offset + i] = static_cast<char>(value >> (8 * i));
  }
}
uint64_t GetLE(std::string_view bytes, size_t offset, size_t width) {
  uint64_t value = 0;
  for (size_t i = 0; i < width; ++i) {
    value |= uint64_t(static_cast<unsigned char>(bytes[offset + i])) << (8 * i);
  }
  return value;
}

struct Result {
  uint64_t index = 0;
  std::array<uint64_t, 22> metrics{};
  std::vector<uint32_t> positions;
  std::exception_ptr error;
  bool ready = false;
};

class Workers {
 public:
  Workers(const Settings& settings, bool recorded)
      : s_(settings),
        recorded_(recorded),
        slots_(std::min(s_.checkpoint, 8 * s_.threads)) {
    try {
      for (uint64_t i = 0; i < std::min(s_.threads, uint64_t(slots_.size()));
           ++i) {
        threads_.emplace_back([this] { Work(); });
      }
    } catch (...) {
      Shutdown();
      throw;
    }
  }
  ~Workers() { Shutdown(); }
  void Begin(uint64_t batch, uint64_t k) {
    std::lock_guard lock(mutex_);
    batch_ = batch;
    k_ = k;
    next_ = cursor_ = completed_ = 0;
    active_ = true;
    halted_ = false;
    changed_.notify_all();
  }
  // The coordinator alone consumes slots, in trial order. Credits include
  // running and finished-but-not-consumed blocks, bounding memory even when
  // trial zero is slow. Workers refill independently, without a wave barrier.
  bool Pop(Result& result, uint64_t& completed, bool& done) {
    std::unique_lock lock(mutex_);
    auto finished = [&] {
      return cursor_ == next_ &&
             (next_ == s_.size || halted_ || product_interrupted());
    };
    changed_.wait_for(lock, std::chrono::milliseconds(50), [&] {
      return finished() ||
             (cursor_ < next_ && slots_[cursor_ % slots_.size()].ready);
    });
    completed = completed_;
    done = finished();
    if (done) {
      active_ = false;
      return false;
    }
    if (cursor_ == next_ || !slots_[cursor_ % slots_.size()].ready) {
      return false;
    }
    auto& slot = slots_[cursor_ % slots_.size()];
    result = std::move(slot);
    slot = Result{};
    ++cursor_;
    changed_.notify_all();
    return true;
  }

 private:
  void Shutdown() noexcept {
    {
      std::lock_guard lock(mutex_);
      shutdown_ = true;
      changed_.notify_all();
    }
    for (auto& t : threads_) {
      if (t.joinable()) {
        t.join();
      }
    }
  }
  void Work() {
    for (;;) {
      uint64_t index;
      {
        std::unique_lock lock(mutex_);
        changed_.wait(lock, [&] {
          return shutdown_ ||
                 (active_ && !halted_ && !product_interrupted() &&
                  next_ < s_.size && next_ - cursor_ < slots_.size());
        });
        if (shutdown_) {
          return;
        }
        // Claiming a block is its start boundary. Once claimed, it runs the
        // complete product decoder even if a signal arrives immediately after.
        if (product_interrupted()) {
          continue;
        }
        index = next_++;
      }
      auto& slot = slots_[index % slots_.size()];
      slot.index = index;
      try {
        if (recorded_) {
          slot.positions.resize(
              std::max<uint64_t>(1, std::min(k_, 524288 - k_)));
        }
#ifdef GF256_MC_TEST_HOOKS
        // Dedicated, noninstalled test executable only.
        if (const char* fail = std::getenv("MC_TEST_FAIL")) {
          if (index == std::stoull(fail)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            throw std::runtime_error("injected worker failure");
          }
        }
        if (const char* slow = std::getenv("MC_TEST_SLOW_FIRST")) {
          if (index == 0) {
            std::this_thread::sleep_for(
                std::chrono::milliseconds(std::stoul(slow)));
          }
        }
#endif
        int status = product_trial_reference(
            s_.seed, batch_, index, k_, s_.passes, s_.anchors, s_.binary,
            slot.metrics.data(), recorded_ ? 1 : 0,
            recorded_ ? slot.positions.data() : nullptr, 0, nullptr);
        Require(status == 0, "native status=" + std::to_string(status));
      } catch (...) {
        slot.error = std::current_exception();
      }
      {
        std::lock_guard lock(mutex_);
        if (slot.error) {
          halted_ = true;
        } else {
          ++completed_;
        }
        slot.ready = true;
        changed_.notify_all();
      }
    }
  }
  Settings s_;
  bool recorded_;
  std::vector<Result> slots_;
  std::vector<std::thread> threads_;
  std::mutex mutex_;
  std::condition_variable changed_;
  uint64_t batch_ = 0, k_ = 0, next_ = 0, cursor_ = 0, completed_ = 0;
  bool active_ = false, halted_ = false, shutdown_ = false;
};

void WriteFlips(File& file, uint64_t batch, uint64_t k, const Result& result) {
  const auto count = std::min(k, 524288 - k);
  std::string bytes(208 + 4 * count, '\0');
  PutLE(bytes, 0, batch, 8);
  PutLE(bytes, 8, result.index, 8);
  PutLE(bytes, 16, k, 4);
  PutLE(bytes, 20, count, 4);
  bytes[24] = k > 262144;
  for (size_t i = 0; i < 22; ++i) {
    PutLE(bytes, 32 + 8 * i, result.metrics[i], 8);
  }
  for (size_t i = 0; i < count; ++i) {
    PutLE(bytes, 208 + 4 * i, result.positions[i], 4);
  }
  file.Write(bytes);
  file.Write(Hash(bytes));
}

void Run(const fs::path& directory, const Settings& s, bool recorded) {
  Require(product_interrupt_install() == 0,
          "could not install signal handlers");
  Json metadata;
  metadata["schema revision"] = 2;
  metadata["created at"] = Timestamp();
  metadata["settings"] = s.ToJson();
  metadata["codeword"] = "zero";
  metadata["code"] = kCode;
  metadata["random algorithm"] = recorded ? kFisherYates : kFloyd;
  AtomicJson(directory / "metadata.json", metadata);
  const auto digest = Hash(Dump(metadata));
  Aggregate aggregate(Hex(digest), 2);
  AtomicJson(directory / "summary.json", aggregate.Summary());
  File journal(directory / "journal.jsonl", O_WRONLY | O_CREAT | O_EXCL);
  File log(directory / "progress.log", O_WRONLY | O_CREAT | O_EXCL);
  std::optional<File> flips;
  if (recorded) {
    flips.emplace(directory / "flips.bin", O_WRONLY | O_CREAT | O_EXCL);
    flips->Write("RSFLIP01" + digest);
    flips->Sync();
  }
  const bool tty = ::isatty(STDERR_FILENO);
  bool bar_visible = false;
  auto progress = [&](const std::string& text, bool console = true) {
    auto line = Timestamp() + " " + text + "\n";
    if (console) {
      if (bar_visible) {
        std::cerr << "\r\033[K";
        bar_visible = false;
      }
      std::cerr << line << std::flush;
    }
    log.Write(line);
  };
  progress("root seed=" + std::to_string(s.seed) +
           " settings persisted before trials");
  journal.Sync();
  log.Sync();
  SyncDirectory(directory);
  SyncDirectory(directory.parent_path());

  const auto start = Clock::now();
  auto last_sync = start, last_report = start, last_bar = start,
       last_checkpoint = start;
  uint64_t completed_total = 0, increment = 0;
  uint64_t batch_completed = 0, batch = 0, k = 0;
  Stats staged;
  uint64_t first = 0, end = 0, flip_start = 0;
  std::string error;
  auto throughput = [&](Clock::time_point now) {
    const double elapsed = std::chrono::duration<double>(now - start).count();
    const double rate =
        elapsed > 0 ? static_cast<double>(completed_total) / elapsed : 0;
    std::ostringstream out;
    out << std::fixed << std::setprecision(3) << "wall blocks/s=" << rate
        << " information MiB/s=" << rate * 56896 / 1048576
        << " elapsed seconds=" << elapsed;
    return out.str();
  };
  auto counts = [&] {
    return "batch=" + std::to_string(batch) + " k=" + std::to_string(k) +
           " trials=" + std::to_string(batch_completed) + "/" +
           std::to_string(s.size) +
           " overall trials=" + std::to_string(completed_total);
  };
  auto checkpoint = [&] {
    if (staged.blocks != 0) {
      // Preflight every aggregate and derived total before publishing a record.
#ifdef GF256_MC_TEST_HOOKS
      if (increment == 1 && std::getenv("MC_TEST_COUNTER_OVERFLOW")) {
        auto exhausted = aggregate.overall;
        exhausted.iterations = UINT64_MAX;
        exhausted.Add(staged);
      }
#endif
      auto next = aggregate.PrepareAdd(k, staged);
      const auto next_increment = CheckedAdd(increment, 1);
      Json r;
      r["schema revision"] = 2;
      r["run identity"] = aggregate.identity;
      r["increment id"] = Number(increment);
      r["batch id"] = batch;
      r["first trial index"] = first;
      r["past last trial index"] = end;
      r["trial count"] = end - first;
      r["flipped bit count"] = k;
      r["statistics"] = staged.ToJson();
      if (flips) {
        flips->Sync();
        r["flip start"] = flip_start;
        r["flip end"] = flips->Offset();
      }
      journal.Write(Dump(r) + "\n");
      aggregate.Commit(std::move(next));
      increment = next_increment;
      staged = Stats{};
    }
    last_checkpoint = Clock::now();
  };
  auto durable = [&] {
    checkpoint();
    journal.Sync();
    log.Sync();
    AtomicJson(directory / "summary.json", aggregate.Summary());
    last_sync = Clock::now();
  };
  auto bar = [&](Clock::time_point now, bool force = false) {
    if (tty && (force || now - last_bar >= std::chrono::milliseconds(200))) {
      double fraction = double(batch_completed) / s.size;
      size_t width = static_cast<size_t>(30 * fraction);
      std::cerr << "\r[" << std::string(width, '#')
                << std::string(30 - width, ' ') << "] batch=" << batch
                << " k=" << k << " trials=" << batch_completed << "/" << s.size
                << " " << std::fixed << std::setprecision(1) << 100 * fraction
                << "% overall trials=" << completed_total << " "
                << throughput(now) << "\033[K" << std::flush;
      bar_visible = true;
      last_bar = now;
    }
  };
  try {
    Workers workers(s, recorded);
    while (!product_interrupted() && error.empty() &&
           (!s.batches || batch < s.batches)) {
      k = product_batch_k(s.seed, batch, s.lo, s.hi);
      batch_completed = 0;
      progress(counts());
      workers.Begin(batch, k);
      bool done = false;
      while (!done) {
        Result result;
        uint64_t observed;
        bool ready = workers.Pop(result, observed, done);
        completed_total =
            CheckedAdd(completed_total, observed - batch_completed);
        batch_completed = observed;
        if (ready) {
          if (result.error) {
            checkpoint();
            try {
              std::rethrow_exception(result.error);
            } catch (const std::exception& e) {
              if (error.empty()) {
                error = "batch=" + std::to_string(batch) +
                        " index=" + std::to_string(result.index) +
                        " failed: " + e.what();
              }
            } catch (...) {
              if (error.empty()) {
                error = "unknown worker exception";
              }
            }
          } else {
            if (staged.blocks != 0 && end != result.index) {
              checkpoint();
            }
            if (staged.blocks == 0) {
              first = result.index;
              if (flips) {
                flip_start = flips->Offset();
              }
            }
            auto next_staged = staged;
            next_staged.Add(result.metrics);
            const auto next_end = CheckedAdd(result.index, 1);
            if (flips) {
              WriteFlips(*flips, batch, k, result);
            }
            end = next_end;
            staged = next_staged;
            if (staged.blocks == s.checkpoint) {
              checkpoint();
            }
          }
        }
        auto now = Clock::now();
        if (now - last_checkpoint >= std::chrono::seconds(1)) {
          checkpoint();
        }
        if (now - last_sync >= std::chrono::seconds(s.sync)) {
          durable();
        }
        if (now - last_report >= std::chrono::seconds(s.report)) {
          progress(counts() + " persisted overall trials=" +
                       std::to_string(aggregate.overall.blocks) + " " +
                       throughput(now),
                   !tty);
          last_report = now;
        }
        bar(now);
      }
      checkpoint();
      bar(Clock::now(), true);
      progress("batch=" + std::to_string(batch) + " k=" + std::to_string(k) +
               " completed trials=" + std::to_string(batch_completed) + "/" +
               std::to_string(s.size) +
               " overall trials=" + std::to_string(completed_total) + " " +
               throughput(Clock::now()));
      Require(batch != UINT64_MAX, "batch identity exhausted; start a new run");
      ++batch;
    }
  } catch (const std::exception& e) {
    // Worker failures are drained above. Coordinator I/O/allocation failures
    // stop and join workers. Never retry a possibly partially written record
    // or replace the last good summary with a partially updated aggregate.
    // Report mode can recover the complete journal prefix after storage repair.
    try {
      journal.Sync();
      log.Sync();
    } catch (...) {
    }
    throw;
  }
  durable();
  progress("finalizing trials=" + std::to_string(completed_total) +
           " interrupted=" + (product_interrupted() ? "True" : "False") +
           " error=" + (error.empty() ? "None" : error) + " " +
           throughput(Clock::now()));
  log.Sync();
  Require(error.empty(), error);
}

void ReadFlips(std::istream& stream,
               const Json& r,
               const Settings& s,
               bool replay,
               bool random,
               unsigned schema) {
  auto start = stream.tellg();
  Require(start >= 0 && static_cast<uint64_t>(start) == U64(r.at("flip start")),
          "noncontiguous flip offsets");
  uint64_t batch = U64(r.at("batch id")),
           first = U64(r.at("first trial index")),
           end = U64(r.at("past last trial index")),
           k = U64(r.at("flipped bit count"));
  Stats stats;
  Json old = LegacyStats();
  const size_t count = std::min(k, 524288 - k);
  std::vector<uint32_t> positions(std::max<size_t>(1, count));
  std::vector<bool> selected(524288);
  for (uint64_t trial = first; trial < end; ++trial) {
    auto body = ReadExact(stream, 208);
    Require(GetLE(body, 0, 8) == batch && GetLE(body, 8, 8) == trial &&
                GetLE(body, 16, 4) == k,
            "flip trial identity mismatch");
    Require(GetLE(body, 20, 4) == count &&
                GetLE(body, 24, 1) == uint64_t(k > 262144),
            "invalid flip count/complement");
    body += ReadExact(stream, count * 4);
    Require(ReadExact(stream, 32) == Hash(body), "corrupt flip checksum");
    std::fill(selected.begin(), selected.end(), false);
    for (size_t i = 0; i < count; ++i) {
      auto p = GetLE(body, 208 + 4 * i, 4);
      Require(p < selected.size() && !selected[p],
              "invalid or duplicate flip position");
      selected[p] = true;
      positions[i] = static_cast<uint32_t>(p);
    }
    std::array<uint64_t, 22> metrics{};
    for (size_t i = 0; i < 22; ++i) {
      metrics[i] = GetLE(body, 32 + 8 * i, 8);
    }
    if (replay) {
      std::array<uint64_t, 22> actual{};
      int status = product_trial_reference(
          s.seed, batch, trial, k, s.passes, s.anchors, s.binary, actual.data(),
          2, positions.data(), random, nullptr);
      Require(status == 0 && actual == metrics,
              "replay mismatch batch=" + std::to_string(batch) + " trial=" +
                  std::to_string(trial) + " status=" + std::to_string(status));
    }
    stats.Add(metrics);
    if (schema == 1) {
      AddLegacyTrial(old, metrics);
    }
  }
  auto end_offset = stream.tellg();
  Require(end_offset >= 0 &&
              static_cast<uint64_t>(end_offset) == U64(r.at("flip end")),
          "flip end offset mismatch");
  Require((schema == 1 ? old : stats.ToJson()) == r.at("statistics"),
          "flip metrics disagree with journal");
}

void Recover(const fs::path& directory, bool replay) {
  std::ifstream meta_file(directory / "metadata.json", std::ios::binary);
  Require(meta_file.good(), "cannot read metadata");
  std::string text;
  char c;
  while (meta_file.get(c)) {
    Require(text.size() < 131072, "metadata exceeds size limit");
    text += c;
  }
  Require(meta_file.eof() && !meta_file.bad(), "metadata read failed");
  Json metadata = Parse(text);
  std::set<std::string> fields{"schema revision", "created at", "settings",
                               "code", "random algorithm"};
  if (metadata.contains("codeword")) {
    fields.insert("codeword");
  }
  Fields(metadata, fields);
  auto schema = U64(metadata.at("schema revision"));
  Require((schema == 1 || schema == 2) && metadata.at("code") == Json(kCode) &&
              (metadata.at("random algorithm") == Json(kFloyd) ||
               metadata.at("random algorithm") == Json(kFisherYates)),
          "incompatible metadata schema/code/random algorithm");
  Require(metadata.at("created at").is_string(), "invalid created at");
  const auto codeword = metadata.contains("codeword")
                            ? metadata.at("codeword").get<std::string>()
                            : "random";
  Require(codeword == "zero" || codeword == "random",
          "incompatible codeword convention");
  auto settings = Settings::FromJson(metadata.at("settings"));
  // Never insert defaults before hashing legacy metadata.
  const auto digest = Hash(Dump(metadata));
  Aggregate aggregate(Hex(digest), static_cast<unsigned>(schema));
  const bool recorded = metadata.at("random algorithm") == Json(kFisherYates);
  Require(!replay || recorded, "this run has no saved flips");
  std::ifstream journal(directory / "journal.jsonl", std::ios::binary), flips;
  Require(journal.good(), "cannot read journal");
  if (recorded) {
    flips.open(directory / "flips.bin", std::ios::binary);
    Require(flips.good() && ReadExact(flips, 40) == "RSFLIP01" + digest,
            "invalid flip file identity/version");
  }
  uint64_t index = 0;
  std::optional<std::pair<uint64_t, uint64_t>> previous;
  bool incomplete = false;
  while (auto line = ReadLine(journal, incomplete)) {
    const auto next_index = CheckedAdd(index, 1);
    try {
      Json r = Parse(*line);
      fields = {"schema revision", "run identity",      "increment id",
                "batch id",        "first trial index", "past last trial index",
                "trial count",     "flipped bit count", "statistics"};
      if (recorded) {
        fields.insert("flip start");
        fields.insert("flip end");
      }
      Fields(r, fields);
      Require(U64(r.at("schema revision")) == schema &&
                  r.at("run identity") == Json(aggregate.identity),
              "incompatible journal identity/schema");
      Require(Natural(r.at("increment id")) == index,
              "duplicate or out-of-order increment id");
      uint64_t batch = U64(r.at("batch id")),
               first = U64(r.at("first trial index")),
               end = U64(r.at("past last trial index")),
               count = U64(r.at("trial count"));
      Require((settings.batches == 0 || batch < settings.batches) &&
                  first < end && end <= settings.size && end - first == count &&
                  count <= settings.checkpoint,
              "invalid trial range/count");
      if (previous) {
        Require(batch >= previous->first &&
                    (batch != previous->first || first >= previous->second),
                "overlapping or out-of-order trial ranges");
      }
      uint64_t k =
          product_batch_k(settings.seed, batch, settings.lo, settings.hi);
      Require(U64(r.at("flipped bit count")) == k,
              "batch k disagrees with seed/settings");
      Stats stats;
      if (schema == 1) {
        const auto& old = r.at("statistics");
        ValidateLegacy(old, count, k, settings.passes);
        stats.blocks = count;
        stats.iterations = Natural(old.at(kMetrics[12]).at("sum"));
        stats.info_raw = Natural(old.at(kMetrics[2]).at("sum"));
        stats.info_post = Natural(old.at(kMetrics[6]).at("sum"));
        stats.full_raw = Natural(old.at(kMetrics[0]).at("sum"));
        stats.full_post = Natural(old.at(kMetrics[4]).at("sum"));
      } else {
        stats = Stats::FromJson(r.at("statistics"), count, k, settings.passes);
      }
      if (recorded) {
        ReadFlips(flips, r, settings, replay, codeword == "random", schema);
      }
      aggregate.Add(k, stats, r.at("statistics"));
      previous = {batch, end};
      index = next_index;
    } catch (const std::exception& e) {
      throw std::runtime_error("journal line " + std::to_string(next_index) +
                               ": " + e.what());
    }
  }
  if (recorded && flips.peek() != std::char_traits<char>::eof()) {
    std::cerr << "unreferenced flip tail ignored (not committed trials)\n";
  }
  AtomicJson(directory / "summary.json", aggregate.Summary());
  std::cerr << Timestamp() << " regenerated " << aggregate.overall.blocks
            << " trials; ignored incomplete tail="
            << (incomplete ? "True" : "False") << "\n";
  if (replay) {
    std::cerr << "verified replay: " << aggregate.overall.blocks
              << " trials, all 22 metrics match\n";
  }
}

int Main(int argc, char** argv) {
  Settings s;
  std::string mode, sampler = "floyd";
  fs::path directory;
  bool seeded = false;
  int options = 0;
  for (int i = 1; i < argc; ++i) {
    std::string flag = argv[i], value;
    auto equals = flag.find('=');
    if (equals != std::string::npos) {
      value = flag.substr(equals + 1);
      flag.resize(equals);
    }
    if (flag == "--help" || flag == "-h") {
      std::cout
          << "Native fixed-weight RS product Monte Carlo (schema 2)\n"
             "Usage: rs-product-monte-carlo --output NEW_DIRECTORY [options]\n"
             "       rs-product-monte-carlo --report DIRECTORY | --replay "
             "DIRECTORY\n"
             "--seed UINT64 (default: generated, printed and persisted before "
             "work)\n"
             "--batch-size 1000 --batches 0 (infinite) --threads 1 (1..1024)\n"
             "--minimum-flipped-bits 2500 --maximum-flipped-bits 2700 "
             "(inclusive, 0..524288)\n"
             "--max-directional-passes 16 (2..1000000)\n"
             "--[no-]anchors --[no-]binary-image (both enabled)\n"
             "--sampler floyd|fisher-yates (default floyd; Fisher-Yates saves "
             "RSFLIP01)\n"
             "--checkpoint-trials 64 (1..4096, also bounds concurrency)\n"
             "--report-seconds 2 --fsync-seconds 5 (1..86400)\n"
             "One uniform k per batch; exactly k flips per all-zero block. No "
             "resume.\n"
             "Counters are exact uint64 JSON integers "
             "(0..18446744073709551615); "
             "overflow is an error.\n"
             "Total bits must also fit uint64; iterations are directional "
             "passes.\n"
             "Legacy schema 1 sums and squared sums exceeding uint64 are "
             "rejected.\n"
             "SIGINT/SIGTERM stop starts, drain whole in-flight blocks, "
             "checkpoint and fsync.\n"
             "Crash loss: bounded worker window/staged increment plus journal "
             "writes since fsync.\n"
             "Ordered journal flush at least once/second subject to in-flight "
             "completion and I/O.\n"
             "Saved flips are fsynced before journal references; replay checks "
             "all 22 private metrics.\n"
             "Legacy absent codeword means random; report ignores only an "
             "incomplete final line.\n";
      return 0;
    }
    ++options;
    if (flag == "--anchors" || flag == "--no-anchors" ||
        flag == "--binary-image" || flag == "--no-binary-image") {
      Require(equals == std::string::npos, "boolean flags take no value");
      if (flag == "--anchors" || flag == "--no-anchors") {
        s.anchors = flag == "--anchors";
      } else {
        s.binary = flag == "--binary-image";
      }
      continue;
    }
    if (equals == std::string::npos) {
      Require(i + 1 < argc, "missing value for " + flag);
      value = argv[++i];
    }
    if (flag == "--output" || flag == "--report" || flag == "--replay") {
      Require(mode.empty() && !value.empty(),
              "specify exactly one output/report/replay directory");
      mode = flag;
      directory = value;
      continue;
    }
    if (flag == "--sampler") {
      sampler = value;
      continue;
    }
    Require(!value.empty() &&
                std::all_of(value.begin(), value.end(),
                            [](char c) { return c >= '0' && c <= '9'; }),
            "expected unsigned decimal integer for " + flag);
    auto v = Decimal(value);
    if (flag == "--seed") {
      s.seed = v;
      seeded = true;
    } else if (flag == "--batch-size") {
      s.size = v;
    } else if (flag == "--batches") {
      s.batches = v;
    } else if (flag == "--threads") {
      s.threads = v;
    } else if (flag == "--minimum-flipped-bits") {
      s.lo = v;
    } else if (flag == "--maximum-flipped-bits") {
      s.hi = v;
    } else if (flag == "--max-directional-passes") {
      s.passes = v;
    } else if (flag == "--checkpoint-trials") {
      s.checkpoint = v;
    } else if (flag == "--report-seconds") {
      s.report = v;
    } else if (flag == "--fsync-seconds") {
      s.sync = v;
    } else {
      throw std::runtime_error("unknown option: " + flag);
    }
  }
  Require(!mode.empty(),
          "specify exactly one of --output, --report or --replay");
  if (mode != "--output") {
    Require(options == 1, mode + " accepts only a directory");
    File lock(directory / "run.lock", O_WRONLY | O_CREAT | O_NOFOLLOW);
    lock.Lock();
    Recover(directory, mode == "--replay");
    return 0;
  }
  s.Validate();
  Require(sampler == "floyd" || sampler == "fisher-yates", "invalid sampler");
  if (!seeded) {
    Require(RAND_bytes(reinterpret_cast<unsigned char*>(&s.seed),
                       sizeof(s.seed)) == 1,
            "OpenSSL random seed generation failed");
  }
  Require(fs::create_directory(directory),
          "output directory must not already exist");
  File lock(directory / "run.lock", O_WRONLY | O_CREAT | O_EXCL);
  lock.Lock();
  Run(directory, s, sampler == "fisher-yates");
  return 0;
}
}  // namespace mc

int main(int argc, char** argv) {
  try {
    return mc::Main(argc, argv);
  } catch (const std::exception& e) {
    std::cerr << mc::Timestamp() << " error: " << e.what() << "\n";
    return 1;
  }
}
