#include <fcntl.h>
#include <getopt.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <memory>
#include <mutex>
#include <queue>
#include <span>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

#include "indicators/indeterminate_progress_bar.hpp"
#include "indicators/progress_bar.hpp"
#include "lin_chung_han/kernels.h"
#include "reed_solomon/lch_decoder.h"
#include "reed_solomon/lch_encoder.h"
#include "xxhash.h"

namespace {

namespace fs = std::filesystem;
using gf2p8::Element;
using gf2p8::lch::Backend;
using gf2p8::lch::Radix;
using gf2p8::lch::Status;
using gf2p8::rs::LCHDecoder;
using gf2p8::rs::LCHEncoder;

constexpr char kMagic[8] = {'L', 'C', 'H', 'R', 'S', '2', '5', '6'};
constexpr uint16_t kVersion = 3;
constexpr uint16_t kHeaderSize = 64;
constexpr uint16_t kMaxFilenameLen = 4096;
constexpr uint32_t kDefaultChunkSize = 16384;
constexpr int kExitOk = 0;
constexpr int kExitRepairable = 1;
constexpr int kExitUnrecoverable = 2;
constexpr int kExitUsage = 3;

enum class ProgressMode { automatic, always, never };

void WriteU16(uint8_t* out, uint16_t value) {
  out[0] = static_cast<uint8_t>(value);
  out[1] = static_cast<uint8_t>(value >> 8);
}

void WriteU32(uint8_t* out, uint32_t value) {
  out[0] = static_cast<uint8_t>(value);
  out[1] = static_cast<uint8_t>(value >> 8);
  out[2] = static_cast<uint8_t>(value >> 16);
  out[3] = static_cast<uint8_t>(value >> 24);
}

void WriteU64(uint8_t* out, uint64_t value) {
  for (int i = 0; i < 8; ++i) {
    out[i] = static_cast<uint8_t>(value >> (8 * i));
  }
}

uint16_t ReadU16(const uint8_t* in) {
  return static_cast<uint16_t>(in[0] | (static_cast<uint16_t>(in[1]) << 8));
}

uint32_t ReadU32(const uint8_t* in) {
  return static_cast<uint32_t>(in[0]) | (static_cast<uint32_t>(in[1]) << 8) |
         (static_cast<uint32_t>(in[2]) << 16) |
         (static_cast<uint32_t>(in[3]) << 24);
}

uint64_t ReadU64(const uint8_t* in) {
  uint64_t value = 0;
  for (int i = 0; i < 8; ++i) {
    value |= static_cast<uint64_t>(in[i]) << (8 * i);
  }
  return value;
}

bool MulOverflow(uint64_t a, uint64_t b, uint64_t* out) {
  if (a != 0 && b > UINT64_MAX / a) {
    return true;
  }
  *out = a * b;
  return false;
}

bool AddOverflow(uint64_t a, uint64_t b, uint64_t* out) {
  if (b > UINT64_MAX - a) {
    return true;
  }
  *out = a + b;
  return false;
}

bool FitsSize(uint64_t value, size_t* out) {
  if (value > static_cast<uint64_t>(static_cast<size_t>(-1))) {
    return false;
  }
  *out = static_cast<size_t>(value);
  return true;
}

const char* BackendName(Backend backend) {
  switch (backend) {
    case Backend::tuned:
      return "tuned";
    case Backend::scalar:
      return "scalar";
    case Backend::ssse3:
      return "ssse3";
    case Backend::avx2:
      return "avx2";
    case Backend::gfni128_affine:
      return "gfni128_affine";
    case Backend::gfni256_affine:
      return "gfni256_affine";
    case Backend::gfni512_affine:
      return "gfni512_affine";
  }
  return "unknown";
}

void PrintHelp(FILE* out) {
  std::fputs(
      "Usage:\n"
      "  lch-rs encode [-] [<file>] -k <K> -r <R> [--chunk-size C] [-o <dir>] "
      "[--name <name>] [--jobs N] [--progress MODE] [-v] [--force]\n"
      "  lch-rs verify [--jobs N] [--progress MODE] [-v] <share> "
      "[<share> ...]\n"
      "  lch-rs decode [-o <dest>|-] [--force] [--jobs N] "
      "[--progress MODE] [-v] <share> [<share> ...]\n"
      "\n"
      "Progress modes: auto (default), always, never.\n"
      "\n"
      "Examples:\n"
      "  lch-rs encode data.bin -k 10 -r 4\n"
      "  lch-rs verify data.bin.lch.14.*\n"
      "  lch-rs decode -o restored.bin data.bin.lch.14.*\n"
      "  tar -cf - dir | lch-rs encode -k 10 -r 4 --name backup\n"
      "  lch-rs decode backup.lch.14.* | tar -xf -\n",
      out);
}

void Error(const char* message) {
  std::fprintf(stderr, "lch-rs: %s\n", message);
}

template <typename... Args>
void Verbose(bool enabled, Args&&... args) {
  if (!enabled) {
    return;
  }
  std::cerr << "lch-rs: ";
  (std::cerr << ... << std::forward<Args>(args));
  std::cerr << '\n';
}

bool ProgressEnabled(ProgressMode mode) {
  return mode == ProgressMode::always ||
         (mode == ProgressMode::automatic && isatty(STDERR_FILENO));
}

std::string FormatBytes(double bytes) {
  constexpr const char* units[] = {"B", "KiB", "MiB", "GiB", "TiB", "PiB"};
  size_t unit = 0;
  while (bytes >= 1024.0 && unit + 1 < std::size(units)) {
    bytes /= 1024.0;
    ++unit;
  }
  char text[64];
  if (unit == 0) {
    std::snprintf(text, sizeof(text), "%.0f %s", bytes, units[unit]);
  } else {
    std::snprintf(text, sizeof(text), "%.1f %s", bytes, units[unit]);
  }
  return text;
}

class ProgressDisplay {
 public:
  ProgressDisplay(std::string label,
                  ProgressMode mode,
                  bool total_known,
                  uint64_t total)
      : total_(total),
        label_(std::move(label)),
        start_(std::chrono::steady_clock::now()) {
    if (!ProgressEnabled(mode)) {
      return;
    }
    if (total_known) {
      bar_ = std::make_unique<indicators::ProgressBar>(
          indicators::option::BarWidth{30},
          indicators::option::PrefixText{label_ + " "},
          indicators::option::Start{"["}, indicators::option::Fill{"="},
          indicators::option::Lead{">"}, indicators::option::Remainder{" "},
          indicators::option::End{"]"},
          indicators::option::ShowPercentage{true},
          indicators::option::ShowElapsedTime{true},
          indicators::option::ShowRemainingTime{true},
          indicators::option::Stream{std::cerr});
    } else {
      spinner_ =
          std::make_unique<indicators::IndeterminateProgressBar>(
              indicators::option::BarWidth{30},
              indicators::option::PrefixText{label_ + " "},
              indicators::option::Start{"["}, indicators::option::Fill{"."},
              indicators::option::Lead{"<=>"}, indicators::option::End{"]"},
              indicators::option::Stream{std::cerr});
    }
    reporter_ = std::thread([this] { Report(); });
  }

  ProgressDisplay(const ProgressDisplay&) = delete;
  ProgressDisplay& operator=(const ProgressDisplay&) = delete;

  ~ProgressDisplay() { Finish(false); }

  void Add(uint64_t bytes) {
    completed_.fetch_add(bytes, std::memory_order_relaxed);
  }

  void Finish(bool success) {
    if (!reporter_.joinable()) {
      return;
    }
    {
      std::lock_guard lock(wait_mutex_);
      stop_ = true;
    }
    wake_.notify_one();
    reporter_.join();
    Render(true, success);
  }

 private:
  std::string Postfix(bool failed) const {
    const uint64_t completed = completed_.load(std::memory_order_relaxed);
    const double elapsed = std::chrono::duration<double>(
                               std::chrono::steady_clock::now() - start_)
                               .count();
    const double rate = elapsed > 0.0 ? static_cast<double>(completed) / elapsed
                                      : 0.0;
    std::string text = FormatBytes(static_cast<double>(completed));
    if (bar_) {
      text += "/" + FormatBytes(static_cast<double>(total_));
    }
    text += " " + FormatBytes(rate) + "/s";
    if (failed) {
      text += " failed";
    }
    return text;
  }

  void Render(bool final, bool success) {
    const bool failed = final && !success;
    if (bar_) {
      bar_->set_option(indicators::option::PostfixText{Postfix(failed)});
      if (final && success) {
        bar_->set_progress(100);
        return;
      }
      const uint64_t completed = completed_.load(std::memory_order_relaxed);
      const size_t percent =
          total_ == 0
              ? 0
              : static_cast<size_t>(std::min<long double>(
                    99.0L, static_cast<long double>(completed) * 100.0L /
                               static_cast<long double>(total_)));
      bar_->set_progress(percent);
      if (final) {
        bar_->mark_as_completed();
      }
      return;
    }
    spinner_->set_option(indicators::option::PostfixText{Postfix(failed)});
    if (final) {
      spinner_->mark_as_completed();
    } else {
      spinner_->tick();
    }
  }

  void Report() {
    Render(false, false);
    std::unique_lock lock(wait_mutex_);
    while (!wake_.wait_for(lock, std::chrono::milliseconds(100),
                           [this] { return stop_; })) {
      lock.unlock();
      Render(false, false);
      lock.lock();
    }
  }

  const uint64_t total_;
  const std::string label_;
  const std::chrono::steady_clock::time_point start_;
  std::unique_ptr<indicators::ProgressBar> bar_;
  std::unique_ptr<indicators::IndeterminateProgressBar> spinner_;
  std::atomic<uint64_t> completed_{0};
  std::thread reporter_;
  std::mutex wait_mutex_;
  std::condition_variable wake_;
  bool stop_ = false;
};

class UniqueFd {
 public:
  UniqueFd() = default;
  explicit UniqueFd(int fd, bool owned = true) : fd_(fd), owned_(owned) {}
  ~UniqueFd() { Reset(); }
  UniqueFd(const UniqueFd&) = delete;
  UniqueFd& operator=(const UniqueFd&) = delete;
  UniqueFd(UniqueFd&& other) noexcept : fd_(other.fd_), owned_(other.owned_) {
    other.fd_ = -1;
    other.owned_ = true;
  }
  UniqueFd& operator=(UniqueFd&& other) noexcept {
    if (this != &other) {
      Reset();
      fd_ = other.fd_;
      owned_ = other.owned_;
      other.fd_ = -1;
      other.owned_ = true;
    }
    return *this;
  }
  int get() const { return fd_; }
  bool valid() const { return fd_ >= 0; }
  void Reset() {
    if (owned_ && fd_ >= 0) {
      close(fd_);
    }
    fd_ = -1;
    owned_ = true;
  }

 private:
  int fd_ = -1;
  bool owned_ = true;
};

class AlignedBuffer {
 public:
  bool Allocate(size_t bytes) {
    data_.reset();
    if (bytes == 0) {
      return true;
    }
    void* pointer = nullptr;
    if (posix_memalign(&pointer, 64, bytes) != 0) {
      return false;
    }
    data_.reset(static_cast<Element*>(pointer));
    std::memset(pointer, 0, bytes);
    return true;
  }
  Element* data() { return data_.get(); }

 private:
  struct Free {
    void operator()(Element* pointer) const { std::free(pointer); }
  };
  std::unique_ptr<Element, Free> data_;
};

class StreamingHash {
 public:
  bool Initialize() {
    state_.reset(XXH3_createState());
    return state_ != nullptr && XXH3_64bits_reset(state_.get()) == XXH_OK;
  }

  bool Update(const void* data, size_t bytes) {
    return XXH3_64bits_update(state_.get(), data, bytes) == XXH_OK;
  }

  uint64_t Digest() const {
    return static_cast<uint64_t>(XXH3_64bits_digest(state_.get()));
  }

 private:
  struct Free {
    void operator()(XXH3_state_t* state) const { XXH3_freeState(state); }
  };
  std::unique_ptr<XXH3_state_t, Free> state_;
};

bool ParseU64(const char* text, uint64_t* value) {
  if (text == nullptr || text[0] == '\0') {
    return false;
  }
  errno = 0;
  char* end = nullptr;
  const unsigned long long parsed = std::strtoull(text, &end, 10);
  if (errno != 0 || end == text || *end != '\0') {
    return false;
  }
  *value = static_cast<uint64_t>(parsed);
  return true;
}

uint32_t DefaultJobs() {
  const unsigned hardware = std::thread::hardware_concurrency();
  return hardware == 0 ? 1 : hardware;
}

bool FilenameValid(std::string_view name) {
  if (name.empty() || name.size() > kMaxFilenameLen) {
    return false;
  }
  return name.find('\0') == std::string_view::npos &&
         name.find('/') == std::string_view::npos &&
         name.find('\\') == std::string_view::npos;
}

std::string ShareFileName(const std::string& name, uint16_t n, uint16_t id) {
  return name + ".lch." + std::to_string(n) + "." + std::to_string(id);
}

bool ParseShareFileName(const std::string& name,
                        std::string* basename,
                        uint16_t* n,
                        uint16_t* id) {
  const auto last = name.rfind('.');
  if (last == std::string::npos || last + 1 >= name.size()) {
    return false;
  }
  uint64_t parsed_id = 0;
  if (!ParseU64(name.c_str() + last + 1, &parsed_id) || parsed_id > 0xffff) {
    return false;
  }
  const std::string rest = name.substr(0, last);
  const auto last2 = rest.rfind('.');
  if (last2 == std::string::npos || last2 + 1 >= rest.size()) {
    return false;
  }
  uint64_t parsed_n = 0;
  if (!ParseU64(rest.c_str() + last2 + 1, &parsed_n) || parsed_n > 0xffff) {
    return false;
  }
  if (last2 < 4 || rest.compare(last2 - 4, 4, ".lch") != 0) {
    return false;
  }
  *basename = rest.substr(0, last2 - 4);
  if (!FilenameValid(*basename)) {
    return false;
  }
  *n = static_cast<uint16_t>(parsed_n);
  *id = static_cast<uint16_t>(parsed_id);
  return true;
}

bool PWriteAll(int fd, const void* buffer, size_t bytes, uint64_t offset) {
  const auto* data = static_cast<const uint8_t*>(buffer);
  size_t written = 0;
  while (written < bytes) {
    const ssize_t result = pwrite(fd, data + written, bytes - written,
                                  static_cast<off_t>(offset + written));
    if (result < 0) {
      if (errno == EINTR) {
        continue;
      }
      return false;
    }
    if (result == 0) {
      return false;
    }
    written += static_cast<size_t>(result);
  }
  return true;
}

bool WriteAll(int fd, const void* buffer, size_t bytes) {
  const auto* data = static_cast<const uint8_t*>(buffer);
  size_t written = 0;
  while (written < bytes) {
    const ssize_t result = write(fd, data + written, bytes - written);
    if (result < 0) {
      if (errno == EINTR) {
        continue;
      }
      return false;
    }
    if (result == 0) {
      return false;
    }
    written += static_cast<size_t>(result);
  }
  return true;
}

bool PReadRange(int fd,
                void* buffer,
                size_t bytes,
                uint64_t offset,
                size_t* got) {
  auto* data = static_cast<uint8_t*>(buffer);
  size_t read_bytes = 0;
  while (read_bytes < bytes) {
    const ssize_t result = pread(fd, data + read_bytes, bytes - read_bytes,
                                 static_cast<off_t>(offset + read_bytes));
    if (result < 0) {
      if (errno == EINTR) {
        continue;
      }
      return false;
    }
    if (result == 0) {
      break;
    }
    read_bytes += static_cast<size_t>(result);
  }
  if (read_bytes < bytes) {
    std::memset(data + read_bytes, 0, bytes - read_bytes);
  }
  if (got != nullptr) {
    *got = read_bytes;
  }
  return true;
}

bool ReadUpTo(int fd, void* buffer, size_t bytes, size_t* got) {
  auto* data = static_cast<uint8_t*>(buffer);
  size_t read_bytes = 0;
  while (read_bytes < bytes) {
    const ssize_t result = read(fd, data + read_bytes, bytes - read_bytes);
    if (result < 0) {
      if (errno == EINTR) {
        continue;
      }
      return false;
    }
    if (result == 0) {
      break;
    }
    read_bytes += static_cast<size_t>(result);
  }
  if (got != nullptr) {
    *got = read_bytes;
  }
  return true;
}

uint64_t HashBytes(const void* data, size_t bytes) {
  return static_cast<uint64_t>(XXH3_64bits(data, bytes));
}

struct Header {
  uint16_t version = 0;
  uint16_t header_size = 0;
  uint16_t filename_len = 0;
  uint16_t share_id = 0;
  uint16_t k = 0;
  uint16_t r = 0;
  uint32_t c = 0;
  uint32_t m = 0;
  uint64_t original_size = 0;
  uint64_t whole_file_hash = 0;
  uint64_t metadata_hash = 0;
};

void PackHeader(uint8_t out[kHeaderSize], const Header& header) {
  std::memset(out, 0, kHeaderSize);
  std::memcpy(out, kMagic, 8);
  WriteU16(out + 8, header.version);
  WriteU16(out + 10, header.header_size);
  WriteU16(out + 12, header.filename_len);
  WriteU16(out + 14, header.share_id);
  WriteU16(out + 16, header.k);
  WriteU16(out + 18, header.r);
  WriteU32(out + 20, header.c);
  WriteU32(out + 24, header.m);
  WriteU64(out + 28, header.original_size);
  WriteU64(out + 36, header.whole_file_hash);
  WriteU64(out + 44, header.metadata_hash);
}

bool UnpackHeader(const uint8_t in[kHeaderSize], Header* header) {
  if (std::memcmp(in, kMagic, 8) != 0) {
    return false;
  }
  header->version = ReadU16(in + 8);
  header->header_size = ReadU16(in + 10);
  header->filename_len = ReadU16(in + 12);
  header->share_id = ReadU16(in + 14);
  header->k = ReadU16(in + 16);
  header->r = ReadU16(in + 18);
  header->c = ReadU32(in + 20);
  header->m = ReadU32(in + 24);
  header->original_size = ReadU64(in + 28);
  header->whole_file_hash = ReadU64(in + 36);
  header->metadata_hash = ReadU64(in + 44);
  for (int i = 52; i < 64; ++i) {
    if (in[i] != 0) {
      return false;
    }
  }
  return header->version == kVersion && header->header_size == kHeaderSize;
}

void WriteHashTable(uint8_t* out, const uint64_t* values, size_t count) {
  for (size_t i = 0; i < count; ++i) {
    WriteU64(out + i * 8, values[i]);
  }
}

void ReadHashTable(const uint8_t* in, uint64_t* values, size_t count) {
  for (size_t i = 0; i < count; ++i) {
    values[i] = ReadU64(in + i * 8);
  }
}

uint64_t MetadataHash(const Header& header,
                      std::string_view filename,
                      const uint64_t* chunk_hashes,
                      size_t chunk_count) {
  Header hashed = header;
  hashed.metadata_hash = 0;
  uint8_t packed[kHeaderSize];
  PackHeader(packed, hashed);
  StreamingHash hash;
  if (!hash.Initialize()) {
    return 0;
  }
  hash.Update(packed, kHeaderSize);
  hash.Update(filename.data(), filename.size());
  std::vector<uint8_t> table(chunk_count * 8);
  if (chunk_count != 0) {
    WriteHashTable(table.data(), chunk_hashes, chunk_count);
    hash.Update(table.data(), table.size());
  }
  return hash.Digest();
}

struct Geometry {
  uint16_t k = 0;
  uint16_t r = 0;
  uint32_t c = 0;
  uint32_t m = 0;
  uint64_t original_size = 0;
  uint64_t whole_file_hash = 0;
  uint64_t stripe_bytes = 0;
  uint64_t payload_offset = 0;
  std::string filename;
};

uint16_t ShareCount(const Geometry& geometry) {
  return static_cast<uint16_t>(geometry.k + geometry.r);
}

struct AcceptedShare {
  UniqueFd fd;
  uint16_t share_id = 0;
  std::vector<uint64_t> chunk_hashes;
};

bool StripeLiveBytes(const Geometry& geometry,
                     uint32_t stripe,
                     uint64_t* live) {
  uint64_t offset = 0;
  if (MulOverflow(stripe, geometry.stripe_bytes, &offset)) {
    return false;
  }
  if (offset >= geometry.original_size) {
    *live = 0;
    return true;
  }
  *live = geometry.original_size - offset;
  if (*live > geometry.stripe_bytes) {
    *live = geometry.stripe_bytes;
  }
  return true;
}

bool ChunkSizeForLiveBytes(const Geometry& geometry,
                           uint64_t live,
                           uint32_t* chunk_size) {
  if (live == 0 || live > geometry.stripe_bytes) {
    return false;
  }
  const uint64_t size = (live - 1) / geometry.k + 1;
  if (size > geometry.c) {
    return false;
  }
  *chunk_size = static_cast<uint32_t>(size);
  return true;
}

bool StripeChunkSize(const Geometry& geometry,
                     uint32_t stripe,
                     uint32_t* chunk_size) {
  if (geometry.m == 0 || stripe >= geometry.m) {
    return false;
  }
  if (stripe + 1 < geometry.m) {
    *chunk_size = geometry.c;
    return true;
  }
  uint64_t live = 0;
  return StripeLiveBytes(geometry, stripe, &live) &&
         ChunkSizeForLiveBytes(geometry, live, chunk_size);
}

bool PayloadBytes(const Geometry& geometry, uint64_t* payload_bytes) {
  if (geometry.m == 0) {
    if (geometry.original_size != 0) {
      return false;
    }
    *payload_bytes = 0;
    return true;
  }
  uint32_t final_chunk_size = 0;
  uint64_t prefix_bytes = 0;
  return StripeChunkSize(geometry, geometry.m - 1, &final_chunk_size) &&
         !MulOverflow(geometry.m - 1, geometry.c, &prefix_bytes) &&
         !AddOverflow(prefix_bytes, final_chunk_size, payload_bytes);
}

bool DimensionsValid(uint16_t k, uint16_t r) {
  LCHEncoder encoder(k, r);
  LCHDecoder decoder(k, r);
  return encoder.Valid() && decoder.Valid();
}

bool GeometryValid(const Geometry& geometry) {
  if (geometry.k == 0 || geometry.r == 0 || geometry.c == 0 ||
      !DimensionsValid(geometry.k, geometry.r)) {
    return false;
  }
  if (geometry.original_size == 0) {
    return geometry.m == 0;
  }
  if (geometry.stripe_bytes == 0) {
    return false;
  }
  const uint64_t expected_m =
      (geometry.original_size - 1) / geometry.stripe_bytes + 1;
  uint64_t payload_bytes = 0;
  return expected_m == geometry.m && PayloadBytes(geometry, &payload_bytes);
}

bool ParseProgressMode(std::string_view text, ProgressMode* mode) {
  if (text == "auto") {
    *mode = ProgressMode::automatic;
  } else if (text == "always") {
    *mode = ProgressMode::always;
  } else if (text == "never") {
    *mode = ProgressMode::never;
  } else {
    return false;
  }
  return true;
}

struct Options {
  std::string command;
  std::string input;
  std::string output;
  std::string name;
  std::vector<std::string> shares;
  uint16_t k = 0;
  uint16_t r = 0;
  uint32_t chunk_size = kDefaultChunkSize;
  uint32_t jobs = DefaultJobs();
  ProgressMode progress_mode = ProgressMode::automatic;
  bool have_output = false;
  bool force = false;
  bool verbose = false;
};

int ParseCommandOptions(int argc, char** argv, Options* options) {
  options->command = argv[0];
  const bool is_encode = options->command == "encode";
  const bool is_verify = options->command == "verify";
  const bool is_decode = options->command == "decode";
  if (!is_encode && !is_verify && !is_decode) {
    Error("unknown command");
    PrintHelp(stderr);
    return kExitUsage;
  }

  const option long_opts[] = {
      {"chunk-size", required_argument, nullptr, 1},
      {"jobs", required_argument, nullptr, 2},
      {"force", no_argument, nullptr, 3},
      {"name", required_argument, nullptr, 4},
      {"progress", required_argument, nullptr, 5},
      {"verbose", no_argument, nullptr, 'v'},
      {"help", no_argument, nullptr, 'h'},
      {nullptr, 0, nullptr, 0},
  };
  const char* short_opts = is_encode ? "k:r:o:hv" : "o:hv";
  optind = 1;
  int flag = 0;
  while ((flag = getopt_long(argc, argv, short_opts, long_opts, nullptr)) !=
         -1) {
    uint64_t parsed = 0;
    switch (flag) {
      case 'k':
        if (!is_encode || !ParseU64(optarg, &parsed) || parsed == 0 ||
            parsed > 0xffff) {
          Error("invalid -k");
          return kExitUsage;
        }
        options->k = static_cast<uint16_t>(parsed);
        break;
      case 'r':
        if (!is_encode || !ParseU64(optarg, &parsed) || parsed == 0 ||
            parsed > 0xffff) {
          Error("invalid -r");
          return kExitUsage;
        }
        options->r = static_cast<uint16_t>(parsed);
        break;
      case 'o':
        if (is_verify) {
          Error("unknown option");
          return kExitUsage;
        }
        options->output = optarg;
        options->have_output = true;
        break;
      case 'h':
        PrintHelp(stdout);
        std::exit(kExitOk);
      case 'v':
        options->verbose = true;
        break;
      case 1:
        if (!is_encode || !ParseU64(optarg, &parsed) || parsed == 0 ||
            parsed > UINT32_MAX) {
          Error("invalid --chunk-size");
          return kExitUsage;
        }
        options->chunk_size = static_cast<uint32_t>(parsed);
        break;
      case 2:
        if (!ParseU64(optarg, &parsed) || parsed == 0 || parsed > UINT32_MAX) {
          Error("invalid --jobs");
          return kExitUsage;
        }
        options->jobs = static_cast<uint32_t>(parsed);
        break;
      case 3:
        if (is_verify) {
          Error("unknown option");
          return kExitUsage;
        }
        options->force = true;
        break;
      case 4:
        if (!is_encode || optarg == nullptr || !FilenameValid(optarg)) {
          Error("invalid --name");
          return kExitUsage;
        }
        options->name = optarg;
        break;
      case 5:
        if (!ParseProgressMode(optarg, &options->progress_mode)) {
          Error("invalid --progress");
          return kExitUsage;
        }
        break;
      default:
        PrintHelp(stderr);
        return kExitUsage;
    }
  }

  for (int i = optind; i < argc; ++i) {
    if (is_encode) {
      if (!options->input.empty()) {
        Error("encode accepts one input");
        return kExitUsage;
      }
      options->input = argv[i];
    } else {
      options->shares.emplace_back(argv[i]);
    }
  }
  if (is_encode) {
    if (options->k == 0 || options->r == 0) {
      Error("encode requires -k and -r");
      return kExitUsage;
    }
    if (options->input.empty() || options->input == "-") {
      if (options->name.empty()) {
        Error("encode from stdin requires --name");
        return kExitUsage;
      }
    }
  } else if (options->shares.empty()) {
    Error("share list required");
    return kExitUsage;
  }
  return kExitOk;
}

class StripeQueue {
 public:
  explicit StripeQueue(size_t capacity) : capacity_(capacity) {}

  bool Push(uint32_t stripe, AlignedBuffer buffer, size_t live) {
    std::unique_lock lock(mutex_);
    not_full_.wait(lock, [&] { return queue_.size() < capacity_ || stop_; });
    if (stop_) {
      return false;
    }
    queue_.push(Item{stripe, std::move(buffer), live});
    not_empty_.notify_one();
    return true;
  }

  bool Pop(uint32_t* stripe, AlignedBuffer* buffer, size_t* live) {
    std::unique_lock lock(mutex_);
    not_empty_.wait(lock, [&] { return !queue_.empty() || done_ || stop_; });
    if (stop_) {
      return false;
    }
    if (queue_.empty()) {
      return false;
    }
    Item item = std::move(queue_.front());
    queue_.pop();
    *stripe = item.stripe;
    *buffer = std::move(item.buffer);
    *live = item.live;
    not_full_.notify_one();
    return true;
  }

  void Done() {
    std::lock_guard lock(mutex_);
    done_ = true;
    not_empty_.notify_all();
  }

  void Stop() {
    std::lock_guard lock(mutex_);
    stop_ = true;
    not_empty_.notify_all();
    not_full_.notify_all();
  }

 private:
  struct Item {
    uint32_t stripe = 0;
    AlignedBuffer buffer;
    size_t live = 0;
  };
  const size_t capacity_;
  std::mutex mutex_;
  std::condition_variable not_empty_;
  std::condition_variable not_full_;
  std::queue<Item> queue_;
  bool done_ = false;
  bool stop_ = false;
};

void SetError(std::atomic<int>* error, int code) {
  int expected = 0;
  error->compare_exchange_strong(expected, code);
}

struct EncodeState {
  const LCHEncoder* encoder = nullptr;
  Geometry geometry;
  std::vector<UniqueFd> share_fds;
  std::vector<uint64_t> chunk_hashes;
  std::mutex hash_mutex;
  std::atomic<int> error{0};
};

bool EncodeStripe(const LCHEncoder& encoder,
                  Element* storage,
                  size_t chunk_size,
                  std::span<Element> workspace) {
  const size_t k = encoder.DataCount();
  const size_t r = encoder.RecoveryCount();
  std::vector<const Element*> data(k);
  std::vector<Element*> recovery(r);
  for (size_t i = 0; i < k; ++i) {
    data[i] = storage + i * chunk_size;
  }
  for (size_t i = 0; i < r; ++i) {
    recovery[i] = storage + (k + i) * chunk_size;
  }
  return encoder.Encode(data, recovery, chunk_size, workspace, Backend::tuned,
                        Radix::radix4) == Status::ok;
}

bool ProcessEncodeStripe(EncodeState* state,
                         uint32_t stripe,
                         size_t live,
                         AlignedBuffer* data,
                         AlignedBuffer* work,
                         std::span<Element> workspace) {
  const uint16_t n = ShareCount(state->geometry);
  const uint32_t c = state->geometry.c;
  uint32_t encoded_chunk_size = 0;
  if (!ChunkSizeForLiveBytes(state->geometry, live, &encoded_chunk_size)) {
    return false;
  }
  const size_t encoded_bytes = static_cast<size_t>(n) * encoded_chunk_size;
  std::memset(work->data(), 0, encoded_bytes);
  std::memcpy(work->data(), data->data(), live);
  if (!EncodeStripe(*state->encoder, work->data(), encoded_chunk_size,
                     workspace)) {
    return false;
  }
  uint64_t payload_index = 0;
  if (MulOverflow(stripe, c, &payload_index)) {
    return false;
  }
  uint64_t offset = 0;
  if (AddOverflow(state->geometry.payload_offset, payload_index, &offset)) {
    return false;
  }
  std::vector<uint64_t> digests(n);
  for (uint16_t i = 0; i < n; ++i) {
    const Element* chunk =
        work->data() + static_cast<size_t>(i) * encoded_chunk_size;
    if (!PWriteAll(state->share_fds[i].get(), chunk, encoded_chunk_size,
                   offset)) {
      return false;
    }
    digests[i] = HashBytes(chunk, encoded_chunk_size);
  }
  {
    std::lock_guard lock(state->hash_mutex);
    std::copy(digests.begin(), digests.end(),
              state->chunk_hashes.begin() + static_cast<size_t>(stripe) * n);
  }
  return true;
}

int FinishEncode(EncodeState* state) {
  const uint16_t n = ShareCount(state->geometry);
  Header header;
  header.version = kVersion;
  header.header_size = kHeaderSize;
  header.filename_len = static_cast<uint16_t>(state->geometry.filename.size());
  header.k = state->geometry.k;
  header.r = state->geometry.r;
  header.c = state->geometry.c;
  header.m = state->geometry.m;
  header.original_size = state->geometry.original_size;
  header.whole_file_hash = state->geometry.whole_file_hash;
  uint64_t payload_bytes = 0;
  if (!PayloadBytes(state->geometry, &payload_bytes)) {
    Error("payload overflow");
    return kExitUsage;
  }
  uint64_t hash_offset = 0;
  if (AddOverflow(state->geometry.payload_offset, payload_bytes,
                  &hash_offset)) {
    Error("share size overflow");
    return kExitUsage;
  }
  uint64_t hash_bytes = 0;
  if (MulOverflow(state->geometry.m, 8, &hash_bytes)) {
    Error("hash table overflow");
    return kExitUsage;
  }
  uint64_t file_size = 0;
  if (AddOverflow(hash_offset, hash_bytes, &file_size)) {
    Error("share size overflow");
    return kExitUsage;
  }
  std::vector<uint8_t> table(static_cast<size_t>(hash_bytes));
  std::vector<uint64_t> share_hashes(state->geometry.m);
  for (uint16_t i = 0; i < n; ++i) {
    header.share_id = i;
    for (uint32_t stripe = 0; stripe < state->geometry.m; ++stripe) {
      share_hashes[stripe] =
          state->chunk_hashes[static_cast<size_t>(stripe) * n + i];
    }
    header.metadata_hash =
        MetadataHash(header, state->geometry.filename, share_hashes.data(),
                     share_hashes.size());
    uint8_t packed[kHeaderSize];
    PackHeader(packed, header);
    const int fd = state->share_fds[i].get();
    if (state->geometry.m != 0) {
      WriteHashTable(table.data(), share_hashes.data(), share_hashes.size());
    }
    if (!PWriteAll(fd, packed, kHeaderSize, 0) ||
        (state->geometry.m != 0 &&
         !PWriteAll(fd, table.data(), table.size(), hash_offset)) ||
        ftruncate(fd, static_cast<off_t>(file_size)) != 0) {
      Error("failed to write share metadata");
      return kExitUsage;
    }
  }
  return kExitOk;
}

int RunEncode(const Options& options) {
  if (!DimensionsValid(options.k, options.r)) {
    Error("unsupported -k/-r");
    return kExitUsage;
  }

  UniqueFd source;
  std::string name = options.name;
  fs::path out_dir = ".";
  const bool read_stdin = options.input.empty() || options.input == "-";
  uint64_t source_size = 0;
  if (read_stdin) {
    source = UniqueFd(STDIN_FILENO, false);
  } else {
    const fs::path input_path(options.input);
    std::error_code exists_error;
    if (!fs::is_regular_file(input_path, exists_error)) {
      Error("input file not found");
      return kExitUsage;
    }
    source_size = fs::file_size(input_path, exists_error);
    if (exists_error) {
      Error("failed to determine input size");
      return kExitUsage;
    }
    if (options.name.empty()) {
      name = input_path.filename().string();
    }
    if (!options.have_output) {
      out_dir = input_path.parent_path();
      if (out_dir.empty()) {
        out_dir = ".";
      }
    }
    source = UniqueFd(open(input_path.c_str(), O_RDONLY));
    if (!source.valid()) {
      Error("failed to open input");
      return kExitUsage;
    }
  }
  if (options.have_output && options.output != "-") {
    out_dir = options.output;
  }
  if (!FilenameValid(name)) {
    Error("invalid share name");
    return kExitUsage;
  }
  std::error_code dir_error;
  if (!fs::is_directory(out_dir, dir_error)) {
    Error("output directory not found");
    return kExitUsage;
  }

  uint64_t stripe_bytes = 0;
  if (MulOverflow(options.k, options.chunk_size, &stripe_bytes)) {
    Error("chunk geometry overflow");
    return kExitUsage;
  }
  size_t stripe_size = 0;
  if (!FitsSize(stripe_bytes, &stripe_size)) {
    Error("chunk geometry overflow");
    return kExitUsage;
  }
  const uint16_t n = static_cast<uint16_t>(options.k + options.r);
  LCHEncoder encoder(options.k, options.r);
  const size_t workspace_bytes = encoder.WorkspaceSize(options.chunk_size);
  if (workspace_bytes == static_cast<size_t>(-1)) {
    Error("workspace overflow");
    return kExitUsage;
  }

  EncodeState state;
  state.encoder = &encoder;
  state.geometry.k = options.k;
  state.geometry.r = options.r;
  state.geometry.c = options.chunk_size;
  state.geometry.stripe_bytes = stripe_bytes;
  state.geometry.payload_offset = kHeaderSize + name.size();
  state.geometry.filename = name;
  state.share_fds.resize(n);
  for (uint16_t i = 0; i < n; ++i) {
    const fs::path path = out_dir / ShareFileName(name, n, i);
    std::error_code exists;
    if (fs::exists(path, exists) && !options.force) {
      Error("share exists (use --force)");
      return kExitUsage;
    }
    UniqueFd fd(open(path.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0644));
    uint8_t placeholder[kHeaderSize]{};
    if (!fd.valid() || !PWriteAll(fd.get(), placeholder, kHeaderSize, 0) ||
        !PWriteAll(fd.get(), name.data(), name.size(), kHeaderSize)) {
      Error("failed to create share");
      return kExitUsage;
    }
    state.share_fds[i] = std::move(fd);
  }

  StreamingHash file_hash;
  if (!file_hash.Initialize()) {
    Error("failed to initialize hash");
    return kExitUsage;
  }
  Verbose(options.verbose, "encode k=", options.k, " r=", options.r,
          " chunk=", options.chunk_size, " jobs=", options.jobs,
          " backend=",
          BackendName(gf2p8::lch::SelectBackend(options.chunk_size)),
          " input=", read_stdin ? "stdin" : options.input.c_str());
  ProgressDisplay progress("Encoding", options.progress_mode, !read_stdin,
                           source_size);

  auto consume = [&](uint32_t stripe, size_t live, AlignedBuffer* data,
                      AlignedBuffer* work, std::span<Element> workspace) {
    if (!ProcessEncodeStripe(&state, stripe, live, data, work, workspace)) {
      SetError(&state.error, kExitUsage);
      return false;
    }
    progress.Add(live);
    return true;
  };

  uint32_t stripe_count = 0;
  uint64_t original_size = 0;
  auto read_next = [&](AlignedBuffer* data, uint64_t* live) {
    if (!data->Allocate(stripe_size)) {
      return false;
    }
    size_t got = 0;
    if (!ReadUpTo(source.get(), data->data(), stripe_size, &got)) {
      return false;
    }
    if (got < stripe_size) {
      std::memset(data->data() + got, 0, stripe_size - got);
    }
    *live = got;
    return true;
  };

  AlignedBuffer work;
  AlignedBuffer workspace_buffer;
  const bool serial = options.jobs == 1;
  if (serial && (!work.Allocate(static_cast<size_t>(n) * options.chunk_size) ||
                  !workspace_buffer.Allocate(workspace_bytes))) {
    progress.Finish(false);
    Error("out of memory");
    return kExitUsage;
  }

  StripeQueue queue(options.jobs);
  std::vector<std::thread> workers;
  if (!serial) {
    workers.reserve(options.jobs);
    for (uint32_t i = 0; i < options.jobs; ++i) {
      workers.emplace_back([&] {
        AlignedBuffer local_work;
        AlignedBuffer local_workspace;
        if (!local_work.Allocate(static_cast<size_t>(n) * options.chunk_size) ||
            !local_workspace.Allocate(workspace_bytes)) {
          SetError(&state.error, kExitUsage);
          queue.Stop();
          return;
        }
        uint32_t stripe = 0;
        AlignedBuffer data;
        size_t live = 0;
        while (queue.Pop(&stripe, &data, &live)) {
          if (!consume(stripe, live, &data, &local_work,
                       std::span<Element>(local_workspace.data(),
                                          workspace_bytes))) {
            queue.Stop();
            return;
          }
        }
      });
    }
  }

  bool ok = true;
  while (ok) {
    AlignedBuffer data;
    uint64_t live = 0;
    if (!read_next(&data, &live)) {
      ok = false;
      break;
    }
    if (live == 0) {
      break;
    }
    if (stripe_count == UINT32_MAX) {
      ok = false;
      break;
    }
    if (!file_hash.Update(data.data(), static_cast<size_t>(live))) {
      ok = false;
      break;
    }
    original_size += live;
    {
      std::lock_guard lock(state.hash_mutex);
      state.chunk_hashes.resize((static_cast<size_t>(stripe_count) + 1) * n);
    }
    if (serial) {
      if (!consume(
              stripe_count, static_cast<size_t>(live), &data, &work,
              std::span<Element>(workspace_buffer.data(), workspace_bytes))) {
        ok = false;
        break;
      }
    } else if (!queue.Push(stripe_count, std::move(data),
                           static_cast<size_t>(live))) {
      ok = false;
      break;
    }
    ++stripe_count;
  }

  if (!serial) {
    if (!ok) {
      SetError(&state.error, kExitUsage);
      queue.Stop();
    } else {
      queue.Done();
    }
    for (auto& worker : workers) {
      worker.join();
    }
  }
  if (!ok || state.error.load() != 0) {
    progress.Finish(false);
    Error("failed to encode");
    return kExitUsage;
  }

  state.geometry.m = stripe_count;
  state.geometry.original_size = original_size;
  state.geometry.whole_file_hash = file_hash.Digest();
  progress.Finish(true);
  return FinishEncode(&state);
}

bool LoadShare(const std::string& path,
               Geometry* geometry,
               std::vector<AcceptedShare>* shares) {
  UniqueFd fd(open(path.c_str(), O_RDONLY));
  if (!fd.valid()) {
    return false;
  }
  uint8_t packed[kHeaderSize];
  size_t got = 0;
  if (!PReadRange(fd.get(), packed, kHeaderSize, 0, &got) ||
      got != kHeaderSize) {
    return false;
  }
  Header header;
  if (!UnpackHeader(packed, &header) || header.k == 0 || header.r == 0 ||
      header.c == 0 || header.filename_len == 0 ||
      header.filename_len > kMaxFilenameLen) {
    return false;
  }
  const uint32_t n = static_cast<uint32_t>(header.k) + header.r;
  if (n > 0xffff || header.share_id >= n) {
    return false;
  }
  std::string filename(header.filename_len, '\0');
  if (!PReadRange(fd.get(), filename.data(), header.filename_len, kHeaderSize,
                  &got) ||
      got != header.filename_len || !FilenameValid(filename)) {
    return false;
  }
  uint64_t stripe_bytes = 0;
  if (MulOverflow(header.k, header.c, &stripe_bytes)) {
    return false;
  }
  Geometry candidate;
  candidate.k = header.k;
  candidate.r = header.r;
  candidate.c = header.c;
  candidate.m = header.m;
  candidate.original_size = header.original_size;
  candidate.whole_file_hash = header.whole_file_hash;
  candidate.stripe_bytes = stripe_bytes;
  candidate.payload_offset = kHeaderSize + header.filename_len;
  if (!GeometryValid(candidate)) {
    return false;
  }
  uint64_t payload_bytes = 0;
  uint64_t hash_bytes = 0;
  if (!PayloadBytes(candidate, &payload_bytes) ||
      MulOverflow(header.m, 8, &hash_bytes)) {
    return false;
  }
  size_t hash_size = 0;
  if (!FitsSize(hash_bytes, &hash_size)) {
    return false;
  }
  uint64_t hash_offset = 0;
  if (AddOverflow(candidate.payload_offset, payload_bytes, &hash_offset)) {
    return false;
  }
  std::vector<uint64_t> chunk_hashes(header.m);
  if (header.m != 0) {
    std::vector<uint8_t> table(hash_size);
    if (!PReadRange(fd.get(), table.data(), hash_size, hash_offset, &got) ||
        got != hash_size) {
      return false;
    }
    ReadHashTable(table.data(), chunk_hashes.data(), chunk_hashes.size());
  }
  if (MetadataHash(header, filename, chunk_hashes.data(),
                   chunk_hashes.size()) != header.metadata_hash) {
    return false;
  }

  std::string parsed_base;
  uint16_t parsed_n = 0;
  uint16_t parsed_id = 0;
  if (!ParseShareFileName(fs::path(path).filename().string(), &parsed_base,
                          &parsed_n, &parsed_id) ||
      parsed_id != header.share_id || parsed_n != static_cast<uint16_t>(n)) {
    return false;
  }

  if (shares->empty()) {
    *geometry = std::move(candidate);
    geometry->filename = std::move(filename);
  } else if (header.k != geometry->k || header.r != geometry->r ||
             header.c != geometry->c || header.m != geometry->m ||
             header.original_size != geometry->original_size ||
             header.whole_file_hash != geometry->whole_file_hash ||
             filename != geometry->filename) {
    return false;
  }

  for (const auto& share : *shares) {
    if (share.share_id == header.share_id) {
      return false;
    }
  }
  AcceptedShare accepted;
  accepted.fd = std::move(fd);
  accepted.share_id = header.share_id;
  accepted.chunk_hashes = std::move(chunk_hashes);
  shares->push_back(std::move(accepted));
  return true;
}

bool IngestShares(const std::vector<std::string>& paths,
                  Geometry* geometry,
                  std::vector<AcceptedShare>* shares) {
  for (const auto& path : paths) {
    LoadShare(path, geometry, shares);
  }
  return !shares->empty();
}

enum class StripeHealth { intact, repairable, unrecoverable };

StripeHealth InspectStripe(const Geometry& geometry,
                           const std::vector<AcceptedShare>& shares,
                           uint32_t stripe,
                           std::vector<uint8_t>* present,
                           uint32_t* chunk_size,
                           Element* storage) {
  const uint16_t n = ShareCount(geometry);
  if (!StripeChunkSize(geometry, stripe, chunk_size)) {
    return StripeHealth::unrecoverable;
  }
  std::memset(storage, 0, static_cast<size_t>(n) * *chunk_size);
  present->assign(n, 0);
  size_t present_count = 0;
  size_t present_data = 0;
  uint64_t payload_index = 0;
  MulOverflow(stripe, geometry.c, &payload_index);
  const uint64_t offset = geometry.payload_offset + payload_index;
  for (const auto& share : shares) {
    Element* dest = storage + static_cast<size_t>(share.share_id) * *chunk_size;
    size_t got = 0;
    if (!PReadRange(share.fd.get(), dest, *chunk_size, offset, &got) ||
        got != *chunk_size) {
      std::memset(dest, 0, *chunk_size);
      continue;
    }
    if (HashBytes(dest, *chunk_size) != share.chunk_hashes[stripe]) {
      std::memset(dest, 0, *chunk_size);
      continue;
    }
    (*present)[share.share_id] = 1;
    ++present_count;
    if (share.share_id < geometry.k) {
      ++present_data;
    }
  }
  if (present_count < geometry.k) {
    return StripeHealth::unrecoverable;
  }
  if (present_data == geometry.k) {
    return StripeHealth::intact;
  }
  return StripeHealth::repairable;
}

bool DecodeStripeData(const LCHDecoder& decoder,
                      Element* storage,
                      const std::vector<uint8_t>& present,
                      std::span<Element> workspace,
                      uint32_t chunk_size) {
  const size_t k = decoder.DataCount();
  const size_t r = decoder.RecoveryCount();
  std::vector<Element*> data(k);
  std::vector<const Element*> recovery(r);
  std::vector<uint8_t> data_present(k);
  std::vector<uint8_t> recovery_present(r);
  for (size_t i = 0; i < k; ++i) {
    data[i] = storage + i * chunk_size;
    data_present[i] = present[i];
  }
  for (size_t i = 0; i < r; ++i) {
    recovery_present[i] = present[k + i];
    recovery[i] =
        recovery_present[i] ? storage + (k + i) * chunk_size : nullptr;
  }
  return decoder.Decode(data, data_present, recovery, recovery_present,
                        chunk_size, workspace, Backend::tuned,
                        Radix::radix4) == Status::ok;
}

int RunVerify(const Options& options) {
  Geometry geometry;
  std::vector<AcceptedShare> shares;
  if (!IngestShares(options.shares, &geometry, &shares)) {
    Error("no usable shares");
    return kExitUsage;
  }
  std::atomic<uint32_t> intact{0};
  std::atomic<uint32_t> repairable{0};
  std::atomic<uint32_t> unrecoverable{0};
  std::atomic<uint32_t> next_stripe{0};
  std::atomic<bool> operation_error{false};
  const uint32_t jobs =
      geometry.m == 0 ? 1 : std::min(options.jobs, geometry.m);
  Verbose(options.verbose, "verify k=", geometry.k, " r=", geometry.r,
          " chunk=", geometry.c, " stripes=", geometry.m, " jobs=", jobs,
          " shares=", shares.size(), " backend=",
          BackendName(gf2p8::lch::SelectBackend(geometry.c)));
  ProgressDisplay progress("Verifying", options.progress_mode, true,
                           geometry.original_size);
  auto worker = [&] {
    AlignedBuffer storage;
    if (!storage.Allocate(static_cast<size_t>(ShareCount(geometry)) *
                           geometry.c)) {
      operation_error.store(true);
      return;
    }
    std::vector<uint8_t> present;
    while (!operation_error.load()) {
      const uint32_t stripe = next_stripe.fetch_add(1);
      if (stripe >= geometry.m) {
        break;
      }
      uint32_t chunk_size = 0;
      const StripeHealth health = InspectStripe(
          geometry, shares, stripe, &present, &chunk_size, storage.data());
      if (health == StripeHealth::intact) {
        intact.fetch_add(1);
      } else if (health == StripeHealth::repairable) {
        repairable.fetch_add(1);
      } else {
        unrecoverable.fetch_add(1);
      }
      uint64_t live = 0;
      if (!StripeLiveBytes(geometry, stripe, &live)) {
        operation_error.store(true);
        return;
      }
      progress.Add(live);
    }
  };
  if (jobs <= 1) {
    worker();
  } else {
    std::vector<std::thread> workers;
    for (uint32_t i = 0; i < jobs; ++i) {
      workers.emplace_back(worker);
    }
    for (auto& thread : workers) {
      thread.join();
    }
  }
  if (operation_error.load()) {
    progress.Finish(false);
    Error("failed to verify shares");
    return kExitUsage;
  }
  progress.Finish(true);
  std::fprintf(stdout, "stripes=%u intact=%u repairable=%u unrecoverable=%u\n",
               geometry.m, intact.load(), repairable.load(),
               unrecoverable.load());
  if (unrecoverable.load() != 0) {
    return kExitUnrecoverable;
  }
  if (repairable.load() != 0) {
    return kExitRepairable;
  }
  return kExitOk;
}

int RunDecode(const Options& options) {
  Geometry geometry;
  std::vector<AcceptedShare> shares;
  if (!IngestShares(options.shares, &geometry, &shares)) {
    Error("no usable shares");
    return kExitUsage;
  }
  LCHDecoder decoder(geometry.k, geometry.r);
  if (!decoder.Valid()) {
    Error("unsupported share geometry");
    return kExitUsage;
  }
  const size_t workspace_bytes = decoder.WorkspaceSize(geometry.c);
  if (workspace_bytes == static_cast<size_t>(-1)) {
    Error("workspace overflow");
    return kExitUsage;
  }

  const bool write_stdout =
      options.output == "-" || (!options.have_output && !isatty(STDOUT_FILENO));
  fs::path dest_path;
  UniqueFd dest;
  bool created_dest = false;
  if (write_stdout) {
    dest = UniqueFd(STDOUT_FILENO, false);
  } else {
    dest_path = options.have_output ? fs::path(options.output)
                                    : fs::path(geometry.filename);
    std::error_code exists_error;
    const bool dest_exists = fs::exists(dest_path, exists_error);
    if (dest_exists && !options.force) {
      Error("destination exists (use --force)");
      return kExitUsage;
    }
    created_dest = !dest_exists;
    dest = UniqueFd(open(dest_path.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0644));
    if (!dest.valid() ||
        ftruncate(dest.get(), static_cast<off_t>(geometry.original_size)) !=
            0) {
      Error("failed to create destination");
      return kExitUsage;
    }
  }

  const bool serial = options.jobs == 1 || write_stdout || geometry.m == 0;
  StreamingHash stdout_hash;
  if (write_stdout && !stdout_hash.Initialize()) {
    Error("failed to initialize hash");
    return kExitUsage;
  }
  const uint32_t jobs = serial ? 1 : std::min(options.jobs, geometry.m);
  Verbose(options.verbose, "decode k=", geometry.k, " r=", geometry.r,
          " chunk=", geometry.c, " stripes=", geometry.m, " jobs=", jobs,
          " shares=", shares.size(), " backend=",
          BackendName(gf2p8::lch::SelectBackend(geometry.c)), " output=",
          write_stdout ? "stdout" : dest_path.c_str());
  ProgressDisplay progress("Decoding", options.progress_mode, true,
                           geometry.original_size);
  std::atomic<int> error{0};
  std::atomic<uint32_t> next_stripe{0};
  auto process = [&](uint32_t stripe, Element* storage,
                     std::span<Element> workspace) {
    std::vector<uint8_t> present;
    uint32_t chunk_size = 0;
    const StripeHealth health =
        InspectStripe(geometry, shares, stripe, &present, &chunk_size, storage);
    if (health == StripeHealth::unrecoverable) {
      SetError(&error, kExitUnrecoverable);
      return false;
    }
    if (health == StripeHealth::repairable &&
        !DecodeStripeData(decoder, storage, present, workspace, chunk_size)) {
      SetError(&error, kExitUnrecoverable);
      return false;
    }
    uint64_t live = 0;
    if (!StripeLiveBytes(geometry, stripe, &live)) {
      SetError(&error, kExitUsage);
      return false;
    }
    const uint64_t offset =
        static_cast<uint64_t>(stripe) * geometry.stripe_bytes;
    if (live == 0) {
      return true;
    }
    if (write_stdout) {
      if (!WriteAll(dest.get(), storage, static_cast<size_t>(live))) {
        SetError(&error, kExitUsage);
        return false;
      }
      if (!stdout_hash.Update(storage, static_cast<size_t>(live))) {
        SetError(&error, kExitUsage);
        return false;
      }
    } else if (!PWriteAll(dest.get(), storage, static_cast<size_t>(live),
                           offset)) {
      SetError(&error, kExitUsage);
      return false;
    }
    progress.Add(live);
    return true;
  };

  auto worker_body = [&] {
    AlignedBuffer storage;
    AlignedBuffer workspace_buffer;
    if (!storage.Allocate(static_cast<size_t>(ShareCount(geometry)) *
                          geometry.c) ||
        !workspace_buffer.Allocate(workspace_bytes)) {
      SetError(&error, kExitUsage);
      return;
    }
    while (error.load() == 0) {
      const uint32_t stripe = next_stripe.fetch_add(1);
      if (stripe >= geometry.m) {
        break;
      }
      if (!process(
              stripe, storage.data(),
              std::span<Element>(workspace_buffer.data(), workspace_bytes))) {
        return;
      }
    }
  };

  if (serial) {
    AlignedBuffer storage;
    AlignedBuffer workspace_buffer;
    if (!storage.Allocate(static_cast<size_t>(ShareCount(geometry)) *
                           geometry.c) ||
        !workspace_buffer.Allocate(workspace_bytes)) {
      progress.Finish(false);
      Error("out of memory");
      if (created_dest) {
        dest.Reset();
        fs::remove(dest_path);
      }
      return kExitUsage;
    }
    for (uint32_t stripe = 0; stripe < geometry.m; ++stripe) {
      if (!process(
              stripe, storage.data(),
              std::span<Element>(workspace_buffer.data(), workspace_bytes))) {
        break;
      }
    }
  } else {
    std::vector<std::thread> workers;
    for (uint32_t i = 0; i < jobs; ++i) {
      workers.emplace_back(worker_body);
    }
    for (auto& thread : workers) {
      thread.join();
    }
  }

  if (error.load() != 0) {
    progress.Finish(false);
    if (created_dest) {
      dest.Reset();
      fs::remove(dest_path);
    }
    if (error.load() == kExitUnrecoverable) {
      Error("unrecoverable");
    } else {
      Error("failed to decode");
    }
    return error.load();
  }

  if (write_stdout) {
    const bool hash_ok = stdout_hash.Digest() == geometry.whole_file_hash;
    progress.Finish(hash_ok);
    if (!hash_ok) {
      Error("reconstructed file hash mismatch");
      return kExitUnrecoverable;
    }
    return kExitOk;
  }
  dest.Reset();
  UniqueFd check(open(dest_path.c_str(), O_RDONLY));
  if (!check.valid()) {
    progress.Finish(false);
    Error("failed to hash reconstructed file");
    if (created_dest) {
      fs::remove(dest_path);
    }
    return kExitUnrecoverable;
  }
  StreamingHash file_hash;
  if (!file_hash.Initialize()) {
    progress.Finish(false);
    Error("failed to hash reconstructed file");
    if (created_dest) {
      check.Reset();
      fs::remove(dest_path);
    }
    return kExitUnrecoverable;
  }
  std::vector<uint8_t> buffer(1 << 16);
  uint64_t remaining = geometry.original_size;
  uint64_t offset = 0;
  while (remaining != 0) {
    const size_t want =
        static_cast<size_t>(std::min<uint64_t>(buffer.size(), remaining));
    size_t got = 0;
    if (!PReadRange(check.get(), buffer.data(), want, offset, &got) ||
        got != want) {
      progress.Finish(false);
      Error("reconstructed file hash mismatch");
      if (created_dest) {
        check.Reset();
        fs::remove(dest_path);
      }
      return kExitUnrecoverable;
    }
    file_hash.Update(buffer.data(), want);
    offset += want;
    remaining -= want;
  }
  const bool hash_ok = file_hash.Digest() == geometry.whole_file_hash;
  progress.Finish(hash_ok);
  if (!hash_ok) {
    Error("reconstructed file hash mismatch");
    if (created_dest) {
      check.Reset();
      fs::remove(dest_path);
    }
    return kExitUnrecoverable;
  }
  return kExitOk;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    PrintHelp(stderr);
    return kExitUsage;
  }
  const std::string_view first = argv[1];
  if (first == "-h" || first == "--help") {
    PrintHelp(stdout);
    return kExitOk;
  }
  Options options;
  const int parsed = ParseCommandOptions(argc - 1, argv + 1, &options);
  if (parsed != kExitOk) {
    return parsed;
  }
  if (options.command == "encode") {
    return RunEncode(options);
  }
  if (options.command == "verify") {
    return RunVerify(options);
  }
  return RunDecode(options);
}
