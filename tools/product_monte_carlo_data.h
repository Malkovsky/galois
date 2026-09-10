#pragma once

#include <array>
#include <bit>
#include <charconv>
#include <cstdint>
#include <map>
#include <nlohmann/json.hpp>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mc {
using Json = nlohmann::json;
inline constexpr std::string_view kCode =
    "RS256,224 x RS256,254 Cantor systematic row major";
inline constexpr std::string_view kSnapshots = "atomic summary v1";
inline constexpr std::string_view kFloyd =
    "splitmix64 domain seeds; mt19937_64; rejection modulo; Floyd complement "
    "v1";
inline constexpr std::string_view kFisherYates =
    "splitmix64 domain seeds; mt19937_64; rejection modulo; persistent "
    "Fisher-Yates complement v1; replay saved flips";
inline constexpr std::array<const char*, 22> kMetrics{
    "initial full block corrupted bits",
    "initial full block corrupted bytes",
    "initial information corrupted bits",
    "initial information corrupted bytes",
    "residual full block bits",
    "residual full block bytes",
    "residual information bits",
    "residual information bytes",
    "message failures",
    "full block failures",
    "zero syndrome outcomes",
    "zero syndrome wrong full blocks",
    "directional passes",
    "accepted bit changes",
    "accepted byte changes",
    "strong accepted bit changes",
    "strong accepted byte changes",
    "weak accepted bit changes",
    "weak accepted byte changes",
    "strong lines visited",
    "weak lines visited",
    "pass limit outcomes"};

inline void Require(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

inline Json Parse(const std::string& text) {
  std::vector<std::set<std::string>> objects;
  try {
    return Json::parse(
        text, [&](int depth, Json::parse_event_t event, Json& value) {
          using E = Json::parse_event_t;
          if (event == E::object_start || event == E::array_start) {
            Require(depth < 32, "JSON nesting exceeds 32");
          }
          if (event == E::object_start) {
            objects.emplace_back();
          }
          if (event == E::object_end) {
            objects.pop_back();
          }
          if (event == E::key) {
            Require(objects.back().insert(value.get<std::string>()).second,
                    "duplicate JSON key");
          }
          // Overflowing integer tokens also become floats in the DOM. Never
          // accept their rounded values, including unused legacy squared sums.
          Require(!value.is_number_float(),
                  "noninteger JSON number or integer exceeds uint64 (including "
                  "legacy moments)");
          return true;
        });
  } catch (const Json::out_of_range&) {
    throw std::runtime_error(
        "JSON number exceeds uint64 (including legacy moments)");
  }
}

inline std::string Dump(const Json& value, bool pretty = false) {
  // The default ordered object map and ASCII escaping match Python's
  // sort_keys=True, ensure_ascii=True, separators=(",", ":") identities.
  return value.dump(pretty ? 2 : -1, ' ', true);
}

inline uint64_t Natural(const Json& value) {
  Require(value.is_number_unsigned() ||
              (value.is_number_integer() && value.get<int64_t>() >= 0),
          "expected unsigned integer");
  return value.get<uint64_t>();
}
inline uint64_t U64(const Json& value) {
  return Natural(value);
}
inline uint64_t Decimal(std::string_view text) {
  Require(!text.empty(), "expected unsigned decimal integer");
  uint64_t value = 0;
  const auto [end, error] =
      std::from_chars(text.data(), text.data() + text.size(), value);
  Require(error != std::errc::result_out_of_range, "integer exceeds uint64");
  Require(error == std::errc{} && end == text.data() + text.size(),
          "expected unsigned decimal integer");
  return value;
}
inline uint64_t CheckedAdd(uint64_t a, uint64_t b) {
  Require(b <= UINT64_MAX - a, "uint64 counter addition overflow");
  return a + b;
}
inline uint64_t CheckedMultiply(uint64_t a, uint64_t b) {
  Require(b == 0 || a <= UINT64_MAX / b,
          "uint64 counter multiplication overflow");
  return a * b;
}
inline Json Number(uint64_t value) {
  return Json(value);
}
inline void Fields(const Json& value, std::set<std::string> expected) {
  Require(value.is_object() && value.size() == expected.size(),
          "incompatible object fields");
  for (const auto& item : value.items()) {
    Require(expected.erase(std::string(item.key())) == 1,
            "incompatible object field");
  }
}

struct Settings {
  uint64_t n1 = 256, k1 = 224, n2 = 256, k2 = 254;
  uint64_t seed = 0, size = 1000, batches = 0, threads = 1, lo = 2500,
           hi = 2700, passes = 16, checkpoint = 64, report = 2, sync = 5;
  bool anchors = true, binary = true;

  /** @brief Transmitted bits per block after dimension validation. */
  uint64_t FullBits() const { return 8 * n1 * n2; }
  /** @brief Information bits per block after dimension validation. */
  uint64_t InfoBits() const { return 8 * k1 * k2; }
  /** @brief Whether legacy implicit dimensions describe this code. */
  bool DefaultDimensions() const {
    return n1 == 256 && k1 == 224 && n2 == 256 && k2 == 254;
  }
  /** @brief Code and coordinate convention, preserving the legacy spelling. */
  std::string Code() const {
    return "RS" + std::to_string(n1) + "," + std::to_string(k1) + " x RS" +
           std::to_string(n2) + "," + std::to_string(k2) +
           " Cantor systematic row major";
  }

  void Validate() const {
    Require(n1 <= 256 && std::has_single_bit(n1) && k1 < n1 && n1 - k1 >= 2 &&
                n1 - k1 <= k1 && std::has_single_bit(n1 - k1) && n2 <= 256 &&
                k2 >= 2 && k2 < n2 &&
                (n2 - k2 == 2 || (n2 == 256 && k2 == 252)),
            "invalid dimensions: strong N,R powers of two, N<=256, 2<=R<=K; "
            "weak N<=256, K>=2, R=2 (shortening supported), or RS(256,252)");
    Require(lo <= hi && hi <= FullBits(),
            "require 0 <= minimum <= maximum <= 8*n1*n2; set smaller explicit "
            "flip bounds for small codes (defaults 2500..2700 are not capped)");
    Require(size > 0, "batch size must be positive");
    Require(threads >= 1 && threads <= 1024, "threads must be in [1,1024]");
    Require(passes >= 2 && passes <= 1000000,
            "pass cap must be in [2,1000000]");
    Require(checkpoint >= 1 && checkpoint <= 4096,
            "checkpoint trials must be in [1,4096]");
    Require(report >= 1 && report <= 86400 && sync >= 1 && sync <= 86400,
            "report/fsync seconds must be in [1,86400]");
  }
  Json ToJson() const {
    Json j;
    j["root seed"] = seed;
    j["batch size"] = size;
    j["batches"] = batches;
    j["threads"] = threads;
    j["minimum flipped bits"] = lo;
    j["maximum flipped bits"] = hi;
    j["maximum directional passes"] = passes;
    j["checkpoint trials"] = checkpoint;
    j["report seconds"] = report;
    j["fsync seconds"] = sync;
    j["anchors"] = anchors;
    j["binary image"] = binary;
    if (!DefaultDimensions()) {
      j["n1"] = n1;
      j["k1"] = k1;
      j["n2"] = n2;
      j["k2"] = k2;
    }
    return j;
  }
  static Settings FromJson(const Json& j) {
    std::set<std::string> fields{"root seed",
                                 "batch size",
                                 "batches",
                                 "threads",
                                 "minimum flipped bits",
                                 "maximum flipped bits",
                                 "maximum directional passes",
                                 "checkpoint trials",
                                 "report seconds",
                                 "fsync seconds",
                                 "anchors",
                                 "binary image"};
    const bool dimensions = j.contains("n1") || j.contains("k1") ||
                            j.contains("n2") || j.contains("k2");
    if (dimensions) {
      fields.insert({"n1", "k1", "n2", "k2"});
    }
    Fields(j, fields);
    Settings s;
    if (dimensions) {
      s.n1 = U64(j.at("n1"));
      s.k1 = U64(j.at("k1"));
      s.n2 = U64(j.at("n2"));
      s.k2 = U64(j.at("k2"));
    }
    s.seed = U64(j.at("root seed"));
    s.size = U64(j.at("batch size"));
    s.batches = U64(j.at("batches"));
    s.threads = U64(j.at("threads"));
    s.lo = U64(j.at("minimum flipped bits"));
    s.hi = U64(j.at("maximum flipped bits"));
    s.passes = U64(j.at("maximum directional passes"));
    s.checkpoint = U64(j.at("checkpoint trials"));
    s.report = U64(j.at("report seconds"));
    s.sync = U64(j.at("fsync seconds"));
    Require(j.at("anchors").is_boolean() && j.at("binary image").is_boolean(),
            "gates must be boolean");
    s.anchors = j.at("anchors").get<bool>();
    s.binary = j.at("binary image").get<bool>();
    s.Validate();
    return s;
  }
};

struct Stats {
  uint64_t blocks = 0, iterations = 0, info_raw = 0, info_post = 0,
           full_raw = 0, full_post = 0;
  void Add(const std::array<uint64_t, 22>& m, uint64_t full_bits = 524288) {
    Add(Stats{1, m[12], m[2], m[6], m[0], m[4]}, full_bits);
  }
  void Add(const Stats& s, uint64_t full_bits = 524288) {
    Stats next{
        CheckedAdd(blocks, s.blocks),     CheckedAdd(iterations, s.iterations),
        CheckedAdd(info_raw, s.info_raw), CheckedAdd(info_post, s.info_post),
        CheckedAdd(full_raw, s.full_raw), CheckedAdd(full_post, s.full_post)};
    CheckedMultiply(next.blocks, full_bits);
    *this = next;
  }
  Json ToJson(const Settings& settings = {}) const {
    Json j;
    j["completed blocks"] = Number(blocks);
    j["total iterations"] = Number(iterations);
    for (bool info : {true, false}) {
      Json bits;
      bits["total bits"] = Number(CheckedMultiply(
          blocks, info ? settings.InfoBits() : settings.FullBits()));
      bits["raw corrupted bits"] = Number(info ? info_raw : full_raw);
      bits["post decoding corrupted bits"] =
          Number(info ? info_post : full_post);
      j[info ? "information bits" : "full-codeword bits"] = std::move(bits);
    }
    return j;
  }
  static Stats FromJson(const Json& j,
                        uint64_t count,
                        uint64_t k,
                        uint64_t passes,
                        const Settings& settings = {}) {
    Fields(j, {"completed blocks", "total iterations", "information bits",
               "full-codeword bits"});
    Stats s;
    s.blocks = Natural(j.at("completed blocks"));
    s.iterations = Natural(j.at("total iterations"));
    // Compare the upper bound by division: count*passes need not fit u64
    // when the actual iteration count does.
    Require(s.blocks == count && s.iterations >= CheckedMultiply(count, 2) &&
                (passes != 0 && s.iterations / passes <= count &&
                 (s.iterations / passes < count || s.iterations % passes == 0)),
            "inconsistent block/iteration counts");
    for (bool info : {true, false}) {
      const auto& bits = j.at(info ? "information bits" : "full-codeword bits");
      Fields(bits, {"total bits", "raw corrupted bits",
                    "post decoding corrupted bits"});
      const auto total = CheckedMultiply(
          count, info ? settings.InfoBits() : settings.FullBits());
      auto raw = Natural(bits.at("raw corrupted bits"));
      auto post = Natural(bits.at("post decoding corrupted bits"));
      Require(Natural(bits.at("total bits")) == total && raw <= total &&
                  post <= total,
              "invalid bit counters");
      (info ? s.info_raw : s.full_raw) = raw;
      (info ? s.info_post : s.full_post) = post;
    }
    Require(s.full_raw == CheckedMultiply(count, k),
            "initial channel is not exact k");
    Require(s.info_raw <= s.full_raw && s.info_post <= s.full_post &&
                s.full_raw - s.info_raw <=
                    CheckedMultiply(
                        count, settings.FullBits() - settings.InfoBits()) &&
                s.full_post - s.info_post <=
                    CheckedMultiply(count,
                                    settings.FullBits() - settings.InfoBits()),
            "inconsistent information/full bit counters");
    return s;
  }
};

inline Json LegacyStats() {
  Json j;
  for (auto name : kMetrics) {
    j[name]["sum"] = uint64_t{0};
    j[name]["squared sum"] = uint64_t{0};
  }
  return j;
}
inline void AddLegacy(Json& to, const Json& from) {
  Json next = to;
  for (auto name : kMetrics) {
    for (auto field : {"sum", "squared sum"}) {
      next[name][field] = Number(CheckedAdd(Natural(to.at(name).at(field)),
                                            Natural(from.at(name).at(field))));
    }
  }
  to.swap(next);
}
inline void AddLegacyTrial(Json& to, const std::array<uint64_t, 22>& m) {
  Json next = to;
  for (size_t i = 0; i < m.size(); ++i) {
    auto& item = next.at(kMetrics[i]);
    item["sum"] = Number(CheckedAdd(Natural(item.at("sum")), m[i]));
    item["squared sum"] = Number(CheckedAdd(Natural(item.at("squared sum")),
                                            CheckedMultiply(m[i], m[i])));
  }
  to.swap(next);
}
inline void ValidateLegacy(const Json& j,
                           uint64_t count,
                           uint64_t k,
                           uint64_t passes,
                           const Settings& settings = {}) {
  Fields(j, std::set<std::string>(kMetrics.begin(), kMetrics.end()));
  Require(passes >= 2 && passes <= 1000000, "invalid legacy pass cap");
  const auto full = settings.FullBits(), info = settings.InfoBits();
  std::array<uint64_t, 22> bounds{full,
                                  full / 8,
                                  info,
                                  info / 8,
                                  full,
                                  full / 8,
                                  info,
                                  info / 8,
                                  1,
                                  1,
                                  1,
                                  1,
                                  passes,
                                  passes * full,
                                  passes * (full / 8),
                                  passes * full,
                                  passes * (full / 8),
                                  passes * full,
                                  passes * (full / 8),
                                  passes * settings.n2,
                                  passes * settings.n1,
                                  1};
  for (size_t i = 0; i < kMetrics.size(); ++i) {
    const auto& item = j.at(kMetrics[i]);
    Fields(item, {"sum", "squared sum"});
    auto sum = Natural(item.at("sum")),
         square = Natural(item.at("squared sum"));
    // Only comparison products are widened, never stored counters. Two u64
    // factors fit exactly; the redundant count*bound*bound bound is implied
    // by sum <= count*bound and square <= bound*sum.
    using Wide = __uint128_t;
    Require(sum <= Wide(count) * bounds[i] &&
                Wide(sum) * sum <= Wide(count) * square && sum <= square &&
                square <= Wide(bounds[i]) * sum,
            "inconsistent legacy moments");
  }
  Require(Natural(j.at(kMetrics[0]).at("sum")) == CheckedMultiply(count, k) &&
              Natural(j.at(kMetrics[0]).at("squared sum")) ==
                  CheckedMultiply(CheckedMultiply(count, k), k),
          "initial channel is not exact k");
  for (size_t i : {13, 14}) {
    Require(Natural(j.at(kMetrics[i]).at("sum")) ==
                CheckedAdd(Natural(j.at(kMetrics[i + 2]).at("sum")),
                           Natural(j.at(kMetrics[i + 4]).at("sum"))),
            "directional accepted totals disagree");
  }
}

struct Aggregate {
  std::string identity;
  unsigned schema;
  Settings settings;
  Stats overall;
  std::map<uint64_t, Stats> by_k;
  Json legacy = LegacyStats();
  std::map<uint64_t, Json> legacy_k;
  /** @brief Initialize an empty aggregate for a run and schema revision. */
  Aggregate(std::string id, unsigned revision, Settings dimensions = {})
      : identity(std::move(id)), schema(revision), settings(dimensions) {}

  struct PreparedAdd {
    Stats overall;
    decltype(by_k)::node_type row;
    Json legacy;
    decltype(legacy_k)::node_type legacy_row;
  };

  /**
   * @brief Check totals and allocate staged nodes without changing this
   * aggregate.
   * @param k Flipped-bit count of the updated stratum.
   * @param stats Counter increment.
   * @param old Legacy moment increment, required for schema 1.
   * @return Prepared update that can be discarded without side effects.
   */
  PreparedAdd PrepareAdd(uint64_t k,
                         const Stats& stats,
                         const Json& old = Json()) const {
    PreparedAdd next;
    next.overall = overall;
    next.overall.Add(stats, settings.FullBits());
    const auto it = by_k.find(k);
    Stats row = it == by_k.end() ? Stats{} : it->second;
    row.Add(stats, settings.FullBits());
    // Stage just one node, with the same allocator as the destination map.
    decltype(by_k) rows;
    rows.emplace(k, row);
    next.row = rows.extract(rows.begin());
    if (schema == 1) {
      next.legacy = legacy;
      AddLegacy(next.legacy, old);
      const auto old_it = legacy_k.find(k);
      Json old_row = old_it == legacy_k.end() ? LegacyStats() : old_it->second;
      AddLegacy(old_row, old);
      decltype(legacy_k) old_rows;
      old_rows.emplace(k, std::move(old_row));
      next.legacy_row = old_rows.extract(old_rows.begin());
    }
    return next;
  }

  /**
   * @brief Publish prepared totals without allocation or exceptions.
   * @param next Update prepared by this aggregate; commit once, with no
   * intervening updates. Node insertion allocates nothing, and the integer
   * comparator cannot throw.
   */
  void Commit(PreparedAdd&& next) noexcept {
    auto row = by_k.insert(std::move(next.row));
    if (!row.inserted) {
      row.position->second = row.node.mapped();
    }
    if (schema == 1) {
      auto old_row = legacy_k.insert(std::move(next.legacy_row));
      if (!old_row.inserted) {
        old_row.position->second.swap(old_row.node.mapped());
      }
      legacy.swap(next.legacy);
    }
    overall = next.overall;
  }

  /** @brief Prepare and commit one increment with strong exception safety. */
  void Add(uint64_t k, const Stats& stats, const Json& old = Json()) {
    Commit(PrepareAdd(k, stats, old));
  }

  /** @brief Serialize the overall totals and all strata in the run's schema. */
  Json Summary() const {
    Json out;
    out["schema revision"] = schema;
    out["run identity"] = identity;
    if (!settings.DefaultDimensions()) {
      out["code parameters"] = {{"n1", settings.n1},
                                {"k1", settings.k1},
                                {"n2", settings.n2},
                                {"k2", settings.k2}};
    }
    Json rows = Json::array();
    for (const auto& [k, s] : by_k) {
      Json row;
      row["flipped bit count"] = k;
      if (schema == 1) {
        row["trial count"] = Number(s.blocks);
        row["statistics"] = legacy_k.at(k);
      } else {
        row["statistics"] = s.ToJson(settings);
      }
      rows.push_back(std::move(row));
    }
    out["by flipped bit count"] = std::move(rows);
    if (schema == 1) {
      out["overall"]["trial count"] = Number(overall.blocks);
      out["overall"]["statistics"] = legacy;
    } else {
      out["overall"]["statistics"] = overall.ToJson(settings);
    }
    return out;
  }
};
}  // namespace mc
