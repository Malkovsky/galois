#pragma once

#include <array>
#include <cstdint>
#include <jsoncons/bigint.hpp>
#include <jsoncons/json.hpp>
#include <jsoncons/json_cursor.hpp>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace mc {
using Json = jsoncons::json;
using Big = jsoncons::bigint;
inline constexpr std::string_view kCode =
    "RS256,224 x RS256,254 Cantor systematic row major";
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
  const auto options = jsoncons::json_options()
                           .lossless_number(true)
                           .max_nesting_depth(32)
                           .err_handler(jsoncons::strict_json_parsing());
  // Validate the library's event stream before DOM construction can discard
  // duplicate keys. Number tokens larger than uint64 retain their exact digits.
  std::vector<std::set<std::string>> objects;
  jsoncons::json_string_cursor cursor(text, options);
  for (; !cursor.done(); cursor.next()) {
    const auto& event = cursor.current();
    using E = jsoncons::staj_event_type;
    if (event.event_type() == E::begin_object) {
      objects.emplace_back();
    }
    if (event.event_type() == E::end_object) {
      objects.pop_back();
    }
    if (event.event_type() == E::key) {
      Require(objects.back().insert(event.get<std::string>()).second,
              "duplicate JSON key");
    }
    Require(event.event_type() != E::double_value &&
                event.tag() != jsoncons::semantic_tag::bigdec,
            "noninteger JSON number");
  }
  return Json::parse(text, options);
}

inline std::string Dump(const Json& value, bool pretty = false) {
  std::string text;
  auto options = jsoncons::json_options()
                     .bigint_format(jsoncons::bigint_chars_format::number)
                     .escape_all_non_ascii(true)
                     .indent_size(2);
  if (!pretty) {
    options.spaces_around_colon(jsoncons::spaces_option::no_spaces)
        .spaces_around_comma(jsoncons::spaces_option::no_spaces);
  }
  value.dump(
      text, options,
      pretty ? jsoncons::indenting::indent : jsoncons::indenting::no_indent);
  return text;
}

inline Big Natural(const Json& value) {
  if (value.is_uint64()) {
    return Big(value.as<uint64_t>());
  }
  Require(value.tag() == jsoncons::semantic_tag::bigint,
          "expected unsigned integer");
  auto n = Big::from_string(value.as<std::string>());
  Require(n >= 0, "expected unsigned integer");
  return n;
}
inline uint64_t U64(const Json& value) {
  auto n = Natural(value);
  Require(n <= Big(UINT64_MAX), "integer exceeds uint64");
  return static_cast<uint64_t>(n);
}
inline Json Number(const Big& value) {
  if (value <= Big(UINT64_MAX)) {
    return Json(static_cast<uint64_t>(value));
  }
  return Json(value.to_string(), jsoncons::semantic_tag::bigint);
}
inline void Fields(const Json& value, std::set<std::string> expected) {
  Require(value.is_object() && value.size() == expected.size(),
          "incompatible object fields");
  for (const auto& item : value.object_range()) {
    Require(expected.erase(std::string(item.key())) == 1,
            "incompatible object field");
  }
}

struct Settings {
  uint64_t seed = 0, size = 1000, batches = 0, threads = 1, lo = 2500,
           hi = 2700, passes = 16, checkpoint = 64, report = 2, sync = 5;
  bool anchors = true, binary = true;

  void Validate() const {
    Require(lo <= hi && hi <= 524288,
            "require 0 <= minimum <= maximum <= 524288");
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
    return j;
  }
  static Settings FromJson(const Json& j) {
    Fields(j, {"root seed", "batch size", "batches", "threads",
               "minimum flipped bits", "maximum flipped bits",
               "maximum directional passes", "checkpoint trials",
               "report seconds", "fsync seconds", "anchors", "binary image"});
    Settings s;
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
    Require(j.at("anchors").is_bool() && j.at("binary image").is_bool(),
            "gates must be boolean");
    s.anchors = j.at("anchors").as<bool>();
    s.binary = j.at("binary image").as<bool>();
    s.Validate();
    return s;
  }
};

struct Stats {
  // All public accumulators are arbitrary precision. Trial metrics remain
  // bounded uint64 values in the private ABI and legacy flip records.
  Big blocks = 0, iterations = 0, info_raw = 0, info_post = 0, full_raw = 0,
      full_post = 0;
  void Add(const std::array<uint64_t, 22>& m) {
    ++blocks;
    iterations += m[12];
    info_raw += m[2];
    info_post += m[6];
    full_raw += m[0];
    full_post += m[4];
  }
  void Add(const Stats& s) {
    blocks += s.blocks;
    iterations += s.iterations;
    info_raw += s.info_raw;
    info_post += s.info_post;
    full_raw += s.full_raw;
    full_post += s.full_post;
  }
  Json ToJson() const {
    Json j;
    j["completed blocks"] = Number(blocks);
    j["total iterations"] = Number(iterations);
    for (bool info : {true, false}) {
      Json bits;
      bits["total bits"] = Number(blocks * (info ? 455168 : 524288));
      bits["raw corrupted bits"] = Number(info ? info_raw : full_raw);
      bits["post decoding corrupted bits"] =
          Number(info ? info_post : full_post);
      j[info ? "information bits" : "full-codeword bits"] = std::move(bits);
    }
    return j;
  }
  static Stats FromJson(const Json& j,
                        const Big& count,
                        uint64_t k,
                        uint64_t passes) {
    Fields(j, {"completed blocks", "total iterations", "information bits",
               "full-codeword bits"});
    Stats s;
    s.blocks = Natural(j.at("completed blocks"));
    s.iterations = Natural(j.at("total iterations"));
    Require(s.blocks == count && s.iterations >= count * 2 &&
                s.iterations <= count * Big(passes),
            "inconsistent block/iteration counts");
    for (bool info : {true, false}) {
      const auto& bits = j.at(info ? "information bits" : "full-codeword bits");
      Fields(bits, {"total bits", "raw corrupted bits",
                    "post decoding corrupted bits"});
      const auto total = count * (info ? 455168 : 524288);
      auto raw = Natural(bits.at("raw corrupted bits"));
      auto post = Natural(bits.at("post decoding corrupted bits"));
      Require(Natural(bits.at("total bits")) == total && raw <= total &&
                  post <= total,
              "invalid bit counters");
      (info ? s.info_raw : s.full_raw) = raw;
      (info ? s.info_post : s.full_post) = post;
    }
    Require(s.full_raw == count * Big(k), "initial channel is not exact k");
    Require(s.info_raw <= s.full_raw && s.info_post <= s.full_post &&
                s.full_raw - s.info_raw <= count * 69120 &&
                s.full_post - s.info_post <= count * 69120,
            "inconsistent information/full bit counters");
    return s;
  }
};

inline Json LegacyStats() {
  Json j;
  for (auto name : kMetrics) {
    j[name]["sum"] = 0;
    j[name]["squared sum"] = 0;
  }
  return j;
}
inline void AddLegacy(Json& to, const Json& from) {
  for (auto name : kMetrics) {
    for (auto field : {"sum", "squared sum"}) {
      to[name][field] = Number(Natural(to.at(name).at(field)) +
                               Natural(from.at(name).at(field)));
    }
  }
}
inline void AddLegacyTrial(Json& to, const std::array<uint64_t, 22>& m) {
  for (size_t i = 0; i < m.size(); ++i) {
    auto& item = to.at(kMetrics[i]);
    item["sum"] = Number(Natural(item.at("sum")) + Big(m[i]));
    item["squared sum"] =
        Number(Natural(item.at("squared sum")) + Big(m[i]) * Big(m[i]));
  }
}
inline void ValidateLegacy(const Json& j,
                           uint64_t count,
                           uint64_t k,
                           uint64_t passes) {
  Fields(j, std::set<std::string>(kMetrics.begin(), kMetrics.end()));
  std::array<uint64_t, 22> bounds{524288,
                                  65536,
                                  455168,
                                  56896,
                                  524288,
                                  65536,
                                  455168,
                                  56896,
                                  1,
                                  1,
                                  1,
                                  1,
                                  passes,
                                  passes * 524288,
                                  passes * 524288,
                                  passes * 524288,
                                  passes * 524288,
                                  passes * 524288,
                                  passes * 524288,
                                  passes * 524288,
                                  passes * 524288,
                                  1};
  for (size_t i = 0; i < kMetrics.size(); ++i) {
    const auto& item = j.at(kMetrics[i]);
    Fields(item, {"sum", "squared sum"});
    auto sum = Natural(item.at("sum")),
         square = Natural(item.at("squared sum"));
    Require(sum <= Big(count) * Big(bounds[i]) &&
                square <= Big(count) * Big(bounds[i]) * Big(bounds[i]) &&
                sum * sum <= Big(count) * square && sum <= square &&
                square <= Big(bounds[i]) * sum,
            "inconsistent legacy moments");
  }
  Require(Natural(j.at(kMetrics[0]).at("sum")) == Big(count) * Big(k) &&
              Natural(j.at(kMetrics[0]).at("squared sum")) ==
                  Big(count) * Big(k) * Big(k),
          "initial channel is not exact k");
  for (size_t i : {13, 14}) {
    Require(Natural(j.at(kMetrics[i]).at("sum")) ==
                Natural(j.at(kMetrics[i + 2]).at("sum")) +
                    Natural(j.at(kMetrics[i + 4]).at("sum")),
            "directional accepted totals disagree");
  }
}
}  // namespace mc
