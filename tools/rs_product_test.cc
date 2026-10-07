#include <algorithm>
#include <bit>
#include <charconv>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "reed_solomon/strong_weak_rs_product_code.h"

namespace {
using gf2p8::Element;
using nlohmann::json;
using namespace gf2p8::rs;
using Bytes = std::vector<Element>;

void Require(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void Keys(const json& value, std::initializer_list<std::string> allowed) {
  Require(value.is_object(), "expected JSON object");
  for (auto it = value.begin(); it != value.end(); ++it) {
    Require(
        std::find(allowed.begin(), allowed.end(), it.key()) != allowed.end(),
        "unknown field: " + it.key());
  }
}

uint64_t Integer(const json& value, uint64_t low, uint64_t high) {
  Require(value.is_number_integer() &&
              (value.is_number_unsigned() || value.get<int64_t>() >= 0),
          "expected nonnegative integer");
  auto n = value.get<uint64_t>();
  Require(n >= low && n <= high, "integer out of bounds");
  return n;
}

size_t Number(const json& j,
              const char* key,
              size_t fallback,
              size_t low,
              size_t high) {
  return j.contains(key) ? Integer(j.at(key), low, high) : fallback;
}

std::string Hex(const Bytes& bytes) {
  const char* digits = "0123456789abcdef";
  std::string text;
  for (auto b : bytes) {
    text += digits[b >> 4];
    text += digits[b & 15];
  }
  return text;
}

Bytes Unhex(const json& value, size_t size) {
  Require(value.is_string(), "hex must be a string");
  const auto text = value.get<std::string>();
  Require(text.size() == 2 * size, "incorrect hex length");
  Bytes bytes(size);
  for (size_t i = 0; i < size; ++i) {
    unsigned n = 0;
    const auto* start = text.data() + 2 * i;
    auto result = std::from_chars(start, start + 2, n, 16);
    Require(result.ec == std::errc{} && result.ptr == start + 2,
            "malformed hex");
    bytes[i] = n;
  }
  return bytes;
}

size_t Bits(const Bytes& word) {
  size_t count = 0;
  for (auto b : word) {
    count += std::popcount(b);
  }
  return count;
}

json Weight(const Bytes& word) {
  return {{"symbols", std::count_if(word.begin(), word.end(),
                                    [](auto b) { return b != 0; })},
          {"bits", Bits(word)}};
}

Bytes Encode(const Bytes& message, size_t n) {
  const size_t k = message.size();
  LCHEncoder encoder(k, n - k);
  Require(encoder.Valid(), "invalid component dimensions");
  Bytes word(n), workspace(encoder.WorkspaceSize(1));
  std::copy(message.begin(), message.end(), word.begin());
  std::vector<const Element*> data(k);
  std::vector<Element*> parity(n - k);
  for (size_t i = 0; i < k; ++i) {
    data[i] = &word[i];
  }
  for (size_t i = k; i < n; ++i) {
    parity[i - k] = &word[i];
  }
  Require(encoder.Encode(data, parity, 1, workspace,
                         gf2p8::lch::Backend::scalar) == gf2p8::lch::Status::ok,
          "component encode failed");
  return word;
}

json Run(const json& input) {
  Keys(input, {"strong_n", "strong_k", "weak_n", "weak_k", "options",
               "transmitted_hex", "errors", "row_masks", "column_masks",
               "expected", "expected_correction"});
  const size_t sn = Number(input, "strong_n", 256, 1, 256);
  const size_t sk = Number(input, "strong_k", 224, 1, 255);
  const size_t wn = Number(input, "weak_n", 256, 1, 256);
  const size_t wk = Number(input, "weak_k", 254, 1, 255);
  StrongWeakRSProductCode code(sn, sk, wn, wk);
  Require(code.Valid(), "unsupported product dimensions");
  ProductDecodeOptions options;
  if (input.contains("options")) {
    const auto& o = input.at("options");
    Keys(o, {"max_directional_passes", "use_anchors", "use_binary_image",
             "use_postprocessing"});
    options.max_directional_passes =
        Number(o, "max_directional_passes", 16, 2, 1024);
    for (const auto* key :
         {"use_anchors", "use_binary_image", "use_postprocessing"}) {
      if (o.contains(key)) {
        Require(o.at(key).is_boolean(), "option must be boolean");
      }
    }
    options.use_anchors = o.value("use_anchors", true);
    options.use_binary_image = o.value("use_binary_image", true);
    options.use_postprocessing = o.value("use_postprocessing", false);
  }
  std::string expected = input.value("expected", std::string());
  Require(!input.contains("expected") || expected == "recovered" ||
              expected == "undetected-failure" ||
              expected == "detected-failure",
          "invalid expected outcome");
  Bytes transmitted(code.BlockSize());
  if (input.contains("transmitted_hex")) {
    transmitted = Unhex(input.at("transmitted_hex"), transmitted.size());
  }
  auto encoded = transmitted;
  Require(code.Encode(encoded, gf2p8::lch::Backend::scalar) ==
                  gf2p8::lch::Status::ok &&
              encoded == transmitted,
          "transmitted block is not a valid systematic codeword");
  auto block = transmitted;
  for (const auto* key : {"errors", "row_masks", "column_masks"}) {
    if (!input.contains(key)) {
      continue;
    }
    const auto& masks = input.at(key);
    Require(masks.is_array() && masks.size() <= 65536, "invalid mask array");
    for (const auto& mask : masks) {
      if (std::string(key) == "errors") {
        Keys(mask, {"row", "column", "xor_hex"});
        size_t row = Integer(mask.at("row"), 0, sn - 1);
        size_t col = Integer(mask.at("column"), 0, wn - 1);
        block[row * wn + col] ^= Unhex(mask.at("xor_hex"), 1)[0];
      } else {
        const bool row = std::string(key) == "row_masks";
        Keys(mask, {row ? "row" : "column", "xor_hex"});
        size_t index =
            Integer(mask.at(row ? "row" : "column"), 0, (row ? sn : wn) - 1);
        auto bytes = Unhex(mask.at("xor_hex"), row ? wn : sn);
        for (size_t i = 0; i < bytes.size(); ++i) {
          block[row ? index * wn + i : i * wn + index] ^= bytes[i];
        }
      }
    }
  }
  auto state = [&](const Bytes& current) {
    Bytes delta(current.size()), information;
    json sparse = json::array(), bad_rows = json::array(),
         bad_columns = json::array();
    for (size_t r = 0; r < sn; ++r) {
      Bytes row(current.begin() + r * wn, current.begin() + (r + 1) * wn);
      if (Encode(Bytes(row.begin(), row.begin() + wk), wn) != row) {
        bad_rows.push_back(r);
      }
      for (size_t c = 0; c < wn; ++c) {
        const size_t p = r * wn + c;
        delta[p] = current[p] ^ transmitted[p];
        if (r < sk && c < wk) {
          information.push_back(delta[p]);
        }
        if (delta[p]) {
          sparse.push_back({{"row", r},
                            {"column", c},
                            {"xor_hex", Hex({delta[p]})},
                            {"actual_hex", Hex({current[p]})},
                            {"transmitted_hex", Hex({transmitted[p]})}});
        }
      }
    }
    for (size_t c = 0; c < wn; ++c) {
      Bytes col(sn);
      for (size_t r = 0; r < sn; ++r) {
        col[r] = current[r * wn + c];
      }
      if (Encode(Bytes(col.begin(), col.begin() + sk), sn) != col) {
        bad_columns.push_back(c);
      }
    }
    return json{{"full_residual", Weight(delta)},
                {"information_residual", Weight(information)},
                {"residual_errors", sparse},
                {"block_hex", Hex(current)},
                {"components",
                 {{"valid", bad_rows.empty() && bad_columns.empty()},
                  {"invalid_rows", bad_rows},
                  {"invalid_columns", bad_columns}}}};
  };
  auto initial = state(block);
  const auto result = code.Correct(block, options);
  auto final = state(block);
  bool counters_match = true;
  if (input.contains("expected_correction")) {
    const auto& expected_correction = input.at("expected_correction");
    Keys(expected_correction,
         {"stall_patterns_corrected", "strong_miscorrections_detected",
          "strong_miscorrections_corrected"});
    for (const auto& [key, actual] :
         std::initializer_list<std::pair<const char*, uint64_t>>{
             {"stall_patterns_corrected", result.stall_patterns_corrected},
             {"strong_miscorrections_detected",
              result.strong_miscorrections_detected},
             {"strong_miscorrections_corrected",
              result.strong_miscorrections_corrected}}) {
      if (expected_correction.contains(key)) {
        counters_match &=
            Integer(expected_correction.at(key), 0, 256) == actual;
      }
    }
  }
  const std::string outcome = block == transmitted ? "recovered"
                              : final["components"]["valid"].get<bool>()
                                  ? "undetected-failure"
                                  : "detected-failure";
  const char* termination =
      result.termination == ProductTermination::no_change ? "no-change"
      : result.termination == ProductTermination::pass_limit
          ? "pass-limit"
          : "invalid-argument";
  return {{"schema_version", 1},
          {"fixture", input},
          {"initial", initial},
          {"final", final},
          {"outcome", outcome},
          {"assertion_passed",
           counters_match && (expected.empty() || expected == outcome)},
          {"postprocessing", options.use_postprocessing ? "enabled" : "none"},
          {"correction",
           {{"termination", termination},
            {"all_zero_syndromes", result.all_zero_syndromes},
            {"directional_passes", result.directional_passes},
            {"strong_lines_visited", result.strong_lines_visited},
            {"weak_lines_visited", result.weak_lines_visited},
            {"changed_symbols", result.changed_symbols},
            {"changed_bits", result.changed_bits},
            {"strong_changed_symbols", result.strong_changed_symbols},
            {"strong_changed_bits", result.strong_changed_bits},
            {"weak_changed_symbols", result.weak_changed_symbols},
            {"weak_changed_bits", result.weak_changed_bits},
            {"stall_patterns_corrected", result.stall_patterns_corrected},
            {"strong_miscorrections_detected",
             result.strong_miscorrections_detected},
            {"strong_miscorrections_corrected",
             result.strong_miscorrections_corrected}}}};
}

json Generate(const json& args) {
  const size_t n = Number(args, "strong_n", 256, 4, 256);
  const size_t k = Number(args, "strong_k", 224, 2, 254);
  Require(k < n && std::has_single_bit(n) && std::has_single_bit(n - k) &&
              n - k >= 2 && n - k <= k,
          "unsupported strong dimensions");
  const size_t iterations = Number(args, "iterations", 2000, 1, 1000000);
  const size_t restart = Number(args, "restart_interval", 128, 1, 1000000);
  const size_t max_bits = Number(args, "max_perturb_bits", 8, 2, 8 * k);
  const size_t retain = Number(args, "retain", 4, 1, 32);
  const uint64_t seed =
      args.contains("seed") ? Integer(args.at("seed"), 0, UINT64_MAX) : 1;
  std::mt19937_64 rng(seed);
  // Linear binary contributions permit arbitrary multi-symbol message changes.
  std::vector<Bytes> contributions;
  for (size_t bit = 0; bit < 8 * k; ++bit) {
    Bytes message(k);
    message[bit / 8] = 1u << (bit % 8);
    contributions.push_back(Encode(message, n));
  }
  std::vector<Bytes> best;
  auto less = [](const Bytes& a, const Bytes& b) {
    return Bits(a) != Bits(b) ? Bits(a) < Bits(b) : a < b;
  };
  auto keep = [&](const Bytes& word) {
    if (Bits(word) == 0 ||
        std::find(best.begin(), best.end(), word) != best.end()) {
      return;
    }
    best.push_back(word);
    std::sort(best.begin(), best.end(), less);
    if (best.size() > retain) {
      best.pop_back();
    }
  };
  for (const auto& word : contributions) {
    keep(word);
  }
  Bytes current = best.front();
  for (size_t step = 0; step < iterations; ++step) {
    if (step % restart == 0) {
      current = contributions[rng() % contributions.size()];
    }
    // Alternate single-bit hill climbing with distinct multi-bit perturbations.
    const size_t count = step % 2 == 0 ? 2 + rng() % (max_bits - 1) : 1;
    std::set<size_t> bits;
    while (bits.size() < count) {
      bits.insert(rng() % (8 * k));
    }
    auto trial = current;
    for (size_t bit : bits) {
      for (size_t i = 0; i < n; ++i) {
        trial[i] ^= contributions[bit][i];
      }
    }
    keep(trial);
    if (Bits(trial) && Bits(trial) <= Bits(current)) {
      current = std::move(trial);
    }
  }
  json codewords = json::array();
  const size_t t = (n - k) / 2;
  LCHDecoder decoder(k, n - k);
  for (const auto& word : best) {
    Bytes message(word.begin(), word.begin() + k);
    Require(Encode(message, n) == word,
            "generated word failed independent reencoding");
    std::vector<size_t> positions;
    for (size_t i = 0; i < n; ++i) {
      if (word[i]) {
        positions.push_back(i);
      }
    }
    Require(positions.size() >= 2 * t + 1,
            "nonzero word violates distance bound");
    std::sort(positions.begin(), positions.end(), [&](size_t a, size_t b) {
      return std::popcount(word[a]) != std::popcount(word[b])
                 ? std::popcount(word[a]) < std::popcount(word[b])
                 : a < b;
    });
    Bytes received(n);
    for (size_t i = 0; i < positions.size() - t; ++i) {
      received[positions[i]] = word[positions[i]];
    }
    auto corrected = received;
    auto result = CorrectCodeword(decoder, corrected);
    Require(result.status == CorrectionStatus::ok && result.error_count == t &&
                corrected == word,
            "miscorrection witness did not decode to exact wrong word");
    codewords.push_back({{"input", Hex(message)}, {"codeword", Hex(word)}});
  }
  return {{"code params", {{"n", n}, {"k", k}}}, {"codewords", codewords}};
}

json ReadJSONStream(std::istream& stream, size_t limit = 4 * 1024 * 1024) {
  std::string text;
  char ch;
  while (stream.get(ch)) {
    Require(text.size() < limit, "JSON input exceeds size limit");
    Require(ch != '\0', "raw NUL byte in fixture");
    text += ch;
  }
  Require(!stream.bad(), "fixture read failed");
  std::vector<std::set<std::string>> keys;
  auto callback = [&](int depth, json::parse_event_t event, json& parsed) {
    Require(depth < 64, "JSON nesting too deep");
    if (event == json::parse_event_t::object_start) {
      keys.emplace_back();
    }
    if (event == json::parse_event_t::key) {
      Require(keys.back().insert(parsed.get<std::string>()).second,
              "duplicate JSON key");
    }
    if (event == json::parse_event_t::object_end) {
      keys.pop_back();
    }
    return true;
  };
  return json::parse(text, callback);
}

json ReadJSON(const std::string& path, size_t limit = 4 * 1024 * 1024) {
  std::ifstream file;
  if (path != "-") {
    file.open(path);
    Require(file.is_open(), "cannot open JSON input: " + path);
  }
  return ReadJSONStream(path == "-" ? std::cin : file, limit);
}

json Pattern(const std::string& name, const json& args, bool decode = true) {
  Require(name == "stall" || name == "miscorrection" || name == "mixed",
          "expected pattern stall|miscorrection|mixed");
  Keys(args, {"strong_n", "strong_k", "weak_n", "weak_k", "seed", "attempts",
              "columns", "codewords", "index"});
  const size_t sn = Number(args, "strong_n", 256, 1, 256);
  const size_t sk = Number(args, "strong_k", 224, 1, 255);
  const size_t wn = Number(args, "weak_n", 175, 4, 256);
  const size_t wk = Number(args, "weak_k", 173, 2, 254);
  Require(StrongWeakRSProductCode(sn, sk, wn, wk).Valid(),
          "unsupported product dimensions");
  Require(wn - wk == 2,
          "pattern construction unsupported for weak R=4; requires R=2");
  Require(sn - sk >= 4, "pattern construction requires strong R>=4");
  const size_t t = (sn - sk) / 2;
  const size_t attempts = Number(args, "attempts", 128, 1, 10000);
  const uint64_t seed =
      args.contains("seed") ? Integer(args.at("seed"), 0, UINT64_MAX) : 1;
  std::mt19937_64 rng(seed);
  std::vector<size_t> columns = name == "miscorrection"
                                    ? std::vector<size_t>{0}
                                    : std::vector<size_t>{0, 1};
  if (args.contains("columns")) {
    const auto text = args.at("columns").get<std::string>();
    columns.clear();
    size_t start = 0;
    do {
      const size_t end = text.find(',', start);
      const auto part = text.substr(start, end - start);
      size_t value;
      const auto parsed =
          std::from_chars(part.data(), part.data() + part.size(), value);
      Require(parsed.ec == std::errc{} &&
                  parsed.ptr == part.data() + part.size() && value < wn,
              "invalid column index");
      columns.push_back(value);
      if (end == std::string::npos) {
        break;
      }
      start = end + 1;
    } while (true);
  }
  Require(
      columns.size() == (name == "miscorrection" ? 1u : 2u) &&
          std::set<size_t>(columns.begin(), columns.end()).size() ==
              columns.size(),
      "columns must be distinct: one for miscorrection, two for stall/mixed");
  Require(!args.contains("index") || args.contains("codewords"),
          "--index requires --codewords");
  Require(name != "stall" || !args.contains("codewords"),
          "stall does not use --codewords");
  LCHDecoder strong(sk, sn - sk);
  Bytes word(sn), witness(sn);
  std::vector<size_t> support;
  json construction = {{"seed", seed},
                       {"attempt_limit", attempts},
                       {"columns", columns},
                       {"strong_radius", t},
                       {"verified", false}};
  if (name != "stall") {
    if (args.contains("codewords")) {
      const auto source = ReadJSON(args.at("codewords").get<std::string>());
      Keys(source, {"code params", "codewords"});
      const auto& params = source.at("code params");
      Keys(params, {"n", "k"});
      Require(Integer(params.at("n"), 1, 256) == sn &&
                  Integer(params.at("k"), 1, 255) == sk,
              "source codeword dimensions do not match strong dimensions");
      const auto& words = source.at("codewords");
      Require(words.is_array() && !words.empty(),
              "codewords must be a nonempty array");
      const size_t index = Number(args, "index", 0, 0, words.size() - 1);
      const auto& selected = words.at(index);
      Keys(selected, {"input", "codeword"});
      word = Unhex(selected.at("codeword"), sn);
      const auto message = Unhex(selected.at("input"), sk);
      Require(std::equal(message.begin(), message.end(), word.begin()) &&
                  Encode(message, sn) == word,
              "source word is not a valid systematic codeword");
      construction["source"] = {{"kind", "compact-codewords"},
                                {"index", index}};
    } else {
      Bytes message(sk);
      const size_t bit = rng() % (8 * sk);
      message[bit / 8] = 1u << (bit % 8);
      word = Encode(message, sn);
      construction["source"] = {{"kind", "seeded-single-information-bit"},
                                {"bit", bit}};
    }
    for (size_t r = 0; r < sn; ++r) {
      if (word[r]) {
        support.push_back(r);
      }
    }
    Require(support.size() >= 2 * t + 1,
            "source word must be nonzero and satisfy the distance bound");
    Require(args.contains("codewords") || support.size() == 2 * t + 1,
            "construction failed: single-information-bit word is not "
            "minimum-symbol");
    std::sort(support.begin(), support.end(), [&](size_t a, size_t b) {
      return std::popcount(word[a]) != std::popcount(word[b])
                 ? std::popcount(word[a]) < std::popcount(word[b])
                 : a < b;
    });
    for (size_t i = 0; i < support.size() - t; ++i) {
      witness[support[i]] = word[support[i]];
    }
    auto corrected = witness;
    const auto result = CorrectCodeword(strong, corrected);
    Require(result.status == CorrectionStatus::ok && result.error_count == t &&
                corrected == word,
            "construction failed: witness does not miscorrect to the exact "
            "source word");
    construction["wrong_codeword_hex"] = Hex(word);
    construction["wrong_codeword_weight"] = Weight(word);
    construction["witness_weight"] = Weight(witness);
    construction["witness_exact_miscorrection"] = true;
    construction["witness_corrections"] = result.error_count;
    support.resize(support.size() - t);
  } else {
    for (size_t r = 0; r < sn; ++r) {
      support.push_back(r);
    }
    for (size_t i = support.size(); i > 1; --i) {
      std::swap(support[i - 1], support[rng() % i]);
    }
    support.resize(t + 1);
  }
  // Selection uses component outcomes, independent weak proposals, and baseline
  // bytes only. Never run postprocessing until one construction is selected.
  for (size_t attempt = 1; attempt <= attempts; ++attempt) {
    Bytes first = witness, second(sn);
    if (name == "stall") {
      for (size_t r : support) {
        first[r] = second[r] = 1 + rng() % 255;
      }
    } else if (name == "mixed") {
      second = witness;
      const size_t r = support[rng() % support.size()];
      // Equal overlapping errors cancel the symbol-sum syndrome (no R=2 BDD
      // proposal). Break one equality to move the second column outside c's
      // radius; only this row can add a spurious clean-column proposal.
      second[r] = static_cast<Element>(1 + (word[r] + rng() % 254) % 255);
    }
    json strong_checks = json::array();
    bool failed = false;
    for (size_t i = 0; i < columns.size(); ++i) {
      const auto& received = i == 0 ? first : second;
      auto corrected = received;
      const auto result = CorrectCodeword(strong, corrected);
      const bool wrong = name != "stall" && i == 0;
      const bool verified =
          wrong ? result.status == CorrectionStatus::ok &&
                      result.error_count == t && corrected == word
                : result.status == CorrectionStatus::uncorrectable &&
                      result.error_count == 0 && corrected == received;
      failed |= !verified;
      strong_checks.push_back(
          {{"column", columns[i]},
           {"status", wrong ? "miscorrected" : "uncorrectable"},
           {"corrections", result.error_count},
           {"verified", verified},
           {"unchanged", corrected == received},
           {"received_weight", Weight(received)}});
    }
    if (failed) {
      continue;
    }
    Bytes after_strong(sn * wn);
    for (size_t r = 0; r < sn; ++r) {
      after_strong[r * wn + columns[0]] = name == "stall" ? first[r] : word[r];
      if (columns.size() == 2) {
        after_strong[r * wn + columns[1]] = second[r];
      }
    }
    std::vector<size_t> proposals(wn);
    bool inconsistent = false;
    const size_t mother_n = std::bit_ceil(wn), mother_k = mother_n - 2;
    LCHDecoder weak(mother_k, 2);
    Bytes unit(wn);
    if (name == "mixed") {
      unit[columns[1]] = 1;
    }
    const auto encoded_unit =
        Encode(Bytes(unit.begin(), unit.begin() + wk), wn);
    for (size_t r = 0; r < sn; ++r) {
      Bytes row(after_strong.begin() + r * wn,
                after_strong.begin() + (r + 1) * wn);
      Bytes padded(mother_n);
      std::copy_n(row.begin(), wk, padded.begin());
      std::copy_n(row.begin() + wk, 2, padded.begin() + mother_k);
      auto candidate = padded;
      const auto result = CorrectCodeword(weak, candidate);
      if (result.status == CorrectionStatus::ok && result.error_count == 1 &&
          std::all_of(candidate.begin() + wk, candidate.begin() + mother_k,
                      [](auto b) { return b == 0; })) {
        for (size_t c = 0; c < wn; ++c) {
          const size_t p = c < wk ? c : mother_k + c - wk;
          if (candidate[p] != padded[p]) {
            ++proposals[c];
          }
        }
      }
      if (name == "mixed") {
        // A one-column erasure can explain the row only if its two parity
        // residuals are a scalar multiple of the erased unit's residuals.
        const auto encoded_row =
            Encode(Bytes(row.begin(), row.begin() + wk), wn);
        const Element a = row[wk] ^ encoded_row[wk],
                      b = row[wk + 1] ^ encoded_row[wk + 1];
        const Element u = unit[wk] ^ encoded_unit[wk],
                      v = unit[wk + 1] ^ encoded_unit[wk + 1];
        inconsistent |=
            gf2p8::MultiplyCantor(a, v) != gf2p8::MultiplyCantor(b, u);
      }
    }
    if (name == "mixed" &&
        (!inconsistent || *std::max_element(proposals.begin(),
                                            proposals.end()) >= sn - sk + 1)) {
      continue;
    }
    json fixture = {
        {"strong_n", sn},
        {"strong_k", sk},
        {"weak_n", wn},
        {"weak_k", wk},
        {"options", {{"use_postprocessing", false}}},
        {"expected", "detected-failure"},
        {"expected_correction",
         {{"stall_patterns_corrected", 0},
          {"strong_miscorrections_detected", 0},
          {"strong_miscorrections_corrected", 0}}},
        {"column_masks",
         json::array({{{"column", columns[0]}, {"xor_hex", Hex(first)}}})}};
    if (columns.size() == 2) {
      fixture["column_masks"].push_back(
          {{"column", columns[1]}, {"xor_hex", Hex(second)}});
    }
    auto baseline = Run(fixture);
    auto invalid_columns = columns;
    if (name != "stall") {
      invalid_columns.erase(invalid_columns.begin());
    }
    std::sort(invalid_columns.begin(), invalid_columns.end());
    if (baseline["outcome"] != "detected-failure" ||
        baseline["final"]["block_hex"] != Hex(after_strong) ||
        baseline["final"]["components"]["invalid_columns"] != invalid_columns) {
      continue;
    }
    construction["verified"] = true;
    construction["attempts_used"] = attempt;
    construction["strong_outcomes"] = strong_checks;
    construction["weak_proposal_support"] = json::object();
    for (size_t c = 0; c < wn; ++c) {
      if (proposals[c] != 0) {
        construction["weak_proposal_support"][std::to_string(c)] = proposals[c];
      }
    }
    construction["consensus_threshold"] = sn - sk + 1;
    construction["baseline_matches_initial_strong_output"] = true;
    if (name == "mixed") {
      construction["known_column_erasure_inconsistent"] = inconsistent;
    }
    if (!decode) {
      fixture.erase("expected");
      fixture.erase("expected_correction");
      return {{"construction", construction}, {"fixture", fixture}};
    }
    fixture["options"]["use_postprocessing"] = true;
    fixture["expected"] = name == "mixed" ? "detected-failure" : "recovered";
    fixture["expected_correction"] = {
        {"stall_patterns_corrected", name == "stall" ? 1 : 0},
        {"strong_miscorrections_detected", name == "stall" ? 0 : 1},
        {"strong_miscorrections_corrected", name == "miscorrection" ? 1 : 0}};
    auto postprocessed = Run(fixture);
    const bool unchanged =
        name != "mixed" ||
        postprocessed["final"]["block_hex"] == baseline["final"]["block_hex"];
    const bool matched = baseline["assertion_passed"].get<bool>() &&
                         postprocessed["assertion_passed"].get<bool>() &&
                         unchanged;
    return {{"schema_version", 1},
            {"scenario", name},
            {"construction", construction},
            {"baseline", baseline},
            {"postprocessed", postprocessed},
            {"expected_matched", matched},
            {"assertion_passed", matched}};
  }
  throw std::runtime_error("construction failed: exhausted " +
                           std::to_string(attempts) +
                           " attempts without verified component/baseline "
                           "conditions; postprocessing was not run");
}

// Rejection sampling avoids modulo bias and library-specific distributions.
uint64_t Uniform(std::mt19937_64& rng, uint64_t bound) {
  const uint64_t threshold = -bound % bound;
  uint64_t value;
  do {
    value = rng();
  } while (value < threshold);
  return value % bound;
}

void Sample(const std::string& scenario, const json& args) {
  Require(
      scenario == "none" || scenario == "stall" ||
          scenario == "miscorrection" || scenario == "mixed" ||
          scenario == "miscorrection-stall" || scenario == "two-miscorrections",
      "expected scenario "
      "none|stall|miscorrection|mixed|miscorrection-stall|two-miscorrections");
  Keys(args, {"strong_n", "strong_k", "weak_n", "weak_k", "seed", "attempts",
              "codewords", "index", "random_bits", "samples", "details"});
  const size_t sn = Number(args, "strong_n", 256, 1, 256);
  const size_t sk = Number(args, "strong_k", 224, 1, 255);
  const size_t wn = Number(args, "weak_n", 175, 4, 256);
  const size_t wk = Number(args, "weak_k", 173, 2, 254);
  Require(StrongWeakRSProductCode(sn, sk, wn, wk).Valid(),
          "unsupported product dimensions");
  Require(scenario == "none" || (wn - wk == 2 && sn - sk >= 4),
          "planted scenarios require weak R=2 and strong R>=4");
  Require(!args.contains("index") || args.contains("codewords"),
          "--index requires --codewords");
  Require((scenario != "none" && scenario != "stall") ||
              (!args.contains("codewords") && !args.contains("index")),
          "codeword sources require a miscorrection scenario");
  const size_t count = Number(args, "samples", 1, 1, 1000000);
  const bool details = Number(args, "details", 0, 0, 1);
  Require(args.contains("random_bits"), "sample requires --random-bits K");
  const size_t protected_columns =
      scenario == "two-miscorrections" ? 2
      : (scenario == "miscorrection" || scenario == "mixed" ||
         scenario == "miscorrection-stall")
          ? 1
          : 0;
  const size_t bit_count = 8 * sn * (wn - protected_columns);
  Require(Integer(args.at("random_bits"), 0, 8 * sn * wn) <= bit_count,
          "random bits exceed eligible space outside miscorrection columns");
  const size_t flips = Number(args, "random_bits", 0, 0, bit_count);
  Number(args, "attempts", 128, 1, 10000);
  const uint64_t seed =
      args.contains("seed") ? Integer(args.at("seed"), 0, UINT64_MAX) : 1;
  std::mt19937_64 seeds(seed);
  json document;
  document["samples"] = json::array();
  for (size_t trial = 0; trial < count; ++trial) {
    const uint64_t planting_seed = seeds(), noise_seed = seeds();
    std::mt19937_64 planting(planting_seed), noise(noise_seed);
    json fixture = {{"strong_n", sn},
                    {"strong_k", sk},
                    {"weak_n", wn},
                    {"weak_k", wk},
                    {"options", {{"use_postprocessing", false}}}};
    json construction = json::array();
    std::vector<size_t> columns(wn);
    for (size_t i = 0; i < wn; ++i) {
      columns[i] = i;
    }
    for (size_t i = 0; i < 3; ++i) {
      std::swap(columns[i], columns[i + Uniform(planting, wn - i)]);
    }
    const auto plant = [&](const std::string& kind, size_t offset) {
      json options = args;
      options.erase("samples");
      options.erase("random_bits");
      options.erase("details");
      options["seed"] = planting();
      options["columns"] = std::to_string(columns[offset]);
      if (kind != "miscorrection") {
        options["columns"] = options["columns"].get<std::string>() + "," +
                             std::to_string(columns[offset + 1]);
      }
      if (kind == "stall") {
        options.erase("codewords");
        options.erase("index");
      }
      auto result = Pattern(kind, options, false);
      construction.push_back({{"scenario", kind},
                              {"seed", options["seed"]},
                              {"checks", result["construction"]}});
      for (const auto& mask : result["fixture"]["column_masks"]) {
        fixture["column_masks"].push_back(mask);
      }
    };
    if (scenario == "miscorrection-stall") {
      plant("miscorrection", 0);
      plant("stall", 1);
    } else if (scenario == "two-miscorrections") {
      plant("miscorrection", 0);
      plant("miscorrection", 1);
    } else if (scenario != "none") {
      plant(scenario, 0);
    }
    std::vector<size_t> eligible_columns;
    for (size_t c = 0; c < wn; ++c) {
      if (std::find(columns.begin(), columns.begin() + protected_columns, c) ==
          columns.begin() + protected_columns) {
        eligible_columns.push_back(c);
      }
    }
    // Uniform K-subset of the bits outside planted miscorrection columns.
    std::set<size_t> positions;
    for (size_t j = bit_count - flips; j < bit_count; ++j) {
      const size_t bit = Uniform(noise, j + 1);
      if (!positions.insert(bit).second) {
        positions.insert(j);
      }
    }
    Bytes random_mask(sn * wn), planted(sn * wn);
    for (size_t bit : positions) {
      const size_t symbol = bit / 8;
      const size_t row = symbol / eligible_columns.size();
      const size_t col = eligible_columns[symbol % eligible_columns.size()];
      random_mask[row * wn + col] ^= 1u << (bit % 8);
    }
    if (fixture.contains("column_masks")) {
      for (const auto& mask : fixture["column_masks"]) {
        const size_t c = mask["column"];
        auto values = Unhex(mask["xor_hex"], sn);
        for (size_t r = 0; r < sn; ++r) {
          planted[r * wn + c] ^= values[r];
        }
      }
    }
    for (size_t r = 0; r < sn; ++r) {
      Bytes row(random_mask.begin() + r * wn,
                random_mask.begin() + (r + 1) * wn);
      if (Bits(row) != 0) {
        fixture["row_masks"].push_back({{"row", r}, {"xor_hex", Hex(row)}});
      }
    }
    size_t cancelled = 0;
    for (size_t i = 0; i < planted.size(); ++i) {
      cancelled += std::popcount(unsigned(planted[i] & random_mask[i]));
    }
    json record = {{"type", "sample"},
                   {"schema_version", 1},
                   {"scenario", scenario},
                   {"sample_index", trial},
                   {"seed", seed},
                   {"planting_seed", planting_seed},
                   {"noise_seed", noise_seed},
                   {"random_bits", flips},
                   {"planted_weight", Weight(planted)},
                   {"cancelled_planted_bits", cancelled},
                   {"construction", construction},
                   {"fixture", fixture}};
    // Emit the complete combined error, without decoding the noisy sample.
    Bytes errors = planted;
    for (size_t i = 0; i < errors.size(); ++i) {
      errors[i] ^= random_mask[i];
    }
    record["received_weight"] = Weight(errors);
    record["error_hex"] = Hex(errors);
    if (trial == 0) {
      json generation = {
          {"scenario", scenario},
          {"random_bits", flips},
          {"seed", seed},
          {"noise",
           "uniform distinct bits, XOR outside planted miscorrection columns"},
          {"excluded_miscorrection_columns", protected_columns},
          {"eligible_noise_bits", bit_count},
          {"columns", "uniform distinct columns, including parity"}};
      if (scenario != "none") {
        generation["attempts"] = Number(args, "attempts", 128, 1, 10000);
      }
      if (args.contains("codewords")) {
        generation["codewords"] = args["codewords"];
        generation["index"] = args.value("index", uint64_t{0});
      }
      document.update(
          json({{"schema_version", 4},
                {"code_params",
                 {{"strong", {{"n", sn}, {"k", sk}}},
                  {"weak", {{"n", wn}, {"k", wk}}}}},
                {"error_generation", generation},
                {"error_format",
                 "row-major Cantor bytes, XOR against all-zero product word"},
                {"sample_count", count}}));
    }
    if (!details) {
      record = {{"type", "sample"},
                {"sample_index", trial},
                {"error_hex", Hex(errors)}};
    }
    document["samples"].push_back(std::move(record));
  }
  std::cout << document.dump(2) << '\n';
  Require(bool(std::cout), "sample output failed");
}

void TestSamples(const std::string& path) {
  const auto header = ReadJSON(path, 512 * 1024 * 1024);
  Keys(header, {"schema_version", "code_params", "error_generation",
                "error_format", "sample_count", "samples"});
  Require(header.at("schema_version") == 4,
          "test requires schema-4 generated samples");
  Require(header.at("error_format") ==
              "row-major Cantor bytes, XOR against all-zero product word",
          "unsupported error format");
  const auto& params = header.at("code_params");
  Keys(params, {"strong", "weak"});
  Keys(params.at("strong"), {"n", "k"});
  Keys(params.at("weak"), {"n", "k"});
  const size_t sn = Integer(params.at("strong").at("n"), 1, 256);
  const size_t sk = Integer(params.at("strong").at("k"), 1, 255);
  const size_t wn = Integer(params.at("weak").at("n"), 1, 256);
  const size_t wk = Integer(params.at("weak").at("k"), 1, 255);
  Require(StrongWeakRSProductCode(sn, sk, wn, wk).Valid(),
          "unsupported product dimensions");
  const size_t count = Integer(header.at("sample_count"), 1, 1000000);
  Require(
      header.at("samples").is_array() && header.at("samples").size() == count,
      "sample count mismatch");
  json results = json::array();
  json totals = {
      {"recovered", 0}, {"undetected-failure", 0}, {"detected-failure", 0}};
  for (size_t i = 0; i < count; ++i) {
    const auto& sample = header.at("samples").at(i);
    Keys(sample, {"type", "sample_index", "error_hex", "schema_version",
                  "scenario", "seed", "planting_seed", "noise_seed",
                  "random_bits", "planted_weight", "cancelled_planted_bits",
                  "construction", "fixture", "received_weight"});
    Require(sample.at("type") == "sample" &&
                Integer(sample.at("sample_index"), 0, count - 1) == i,
            "sample index/type mismatch");
    const auto errors = Unhex(sample.at("error_hex"), sn * wn);
    json fixture = {
        {"strong_n", sn}, {"strong_k", sk}, {"weak_n", wn}, {"weak_k", wk}};
    for (size_t r = 0; r < sn; ++r) {
      fixture["row_masks"].push_back(
          {{"row", r},
           {"xor_hex", Hex(Bytes(errors.begin() + r * wn,
                                 errors.begin() + (r + 1) * wn))}});
    }
    fixture["options"]["use_postprocessing"] = true;
    const auto decoded = Run(fixture);
    const std::string outcome = decoded["outcome"];
    totals[outcome] = totals[outcome].get<uint64_t>() + 1;
    const auto& correction = decoded["correction"];
    json patterns = json::array();
    if (correction["stall_patterns_corrected"].get<uint64_t>() != 0) {
      patterns.push_back("stall-repaired");
    }
    if (correction["strong_miscorrections_detected"].get<uint64_t>() != 0) {
      patterns.push_back("strong-miscorrection-detected");
    }
    if (correction["strong_miscorrections_corrected"].get<uint64_t>() != 0) {
      patterns.push_back("strong-miscorrection-repaired");
    }
    json result = {
        {"sample_index", i}, {"outcome", outcome}, {"patterns", patterns}};
    results.push_back(std::move(result));
  }
  std::cout << json({{"schema_version", 5},
                     {"sample_count", count},
                     {"results", results},
                     {"outcomes", totals}})
                   .dump(2)
            << '\n';
  Require(bool(std::cout), "test summary output failed");
}

constexpr const char* kHelp =
    R"(rs-product-test: deterministic product-decoder fixtures

Usage:
  rs-product-test generate [options]
  rs-product-test pattern stall|miscorrection|mixed [options]
  rs-product-test sample SCENARIO --random-bits K [options]
  rs-product-test test FILE.json    Decode generated samples (or - for stdin)
  rs-product-test run FILE.json
  rs-product-test run -             Read JSON from stdin
  rs-product-test --help

Sampling:
  SCENARIO: none|stall|miscorrection|mixed|miscorrection-stall
            |two-miscorrections
  --samples N --seed UINT64      Default: 1, 1; samples 1..1000000
  --random-bits K                Exactly K distinct uniformly sampled bits
  --details 0|1                 Full construction/replay records (default: 0)
  Uses pattern dimensions, attempts and optional codewords/index below.
  Planted columns are uniformly selected, distinct, including parity columns.
  miscorrection-stall plants a miscorrection plus a separate two-column stall.
  two-miscorrections plants two independently verified strong miscorrections.
  Both use the selected source word if --codewords/--index is supplied;
  otherwise each uses its own seeded word.
  Noise excludes every planted strong-miscorrection column, preserving its
  initial BDD miscorrection. Elsewhere it may cancel planted errors. K must fit
  the eligible bits. Subsequent product decoding may still repair the sample.
  none/stall allow noise across the full block. No outcome filtering.
  Indented JSON: code params, error-generation scheme, sample_count and samples
  array. Each sample contains its index and full combined error_hex.
  error_hex: 2*Nstrong*Nweak hex digits, row-major Cantor-coordinate bytes.
  These XOR masks target an all-zero transmitted product word.
  Sample does not decode noisy samples. Planting still verifies clean scenarios.
  --details 1 adds construction diagnostics and equivalent replay fixtures.
  Replay detailed fixture with run; enable options.use_postprocessing for ON.
  Construction failure exits 2 without a JSON document.
  These are conditional experiments, not weighted channel error-rate estimates.

Testing:
  test reads schema-4 sample JSON (up to 512 MiB) and decodes error_hex with
  postprocessing enabled. Schema-5 results contain one final outcome and patterns:
  stall-repaired, strong-miscorrection-detected, strong-miscorrection-repaired.
  Markers reflect decoder events, not planted labels; empty is not proof that
  no rare pattern occurred. A stall repair can still produce undetected-failure.
  Output is indented JSON with per-sample results and aggregate outcome totals.
  It does not regenerate errors or read source codeword files. Malformed input
  exits 2 without JSON output. Split larger experiments into smaller files.

Patterns:
  --strong-n N --strong-k K      Default: 256, 224
  --weak-n N --weak-k K          Default: 175, 173
  --seed UINT64 --attempts N     Default: 1, 128; attempts 1..10000
  --columns I[,J]                Zero-based, including parity; default: 0[,1]
  --codewords FILE --index N     Compact generator source; default index: 0
  Stall: two verified failed columns, equal errors on the same >t rows.
  Miscorrection: retain w-t lowest-popcount entries of a valid wrong word.
  Mixed: that witness plus an overlapping, verified failed second column.
  Source is for miscorrection/mixed only; selected word must be nonzero,
  systematic and valid with matching strong dimensions. Without a source,
  encode one seeded information bit: minimum symbols, NOT minimum bits.
  Requires supported strong dimensions with R>=4 and weak R=2 (not R=4).
  Select using strong outcomes, weak proposals/equations and baseline bytes;
  then run postprocessing ONCE, never resampling on an expected mismatch.
  Paired reports contain construction checks and replayable branch fixtures.

Generator:
  --strong-n N --strong-k K      Code dimensions (default: 256, 224)
  --seed UINT64                  Default: 1
  --iterations N                 Search trials (default: 2000; 1..1000000)
  --restart-interval N           Restart interval (default: 128; 1..1000000)
  --max-perturb-bits N           Message-bit perturbation cap (default: 8; 2..8*K)
  --retain N                     Words to retain (default: 4; 1..32)
  N and N-K must be powers of two; N<=256, 2<=N-K<=K (unshortened).
  Searches low-bitweight nonzero words; NOT a global minimum search.
  Each word is reencoded and verified to admit a miscorrection witness.

Fixture:
  Example (shown dimensions and decoder options are defaults):
  {"strong_n":256,"strong_k":224,"weak_n":256,"weak_k":254,
   "options":{"max_directional_passes":16,"use_anchors":true,
              "use_binary_image":true,"use_postprocessing":false},
   "errors":[{"row":0,"column":0,"xor_hex":"03"}],"expected":"recovered"}
  Optional row_masks:    [{"row":1,"xor_hex":"<weak_n-byte hex>"}]
  Optional column_masks: [{"column":0,"xor_hex":"<strong_n-byte hex>"}]
  Optional transmitted_hex: valid row-major systematic product word
  (strong_n*weak_n bytes); default: all zero. Missing mask arrays are empty.
  Coordinates are zero-based, including parity; all masks compose by XOR.
  Hex: exactly two digits per byte, no whitespace or prefix.
  Dimensions must be supported by the product API; pass cap: 2..1024.
  Postprocessing is off by default; when enabled, runs once after stall/cap.
  Optional expected: recovered (exact match), undetected-failure (valid but wrong),
  or detected-failure (invalid final component, not a diagnosis of a failed pass).
  Component validity alone cannot prove equality to the transmitted word.
  Optional expected_correction: object asserting stall_patterns_corrected,
  strong_miscorrections_detected and/or strong_miscorrections_corrected.
  Limits: 4 MiB input, 65536 masks per array; unknown/duplicate keys rejected.

Output:
  Indented JSON on stdout; no files created or overwritten.
  Generate emits only:
  {"code params":{"n":N,"k":K},
   "codewords":[{"input":"<K-byte hex>","codeword":"<N-byte hex>"}]}
  No metadata, weights, witnesses, or fixtures are embedded in generation output.
  Run reports outcome, residual errors and correction work.
  Exit: 0 success; 1 expected-outcome mismatch (report emitted);
  2 invalid input/construction failure (stderr diagnostic, no JSON report).
)";
}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 2 && std::string(argv[1]) == "--help") {
      std::cout << kHelp;
      return 0;
    }
    Require(argc >= 2, "missing subcommand; use --help");
    if (std::string(argv[1]) == "test") {
      Require(argc == 3, "expected test FILE.jsonl");
      TestSamples(argv[2]);
      return 0;
    }
    json report;
    const bool pattern = std::string(argv[1]) == "pattern";
    const bool sample = std::string(argv[1]) == "sample";
    if (std::string(argv[1]) == "generate" || pattern || sample) {
      Require(!(pattern || sample) || argc >= 3, "missing scenario name");
      json args = json::object();
      for (int i = (pattern || sample) ? 3 : 2; i < argc; i += 2) {
        Require(i + 1 < argc, "missing option value");
        std::string key = argv[i];
        Require(key.starts_with("--"), "expected --option");
        key.erase(0, 2);
        std::replace(key.begin(), key.end(), '-', '_');
        Require(!args.contains(key), "duplicate option");
        std::string text = argv[i + 1];
        if ((pattern || sample) && (key == "columns" || key == "codewords")) {
          args[key] = text;
          continue;
        }
        uint64_t number;
        auto parsed =
            std::from_chars(text.data(), text.data() + text.size(), number);
        Require(
            parsed.ec == std::errc{} && parsed.ptr == text.data() + text.size(),
            "invalid unsigned option value");
        args[key] = number;
      }
      if (sample) {
        Sample(argv[2], args);
        return 0;
      } else if (pattern) {
        report = Pattern(argv[2], args);
      } else {
        Keys(args, {"strong_n", "strong_k", "seed", "iterations",
                    "restart_interval", "max_perturb_bits", "retain"});
        report = Generate(args);
      }
    } else {
      Require(std::string(argv[1]) == "run" && argc == 3,
              "expected run FILE.json; use --help");
      report = Run(ReadJSON(argv[2]));
    }
    if (pattern) {
      nlohmann::ordered_json output;
      for (const auto* key :
           {"schema_version", "scenario", "construction", "expected_matched",
            "assertion_passed", "baseline", "postprocessed"}) {
        output[key] = report.at(key);
      }
      std::cout << output.dump(2) << '\n';
    } else {
      std::cout << report.dump(2) << '\n';
    }
    return report.value("assertion_passed", true) ? 0 : 1;
  } catch (const std::exception& error) {
    std::cerr << "rs-product-test: " << error.what() << '\n';
    return 2;
  }
}
