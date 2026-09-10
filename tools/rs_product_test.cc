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
  Keys(input,
       {"strong_n", "strong_k", "weak_n", "weak_k", "options",
        "transmitted_hex", "errors", "row_masks", "column_masks", "expected"});
  const size_t sn = Number(input, "strong_n", 256, 1, 256);
  const size_t sk = Number(input, "strong_k", 224, 1, 255);
  const size_t wn = Number(input, "weak_n", 256, 1, 256);
  const size_t wk = Number(input, "weak_k", 254, 1, 255);
  StrongWeakRSProductCode code(sn, sk, wn, wk);
  Require(code.Valid(), "unsupported product dimensions");
  ProductDecodeOptions options;
  if (input.contains("options")) {
    const auto& o = input.at("options");
    Keys(o, {"max_directional_passes", "use_anchors", "use_binary_image"});
    options.max_directional_passes =
        Number(o, "max_directional_passes", 16, 2, 1024);
    for (const auto* key : {"use_anchors", "use_binary_image"}) {
      if (o.contains(key)) {
        Require(o.at(key).is_boolean(), "option must be boolean");
      }
    }
    options.use_anchors = o.value("use_anchors", true);
    options.use_binary_image = o.value("use_binary_image", true);
  }
  std::string expected = input.value("expected", std::string());
  Require(!input.contains("expected") || expected == "recovered" ||
              expected == "valid-wrong" || expected == "detected-failure",
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
  const std::string outcome = block == transmitted ? "recovered"
                              : final["components"]["valid"].get<bool>()
                                  ? "valid-wrong"
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
          {"assertion_passed", expected.empty() || expected == outcome},
          {"postprocessing", "none"},
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
            {"weak_changed_bits", result.weak_changed_bits}}}};
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

constexpr const char* kHelp =
    R"(rs-product-test: testing fixtures only; no production postprocessor.
Usage:
  rs-product-test generate [--strong-n 256 --strong-k 224 --seed 1]
      [--iterations 2000 --restart-interval 128 --max-perturb-bits 8 --retain 4]
  rs-product-test run FILE.json
  rs-product-test run -                 (read JSON from stdin)
All results are indented JSON on stdout. No files are created or overwritten.
Exit codes: 0 success, 1 expected-outcome mismatch (report emitted), 2 invalid
input/internal verification failure (diagnostic on stderr, no JSON report).

Fixture schema (unknown fields rejected; coordinates zero-based):
 {"strong_n":256,"strong_k":224,"weak_n":256,"weak_k":254,
  "options":{"max_directional_passes":16,"use_anchors":true,"use_binary_image":true},
  "errors":[{"row":0,"column":0,"xor_hex":"03"}],
  "row_masks":[{"row":1,"xor_hex":"<weak_n bytes>"}],
  "column_masks":[{"column":0,"xor_hex":"<strong_n bytes>"}],
  "expected":"recovered"}
Omitted transmitted_hex means zero. Otherwise supply a full row-major valid
systematic product word (strong_n*weak_n bytes). Hex has exactly two digits per
byte, no whitespace or prefix. Masks compose by XOR, including duplicates and
parity coordinates. Expected is optional: recovered, valid-wrong, detected-failure.
Recovered means exact full equality; valid-wrong means independently valid but
different; detected-failure means an invalid final component, not a proof of
which decoder pass failed. Reports include full/information residual weights,
sparse coordinates, independent scalar component reencoding, and correction work.
Dimensions must be supported by the public product API; pass cap is 2..1024.
Input is limited to 4 MiB, arrays to 65536 masks; malformed/duplicate keys rejected.

Generation supports unshortened power-of-two strong N<=256, power-of-two
2<=N-K<=K. All K*8 one-bit seeds are encoded. Seeded mt19937_64 modulo sampling
alternates one-bit trials and 2..max-perturb-bits distinct message-bit trials,
XORing encoded contributions across arbitrary information symbols. Restart at a
random one-bit seed every restart-interval trials; accept nonzero non-increasing
bitweight moves, retain the best distinct nonzero words from ALL evaluated trials.
Ties use lexicographic bytes. This bounded heuristic is NOT a global minimum
search. Iterations: 1..1000000; retain: 1..32; perturb bits: 2..8*K.
Each retained word is independently reencoded. Its genuine witness retains the
w-t lowest-popcount nonzero symbols (index breaks ties), removing t highest.
Thus its injected symbolweight is w-t>t, distance to the wrong word is t, and
CorrectCodeword is verified to return that exact nonzero wrong word.
Generation output contains only:
 {"code params":{"n":N,"k":K},
  "codewords":[{"input":"<K-byte hex>","codeword":"<N-byte hex>"},...]}
No metadata, weights, witnesses, or fixtures are embedded in generation output.
)";
}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 2 && std::string(argv[1]) == "--help") {
      std::cout << kHelp;
      return 0;
    }
    Require(argc >= 2, "missing subcommand; use --help");
    json report;
    if (std::string(argv[1]) == "generate") {
      json args = json::object();
      for (int i = 2; i < argc; i += 2) {
        Require(i + 1 < argc, "missing option value");
        std::string key = argv[i];
        Require(key.starts_with("--"), "expected --option");
        key.erase(0, 2);
        std::replace(key.begin(), key.end(), '-', '_');
        Require(!args.contains(key), "duplicate option");
        std::string text = argv[i + 1];
        uint64_t number;
        auto parsed =
            std::from_chars(text.data(), text.data() + text.size(), number);
        Require(
            parsed.ec == std::errc{} && parsed.ptr == text.data() + text.size(),
            "invalid unsigned option value");
        args[key] = number;
      }
      Keys(args, {"strong_n", "strong_k", "seed", "iterations",
                  "restart_interval", "max_perturb_bits", "retain"});
      report = Generate(args);
    } else {
      Require(std::string(argv[1]) == "run" && argc == 3,
              "expected run FILE.json; use --help");
      std::ifstream file;
      if (std::string(argv[2]) != "-") {
        file.open(argv[2]);
        Require(file.is_open(), "cannot open fixture");
      }
      std::istream& stream = std::string(argv[2]) == "-" ? std::cin : file;
      std::string text;
      char ch;
      while (stream.get(ch)) {
        Require(text.size() < 4 * 1024 * 1024, "fixture exceeds 4 MiB");
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
      report = Run(json::parse(text, callback));
    }
    std::cout << report.dump(2) << '\n';
    return report.value("assertion_passed", true) ? 0 : 1;
  } catch (const std::exception& error) {
    std::cerr << "rs-product-test: " << error.what() << '\n';
    return 2;
  }
}
