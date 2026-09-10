#include <algorithm>
#include <bit>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "benchmark/benchmark.h"
#include "reed_solomon/product_code_internal.h"
#include "reed_solomon/strong_weak_rs_product_code.h"

namespace {

void BenchmarkProductCorrectionBSC(benchmark::State& state,
                                   int batch_passes = -1,
                                   size_t weak_n = 256,
                                   int optimizations = -1) {
  using gf2p8::Element;
  using gf2p8::rs::ProductCorrectionResult;
  using gf2p8::rs::ProductTermination;
  constexpr size_t kN = 256;
  constexpr size_t kStrongK = 224;
  const size_t kWeakK = weak_n - 2;
  const size_t kInformationBytes = kStrongK * kWeakK;
  const size_t kBlockBytes = kN * weak_n;
  constexpr size_t kCorpusCount = 64;
  constexpr size_t kPassLimit = 16;
  constexpr uint32_t kSeed = 0x5b5c0224;
  gf2p8::rs::StrongWeakRSProductCode code(kN, kStrongK, weak_n, kWeakK);
  const auto correct = [&](std::vector<Element>& block) {
    if (optimizations >= 0) {
      return gf2p8::rs::detail::ProductCorrectionAccess::Experiment(
          code, block, gf2p8::rs::ProductDecodeOptions{kPassLimit},
          optimizations);
    }
    return batch_passes < 0
               ? code.Correct(block, kPassLimit)
               : gf2p8::rs::detail::ProductCorrectionAccess::Correct(
                     code, block, kPassLimit, batch_passes);
  };
  if (!code.Valid() || code.BlockSize() != kBlockBytes) {
    state.SkipWithError("invalid product-code dimensions");
    return;
  }

  std::mt19937 messages(kSeed);
  std::mt19937 channel(kSeed ^ 0x9e3779b9U);
  std::vector<std::vector<Element>> original(kCorpusCount,
                                             std::vector<Element>(kBlockBytes));
  auto corrupted = original;
  auto work = original;
  std::vector<ProductCorrectionResult> results(kCorpusCount);
  uint64_t channel_bits = 0, channel_symbols = 0;
  for (size_t sample = 0; sample < kCorpusCount; ++sample) {
    auto& block = original[sample];
    for (size_t row = 0; row < kStrongK; ++row) {
      for (size_t col = 0; col < kWeakK; ++col) {
        block[row * weak_n + col] = static_cast<Element>(messages());
      }
    }
    if (code.Encode(block) != gf2p8::lch::Status::ok) {
      state.SkipWithError("product corpus encoding failed");
      return;
    }
    work[sample] = block;
    const auto clean = code.Correct(work[sample], kPassLimit);
    if (!clean.all_zero_syndromes || clean.changed_symbols != 0 ||
        work[sample] != block) {
      state.SkipWithError("encoded product corpus is not valid");
      return;
    }
    corrupted[sample] = block;
    for (auto& value : corrupted[sample]) {
      const auto before = value;
      for (unsigned bit = 0; bit < 8; ++bit) {
        // Reject the incomplete residue range: exact Bernoulli(1/200),
        // reproducible across standard libraries, with no fixed error count.
        uint32_t draw;
        do {
          draw = static_cast<uint32_t>(channel());
        } while (draw >= 4294967200U);
        if (draw % 200 == 0) {
          value ^= static_cast<Element>(1U << bit);
          ++channel_bits;
        }
      }
      channel_symbols += value != before;
    }
  }

  uint64_t data_residual_bits = 0, codeword_residual_bits = 0;
  uint64_t message_failures = 0, block_failures = 0, validity_failures = 0;
  uint64_t valid_wrong_blocks = 0, pass_limits = 0, passes = 0;
  uint64_t strong_lines = 0, weak_lines = 0;
  uint64_t accepted_bits = 0, accepted_symbols = 0;
  uint64_t strong_bits = 0, strong_symbols = 0, weak_bits = 0, weak_symbols = 0;
  size_t max_passes = 0;
  // Validate this exact channel corpus against the retained single path outside
  // timing, including scheduling outcomes, not just final message recovery.
  for (size_t sample = 0; sample < kCorpusCount; ++sample) {
    auto reference = corrupted[sample];
    const auto expected = gf2p8::rs::detail::ProductCorrectionAccess::Correct(
        code, reference, kPassLimit, 0, false);
    work[sample] = corrupted[sample];
    const auto actual = correct(work[sample]);
    if (work[sample] != reference ||
        actual.termination != expected.termination ||
        actual.all_zero_syndromes != expected.all_zero_syndromes ||
        actual.directional_passes != expected.directional_passes ||
        actual.strong_lines_visited != expected.strong_lines_visited ||
        actual.weak_lines_visited != expected.weak_lines_visited ||
        actual.changed_symbols != expected.changed_symbols ||
        actual.changed_bits != expected.changed_bits ||
        actual.strong_changed_bits != expected.strong_changed_bits ||
        actual.weak_changed_bits != expected.weak_changed_bits ||
        actual.strong_changed_symbols != expected.strong_changed_symbols ||
        actual.weak_changed_symbols != expected.weak_changed_symbols) {
      state.SkipWithError("batch/single corpus differential mismatch");
      return;
    }
  }
  for (auto _ : state) {
    state.PauseTiming();
    for (size_t sample = 0; sample < kCorpusCount; ++sample) {
      std::copy(corrupted[sample].begin(), corrupted[sample].end(),
                work[sample].begin());
    }
    state.ResumeTiming();
    for (size_t sample = 0; sample < kCorpusCount; ++sample) {
      results[sample] = correct(work[sample]);
      benchmark::DoNotOptimize(results[sample]);
      benchmark::ClobberMemory();
    }
    state.PauseTiming();
    for (size_t sample = 0; sample < kCorpusCount; ++sample) {
      const auto& result = results[sample];
      if (result.termination == ProductTermination::invalid_argument) {
        state.SkipWithError("product correction rejected corpus input");
        break;
      }
      uint64_t data_bits = 0, block_bits = 0;
      for (size_t pos = 0; pos < kBlockBytes; ++pos) {
        const auto bits = std::popcount(
            static_cast<unsigned>(work[sample][pos] ^ original[sample][pos]));
        block_bits += bits;
        if (pos / weak_n < kStrongK && pos % weak_n < kWeakK) {
          data_bits += bits;
        }
      }
      data_residual_bits += data_bits;
      codeword_residual_bits += block_bits;
      message_failures += data_bits != 0;
      block_failures += block_bits != 0;
      validity_failures += !result.all_zero_syndromes;
      valid_wrong_blocks += result.all_zero_syndromes && block_bits != 0;
      pass_limits += result.termination == ProductTermination::pass_limit;
      passes += result.directional_passes;
      max_passes = std::max(max_passes, result.directional_passes);
      strong_lines += result.strong_lines_visited;
      weak_lines += result.weak_lines_visited;
      accepted_bits += result.changed_bits;
      accepted_symbols += result.changed_symbols;
      strong_bits += result.strong_changed_bits;
      strong_symbols += result.strong_changed_symbols;
      weak_bits += result.weak_changed_bits;
      weak_symbols += result.weak_changed_symbols;
    }
    state.ResumeTiming();
    if (state.skipped()) {
      break;
    }
  }

  const double sweeps = static_cast<double>(state.iterations());
  if (sweeps == 0 || state.skipped()) {
    return;
  }
  const double blocks = sweeps * kCorpusCount;
  state.counters["corpus_blocks"] = kCorpusCount;
  state.counters["timed_blocks"] = blocks;
  state.counters["seed"] = kSeed;
  state.counters["information_bytes_per_block"] = kInformationBytes;
  state.counters["target_bit_probability"] = 0.005;
  state.counters["channel_codeword_BER"] =
      static_cast<double>(channel_bits) / (kCorpusCount * kBlockBytes * 8);
  state.counters["channel_flipped_bits"] = static_cast<double>(channel_bits);
  state.counters["channel_corrupted_bytes"] =
      static_cast<double>(channel_symbols);
  state.counters["corpus_data_residual_bits"] = data_residual_bits / sweeps;
  state.counters["corpus_codeword_residual_bits"] =
      codeword_residual_bits / sweeps;
  state.counters["data_residual_BER"] =
      data_residual_bits / (blocks * kInformationBytes * 8);
  state.counters["codeword_residual_BER"] =
      codeword_residual_bits / (blocks * kBlockBytes * 8);
  // Counts per unique corpus, not inflated by timed replays of the same noise.
  state.counters["corpus_message_failures"] = message_failures / sweeps;
  state.counters["corpus_block_failures"] = block_failures / sweeps;
  state.counters["corpus_validity_failures"] = validity_failures / sweeps;
  state.counters["corpus_valid_wrong_blocks"] = valid_wrong_blocks / sweeps;
  state.counters["message_recovery_fraction"] = 1 - message_failures / blocks;
  state.counters["block_recovery_fraction"] = 1 - block_failures / blocks;
  state.counters["corpus_pass_limits"] = pass_limits / sweeps;
  state.counters["mean_directional_passes"] = passes / blocks;
  state.counters["max_directional_passes"] = static_cast<double>(max_passes);
  state.counters["mean_strong_lines"] = strong_lines / blocks;
  state.counters["mean_weak_lines"] = weak_lines / blocks;
  state.counters["mean_accepted_bit_changes"] = accepted_bits / blocks;
  state.counters["mean_accepted_byte_changes"] = accepted_symbols / blocks;
  state.counters["mean_strong_accepted_bit_changes"] = strong_bits / blocks;
  state.counters["mean_strong_accepted_byte_changes"] = strong_symbols / blocks;
  state.counters["mean_weak_accepted_bit_changes"] = weak_bits / blocks;
  state.counters["mean_weak_accepted_byte_changes"] = weak_symbols / blocks;
  state.counters["pass_limit"] = kPassLimit;
  state.SetItemsProcessed(state.iterations() * kCorpusCount);
  state.SetBytesProcessed(state.iterations() * kCorpusCount *
                          kInformationBytes);
  state.SetLabel(
      "64 blocks/iteration; information bytes=Kstrong*Kweak; "
      "Correct includes exact final validity; tuned backend");
}

const auto* kProductCorrectionBSC = benchmark::RegisterBenchmark(
    "LCH/Owned/StrongWeakRSProductCode/Correct/BSC005/"
    "Nstrong:256/Kstrong:224/Nweak:256/Kweak:254",
    [](benchmark::State& state) { BenchmarkProductCorrectionBSC(state); });

const auto kProductExperiments = [] {
  for (int variant : {0, 1, 2, 3, 4, 6, 7}) {
    for (size_t n : {256u, 175u}) {
      const auto name = "LCH/Owned/StrongWeakRSProductCode/Experiment/" +
                        std::to_string(variant) + "/" + std::to_string(n);
      benchmark::RegisterBenchmark(name.c_str(), [=](benchmark::State& state) {
        BenchmarkProductCorrectionBSC(state, -1, n, variant);
      });
    }
  }
  return true;
}();

const auto* kProductCorrectionSingle = benchmark::RegisterBenchmark(
    "LCH/Owned/StrongWeakRSProductCode/Correct/BSC005/Single",
    [](benchmark::State& state) { BenchmarkProductCorrectionBSC(state, 0); });
const auto* kProductCorrectionShortened = benchmark::RegisterBenchmark(
    "LCH/Owned/StrongWeakRSProductCode/Correct/BSC005/"
    "Nstrong:256/Kstrong:224/Nweak:175/Kweak:173",
    [](benchmark::State& state) {
      BenchmarkProductCorrectionBSC(state, -1, 175);
    });
const auto* kProductCorrectionStrongBatch = benchmark::RegisterBenchmark(
    "LCH/Owned/StrongWeakRSProductCode/Correct/BSC005/StrongBatch",
    [](benchmark::State& state) { BenchmarkProductCorrectionBSC(state, 1); });
const auto* kProductCorrectionBothBatch = benchmark::RegisterBenchmark(
    "LCH/Owned/StrongWeakRSProductCode/Correct/BSC005/BothBatch",
    [](benchmark::State& state) { BenchmarkProductCorrectionBSC(state, 2); });

const auto* kProductCorrectionGenericShortened = benchmark::RegisterBenchmark(
    "LCH/Owned/StrongWeakRSProductCode/Correct/BSC005/Generic175",
    [](benchmark::State& state) {
      BenchmarkProductCorrectionBSC(state, 2, 175);
    });

void BenchmarkProductEncode(benchmark::State& state) {
  const size_t n = state.range(0);
  gf2p8::rs::StrongWeakRSProductCode code(256, 224, n, n - 2);
  std::vector<gf2p8::Element> block(code.BlockSize());
  std::mt19937 random(0x5b5c0224);
  for (auto& value : block) {
    value = static_cast<gf2p8::Element>(random());
  }
  if (code.Encode(block) != gf2p8::lch::Status::ok) {
    state.SkipWithError("encoding failed");
    return;
  }
  const auto original = block;
  for (auto _ : state) {
    benchmark::DoNotOptimize(code.Encode(block));
    benchmark::ClobberMemory();
  }
  if (block != original || !code.Correct(block).all_zero_syndromes) {
    state.SkipWithError("encoding mismatch");
  }
  state.SetBytesProcessed(state.iterations() * 224 * (n - 2));
}

const auto* kProductEncode =
    benchmark::RegisterBenchmark("LCH/Owned/StrongWeakRSProductCode/Encode",
                                 BenchmarkProductEncode)
        ->Arg(256)
        ->Arg(175);

}  // namespace
