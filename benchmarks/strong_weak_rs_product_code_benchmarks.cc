#include <algorithm>
#include <bit>
#include <cstdint>
#include <random>
#include <vector>

#include "benchmark/benchmark.h"
#include "reed_solomon/strong_weak_rs_product_code.h"
#include "reed_solomon/product_code_internal.h"

namespace {

void BenchmarkProductCorrectionBSC(benchmark::State& state, int batch_passes = -1) {
  using gf2p8::Element;
  using gf2p8::rs::ProductCorrectionResult;
  using gf2p8::rs::ProductTermination;
  constexpr size_t kN = 256;
  constexpr size_t kStrongK = 224;
  constexpr size_t kWeakK = 254;
  constexpr size_t kInformationBytes = kStrongK * kWeakK;
  constexpr size_t kBlockBytes = kN * kN;
  constexpr size_t kCorpusCount = 64;
  constexpr size_t kPassLimit = 16;
  constexpr uint32_t kSeed = 0x5b5c0224;
  gf2p8::rs::StrongWeakRSProductCode code(kN, kStrongK, kN, kWeakK);
  if (!code.Valid() || code.BlockSize() != kBlockBytes) {
    state.SkipWithError("invalid product-code dimensions");
    return;
  }

  std::mt19937 messages(kSeed);
  std::mt19937 channel(kSeed ^ 0x9e3779b9U);
  std::vector<std::vector<Element>> original(
      kCorpusCount, std::vector<Element>(kBlockBytes));
  auto corrupted = original;
  auto work = original;
  std::vector<ProductCorrectionResult> results(kCorpusCount);
  uint64_t channel_bits = 0, channel_symbols = 0;
  for (size_t sample = 0; sample < kCorpusCount; ++sample) {
    auto& block = original[sample];
    for (size_t row = 0; row < kStrongK; ++row) {
      for (size_t col = 0; col < kWeakK; ++col) {
        block[row * kN + col] = static_cast<Element>(messages());
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
  size_t max_passes = 0;
  // Validate this exact channel corpus against the retained single path outside
  // timing, including scheduling outcomes, not just final message recovery.
  for (size_t sample = 0; sample < kCorpusCount; ++sample) {
    auto reference = corrupted[sample];
    const auto expected = gf2p8::rs::detail::ProductCorrectionAccess::Correct(
        code, reference, kPassLimit, 0, false);
    work[sample] = corrupted[sample];
    const auto actual = batch_passes < 0 ? code.Correct(work[sample], kPassLimit)
        : gf2p8::rs::detail::ProductCorrectionAccess::Correct(
              code, work[sample], kPassLimit, batch_passes);
    if (work[sample] != reference || actual.termination != expected.termination ||
        actual.all_zero_syndromes != expected.all_zero_syndromes ||
        actual.directional_passes != expected.directional_passes ||
        actual.strong_lines_visited != expected.strong_lines_visited ||
        actual.weak_lines_visited != expected.weak_lines_visited ||
        actual.changed_symbols != expected.changed_symbols) {
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
      results[sample] = batch_passes < 0 ? code.Correct(work[sample], kPassLimit)
          : gf2p8::rs::detail::ProductCorrectionAccess::Correct(
                code, work[sample], kPassLimit, batch_passes);
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
        const auto bits = std::popcount(static_cast<unsigned>(
            work[sample][pos] ^ original[sample][pos]));
        block_bits += bits;
        if (pos / kN < kStrongK && pos % kN < kWeakK) {
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
    }
    state.ResumeTiming();
    if (state.skipped()) break;
  }

  const double sweeps = static_cast<double>(state.iterations());
  if (sweeps == 0 || state.skipped()) return;
  const double blocks = sweeps * kCorpusCount;
  state.counters["corpus_blocks"] = kCorpusCount;
  state.counters["timed_blocks"] = blocks;
  state.counters["seed"] = kSeed;
  state.counters["information_bytes_per_block"] = kInformationBytes;
  state.counters["target_bit_probability"] = 0.005;
  state.counters["channel_codeword_BER"] =
      static_cast<double>(channel_bits) / (kCorpusCount * kBlockBytes * 8);
  state.counters["channel_flipped_bits"] = static_cast<double>(channel_bits);
  state.counters["channel_corrupted_bytes"] = static_cast<double>(channel_symbols);
  state.counters["corpus_data_residual_bits"] = data_residual_bits / sweeps;
  state.counters["corpus_codeword_residual_bits"] = codeword_residual_bits / sweeps;
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
  state.counters["pass_limit"] = kPassLimit;
  state.SetItemsProcessed(state.iterations() * kCorpusCount);
  state.SetBytesProcessed(state.iterations() * kCorpusCount * kInformationBytes);
  state.SetLabel("64 blocks/iteration; information bytes=224*254; "
                  "Correct includes exact final validity; tuned backend");
}

const auto* kProductCorrectionBSC = benchmark::RegisterBenchmark(
    "LCH/Owned/StrongWeakRSProductCode/Correct/BSC005/"
    "Nstrong:256/Kstrong:224/Nweak:256/Kweak:254",
    [](benchmark::State& state) { BenchmarkProductCorrectionBSC(state); });

const auto* kProductCorrectionSingle = benchmark::RegisterBenchmark(
    "LCH/Owned/StrongWeakRSProductCode/Correct/BSC005/Single",
    [](benchmark::State& state) { BenchmarkProductCorrectionBSC(state, 0); });
const auto* kProductCorrectionStrongBatch = benchmark::RegisterBenchmark(
    "LCH/Owned/StrongWeakRSProductCode/Correct/BSC005/StrongBatch",
    [](benchmark::State& state) { BenchmarkProductCorrectionBSC(state, 1); });
const auto* kProductCorrectionBothBatch = benchmark::RegisterBenchmark(
    "LCH/Owned/StrongWeakRSProductCode/Correct/BSC005/BothBatch",
    [](benchmark::State& state) { BenchmarkProductCorrectionBSC(state, 2); });

}  // namespace
