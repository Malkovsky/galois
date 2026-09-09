// 2026-09-08, Ryzen 8845HS/WSL, /check cli Release -march=native.
// CLI: seed=42, k=2600, 2 batches x 4000 trials, default gates/checkpoints;
// all workers verified with affinity 0-15. Whole-simulation wall blocks/s:
// threads       1        2        4        8
// before     634.5   1132.6   2089.3    821.9
// after     1124.4   1907.5   3191.0   1025.0
// Single-worker native phase wall us/block: encode 383 -> 114,
// both difference scans 375 -> 30. Timers removed; seeded metrics unchanged.
// Eight-worker CPU-time inflation remains unresolved on this host.
//
// 2026-09-08 persistent Fisher-Yates experiment, same host and /check cli
// build. seed=42, k=2600, 2x4000 trials, affinity 0-15, default
// gates/checkpoints. Median of 3 rotated-order runs; whole-process wall
// blocks/s incl. final fsync: threads                  1        4        8
// Floyd                 1081.5   2782.4    839.5
// FY diagnostic/no disk   945.9   2662.9    956.3
// FY saved flips          806.8   1959.8    785.8
// Each saved run: 85,120,040 bytes; eight-worker measurements were noisy.
// Retained opt-in only: --sampler fisher-yates always saves replayable flips;
// default Floyd remains unchanged. No-disk FY is not exposed by the CLI.

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <csignal>
#include <cstdint>
#include <exception>
#include <numeric>
#include <random>
#include <vector>

#include "reed_solomon/strong_weak_rs_product_code.h"

namespace {
static_assert(std::atomic<bool>::is_always_lock_free);
std::atomic<bool> interrupted{false};
void Interrupt(int) {
  interrupted = 1;
}

uint64_t Mix(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}
uint64_t Uniform(std::mt19937_64& rng, uint64_t n) {
  const uint64_t threshold = -n % n;
  uint64_t x;
  do {
    x = rng();
  } while (x < threshold);
  return x % n;
}

gf2p8::lch::Status Encode(std::span<gf2p8::Element> block) {
  using gf2p8::Element;
  thread_local gf2p8::rs::LCHEncoder weak(254, 2), strong(224, 32);
  thread_local std::vector<Element> columns(256 * 224);
  thread_local std::vector<Element> workspace(
      std::max(weak.WorkspaceSize(224), strong.WorkspaceSize(256)));
  std::array<const Element*, 254> data{};
  std::array<Element*, 32> recovery{};
  // Independent weak rows become SIMD byte lanes of one shard encode.
  for (size_t col = 0; col < 254; ++col) {
    data[col] = columns.data() + col * 224;
    for (size_t row = 0; row < 224; ++row) {
      columns[col * 224 + row] = block[row * 256 + col];
    }
  }
  recovery[0] = columns.data() + 254 * 224;
  recovery[1] = columns.data() + 255 * 224;
  auto status = weak.Encode(data, std::span(recovery).first(2), 224, workspace);
  if (status != gf2p8::lch::Status::ok) {
    return status;
  }
  for (size_t row = 0; row < 224; ++row) {
    block[row * 256 + 254] = recovery[0][row];
    block[row * 256 + 255] = recovery[1][row];
    data[row] = block.data() + row * 256;
  }
  for (size_t row = 0; row < 32; ++row) {
    recovery[row] = block.data() + (224 + row) * 256;
  }
  return strong.Encode(std::span(data).first(224), recovery, 256, workspace);
}

template <bool Random>
void CountDifferences(const std::vector<gf2p8::Element>& block,
                      const gf2p8::Element* original,
                      uint64_t* out) {
  // Keep the contiguous reduction free of per-byte information-region branches.
  for (size_t row = 0; row < 256; ++row) {
    uint64_t bits = 0, bytes = 0;
    for (size_t col = 0; col < 254; ++col) {
      const unsigned d =
          block[row * 256 + col] ^ (Random ? original[row * 256 + col] : 0);
      bits += std::popcount(d);
      bytes += d != 0;
    }
    out[0] += bits;
    out[1] += bytes;
    if (row < 224) {
      out[2] += bits;
      out[3] += bytes;
    }
    for (size_t col = 254; col < 256; ++col) {
      const unsigned d =
          block[row * 256 + col] ^ (Random ? original[row * 256 + col] : 0);
      out[0] += std::popcount(d);
      out[1] += d != 0;
    }
  }
}
template <bool Random>
int Trial(uint64_t seed,
          uint64_t batch,
          uint64_t trial,
          uint64_t k,
          uint64_t passes,
          int anchors,
          int binary,
          uint64_t* output,
          int sampler,
          uint32_t* positions,
          uint8_t* residual) {
  try {
    // Local counters cannot alias the byte buffers and are published only once.
    uint64_t out[22]{};
    if (!output || k > 524288 || passes < 2 || passes > 1000000 ||
        sampler < 0 || sampler > 2 || (sampler == 2 && !positions)) {
      return 1;
    }
    thread_local gf2p8::rs::StrongWeakRSProductCode code;
    thread_local std::vector<gf2p8::Element> block(65536);
    thread_local std::vector<uint8_t> selected;
    thread_local std::vector<uint32_t> permutation;
    const auto key = Mix(seed ^ Mix(batch) ^ Mix(trial ^ 0x545249414cULL));
    std::mt19937_64 noise(Mix(key ^ 0x4e4f495345ULL));
    const gf2p8::Element* original = nullptr;
    if constexpr (Random) {
      // Legacy replay/reference only. Zero trials neither allocate this buffer
      // nor construct a message PRNG or encoder, even on first worker use.
      thread_local std::vector<gf2p8::Element> reference(65536);
      std::fill(reference.begin(), reference.end(), 0);
      std::mt19937_64 message(Mix(key ^ 0x4d455353414745ULL));
      for (size_t r = 0; r < 224; ++r) {
        for (size_t c = 0; c < 254; ++c) {
          reference[r * 256 + c] = static_cast<uint8_t>(message());
        }
      }
      if (Encode(reference) != gf2p8::lch::Status::ok) {
        return 2;
      }
      original = reference.data();
      block = reference;
    } else {
      // Linearity gives S(c XOR e) = S(e). BDD deltas, anchor decisions and
      // binary-image delta gates therefore depend on e, not the sent codeword.
      std::fill(block.begin(), block.end(), 0);
    }
    // Sample the smaller of the flip set and its complement. Floyd uses a
    // cleared bitmap; opt-in Fisher-Yates retains a 2 MiB worker permutation.
    const bool complement = k > 262144;
    const size_t count = complement ? 524288 - k : k;
    if (sampler != 1) {
      selected.assign(524288, 0);
    } else if (permutation.empty()) {
      permutation.resize(524288);
      std::iota(permutation.begin(), permutation.end(), 0u);
    }
    if (complement) {
      for (auto& byte : block) {
        byte ^= 255;
      }
    }
    for (size_t i = 0; i < count; ++i) {
      size_t pos;
      if (sampler == 1) {
        // Any prior permutation gives a uniform subset. Do not reset it;
        // replay requires saved positions, not just schedule-dependent seeds.
        const size_t draw = i + Uniform(noise, 524288 - i);
        std::swap(permutation[i], permutation[draw]);
        pos = permutation[i];
      } else if (sampler == 2) {
        pos = positions[i];
        if (pos >= 524288 || selected[pos]) {
          return 6;
        }
        selected[pos] = 1;
      } else {
        const size_t j = 524288 - count + i;
        const size_t draw = Uniform(noise, j + 1);
        pos = selected[draw] ? j : draw;
        selected[pos] = 1;
      }
      if (positions && sampler != 2) {
        positions[i] = static_cast<uint32_t>(pos);
      }
      block[pos / 8] ^= static_cast<uint8_t>(1u << (pos % 8));
    }
    CountDifferences<Random>(block, original, out);
    if (out[0] != k) {
      return 4;
    }
    const auto result = code.Correct(
        block, {static_cast<size_t>(passes), anchors != 0, binary != 0});
    if (result.termination == gf2p8::rs::ProductTermination::invalid_argument) {
      return 3;
    }
    CountDifferences<Random>(block, original, out + 4);
    out[8] = out[6] != 0;
    out[9] = out[4] != 0;
    out[10] = result.all_zero_syndromes;
    out[11] = result.all_zero_syndromes && out[9];
    out[12] = result.directional_passes;
    out[13] = result.changed_bits;
    out[14] = result.changed_symbols;
    out[15] = result.strong_changed_bits;
    out[16] = result.strong_changed_symbols;
    out[17] = result.weak_changed_bits;
    out[18] = result.weak_changed_symbols;
    out[19] = result.strong_lines_visited;
    out[20] = result.weak_lines_visited;
    out[21] = result.termination == gf2p8::rs::ProductTermination::pass_limit;
    if (residual) {
      for (size_t i = 0; i < block.size(); ++i) {
        residual[i] = block[i] ^ (Random ? original[i] : 0);
      }
    }
    std::copy(out, out + 22, output);
    return 0;
  } catch (...) {
    return 5;
  }
}
}  // namespace

// Private trial ABI for the native CLI and test-only reference. No C++
// exception crosses the boundary. Counters/residuals publish only after
// successful trials; sampler output positions are scratch and must be ignored
// on failure.
extern "C" {
int product_interrupt_install() {
  interrupted = 0;
  struct sigaction action {};
  action.sa_handler = Interrupt;
  sigemptyset(&action.sa_mask);
  return sigaction(SIGINT, &action, nullptr) ||
         sigaction(SIGTERM, &action, nullptr);
}
int product_interrupted() {
  return interrupted;
}
uint64_t product_batch_k(uint64_t seed,
                         uint64_t batch,
                         uint64_t lo,
                         uint64_t hi) {
  std::mt19937_64 rng(Mix(seed ^ Mix(batch) ^ 0x4241544348ULL));
  return lo + Uniform(rng, hi - lo + 1);
}
// Explicit random reference for old saved records and differential tests.
// Optional residual is the complete decoded block XOR the original codeword.
int product_trial_reference(uint64_t seed,
                            uint64_t batch,
                            uint64_t trial,
                            uint64_t k,
                            uint64_t passes,
                            int anchors,
                            int binary,
                            uint64_t* output,
                            int sampler,
                            uint32_t* positions,
                            int random,
                            uint8_t* residual) {
  if (random != 0 && random != 1) {
    return 1;
  }
  return random ? Trial<true>(seed, batch, trial, k, passes, anchors, binary,
                              output, sampler, positions, residual)
                : Trial<false>(seed, batch, trial, k, passes, anchors, binary,
                               output, sampler, positions, residual);
}
int product_trial_flips(uint64_t seed,
                        uint64_t batch,
                        uint64_t trial,
                        uint64_t k,
                        uint64_t passes,
                        int anchors,
                        int binary,
                        uint64_t* output,
                        int sampler,
                        uint32_t* positions) {
  return Trial<false>(seed, batch, trial, k, passes, anchors, binary, output,
                      sampler, positions, nullptr);
}
int product_trial(uint64_t seed,
                  uint64_t batch,
                  uint64_t trial,
                  uint64_t k,
                  uint64_t passes,
                  int anchors,
                  int binary,
                  uint64_t* output) {
  return product_trial_flips(seed, batch, trial, k, passes, anchors, binary,
                             output, 0, nullptr);
}
// Bounded coordinator task. Publish each successful prefix trial separately;
// check signals between blocks rather than delaying them for the whole task.
int product_trials(uint64_t seed,
                   uint64_t batch,
                   uint64_t first,
                   uint64_t k,
                   uint64_t passes,
                   int anchors,
                   int binary,
                   uint64_t* output,
                   uint64_t count,
                   int sampler,
                   uint32_t* positions,
                   uint64_t* completed) {
  if (!completed) {
    return 1;
  }
  *completed = 0;
  if (!output || count == 0 || count > 4 || first > UINT64_MAX - (count - 1) ||
      k > 524288 || sampler < 0 || sampler > 1 ||
      (sampler == 1 && !positions)) {
    return 1;
  }
  const uint64_t stride = std::min(k, 524288 - k);
  for (uint64_t i = 0; i < count && !interrupted; ++i) {
    const int status = product_trial_flips(
        seed, batch, first + i, k, passes, anchors, binary, output + 22 * i,
        sampler, positions ? positions + stride * i : nullptr);
    if (status) {
      return status;
    }
    ++*completed;
  }
  return 0;
}
}
