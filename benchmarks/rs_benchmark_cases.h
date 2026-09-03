#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace gf256_benchmarks {

struct RsBenchmarkCase {
  int64_t data_count;
  int64_t recovery_count;
  int64_t bytes;
};

// Full-length low-/high-rate comparison grid used by the owned LCH codec and
// the native XDRS and ISA-L references.
inline constexpr std::array<RsBenchmarkCase, 9> kLCHComparisonCases{{
    {8, 248, 1024},
    {16, 240, 1024},
    {32, 224, 1024},
    {64, 192, 1024},
    {128, 128, 1024},
    {192, 64, 1024},
    {224, 32, 1024},
    {240, 16, 1024},
    {248, 8, 1024},
}};

// Native Leopard supports only R<=K, so its rows use the high-rate subset.
inline constexpr std::array<RsBenchmarkCase, 5> kNativeLeopardCases{{
    {128, 128, 1024},
    {192, 64, 1024},
    {224, 32, 1024},
    {240, 16, 1024},
    {248, 8, 1024},
}};

constexpr bool NativeXDRSLowRate(int64_t data_count) {
  return data_count <= 128;
}

using ErasurePattern = std::array<uint8_t, 256>;

// Canonical logical order is [recovery][data]. This deterministically
// reproduces the upstream benchmark's shuffled maximum-erasure workload
// without timing random-number generation.
inline ErasurePattern MaxErasurePattern(size_t data_count,
                                        size_t recovery_count) {
  ErasurePattern erased{};
  for (size_t i = 0; i < recovery_count; ++i) {
    erased[i] = 1;
  }

  uint32_t random = 0x6d2b79f5U ^ static_cast<uint32_t>(data_count * 257U);
  const auto next_random = [&random] {
    random ^= random << 13;
    random ^= random >> 17;
    random ^= random << 5;
    return random;
  };
  for (size_t i = data_count + recovery_count - 1; i > 0; --i) {
    const size_t position = next_random() % (i + 1);
    std::swap(erased[i], erased[position]);
  }
  return erased;
}

inline ErasurePattern XDRSMaxErasurePattern(size_t data_count,
                                            size_t recovery_count) {
  const ErasurePattern canonical =
      MaxErasurePattern(data_count, recovery_count);
  if (!NativeXDRSLowRate(static_cast<int64_t>(data_count))) {
    return canonical;
  }

  ErasurePattern active_layout{};
  for (size_t i = 0; i < data_count; ++i) {
    active_layout[i] = canonical[recovery_count + i];
  }
  for (size_t i = 0; i < recovery_count; ++i) {
    active_layout[data_count + i] = canonical[i];
  }
  return active_layout;
}

inline size_t DecodeInputBytesPerIteration(size_t data_count, size_t bytes) {
  return data_count * bytes;
}

}  // namespace gf256_benchmarks
