#include <benchmark/benchmark.h>

#include <cstdint>
#include <random>
#include <vector>

#include "matrix.h"

template <typename T, typename R>
void FillRandom(std::vector<T>& v, R& rng) {
  size_t step = 8 / sizeof(T);
  size_t length = v.size() / step;
  for (size_t i = 0; i < length; ++i) {
    uint64_t r = rng();
    *(uint64_t*)(v.data() + i * step) = r;
  }

  for (size_t i = length * step; i < v.size(); ++i) {
    v[i] = rng();
  }
}

template <typename F>
static void RunMatMulBenchmark(benchmark::State& state, F matmul) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<gf2p8::Element> left(n * n);
  std::vector<gf2p8::Element> right(n * n);
  std::vector<gf2p8::Element> result(n * n);

  FillRandom(left, rng);
  FillRandom(right, rng);

  for (auto _ : state) {
    matmul(left.data(), right.data(), n, result.data());
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }

  state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(n) * static_cast<int64_t>(n) *
                          static_cast<int64_t>(n));
}

static void BM_MatMulBase(benchmark::State& state) {
  RunMatMulBenchmark(
      state, [](const gf2p8::Element* left, const gf2p8::Element* right,
                size_t n, gf2p8::Element* result) {
        gf2p8::MatMul(left, right, n, n, n, gf2p8::AddScaledRow, result);
      });
}

static void BM_MatMulBlockedLUT(benchmark::State& state) {
  RunMatMulBenchmark(state,
                     [](const gf2p8::Element* left, const gf2p8::Element* right,
                        size_t n, gf2p8::Element* result) {
                       gf2p8::MatMulBlockedLUT(left, right, n, n, n, result);
                     });
}

static void BM_MatMulBlockedGFNI(benchmark::State& state) {
  RunMatMulBenchmark(state,
                     [](const gf2p8::Element* left, const gf2p8::Element* right,
                        size_t n, gf2p8::Element* result) {
                       gf2p8::MatMulBlockedGFNI(left, right, n, n, n, result);
                     });
}

BENCHMARK(BM_MatMulBase)
    ->Name("SharedShuffleRows")
    ->ArgNames({"n"})
    ->RangeMultiplier(2)
    ->Range(16, 2048);

BENCHMARK(BM_MatMulBlockedLUT)
    ->Name("BlockedLowHighSIMD")
    ->ArgNames({"n"})
    ->RangeMultiplier(2)
    ->Range(16, 2048);

BENCHMARK(BM_MatMulBlockedGFNI)
    ->Name("BlockedGFNIAffine")
    ->ArgNames({"n"})
    ->RangeMultiplier(2)
    ->Range(16, 2048);
