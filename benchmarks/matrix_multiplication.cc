#include "matrix.h"

#include <benchmark/benchmark.h>
#include <cstdint>
#include <random>
#include <vector>

template <typename T, typename R> void FillRandom(std::vector<T> &v, R &rng) {
  size_t step = 8 / sizeof(T);
  size_t length = v.size() / step;
  for (size_t i = 0; i < length; ++i) {
    uint64_t r = rng();
    *(uint64_t *)(v.data() + i * step) = r;
  }

  for (size_t i = length * step; i < v.size(); ++i) {
    v[i] = rng();
  }
}

template <typename F>
static void RunMatMulBenchmark(benchmark::State &state, F matmul) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<gf_2_8::element_t> left(n * n);
  std::vector<gf_2_8::element_t> right(n * n);
  std::vector<gf_2_8::element_t> result(n * n);

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

static void BM_MatMulBase(benchmark::State &state) {
  gf_2_8::Init();
  RunMatMulBenchmark(state,
                     [](const gf_2_8::element_t *left,
                        const gf_2_8::element_t *right, size_t n,
                        gf_2_8::element_t *result) {
                       gf_2_8::MatMul(left, right, n, n, n,
                                     gf_2_8::AddScaledRowBase, result);
                     });
}

static void BM_MatMulBlockedLUT(benchmark::State &state) {
  gf_2_8::Init();
  RunMatMulBenchmark(state,
                     [](const gf_2_8::element_t *left,
                        const gf_2_8::element_t *right, size_t n,
                        gf_2_8::element_t *result) {
                       gf_2_8::MatMulBlockedLUT(left, right, n, n, n, result);
                     });
}

static void BM_MatMulSIMD(benchmark::State &state) {
  gf_2_8::Init();
  RunMatMulBenchmark(state,
                     [](const gf_2_8::element_t *left,
                        const gf_2_8::element_t *right, size_t n,
                        gf_2_8::element_t *result) {
                       gf_2_8::MatMul(left, right, n, n, n,
                                      gf_2_8::AddScaledRowSIMD, result);
                     });
}

static void BM_MatMulGFNIGeneral(benchmark::State &state) {
  gf_2_8::Init();
  gf_2_8::InitGFNI();
  RunMatMulBenchmark(state,
                     [](const gf_2_8::element_t *left,
                        const gf_2_8::element_t *right, size_t n,
                        gf_2_8::element_t *result) {
                       gf_2_8::MatMul(left, right, n, n, n,
                                      gf_2_8::AddScaledRowGFNIGeneral, result);
                     });
}

static void BM_MatMulGFNIDedicated(benchmark::State &state) {
  gf_2_8::Init();
  gf_2_8::InitGFNI();
  RunMatMulBenchmark(state,
                     [](const gf_2_8::element_t *left,
                        const gf_2_8::element_t *right, size_t n,
                        gf_2_8::element_t *result) {
                       gf_2_8::MatMul(left, right, n, n, n,
                                      gf_2_8::AddScaledRowGFNIDedicated,
                                      result);
                     });
}

static void BM_MatMulBlockedGFNI(benchmark::State &state) {
  gf_2_8::Init();
  gf_2_8::InitGFNI();
  RunMatMulBenchmark(state,
                     [](const gf_2_8::element_t *left,
                        const gf_2_8::element_t *right, size_t n,
                        gf_2_8::element_t *result) {
                       gf_2_8::MatMulBlockedGFNI(left, right, n, n, n, result);
                     });
}

BENCHMARK(BM_MatMulBase)
    ->Name("BinaryTable")
    ->ArgNames({"n"})
    ->RangeMultiplier(2)
    ->Range(16, 2048);

BENCHMARK(BM_MatMulBlockedLUT)
    ->Name("BlockedLowHighSIMD")
    ->ArgNames({"n"})
    ->RangeMultiplier(2)
    ->Range(16, 2048);

BENCHMARK(BM_MatMulSIMD)
    ->Name("LowHighSIMDTables")
    ->ArgNames({"n"})
    ->RangeMultiplier(2)
    ->Range(16, 2048);

BENCHMARK(BM_MatMulGFNIGeneral)
    ->Name("GFNIAffine")
    ->ArgNames({"n"})
    ->RangeMultiplier(2)
    ->Range(16, 2048);

BENCHMARK(BM_MatMulGFNIDedicated)
    ->Name("GFNIMul")
    ->ArgNames({"n"})
    ->RangeMultiplier(2)
    ->Range(16, 2048);

BENCHMARK(BM_MatMulBlockedGFNI)
    ->Name("BlockedGFNI")
    ->ArgNames({"n"})
    ->RangeMultiplier(2)
    ->Range(16, 2048);
