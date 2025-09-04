#include "field.h"
#include "gf256/gf256.h"

#include <benchmark/benchmark.h>
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

static void BM_MatMulBase(benchmark::State &state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);
  gf_2_8::Init();

  std::vector<gf_2_8::element_t> left(n * n);
  std::vector<gf_2_8::element_t> right(n * n);
  std::vector<gf_2_8::element_t> result(n * n);

  for (auto _ : state) {
    FillRandom(left, rng);
    FillRandom(right, rng);
    gf_2_8::MatMul(left.data(), right.data(), n, n, n, gf_2_8::AddScaledRowBase,
                   result.data());
  }
}

static void BM_MatMulSIMD(benchmark::State &state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);
  gf256_init();
  

  std::vector<gf_2_8::element_t> left(n * n);
  std::vector<gf_2_8::element_t> right(n * n);
  std::vector<gf_2_8::element_t> result(n * n);

  for (auto _ : state) {
    FillRandom(left, rng);
    FillRandom(right, rng);
    gf_2_8::MatMul(left.data(), right.data(), n, n, n, gf_2_8::AddScaledRowSIMD,
                   result.data());
  }
}

static void BM_MatMulGFNIGeneral(benchmark::State &state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);
  gf_2_8::InitGFNI();

  std::vector<gf_2_8::element_t> left(n * n);
  std::vector<gf_2_8::element_t> right(n * n);
  std::vector<gf_2_8::element_t> result(n * n);

  for (auto _ : state) {
    FillRandom(left, rng);
    FillRandom(right, rng);
    gf_2_8::MatMul(left.data(), right.data(), n, n, n,
                   gf_2_8::AddScaledRowGFNIGeneral, result.data());
  }
}

static void BM_MatMulGFNIDedicated(benchmark::State &state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<gf_2_8::element_t> left(n * n);
  std::vector<gf_2_8::element_t> right(n * n);
  std::vector<gf_2_8::element_t> result(n * n);

  for (auto _ : state) {
    FillRandom(left, rng);
    FillRandom(right, rng);
    gf_2_8::MatMul(left.data(), right.data(), n, n, n,
                   gf_2_8::AddScaledRowGFNIDedicated, result.data());
  }
}

BENCHMARK(BM_MatMulBase)
    ->Name("BinaryTable")
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
