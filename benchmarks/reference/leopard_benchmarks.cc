#include <benchmark/benchmark.h>
#include <leopard.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "../rs_benchmark_cases.h"
#include "reference_benchmarks.h"

namespace gf256_benchmarks {
namespace {

std::vector<std::vector<uint8_t>> MakeBuffers(size_t count, size_t bytes) {
  std::mt19937 random(42);
  std::vector<std::vector<uint8_t>> result(count, std::vector<uint8_t>(bytes));
  for (auto& buffer : result) {
    std::generate(buffer.begin(), buffer.end(),
                  [&random] { return static_cast<uint8_t>(random()); });
  }
  return result;
}

bool VerifyCodeword(unsigned k,
                    unsigned r,
                    size_t bytes,
                    const std::vector<std::vector<uint8_t>>& original,
                    const std::vector<std::vector<uint8_t>>& encode_work) {
  std::vector<const void*> original_pointers(k);
  for (size_t i = 1; i < original.size(); ++i) {
    original_pointers[i] = original[i].data();
  }
  std::vector<const void*> recovery_pointers(r);
  for (size_t i = 0; i < recovery_pointers.size(); ++i) {
    recovery_pointers[i] = encode_work[i].data();
  }
  const unsigned decode_work_count = leo_decode_work_count(k, r);
  auto decode_work = MakeBuffers(decode_work_count, bytes);
  std::vector<void*> decode_work_pointers(decode_work_count);
  for (size_t i = 0; i < decode_work.size(); ++i) {
    decode_work_pointers[i] = decode_work[i].data();
  }
  return leo_decode(bytes, k, r, decode_work_count, original_pointers.data(),
                    recovery_pointers.data(),
                    decode_work_pointers.data()) == Leopard_Success &&
         decode_work[0] == original[0];
}

void NativeLeopardEncode(benchmark::State& state) {
  const unsigned k = static_cast<unsigned>(state.range(0));
  const unsigned r = static_cast<unsigned>(state.range(1));
  const size_t bytes = static_cast<size_t>(state.range(2));
  const unsigned work_count = leo_encode_work_count(k, r);
  auto original = MakeBuffers(k, bytes);
  auto work = MakeBuffers(work_count, bytes);
  std::vector<const void*> original_pointers(k);
  std::vector<void*> work_pointers(work_count);
  for (size_t i = 0; i < original.size(); ++i) {
    original_pointers[i] = original[i].data();
  }
  for (size_t i = 0; i < work.size(); ++i) {
    work_pointers[i] = work[i].data();
  }
  if (leo_encode(bytes, k, r, work_count, original_pointers.data(),
                 work_pointers.data()) != Leopard_Success ||
      !VerifyCodeword(k, r, bytes, original, work)) {
    state.SkipWithError("native Leopard encode verification failed");
    return;
  }

  for (auto _ : state) {
    LeopardResult result =
        leo_encode(bytes, k, r, work_count, original_pointers.data(),
                   work_pointers.data());
    benchmark::DoNotOptimize(result);
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * k *
                          static_cast<int64_t>(bytes));
}

void NativeLeopardDecode(benchmark::State& state) {
  const unsigned k = static_cast<unsigned>(state.range(0));
  const unsigned r = static_cast<unsigned>(state.range(1));
  const size_t bytes = static_cast<size_t>(state.range(2));
  const unsigned encode_work_count = leo_encode_work_count(k, r);
  auto original = MakeBuffers(k, bytes);
  auto encode_work = MakeBuffers(encode_work_count, bytes);
  std::vector<const void*> original_pointers(k);
  std::vector<void*> encode_work_pointers(encode_work_count);
  for (size_t i = 0; i < original.size(); ++i) {
    original_pointers[i] = original[i].data();
  }
  for (size_t i = 0; i < encode_work.size(); ++i) {
    encode_work_pointers[i] = encode_work[i].data();
  }
  if (leo_encode(bytes, k, r, encode_work_count, original_pointers.data(),
                 encode_work_pointers.data()) != Leopard_Success) {
    state.SkipWithError("native Leopard encode setup failed");
    return;
  }

  const auto erased = MaxErasurePattern(k, r);
  size_t missing_data_count = 0;
  for (size_t i = 0; i < original.size(); ++i) {
    if (erased[r + i]) {
      original_pointers[i] = nullptr;
      ++missing_data_count;
    }
  }
  std::vector<const void*> recovery_pointers(r);
  for (size_t i = 0; i < recovery_pointers.size(); ++i) {
    recovery_pointers[i] = erased[i] ? nullptr : encode_work[i].data();
  }
  const unsigned decode_work_count = leo_decode_work_count(k, r);
  auto decode_work = MakeBuffers(decode_work_count, bytes);
  std::vector<void*> decode_work_pointers(decode_work_count);
  for (size_t i = 0; i < decode_work.size(); ++i) {
    decode_work_pointers[i] = decode_work[i].data();
  }
  if (leo_decode(bytes, k, r, decode_work_count, original_pointers.data(),
                 recovery_pointers.data(),
                 decode_work_pointers.data()) != Leopard_Success ||
      std::count(erased.begin(), erased.end(), uint8_t{1}) != r ||
      missing_data_count == 0) {
    state.SkipWithError("native Leopard decode verification failed");
    return;
  }
  for (size_t i = 0; i < original.size(); ++i) {
    if (erased[r + i] && decode_work[i] != original[i]) {
      state.SkipWithError("native Leopard decode verification failed");
      return;
    }
  }

  for (auto _ : state) {
    LeopardResult result =
        leo_decode(bytes, k, r, decode_work_count, original_pointers.data(),
                   recovery_pointers.data(), decode_work_pointers.data());
    benchmark::DoNotOptimize(result);
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(
      static_cast<int64_t>(state.iterations()) *
      static_cast<int64_t>(DecodeInputBytesPerIteration(k, bytes)));
}

}  // namespace

void RegisterLeopardReferenceBenchmarks() {
  if (leo_init() != 0) {
    return;
  }
#if defined(GF256_HIGH_LEVEL_BENCHMARKS_ONLY)
  for (const auto& [name, function] :
       std::array<std::pair<const char*, void (*)(benchmark::State&)>, 2>{{
           {"RS/Leopard/Native/Encode", &NativeLeopardEncode},
           {"RS/Leopard/Native/DecodeMax", &NativeLeopardDecode},
       }}) {
    auto* registered = benchmark::RegisterBenchmark(name, function);
    for (const auto& test_case : kNativeLeopardCases) {
      registered->Args(
          {test_case.data_count, test_case.recovery_count, test_case.bytes});
    }
    registered->ArgNames({"K", "R", "bytes"});
  }
#else
  for (const auto& [name, function] :
       std::array<std::pair<const char*, void (*)(benchmark::State&)>, 2>{
           {{"Leopard/Native/Encode", &NativeLeopardEncode},
            {"Leopard/Native/DecodeMax", &NativeLeopardDecode}}}) {
    auto* registered = benchmark::RegisterBenchmark(name, function);
    registered->Args({32, 16, 1024})
        ->Args({129, 64, 1024})
        ->Args({128, 128, 1024})
        ->ArgNames({"K", "R", "bytes"});
  }
#endif
}

}  // namespace gf256_benchmarks
