#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "reed_solomon/leopard.h"
#include "reed_solomon/lin_chung_han/transform.h"
#include "reed_solomon/xdrs.h"
#include "rs_benchmark_cases.h"

#if defined(__linux__)
#include <sched.h>
#endif

#if defined(GF256_HAVE_REFERENCE_BENCHMARKS)
#include "reference/reference_benchmarks.h"
#endif

#ifndef GF256_GIT_REVISION
#define GF256_GIT_REVISION "unknown"
#endif
#ifndef GF256_BUILD_TYPE
#define GF256_BUILD_TYPE "unknown"
#endif
#ifndef GF256_COMPILER_ID
#define GF256_COMPILER_ID "unknown"
#endif
#ifndef GF256_COMPILER_VERSION
#define GF256_COMPILER_VERSION "unknown"
#endif
#ifndef GF256_NATIVE_ISA
#define GF256_NATIVE_ISA "unknown"
#endif

namespace {

using gf2p8::Element;
using gf2p8::lch::Backend;
using gf2p8::lch::Context;
using gf2p8::lch::Radix;
using gf2p8::lch::Status;

const char* BackendName(Backend backend) {
  switch (backend) {
    case Backend::tuned:
      return "Tuned";
    case Backend::scalar:
      return "Scalar";
    case Backend::ssse3:
      return "SSSE3";
    case Backend::avx2:
      return "AVX2";
    case Backend::gfni128_affine:
      return "GFNI128Affine";
    case Backend::gfni256_affine:
      return "GFNI256Affine";
    case Backend::gfni512_affine:
      return "GFNI512Affine";
  }
  return "Unknown";
}

std::vector<Element*> MutablePointers(
    std::vector<std::vector<Element>>& shards) {
  std::vector<Element*> result;
  result.reserve(shards.size());
  for (auto& shard : shards) {
    result.push_back(shard.data());
  }
  return result;
}

std::vector<const Element*> ConstPointers(
    const std::vector<std::vector<Element>>& shards) {
  std::vector<const Element*> result;
  result.reserve(shards.size());
  for (const auto& shard : shards) {
    result.push_back(shard.data());
  }
  return result;
}

std::vector<std::vector<Element>> MakeShards(size_t count, size_t bytes) {
  std::mt19937 random(42);
  std::vector<std::vector<Element>> result(count, std::vector<Element>(bytes));
  for (auto& shard : result) {
    std::generate(shard.begin(), shard.end(),
                  [&random] { return static_cast<Element>(random()); });
  }
  return result;
}

void LCHBenchmark(benchmark::State& state,
                  Backend backend,
                  Radix radix,
                  bool inverse) {
  if (!gf2p8::lch::BackendAvailable(backend)) {
    state.SkipWithError("backend was not compiled");
    return;
  }
  const size_t transform_size = static_cast<size_t>(state.range(0));
  const size_t bytes = static_cast<size_t>(state.range(1));
  const Context& context = Context::Shared();
  auto shards = MakeShards(transform_size, bytes);
  auto pointers = MutablePointers(shards);
  const gf2p8::lch::TransformOptions options{backend, radix};

  for (auto _ : state) {
    Status status = inverse
                        ? gf2p8::lch::IFFT(context, pointers, bytes, 0, options)
                        : gf2p8::lch::FFT(context, pointers, bytes, 0, options);
    benchmark::DoNotOptimize(status);
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(transform_size) *
                          static_cast<int64_t>(bytes));
}

void LeopardEncodeBenchmark(benchmark::State& state,
                            Backend backend,
                            Radix radix) {
  if (!gf2p8::lch::BackendAvailable(backend)) {
    state.SkipWithError("backend was not compiled");
    return;
  }
  const size_t original_count = static_cast<size_t>(state.range(0));
  const size_t recovery_count = static_cast<size_t>(state.range(1));
  const size_t bytes = static_cast<size_t>(state.range(2));
  const gf2p8::rs::Leopard code(original_count, recovery_count);
  auto original = MakeShards(original_count, bytes);
  auto recovery = MakeShards(recovery_count, bytes);
  const auto original_pointers = ConstPointers(original);
  auto recovery_pointers = MutablePointers(recovery);
  std::vector<Element> workspace(code.EncodeWorkspaceSize(bytes));
  auto expected_recovery = MakeShards(recovery_count, bytes);
  auto expected_pointers = MutablePointers(expected_recovery);
  std::vector<Element> expected_workspace(code.EncodeWorkspaceSize(bytes));
  if (code.Encode(original_pointers, expected_pointers, bytes,
                  expected_workspace, Backend::scalar,
                  Radix::radix2) != Status::ok ||
      code.Encode(original_pointers, recovery_pointers, bytes, workspace,
                  backend, radix) != Status::ok ||
      recovery != expected_recovery) {
    state.SkipWithError("Leopard encode verification failed");
    return;
  }

  for (auto _ : state) {
    Status status = code.Encode(original_pointers, recovery_pointers, bytes,
                                workspace, backend, radix);
    benchmark::DoNotOptimize(status);
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(original_count) *
                          static_cast<int64_t>(bytes));
}

void LeopardDecodeBenchmark(benchmark::State& state,
                            Backend backend,
                            Radix radix) {
  if (!gf2p8::lch::BackendAvailable(backend)) {
    state.SkipWithError("backend was not compiled");
    return;
  }
  const size_t original_count = static_cast<size_t>(state.range(0));
  const size_t recovery_count = static_cast<size_t>(state.range(1));
  const size_t bytes = static_cast<size_t>(state.range(2));
  const gf2p8::rs::Leopard code(original_count, recovery_count);
  auto original = MakeShards(original_count, bytes);
  auto recovery = MakeShards(recovery_count, bytes);
  const auto original_pointers = ConstPointers(original);
  auto recovery_pointers = MutablePointers(recovery);
  std::vector<Element> encode_workspace(code.EncodeWorkspaceSize(bytes));
  if (code.Encode(original_pointers, recovery_pointers, bytes, encode_workspace,
                  backend, radix) != Status::ok) {
    state.SkipWithError("encode setup failed");
    return;
  }

  std::vector<Element*> codeword;
  codeword.reserve(original_count + recovery_count);
  codeword.insert(codeword.end(), recovery_pointers.begin(),
                  recovery_pointers.end());
  auto mutable_original = MutablePointers(original);
  codeword.insert(codeword.end(), mutable_original.begin(),
                  mutable_original.end());
  const auto erased =
      gf256_benchmarks::MaxErasurePattern(original_count, recovery_count);
  std::vector<uint8_t> present(codeword.size());
  for (size_t i = 0; i < present.size(); ++i) {
    present[i] = static_cast<uint8_t>(!erased[i]);
  }
  std::vector<Element> workspace(code.DecodeWorkspaceSize(bytes));
  const auto expected_original = original;
  size_t missing_data_count = 0;
  for (size_t i = 0; i < original_count; ++i) {
    if (erased[recovery_count + i]) {
      std::fill(original[i].begin(), original[i].end(), 0);
      ++missing_data_count;
    }
  }
  if (code.Decode(codeword, present, bytes, workspace, backend, radix) !=
          Status::ok ||
      std::count(erased.begin(), erased.end(), uint8_t{1}) != recovery_count ||
      missing_data_count == 0) {
    state.SkipWithError("Leopard decode verification failed");
    return;
  }
  for (size_t i = 0; i < original_count; ++i) {
    if (erased[recovery_count + i] && original[i] != expected_original[i]) {
      state.SkipWithError("Leopard decode verification failed");
      return;
    }
  }

  for (auto _ : state) {
    Status status =
        code.Decode(codeword, present, bytes, workspace, backend, radix);
    benchmark::DoNotOptimize(status);
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(
      static_cast<int64_t>(state.iterations()) *
      static_cast<int64_t>(gf256_benchmarks::DecodeInputBytesPerIteration(
          original_count, bytes)));
}

void XDRSEncodeBenchmark(benchmark::State& state,
                         Backend backend,
                         gf2p8::rs::XDRSRate rate,
                         Radix radix,
                         size_t bytes_range = 1) {
  if (!gf2p8::lch::BackendAvailable(backend)) {
    state.SkipWithError("backend was not compiled");
    return;
  }
  const size_t data_count = static_cast<size_t>(state.range(0));
  const size_t bytes = static_cast<size_t>(state.range(bytes_range));
  const gf2p8::rs::XDRS code(data_count, rate);
  auto data = MakeShards(data_count, bytes);
  auto recovery = MakeShards(code.RecoveryCount(), bytes);
  const auto data_pointers = ConstPointers(data);
  auto recovery_pointers = MutablePointers(recovery);
  std::vector<Element> workspace(code.EncodeWorkspaceSize(bytes));

  for (auto _ : state) {
    Status status = code.Encode(data_pointers, recovery_pointers, bytes,
                                workspace, backend, radix);
    benchmark::DoNotOptimize(status);
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(data_count) *
                          static_cast<int64_t>(bytes));
}

void XDRSDecodeBenchmark(benchmark::State& state,
                         Backend backend,
                         gf2p8::rs::XDRSRate rate,
                         Radix radix,
                         size_t bytes_range = 1) {
  if (!gf2p8::lch::BackendAvailable(backend)) {
    state.SkipWithError("backend was not compiled");
    return;
  }
  const size_t data_count = static_cast<size_t>(state.range(0));
  const size_t bytes = static_cast<size_t>(state.range(bytes_range));
  const gf2p8::rs::XDRS code(data_count, rate);
  auto data = MakeShards(data_count, bytes);
  auto recovery = MakeShards(code.RecoveryCount(), bytes);
  const auto data_pointers = ConstPointers(data);
  auto recovery_pointers = MutablePointers(recovery);
  std::vector<Element> encode_workspace(code.EncodeWorkspaceSize(bytes));
  if (code.Encode(data_pointers, recovery_pointers, bytes, encode_workspace,
                  backend, radix) != Status::ok) {
    state.SkipWithError("encode setup failed");
    return;
  }

  auto mutable_data = MutablePointers(data);
  std::vector<Element*> codeword;
  codeword.reserve(Context::kFieldSize);
  if (rate == gf2p8::rs::XDRSRate::low) {
    codeword.insert(codeword.end(), mutable_data.begin(), mutable_data.end());
    codeword.insert(codeword.end(), recovery_pointers.begin(),
                    recovery_pointers.end());
  } else {
    codeword.insert(codeword.end(), recovery_pointers.begin(),
                    recovery_pointers.end());
    codeword.insert(codeword.end(), mutable_data.begin(), mutable_data.end());
  }
  const size_t recovery_count = code.RecoveryCount();
  const auto erased =
      gf256_benchmarks::XDRSMaxErasurePattern(data_count, recovery_count);
  std::vector<uint8_t> present(codeword.size());
  for (size_t i = 0; i < present.size(); ++i) {
    present[i] = static_cast<uint8_t>(!erased[i]);
  }
  std::vector<Element> workspace(code.DecodeWorkspaceSize(bytes));
  const auto expected_data = data;
  size_t missing_data_count = 0;
  for (size_t i = 0; i < data_count; ++i) {
    const size_t codeword_index =
        rate == gf2p8::rs::XDRSRate::low ? i : recovery_count + i;
    if (erased[codeword_index]) {
      std::fill(data[i].begin(), data[i].end(), 0);
      ++missing_data_count;
    }
  }
  if (code.Decode(codeword, present, bytes, workspace, backend, radix) !=
          Status::ok ||
      std::count(erased.begin(), erased.end(), uint8_t{1}) != recovery_count ||
      missing_data_count == 0) {
    state.SkipWithError("XDRS decode verification failed");
    return;
  }
  for (size_t i = 0; i < data_count; ++i) {
    const size_t codeword_index =
        rate == gf2p8::rs::XDRSRate::low ? i : recovery_count + i;
    if (erased[codeword_index] && data[i] != expected_data[i]) {
      state.SkipWithError("XDRS decode verification failed");
      return;
    }
  }

  for (auto _ : state) {
    Status status =
        code.Decode(codeword, present, bytes, workspace, backend, radix);
    benchmark::DoNotOptimize(status);
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(
      static_cast<int64_t>(state.iterations()) *
      static_cast<int64_t>(
          gf256_benchmarks::DecodeInputBytesPerIteration(data_count, bytes)));
}

void RegisterLCHBenchmarks() {
  const Backend backends[] = {
      Backend::scalar,         Backend::ssse3,          Backend::avx2,
      Backend::gfni128_affine, Backend::gfni256_affine, Backend::gfni512_affine,
  };
  const int64_t lengths[] = {1,    8,    15,   16,   17,    31,    32,
                             33,   63,   64,   65,   128,   256,   512,
                             1024, 2048, 4096, 8192, 16384, 32768, 65536};
  for (const Backend backend : backends) {
    for (const Radix radix : {Radix::radix2, Radix::radix4}) {
      for (const bool inverse : {false, true}) {
        const std::string name =
            std::string("LCH/") + (inverse ? "IFFT/" : "FFT/") +
            BackendName(backend) +
            (radix == Radix::radix2 ? "/Radix2" : "/Radix4");
        auto* benchmark = benchmark::RegisterBenchmark(
            name, [backend, radix, inverse](benchmark::State& state) {
              LCHBenchmark(state, backend, radix, inverse);
            });
        for (int64_t size = 2; size <= 256; size *= 2) {
          for (const int64_t length : lengths) {
            benchmark->Args({size, length});
          }
        }
        benchmark->ArgNames({"transform", "bytes"});
      }
    }
  }
}

void RegisterRsBenchmarks() {
  const Backend backends[] = {Backend::scalar,         Backend::ssse3,
                              Backend::avx2,           Backend::gfni128_affine,
                              Backend::gfni256_affine, Backend::gfni512_affine};
  for (const Backend backend : backends) {
    for (const Radix radix : {Radix::radix2, Radix::radix4}) {
      const char* const radix_name =
          radix == Radix::radix2 ? "Radix2" : "Radix4";
      auto* leopard = benchmark::RegisterBenchmark(
          std::string("Leopard/Polished/Encode/") + BackendName(backend) + "/" +
              radix_name,
          [backend, radix](benchmark::State& state) {
            LeopardEncodeBenchmark(state, backend, radix);
          });
      leopard->Args({32, 16, 1024})
          ->Args({129, 64, 1024})
          ->Args({128, 128, 1024})
          ->ArgNames({"K", "R", "bytes"});

      auto* leopard_decode = benchmark::RegisterBenchmark(
          std::string("Leopard/Polished/DecodeMax/") + BackendName(backend) +
              "/" + radix_name,
          [backend, radix](benchmark::State& state) {
            LeopardDecodeBenchmark(state, backend, radix);
          });
      leopard_decode->Args({32, 16, 1024})
          ->Args({129, 64, 1024})
          ->Args({128, 128, 1024})
          ->ArgNames({"K", "R", "bytes"});

      auto* low = benchmark::RegisterBenchmark(
          std::string("XDRS/Polished/Low/Encode/") + BackendName(backend) +
              "/" + radix_name,
          [backend, radix](benchmark::State& state) {
            XDRSEncodeBenchmark(state, backend, gf2p8::rs::XDRSRate::low,
                                radix);
          });
      for (const int64_t k : {8, 16, 32, 64, 128}) {
        low->Args({k, 1024});
      }
      low->ArgNames({"K", "bytes"});

      auto* low_decode = benchmark::RegisterBenchmark(
          std::string("XDRS/Polished/Low/DecodeMax/") + BackendName(backend) +
              "/" + radix_name,
          [backend, radix](benchmark::State& state) {
            XDRSDecodeBenchmark(state, backend, gf2p8::rs::XDRSRate::low,
                                radix);
          });
      for (const int64_t k : {8, 16, 32, 64, 128}) {
        low_decode->Args({k, 1024});
      }
      low_decode->ArgNames({"K", "bytes"});

      auto* high = benchmark::RegisterBenchmark(
          std::string("XDRS/Polished/High/Encode/") + BackendName(backend) +
              "/" + radix_name,
          [backend, radix](benchmark::State& state) {
            XDRSEncodeBenchmark(state, backend, gf2p8::rs::XDRSRate::high,
                                radix);
          });
      for (const int64_t recovery : {128, 64, 32, 16, 8}) {
        high->Args({256 - recovery, 1024});
      }
      high->ArgNames({"K", "bytes"});

      auto* high_decode = benchmark::RegisterBenchmark(
          std::string("XDRS/Polished/High/DecodeMax/") + BackendName(backend) +
              "/" + radix_name,
          [backend, radix](benchmark::State& state) {
            XDRSDecodeBenchmark(state, backend, gf2p8::rs::XDRSRate::high,
                                radix);
          });
      for (const int64_t recovery : {128, 64, 32, 16, 8}) {
        high_decode->Args({256 - recovery, 1024});
      }
      high_decode->ArgNames({"K", "bytes"});
    }
  }
}

void RegisterHighLevelRsBenchmarks() {
  constexpr Radix radix = Radix::radix4;

  const auto register_leopard = [](const char* name, Backend selected_backend,
                                   bool decode) {
    auto* registered = benchmark::RegisterBenchmark(
        name, [selected_backend, decode](benchmark::State& state) {
          if (decode) {
            LeopardDecodeBenchmark(state, selected_backend, radix);
          } else {
            LeopardEncodeBenchmark(state, selected_backend, radix);
          }
        });
    for (const auto& test_case : gf256_benchmarks::kLeopardComparisonCases) {
      registered->Args(
          {test_case.data_count, test_case.recovery_count, test_case.bytes});
    }
    registered->ArgNames({"K", "R", "bytes"});
  };

  register_leopard("RS/Leopard/Owned/AVX2/Encode", Backend::avx2, false);
  register_leopard("RS/Leopard/Owned/AVX2/DecodeMax", Backend::avx2, true);
  register_leopard("RS/Leopard/Owned/GFNI256/Encode", Backend::gfni256_affine,
                   false);
  register_leopard("RS/Leopard/Owned/GFNI256/DecodeMax",
                   Backend::gfni256_affine, true);

  const auto register_xdrs = [](const char* name, Backend selected_backend,
                                bool decode) {
    auto* registered = benchmark::RegisterBenchmark(
        name, [selected_backend, decode](benchmark::State& state) {
          const gf2p8::rs::XDRSRate rate =
              gf256_benchmarks::XDRSLowRate(state.range(0))
                  ? gf2p8::rs::XDRSRate::low
                  : gf2p8::rs::XDRSRate::high;
          if (decode) {
            XDRSDecodeBenchmark(state, selected_backend, rate, radix, 2);
          } else {
            XDRSEncodeBenchmark(state, selected_backend, rate, radix, 2);
          }
        });
    for (const auto& test_case : gf256_benchmarks::kXDRSLogCases) {
      registered->Args(
          {test_case.data_count, test_case.recovery_count, test_case.bytes});
    }
    registered->ArgNames({"K", "R", "bytes"});
  };

  register_xdrs("RS/XDRS/Owned/AVX2/Encode", Backend::avx2, false);
  register_xdrs("RS/XDRS/Owned/AVX2/DecodeMax", Backend::avx2, true);
  register_xdrs("RS/XDRS/Owned/GFNI256/Encode", Backend::gfni256_affine, false);
  register_xdrs("RS/XDRS/Owned/GFNI256/DecodeMax", Backend::gfni256_affine,
                true);
}

std::string CpuPinningStatus() {
#if defined(__linux__)
  cpu_set_t affinity;
  CPU_ZERO(&affinity);
  if (sched_getaffinity(0, sizeof(affinity), &affinity) == 0) {
    return CPU_COUNT(&affinity) == 1 ? "pinned" : "not_pinned";
  }
#endif
  return "unknown";
}

}  // namespace

int main(int argc, char** argv) {
  benchmark::AddCustomContext("git_revision", GF256_GIT_REVISION);
  benchmark::AddCustomContext("build_type", GF256_BUILD_TYPE);
  benchmark::AddCustomContext("compiler", std::string(GF256_COMPILER_ID) + " " +
                                              GF256_COMPILER_VERSION);
  benchmark::AddCustomContext("native_isa", GF256_NATIVE_ISA);
  benchmark::AddCustomContext("cpu_pinning", CpuPinningStatus());
#if defined(__GFNI__)
  benchmark::AddCustomContext("compiled_gfni", "true");
#else
  benchmark::AddCustomContext("compiled_gfni", "false");
#endif
#if defined(GF256_HIGH_LEVEL_BENCHMARKS_ONLY)
  RegisterHighLevelRsBenchmarks();
#else
  RegisterLCHBenchmarks();
  RegisterRsBenchmarks();
#endif
#if defined(GF256_HAVE_REFERENCE_BENCHMARKS)
  gf256_benchmarks::RegisterLeopardReferenceBenchmarks();
  gf256_benchmarks::RegisterXDRSReferenceBenchmarks();
#endif
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv)) {
    return 1;
  }
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return 0;
}
