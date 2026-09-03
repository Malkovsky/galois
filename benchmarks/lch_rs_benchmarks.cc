#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "lin_chung_han/transform.h"
#include "reed_solomon/lch_decoder.h"
#include "reed_solomon/lch_encoder.h"
#include "rs_benchmark_cases.h"

#if defined(GF256_ENABLE_GFNI512_RADIX8_EXPERIMENT)
#include "lin_chung_han/experiment/gfni512_radix8.h"
#include "reed_solomon/experiment/gfni512_radix8.h"
#endif

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

void LCHEncodeBenchmark(benchmark::State& state, Backend backend, Radix radix) {
  if (!gf2p8::lch::BackendAvailable(backend)) {
    state.SkipWithError("backend was not compiled");
    return;
  }
  const size_t data_count = static_cast<size_t>(state.range(0));
  const size_t recovery_count = static_cast<size_t>(state.range(1));
  const size_t bytes = static_cast<size_t>(state.range(2));
  const gf2p8::rs::LCHEncoder encoder(data_count, recovery_count);
  auto data = MakeShards(data_count, bytes);
  auto recovery = MakeShards(recovery_count, bytes);
  const auto data_pointers = ConstPointers(data);
  auto recovery_pointers = MutablePointers(recovery);
  std::vector<Element> workspace(encoder.WorkspaceSize(bytes));
  auto expected_recovery = MakeShards(recovery_count, bytes);
  auto expected_pointers = MutablePointers(expected_recovery);
  std::vector<Element> expected_workspace(encoder.WorkspaceSize(bytes));
  if (encoder.Encode(data_pointers, expected_pointers, bytes,
                     expected_workspace, Backend::scalar,
                     Radix::radix2) != Status::ok ||
      encoder.Encode(data_pointers, recovery_pointers, bytes, workspace,
                     backend, radix) != Status::ok ||
      recovery != expected_recovery) {
    state.SkipWithError("LCH encode verification failed");
    return;
  }

  for (auto _ : state) {
    Status status = encoder.Encode(data_pointers, recovery_pointers, bytes,
                                   workspace, backend, radix);
    benchmark::DoNotOptimize(status);
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(data_count) *
                          static_cast<int64_t>(bytes));
}

#if defined(GF256_ENABLE_GFNI512_RADIX8_EXPERIMENT)

void Radix8Control(std::vector<std::vector<Element>>& shards,
                   size_t bytes,
                   const gf2p8::lch::detail::ResolvedKernels& kernels,
                   const Context& context) {
  auto pointers = MutablePointers(shards);
  constexpr size_t kOffset = 8;
  const Element top = context.Skew(2, kOffset);
  const Element middle0 = context.Skew(1, kOffset);
  const Element middle1 = context.Skew(1, kOffset ^ 4);
  const Element low0 = context.Skew(0, kOffset);
  const Element low1 = context.Skew(0, kOffset ^ 2);
  const Element low2 = context.Skew(0, kOffset ^ 4);
  const Element low3 = context.Skew(0, kOffset ^ 6);
  kernels.ifft_radix4(pointers[0], pointers[1], pointers[2], pointers[3], bytes,
                      middle0, low0, low1, context.Tables());
  kernels.ifft_radix4(pointers[4], pointers[5], pointers[6], pointers[7], bytes,
                      middle1, low2, low3, context.Tables());
  for (size_t i = 0; i < 4; ++i) {
    kernels.ifft_radix2(pointers[i], pointers[i + 4], bytes, top,
                        context.Tables());
  }
}

void Radix8Candidate(
    std::vector<std::vector<Element>>& shards,
    size_t bytes,
    const gf2p8::lch::detail::experiment::radix8::Kernels& kernels,
    const Context& context) {
  auto pointers = MutablePointers(shards);
  constexpr size_t kOffset = 8;
  kernels.ifft_radix8(pointers[0], pointers[1], pointers[2], pointers[3],
                      pointers[4], pointers[5], pointers[6], pointers[7], bytes,
                      context.Skew(2, kOffset), context.Skew(1, kOffset),
                      context.Skew(1, kOffset ^ 4), context.Skew(0, kOffset),
                      context.Skew(0, kOffset ^ 2),
                      context.Skew(0, kOffset ^ 4),
                      context.Skew(0, kOffset ^ 6), context.Tables());
}

void Radix8LeafBenchmark(benchmark::State& state, bool experiment) {
  const auto* radix8 = gf2p8::lch::detail::experiment::radix8::ResolveKernels();
  const size_t bytes = static_cast<size_t>(state.range(0));
  const auto* base =
      gf2p8::lch::detail::ResolveKernels(Backend::gfni512_affine, bytes);
  if (radix8 == nullptr || base == nullptr) {
    state.SkipWithError("GFNI512 radix-8 experiment was not compiled");
    return;
  }
  const Context& context = Context::Shared();
  const auto input = MakeShards(8, bytes);
  auto expected = input;
  auto expected_pointers = MutablePointers(expected);
  if (gf2p8::lch::IFFT(context, expected_pointers, bytes, 8,
                       {Backend::scalar, Radix::radix2}) != Status::ok) {
    state.SkipWithError("scalar leaf reference failed");
    return;
  }
  auto shards = input;
  if (experiment) {
    Radix8Candidate(shards, bytes, *radix8, context);
  } else {
    Radix8Control(shards, bytes, *base, context);
  }
  if (shards != expected) {
    state.SkipWithError("radix-8 leaf verification failed");
    return;
  }

  for (auto _ : state) {
    if (experiment) {
      Radix8Candidate(shards, bytes, *radix8, context);
    } else {
      Radix8Control(shards, bytes, *base, context);
    }
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * 8 *
                          static_cast<int64_t>(bytes));
}

void Radix8IFFTBenchmark(benchmark::State& state, bool experiment) {
  const auto* radix8 = gf2p8::lch::detail::experiment::radix8::ResolveKernels();
  const size_t transform_size = static_cast<size_t>(state.range(0));
  const size_t bytes = static_cast<size_t>(state.range(1));
  const auto* base =
      gf2p8::lch::detail::ResolveKernels(Backend::gfni512_affine, bytes);
  if (radix8 == nullptr || base == nullptr) {
    state.SkipWithError("GFNI512 radix-8 experiment was not compiled");
    return;
  }
  const Context& context = Context::Shared();
  const auto input = MakeShards(transform_size, bytes);
  auto expected = input;
  auto expected_pointers = MutablePointers(expected);
  if (gf2p8::lch::IFFT(context, expected_pointers, bytes, 0,
                       {Backend::scalar, Radix::radix2}) != Status::ok) {
    state.SkipWithError("scalar IFFT reference failed");
    return;
  }
  auto shards = input;
  auto pointers = MutablePointers(shards);
  const Status verification =
      experiment
          ? gf2p8::lch::detail::experiment::radix8::IFFTResolved(
                context, pointers, bytes, 0, transform_size, *base, *radix8)
          : gf2p8::lch::detail::IFFTResolved(context, pointers, bytes, 0,
                                             transform_size, *base,
                                             Radix::radix4);
  if (verification != Status::ok || shards != expected) {
    state.SkipWithError("radix-8 IFFT verification failed");
    return;
  }

  for (auto _ : state) {
    Status status =
        experiment
            ? gf2p8::lch::detail::experiment::radix8::IFFTResolved(
                  context, pointers, bytes, 0, transform_size, *base, *radix8)
            : gf2p8::lch::detail::IFFTResolved(context, pointers, bytes, 0,
                                               transform_size, *base,
                                               Radix::radix4);
    benchmark::DoNotOptimize(status);
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(transform_size) *
                          static_cast<int64_t>(bytes));
}

size_t NextPowerOfTwo(size_t value) {
  size_t result = 1;
  while (result < value) {
    result *= 2;
  }
  return result;
}

void Radix8EncodeBenchmark(benchmark::State& state, bool experiment) {
  const auto* radix8 = gf2p8::lch::detail::experiment::radix8::ResolveKernels();
  const size_t data_count = static_cast<size_t>(state.range(0));
  const size_t recovery_count = static_cast<size_t>(state.range(1));
  const size_t bytes = static_cast<size_t>(state.range(2));
  const auto* base =
      gf2p8::lch::detail::ResolveKernels(Backend::gfni512_affine, bytes);
  if (radix8 == nullptr || base == nullptr) {
    state.SkipWithError("GFNI512 radix-8 experiment was not compiled");
    return;
  }
  const gf2p8::rs::LCHEncoder encoder(data_count, recovery_count);
  const size_t transform_size = NextPowerOfTwo(recovery_count);
  const auto data = MakeShards(data_count, bytes);
  const auto data_pointers = ConstPointers(data);
  auto expected = MakeShards(recovery_count, bytes);
  auto expected_pointers = MutablePointers(expected);
  std::vector<Element> expected_workspace(encoder.WorkspaceSize(bytes));
  if (encoder.Encode(data_pointers, expected_pointers, bytes,
                     expected_workspace, Backend::scalar,
                     Radix::radix2) != Status::ok) {
    state.SkipWithError("scalar encode reference failed");
    return;
  }
  auto recovery = MakeShards(recovery_count, bytes);
  auto recovery_pointers = MutablePointers(recovery);
  std::vector<Element> workspace(encoder.WorkspaceSize(bytes));
  const auto run = [&] {
    return experiment ? gf2p8::rs::detail::experiment::radix8::EncodeLCH(
                            Context::Shared(), data_pointers, recovery_pointers,
                            bytes, transform_size, workspace, *base, *radix8)
                      : encoder.Encode(data_pointers, recovery_pointers, bytes,
                                       workspace, Backend::gfni512_affine,
                                       Radix::radix4);
  };
  if (run() != Status::ok || recovery != expected) {
    state.SkipWithError("radix-8 encode verification failed");
    return;
  }

  for (auto _ : state) {
    Status status = run();
    benchmark::DoNotOptimize(status);
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(data_count) *
                          static_cast<int64_t>(bytes));
}

#endif

void LCHDecodeBenchmark(benchmark::State& state, Backend backend, Radix radix) {
  if (!gf2p8::lch::BackendAvailable(backend)) {
    state.SkipWithError("backend was not compiled");
    return;
  }
  const size_t data_count = static_cast<size_t>(state.range(0));
  const size_t recovery_count = static_cast<size_t>(state.range(1));
  const size_t bytes = static_cast<size_t>(state.range(2));
  const gf2p8::rs::LCHEncoder encoder(data_count, recovery_count);
  const gf2p8::rs::LCHDecoder decoder(data_count, recovery_count);
  auto data = MakeShards(data_count, bytes);
  auto recovery = MakeShards(recovery_count, bytes);
  const auto data_input = ConstPointers(data);
  auto recovery_pointers = MutablePointers(recovery);
  std::vector<Element> encode_workspace(encoder.WorkspaceSize(bytes));
  if (encoder.Encode(data_input, recovery_pointers, bytes, encode_workspace,
                     backend, radix) != Status::ok) {
    state.SkipWithError("encode setup failed");
    return;
  }

  auto data_pointers = MutablePointers(data);
  auto recovery_input = ConstPointers(recovery);
  const auto erased =
      gf256_benchmarks::MaxErasurePattern(data_count, recovery_count);
  std::vector<uint8_t> data_present(data_count);
  std::vector<uint8_t> recovery_present(recovery_count);
  for (size_t i = 0; i < recovery_count; ++i) {
    recovery_present[i] = static_cast<uint8_t>(!erased[i]);
    if (erased[i] != 0) {
      recovery_input[i] = nullptr;
    }
  }
  const auto expected_data = data;
  size_t missing_data_count = 0;
  for (size_t i = 0; i < data_count; ++i) {
    data_present[i] = static_cast<uint8_t>(!erased[recovery_count + i]);
    if (erased[recovery_count + i]) {
      std::fill(data[i].begin(), data[i].end(), 0);
      ++missing_data_count;
    }
  }
  std::vector<Element> workspace(decoder.WorkspaceSize(bytes));
  if (decoder.Decode(data_pointers, data_present, recovery_input,
                     recovery_present, bytes, workspace, backend,
                     radix) != Status::ok ||
      std::count(erased.begin(), erased.end(), uint8_t{1}) != recovery_count ||
      missing_data_count == 0) {
    state.SkipWithError("LCH decode verification failed");
    return;
  }
  for (size_t i = 0; i < data_count; ++i) {
    if (erased[recovery_count + i] && data[i] != expected_data[i]) {
      state.SkipWithError("LCH decode verification failed");
      return;
    }
  }

  for (auto _ : state) {
    Status status =
        decoder.Decode(data_pointers, data_present, recovery_input,
                       recovery_present, bytes, workspace, backend, radix);
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

#if defined(GF256_ENABLE_GFNI512_RADIX8_EXPERIMENT)
void RegisterRadix8ExperimentBenchmarks() {
  const int64_t leaf_lengths[] = {1,   8,   63,   64,   65,
                                  128, 256, 1024, 4096, 65536};
  for (const bool experiment : {false, true}) {
    auto* leaf = benchmark::RegisterBenchmark(
        std::string("LCH/IFFTRadix8Leaf/GFNI512Affine/") +
            (experiment ? "Radix8Experiment" : "Radix4Control"),
        [experiment](benchmark::State& state) {
          Radix8LeafBenchmark(state, experiment);
        });
    for (const int64_t bytes : leaf_lengths) {
      leaf->Arg(bytes);
    }
    leaf->ArgName("bytes");

    auto* ifft = benchmark::RegisterBenchmark(
        std::string("LCH/IFFT/GFNI512Affine/") +
            (experiment ? "Radix8Experiment" : "Radix4Control"),
        [experiment](benchmark::State& state) {
          Radix8IFFTBenchmark(state, experiment);
        });
    for (int64_t transform = 8; transform <= 256; transform *= 2) {
      for (const int64_t bytes : {1, 63, 64, 65, 256, 1024, 4096, 65536}) {
        ifft->Args({transform, bytes});
      }
    }
    ifft->ArgNames({"transform", "bytes"});

    auto* encode = benchmark::RegisterBenchmark(
        std::string("LCH/Experiment/Encode/GFNI512Affine/") +
            (experiment ? "Radix8" : "Radix4Control"),
        [experiment](benchmark::State& state) {
          Radix8EncodeBenchmark(state, experiment);
        });
    encode->Args({32, 16, 1024})
        ->Args({224, 32, 1024})
        ->Args({129, 64, 1024})
        ->Args({128, 128, 1024})
        ->ArgNames({"K", "R", "bytes"});
  }
}
#endif

void RegisterRsBenchmarks() {
  const Backend backends[] = {Backend::scalar,         Backend::ssse3,
                              Backend::avx2,           Backend::gfni128_affine,
                              Backend::gfni256_affine, Backend::gfni512_affine};
  for (const Backend backend : backends) {
    for (const Radix radix : {Radix::radix2, Radix::radix4}) {
      const char* const radix_name =
          radix == Radix::radix2 ? "Radix2" : "Radix4";
      auto* encode = benchmark::RegisterBenchmark(
          std::string("LCH/Owned/Encode/") + BackendName(backend) + "/" +
              radix_name,
          [backend, radix](benchmark::State& state) {
            LCHEncodeBenchmark(state, backend, radix);
          });
      auto* decode = benchmark::RegisterBenchmark(
          std::string("LCH/Owned/DecodeMax/") + BackendName(backend) + "/" +
              radix_name,
          [backend, radix](benchmark::State& state) {
            LCHDecodeBenchmark(state, backend, radix);
          });
      for (const auto& test_case : gf256_benchmarks::kLCHComparisonCases) {
        encode->Args(
            {test_case.data_count, test_case.recovery_count, test_case.bytes});
        decode->Args(
            {test_case.data_count, test_case.recovery_count, test_case.bytes});
      }
      encode->Args({5, 3, 1024})->Args({5, 6, 1024});
      decode->Args({5, 3, 1024})->Args({5, 6, 1024});
      encode->ArgNames({"K", "R", "bytes"});
      decode->ArgNames({"K", "R", "bytes"});
    }
  }
}

void RegisterHighLevelRsBenchmarks() {
  constexpr Radix radix = Radix::radix4;

  const auto register_lch = [](const char* name, Backend selected_backend,
                               bool decode) {
    auto* registered = benchmark::RegisterBenchmark(
        name, [selected_backend, decode](benchmark::State& state) {
          if (decode) {
            LCHDecodeBenchmark(state, selected_backend, radix);
          } else {
            LCHEncodeBenchmark(state, selected_backend, radix);
          }
        });
    for (const auto& test_case : gf256_benchmarks::kLCHComparisonCases) {
      registered->Args(
          {test_case.data_count, test_case.recovery_count, test_case.bytes});
    }
    registered->ArgNames({"K", "R", "bytes"});
  };

  register_lch("RS/LCH/Owned/AVX2/Encode", Backend::avx2, false);
  register_lch("RS/LCH/Owned/AVX2/DecodeMax", Backend::avx2, true);
  register_lch("RS/LCH/Owned/GFNI256/Encode", Backend::gfni256_affine, false);
  register_lch("RS/LCH/Owned/GFNI256/DecodeMax", Backend::gfni256_affine, true);
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
#if defined(GF256_ENABLE_GFNI512_RADIX8_EXPERIMENT)
  RegisterRadix8ExperimentBenchmarks();
#endif
#endif
#if defined(GF256_HAVE_REFERENCE_BENCHMARKS)
  gf256_benchmarks::RegisterISALReferenceBenchmarks();
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
