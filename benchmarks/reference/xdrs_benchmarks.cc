#include <P_function.h>
#include <benchmark/benchmark.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "../rs_benchmark_cases.h"
#include "reference_benchmarks.h"

namespace gf256_benchmarks {
namespace {

struct XDRSTables {
  std::array<GFSymbol, Size> index_of{};
  std::array<GFSymbol, Size> alpha_to{};
  // Upstream init_dec() accesses the sentinel entry at index mod.
  std::array<GFSymbol, Size> skew{};
  std::array<GFSymbol, Size> derivative{};
  std::array<GFSymbol, Size> log_walsh{};
  std::array<GFSymbol, len - 1> base{};
  std::array<GFSymbol, Size> coefficient_low{};
  std::array<GFSymbol, Size> coefficient_high{};
  std::array<int, len> subspace{};
};

std::vector<std::vector<GFSymbol>> MakeBuffers(size_t count, size_t bytes) {
  std::mt19937 random(42);
  std::vector<std::vector<GFSymbol>> result(count,
                                            std::vector<GFSymbol>(bytes));
  for (auto& buffer : result) {
    std::generate(buffer.begin(), buffer.end(),
                  [&random] { return static_cast<GFSymbol>(random()); });
  }
  return result;
}

std::vector<GFSymbol*> MutablePointers(
    std::vector<std::vector<GFSymbol>>& buffers) {
  std::vector<GFSymbol*> result;
  result.reserve(buffers.size());
  for (auto& buffer : buffers) {
    result.push_back(buffer.data());
  }
  return result;
}

std::vector<const GFSymbol*> ConstPointers(
    const std::vector<std::vector<GFSymbol>>& buffers) {
  std::vector<const GFSymbol*> result;
  result.reserve(buffers.size());
  for (const auto& buffer : buffers) {
    result.push_back(buffer.data());
  }
  return result;
}

void InitializeXDRS(XDRSTables& tables, unsigned bytes, unsigned k) {
  ::function::main_params(bytes, 1, k);
  ::function::array_params(
      tables.index_of.data(), tables.alpha_to.data(), tables.skew.data(),
      tables.derivative.data(), tables.log_walsh.data(), tables.base.data(),
      tables.coefficient_low.data(), tables.coefficient_high.data(),
      tables.subspace.data());
  ::function::init();
  ::function::init_dec();
}

bool VerifyXDRSRecovery(
    unsigned k, unsigned bytes, bool low_rate,
    const std::vector<std::vector<GFSymbol>>& data,
    const std::vector<std::vector<GFSymbol>>& recovery) {
  const unsigned recovery_count = Size - k;
  auto data_pointers = ConstPointers(data);
  auto recovery_pointers = ConstPointers(recovery);
  auto codeword = MakeBuffers(Size, bytes);
  auto codeword_pointers = MutablePointers(codeword);
  std::unique_ptr<bool[]> erasure(new bool[Size]{});
  const auto erased = XDRSMaxErasurePattern(k, recovery_count);
  size_t missing_data_count = 0;
  for (size_t i = 0; i < erased.size(); ++i) {
    erasure[i] = erased[i] != 0;
  }
  for (size_t i = 0; i < k; ++i) {
    const size_t codeword_index = low_rate ? i : recovery_count + i;
    missing_data_count += erased[codeword_index] != 0;
  }
  std::array<GFSymbol, Size> locator{};
  if (low_rate) {
    ::function::Algorithm2::ReedSolomondecodeL(
        data_pointers.data(), recovery_pointers.data(),
        codeword_pointers.data(), erasure.get(), locator.data());
  } else {
    ::function::Algorithm3::ReedSolomondecodeH(
        data_pointers.data(), recovery_pointers.data(),
        codeword_pointers.data(), erasure.get(), locator.data());
  }
  if (std::count(erased.begin(), erased.end(), uint8_t{1}) != recovery_count ||
      missing_data_count == 0) {
    return false;
  }
  for (size_t i = 0; i < k; ++i) {
    const size_t codeword_index = low_rate ? i : recovery_count + i;
    if (erased[codeword_index] && codeword[codeword_index] != data[i]) {
      return false;
    }
  }
  return true;
}

void NativeXDRSEncode(benchmark::State& state,
                      bool low_rate,
                      size_t bytes_range = 1) {
  const unsigned k = static_cast<unsigned>(state.range(0));
  const unsigned bytes = static_cast<unsigned>(state.range(bytes_range));
  const unsigned recovery_count = Size - k;
  const unsigned parity_workspace_count =
      low_rate ? recovery_count : 2 * recovery_count;
  XDRSTables tables;
  InitializeXDRS(tables, bytes, k);
  auto data = MakeBuffers(k, bytes);
  auto recovery = MakeBuffers(parity_workspace_count, bytes);
  auto data_pointers = ConstPointers(data);
  auto recovery_pointers = MutablePointers(recovery);
  const auto encode = [&] {
    if (low_rate) {
      ::function::ReedSolomonEncodeL(data_pointers.data(),
                                     recovery_pointers.data());
    } else {
      ::function::ReedSolomonEncodeH(data_pointers.data(),
                                     recovery_pointers.data());
    }
  };
  encode();
  if (!VerifyXDRSRecovery(k, bytes, low_rate, data, recovery)) {
    state.SkipWithError("native XDRS encode verification failed");
    return;
  }

  for (auto _ : state) {
    encode();
    benchmark::ClobberMemory();
  }
  if (!VerifyXDRSRecovery(k, bytes, low_rate, data, recovery)) {
    state.SkipWithError("native XDRS encode verification failed");
    return;
  }
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * k * bytes);
}

void NativeXDRSDecode(benchmark::State& state,
                      bool low_rate,
                      size_t bytes_range = 1) {
  const unsigned k = static_cast<unsigned>(state.range(0));
  const unsigned bytes = static_cast<unsigned>(state.range(bytes_range));
  const unsigned recovery_count = Size - k;
  const unsigned parity_workspace_count =
      low_rate ? recovery_count : 2 * recovery_count;
  XDRSTables tables;
  InitializeXDRS(tables, bytes, k);
  auto data = MakeBuffers(k, bytes);
  auto recovery = MakeBuffers(parity_workspace_count, bytes);
  auto data_pointers = ConstPointers(data);
  auto recovery_pointers = MutablePointers(recovery);
  if (low_rate) {
    ::function::ReedSolomonEncodeL(data_pointers.data(),
                                   recovery_pointers.data());
  } else {
    ::function::ReedSolomonEncodeH(data_pointers.data(),
                                   recovery_pointers.data());
  }
  auto recovery_const = ConstPointers(recovery);
  auto codeword = MakeBuffers(Size, bytes);
  auto codeword_pointers = MutablePointers(codeword);
  std::unique_ptr<bool[]> erasure(new bool[Size]{});
  const auto erased = XDRSMaxErasurePattern(k, recovery_count);
  size_t missing_data_count = 0;
  for (size_t i = 0; i < erased.size(); ++i) {
    erasure[i] = erased[i] != 0;
  }
  for (size_t i = 0; i < k; ++i) {
    const size_t codeword_index = low_rate ? i : recovery_count + i;
    missing_data_count += erased[codeword_index] != 0;
  }
  std::array<GFSymbol, Size> locator{};

  const auto decode = [&] {
    if (low_rate) {
      ::function::Algorithm2::ReedSolomondecodeL(
          data_pointers.data(), recovery_const.data(), codeword_pointers.data(),
          erasure.get(), locator.data());
    } else {
      ::function::Algorithm3::ReedSolomondecodeH(
          data_pointers.data(), recovery_const.data(), codeword_pointers.data(),
          erasure.get(), locator.data());
    }
  };
  decode();
  if (std::count(erased.begin(), erased.end(), uint8_t{1}) != recovery_count ||
      missing_data_count == 0) {
    state.SkipWithError("native XDRS decode verification failed");
    return;
  }
  for (size_t i = 0; i < k; ++i) {
    const size_t codeword_index = low_rate ? i : recovery_count + i;
    if (erased[codeword_index] && codeword[codeword_index] != data[i]) {
      state.SkipWithError("native XDRS decode verification failed");
      return;
    }
  }

  for (auto _ : state) {
    decode();
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(
      static_cast<int64_t>(state.iterations()) *
      static_cast<int64_t>(DecodeInputBytesPerIteration(k, bytes)));
}

void RegisterRate(const char* rate_name,
                  bool low_rate,
                  std::initializer_list<int64_t> data_counts) {
  auto* encode = benchmark::RegisterBenchmark(
      std::string("XDRS/Native/") + rate_name + "/Encode",
      [low_rate](benchmark::State& state) {
        NativeXDRSEncode(state, low_rate);
      });
  auto* decode = benchmark::RegisterBenchmark(
      std::string("XDRS/Native/") + rate_name + "/DecodeMax",
      [low_rate](benchmark::State& state) {
        NativeXDRSDecode(state, low_rate);
      });
  for (const int64_t k : data_counts) {
    encode->Args({k, 1024});
    decode->Args({k, 1024});
  }
  encode->ArgNames({"K", "bytes"});
  decode->ArgNames({"K", "bytes"});
}

}  // namespace

void RegisterXDRSReferenceBenchmarks() {
#if defined(GF256_HIGH_LEVEL_BENCHMARKS_ONLY)
  const auto register_operation = [](const char* name, bool decode) {
    auto* registered =
        benchmark::RegisterBenchmark(name, [decode](benchmark::State& state) {
          const bool low_rate = NativeXDRSLowRate(state.range(0));
          if (decode) {
            NativeXDRSDecode(state, low_rate, 2);
          } else {
            NativeXDRSEncode(state, low_rate, 2);
          }
        });
    for (const auto& test_case : kLCHComparisonCases) {
      registered->Args(
          {test_case.data_count, test_case.recovery_count, test_case.bytes});
    }
    registered->ArgNames({"K", "R", "bytes"});
  };
  register_operation("RS/XDRS/Native/Encode", false);
  register_operation("RS/XDRS/Native/DecodeMax", true);
#else
  RegisterRate("Low", true, {8, 16, 32, 64, 128});
  RegisterRate("High", false, {128, 192, 224, 240, 248});
#endif
}

}  // namespace gf256_benchmarks
