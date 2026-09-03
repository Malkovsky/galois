#include <benchmark/benchmark.h>
#include <galois.h>
#include <jerasure.h>
#include <reed_sol.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../rs_benchmark_cases.h"
#include "field.h"
#include "reference_benchmarks.h"

namespace gf256_benchmarks {
namespace {

using gf2p8::Element;

struct FreeAllocation {
  template <typename T>
  void operator()(T* pointer) const {
    std::free(pointer);
  }
};

using ByteBuffer = std::unique_ptr<Element, FreeAllocation>;
using Matrix = std::unique_ptr<int, FreeAllocation>;

void InitializeJerasure() {
  static const bool initialized = [] {
    if (galois_init_default_field(8) != 0) {
      throw std::runtime_error("failed to initialize Jerasure GF(256)");
    }
    return true;
  }();
  static_cast<void>(initialized);
}

std::vector<ByteBuffer> MakeBuffers(size_t count,
                                    size_t bytes,
                                    bool initialize) {
  std::mt19937 random(42);
  std::vector<ByteBuffer> result;
  result.reserve(count);
  for (size_t buffer = 0; buffer < count; ++buffer) {
    auto* raw = static_cast<Element*>(std::aligned_alloc(64, bytes));
    if (raw == nullptr) {
      throw std::bad_alloc();
    }
    if (initialize) {
      std::generate_n(raw, bytes,
                      [&random] { return static_cast<Element>(random()); });
    } else {
      std::fill_n(raw, bytes, Element{0});
    }
    result.emplace_back(raw);
  }
  return result;
}

std::vector<char*> MutablePointers(std::vector<ByteBuffer>& buffers) {
  std::vector<char*> result;
  result.reserve(buffers.size());
  for (auto& buffer : buffers) {
    result.push_back(reinterpret_cast<char*>(buffer.get()));
  }
  return result;
}

class JerasureCode {
 public:
  JerasureCode(size_t k, size_t r, size_t bytes)
      : k_(static_cast<int>(k)),
        r_(static_cast<int>(r)),
        bytes_(static_cast<int>(bytes)),
        data_(MakeBuffers(k, bytes, true)),
        recovery_(MakeBuffers(r, bytes, false)),
        data_pointers_(MutablePointers(data_)),
        recovery_pointers_(MutablePointers(recovery_)),
        original_data_(k, std::vector<Element>(bytes)),
        erased_(k + r),
        decoding_matrix_(k * k),
        survivor_ids_(k),
        missing_data_(k) {
    InitializeJerasure();
    matrix_.reset(reed_sol_vandermonde_coding_matrix(k_, r_, 8));
    if (matrix_ == nullptr) {
      throw std::runtime_error(
          "failed to construct Jerasure Vandermonde matrix");
    }

    for (size_t i = 0; i < k; ++i) {
      std::copy_n(data_[i].get(), bytes, original_data_[i].begin());
    }
    Encode();

    const auto canonical = MaxErasurePattern(k, r);
    for (size_t i = 0; i < k; ++i) {
      erased_[i] = canonical[r + i] != 0;
      if (erased_[i]) {
        missing_data_[missing_data_count_++] = static_cast<int>(i);
      }
    }
    for (size_t i = 0; i < r; ++i) {
      erased_[k + i] = canonical[i] != 0;
    }
  }

  void Encode() {
    jerasure_matrix_encode(k_, r_, 8, matrix_.get(), data_pointers_.data(),
                           recovery_pointers_.data(), bytes_);
  }

  bool PrepareDecode() {
    return jerasure_make_decoding_matrix(
               k_, r_, 8, matrix_.get(), erased_.data(),
               decoding_matrix_.data(), survivor_ids_.data()) == 0;
  }

  void RepairPrepared() {
    for (size_t output = 0; output < missing_data_count_; ++output) {
      const int data_id = missing_data_[output];
      jerasure_matrix_dotprod(
          k_, 8, decoding_matrix_.data() + data_id * k_, survivor_ids_.data(),
          data_id, data_pointers_.data(), recovery_pointers_.data(), bytes_);
    }
  }

  bool Decode() {
    if (!PrepareDecode()) {
      return false;
    }
    RepairPrepared();
    return true;
  }

  bool VerifyEncode() const {
    for (int output = 0; output < r_; ++output) {
      const int* coefficients = matrix_.get() + output * k_;
      for (int byte = 0; byte < bytes_; ++byte) {
        Element expected = 0;
        for (int input = 0; input < k_; ++input) {
          expected ^=
              gf2p8::MultiplyStandard(static_cast<Element>(coefficients[input]),
                                      data_[input].get()[byte]);
        }
        if (recovery_[output].get()[byte] != expected) {
          return false;
        }
      }
    }
    return true;
  }

  bool VerifyDecode() {
    const size_t erasure_count = std::count(erased_.begin(), erased_.end(), 1);
    if (erasure_count != static_cast<size_t>(r_) || missing_data_count_ == 0) {
      return false;
    }
    for (size_t output = 0; output < missing_data_count_; ++output) {
      std::fill_n(data_[missing_data_[output]].get(), bytes_, Element{0xa5});
    }
    for (int recovery = 0; recovery < r_; ++recovery) {
      if (erased_[static_cast<size_t>(k_ + recovery)] != 0) {
        std::fill_n(recovery_[recovery].get(), bytes_, Element{0x5a});
      }
    }
    if (!Decode()) {
      return false;
    }
    for (size_t output = 0; output < missing_data_count_; ++output) {
      const size_t data_id = static_cast<size_t>(missing_data_[output]);
      if (!std::equal(original_data_[data_id].begin(),
                      original_data_[data_id].end(), data_[data_id].get())) {
        return false;
      }
    }
    return true;
  }

 private:
  int k_;
  int r_;
  int bytes_;
  Matrix matrix_;
  std::vector<ByteBuffer> data_;
  std::vector<ByteBuffer> recovery_;
  std::vector<char*> data_pointers_;
  std::vector<char*> recovery_pointers_;
  std::vector<std::vector<Element>> original_data_;
  std::vector<int> erased_;
  std::vector<int> decoding_matrix_;
  std::vector<int> survivor_ids_;
  std::vector<int> missing_data_;
  size_t missing_data_count_ = 0;
};

void JerasureEncode(benchmark::State& state) {
  const size_t k = static_cast<size_t>(state.range(0));
  const size_t r = static_cast<size_t>(state.range(1));
  const size_t bytes = static_cast<size_t>(state.range(2));
  JerasureCode code(k, r, bytes);
  if (!code.VerifyEncode()) {
    state.SkipWithError("Jerasure encode verification failed");
    return;
  }

  for (auto _ : state) {
    code.Encode();
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(k * bytes));
}

void JerasureDecodeMax(benchmark::State& state) {
  const size_t k = static_cast<size_t>(state.range(0));
  const size_t r = static_cast<size_t>(state.range(1));
  const size_t bytes = static_cast<size_t>(state.range(2));
  JerasureCode code(k, r, bytes);
  if (!code.VerifyEncode() || !code.VerifyDecode()) {
    state.SkipWithError("Jerasure DecodeMax verification failed");
    return;
  }

  for (auto _ : state) {
    bool ok = code.Decode();
    benchmark::DoNotOptimize(ok);
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(
      static_cast<int64_t>(state.iterations()) *
      static_cast<int64_t>(DecodeInputBytesPerIteration(k, bytes)));
}

void JerasureRepairPrepared(benchmark::State& state) {
  const size_t k = static_cast<size_t>(state.range(0));
  const size_t r = static_cast<size_t>(state.range(1));
  const size_t bytes = static_cast<size_t>(state.range(2));
  JerasureCode code(k, r, bytes);
  if (!code.VerifyEncode() || !code.VerifyDecode() || !code.PrepareDecode()) {
    state.SkipWithError("Jerasure prepared repair verification failed");
    return;
  }

  for (auto _ : state) {
    code.RepairPrepared();
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(
      static_cast<int64_t>(state.iterations()) *
      static_cast<int64_t>(DecodeInputBytesPerIteration(k, bytes)));
}

void RegisterTopLevel(const char* prefix) {
  for (const auto& [suffix, function] :
       std::array<std::pair<const char*, void (*)(benchmark::State&)>, 2>{{
           {"Encode", &JerasureEncode},
           {"DecodeMax", &JerasureDecodeMax},
       }}) {
    auto* registered =
        benchmark::RegisterBenchmark(std::string(prefix) + suffix, function);
    for (const auto& test_case : kLCHComparisonCases) {
      registered->Args(
          {test_case.data_count, test_case.recovery_count, test_case.bytes});
    }
    registered->ArgNames({"K", "R", "bytes"});
  }
}

}  // namespace

void RegisterJerasureReferenceBenchmarks() {
#if defined(GF256_HIGH_LEVEL_BENCHMARKS_ONLY)
  RegisterTopLevel("RS/Jerasure/Native/");
#else
  RegisterTopLevel("Jerasure/Native/");
  auto* repair = benchmark::RegisterBenchmark("Jerasure/Native/RepairPrepared",
                                              &JerasureRepairPrepared);
  for (const auto& test_case : kLCHComparisonCases) {
    repair->Args(
        {test_case.data_count, test_case.recovery_count, test_case.bytes});
  }
  repair->ArgNames({"K", "R", "bytes"});
#endif
}

}  // namespace gf256_benchmarks
