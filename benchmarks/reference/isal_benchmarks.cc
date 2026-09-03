#include <benchmark/benchmark.h>
#define ISAL_DEPRECATED_INTERNAL
#include <isa-l/erasure_code.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "../rs_benchmark_cases.h"
#include "field.h"
#include "lin_chung_han/kernels.h"
#include "lin_chung_han/kernels_internal.h"
#include "reference_benchmarks.h"

namespace gf256_benchmarks {
namespace {

using gf2p8::Element;
namespace lch = gf2p8::lch;

std::vector<std::vector<Element>> MakeBuffers(size_t count, size_t bytes) {
  std::mt19937 random(42);
  std::vector<std::vector<Element>> result(count, std::vector<Element>(bytes));
  for (auto& buffer : result) {
    std::generate(buffer.begin(), buffer.end(),
                  [&random] { return static_cast<Element>(random()); });
  }
  return result;
}

std::vector<Element*> MutablePointers(
    std::vector<std::vector<Element>>& buffers) {
  std::vector<Element*> result;
  result.reserve(buffers.size());
  for (auto& buffer : buffers) {
    result.push_back(buffer.data());
  }
  return result;
}

bool VerifyParity(size_t k,
                  size_t r,
                  size_t bytes,
                  const std::vector<Element>& generator,
                  const std::vector<std::vector<Element>>& data,
                  const std::vector<std::vector<Element>>& recovery) {
  for (size_t output = 0; output < r; ++output) {
    const Element* coefficients = generator.data() + (k + output) * k;
    for (size_t byte = 0; byte < bytes; ++byte) {
      Element expected = 0;
      for (size_t input = 0; input < k; ++input) {
        expected ^=
            gf2p8::MultiplyStandard(coefficients[input], data[input][byte]);
      }
      if (recovery[output][byte] != expected) {
        return false;
      }
    }
  }
  return true;
}

class ISALCode {
 public:
  ISALCode(size_t k, size_t r, size_t bytes)
      : k_(k),
        r_(r),
        bytes_(bytes),
        generator_((k + r) * k),
        encode_tables_(32 * k * r),
        data_(MakeBuffers(k, bytes)),
        recovery_(MakeBuffers(r, bytes)),
        data_pointers_(MutablePointers(data_)),
        recovery_pointers_(MutablePointers(recovery_)),
        erased_(k + r),
        decode_sources_(k),
        decode_outputs_(k),
        decode_output_pointers_(k),
        matrix_work_(k * k),
        matrix_inverse_(k * k),
        decode_matrix_(k * k),
        decode_tables_(32 * k * k),
        missing_data_indices_(k) {
    gf_gen_cauchy1_matrix(generator_.data(), static_cast<int>(k + r),
                          static_cast<int>(k));
    ec_init_tables(static_cast<int>(k), static_cast<int>(r),
                   generator_.data() + k * k, encode_tables_.data());
    Encode();

    const auto canonical = MaxErasurePattern(k, r);
    for (size_t i = 0; i < k; ++i) {
      erased_[i] = canonical[r + i];
    }
    for (size_t i = 0; i < r; ++i) {
      erased_[k + i] = canonical[i];
    }
    for (size_t i = 0; i < decode_outputs_.size(); ++i) {
      decode_outputs_[i].resize(bytes);
      decode_output_pointers_[i] = decode_outputs_[i].data();
    }
  }

  void Encode() {
    ec_encode_data(static_cast<int>(bytes_), static_cast<int>(k_),
                   static_cast<int>(r_), encode_tables_.data(),
                   data_pointers_.data(), recovery_pointers_.data());
  }

  bool PrepareDecode() {
    size_t source_count = 0;
    missing_data_count_ = 0;
    for (size_t row = 0; row < k_ + r_; ++row) {
      if (erased_[row]) {
        if (row < k_) {
          missing_data_indices_[missing_data_count_++] = row;
        }
        continue;
      }

      std::copy_n(generator_.data() + row * k_, k_,
                  matrix_work_.data() + source_count * k_);
      decode_sources_[source_count] =
          row < k_ ? data_[row].data() : recovery_[row - k_].data();
      ++source_count;
    }
    if (source_count != k_ || missing_data_count_ == 0 ||
        gf_invert_matrix(matrix_work_.data(), matrix_inverse_.data(),
                         static_cast<int>(k_)) != 0) {
      return false;
    }

    for (size_t output = 0; output < missing_data_count_; ++output) {
      std::copy_n(matrix_inverse_.data() + missing_data_indices_[output] * k_,
                  k_, decode_matrix_.data() + output * k_);
    }
    ec_init_tables(static_cast<int>(k_), static_cast<int>(missing_data_count_),
                   decode_matrix_.data(), decode_tables_.data());
    return true;
  }

  void RepairPrepared() {
    ec_encode_data(static_cast<int>(bytes_), static_cast<int>(k_),
                   static_cast<int>(missing_data_count_), decode_tables_.data(),
                   decode_sources_.data(), decode_output_pointers_.data());
  }

  bool Decode() {
    if (!PrepareDecode()) {
      return false;
    }
    RepairPrepared();
    return true;
  }

  bool VerifyEncode() const {
    return VerifyParity(k_, r_, bytes_, generator_, data_, recovery_);
  }

  bool VerifyDecode() {
    if (std::count(erased_.begin(), erased_.end(), uint8_t{1}) != r_ ||
        !Decode()) {
      return false;
    }
    size_t output = 0;
    for (size_t i = 0; i < k_; ++i) {
      if (erased_[i] && decode_outputs_[output++] != data_[i]) {
        return false;
      }
    }
    return output != 0;
  }

 private:
  size_t k_;
  size_t r_;
  size_t bytes_;
  std::vector<Element> generator_;
  std::vector<Element> encode_tables_;
  std::vector<std::vector<Element>> data_;
  std::vector<std::vector<Element>> recovery_;
  std::vector<Element*> data_pointers_;
  std::vector<Element*> recovery_pointers_;
  std::vector<uint8_t> erased_;
  std::vector<Element*> decode_sources_;
  std::vector<std::vector<Element>> decode_outputs_;
  std::vector<Element*> decode_output_pointers_;
  std::vector<Element> matrix_work_;
  std::vector<Element> matrix_inverse_;
  std::vector<Element> decode_matrix_;
  std::vector<Element> decode_tables_;
  std::vector<size_t> missing_data_indices_;
  size_t missing_data_count_ = 0;
};

void ISALEncode(benchmark::State& state) {
  const size_t k = static_cast<size_t>(state.range(0));
  const size_t r = static_cast<size_t>(state.range(1));
  const size_t bytes = static_cast<size_t>(state.range(2));
  ISALCode code(k, r, bytes);
  if (!code.VerifyEncode()) {
    state.SkipWithError("ISA-L encode verification failed");
    return;
  }

  for (auto _ : state) {
    code.Encode();
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * k *
                          static_cast<int64_t>(bytes));
}

void ISALDecodeMax(benchmark::State& state) {
  const size_t k = static_cast<size_t>(state.range(0));
  const size_t r = static_cast<size_t>(state.range(1));
  const size_t bytes = static_cast<size_t>(state.range(2));
  ISALCode code(k, r, bytes);
  if (!code.VerifyEncode() || !code.VerifyDecode()) {
    state.SkipWithError("ISA-L DecodeMax verification failed");
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

void ISALRepairPrepared(benchmark::State& state) {
  const size_t k = static_cast<size_t>(state.range(0));
  const size_t r = static_cast<size_t>(state.range(1));
  const size_t bytes = static_cast<size_t>(state.range(2));
  ISALCode code(k, r, bytes);
  if (!code.VerifyEncode() || !code.VerifyDecode() || !code.PrepareDecode()) {
    state.SkipWithError("ISA-L prepared repair verification failed");
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

std::array<Element, 32> ISALCantorTable(Element coefficient) {
  std::array<Element, 32> result{};
  const auto& row = gf2p8::Tables().shuffle[coefficient];
  std::copy_n(row.data(), 16, result.data());
  std::copy_n(row.data() + 32, 16, result.data() + 16);
  return result;
}

enum class PrimitiveOperation { scale, add_scaled };
enum class PrimitiveImplementation { owned, isal };

void PrimitiveBenchmark(benchmark::State& state,
                        PrimitiveOperation operation,
                        PrimitiveImplementation implementation) {
  constexpr Element kCoefficient = 0x53;
  const size_t bytes = static_cast<size_t>(state.range(0));
  auto buffers = MakeBuffers(3, bytes);
  auto expected = buffers[1];
  auto candidate = buffers[1];
  auto source_pointers = MutablePointers(buffers);
  const auto isal_table = ISALCantorTable(kCoefficient);
  const auto& tables = gf2p8::Tables();
  const lch::detail::ResolvedKernels* kernels =
      lch::detail::ResolveKernels(lch::Backend::avx2, bytes);
  const lch::Status expected_status =
      operation == PrimitiveOperation::scale
          ? lch::Scale(expected.data(), buffers[0].data(), bytes, kCoefficient,
                       lch::Backend::avx2, tables)
          : lch::AddScaled(expected.data(), buffers[0].data(), bytes,
                           kCoefficient, lch::Backend::avx2, tables);
  if (implementation == PrimitiveImplementation::owned) {
    candidate = buffers[1];
    if (operation == PrimitiveOperation::scale) {
      kernels->scale(candidate.data(), buffers[0].data(), bytes, kCoefficient,
                     tables);
    } else {
      kernels->add_scaled(candidate.data(), buffers[0].data(), bytes,
                          kCoefficient, tables);
    }
  } else if (operation == PrimitiveOperation::scale) {
    gf_vect_dot_prod_avx2(static_cast<int>(bytes), 1,
                          const_cast<Element*>(isal_table.data()),
                          source_pointers.data(), candidate.data());
  } else {
    gf_vect_mad_avx2(static_cast<int>(bytes), 1, 0,
                     const_cast<Element*>(isal_table.data()), buffers[0].data(),
                     candidate.data());
  }
  if (expected_status != lch::Status::ok || candidate != expected) {
    state.SkipWithError("ISA-L Cantor primitive verification failed");
    return;
  }

  candidate = buffers[1];
  for (auto _ : state) {
    if (implementation == PrimitiveImplementation::owned) {
      if (operation == PrimitiveOperation::scale) {
        kernels->scale(candidate.data(), buffers[0].data(), bytes, kCoefficient,
                       tables);
      } else {
        kernels->add_scaled(candidate.data(), buffers[0].data(), bytes,
                            kCoefficient, tables);
      }
    } else if (operation == PrimitiveOperation::scale) {
      gf_vect_dot_prod_avx2(static_cast<int>(bytes), 1,
                            const_cast<Element*>(isal_table.data()),
                            source_pointers.data(), candidate.data());
    } else {
      gf_vect_mad_avx2(static_cast<int>(bytes), 1, 0,
                       const_cast<Element*>(isal_table.data()),
                       buffers[0].data(), candidate.data());
    }
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(bytes));
}

void RegisterTopLevel(const char* prefix) {
  for (const auto& [suffix, function] :
       std::array<std::pair<const char*, void (*)(benchmark::State&)>, 2>{{
           {"Encode", &ISALEncode},
           {"DecodeMax", &ISALDecodeMax},
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

void RegisterISALReferenceBenchmarks() {
#if defined(GF256_HIGH_LEVEL_BENCHMARKS_ONLY)
  RegisterTopLevel("RS/ISA-L/Native/");
#else
  RegisterTopLevel("ISA-L/Native/");
  auto* repair = benchmark::RegisterBenchmark("ISA-L/Native/RepairPrepared",
                                              &ISALRepairPrepared);
  for (const auto& test_case : kLCHComparisonCases) {
    repair->Args(
        {test_case.data_count, test_case.recovery_count, test_case.bytes});
  }
  repair->ArgNames({"K", "R", "bytes"});
  if (!lch::BackendAvailable(lch::Backend::avx2)) {
    return;
  }
  for (const auto& [name, operation, implementation] : std::array<
           std::tuple<const char*, PrimitiveOperation, PrimitiveImplementation>,
           4>{{
           {"Kernel/Scale/Owned/AVX2", PrimitiveOperation::scale,
            PrimitiveImplementation::owned},
           {"Kernel/Scale/ISA-L/AVX2", PrimitiveOperation::scale,
            PrimitiveImplementation::isal},
           {"Kernel/AddScaled/Owned/AVX2", PrimitiveOperation::add_scaled,
            PrimitiveImplementation::owned},
           {"Kernel/AddScaled/ISA-L/AVX2", PrimitiveOperation::add_scaled,
            PrimitiveImplementation::isal},
       }}) {
    auto* registered = benchmark::RegisterBenchmark(
        name, [operation, implementation](benchmark::State& state) {
          PrimitiveBenchmark(state, operation, implementation);
        });
    for (const int64_t bytes : {64, 128, 1024, 4096, 65536}) {
      registered->Arg(bytes);
    }
    registered->ArgName("bytes");
  }
#endif
}

}  // namespace gf256_benchmarks
