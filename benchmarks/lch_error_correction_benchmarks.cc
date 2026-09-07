#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <random>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "benchmark/benchmark.h"
#include "lin_chung_han/codeword_transform_internal.h"
#include "lin_chung_han/transform.h"
#include "reed_solomon/error_correction/internal.h"
#include "reed_solomon/lch_encoder.h"
#include "rs_benchmark_cases.h"

namespace {

using gf2p8::Element;
using gf2p8::lch::Backend;
using gf2p8::lch::Radix;
using gf2p8::lch::Status;
using gf2p8::lch::detail::FFTCodewordBlocks;
#if defined(GF256_ENABLE_CODEWORD_CANTOR_AFFINE_EXPERIMENT)
using gf2p8::lch::detail::FFTCodewordBlocksCantorAffine;
#endif
using gf2p8::lch::detail::IFFTCodewordBlocks;
#if defined(GF256_ENABLE_CODEWORD_CANTOR_AFFINE_EXPERIMENT)
using gf2p8::lch::detail::IFFTCodewordBlocksCantorAffine;
#endif
using gf2p8::rs::LCHDecoder;
using gf2p8::rs::LCHEncoder;
using gf2p8::rs::detail::error_correction::CorrectBatch;
using gf2p8::rs::detail::error_correction::CorrectionStatus;
using gf2p8::rs::detail::error_correction::CorrectOne;

struct FullCodeCase {
  size_t data_count;
  size_t recovery_count;
};

enum class ErrorProfile {
  clean,
  data_one,
  recovery_one,
  mixed_two,
  data_max,
  recovery_max,
  max_radius,
};

constexpr std::array<ErrorProfile, 7> kErrorProfiles = {
    ErrorProfile::clean,        ErrorProfile::data_one,
    ErrorProfile::recovery_one, ErrorProfile::mixed_two,
    ErrorProfile::data_max,     ErrorProfile::recovery_max,
    ErrorProfile::max_radius,
};

struct Error {
  size_t position;
  Element magnitude;
};

std::vector<Element*> MutablePointers(
    std::vector<std::vector<Element>>& shards) {
  std::vector<Element*> pointers;
  pointers.reserve(shards.size());
  for (auto& shard : shards) {
    pointers.push_back(shard.data());
  }
  return pointers;
}

std::vector<const Element*> ConstPointers(
    const std::vector<std::vector<Element>>& shards) {
  std::vector<const Element*> pointers;
  pointers.reserve(shards.size());
  for (const auto& shard : shards) {
    pointers.push_back(shard.data());
  }
  return pointers;
}

const char* ProfileName(ErrorProfile profile) {
  switch (profile) {
    case ErrorProfile::clean:
      return "Clean";
    case ErrorProfile::data_one:
      return "DataOne";
    case ErrorProfile::recovery_one:
      return "RecoveryOne";
    case ErrorProfile::mixed_two:
      return "MixedTwo";
    case ErrorProfile::data_max:
      return "DataMax";
    case ErrorProfile::recovery_max:
      return "RecoveryMax";
    case ErrorProfile::max_radius:
      return "MaxRadius";
  }
  return "Unknown";
}

std::vector<Error> MakeErrors(FullCodeCase code, ErrorProfile profile) {
  if (profile == ErrorProfile::clean) {
    return {};
  }

  if (profile == ErrorProfile::data_one) {
    return {{code.data_count / 3, Element{0x5b}}};
  }
  if (profile == ErrorProfile::recovery_one) {
    return {{code.data_count + code.recovery_count / 3, Element{0xa7}}};
  }
  if (profile == ErrorProfile::mixed_two) {
    return {
        {code.data_count / 3, Element{0x5b}},
        {code.data_count + code.recovery_count / 3, Element{0xa7}},
    };
  }

  const size_t error_count = code.recovery_count / 2;
  const size_t data_errors = profile == ErrorProfile::data_max ? error_count
                             : profile == ErrorProfile::recovery_max
                                 ? 0
                                 : (error_count + 1) / 2;
  const size_t recovery_errors = error_count - data_errors;
  std::vector<size_t> data_positions(code.data_count);
  std::vector<size_t> recovery_positions(code.recovery_count);
  std::iota(data_positions.begin(), data_positions.end(), size_t{0});
  std::iota(recovery_positions.begin(), recovery_positions.end(),
            code.data_count);

  std::mt19937 random(static_cast<uint32_t>(
      0x9e3779b9U ^ (code.data_count << 8U) ^ code.recovery_count));
  std::shuffle(data_positions.begin(), data_positions.end(), random);
  std::shuffle(recovery_positions.begin(), recovery_positions.end(), random);

  std::vector<Error> errors;
  errors.reserve(error_count);
  for (size_t i = 0; i < data_errors; ++i) {
    const Element magnitude =
        static_cast<Element>(((29 * data_positions[i] + 71 * i) % 255) + 1);
    errors.push_back({data_positions[i], magnitude});
  }
  for (size_t i = 0; i < recovery_errors; ++i) {
    const Element magnitude =
        static_cast<Element>(((43 * recovery_positions[i] + 97 * i) % 255) + 1);
    errors.push_back({recovery_positions[i], magnitude});
  }
  return errors;
}

class ScalarCorrectionInput {
 public:
  ScalarCorrectionInput(FullCodeCase code, ErrorProfile profile)
      : code_(code),
        encoder_(code.data_count, code.recovery_count),
        decoder_(code.data_count, code.recovery_count),
        clean_data_(code.data_count),
        data_(code.data_count),
        recovery_(code.recovery_count),
        mask_(code.data_count + code.recovery_count),
        expected_mask_(code.data_count + code.recovery_count),
        errors_(MakeErrors(code, profile)) {
    if (!encoder_.Valid() || !decoder_.Valid()) {
      error_ = "invalid benchmark dimensions";
      return;
    }

    for (size_t i = 0; i < clean_data_.size(); ++i) {
      clean_data_[i] = static_cast<Element>(
          (131 * i + 17 * code.data_count + code.recovery_count) & 0xffU);
    }
    if (!Encode()) {
      error_ = "failed to encode benchmark codeword";
      return;
    }

    data_ = clean_data_;
    for (const Error& error : errors_) {
      expected_mask_[error.position] = 1;
      if (error.position < code_.data_count) {
        data_[error.position] ^= error.magnitude;
        has_data_errors_ = true;
      } else {
        recovery_[error.position - code_.data_count] ^= error.magnitude;
      }
    }

    if (!Validate()) {
      error_ = "error-correction benchmark fixture failed validation";
    }
  }

  bool Valid() const { return error_.empty(); }
  const char* ErrorMessage() const { return error_.c_str(); }
  size_t ErrorCount() const { return errors_.size(); }
  bool HasDataErrors() const { return has_data_errors_; }

  gf2p8::rs::detail::error_correction::CorrectionResult Correct() {
    return CorrectOne(decoder_, data_, recovery_, mask_);
  }

  void RestoreDataErrors() {
    if (!has_data_errors_) {
      return;
    }
    for (const Error& error : errors_) {
      if (error.position < code_.data_count) {
        data_[error.position] ^= error.magnitude;
      }
    }
  }

 private:
  bool Encode() {
    std::vector<const Element*> data_pointers(code_.data_count);
    std::vector<Element*> recovery_pointers(code_.recovery_count);
    for (size_t i = 0; i < code_.data_count; ++i) {
      data_pointers[i] = &clean_data_[i];
    }
    for (size_t i = 0; i < code_.recovery_count; ++i) {
      recovery_pointers[i] = &recovery_[i];
    }
    std::vector<Element> workspace(encoder_.WorkspaceSize(1));
    return encoder_.Encode(data_pointers, recovery_pointers, 1, workspace,
                           Backend::scalar, Radix::radix2) == Status::ok;
  }

  bool Validate() const {
    std::vector<Element> data = data_;
    std::vector<Element> recovery = recovery_;
    std::vector<uint8_t> mask(mask_.size(), 0xa5);
    const auto result = CorrectOne(decoder_, data, recovery, mask);
    return result.status == CorrectionStatus::ok &&
           result.error_count == errors_.size() && data == clean_data_ &&
           recovery == recovery_ && mask == expected_mask_;
  }

  FullCodeCase code_;
  LCHEncoder encoder_;
  LCHDecoder decoder_;
  std::vector<Element> clean_data_;
  std::vector<Element> data_;
  std::vector<Element> recovery_;
  std::vector<uint8_t> mask_;
  std::vector<uint8_t> expected_mask_;
  std::vector<Error> errors_;
  bool has_data_errors_ = false;
  std::string error_;
};

struct BatchError {
  size_t byte;
  size_t position;
  Element magnitude;
};

std::vector<BatchError> MakeBatchErrors(FullCodeCase code,
                                        ErrorProfile profile,
                                        size_t byte_count) {
  std::vector<BatchError> errors;
  if (profile == ErrorProfile::clean) {
    return errors;
  }
  const size_t radius = code.recovery_count / 2;
  const size_t error_count =
      profile == ErrorProfile::data_one || profile == ErrorProfile::recovery_one
          ? 1
      : profile == ErrorProfile::mixed_two ? 2
                                           : radius;
  errors.reserve(byte_count * error_count);

  for (size_t byte = 0; byte < byte_count; ++byte) {
    std::vector<size_t> candidates;
    if (profile == ErrorProfile::data_one ||
        profile == ErrorProfile::data_max) {
      candidates.resize(code.data_count);
      std::iota(candidates.begin(), candidates.end(), size_t{0});
    } else if (profile == ErrorProfile::recovery_one ||
               profile == ErrorProfile::recovery_max) {
      candidates.resize(code.recovery_count);
      std::iota(candidates.begin(), candidates.end(), code.data_count);
    } else {
      candidates.resize(code.data_count + code.recovery_count);
      std::iota(candidates.begin(), candidates.end(), size_t{0});
    }
    std::mt19937 random(
        static_cast<uint32_t>(0xd1b50000U ^ (code.data_count << 12U) ^
                              (code.recovery_count << 5U) ^ byte));
    std::shuffle(candidates.begin(), candidates.end(), random);
    if (profile == ErrorProfile::mixed_two) {
      candidates[0] = byte % code.data_count;
      candidates[1] = code.data_count + byte % code.recovery_count;
    }
    for (size_t i = 0; i < error_count; ++i) {
      errors.push_back({
          .byte = byte,
          .position = candidates[i],
          .magnitude = static_cast<Element>(
              1 + ((61 * byte + 47 * i + code.recovery_count) % 255)),
      });
    }
  }
  return errors;
}

class BatchCorrectionInput {
 public:
  static constexpr size_t kByteCount = 32;

  BatchCorrectionInput(FullCodeCase code, ErrorProfile profile)
      : code_(code),
        encoder_(code.data_count, code.recovery_count),
        decoder_(code.data_count, code.recovery_count),
        clean_data_(code.data_count, std::vector<Element>(kByteCount)),
        data_(code.data_count, std::vector<Element>(kByteCount)),
        recovery_(code.recovery_count, std::vector<Element>(kByteCount)),
        mask_((code.data_count + code.recovery_count) * kByteCount),
        expected_mask_(mask_.size()),
        results_(kByteCount),
        errors_(MakeBatchErrors(code, profile, kByteCount)),
        errors_per_codeword_(errors_.size() / kByteCount) {
    if (!encoder_.Valid() || !decoder_.Valid() ||
        !gf2p8::lch::BackendAvailable(Backend::gfni256_affine)) {
      error_ = "GFNI batch backend is unavailable";
      return;
    }
    for (size_t shard = 0; shard < code.data_count; ++shard) {
      for (size_t byte = 0; byte < kByteCount; ++byte) {
        clean_data_[shard][byte] = static_cast<Element>(
            (131 * shard + 73 * byte + code.data_count) & 0xffU);
      }
    }
    data_ = clean_data_;
    auto data_input = ConstPointers(clean_data_);
    auto recovery_output = MutablePointers(recovery_);
    std::vector<Element> workspace(encoder_.WorkspaceSize(kByteCount));
    if (encoder_.Encode(data_input, recovery_output, kByteCount, workspace,
                        Backend::tuned, Radix::radix4) != Status::ok) {
      error_ = "failed to encode batch benchmark input";
      return;
    }
    for (const BatchError& error : errors_) {
      expected_mask_[error.position * kByteCount + error.byte] = 1;
      if (error.position < code.data_count) {
        data_[error.position][error.byte] ^= error.magnitude;
        has_data_errors_ = true;
      } else {
        recovery_[error.position - code.data_count][error.byte] ^=
            error.magnitude;
      }
    }
    data_pointers_ = MutablePointers(data_);
    recovery_pointers_ = ConstPointers(recovery_);
    if (!Validate()) {
      error_ = "batch error-correction fixture failed validation";
    }
  }

  bool Valid() const { return error_.empty(); }
  const char* ErrorMessage() const { return error_.c_str(); }
  size_t ErrorsPerCodeword() const { return errors_per_codeword_; }
  bool HasDataErrors() const { return has_data_errors_; }

  CorrectionStatus Correct() {
    return CorrectBatch(decoder_, data_pointers_, recovery_pointers_,
                        kByteCount, results_, mask_);
  }

  void RestoreDataErrors() {
    for (const BatchError& error : errors_) {
      if (error.position < code_.data_count) {
        data_[error.position][error.byte] ^= error.magnitude;
      }
    }
  }

  const auto& Results() const { return results_; }

  bool OutcomesMatch() const {
    return std::all_of(results_.begin(), results_.end(),
                       [&](const auto& result) {
                         return result.status == CorrectionStatus::ok &&
                                result.error_count == errors_per_codeword_;
                       });
  }

 private:
  bool Validate() const {
    auto data = data_;
    auto recovery = recovery_;
    auto data_pointers = MutablePointers(data);
    auto recovery_pointers = ConstPointers(recovery);
    std::vector<gf2p8::rs::detail::error_correction::CorrectionResult> results(
        kByteCount);
    std::vector<uint8_t> mask(mask_.size(), 0xa5);
    if (CorrectBatch(decoder_, data_pointers, recovery_pointers, kByteCount,
                     results, mask) != CorrectionStatus::ok ||
        data != clean_data_ || recovery != recovery_ ||
        mask != expected_mask_) {
      return false;
    }
    return std::all_of(results.begin(), results.end(), [&](const auto& result) {
      return result.status == CorrectionStatus::ok &&
             result.error_count == errors_per_codeword_;
    });
  }

  FullCodeCase code_;
  LCHEncoder encoder_;
  LCHDecoder decoder_;
  std::vector<std::vector<Element>> clean_data_;
  std::vector<std::vector<Element>> data_;
  std::vector<std::vector<Element>> recovery_;
  std::vector<Element*> data_pointers_;
  std::vector<const Element*> recovery_pointers_;
  std::vector<uint8_t> mask_;
  std::vector<uint8_t> expected_mask_;
  std::vector<gf2p8::rs::detail::error_correction::CorrectionResult> results_;
  std::vector<BatchError> errors_;
  size_t errors_per_codeword_ = 0;
  bool has_data_errors_ = false;
  std::string error_;
};

void BenchmarkScalarCorrection(benchmark::State& state,
                               FullCodeCase code,
                               ErrorProfile profile) {
  ScalarCorrectionInput input(code, profile);
  if (!input.Valid()) {
    state.SkipWithError(input.ErrorMessage());
    return;
  }

  for (auto _ : state) {
    auto result = input.Correct();
    benchmark::DoNotOptimize(result.status);
    benchmark::DoNotOptimize(result.error_count);
    benchmark::ClobberMemory();
    if (result.status != CorrectionStatus::ok ||
        result.error_count != input.ErrorCount()) {
      state.SkipWithError("error correction failed during benchmark");
      break;
    }

    if (input.HasDataErrors()) {
      state.PauseTiming();
      input.RestoreDataErrors();
      state.ResumeTiming();
    }
  }

  state.counters["K"] = static_cast<double>(code.data_count);
  state.counters["R"] = static_cast<double>(code.recovery_count);
  state.counters["errors"] = static_cast<double>(input.ErrorCount());
  state.SetItemsProcessed(state.iterations());
  state.SetBytesProcessed(state.iterations() *
                          static_cast<int64_t>(code.data_count));
}

void BenchmarkBatchCorrection(benchmark::State& state,
                              FullCodeCase code,
                              ErrorProfile profile) {
  BatchCorrectionInput input(code, profile);
  if (!input.Valid()) {
    state.SkipWithError(input.ErrorMessage());
    return;
  }

  for (auto _ : state) {
    CorrectionStatus status = input.Correct();
    benchmark::DoNotOptimize(status);
    benchmark::DoNotOptimize(input.Results().data());
    benchmark::ClobberMemory();
    if (status != CorrectionStatus::ok || !input.OutcomesMatch()) {
      state.SkipWithError("batch error correction failed during benchmark");
      break;
    }
    if (input.HasDataErrors()) {
      state.PauseTiming();
      input.RestoreDataErrors();
      state.ResumeTiming();
    }
  }

  state.counters["K"] = static_cast<double>(code.data_count);
  state.counters["R"] = static_cast<double>(code.recovery_count);
  state.counters["errors"] = static_cast<double>(input.ErrorsPerCodeword());
  state.counters["batch_bytes"] =
      static_cast<double>(BatchCorrectionInput::kByteCount);
  state.SetItemsProcessed(state.iterations() *
                          BatchCorrectionInput::kByteCount);
  state.SetBytesProcessed(
      state.iterations() * static_cast<int64_t>(code.data_count) *
      static_cast<int64_t>(BatchCorrectionInput::kByteCount));
}

void BenchmarkCodewordTransform(benchmark::State& state,
                                bool inverse,
                                Backend backend,
                                bool cantor_affine) {
  if ((backend != Backend::scalar || cantor_affine) &&
      !gf2p8::lch::BackendAvailable(Backend::gfni256_affine)) {
    state.SkipWithError("backend was not compiled");
    return;
  }
  const size_t value_count = static_cast<size_t>(state.range(0));
  const size_t block_size = static_cast<size_t>(state.range(1));
  std::mt19937 random(42);
  std::vector<Element> values(value_count);
  std::generate(values.begin(), values.end(),
                [&random] { return static_cast<Element>(random()); });
  const auto run = [&] {
#if defined(GF256_ENABLE_CODEWORD_CANTOR_AFFINE_EXPERIMENT)
    if (cantor_affine) {
      return inverse
                 ? IFFTCodewordBlocksCantorAffine(gf2p8::lch::Context::Shared(),
                                                  values, block_size, 0)
                 : FFTCodewordBlocksCantorAffine(gf2p8::lch::Context::Shared(),
                                                 values, block_size, 0);
    }
#else
    (void)cantor_affine;
#endif
    return inverse ? IFFTCodewordBlocks(gf2p8::lch::Context::Shared(), values,
                                        block_size, 0, backend)
                   : FFTCodewordBlocks(gf2p8::lch::Context::Shared(), values,
                                       block_size, 0, backend);
  };
  if (run() != Status::ok) {
    state.SkipWithError("transform validation failed");
    return;
  }

  for (auto _ : state) {
    auto status = run();
    benchmark::DoNotOptimize(status);
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) *
                          static_cast<int64_t>(value_count));
}

void RegisterScalarCorrectionBenchmarks() {
  for (const ErrorProfile profile : kErrorProfiles) {
    auto* registered = benchmark::RegisterBenchmark(
        std::string("LCH/Owned/ErrorCorrection/FDMA/TunedBlocks/") +
            ProfileName(profile),
        [profile](benchmark::State& state) {
          BenchmarkScalarCorrection(state,
                                    {static_cast<size_t>(state.range(0)),
                                     static_cast<size_t>(state.range(1))},
                                    profile);
        });
    for (const auto code : gf256_benchmarks::kFullCodeErrorCorrectionCases) {
      if (profile == ErrorProfile::mixed_two && code.recovery_count < 4) {
        continue;
      }
      registered->Args({code.data_count, code.recovery_count});
    }
    registered->ArgNames({"K", "R"});
  }
}

void RegisterBatchCorrectionBenchmarks() {
  for (const ErrorProfile profile : kErrorProfiles) {
    auto* registered = benchmark::RegisterBenchmark(
        std::string("LCH/Owned/ErrorCorrection/FDMA/BatchGFNI32/") +
            ProfileName(profile),
        [profile](benchmark::State& state) {
          BenchmarkBatchCorrection(state,
                                   {static_cast<size_t>(state.range(0)),
                                    static_cast<size_t>(state.range(1))},
                                   profile);
        });
    for (const auto code : gf256_benchmarks::kFullCodeErrorCorrectionCases) {
      if (profile == ErrorProfile::mixed_two && code.recovery_count < 4) {
        continue;
      }
      registered->Args({code.data_count, code.recovery_count});
    }
    registered->ArgNames({"K", "R"});
  }
}

void RegisterCodewordTransformBenchmarks() {
  struct Variant {
    Backend backend;
    const char* name;
    bool cantor_affine;
  };
  std::vector<Variant> variants = {
      {Backend::scalar, "Scalar", false},
      {Backend::gfni256_affine, "GFNI256Mul", false},
  };
#if defined(GF256_ENABLE_CODEWORD_CANTOR_AFFINE_EXPERIMENT)
  variants.push_back({Backend::gfni256_affine, "GFNI256CantorAffine", true});
#endif
  for (const Variant variant : variants) {
    for (const bool inverse : {false, true}) {
      auto* registered = benchmark::RegisterBenchmark(
          std::string("LCH/Codeword/") +
              (inverse ? "IFFTBlocks/" : "FFTBlocks/") + variant.name,
          [inverse, variant](benchmark::State& state) {
            BenchmarkCodewordTransform(state, inverse, variant.backend,
                                       variant.cantor_affine);
          });
      for (const auto code : gf256_benchmarks::kFullCodeErrorCorrectionCases) {
        if (variant.backend != Backend::scalar &&
            code.data_count + code.recovery_count < 32) {
          continue;
        }
        registered->Args(
            {code.data_count + code.recovery_count, code.recovery_count});
      }
      registered->ArgNames({"N", "R"});
    }
  }
}

const bool kScalarCorrectionBenchmarksRegistered = [] {
  RegisterScalarCorrectionBenchmarks();
  RegisterBatchCorrectionBenchmarks();
  RegisterCodewordTransformBenchmarks();
  return true;
}();

}  // namespace
