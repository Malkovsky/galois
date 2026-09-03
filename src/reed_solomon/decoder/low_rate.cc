#include <algorithm>
#include <array>

#include "reed_solomon/decoder/internal.h"

namespace gf2p8::rs::detail::decoder {

std::vector<Element> MakeLowRateCoefficients(const CodeParameters& parameters) {
  const size_t transform_size = parameters.transform_size;
  const size_t block_count =
      1 + (parameters.recovery_count + transform_size - 1) / transform_size;
  std::vector<Element> coefficients(block_count);
  const lch::Context& context = lch::Context::Shared();

  Element numerator = 1;
  for (size_t i = 1; i < transform_size; ++i) {
    numerator = MultiplyCantor(numerator, context.EvaluationPoint(i));
  }
  for (size_t block = 1; block < block_count; ++block) {
    Element denominator = 1;
    for (size_t i = 0; i < transform_size; ++i) {
      denominator = MultiplyCantor(
          denominator, context.EvaluationPoint(block * transform_size) ^
                           context.EvaluationPoint(i));
    }
    coefficients[block] =
        MultiplyCantor(numerator, Inverse(context, denominator));
  }
  return coefficients;
}

lch::Status DecodeLowRate(const CodeParameters& parameters,
                          std::span<const Element> coefficients,
                          const Inputs& inputs,
                          const lch::detail::ResolvedKernels& kernels,
                          lch::Radix radix) {
  const size_t k = parameters.data_count;
  const size_t r = parameters.recovery_count;
  const size_t d = parameters.transform_size;
  const size_t n = parameters.mother_size;
  const lch::Context& context = lch::Context::Shared();

  Erased erased{};
  for (size_t i = 0; i < k; ++i) {
    erased[i] = inputs.data_present[i] == 0;
  }
  for (size_t i = 0; i < r; ++i) {
    erased[d + i] = inputs.recovery_present[i] == 0;
  }
  for (size_t position = d + r; position < n; ++position) {
    erased[position] = 1;
  }

  auto locator = std::span<Element, lch::Context::kFieldSize>(
      inputs.workspace.data(), 256);
  ErrorLocator(context, erased, n, locator);
  std::array<Element*, lch::Context::kFieldSize> scratch{};
  Element zero_byte = 0;
  Element* const base =
      inputs.byte_count == 0 ? &zero_byte : inputs.workspace.data() + 256;
  for (size_t i = 0; i < n; ++i) {
    scratch[i] = base + i * inputs.byte_count;
  }

  for (size_t position = 0; position < n; ++position) {
    const Element* source = nullptr;
    if (position < k && inputs.data_present[position] != 0) {
      source = inputs.data[position];
    } else if (position >= d && position < d + r) {
      const size_t recovery_index = position - d;
      if (inputs.recovery_present[recovery_index] != 0) {
        source = inputs.recovery[recovery_index];
      }
    }
    if (source == nullptr) {
      std::fill_n(scratch[position], inputs.byte_count, 0);
    } else {
      kernels.scale(scratch[position], source, inputs.byte_count,
                    locator[position], context.Tables());
    }
  }

  // Transform only public prefixes; fully punctured blocks contribute zero.
  lch::Status status = lch::detail::IFFTResolved(
      context, std::span<Element* const>(scratch).first(d), inputs.byte_count,
      0, k, kernels, radix);
  if (status != lch::Status::ok) {
    return status;
  }
  for (size_t recovery_offset = 0; recovery_offset < r; recovery_offset += d) {
    const size_t live_count = std::min(d, r - recovery_offset);
    status = lch::detail::IFFTResolved(
        context,
        std::span<Element* const>(scratch).subspan(d + recovery_offset, d),
        inputs.byte_count, d + recovery_offset, live_count, kernels, radix);
    if (status != lch::Status::ok) {
      return status;
    }
  }
  status = LowRateFormalDerivative(std::span<Element* const>(scratch).first(d),
                                   inputs.byte_count, kernels);
  if (status != lch::Status::ok) {
    return status;
  }
  for (size_t block = 1; block < coefficients.size(); ++block) {
    for (size_t i = 0; i < d; ++i) {
      kernels.add_scaled(scratch[i], scratch[block * d + i], inputs.byte_count,
                         coefficients[block], context.Tables());
    }
  }
  std::array<uint8_t, lch::Context::kFieldSize> requested{};
  for (size_t i = 0; i < k; ++i) {
    requested[i] = inputs.data_present[i] == 0;
  }
  status = lch::detail::FFTResolved(
      context, std::span<Element* const>(scratch).first(d), inputs.byte_count,
      0, std::span<const uint8_t>(requested).first(d), kernels, radix);
  if (status != lch::Status::ok) {
    return status;
  }
  for (size_t i = 0; i < k; ++i) {
    if (inputs.data_present[i] != 0) {
      continue;
    }
    status = MultiplyCopy(inputs.data[i], scratch[i], inputs.byte_count,
                          Inverse(context, locator[i]), kernels, context);
    if (status != lch::Status::ok) {
      return status;
    }
  }
  return lch::Status::ok;
}

}  // namespace gf2p8::rs::detail::decoder
