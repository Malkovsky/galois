#include <algorithm>
#include <array>

#include "reed_solomon/decoder/internal.h"

namespace gf2p8::rs::detail::decoder {

std::vector<Element> MakeHighRateCoefficients(
    const CodeParameters& parameters) {
  const size_t m = parameters.transform_size;
  const size_t n = parameters.mother_size;
  std::vector<Element> coefficients(parameters.data_count);
  const lch::Context& context = lch::Context::Shared();

  // Algorithm 3 normalization for the padded [N, N-M] mother code.
  Element common = 1;
  for (size_t block = m; block < n; block *= 2) {
    for (size_t i = 0; i < block; ++i) {
      common = MultiplyCantor(
          common, context.EvaluationPoint(block) ^ context.EvaluationPoint(i));
    }
  }
  for (size_t data_index = 0; data_index < parameters.data_count;
       ++data_index) {
    const size_t position = m + data_index;
    Element denominator = common;
    for (size_t i = 0; i < m; ++i) {
      denominator =
          MultiplyCantor(denominator, context.EvaluationPoint(position) ^
                                          context.EvaluationPoint(i));
    }
    coefficients[data_index] = Inverse(context, denominator);
  }
  return coefficients;
}

lch::Status DecodeHighRate(const CodeParameters& parameters,
                           std::span<const Element> coefficients,
                           const Inputs& inputs,
                           const lch::detail::ResolvedKernels& kernels,
                           lch::Radix radix) {
  const size_t k = parameters.data_count;
  const size_t r = parameters.recovery_count;
  const size_t m = parameters.transform_size;
  const size_t n = parameters.mother_size;
  const lch::Context& context = lch::Context::Shared();

  // Virtual layout: [R recovery][M-R punctured][K data][known-zero tail].
  Erased erased{};
  for (size_t i = 0; i < r; ++i) {
    erased[i] = inputs.recovery_present[i] == 0;
  }
  for (size_t i = r; i < m; ++i) {
    erased[i] = 1;
  }
  Element zero_byte = 0;
  Element* zero_source = inputs.byte_count == 0 ? &zero_byte : nullptr;
  for (size_t i = 0; i < k; ++i) {
    erased[m + i] = inputs.data_present[i] == 0;
    if (erased[m + i] != 0 && zero_source == nullptr) {
      zero_source = inputs.data[i];
    }
  }
  std::fill_n(zero_source, inputs.byte_count, 0);

  auto locator = std::span<Element, lch::Context::kFieldSize>(
      inputs.workspace.data(), 256);
  ErrorLocator(context, erased, m + k, locator);
  std::array<Element*, lch::Context::kFieldSize> scratch{};
  Element* const base =
      inputs.byte_count == 0 ? &zero_byte : inputs.workspace.data() + 256;
  for (size_t i = 0; i < n; ++i) {
    scratch[i] = base + i * inputs.byte_count;
  }

  std::array<const Element*, lch::Context::kFieldSize> input{};
  auto accumulator = std::span<Element* const>(scratch).first(m);
  lch::Status status = lch::Status::ok;

  // Fold only the logical prefixes; punctured and shortened suffixes are zero.
  for (size_t i = 0; i < r; ++i) {
    input[i] = erased[i] != 0 ? zero_source : inputs.recovery[i];
  }
  status = lch::detail::IFFTResolved(
      context, std::span<const Element* const>(input).first(r), accumulator,
      inputs.byte_count, 0, kernels, radix);
  if (status != lch::Status::ok) {
    return status;
  }

  for (size_t data_offset = 0; data_offset < k; data_offset += m) {
    const size_t live_count = std::min(m, k - data_offset);
    bool has_present_input = false;
    for (size_t i = 0; i < live_count; ++i) {
      if (inputs.data_present[data_offset + i] == 0) {
        input[i] = zero_source;
      } else {
        input[i] = inputs.data[data_offset + i];
        has_present_input = true;
      }
    }
    if (!has_present_input) {
      continue;
    }
    const size_t position = m + data_offset;
    const auto source =
        std::span<const Element* const>(input).first(live_count);
    auto block = std::span<Element* const>(scratch).subspan(position, m);
    status =
        lch::detail::IFFTResolved(context, source, block, accumulator,
                                  inputs.byte_count, position, kernels, radix);
    if (status != lch::Status::ok) {
      return status;
    }
  }
  std::array<uint8_t, lch::Context::kFieldSize> requested{};
  for (size_t i = 0; i < r; ++i) {
    requested[i] = inputs.recovery_present[i] != 0;
  }
  status = lch::detail::FFTResolved(
      context, accumulator, inputs.byte_count, 0,
      std::span<const uint8_t>(requested).first(m), kernels, radix);
  if (status != lch::Status::ok) {
    return status;
  }

  for (size_t i = 0; i < m; ++i) {
    if (i >= r || erased[i] != 0) {
      std::fill_n(scratch[m + i], inputs.byte_count, 0);
    } else {
      kernels.scale(scratch[m + i], scratch[i], inputs.byte_count, locator[i],
                    context.Tables());
    }
  }
  auto polynomial = std::span<Element* const>(scratch).subspan(m, m);
  status = lch::detail::IFFTResolved(context, polynomial, inputs.byte_count, 0,
                                     r, kernels, radix);
  if (status != lch::Status::ok) {
    return status;
  }

  for (size_t data_offset = 0; data_offset < k; data_offset += m) {
    const size_t live_count = std::min(m, k - data_offset);
    std::fill_n(requested.begin(), m, 0);
    bool has_requested_output = false;
    for (size_t i = 0; i < live_count; ++i) {
      requested[i] = inputs.data_present[data_offset + i] == 0;
      has_requested_output |= requested[i] != 0;
    }
    if (!has_requested_output) {
      continue;
    }
    const size_t position = m + data_offset;
    auto output = data_offset == 0
                      ? std::span<Element* const>(scratch).first(m)
                      : std::span<Element* const>(scratch).subspan(position, m);
    for (size_t i = 0; i < m; ++i) {
      std::copy_n(polynomial[i], inputs.byte_count, output[i]);
    }
    status = lch::detail::FFTResolved(
        context, output, inputs.byte_count, position,
        std::span<const uint8_t>(requested).first(m), kernels, radix);
    if (status != lch::Status::ok) {
      return status;
    }
    for (size_t i = 0; i < live_count; ++i) {
      if (requested[i] == 0) {
        continue;
      }
      const size_t data_index = data_offset + i;
      const Element factor = MultiplyCantor(
          coefficients[data_index], Inverse(context, locator[m + data_index]));
      status = MultiplyCopy(inputs.data[data_index], output[i],
                            inputs.byte_count, factor, kernels, context);
      if (status != lch::Status::ok) {
        return status;
      }
    }
  }
  return lch::Status::ok;
}

}  // namespace gf2p8::rs::detail::decoder
