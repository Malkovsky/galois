#pragma once

#include <array>
#include <span>
#include <vector>

#include "lin_chung_han/transform_internal.h"
#include "reed_solomon/code_parameters.h"

namespace gf2p8::rs::detail::decoder {

using Erased = std::array<uint8_t, lch::Context::kFieldSize>;

struct Inputs {
  std::span<Element* const> data;
  std::span<const uint8_t> data_present;
  std::span<const Element* const> recovery;
  std::span<const uint8_t> recovery_present;
  size_t byte_count;
  std::span<Element> workspace;
};

std::vector<Element> MakeLowRateCoefficients(const CodeParameters& parameters);
std::vector<Element> MakeHighRateCoefficients(const CodeParameters& parameters);

lch::Status DecodeLowRate(const CodeParameters& parameters,
                          std::span<const Element> coefficients,
                          const Inputs& inputs,
                          const lch::detail::ResolvedKernels& kernels,
                          lch::Radix radix);
lch::Status DecodeHighRate(const CodeParameters& parameters,
                           std::span<const Element> coefficients,
                           const Inputs& inputs,
                           const lch::detail::ResolvedKernels& kernels,
                           lch::Radix radix);

lch::Status MultiplyCopy(Element* destination,
                         const Element* source,
                         size_t byte_count,
                         Element coefficient,
                         const lch::detail::ResolvedKernels& kernels,
                         const lch::Context& context);
lch::Status LowRateFormalDerivative(
    std::span<Element* const> coefficients,
    size_t byte_count,
    const lch::detail::ResolvedKernels& kernels);
void ErrorLocator(const lch::Context& context,
                  const Erased& erased,
                  size_t live_count,
                  std::span<Element, lch::Context::kFieldSize> values);
Element Inverse(const lch::Context& context, Element value);

}  // namespace gf2p8::rs::detail::decoder
