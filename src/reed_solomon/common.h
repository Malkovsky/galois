#pragma once

#include <array>
#include <cstddef>
#include <span>

#include "reed_solomon/lin_chung_han/transform_internal.h"

namespace gf2p8::rs::detail {

size_t ShardWorkspaceSize(size_t shard_count,
                          size_t byte_count,
                          size_t prefix_bytes = 0);
lch::Status ValidateBackend(lch::Backend backend);
lch::Backend ResolveBackend(lch::Backend backend, size_t byte_count);

lch::Status EncodeLCH(const lch::Context& context,
                      std::span<const Element* const> data,
                      std::span<Element* const> recovery,
                      size_t byte_count,
                      size_t transform_size,
                      std::span<Element> workspace,
                      const lch::detail::ResolvedKernels& kernels,
                      lch::Radix radix);

lch::Status MultiplyCopy(Element* destination,
                         const Element* source,
                         size_t byte_count,
                         Element coefficient,
                         const lch::detail::ResolvedKernels& kernels,
                         const lch::Context& context);
lch::Status FormalDerivative(std::span<Element* const> coefficients,
                             size_t byte_count,
                             const lch::detail::ResolvedKernels& kernels);
void ErrorLocator(const lch::Context& context,
                  const std::array<uint8_t, lch::Context::kFieldSize>& erased,
                  std::span<Element, lch::Context::kFieldSize> values);
void ErrorLocatorLogs(
    const lch::Context& context,
    const std::array<uint8_t, lch::Context::kFieldSize>& erased,
    size_t live_count,
    std::span<Element, lch::Context::kFieldSize> values);
Element Inverse(const lch::Context& context, Element value);

}  // namespace gf2p8::rs::detail
