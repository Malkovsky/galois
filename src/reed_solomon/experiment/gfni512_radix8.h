#pragma once

#include <span>

#include "lin_chung_han/experiment/gfni512_radix8.h"

namespace gf2p8::rs::detail::experiment::radix8 {

lch::Status EncodeLCH(const lch::Context& context,
                      std::span<const Element* const> data,
                      std::span<Element* const> recovery,
                      size_t byte_count,
                      size_t transform_size,
                      std::span<Element> workspace,
                      const lch::detail::ResolvedKernels& base,
                      const lch::detail::experiment::radix8::Kernels& radix8);

}  // namespace gf2p8::rs::detail::experiment::radix8
