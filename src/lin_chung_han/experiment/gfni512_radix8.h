#pragma once

#include <span>

#include "lin_chung_han/transform_internal.h"

namespace gf2p8::lch::detail::experiment::radix8 {

struct Kernels {
  using IFFTFunction = void (*)(Element*,
                                Element*,
                                Element*,
                                Element*,
                                Element*,
                                Element*,
                                Element*,
                                Element*,
                                size_t,
                                Element,
                                Element,
                                Element,
                                Element,
                                Element,
                                Element,
                                Element,
                                const MultiplicationTables&);
  using IFFTXorFunction = void (*)(const Element*,
                                   const Element*,
                                   const Element*,
                                   const Element*,
                                   const Element*,
                                   const Element*,
                                   const Element*,
                                   const Element*,
                                   Element*,
                                   Element*,
                                   Element*,
                                   Element*,
                                   Element*,
                                   Element*,
                                   Element*,
                                   Element*,
                                   size_t,
                                   Element,
                                   Element,
                                   Element,
                                   Element,
                                   Element,
                                   Element,
                                   Element,
                                   const MultiplicationTables&);

  IFFTFunction ifft_radix8;
  IFFTXorFunction ifft_radix8_xor;
};

const Kernels* ResolveKernels();

Status IFFTResolved(const Context& context,
                    std::span<Element* const> shards,
                    size_t byte_count,
                    size_t evaluation_offset,
                    size_t input_count,
                    const detail::ResolvedKernels& base,
                    const Kernels& radix8);

Status IFFTResolved(const Context& context,
                    std::span<const Element* const> input,
                    std::span<Element* const> work,
                    size_t byte_count,
                    size_t evaluation_offset,
                    const detail::ResolvedKernels& base,
                    const Kernels& radix8);

Status IFFTResolved(const Context& context,
                    std::span<const Element* const> input,
                    std::span<Element* const> work,
                    std::span<Element* const> xor_accumulator,
                    size_t byte_count,
                    size_t evaluation_offset,
                    const detail::ResolvedKernels& base,
                    const Kernels& radix8);

}  // namespace gf2p8::lch::detail::experiment::radix8
