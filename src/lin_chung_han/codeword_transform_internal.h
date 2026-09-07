#pragma once

#include <array>
#include <cstddef>
#include <span>

#include "lin_chung_han/transform.h"

namespace gf2p8::lch::detail {

// Applies independent, adjacent block transforms to one contiguous scalar
// codeword. Block j is evaluated at evaluation_offset + j * block_size.
Status FFTCodewordBlocks(const Context& context,
                         std::span<Element> values,
                         size_t block_size,
                         size_t evaluation_offset,
                         Backend backend);

Status IFFTCodewordBlocks(const Context& context,
                          std::span<Element> values,
                          size_t block_size,
                          size_t evaluation_offset,
                          Backend backend);

#if defined(__GFNI__) && defined(__AVX2__)
// Private field-isomorphism helpers shared by position-vectorized kernels.
const std::array<Element, Context::kFieldSize>& CantorToAESMap();
const std::array<Element, Context::kFieldSize>& AESToCantorMap();
void ConvertCantorToAES(std::span<Element> values);
void ConvertAESToCantor(std::span<Element> values);
#endif

#if defined(GF256_ENABLE_CODEWORD_CANTOR_AFFINE_EXPERIMENT)
// Opt-in native-Cantor GFNI experiment. It applies packed fixed-factor affine
// matrices with a shift/mask butterfly network.
Status FFTCodewordBlocksCantorAffine(const Context& context,
                                     std::span<Element> values,
                                     size_t block_size,
                                     size_t evaluation_offset);

Status IFFTCodewordBlocksCantorAffine(const Context& context,
                                      std::span<Element> values,
                                      size_t block_size,
                                      size_t evaluation_offset);
#endif

}  // namespace gf2p8::lch::detail
