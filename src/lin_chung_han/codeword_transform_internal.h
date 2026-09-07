#pragma once

#include <array>
#include <cstddef>
#include <span>

#include "lin_chung_han/transform.h"

namespace gf2p8::lch::detail {

/**
 * @brief Applies forward LCH transforms to adjacent scalar-codeword blocks.
 *
 * @param context LCH field context and skew tables.
 * @param values Contiguous values transformed in place.
 * @param block_size Power-of-two size of each independent transform.
 * @param evaluation_offset Aligned evaluation offset of the first block.
 * @param backend Scalar or tuned implementation to use.
 * @return Transform status; `Status::ok` on success.
 *
 * @note Block j is evaluated at
 * `evaluation_offset + j * block_size`.
 */
Status FFTCodewordBlocks(const Context& context,
                         std::span<Element> values,
                         size_t block_size,
                         size_t evaluation_offset,
                         Backend backend);

/**
 * @brief Applies inverse LCH transforms to adjacent scalar-codeword blocks.
 *
 * @param context LCH field context and skew tables.
 * @param values Contiguous values transformed in place.
 * @param block_size Power-of-two size of each independent transform.
 * @param evaluation_offset Aligned evaluation offset of the first block.
 * @param backend Scalar or tuned implementation to use.
 * @return Transform status; `Status::ok` on success.
 *
 * @note Block j is evaluated at
 * `evaluation_offset + j * block_size`.
 */
Status IFFTCodewordBlocks(const Context& context,
                          std::span<Element> values,
                          size_t block_size,
                          size_t evaluation_offset,
                          Backend backend);

#if defined(__GFNI__) && defined(__AVX2__)
/**
 * @brief Returns the native-Cantor to AES-coordinate byte map.
 * @return Immutable 256-entry field-isomorphism map.
 */
const std::array<Element, Context::kFieldSize>& CantorToAESMap();

/**
 * @brief Returns the AES-coordinate to native-Cantor byte map.
 * @return Immutable 256-entry inverse field-isomorphism map.
 */
const std::array<Element, Context::kFieldSize>& AESToCantorMap();

/**
 * @brief Converts native-Cantor bytes to AES-isomorphic coordinates in place.
 * @param values Field elements to convert.
 */
void ConvertCantorToAES(std::span<Element> values);

/**
 * @brief Converts AES-isomorphic bytes to native-Cantor coordinates in place.
 * @param values Field elements to convert.
 */
void ConvertAESToCantor(std::span<Element> values);
#endif

#if defined(GF256_ENABLE_CODEWORD_CANTOR_AFFINE_EXPERIMENT)
/**
 * @brief Applies the opt-in native-Cantor affine forward block transform.
 *
 * @param context LCH field context and skew tables.
 * @param values Contiguous values transformed in place.
 * @param block_size Power-of-two size of each independent transform.
 * @param evaluation_offset Aligned evaluation offset of the first block.
 * @return Transform status; `Status::ok` on success.
 *
 * @note This experiment uses packed fixed-factor GFNI affine matrices and a
 * shift/mask butterfly network.
 */
Status FFTCodewordBlocksCantorAffine(const Context& context,
                                     std::span<Element> values,
                                     size_t block_size,
                                     size_t evaluation_offset);

/**
 * @brief Applies the opt-in native-Cantor affine inverse block transform.
 *
 * @param context LCH field context and skew tables.
 * @param values Contiguous values transformed in place.
 * @param block_size Power-of-two size of each independent transform.
 * @param evaluation_offset Aligned evaluation offset of the first block.
 * @return Transform status; `Status::ok` on success.
 *
 * @note This experiment uses packed fixed-factor GFNI affine matrices and a
 * shift/mask butterfly network.
 */
Status IFFTCodewordBlocksCantorAffine(const Context& context,
                                      std::span<Element> values,
                                      size_t block_size,
                                      size_t evaluation_offset);
#endif

}  // namespace gf2p8::lch::detail
