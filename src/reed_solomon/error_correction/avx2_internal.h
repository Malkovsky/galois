#pragma once

#if defined(__AVX2__)

#include <immintrin.h>

#include <array>
#include <cstddef>

#include "field.h"

namespace gf2p8::rs::detail::error_correction::avx2 {

struct PreparedMultiplier {
  alignas(32) __m256i multiples[8];
};

/**
 * @brief Multiplies 32 Cantor-coordinate lanes by one field coefficient.
 * @param values Input field elements.
 * @param coefficient Shared Cantor-coordinate multiplier.
 * @param tables Shared field multiplication tables.
 * @return The 32 lane-wise field products.
 */
inline __m256i MultiplyFixed(__m256i values,
                             Element coefficient,
                             const MultiplicationTables& tables) {
  if (coefficient == 0) {
    return _mm256_setzero_si256();
  }
  if (coefficient == 1) {
    return values;
  }
  const Element* row = tables.shuffle[coefficient].data();
  const __m256i low = _mm256_load_si256(reinterpret_cast<const __m256i*>(row));
  const __m256i high =
      _mm256_load_si256(reinterpret_cast<const __m256i*>(row + 32));
  const __m256i nibble_mask = _mm256_set1_epi8(0x0f);
  return _mm256_xor_si256(
      _mm256_shuffle_epi8(low, _mm256_and_si256(values, nibble_mask)),
      _mm256_shuffle_epi8(
          high, _mm256_and_si256(_mm256_srli_epi64(values, 4), nibble_mask)));
}

/**
 * @brief Precomputes the eight Cantor-basis multiples of 32 multipliers.
 * @param multiplier Lane-varying Cantor-coordinate multipliers.
 * @param tables Shared field multiplication tables.
 * @return Products of every multiplier lane by each Cantor basis vector.
 */
inline PreparedMultiplier PrepareMultiplier(
    __m256i multiplier,
    const MultiplicationTables& tables) {
  PreparedMultiplier prepared;
  prepared.multiples[0] = multiplier;
  const __m256i nibble_mask = _mm256_set1_epi8(0x0f);
  const __m256i low_indices = _mm256_and_si256(multiplier, nibble_mask);
  const __m256i high_indices =
      _mm256_and_si256(_mm256_srli_epi64(multiplier, 4), nibble_mask);
  for (size_t bit = 1; bit < 8; ++bit) {
    const Element* row = tables.shuffle[size_t{1} << bit].data();
    const __m256i low =
        _mm256_load_si256(reinterpret_cast<const __m256i*>(row));
    const __m256i high =
        _mm256_load_si256(reinterpret_cast<const __m256i*>(row + 32));
    prepared.multiples[bit] =
        _mm256_xor_si256(_mm256_shuffle_epi8(low, low_indices),
                         _mm256_shuffle_epi8(high, high_indices));
  }
  return prepared;
}

/**
 * @brief Multiplies lane selectors by prepared lane-varying multipliers.
 * @param selector Lane-varying Cantor-coordinate multiplicands.
 * @param prepared Precomputed basis multiples of the other operands.
 * @return The 32 lane-wise field products.
 */
inline __m256i MultiplyPrepared(__m256i selector,
                                const PreparedMultiplier& prepared) {
  const __m256i zero = _mm256_setzero_si256();
  __m256i result = _mm256_blendv_epi8(zero, prepared.multiples[7], selector);
  for (size_t bit = 7; bit-- > 0;) {
    selector = _mm256_add_epi8(selector, selector);
    result = _mm256_xor_si256(
        result, _mm256_blendv_epi8(zero, prepared.multiples[bit], selector));
  }
  return result;
}

/**
 * @brief Multiplies two vectors of lane-varying Cantor-coordinate elements.
 * @param first First vector of field elements.
 * @param second Second vector of field elements.
 * @param tables Shared field multiplication tables.
 * @return The 32 lane-wise field products.
 */
inline __m256i MultiplyVariable(__m256i first,
                                __m256i second,
                                const MultiplicationTables& tables) {
  return MultiplyPrepared(first, PrepareMultiplier(second, tables));
}

}  // namespace gf2p8::rs::detail::error_correction::avx2

#endif
