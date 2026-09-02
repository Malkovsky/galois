#pragma once

#include <cstddef>
#include <functional>

#include "field.h"

namespace gf2p8 {

/**
 * @brief Multiplies matrices with a caller-provided row FMA.
 * @details Performs textbook matrix multiplication with `ikj` loop order over
 * row-major matrices. Coordinate semantics are determined by the supplied FMA.
 * @param left Left matrix with dimensions @p m_i by @p m_k.
 * @param right Right matrix with dimensions @p m_k by @p m_j.
 * @param m_i Number of rows in @p left and @p result.
 * @param m_k Shared inner dimension.
 * @param m_j Number of columns in @p right and @p result.
 * @param fma Row scale-and-add operation.
 * @param result Output matrix with dimensions @p m_i by @p m_j.
 */
void MatMul(const Element* left,
            const Element* right,
            size_t m_i,
            size_t m_k,
            size_t m_j,
            std::function<void(Element*, const Element*, Element, size_t)> fma,
            Element* result);

/**
 * @brief Blocked matrix multiplication using lookup tables
 * @details
 * Performs row-major matrix multiplication over GF(256) using the same
 * packed-nibble SIMD lookup-table FMA strategy as AddScaledRow. The wide
 * AVX2 path repacks right-hand panels into 64-byte tile-major order to improve
 * cache locality. All matrix elements use Cantor coordinates.
 * @param left Left matrix with dimensions @p m_i by @p m_k.
 * @param right Right matrix with dimensions @p m_k by @p m_j.
 * @param m_i Number of rows in @p left and @p result.
 * @param m_k Shared inner dimension.
 * @param m_j Number of columns in @p right and @p result.
 * @param result Output matrix with dimensions @p m_i by @p m_j.
 */
void MatMulBlockedLUT(const Element* left,
                      const Element* right,
                      size_t m_i,
                      size_t m_k,
                      size_t m_j,
                      Element* result);

/**
 * @brief Blocked matrix multiplication optimized for GFNI affine maps.
 * @details
 * Performs row-major matrix multiplication over GF(256). Uses a small
 * register-blocked fixed-coefficient affine micro-kernel for native Cantor
 * coordinates when AVX-512 GFNI is available at compile time and falls back to
 * MatMulBlockedLUT otherwise.
 * @param left Left matrix with dimensions @p m_i by @p m_k.
 * @param right Right matrix with dimensions @p m_k by @p m_j.
 * @param m_i Number of rows in @p left and @p result.
 * @param m_k Shared inner dimension.
 * @param m_j Number of columns in @p right and @p result.
 * @param result Output matrix with dimensions @p m_i by @p m_j.
 */
void MatMulBlockedGFNI(const Element* left,
                       const Element* right,
                       size_t m_i,
                       size_t m_k,
                       size_t m_j,
                       Element* result);

}  // namespace gf2p8
