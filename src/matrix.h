#pragma once

#include "field.h"

#include <cstddef>
#include <functional>

namespace gf_2_8 {

/**
 * @brief baseline
 * @details
 * Performs textbook matrix multiplication with ikj loop order over row-major
 * matrices @p left and @p right of sizes m_i*m_k and m_k*m_j respectively using
 * @p fma to perform row scaling and addition. Multiplication result is put into
 * @p result.
 */
void MatMul(
    const element_t *left, const element_t *right, size_t m_i, size_t m_k,
    size_t m_j,
    std::function<void(element_t *, const element_t *, element_t, size_t)> fma,
    element_t *result);

/**
 * @brief Blocked matrix multiplication using lookup tables
 * @details
 * Performs row-major matrix multiplication over GF(256) using the same
 * low/high-nibble SIMD lookup-table FMA strategy as AddScaledRowSIMD. The wide
 * AVX2 path repacks right-hand panels into 64-byte tile-major order to improve
 * cache locality.
 */
void MatMulBlockedLUT(const element_t *left, const element_t *right,
                      size_t m_i, size_t m_k, size_t m_j, element_t *result);

/**
 * @brief Blocked matrix multiplication optimized for GFNI
 * @details
 * Performs row-major matrix multiplication over GF(256). Uses a small
 * register-blocked micro-kernel when AVX-512 GFNI is available at compile time
 * and falls back to MatMulBlockedLUT otherwise.
 */
void MatMulBlockedGFNI(const element_t *left, const element_t *right,
                       size_t m_i, size_t m_k, size_t m_j, element_t *result);

} // namespace gf_2_8
