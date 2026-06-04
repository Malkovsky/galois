#include "matrix.h"

#include <cstring>
#include <immintrin.h>
#include <memory>

namespace gf_2_8 {

#if defined(__AVX2__)
extern __m256i simd256_table_lo[256];
extern __m256i simd256_table_hi[256];
#endif

void MatMul(
    const element_t *left, const element_t *right, size_t m_i, size_t m_k,
    size_t m_j,
    std::function<void(element_t *, const element_t *, element_t, size_t)> fma,
    element_t *result) {
  std::memset(result, 0, m_i * m_j);
  for (size_t i = 0; i < m_i; ++i) {
    auto right_row = right;
    for (size_t k = 0; k < m_k; ++k, left++, right_row += m_j) {
      fma(result, right_row, *left, m_j);
    }
    result += m_j;
  }
}

void MatMulBlockedLUT(const element_t *left, const element_t *right,
                      size_t m_i, size_t m_k, size_t m_j, element_t *result) {
  constexpr size_t kCols = 2048;

  std::memset(result, 0, m_i * m_j);

#if defined(__AVX2__)
  const __m256i mask = _mm256_set1_epi8(0x0f);
#endif

  for (size_t j = 0; j < m_j; j += kCols) {
    const size_t cols = (m_j - j < kCols) ? (m_j - j) : kCols;
#if defined(__AVX2__)
    // Narrow panels benefit from loading each right-hand vector once and
    // applying it to four result rows.
    if (cols <= 256) {
      constexpr size_t kRows = 4;
      for (size_t i = 0; i < m_i; i += kRows) {
        const size_t rows = (m_i - i < kRows) ? (m_i - i) : kRows;
        for (size_t k = 0; k < m_k; ++k) {
          const element_t *right_panel = right + k * m_j + j;
          const element_t a0 = rows > 0 ? left[(i + 0) * m_k + k] : 0;
          const element_t a1 = rows > 1 ? left[(i + 1) * m_k + k] : 0;
          const element_t a2 = rows > 2 ? left[(i + 2) * m_k + k] : 0;
          const element_t a3 = rows > 3 ? left[(i + 3) * m_k + k] : 0;

          const __m256i lo0 = simd256_table_lo[a0];
          const __m256i hi0 = simd256_table_hi[a0];
          const __m256i lo1 = simd256_table_lo[a1];
          const __m256i hi1 = simd256_table_hi[a1];
          const __m256i lo2 = simd256_table_lo[a2];
          const __m256i hi2 = simd256_table_hi[a2];
          const __m256i lo3 = simd256_table_lo[a3];
          const __m256i hi3 = simd256_table_hi[a3];

          size_t c = 0;
          for (; c + 32 <= cols; c += 32) {
            const __m256i y_reg = _mm256_loadu_si256(
                reinterpret_cast<const __m256i *>(right_panel + c));
            const __m256i y_lo = _mm256_and_si256(y_reg, mask);
            const __m256i y_hi =
                _mm256_and_si256(_mm256_srli_epi64(y_reg, 4), mask);

#define GF256_APPLY_AVX2_ROW(row)                                             \
  if (a##row != 0) {                                                          \
    element_t *out = result + (i + row) * m_j + j + c;                        \
    const __m256i product = _mm256_xor_si256(                                 \
        _mm256_shuffle_epi8(lo##row, y_lo),                                   \
        _mm256_shuffle_epi8(hi##row, y_hi));                                  \
    const __m256i old =                                                       \
        _mm256_loadu_si256(reinterpret_cast<const __m256i *>(out));           \
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(out),                     \
                        _mm256_xor_si256(old, product));                      \
  }

            GF256_APPLY_AVX2_ROW(0)
            GF256_APPLY_AVX2_ROW(1)
            GF256_APPLY_AVX2_ROW(2)
            GF256_APPLY_AVX2_ROW(3)

#undef GF256_APPLY_AVX2_ROW
          }

          if (c < cols) {
            if (a0 != 0) {
              AddScaledRowBase(result + (i + 0) * m_j + j + c, right_panel + c,
                               a0, cols - c);
            }
            if (a1 != 0) {
              AddScaledRowBase(result + (i + 1) * m_j + j + c, right_panel + c,
                               a1, cols - c);
            }
            if (a2 != 0) {
              AddScaledRowBase(result + (i + 2) * m_j + j + c, right_panel + c,
                               a2, cols - c);
            }
            if (a3 != 0) {
              AddScaledRowBase(result + (i + 3) * m_j + j + c, right_panel + c,
                               a3, cols - c);
            }
          }
        }
      }
      continue;
    }

    // Wide AVX2 panels use a small dot-product tile. The eight accumulators
    // keep 4x64 result bytes in registers across the k loop, avoiding repeated
    // result read/modify/write traffic while each right-hand tile feeds four
    // output rows. For large panels, the right-hand matrix is repacked into
    // 64-byte tile-major order so the inner k loop reads B linearly instead of
    // taking one cache line from every row with a wide stride.
    if (cols >= 512) {
      constexpr size_t kRows = 4;
      constexpr size_t kTileCols = 64;
      const size_t vector_cols = cols & ~(kTileCols - 1);
      const size_t tiles = vector_cols / kTileCols;
      const bool pack_right = m_i > kRows;
      std::unique_ptr<element_t[]> packed_right;

      if (pack_right) {
        packed_right.reset(new element_t[tiles * m_k * kTileCols]);
        for (size_t k = 0; k < m_k; ++k) {
          const element_t *right_row = right + k * m_j + j;
          for (size_t tile = 0; tile < tiles; ++tile) {
            std::memcpy(packed_right.get() + (tile * m_k + k) * kTileCols,
                        right_row + tile * kTileCols, kTileCols);
          }
        }
      }

      size_t c = 0;
      for (; c + kTileCols <= cols; c += kTileCols) {
        const size_t tile = c / kTileCols;
        const element_t *right_tile_base =
            pack_right ? packed_right.get() + tile * m_k * kTileCols
                       : right + j + c;
        const size_t right_stride = pack_right ? kTileCols : m_j;

        for (size_t i = 0; i < m_i; i += kRows) {
          const size_t rows = (m_i - i < kRows) ? (m_i - i) : kRows;
          __m256i acc00 = _mm256_setzero_si256();
          __m256i acc01 = _mm256_setzero_si256();
          __m256i acc10 = _mm256_setzero_si256();
          __m256i acc11 = _mm256_setzero_si256();
          __m256i acc20 = _mm256_setzero_si256();
          __m256i acc21 = _mm256_setzero_si256();
          __m256i acc30 = _mm256_setzero_si256();
          __m256i acc31 = _mm256_setzero_si256();

          for (size_t k = 0; k < m_k; ++k) {
            const element_t *right_tile = right_tile_base + k * right_stride;
            const __m256i y0 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i *>(right_tile));
            const __m256i y1 = _mm256_loadu_si256(
                reinterpret_cast<const __m256i *>(right_tile + 32));
            const __m256i y0_lo = _mm256_and_si256(y0, mask);
            const __m256i y0_hi =
                _mm256_and_si256(_mm256_srli_epi64(y0, 4), mask);
            const __m256i y1_lo = _mm256_and_si256(y1, mask);
            const __m256i y1_hi =
                _mm256_and_si256(_mm256_srli_epi64(y1, 4), mask);

#define GF256_ACCUM_64_AVX2_ROW(row)                                          \
  if (rows > row) {                                                           \
    const element_t a = left[(i + row) * m_k + k];                            \
    if (a != 0) {                                                             \
      const __m256i table_lo = simd256_table_lo[a];                           \
      const __m256i table_hi = simd256_table_hi[a];                           \
      acc##row##0 = _mm256_xor_si256(                                         \
          acc##row##0, _mm256_xor_si256(_mm256_shuffle_epi8(table_lo, y0_lo), \
                                        _mm256_shuffle_epi8(table_hi, y0_hi))); \
      acc##row##1 = _mm256_xor_si256(                                         \
          acc##row##1, _mm256_xor_si256(_mm256_shuffle_epi8(table_lo, y1_lo), \
                                        _mm256_shuffle_epi8(table_hi, y1_hi))); \
    }                                                                         \
  }

            GF256_ACCUM_64_AVX2_ROW(0)
            GF256_ACCUM_64_AVX2_ROW(1)
            GF256_ACCUM_64_AVX2_ROW(2)
            GF256_ACCUM_64_AVX2_ROW(3)

#undef GF256_ACCUM_64_AVX2_ROW
          }

          if (rows > 0) {
            element_t *out = result + (i + 0) * m_j + j + c;
            _mm256_storeu_si256(reinterpret_cast<__m256i *>(out), acc00);
            _mm256_storeu_si256(reinterpret_cast<__m256i *>(out + 32), acc01);
          }
          if (rows > 1) {
            element_t *out = result + (i + 1) * m_j + j + c;
            _mm256_storeu_si256(reinterpret_cast<__m256i *>(out), acc10);
            _mm256_storeu_si256(reinterpret_cast<__m256i *>(out + 32), acc11);
          }
          if (rows > 2) {
            element_t *out = result + (i + 2) * m_j + j + c;
            _mm256_storeu_si256(reinterpret_cast<__m256i *>(out), acc20);
            _mm256_storeu_si256(reinterpret_cast<__m256i *>(out + 32), acc21);
          }
          if (rows > 3) {
            element_t *out = result + (i + 3) * m_j + j + c;
            _mm256_storeu_si256(reinterpret_cast<__m256i *>(out), acc30);
            _mm256_storeu_si256(reinterpret_cast<__m256i *>(out + 32), acc31);
          }
        }
      }

      if (c < cols) {
        for (size_t i = 0; i < m_i; i += kRows) {
          const size_t rows = (m_i - i < kRows) ? (m_i - i) : kRows;
          for (size_t k = 0; k < m_k; ++k) {
            const element_t *right_tail = right + k * m_j + j + c;
            for (size_t r = 0; r < rows; ++r) {
              const element_t a = left[(i + r) * m_k + k];
              if (a != 0) {
                AddScaledRowSIMD(result + (i + r) * m_j + j + c, right_tail, a,
                                 cols - c);
              }
            }
          }
        }
      }
      continue;
    }
#endif

    // For 2048-byte rows, eight result rows keep about 16 KiB of result data
    // hot while each right-hand row is reused across the row block.
    constexpr size_t kRows = 8;
    for (size_t i = 0; i < m_i; i += kRows) {
      const size_t rows = (m_i - i < kRows) ? (m_i - i) : kRows;
      for (size_t k = 0; k < m_k; ++k) {
        const element_t *right_panel = right + k * m_j + j;
        for (size_t r = 0; r < rows; ++r) {
          const element_t a = left[(i + r) * m_k + k];
          if (a != 0) {
            AddScaledRowSIMD(result + (i + r) * m_j + j, right_panel, a, cols);
          }
        }
      }
    }
  }
}

void MatMulBlockedGFNI(const element_t *left, const element_t *right,
                       size_t m_i, size_t m_k, size_t m_j,
                       element_t *result) {
#if defined(__GFNI__) && defined(__AVX512F__) && defined(__AVX512BW__)
  // Eight AVX-512 accumulators still fit comfortably in registers and halve
  // right-panel traffic compared with the original four-row micro-kernel.
  constexpr size_t kRows = 8;
  constexpr size_t kCols = 64;

  std::memset(result, 0, m_i * m_j);

  for (size_t i = 0; i < m_i; i += kRows) {
    const size_t rows = (m_i - i < kRows) ? (m_i - i) : kRows;
    for (size_t j = 0; j < m_j; j += kCols) {
      const size_t cols = (m_j - j < kCols) ? (m_j - j) : kCols;
      const uint64_t mask_bits = (cols == kCols) ? ~0ULL : ((1ULL << cols) - 1);
      const __mmask64 mask = static_cast<__mmask64>(mask_bits);

      __m512i acc0 = _mm512_setzero_si512();
      __m512i acc1 = _mm512_setzero_si512();
      __m512i acc2 = _mm512_setzero_si512();
      __m512i acc3 = _mm512_setzero_si512();
      __m512i acc4 = _mm512_setzero_si512();
      __m512i acc5 = _mm512_setzero_si512();
      __m512i acc6 = _mm512_setzero_si512();
      __m512i acc7 = _mm512_setzero_si512();

      for (size_t k = 0; k < m_k; ++k) {
        const __m512i b =
            _mm512_maskz_loadu_epi8(mask, right + k * m_j + j);

        if (rows > 0) {
          const __m512i a =
              _mm512_set1_epi8(static_cast<char>(left[(i + 0) * m_k + k]));
          acc0 = _mm512_xor_si512(acc0, _mm512_gf2p8mul_epi8(a, b));
        }
        if (rows > 1) {
          const __m512i a =
              _mm512_set1_epi8(static_cast<char>(left[(i + 1) * m_k + k]));
          acc1 = _mm512_xor_si512(acc1, _mm512_gf2p8mul_epi8(a, b));
        }
        if (rows > 2) {
          const __m512i a =
              _mm512_set1_epi8(static_cast<char>(left[(i + 2) * m_k + k]));
          acc2 = _mm512_xor_si512(acc2, _mm512_gf2p8mul_epi8(a, b));
        }
        if (rows > 3) {
          const __m512i a =
              _mm512_set1_epi8(static_cast<char>(left[(i + 3) * m_k + k]));
          acc3 = _mm512_xor_si512(acc3, _mm512_gf2p8mul_epi8(a, b));
        }
        if (rows > 4) {
          const __m512i a =
              _mm512_set1_epi8(static_cast<char>(left[(i + 4) * m_k + k]));
          acc4 = _mm512_xor_si512(acc4, _mm512_gf2p8mul_epi8(a, b));
        }
        if (rows > 5) {
          const __m512i a =
              _mm512_set1_epi8(static_cast<char>(left[(i + 5) * m_k + k]));
          acc5 = _mm512_xor_si512(acc5, _mm512_gf2p8mul_epi8(a, b));
        }
        if (rows > 6) {
          const __m512i a =
              _mm512_set1_epi8(static_cast<char>(left[(i + 6) * m_k + k]));
          acc6 = _mm512_xor_si512(acc6, _mm512_gf2p8mul_epi8(a, b));
        }
        if (rows > 7) {
          const __m512i a =
              _mm512_set1_epi8(static_cast<char>(left[(i + 7) * m_k + k]));
          acc7 = _mm512_xor_si512(acc7, _mm512_gf2p8mul_epi8(a, b));
        }
      }

      _mm512_mask_storeu_epi8(result + (i + 0) * m_j + j, mask, acc0);
      if (rows > 1) {
        _mm512_mask_storeu_epi8(result + (i + 1) * m_j + j, mask, acc1);
      }
      if (rows > 2) {
        _mm512_mask_storeu_epi8(result + (i + 2) * m_j + j, mask, acc2);
      }
      if (rows > 3) {
        _mm512_mask_storeu_epi8(result + (i + 3) * m_j + j, mask, acc3);
      }
      if (rows > 4) {
        _mm512_mask_storeu_epi8(result + (i + 4) * m_j + j, mask, acc4);
      }
      if (rows > 5) {
        _mm512_mask_storeu_epi8(result + (i + 5) * m_j + j, mask, acc5);
      }
      if (rows > 6) {
        _mm512_mask_storeu_epi8(result + (i + 6) * m_j + j, mask, acc6);
      }
      if (rows > 7) {
        _mm512_mask_storeu_epi8(result + (i + 7) * m_j + j, mask, acc7);
      }
    }
  }
#else
  MatMulBlockedLUT(left, right, m_i, m_k, m_j, result);
#endif
}

} // namespace gf_2_8
