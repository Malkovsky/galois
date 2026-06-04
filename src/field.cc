#include "field.h"

#include <cstring>
#include <immintrin.h>

namespace gf_2_8 {

/**
 * Irreducible polynomial defining the field, compatible with VAES/GFNI
 */
const element_t irreducible_poly = 0x1b; // x^4+x^3+x+1

/**
 * Primitive element α used in the operations. Note that for
 * irreducible polynomial 0x1b the minimum primitive element is
 * x + 1, i.e. `3` in binary.
 */
const element_t primitive_element = 3; // x

/**
 * Binary multiplication tables
 */
static element_t binary_table[256 * 256];

/* GF(256) tables */
static element_t exp[256]; /* α^i */
static element_t log[256]; /* log_α(i) */

/* 8x8 Matrices for multiplying on particular element */
static uint64_t gfni_matrix[256];

#if defined(__SSSE3__)
static __m128i simd_table_lo[256];
static __m128i simd_table_hi[256];
#if defined(__AVX2__)
__m256i simd256_table_lo[256];
__m256i simd256_table_hi[256];
#endif
#endif

void Init(void) {
  element_t x = 1;

  /* Generate exponential table and logarithm table */
  for (size_t i = 0; i < 255; i++) {
    exp[i] = x;
    log[x] = i;

    x = Multiply(x, primitive_element);
  }

  exp[255] = 1;

  for (size_t i = 0; i < 256; ++i) {
    for (size_t j = 0; j < 256; ++j) {
      binary_table[256 * i + j] = Multiply(i, j);
    }
  }

#if defined(__SSSE3__)
  for (size_t z = 0; z < 256; ++z) {
    alignas(16) element_t lo[16];
    alignas(16) element_t hi[16];
    const element_t *z_table = binary_table + 256 * z;
    for (size_t nibble = 0; nibble < 16; ++nibble) {
      lo[nibble] = z_table[nibble];
      hi[nibble] = z_table[nibble << 4];
    }

    simd_table_lo[z] =
        _mm_load_si128(reinterpret_cast<const __m128i *>(lo));
    simd_table_hi[z] =
        _mm_load_si128(reinterpret_cast<const __m128i *>(hi));
#if defined(__AVX2__)
    simd256_table_lo[z] = _mm256_broadcastsi128_si256(simd_table_lo[z]);
    simd256_table_hi[z] = _mm256_broadcastsi128_si256(simd_table_hi[z]);
#endif
  }
#endif
}

void InitGFNI(void) {
  for (int16_t y = 0; y < 256; ++y) {
    uint64_t mt = 0;
    element_t row = y;
    for (size_t i = 0, shift = 0; i < 8; ++i, shift += 8) {
      mt |= ((uint64_t)row << shift);
      row = (row << 1) ^ ((row >> 7) * irreducible_poly);
    }
    gfni_matrix[y] = 0;

    for (size_t i = 0; i < 8; ++i) {
      for (size_t j = 0; j < 8; ++j) {
        gfni_matrix[y] |= ((mt >> (8 * i + j)) & 1) << (8 * (7 - j) + i);
      }
    }
  }
}

element_t Zero() { return 0; }

element_t One() { return 1; }

element_t Add(element_t a, element_t b) {
  /* Addition in GF(256) is just XOR */
  return a ^ b;
}

element_t Sub(element_t a, element_t b) {
  /* Subtraction is the same as addition in GF(256) */
  return Add(a, b);
}

element_t MultiplyLUT(element_t a, element_t b) {
  /* Handle special case of multiplication by 0 */
  if (a == 0 || b == 0)
    return 0;
  auto p = log[a] + log[b];
  // These manipulations substracts 255 when p > 255
  return exp[(p & 255) + (p >> 8)];
}

element_t Multiply(element_t a, element_t b) {
  element_t result = 0;
  while (a) {
    result ^= b * (a & 1);
    a >>= 1;
    // b << 1 is polynomial multiplication by x
    // (b >> 7) checks whether b is polynomial of degree 7
    // ^ (irreducible_poly * (b >> 7)) does polynomial mod
    // that case
    b = (b << 1) ^ (irreducible_poly * (b >> 7));
  }
  return result;
}

// TODO: implement properly for testing/demo purpose
// element_t MultiplyGFNI(element_t a, element_t b) {
//     uint64_t x = 0x0101010101010101ULL * b & gfni_matrix[a];
//     x ^= x >> 4;
//     x ^= x >> 2;
//     x ^= x >> 1;
//     x &= 0x0101010101010101ULL;
//     return (uint8_t)((x * 0x0102040810204080ULL) >> 56);
// }

element_t Div(element_t a, element_t b) {
  if (a == 0) {
    return 0; /* 0 / b = 0 */
  }

  /* a / b = a * b^(-1) = α^(log(a) - log(b)) */
  int log_diff = log[a] - log[b];
  if (log_diff < 0) {
    log_diff += 255;
  }

  return exp[log_diff];
}

element_t Inv(element_t a) {
  if (a == 0) {
    return 0;
  }
  /* a^(-1) = α^(255 - log(a)) */
  return exp[255 - log[a]];
}

element_t Pow(element_t a, int n) {
  if (a == 0) {
    return (n == 0) ? 1 : 0; /* 0^0 = 1, 0^n = 0 for n > 0 */
  }

  /* If n is negative, we compute the inverse raised to |n| */
  if (n < 0) {
    a = Inv(a);
    n = -n;
  }

  /* a^n = α^(log(a) * n mod 255) */
  int log_res = (log[a] * n) % 255;
  return exp[log_res];
}

void AddScaledRowBase(element_t *x, const element_t *y, element_t z,
                      size_t length) {
  if (z == 0) {
    return;
  }
  element_t *z_table = binary_table + 256 * z;
  for (size_t i = 0; i < length; ++i) {
    x[i] ^= z_table[y[i]];
  }
}

void AddScaledRowSIMD(element_t *x, const element_t *y, element_t z,
                      size_t length) {
  if (z == 0) {
    return;
  }

  size_t processed = 0;
#if defined(__AVX2__)
  const __m256i mask = _mm256_set1_epi8(0x0f);
  const __m256i table_lo = simd256_table_lo[z];
  const __m256i table_hi = simd256_table_hi[z];
  while (processed + 32 <= length) {
    const __m256i y_reg =
        _mm256_loadu_si256(reinterpret_cast<const __m256i *>(y + processed));
    const __m256i lo = _mm256_shuffle_epi8(table_lo,
                                           _mm256_and_si256(y_reg, mask));
    const __m256i hi = _mm256_shuffle_epi8(
        table_hi, _mm256_and_si256(_mm256_srli_epi64(y_reg, 4), mask));
    const __m256i x_reg =
        _mm256_loadu_si256(reinterpret_cast<const __m256i *>(x + processed));
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(x + processed),
                        _mm256_xor_si256(x_reg, _mm256_xor_si256(lo, hi)));
    processed += 32;
  }
#elif defined(__SSSE3__)
  const __m128i mask = _mm_set1_epi8(0x0f);
  const __m128i table_lo = simd_table_lo[z];
  const __m128i table_hi = simd_table_hi[z];
  while (processed + 16 <= length) {
    const __m128i y_reg =
        _mm_loadu_si128(reinterpret_cast<const __m128i *>(y + processed));
    const __m128i lo =
        _mm_shuffle_epi8(table_lo, _mm_and_si128(y_reg, mask));
    const __m128i hi = _mm_shuffle_epi8(
        table_hi, _mm_and_si128(_mm_srli_epi64(y_reg, 4), mask));
    const __m128i x_reg =
        _mm_loadu_si128(reinterpret_cast<const __m128i *>(x + processed));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(x + processed),
                     _mm_xor_si128(x_reg, _mm_xor_si128(lo, hi)));
    processed += 16;
  }
#endif

  AddScaledRowBase(x + processed, y + processed, z, length - processed);
}

void AddScaledRowGFNIGeneral(element_t *x, const element_t *y, element_t z,
                             size_t length) {
  if (z == 0) {
    return;
  }
  size_t processed = 0;
#if defined(__GFNI__) and defined(__AVX512F__)
  __m512i z_matrix = _mm512_set1_epi64(gfni_matrix[z]);
  while (processed + 64 <= length) {
    auto x_reg = _mm512_loadu_epi8(x);
    auto y_reg = _mm512_loadu_epi8(y);
    x_reg = _mm512_xor_epi64(x_reg,
                             _mm512_gf2p8affine_epi64_epi8(y_reg, z_matrix, 0));
    _mm512_storeu_epi8(x, x_reg);
    x += 64;
    y += 64;
    processed += 64;
  }
#endif
  AddScaledRowBase(x, y, z, length - processed);
}

void AddScaledRowGFNIDedicated(element_t *x, const element_t *y, element_t z,
                               size_t length) {
  if (z == 0) {
    return;
  }
  size_t processed = 0;
#if defined(__GFNI__) and defined(__AVX512F__)
  __m512i z_reg = _mm512_set1_epi8(z);
  while (processed + 64 <= length) {
    auto x_reg = _mm512_loadu_epi8(x);
    auto y_reg = _mm512_loadu_epi8(y);
    x_reg = _mm512_xor_epi64(x_reg, _mm512_gf2p8mul_epi8(y_reg, z_reg));
    _mm512_storeu_epi8(x, x_reg);
    x += 64;
    y += 64;
    processed += 64;
  }
#endif
  AddScaledRowBase(x, y, z, length - processed);
}

} // namespace gf_2_8

namespace gf_2_16 {

gf_2_8::element_t delta = 0x20;

element_t Zero() { return 0; }

element_t One() { return 1; }

element_t Add(element_t a, element_t b) { return a ^ b; }

element_t Sub(element_t a, element_t b) { return a ^ b; }

element_t Multiply(element_t a, element_t b) {
  // a = a_0 + a_1x, b = b_0 + b_1x
  // all four are from GF(256)
  gf_2_8::element_t a_0 = a & 255;
  gf_2_8::element_t a_1 = a >> 8;
  gf_2_8::element_t b_0 = b & 255;
  gf_2_8::element_t b_1 = b >> 8;

  auto t = gf_2_8::MultiplyLUT(a_1, b_1);
  auto low_bits =
      gf_2_8::Add(gf_2_8::MultiplyLUT(a_0, b_0), gf_2_8::MultiplyLUT(t, delta));
  auto high_bits = gf_2_8::Add(
      gf_2_8::Add(gf_2_8::MultiplyLUT(a_0, b_1), gf_2_8::MultiplyLUT(a_1, b_0)),
      t);
  return low_bits + (high_bits << 8);
}

element_t Inv(element_t a) {
  element_t result = One();
  element_t b = Multiply(a, a);
  for (size_t i = 1; i < 16; ++i) {
    result = Multiply(result, b);
    b = Multiply(b, b);
  }
  return result;
}

element_t Div(element_t a, element_t b) { return Multiply(a, Inv(b)); }

element_t Pow(element_t a, size_t n) {
  element_t result = One();
  while (n) {
    if (n & 1) {
      result = Multiply(result, a);
    }
    a = Multiply(a, a);
    n >>= 1;
  }
  return result;
}

element_t InvIT(element_t a) {
  element_t a_r = a;
  for (size_t i = 0; i < 8; ++i) {
    a_r = Multiply(a_r, a_r);
  }
  gf_2_8::element_t a_r1 = Multiply(a_r, a);
  return Multiply(a_r, gf_2_8::Inv(a_r1));
}

} // namespace gf_2_16
