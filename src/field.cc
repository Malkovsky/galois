#include "field.h"

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

/* GF(256) tables */
static element_t exp[256]; /* α^i */
static element_t log[256]; /* log_α(i) */

/* 8x8 Matrices for multiplying on particular element */
static uint64_t gfni_matrix[256];

void Init(void) {
  element_t x = 1;

  /* Generate exponential table and logarithm table */
  for (size_t i = 0; i < 255; i++) {
    exp[i] = x;
    log[x] = i;

    x = Multiply(x, primitive_element);
  }
  exp[255] = 1;
}

void InitGFNI(void) {
  for (int16_t y = 0; y < 256; ++y) {
    gfni_matrix[y] = 0;
    element_t row = y;
    for (size_t i = 0, shift = 0; i < 8; ++i, shift += 8) {
      gfni_matrix[y] |= ((uint64_t)row << shift);
      row = (row << 1) ^ ((row >> 7) * irreducible_poly);
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

element_t MultiplyGFNI(element_t a, element_t b) {
  return ((a & 1) * (gfni_matrix[b] & 255)) ^
         (((a >> 1) & 1) * ((gfni_matrix[b] >> 8) & 255)) ^
         (((a >> 2) & 1) * ((gfni_matrix[b] >> 16) & 255)) ^
         (((a >> 3) & 1) * ((gfni_matrix[b] >> 24) & 255)) ^
         (((a >> 4) & 1) * ((gfni_matrix[b] >> 32) & 255)) ^
         (((a >> 5) & 1) * ((gfni_matrix[b] >> 40) & 255)) ^
         (((a >> 6) & 1) * ((gfni_matrix[b] >> 48) & 255)) ^
         (((a >> 7) & 1) * ((gfni_matrix[b] >> 56) & 255));
}

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

element_t Div(element_t a, element_t b) { 
  return Multiply(a, Inv(b));
}

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
