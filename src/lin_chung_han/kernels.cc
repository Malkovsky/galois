#include "lin_chung_han/kernels_internal.h"

// The radix-2/radix-4 butterfly scheduling follows Christopher A. Taylor's
// Leopard-RS kernels (BSD-3-Clause); this is an independent Cantor-coordinate,
// tail-safe implementation with explicit ISA variants.

#if defined(__i386__) || defined(__x86_64__) || defined(_M_IX86) || \
    defined(_M_X64)
#include <immintrin.h>
#endif

#include <algorithm>

namespace gf2p8::lch {
namespace detail {

Status GFNIMultiplyAdd(Element* x,
                       const Element* y,
                       size_t byte_count,
                       Element coefficient,
                       Backend backend,
                       const MultiplicationTables& tables);
Status GFNIMultiply(Element* destination,
                    const Element* source,
                    size_t byte_count,
                    Element coefficient,
                    Backend backend,
                    const MultiplicationTables& tables);
Status GFNIFFTRadix4(Element* x0,
                     Element* x1,
                     Element* x2,
                     Element* x3,
                     size_t byte_count,
                     Element top,
                     Element low,
                     Element high,
                     Backend backend,
                     const MultiplicationTables& tables);
Status GFNIIFFTRadix4(Element* x0,
                      Element* x1,
                      Element* x2,
                      Element* x3,
                      size_t byte_count,
                      Element top,
                      Element low,
                      Element high,
                      Backend backend,
                      const MultiplicationTables& tables);
Status GFNIIFFTRadix2Xor(const Element* x,
                         const Element* y,
                         Element* output_x,
                         Element* output_y,
                         size_t byte_count,
                         Element coefficient,
                         Backend backend,
                         const MultiplicationTables& tables);
Status GFNIIFFTRadix2Copy(const Element* x,
                          const Element* y,
                          Element* output_x,
                          Element* output_y,
                          size_t byte_count,
                          Element coefficient,
                          Backend backend,
                          const MultiplicationTables& tables);
Status GFNIIFFTRadix4Copy(const Element* x0,
                          const Element* x1,
                          const Element* x2,
                          const Element* x3,
                          Element* output0,
                          Element* output1,
                          Element* output2,
                          Element* output3,
                          size_t byte_count,
                          Element top,
                          Element low,
                          Element high,
                          Backend backend,
                          const MultiplicationTables& tables);
Status GFNIIFFTRadix4Xor(const Element* x0,
                         const Element* x1,
                         const Element* x2,
                         const Element* x3,
                         Element* output0,
                         Element* output1,
                         Element* output2,
                         Element* output3,
                         size_t byte_count,
                         Element top,
                         Element low,
                         Element high,
                         Backend backend,
                         const MultiplicationTables& tables);

}  // namespace detail
namespace {

bool IsGFNIBackend(Backend backend) {
  switch (backend) {
    case Backend::gfni128_affine:
    case Backend::gfni256_affine:
    case Backend::gfni512_affine:
      return true;
    default:
      return false;
  }
}

void ScalarMultiplyAdd(Element* x,
                       const Element* y,
                       size_t byte_count,
                       Element coefficient,
                       const MultiplicationTables& tables) {
  const auto& row = tables.shuffle[coefficient];
  for (size_t i = 0; i < byte_count; ++i) {
    x[i] ^= row[y[i] & 0x0f] ^ row[32 + (y[i] >> 4)];
  }
}

void ScalarMultiply(Element* destination,
                    const Element* source,
                    size_t byte_count,
                    Element coefficient,
                    const MultiplicationTables& tables) {
  const auto& row = tables.shuffle[coefficient];
  for (size_t i = 0; i < byte_count; ++i) {
    destination[i] = row[source[i] & 0x0f] ^ row[32 + (source[i] >> 4)];
  }
}

Element ScalarProduct(Element value,
                      Element coefficient,
                      const MultiplicationTables& tables) {
  const auto& row = tables.shuffle[coefficient];
  return row[value & 0x0f] ^ row[32 + (value >> 4)];
}

#if defined(__SSSE3__)
struct Factor128 {
  __m128i low;
  __m128i high;
};

Factor128 LoadFactor128(Element coefficient,
                        const MultiplicationTables& tables) {
  const Element* const row = tables.shuffle[coefficient].data();
  return {_mm_load_si128(reinterpret_cast<const __m128i*>(row)),
          _mm_load_si128(reinterpret_cast<const __m128i*>(row + 32))};
}

__m128i Product128(__m128i input, const Factor128& factor) {
  const __m128i mask = _mm_set1_epi8(0x0f);
  return _mm_xor_si128(
      _mm_shuffle_epi8(factor.low, _mm_and_si128(input, mask)),
      _mm_shuffle_epi8(factor.high,
                       _mm_and_si128(_mm_srli_epi64(input, 4), mask)));
}

__m128i AddProduct128(__m128i accumulator,
                      __m128i input,
                      const Factor128& factor,
                      Element coefficient) {
  if (coefficient == 0) {
    return accumulator;
  }
  return _mm_xor_si128(accumulator, Product128(input, factor));
}
#endif

#if defined(__SSE2__)
inline void Xor64Sse(Element* destination,
                     const Element* source,
                     size_t offset) {
  for (size_t lane = 0; lane < 4; ++lane) {
    const size_t lane_offset = offset + 16 * lane;
    const __m128i x = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(destination + lane_offset));
    const __m128i y =
        _mm_loadu_si128(reinterpret_cast<const __m128i*>(source + lane_offset));
    _mm_storeu_si128(reinterpret_cast<__m128i*>(destination + lane_offset),
                     _mm_xor_si128(x, y));
  }
}
#endif

#if defined(__AVX2__)
struct Factor256 {
  __m256i low;
  __m256i high;
};

Factor256 LoadFactor256(Element coefficient,
                        const MultiplicationTables& tables) {
  const Element* const row = tables.shuffle[coefficient].data();
  return {_mm256_load_si256(reinterpret_cast<const __m256i*>(row)),
          _mm256_load_si256(reinterpret_cast<const __m256i*>(row + 32))};
}

__m256i Product256(__m256i input, const Factor256& factor) {
  const __m256i mask = _mm256_set1_epi8(0x0f);
  return _mm256_xor_si256(
      _mm256_shuffle_epi8(factor.low, _mm256_and_si256(input, mask)),
      _mm256_shuffle_epi8(factor.high,
                          _mm256_and_si256(_mm256_srli_epi64(input, 4), mask)));
}

__m256i AddProduct256(__m256i accumulator,
                      __m256i input,
                      const Factor256& factor,
                      Element coefficient) {
  if (coefficient == 0) {
    return accumulator;
  }
  return _mm256_xor_si256(accumulator, Product256(input, factor));
}

inline void Xor128Avx2(Element* destination,
                       const Element* source,
                       size_t offset) {
  const __m256i x0 = _mm256_loadu_si256(
      reinterpret_cast<const __m256i*>(destination + offset));
  const __m256i x1 = _mm256_loadu_si256(
      reinterpret_cast<const __m256i*>(destination + offset + 32));
  const __m256i x2 = _mm256_loadu_si256(
      reinterpret_cast<const __m256i*>(destination + offset + 64));
  const __m256i x3 = _mm256_loadu_si256(
      reinterpret_cast<const __m256i*>(destination + offset + 96));
  const __m256i y0 =
      _mm256_loadu_si256(reinterpret_cast<const __m256i*>(source + offset));
  const __m256i y1 = _mm256_loadu_si256(
      reinterpret_cast<const __m256i*>(source + offset + 32));
  const __m256i y2 = _mm256_loadu_si256(
      reinterpret_cast<const __m256i*>(source + offset + 64));
  const __m256i y3 = _mm256_loadu_si256(
      reinterpret_cast<const __m256i*>(source + offset + 96));
  _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination + offset),
                      _mm256_xor_si256(x0, y0));
  _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination + offset + 32),
                      _mm256_xor_si256(x1, y1));
  _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination + offset + 64),
                      _mm256_xor_si256(x2, y2));
  _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination + offset + 96),
                      _mm256_xor_si256(x3, y3));
}
#endif

void FFTRadix4ScalarTail(Element* x0,
                         Element* x1,
                         Element* x2,
                         Element* x3,
                         size_t begin,
                         size_t byte_count,
                         Element top,
                         Element low,
                         Element high,
                         const MultiplicationTables& tables) {
  for (size_t i = begin; i < byte_count; ++i) {
    x0[i] ^= ScalarProduct(x2[i], top, tables);
    x2[i] ^= x0[i];
    x1[i] ^= ScalarProduct(x3[i], top, tables);
    x3[i] ^= x1[i];
    x0[i] ^= ScalarProduct(x1[i], low, tables);
    x1[i] ^= x0[i];
    x2[i] ^= ScalarProduct(x3[i], high, tables);
    x3[i] ^= x2[i];
  }
}

void IFFTRadix4ScalarTail(Element* x0,
                          Element* x1,
                          Element* x2,
                          Element* x3,
                          size_t begin,
                          size_t byte_count,
                          Element top,
                          Element low,
                          Element high,
                          const MultiplicationTables& tables) {
  for (size_t i = begin; i < byte_count; ++i) {
    x1[i] ^= x0[i];
    x0[i] ^= ScalarProduct(x1[i], low, tables);
    x3[i] ^= x2[i];
    x2[i] ^= ScalarProduct(x3[i], high, tables);
    x2[i] ^= x0[i];
    x0[i] ^= ScalarProduct(x2[i], top, tables);
    x3[i] ^= x1[i];
    x1[i] ^= ScalarProduct(x3[i], top, tables);
  }
}

void IFFTRadix4XorScalarTail(const Element* x0,
                             const Element* x1,
                             const Element* x2,
                             const Element* x3,
                             Element* output0,
                             Element* output1,
                             Element* output2,
                             Element* output3,
                             size_t begin,
                             size_t byte_count,
                             Element top,
                             Element low,
                             Element high,
                             const MultiplicationTables& tables) {
  for (size_t i = begin; i < byte_count; ++i) {
    Element a = x0[i];
    Element b = x1[i] ^ a;
    a ^= ScalarProduct(b, low, tables);
    Element c = x2[i];
    Element d = x3[i] ^ c;
    c ^= ScalarProduct(d, high, tables);
    c ^= a;
    a ^= ScalarProduct(c, top, tables);
    d ^= b;
    b ^= ScalarProduct(d, top, tables);
    output0[i] ^= a;
    output1[i] ^= b;
    output2[i] ^= c;
    output3[i] ^= d;
  }
}

#if defined(__AVX2__)
void ScaleAvx2Impl(Element* destination,
                   const Element* source,
                   size_t byte_count,
                   Element coefficient,
                   const MultiplicationTables& tables) {
  const Factor256 factor = LoadFactor256(coefficient, tables);
  size_t processed = 0;
  if ((byte_count & 63U) == 0) {
    for (; processed < byte_count; processed += 64) {
      const __m256i input0 = _mm256_loadu_si256(
          reinterpret_cast<const __m256i*>(source + processed));
      const __m256i input1 = _mm256_loadu_si256(
          reinterpret_cast<const __m256i*>(source + processed + 32));
      _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination + processed),
                          Product256(input0, factor));
      _mm256_storeu_si256(
          reinterpret_cast<__m256i*>(destination + processed + 32),
          Product256(input1, factor));
    }
    return;
  }
  for (; processed + 64 <= byte_count; processed += 64) {
    const __m256i input0 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(source + processed));
    const __m256i input1 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(source + processed + 32));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination + processed),
                        Product256(input0, factor));
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(destination + processed + 32),
        Product256(input1, factor));
  }
  for (; processed + 32 <= byte_count; processed += 32) {
    const __m256i input = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(source + processed));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination + processed),
                        Product256(input, factor));
  }
  ScalarMultiply(destination + processed, source + processed,
                 byte_count - processed, coefficient, tables);
}

void XorAvx2Impl(Element* destination,
                 const Element* source,
                 size_t byte_count) {
  size_t processed = 0;
  if ((byte_count & 127U) == 0) {
    for (; processed < byte_count; processed += 128) {
      Xor128Avx2(destination, source, processed);
    }
    return;
  }
  for (; processed + 128 <= byte_count; processed += 128) {
    Xor128Avx2(destination, source, processed);
  }
  for (; processed + 32 <= byte_count; processed += 32) {
    const __m256i x = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(destination + processed));
    const __m256i y = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(source + processed));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination + processed),
                        _mm256_xor_si256(x, y));
  }
  for (; processed < byte_count; ++processed) {
    destination[processed] ^= source[processed];
  }
}

void Xor4Avx2Impl(Element* destination0,
                  const Element* source0,
                  Element* destination1,
                  const Element* source1,
                  Element* destination2,
                  const Element* source2,
                  Element* destination3,
                  const Element* source3,
                  size_t byte_count) {
  size_t processed = 0;
  if ((byte_count & 127U) == 0) {
    for (; processed < byte_count; processed += 128) {
      Xor128Avx2(destination0, source0, processed);
      Xor128Avx2(destination1, source1, processed);
      Xor128Avx2(destination2, source2, processed);
      Xor128Avx2(destination3, source3, processed);
    }
    return;
  }
  for (; processed + 128 <= byte_count; processed += 128) {
    Xor128Avx2(destination0, source0, processed);
    Xor128Avx2(destination1, source1, processed);
    Xor128Avx2(destination2, source2, processed);
    Xor128Avx2(destination3, source3, processed);
  }
  for (; processed + 32 <= byte_count; processed += 32) {
    const __m256i x0 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(destination0 + processed));
    const __m256i x1 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(destination1 + processed));
    const __m256i x2 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(destination2 + processed));
    const __m256i x3 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(destination3 + processed));
    const __m256i y0 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(source0 + processed));
    const __m256i y1 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(source1 + processed));
    const __m256i y2 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(source2 + processed));
    const __m256i y3 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(source3 + processed));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination0 + processed),
                        _mm256_xor_si256(x0, y0));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination1 + processed),
                        _mm256_xor_si256(x1, y1));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination2 + processed),
                        _mm256_xor_si256(x2, y2));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination3 + processed),
                        _mm256_xor_si256(x3, y3));
  }
  for (; processed < byte_count; ++processed) {
    destination0[processed] ^= source0[processed];
    destination1[processed] ^= source1[processed];
    destination2[processed] ^= source2[processed];
    destination3[processed] ^= source3[processed];
  }
}

template <bool ExactVectors = false>
void FFTRadix2Avx2Impl(Element* x,
                       Element* y,
                       size_t byte_count,
                       Element coefficient,
                       const MultiplicationTables& tables) {
  const Factor256 factor = LoadFactor256(coefficient, tables);
  size_t processed = 0;
  for (; processed + 32 <= byte_count; processed += 32) {
    __m256i a =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x + processed));
    __m256i b =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(y + processed));
    a = AddProduct256(a, b, factor, coefficient);
    b = _mm256_xor_si256(b, a);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(x + processed), a);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(y + processed), b);
  }
  if constexpr (!ExactVectors) {
    for (; processed < byte_count; ++processed) {
      x[processed] ^= ScalarProduct(y[processed], coefficient, tables);
      y[processed] ^= x[processed];
    }
  }
}

template <bool ExactVectors = false>
void IFFTRadix2Avx2Impl(Element* x,
                        Element* y,
                        size_t byte_count,
                        Element coefficient,
                        const MultiplicationTables& tables) {
  const Factor256 factor = LoadFactor256(coefficient, tables);
  size_t processed = 0;
  for (; processed + 32 <= byte_count; processed += 32) {
    __m256i a =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x + processed));
    __m256i b =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(y + processed));
    b = _mm256_xor_si256(b, a);
    a = AddProduct256(a, b, factor, coefficient);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(x + processed), a);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(y + processed), b);
  }
  if constexpr (!ExactVectors) {
    for (; processed < byte_count; ++processed) {
      y[processed] ^= x[processed];
      x[processed] ^= ScalarProduct(y[processed], coefficient, tables);
    }
  }
}

template <bool ExactVectors = false>
void IFFTRadix2XorAvx2Impl(const Element* x,
                           const Element* y,
                           Element* output_x,
                           Element* output_y,
                           size_t byte_count,
                           Element coefficient,
                           const MultiplicationTables& tables) {
  const Factor256 factor = LoadFactor256(coefficient, tables);
  size_t processed = 0;
  for (; processed + 32 <= byte_count; processed += 32) {
    __m256i a =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x + processed));
    __m256i b =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(y + processed));
    b = _mm256_xor_si256(b, a);
    a = AddProduct256(a, b, factor, coefficient);
    const __m256i ox = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(output_x + processed));
    const __m256i oy = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(output_y + processed));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(output_x + processed),
                        _mm256_xor_si256(ox, a));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(output_y + processed),
                        _mm256_xor_si256(oy, b));
  }
  if constexpr (!ExactVectors) {
    for (; processed < byte_count; ++processed) {
      const Element b = y[processed] ^ x[processed];
      const Element a = x[processed] ^ ScalarProduct(b, coefficient, tables);
      output_x[processed] ^= a;
      output_y[processed] ^= b;
    }
  }
}

template <bool ExactVectors = false>
void IFFTRadix2CopyAvx2Impl(const Element* x,
                            const Element* y,
                            Element* output_x,
                            Element* output_y,
                            size_t byte_count,
                            Element coefficient,
                            const MultiplicationTables& tables) {
  const Factor256 factor = LoadFactor256(coefficient, tables);
  size_t processed = 0;
  for (; processed + 32 <= byte_count; processed += 32) {
    __m256i a =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x + processed));
    __m256i b =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(y + processed));
    b = _mm256_xor_si256(b, a);
    a = AddProduct256(a, b, factor, coefficient);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(output_x + processed), a);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(output_y + processed), b);
  }
  if constexpr (!ExactVectors) {
    for (; processed < byte_count; ++processed) {
      const Element b = y[processed] ^ x[processed];
      output_x[processed] =
          x[processed] ^ ScalarProduct(b, coefficient, tables);
      output_y[processed] = b;
    }
  }
}

template <bool ExactVectors = false>
void FFTRadix4Avx2Impl(Element* x0,
                       Element* x1,
                       Element* x2,
                       Element* x3,
                       size_t byte_count,
                       Element top,
                       Element low,
                       Element high,
                       const MultiplicationTables& tables) {
  const Factor256 top_factor = LoadFactor256(top, tables);
  const Factor256 low_factor = LoadFactor256(low, tables);
  const Factor256 high_factor = LoadFactor256(high, tables);
  size_t processed = 0;
  for (; processed + 32 <= byte_count; processed += 32) {
    __m256i a =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x0 + processed));
    __m256i b =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x1 + processed));
    __m256i c =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x2 + processed));
    __m256i d =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x3 + processed));
    if (top != 0) {
      a = _mm256_xor_si256(a, Product256(c, top_factor));
      b = _mm256_xor_si256(b, Product256(d, top_factor));
    }
    c = _mm256_xor_si256(c, a);
    d = _mm256_xor_si256(d, b);
    if (low != 0) {
      a = _mm256_xor_si256(a, Product256(b, low_factor));
    }
    b = _mm256_xor_si256(b, a);
    if (high != 0) {
      c = _mm256_xor_si256(c, Product256(d, high_factor));
    }
    d = _mm256_xor_si256(d, c);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(x0 + processed), a);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(x1 + processed), b);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(x2 + processed), c);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(x3 + processed), d);
  }
  if constexpr (!ExactVectors) {
    FFTRadix4ScalarTail(x0, x1, x2, x3, processed, byte_count, top, low, high,
                        tables);
  }
}

template <bool ExactVectors = false>
void IFFTRadix4Avx2Impl(Element* x0,
                        Element* x1,
                        Element* x2,
                        Element* x3,
                        size_t byte_count,
                        Element top,
                        Element low,
                        Element high,
                        const MultiplicationTables& tables) {
  const Factor256 top_factor = LoadFactor256(top, tables);
  const Factor256 low_factor = LoadFactor256(low, tables);
  const Factor256 high_factor = LoadFactor256(high, tables);
  size_t processed = 0;
  for (; processed + 32 <= byte_count; processed += 32) {
    __m256i a =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x0 + processed));
    __m256i b =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x1 + processed));
    __m256i c =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x2 + processed));
    __m256i d =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x3 + processed));
    b = _mm256_xor_si256(b, a);
    if (low != 0) {
      a = _mm256_xor_si256(a, Product256(b, low_factor));
    }
    d = _mm256_xor_si256(d, c);
    if (high != 0) {
      c = _mm256_xor_si256(c, Product256(d, high_factor));
    }
    c = _mm256_xor_si256(c, a);
    d = _mm256_xor_si256(d, b);
    if (top != 0) {
      a = _mm256_xor_si256(a, Product256(c, top_factor));
      b = _mm256_xor_si256(b, Product256(d, top_factor));
    }
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(x0 + processed), a);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(x1 + processed), b);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(x2 + processed), c);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(x3 + processed), d);
  }
  if constexpr (!ExactVectors) {
    IFFTRadix4ScalarTail(x0, x1, x2, x3, processed, byte_count, top, low, high,
                         tables);
  }
}

template <bool ExactVectors = false>
void IFFTRadix4XorAvx2Impl(const Element* x0,
                           const Element* x1,
                           const Element* x2,
                           const Element* x3,
                           Element* output0,
                           Element* output1,
                           Element* output2,
                           Element* output3,
                           size_t byte_count,
                           Element top,
                           Element low,
                           Element high,
                           const MultiplicationTables& tables) {
  const Factor256 top_factor = LoadFactor256(top, tables);
  const Factor256 low_factor = LoadFactor256(low, tables);
  const Factor256 high_factor = LoadFactor256(high, tables);
  size_t processed = 0;
  for (; processed + 32 <= byte_count; processed += 32) {
    __m256i a =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x0 + processed));
    __m256i b =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x1 + processed));
    __m256i c =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x2 + processed));
    __m256i d =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x3 + processed));
    b = _mm256_xor_si256(b, a);
    if (low != 0) {
      a = _mm256_xor_si256(a, Product256(b, low_factor));
    }
    d = _mm256_xor_si256(d, c);
    if (high != 0) {
      c = _mm256_xor_si256(c, Product256(d, high_factor));
    }
    c = _mm256_xor_si256(c, a);
    d = _mm256_xor_si256(d, b);
    if (top != 0) {
      a = _mm256_xor_si256(a, Product256(c, top_factor));
      b = _mm256_xor_si256(b, Product256(d, top_factor));
    }
    const __m256i o0 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(output0 + processed));
    const __m256i o1 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(output1 + processed));
    const __m256i o2 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(output2 + processed));
    const __m256i o3 = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(output3 + processed));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(output0 + processed),
                        _mm256_xor_si256(o0, a));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(output1 + processed),
                        _mm256_xor_si256(o1, b));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(output2 + processed),
                        _mm256_xor_si256(o2, c));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(output3 + processed),
                        _mm256_xor_si256(o3, d));
  }
  if constexpr (!ExactVectors) {
    IFFTRadix4XorScalarTail(x0, x1, x2, x3, output0, output1, output2, output3,
                            processed, byte_count, top, low, high, tables);
  }
}

template <bool ExactVectors = false>
void IFFTRadix4CopyAvx2Impl(const Element* x0,
                            const Element* x1,
                            const Element* x2,
                            const Element* x3,
                            Element* output0,
                            Element* output1,
                            Element* output2,
                            Element* output3,
                            size_t byte_count,
                            Element top,
                            Element low,
                            Element high,
                            const MultiplicationTables& tables) {
  const Factor256 top_factor = LoadFactor256(top, tables);
  const Factor256 low_factor = LoadFactor256(low, tables);
  const Factor256 high_factor = LoadFactor256(high, tables);
  size_t processed = 0;
  for (; processed + 32 <= byte_count; processed += 32) {
    __m256i a =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x0 + processed));
    __m256i b =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x1 + processed));
    __m256i c =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x2 + processed));
    __m256i d =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x3 + processed));
    b = _mm256_xor_si256(b, a);
    if (low != 0) {
      a = _mm256_xor_si256(a, Product256(b, low_factor));
    }
    d = _mm256_xor_si256(d, c);
    if (high != 0) {
      c = _mm256_xor_si256(c, Product256(d, high_factor));
    }
    c = _mm256_xor_si256(c, a);
    d = _mm256_xor_si256(d, b);
    if (top != 0) {
      a = _mm256_xor_si256(a, Product256(c, top_factor));
      b = _mm256_xor_si256(b, Product256(d, top_factor));
    }
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(output0 + processed), a);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(output1 + processed), b);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(output2 + processed), c);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(output3 + processed), d);
  }
  if constexpr (!ExactVectors) {
    for (; processed < byte_count; ++processed) {
      Element a = x0[processed];
      Element b = x1[processed] ^ a;
      a ^= ScalarProduct(b, low, tables);
      Element c = x2[processed];
      Element d = x3[processed] ^ c;
      c ^= ScalarProduct(d, high, tables);
      c ^= a;
      d ^= b;
      a ^= ScalarProduct(c, top, tables);
      b ^= ScalarProduct(d, top, tables);
      output0[processed] = a;
      output1[processed] = b;
      output2[processed] = c;
      output3[processed] = d;
    }
  }
}
#endif

Status AddScaledImpl(Element* x,
                     const Element* y,
                     size_t byte_count,
                     Element coefficient,
                     Backend backend,
                     const MultiplicationTables& tables) {
  if (coefficient == 0 || byte_count == 0) {
    return Status::ok;
  }
  if (IsGFNIBackend(backend)) {
    return detail::GFNIMultiplyAdd(x, y, byte_count, coefficient, backend,
                                   tables);
  }

  size_t processed = 0;
#if defined(__AVX2__)
  if (backend == Backend::avx2) {
    const __m256i mask = _mm256_set1_epi8(0x0f);
    const Element* const row = tables.shuffle[coefficient].data();
    const __m256i low =
        _mm256_load_si256(reinterpret_cast<const __m256i*>(row));
    const __m256i high =
        _mm256_load_si256(reinterpret_cast<const __m256i*>(row + 32));
    for (; processed + 32 <= byte_count; processed += 32) {
      const __m256i input =
          _mm256_loadu_si256(reinterpret_cast<const __m256i*>(y + processed));
      const __m256i product = _mm256_xor_si256(
          _mm256_shuffle_epi8(low, _mm256_and_si256(input, mask)),
          _mm256_shuffle_epi8(
              high, _mm256_and_si256(_mm256_srli_epi64(input, 4), mask)));
      const __m256i output =
          _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x + processed));
      _mm256_storeu_si256(reinterpret_cast<__m256i*>(x + processed),
                          _mm256_xor_si256(output, product));
    }
    ScalarMultiplyAdd(x + processed, y + processed, byte_count - processed,
                      coefficient, tables);
    return Status::ok;
  }
#endif
#if defined(__SSSE3__)
  if (backend == Backend::ssse3) {
    const __m128i mask = _mm_set1_epi8(0x0f);
    const Element* const row = tables.shuffle[coefficient].data();
    const __m128i low = _mm_load_si128(reinterpret_cast<const __m128i*>(row));
    const __m128i high =
        _mm_load_si128(reinterpret_cast<const __m128i*>(row + 32));
    for (; processed + 16 <= byte_count; processed += 16) {
      const __m128i input =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(y + processed));
      const __m128i product = _mm_xor_si128(
          _mm_shuffle_epi8(low, _mm_and_si128(input, mask)),
          _mm_shuffle_epi8(high,
                           _mm_and_si128(_mm_srli_epi64(input, 4), mask)));
      const __m128i output =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(x + processed));
      _mm_storeu_si128(reinterpret_cast<__m128i*>(x + processed),
                       _mm_xor_si128(output, product));
    }
    ScalarMultiplyAdd(x + processed, y + processed, byte_count - processed,
                      coefficient, tables);
    return Status::ok;
  }
#endif
  ScalarMultiplyAdd(x + processed, y + processed, byte_count - processed,
                    coefficient, tables);
  return Status::ok;
}

Status MultiplyImpl(Element* destination,
                    const Element* source,
                    size_t byte_count,
                    Element coefficient,
                    Backend backend,
                    const MultiplicationTables& tables) {
  if (coefficient == 0) {
    std::fill_n(destination, byte_count, 0);
    return Status::ok;
  }
  if (coefficient == 1) {
    std::copy_n(source, byte_count, destination);
    return Status::ok;
  }
  if (byte_count == 0) {
    return Status::ok;
  }
  if (IsGFNIBackend(backend)) {
    return detail::GFNIMultiply(destination, source, byte_count, coefficient,
                                backend, tables);
  }

  size_t processed = 0;
#if defined(__AVX2__)
  if (backend == Backend::avx2) {
    ScaleAvx2Impl(destination, source, byte_count, coefficient, tables);
    return Status::ok;
  }
#endif
#if defined(__SSSE3__)
  if (backend == Backend::ssse3) {
    const Factor128 factor = LoadFactor128(coefficient, tables);
    for (; processed + 64 <= byte_count; processed += 64) {
      for (size_t lane = 0; lane < 4; ++lane) {
        const __m128i input = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(source + processed + 16 * lane));
        _mm_storeu_si128(
            reinterpret_cast<__m128i*>(destination + processed + 16 * lane),
            Product128(input, factor));
      }
    }
    for (; processed + 16 <= byte_count; processed += 16) {
      const __m128i input =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(source + processed));
      _mm_storeu_si128(reinterpret_cast<__m128i*>(destination + processed),
                       Product128(input, factor));
    }
    ScalarMultiply(destination + processed, source + processed,
                   byte_count - processed, coefficient, tables);
    return Status::ok;
  }
#endif
  ScalarMultiply(destination + processed, source + processed,
                 byte_count - processed, coefficient, tables);
  return Status::ok;
}

void XorImpl(Element* destination,
             const Element* source,
             size_t byte_count,
             Backend backend) {
  size_t processed = 0;
#if defined(__AVX512F__) && defined(__AVX512BW__)
  if (backend == Backend::gfni512_affine) {
    for (; processed + 64 <= byte_count; processed += 64) {
      const __m512i x = _mm512_loadu_si512(destination + processed);
      const __m512i y = _mm512_loadu_si512(source + processed);
      _mm512_storeu_si512(destination + processed, _mm512_xor_si512(x, y));
    }
    for (; processed < byte_count; ++processed) {
      destination[processed] ^= source[processed];
    }
    return;
  }
#endif
#if defined(__AVX2__)
  if (backend == Backend::avx2 || backend == Backend::gfni256_affine) {
    XorAvx2Impl(destination, source, byte_count);
    return;
  }
#endif
#if defined(__SSE2__)
  if (backend == Backend::ssse3 || backend == Backend::gfni128_affine) {
    for (; processed + 64 <= byte_count; processed += 64) {
      Xor64Sse(destination, source, processed);
    }
    for (; processed + 16 <= byte_count; processed += 16) {
      const __m128i x = _mm_loadu_si128(
          reinterpret_cast<const __m128i*>(destination + processed));
      const __m128i y =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(source + processed));
      _mm_storeu_si128(reinterpret_cast<__m128i*>(destination + processed),
                       _mm_xor_si128(x, y));
    }
    for (; processed < byte_count; ++processed) {
      destination[processed] ^= source[processed];
    }
    return;
  }
#endif
  for (; processed < byte_count; ++processed) {
    destination[processed] ^= source[processed];
  }
}

void Xor4Impl(Element* destination0,
              const Element* source0,
              Element* destination1,
              const Element* source1,
              Element* destination2,
              const Element* source2,
              Element* destination3,
              const Element* source3,
              size_t byte_count,
              Backend backend) {
  size_t processed = 0;
#if defined(__AVX512F__) && defined(__AVX512BW__)
  if (backend == Backend::gfni512_affine) {
    for (; processed + 64 <= byte_count; processed += 64) {
      const __m512i x0 = _mm512_loadu_si512(destination0 + processed);
      const __m512i x1 = _mm512_loadu_si512(destination1 + processed);
      const __m512i x2 = _mm512_loadu_si512(destination2 + processed);
      const __m512i x3 = _mm512_loadu_si512(destination3 + processed);
      _mm512_storeu_si512(
          destination0 + processed,
          _mm512_xor_si512(x0, _mm512_loadu_si512(source0 + processed)));
      _mm512_storeu_si512(
          destination1 + processed,
          _mm512_xor_si512(x1, _mm512_loadu_si512(source1 + processed)));
      _mm512_storeu_si512(
          destination2 + processed,
          _mm512_xor_si512(x2, _mm512_loadu_si512(source2 + processed)));
      _mm512_storeu_si512(
          destination3 + processed,
          _mm512_xor_si512(x3, _mm512_loadu_si512(source3 + processed)));
    }
    for (; processed < byte_count; ++processed) {
      destination0[processed] ^= source0[processed];
      destination1[processed] ^= source1[processed];
      destination2[processed] ^= source2[processed];
      destination3[processed] ^= source3[processed];
    }
    return;
  }
#endif
#if defined(__AVX2__)
  if (backend == Backend::avx2 || backend == Backend::gfni256_affine) {
    Xor4Avx2Impl(destination0, source0, destination1, source1, destination2,
                 source2, destination3, source3, byte_count);
    return;
  }
#endif
#if defined(__SSE2__)
  if (backend == Backend::ssse3 || backend == Backend::gfni128_affine) {
    for (; processed + 64 <= byte_count; processed += 64) {
      Xor64Sse(destination0, source0, processed);
      Xor64Sse(destination1, source1, processed);
      Xor64Sse(destination2, source2, processed);
      Xor64Sse(destination3, source3, processed);
    }
    for (; processed + 16 <= byte_count; processed += 16) {
      const __m128i x0 = _mm_loadu_si128(
          reinterpret_cast<const __m128i*>(destination0 + processed));
      const __m128i x1 = _mm_loadu_si128(
          reinterpret_cast<const __m128i*>(destination1 + processed));
      const __m128i x2 = _mm_loadu_si128(
          reinterpret_cast<const __m128i*>(destination2 + processed));
      const __m128i x3 = _mm_loadu_si128(
          reinterpret_cast<const __m128i*>(destination3 + processed));
      _mm_storeu_si128(
          reinterpret_cast<__m128i*>(destination0 + processed),
          _mm_xor_si128(x0, _mm_loadu_si128(reinterpret_cast<const __m128i*>(
                                source0 + processed))));
      _mm_storeu_si128(
          reinterpret_cast<__m128i*>(destination1 + processed),
          _mm_xor_si128(x1, _mm_loadu_si128(reinterpret_cast<const __m128i*>(
                                source1 + processed))));
      _mm_storeu_si128(
          reinterpret_cast<__m128i*>(destination2 + processed),
          _mm_xor_si128(x2, _mm_loadu_si128(reinterpret_cast<const __m128i*>(
                                source2 + processed))));
      _mm_storeu_si128(
          reinterpret_cast<__m128i*>(destination3 + processed),
          _mm_xor_si128(x3, _mm_loadu_si128(reinterpret_cast<const __m128i*>(
                                source3 + processed))));
    }
    for (; processed < byte_count; ++processed) {
      destination0[processed] ^= source0[processed];
      destination1[processed] ^= source1[processed];
      destination2[processed] ^= source2[processed];
      destination3[processed] ^= source3[processed];
    }
    return;
  }
#endif
  for (; processed < byte_count; ++processed) {
    destination0[processed] ^= source0[processed];
    destination1[processed] ^= source1[processed];
    destination2[processed] ^= source2[processed];
    destination3[processed] ^= source3[processed];
  }
}

}  // namespace

bool BackendAvailable(Backend backend) {
  switch (backend) {
    case Backend::tuned:
    case Backend::scalar:
      return true;
    case Backend::ssse3:
#if defined(__SSSE3__)
      return true;
#else
      return false;
#endif
    case Backend::avx2:
#if defined(__AVX2__)
      return true;
#else
      return false;
#endif
    case Backend::gfni128_affine:
#if defined(__GFNI__) && defined(__SSE2__)
      return true;
#else
      return false;
#endif
    case Backend::gfni256_affine:
#if defined(__GFNI__) && defined(__AVX2__)
      return true;
#else
      return false;
#endif
    case Backend::gfni512_affine:
#if defined(__GFNI__) && defined(__AVX512F__) && defined(__AVX512BW__)
      return true;
#else
      return false;
#endif
  }
  return false;
}

Backend SelectBackend(size_t byte_count) {
  static_cast<void>(byte_count);
#if defined(__GFNI__) && defined(__AVX2__)
  if (byte_count >= 128) {
    return Backend::gfni256_affine;
  }
#endif
#if defined(__AVX2__)
  if (byte_count >= 32) {
    return Backend::avx2;
  }
#endif
#if defined(__SSSE3__)
  if (byte_count >= 16) {
    return Backend::ssse3;
  }
#endif
  return Backend::scalar;
}

Status AddScaled(Element* destination,
                 const Element* source,
                 size_t byte_count,
                 Element coefficient,
                 Backend backend,
                 const MultiplicationTables& tables) {
  if (backend == Backend::tuned) {
    backend = SelectBackend(byte_count);
  }
  if (!BackendAvailable(backend)) {
    return Status::unsupported_backend;
  }
  return AddScaledImpl(destination, source, byte_count, coefficient, backend,
                       tables);
}

Status Scale(Element* destination,
             const Element* source,
             size_t byte_count,
             Element coefficient,
             Backend backend,
             const MultiplicationTables& tables) {
  if (backend == Backend::tuned) {
    backend = SelectBackend(byte_count);
  }
  if (!BackendAvailable(backend)) {
    return Status::unsupported_backend;
  }
  return detail::ScaleUnchecked(destination, source, byte_count, coefficient,
                                backend, tables);
}

Status Xor(Element* destination,
           const Element* source,
           size_t byte_count,
           Backend backend) {
  if (backend == Backend::tuned) {
    backend = SelectBackend(byte_count);
  }
  if (!BackendAvailable(backend)) {
    return Status::unsupported_backend;
  }
  detail::XorUnchecked(destination, source, byte_count, backend);
  return Status::ok;
}

Status Xor4(Element* destination0,
            const Element* source0,
            Element* destination1,
            const Element* source1,
            Element* destination2,
            const Element* source2,
            Element* destination3,
            const Element* source3,
            size_t byte_count,
            Backend backend) {
  if (backend == Backend::tuned) {
    backend = SelectBackend(byte_count);
  }
  if (!BackendAvailable(backend)) {
    return Status::unsupported_backend;
  }
  detail::Xor4Unchecked(destination0, source0, destination1, source1,
                        destination2, source2, destination3, source3,
                        byte_count, backend);
  return Status::ok;
}

namespace detail {

Status ScaleUnchecked(Element* destination,
                      const Element* source,
                      size_t byte_count,
                      Element coefficient,
                      Backend backend,
                      const MultiplicationTables& tables) {
#if defined(__AVX2__)
  if (backend == Backend::avx2) {
    ScaleAvx2Impl(destination, source, byte_count, coefficient, tables);
    return Status::ok;
  }
#endif
  return MultiplyImpl(destination, source, byte_count, coefficient, backend,
                      tables);
}

void XorUnchecked(Element* destination,
                  const Element* source,
                  size_t byte_count,
                  Backend backend) {
#if defined(__AVX2__)
  if (backend == Backend::avx2) {
    XorAvx2Impl(destination, source, byte_count);
    return;
  }
#endif
  XorImpl(destination, source, byte_count, backend);
}

void Xor4Unchecked(Element* destination0,
                   const Element* source0,
                   Element* destination1,
                   const Element* source1,
                   Element* destination2,
                   const Element* source2,
                   Element* destination3,
                   const Element* source3,
                   size_t byte_count,
                   Backend backend) {
#if defined(__AVX2__)
  if (backend == Backend::avx2) {
    Xor4Avx2Impl(destination0, source0, destination1, source1, destination2,
                 source2, destination3, source3, byte_count);
    return;
  }
#endif
  Xor4Impl(destination0, source0, destination1, source1, destination2, source2,
           destination3, source3, byte_count, backend);
}

}  // namespace detail

void detail::FFTRadix2Unchecked(Element* x,
                                Element* y,
                                size_t byte_count,
                                Element coefficient,
                                Backend backend,
                                const MultiplicationTables& tables) {
  if (IsGFNIBackend(backend)) {
    static_cast<void>(
        AddScaledImpl(x, y, byte_count, coefficient, backend, tables));
    XorImpl(y, x, byte_count, backend);
    return;
  }

  size_t processed = 0;
#if defined(__AVX2__)
  if (backend == Backend::avx2) {
    FFTRadix2Avx2Impl(x, y, byte_count, coefficient, tables);
    return;
  }
#endif
#if defined(__SSSE3__)
  if (backend == Backend::ssse3) {
    const Factor128 factor = LoadFactor128(coefficient, tables);
    for (; processed + 16 <= byte_count; processed += 16) {
      __m128i a =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(x + processed));
      __m128i b =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(y + processed));
      a = AddProduct128(a, b, factor, coefficient);
      b = _mm_xor_si128(b, a);
      _mm_storeu_si128(reinterpret_cast<__m128i*>(x + processed), a);
      _mm_storeu_si128(reinterpret_cast<__m128i*>(y + processed), b);
    }
    for (; processed < byte_count; ++processed) {
      x[processed] ^= ScalarProduct(y[processed], coefficient, tables);
      y[processed] ^= x[processed];
    }
    return;
  }
#endif
  for (; processed < byte_count; ++processed) {
    x[processed] ^= ScalarProduct(y[processed], coefficient, tables);
    y[processed] ^= x[processed];
  }
}

void detail::IFFTRadix2Unchecked(Element* x,
                                 Element* y,
                                 size_t byte_count,
                                 Element coefficient,
                                 Backend backend,
                                 const MultiplicationTables& tables) {
  if (IsGFNIBackend(backend)) {
    XorImpl(y, x, byte_count, backend);
    static_cast<void>(
        AddScaledImpl(x, y, byte_count, coefficient, backend, tables));
    return;
  }

  size_t processed = 0;
#if defined(__AVX2__)
  if (backend == Backend::avx2) {
    IFFTRadix2Avx2Impl(x, y, byte_count, coefficient, tables);
    return;
  }
#endif
#if defined(__SSSE3__)
  if (backend == Backend::ssse3) {
    const Factor128 factor = LoadFactor128(coefficient, tables);
    for (; processed + 16 <= byte_count; processed += 16) {
      __m128i a =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(x + processed));
      __m128i b =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(y + processed));
      b = _mm_xor_si128(b, a);
      a = AddProduct128(a, b, factor, coefficient);
      _mm_storeu_si128(reinterpret_cast<__m128i*>(x + processed), a);
      _mm_storeu_si128(reinterpret_cast<__m128i*>(y + processed), b);
    }
    for (; processed < byte_count; ++processed) {
      y[processed] ^= x[processed];
      x[processed] ^= ScalarProduct(y[processed], coefficient, tables);
    }
    return;
  }
#endif
  for (; processed < byte_count; ++processed) {
    y[processed] ^= x[processed];
    x[processed] ^= ScalarProduct(y[processed], coefficient, tables);
  }
}

void detail::IFFTRadix2XorUnchecked(const Element* x,
                                    const Element* y,
                                    Element* output_x,
                                    Element* output_y,
                                    size_t byte_count,
                                    Element coefficient,
                                    Backend backend,
                                    const MultiplicationTables& tables) {
  if (IsGFNIBackend(backend)) {
    static_cast<void>(detail::GFNIIFFTRadix2Xor(
        x, y, output_x, output_y, byte_count, coefficient, backend, tables));
    return;
  }

  size_t processed = 0;
#if defined(__AVX2__)
  if (backend == Backend::avx2) {
    IFFTRadix2XorAvx2Impl(x, y, output_x, output_y, byte_count, coefficient,
                          tables);
    return;
  }
#endif
#if defined(__SSSE3__)
  if (backend == Backend::ssse3) {
    const Factor128 factor = LoadFactor128(coefficient, tables);
    for (; processed + 16 <= byte_count; processed += 16) {
      __m128i a =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(x + processed));
      __m128i b =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(y + processed));
      b = _mm_xor_si128(b, a);
      a = AddProduct128(a, b, factor, coefficient);
      const __m128i ox = _mm_loadu_si128(
          reinterpret_cast<const __m128i*>(output_x + processed));
      const __m128i oy = _mm_loadu_si128(
          reinterpret_cast<const __m128i*>(output_y + processed));
      _mm_storeu_si128(reinterpret_cast<__m128i*>(output_x + processed),
                       _mm_xor_si128(ox, a));
      _mm_storeu_si128(reinterpret_cast<__m128i*>(output_y + processed),
                       _mm_xor_si128(oy, b));
    }
    for (; processed < byte_count; ++processed) {
      const Element b_value = y[processed] ^ x[processed];
      const Element a_value =
          x[processed] ^ ScalarProduct(b_value, coefficient, tables);
      output_x[processed] ^= a_value;
      output_y[processed] ^= b_value;
    }
    return;
  }
#endif
  for (; processed < byte_count; ++processed) {
    const Element b_value = y[processed] ^ x[processed];
    const Element a_value =
        x[processed] ^ ScalarProduct(b_value, coefficient, tables);
    output_x[processed] ^= a_value;
    output_y[processed] ^= b_value;
  }
}

void detail::FFTRadix4Unchecked(Element* x0,
                                Element* x1,
                                Element* x2,
                                Element* x3,
                                size_t byte_count,
                                Element top,
                                Element low,
                                Element high,
                                Backend backend,
                                const MultiplicationTables& tables) {
  if (IsGFNIBackend(backend)) {
    static_cast<void>(detail::GFNIFFTRadix4(x0, x1, x2, x3, byte_count, top,
                                            low, high, backend, tables));
    return;
  }

  size_t processed = 0;
#if defined(__AVX2__)
  if (backend == Backend::avx2) {
    FFTRadix4Avx2Impl(x0, x1, x2, x3, byte_count, top, low, high, tables);
    return;
  }
#endif
#if defined(__SSSE3__)
  if (backend == Backend::ssse3) {
    const Factor128 top_factor = LoadFactor128(top, tables);
    const Factor128 low_factor = LoadFactor128(low, tables);
    const Factor128 high_factor = LoadFactor128(high, tables);
    for (; processed + 16 <= byte_count; processed += 16) {
      __m128i a =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(x0 + processed));
      __m128i b =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(x1 + processed));
      __m128i c =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(x2 + processed));
      __m128i d =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(x3 + processed));
      if (top != 0) {
        a = _mm_xor_si128(a, Product128(c, top_factor));
        b = _mm_xor_si128(b, Product128(d, top_factor));
      }
      c = _mm_xor_si128(c, a);
      d = _mm_xor_si128(d, b);
      if (low != 0) {
        a = _mm_xor_si128(a, Product128(b, low_factor));
      }
      b = _mm_xor_si128(b, a);
      if (high != 0) {
        c = _mm_xor_si128(c, Product128(d, high_factor));
      }
      d = _mm_xor_si128(d, c);
      _mm_storeu_si128(reinterpret_cast<__m128i*>(x0 + processed), a);
      _mm_storeu_si128(reinterpret_cast<__m128i*>(x1 + processed), b);
      _mm_storeu_si128(reinterpret_cast<__m128i*>(x2 + processed), c);
      _mm_storeu_si128(reinterpret_cast<__m128i*>(x3 + processed), d);
    }
    FFTRadix4ScalarTail(x0, x1, x2, x3, processed, byte_count, top, low, high,
                        tables);
    return;
  }
#endif
  FFTRadix4ScalarTail(x0, x1, x2, x3, processed, byte_count, top, low, high,
                      tables);
}

void detail::IFFTRadix4Unchecked(Element* x0,
                                 Element* x1,
                                 Element* x2,
                                 Element* x3,
                                 size_t byte_count,
                                 Element top,
                                 Element low,
                                 Element high,
                                 Backend backend,
                                 const MultiplicationTables& tables) {
  if (IsGFNIBackend(backend)) {
    static_cast<void>(detail::GFNIIFFTRadix4(x0, x1, x2, x3, byte_count, top,
                                             low, high, backend, tables));
    return;
  }

  size_t processed = 0;
#if defined(__AVX2__)
  if (backend == Backend::avx2) {
    IFFTRadix4Avx2Impl(x0, x1, x2, x3, byte_count, top, low, high, tables);
    return;
  }
#endif
#if defined(__SSSE3__)
  if (backend == Backend::ssse3) {
    const Factor128 top_factor = LoadFactor128(top, tables);
    const Factor128 low_factor = LoadFactor128(low, tables);
    const Factor128 high_factor = LoadFactor128(high, tables);
    for (; processed + 16 <= byte_count; processed += 16) {
      __m128i a =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(x0 + processed));
      __m128i b =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(x1 + processed));
      __m128i c =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(x2 + processed));
      __m128i d =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(x3 + processed));
      b = _mm_xor_si128(b, a);
      if (low != 0) {
        a = _mm_xor_si128(a, Product128(b, low_factor));
      }
      d = _mm_xor_si128(d, c);
      if (high != 0) {
        c = _mm_xor_si128(c, Product128(d, high_factor));
      }
      c = _mm_xor_si128(c, a);
      d = _mm_xor_si128(d, b);
      if (top != 0) {
        a = _mm_xor_si128(a, Product128(c, top_factor));
        b = _mm_xor_si128(b, Product128(d, top_factor));
      }
      _mm_storeu_si128(reinterpret_cast<__m128i*>(x0 + processed), a);
      _mm_storeu_si128(reinterpret_cast<__m128i*>(x1 + processed), b);
      _mm_storeu_si128(reinterpret_cast<__m128i*>(x2 + processed), c);
      _mm_storeu_si128(reinterpret_cast<__m128i*>(x3 + processed), d);
    }
    IFFTRadix4ScalarTail(x0, x1, x2, x3, processed, byte_count, top, low, high,
                         tables);
    return;
  }
#endif
  IFFTRadix4ScalarTail(x0, x1, x2, x3, processed, byte_count, top, low, high,
                       tables);
}

void detail::IFFTRadix4XorUnchecked(const Element* x0,
                                    const Element* x1,
                                    const Element* x2,
                                    const Element* x3,
                                    Element* output0,
                                    Element* output1,
                                    Element* output2,
                                    Element* output3,
                                    size_t byte_count,
                                    Element top,
                                    Element low,
                                    Element high,
                                    Backend backend,
                                    const MultiplicationTables& tables) {
  if (IsGFNIBackend(backend)) {
    static_cast<void>(detail::GFNIIFFTRadix4Xor(
        x0, x1, x2, x3, output0, output1, output2, output3, byte_count, top,
        low, high, backend, tables));
    return;
  }

  size_t processed = 0;
#if defined(__AVX2__)
  if (backend == Backend::avx2) {
    IFFTRadix4XorAvx2Impl(x0, x1, x2, x3, output0, output1, output2, output3,
                          byte_count, top, low, high, tables);
    return;
  }
#endif
#if defined(__SSSE3__)
  if (backend == Backend::ssse3) {
    const Factor128 top_factor = LoadFactor128(top, tables);
    const Factor128 low_factor = LoadFactor128(low, tables);
    const Factor128 high_factor = LoadFactor128(high, tables);
    for (; processed + 16 <= byte_count; processed += 16) {
      __m128i a =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(x0 + processed));
      __m128i b =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(x1 + processed));
      __m128i c =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(x2 + processed));
      __m128i d =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(x3 + processed));
      b = _mm_xor_si128(b, a);
      if (low != 0) {
        a = _mm_xor_si128(a, Product128(b, low_factor));
      }
      d = _mm_xor_si128(d, c);
      if (high != 0) {
        c = _mm_xor_si128(c, Product128(d, high_factor));
      }
      c = _mm_xor_si128(c, a);
      d = _mm_xor_si128(d, b);
      if (top != 0) {
        a = _mm_xor_si128(a, Product128(c, top_factor));
        b = _mm_xor_si128(b, Product128(d, top_factor));
      }
      const __m128i o0 = _mm_loadu_si128(
          reinterpret_cast<const __m128i*>(output0 + processed));
      const __m128i o1 = _mm_loadu_si128(
          reinterpret_cast<const __m128i*>(output1 + processed));
      const __m128i o2 = _mm_loadu_si128(
          reinterpret_cast<const __m128i*>(output2 + processed));
      const __m128i o3 = _mm_loadu_si128(
          reinterpret_cast<const __m128i*>(output3 + processed));
      _mm_storeu_si128(reinterpret_cast<__m128i*>(output0 + processed),
                       _mm_xor_si128(o0, a));
      _mm_storeu_si128(reinterpret_cast<__m128i*>(output1 + processed),
                       _mm_xor_si128(o1, b));
      _mm_storeu_si128(reinterpret_cast<__m128i*>(output2 + processed),
                       _mm_xor_si128(o2, c));
      _mm_storeu_si128(reinterpret_cast<__m128i*>(output3 + processed),
                       _mm_xor_si128(o3, d));
    }
    IFFTRadix4XorScalarTail(x0, x1, x2, x3, output0, output1, output2, output3,
                            processed, byte_count, top, low, high, tables);
    return;
  }
#endif
  IFFTRadix4XorScalarTail(x0, x1, x2, x3, output0, output1, output2, output3,
                          processed, byte_count, top, low, high, tables);
}

namespace {

template <Backend backend>
void AddScaledResolved(Element* destination,
                       const Element* source,
                       size_t byte_count,
                       Element coefficient,
                       const MultiplicationTables& tables) {
  static_cast<void>(AddScaledImpl(destination, source, byte_count, coefficient,
                                  backend, tables));
}

template <Backend backend>
void ScaleResolved(Element* destination,
                   const Element* source,
                   size_t byte_count,
                   Element coefficient,
                   const MultiplicationTables& tables) {
  static_cast<void>(MultiplyImpl(destination, source, byte_count, coefficient,
                                 backend, tables));
}

template <Backend backend>
void XorResolved(Element* destination,
                 const Element* source,
                 size_t byte_count) {
  XorImpl(destination, source, byte_count, backend);
}

template <Backend backend>
void Xor4Resolved(Element* destination0,
                  const Element* source0,
                  Element* destination1,
                  const Element* source1,
                  Element* destination2,
                  const Element* source2,
                  Element* destination3,
                  const Element* source3,
                  size_t byte_count) {
  Xor4Impl(destination0, source0, destination1, source1, destination2, source2,
           destination3, source3, byte_count, backend);
}

template <Backend backend>
void FFTRadix2Resolved(Element* x,
                       Element* y,
                       size_t byte_count,
                       Element coefficient,
                       const MultiplicationTables& tables) {
  detail::FFTRadix2Unchecked(x, y, byte_count, coefficient, backend, tables);
}

template <Backend backend>
void IFFTRadix2Resolved(Element* x,
                        Element* y,
                        size_t byte_count,
                        Element coefficient,
                        const MultiplicationTables& tables) {
  detail::IFFTRadix2Unchecked(x, y, byte_count, coefficient, backend, tables);
}

template <Backend backend>
void IFFTRadix2XorResolved(const Element* x,
                           const Element* y,
                           Element* output_x,
                           Element* output_y,
                           size_t byte_count,
                           Element coefficient,
                           const MultiplicationTables& tables) {
  detail::IFFTRadix2XorUnchecked(x, y, output_x, output_y, byte_count,
                                 coefficient, backend, tables);
}

template <Backend backend>
void FFTRadix4Resolved(Element* x0,
                       Element* x1,
                       Element* x2,
                       Element* x3,
                       size_t byte_count,
                       Element top,
                       Element low,
                       Element high,
                       const MultiplicationTables& tables) {
  detail::FFTRadix4Unchecked(x0, x1, x2, x3, byte_count, top, low, high,
                             backend, tables);
}

template <Backend backend>
void IFFTRadix4Resolved(Element* x0,
                        Element* x1,
                        Element* x2,
                        Element* x3,
                        size_t byte_count,
                        Element top,
                        Element low,
                        Element high,
                        const MultiplicationTables& tables) {
  detail::IFFTRadix4Unchecked(x0, x1, x2, x3, byte_count, top, low, high,
                              backend, tables);
}

template <Backend backend>
Status IFFTRadix2CopyGFNIResolved(const Element* x,
                                  const Element* y,
                                  Element* output_x,
                                  Element* output_y,
                                  size_t byte_count,
                                  Element coefficient,
                                  const MultiplicationTables& tables) {
  return detail::GFNIIFFTRadix2Copy(x, y, output_x, output_y, byte_count,
                                    coefficient, backend, tables);
}

template <Backend backend>
Status IFFTRadix4CopyGFNIResolved(const Element* x0,
                                  const Element* x1,
                                  const Element* x2,
                                  const Element* x3,
                                  Element* output0,
                                  Element* output1,
                                  Element* output2,
                                  Element* output3,
                                  size_t byte_count,
                                  Element top,
                                  Element low,
                                  Element high,
                                  const MultiplicationTables& tables) {
  return detail::GFNIIFFTRadix4Copy(x0, x1, x2, x3, output0, output1, output2,
                                    output3, byte_count, top, low, high,
                                    backend, tables);
}

template <Backend backend>
void IFFTRadix4XorResolved(const Element* x0,
                           const Element* x1,
                           const Element* x2,
                           const Element* x3,
                           Element* output0,
                           Element* output1,
                           Element* output2,
                           Element* output3,
                           size_t byte_count,
                           Element top,
                           Element low,
                           Element high,
                           const MultiplicationTables& tables) {
  detail::IFFTRadix4XorUnchecked(x0, x1, x2, x3, output0, output1, output2,
                                 output3, byte_count, top, low, high, backend,
                                 tables);
}

#if defined(__AVX2__)
template <>
void ScaleResolved<Backend::avx2>(Element* destination,
                                  const Element* source,
                                  size_t byte_count,
                                  Element coefficient,
                                  const MultiplicationTables& tables) {
  ScaleAvx2Impl(destination, source, byte_count, coefficient, tables);
}

template <>
void XorResolved<Backend::avx2>(Element* destination,
                                const Element* source,
                                size_t byte_count) {
  XorAvx2Impl(destination, source, byte_count);
}

template <>
void Xor4Resolved<Backend::avx2>(Element* destination0,
                                 const Element* source0,
                                 Element* destination1,
                                 const Element* source1,
                                 Element* destination2,
                                 const Element* source2,
                                 Element* destination3,
                                 const Element* source3,
                                 size_t byte_count) {
  Xor4Avx2Impl(destination0, source0, destination1, source1, destination2,
               source2, destination3, source3, byte_count);
}

template <>
void FFTRadix2Resolved<Backend::avx2>(Element* x,
                                      Element* y,
                                      size_t byte_count,
                                      Element coefficient,
                                      const MultiplicationTables& tables) {
  FFTRadix2Avx2Impl(x, y, byte_count, coefficient, tables);
}

template <>
void IFFTRadix2Resolved<Backend::avx2>(Element* x,
                                       Element* y,
                                       size_t byte_count,
                                       Element coefficient,
                                       const MultiplicationTables& tables) {
  IFFTRadix2Avx2Impl(x, y, byte_count, coefficient, tables);
}

template <>
void IFFTRadix2XorResolved<Backend::avx2>(const Element* x,
                                          const Element* y,
                                          Element* output_x,
                                          Element* output_y,
                                          size_t byte_count,
                                          Element coefficient,
                                          const MultiplicationTables& tables) {
  IFFTRadix2XorAvx2Impl(x, y, output_x, output_y, byte_count, coefficient,
                        tables);
}

template <>
void FFTRadix4Resolved<Backend::avx2>(Element* x0,
                                      Element* x1,
                                      Element* x2,
                                      Element* x3,
                                      size_t byte_count,
                                      Element top,
                                      Element low,
                                      Element high,
                                      const MultiplicationTables& tables) {
  FFTRadix4Avx2Impl(x0, x1, x2, x3, byte_count, top, low, high, tables);
}

template <>
void IFFTRadix4Resolved<Backend::avx2>(Element* x0,
                                       Element* x1,
                                       Element* x2,
                                       Element* x3,
                                       size_t byte_count,
                                       Element top,
                                       Element low,
                                       Element high,
                                       const MultiplicationTables& tables) {
  IFFTRadix4Avx2Impl(x0, x1, x2, x3, byte_count, top, low, high, tables);
}

template <>
void IFFTRadix4XorResolved<Backend::avx2>(const Element* x0,
                                          const Element* x1,
                                          const Element* x2,
                                          const Element* x3,
                                          Element* output0,
                                          Element* output1,
                                          Element* output2,
                                          Element* output3,
                                          size_t byte_count,
                                          Element top,
                                          Element low,
                                          Element high,
                                          const MultiplicationTables& tables) {
  IFFTRadix4XorAvx2Impl(x0, x1, x2, x3, output0, output1, output2, output3,
                        byte_count, top, low, high, tables);
}

void FFTRadix4Avx2ExactResolved(Element* x0,
                                Element* x1,
                                Element* x2,
                                Element* x3,
                                size_t byte_count,
                                Element top,
                                Element low,
                                Element high,
                                const MultiplicationTables& tables) {
  FFTRadix4Avx2Impl<true>(x0, x1, x2, x3, byte_count, top, low, high, tables);
}

void IFFTRadix4Avx2ExactResolved(Element* x0,
                                 Element* x1,
                                 Element* x2,
                                 Element* x3,
                                 size_t byte_count,
                                 Element top,
                                 Element low,
                                 Element high,
                                 const MultiplicationTables& tables) {
  IFFTRadix4Avx2Impl<true>(x0, x1, x2, x3, byte_count, top, low, high, tables);
}

void IFFTRadix4XorAvx2ExactResolved(const Element* x0,
                                    const Element* x1,
                                    const Element* x2,
                                    const Element* x3,
                                    Element* output0,
                                    Element* output1,
                                    Element* output2,
                                    Element* output3,
                                    size_t byte_count,
                                    Element top,
                                    Element low,
                                    Element high,
                                    const MultiplicationTables& tables) {
  IFFTRadix4XorAvx2Impl<true>(x0, x1, x2, x3, output0, output1, output2,
                              output3, byte_count, top, low, high, tables);
}

void FFTRadix2Avx2ExactResolved(Element* x,
                                Element* y,
                                size_t byte_count,
                                Element coefficient,
                                const MultiplicationTables& tables) {
  FFTRadix2Avx2Impl<true>(x, y, byte_count, coefficient, tables);
}

void IFFTRadix2Avx2ExactResolved(Element* x,
                                 Element* y,
                                 size_t byte_count,
                                 Element coefficient,
                                 const MultiplicationTables& tables) {
  IFFTRadix2Avx2Impl<true>(x, y, byte_count, coefficient, tables);
}

void IFFTRadix2XorAvx2ExactResolved(const Element* x,
                                    const Element* y,
                                    Element* output_x,
                                    Element* output_y,
                                    size_t byte_count,
                                    Element coefficient,
                                    const MultiplicationTables& tables) {
  IFFTRadix2XorAvx2Impl<true>(x, y, output_x, output_y, byte_count, coefficient,
                              tables);
}

Status IFFTRadix2CopyAvx2Resolved(const Element* x,
                                  const Element* y,
                                  Element* output_x,
                                  Element* output_y,
                                  size_t byte_count,
                                  Element coefficient,
                                  const MultiplicationTables& tables) {
  IFFTRadix2CopyAvx2Impl(x, y, output_x, output_y, byte_count, coefficient,
                         tables);
  return Status::ok;
}

Status IFFTRadix2CopyAvx2ExactResolved(const Element* x,
                                       const Element* y,
                                       Element* output_x,
                                       Element* output_y,
                                       size_t byte_count,
                                       Element coefficient,
                                       const MultiplicationTables& tables) {
  IFFTRadix2CopyAvx2Impl<true>(x, y, output_x, output_y, byte_count,
                               coefficient, tables);
  return Status::ok;
}

Status IFFTRadix4CopyAvx2Resolved(const Element* x0,
                                  const Element* x1,
                                  const Element* x2,
                                  const Element* x3,
                                  Element* output0,
                                  Element* output1,
                                  Element* output2,
                                  Element* output3,
                                  size_t byte_count,
                                  Element top,
                                  Element low,
                                  Element high,
                                  const MultiplicationTables& tables) {
  IFFTRadix4CopyAvx2Impl(x0, x1, x2, x3, output0, output1, output2, output3,
                         byte_count, top, low, high, tables);
  return Status::ok;
}

Status IFFTRadix4CopyAvx2ExactResolved(const Element* x0,
                                       const Element* x1,
                                       const Element* x2,
                                       const Element* x3,
                                       Element* output0,
                                       Element* output1,
                                       Element* output2,
                                       Element* output3,
                                       size_t byte_count,
                                       Element top,
                                       Element low,
                                       Element high,
                                       const MultiplicationTables& tables) {
  IFFTRadix4CopyAvx2Impl<true>(x0, x1, x2, x3, output0, output1, output2,
                               output3, byte_count, top, low, high, tables);
  return Status::ok;
}

const detail::ResolvedKernels& Avx2ExactKernelSet() {
  static const detail::ResolvedKernels kernels{
      AddScaledResolved<Backend::avx2>, ScaleResolved<Backend::avx2>,
      XorResolved<Backend::avx2>,       Xor4Resolved<Backend::avx2>,
      IFFTRadix2CopyAvx2ExactResolved,  IFFTRadix4CopyAvx2ExactResolved,
      FFTRadix2Avx2ExactResolved,       IFFTRadix2Avx2ExactResolved,
      IFFTRadix2XorAvx2ExactResolved,   FFTRadix4Avx2ExactResolved,
      IFFTRadix4Avx2ExactResolved,      IFFTRadix4XorAvx2ExactResolved};
  return kernels;
}
#endif

template <Backend backend>
const detail::ResolvedKernels& KernelSet() {
#if defined(__AVX2__)
  if constexpr (backend == Backend::avx2) {
    static const detail::ResolvedKernels kernels{
        AddScaledResolved<backend>,     ScaleResolved<backend>,
        XorResolved<backend>,           Xor4Resolved<backend>,
        IFFTRadix2CopyAvx2Resolved,     IFFTRadix4CopyAvx2Resolved,
        FFTRadix2Resolved<backend>,     IFFTRadix2Resolved<backend>,
        IFFTRadix2XorResolved<backend>, FFTRadix4Resolved<backend>,
        IFFTRadix4Resolved<backend>,    IFFTRadix4XorResolved<backend>};
    return kernels;
  }
#endif
  if constexpr (backend == Backend::gfni128_affine ||
                backend == Backend::gfni256_affine ||
                backend == Backend::gfni512_affine) {
    static const detail::ResolvedKernels kernels{
        AddScaledResolved<backend>,
        ScaleResolved<backend>,
        XorResolved<backend>,
        Xor4Resolved<backend>,
        IFFTRadix2CopyGFNIResolved<backend>,
        IFFTRadix4CopyGFNIResolved<backend>,
        FFTRadix2Resolved<backend>,
        IFFTRadix2Resolved<backend>,
        IFFTRadix2XorResolved<backend>,
        FFTRadix4Resolved<backend>,
        IFFTRadix4Resolved<backend>,
        IFFTRadix4XorResolved<backend>};
    return kernels;
  }
  static const detail::ResolvedKernels kernels{AddScaledResolved<backend>,
                                               ScaleResolved<backend>,
                                               XorResolved<backend>,
                                               Xor4Resolved<backend>,
                                               nullptr,
                                               nullptr,
                                               FFTRadix2Resolved<backend>,
                                               IFFTRadix2Resolved<backend>,
                                               IFFTRadix2XorResolved<backend>,
                                               FFTRadix4Resolved<backend>,
                                               IFFTRadix4Resolved<backend>,
                                               IFFTRadix4XorResolved<backend>};
  return kernels;
}

}  // namespace

const detail::ResolvedKernels* detail::ResolveKernels(Backend backend) {
  if (backend == Backend::tuned || !BackendAvailable(backend)) {
    return nullptr;
  }
  switch (backend) {
    case Backend::scalar:
      return &KernelSet<Backend::scalar>();
    case Backend::ssse3:
      return &KernelSet<Backend::ssse3>();
    case Backend::avx2:
      return &KernelSet<Backend::avx2>();
    case Backend::gfni128_affine:
      return &KernelSet<Backend::gfni128_affine>();
    case Backend::gfni256_affine:
      return &KernelSet<Backend::gfni256_affine>();
    case Backend::gfni512_affine:
      return &KernelSet<Backend::gfni512_affine>();
    case Backend::tuned:
      break;
  }
  return nullptr;
}

const detail::ResolvedKernels* detail::ResolveKernels(Backend backend,
                                                      size_t byte_count) {
#if defined(__AVX2__)
  if (backend == Backend::avx2 && (byte_count & 31U) == 0) {
    return &Avx2ExactKernelSet();
  }
#else
  static_cast<void>(byte_count);
#endif
  return ResolveKernels(backend);
}

Status FFTRadix2(Element* x,
                 Element* y,
                 size_t byte_count,
                 Element coefficient,
                 Backend backend,
                 const MultiplicationTables& tables) {
  if (backend == Backend::tuned) {
    backend = SelectBackend(byte_count);
  }
  if (!BackendAvailable(backend)) {
    return Status::unsupported_backend;
  }
  detail::FFTRadix2Unchecked(x, y, byte_count, coefficient, backend, tables);
  return Status::ok;
}

Status IFFTRadix2(Element* x,
                  Element* y,
                  size_t byte_count,
                  Element coefficient,
                  Backend backend,
                  const MultiplicationTables& tables) {
  if (backend == Backend::tuned) {
    backend = SelectBackend(byte_count);
  }
  if (!BackendAvailable(backend)) {
    return Status::unsupported_backend;
  }
  detail::IFFTRadix2Unchecked(x, y, byte_count, coefficient, backend, tables);
  return Status::ok;
}

Status IFFTRadix2Xor(const Element* x,
                     const Element* y,
                     Element* output_x,
                     Element* output_y,
                     size_t byte_count,
                     Element coefficient,
                     Backend backend,
                     const MultiplicationTables& tables) {
  if (backend == Backend::tuned) {
    backend = SelectBackend(byte_count);
  }
  if (!BackendAvailable(backend)) {
    return Status::unsupported_backend;
  }
  detail::IFFTRadix2XorUnchecked(x, y, output_x, output_y, byte_count,
                                 coefficient, backend, tables);
  return Status::ok;
}

Status FFTRadix4(Element* x0,
                 Element* x1,
                 Element* x2,
                 Element* x3,
                 size_t byte_count,
                 Element top,
                 Element low,
                 Element high,
                 Backend backend,
                 const MultiplicationTables& tables) {
  if (backend == Backend::tuned) {
    backend = SelectBackend(byte_count);
  }
  if (!BackendAvailable(backend)) {
    return Status::unsupported_backend;
  }
  detail::FFTRadix4Unchecked(x0, x1, x2, x3, byte_count, top, low, high,
                             backend, tables);
  return Status::ok;
}

Status IFFTRadix4(Element* x0,
                  Element* x1,
                  Element* x2,
                  Element* x3,
                  size_t byte_count,
                  Element top,
                  Element low,
                  Element high,
                  Backend backend,
                  const MultiplicationTables& tables) {
  if (backend == Backend::tuned) {
    backend = SelectBackend(byte_count);
  }
  if (!BackendAvailable(backend)) {
    return Status::unsupported_backend;
  }
  detail::IFFTRadix4Unchecked(x0, x1, x2, x3, byte_count, top, low, high,
                              backend, tables);
  return Status::ok;
}

Status IFFTRadix4Xor(const Element* x0,
                     const Element* x1,
                     const Element* x2,
                     const Element* x3,
                     Element* output0,
                     Element* output1,
                     Element* output2,
                     Element* output3,
                     size_t byte_count,
                     Element top,
                     Element low,
                     Element high,
                     Backend backend,
                     const MultiplicationTables& tables) {
  if (backend == Backend::tuned) {
    backend = SelectBackend(byte_count);
  }
  if (!BackendAvailable(backend)) {
    return Status::unsupported_backend;
  }
  detail::IFFTRadix4XorUnchecked(x0, x1, x2, x3, output0, output1, output2,
                                 output3, byte_count, top, low, high, backend,
                                 tables);
  return Status::ok;
}

}  // namespace gf2p8::lch
