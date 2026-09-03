#include "lin_chung_han/kernels.h"

/**
 * GFNI realization of the Leopard-style LCH butterfly schedule using
 * affine transform.
 */
#if defined(__i386__) || defined(__x86_64__) || defined(_M_IX86) || \
    defined(_M_X64)
#include <immintrin.h>
#endif

namespace gf2p8::lch::detail {
namespace {

void ScalarTail(Element* x,
                const Element* y,
                size_t begin,
                size_t byte_count,
                Element coefficient,
                const MultiplicationTables& tables) {
  const auto& row = tables.shuffle[coefficient];
  for (size_t i = begin; i < byte_count; ++i) {
    x[i] ^= row[y[i] & 0x0f] ^ row[32 + (y[i] >> 4)];
  }
}

void ScalarMultiplyTail(Element* destination,
                        const Element* source,
                        size_t begin,
                        size_t byte_count,
                        Element coefficient,
                        const MultiplicationTables& tables) {
  const auto& row = tables.shuffle[coefficient];
  for (size_t i = begin; i < byte_count; ++i) {
    destination[i] = row[source[i] & 0x0f] ^ row[32 + (source[i] >> 4)];
  }
}

Element ScalarProduct(Element value,
                      Element coefficient,
                      const MultiplicationTables& tables) {
  const auto& row = tables.shuffle[coefficient];
  return row[value & 0x0f] ^ row[32 + (value >> 4)];
}

void FFTScalarTail(Element* x0,
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

void IFFTScalarTail(Element* x0,
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

template <typename Ops>
typename Ops::Vector AddProduct(typename Ops::Vector accumulator,
                                typename Ops::Vector input,
                                typename Ops::Vector factor,
                                Element coefficient,
                                bool affine) {
  if (coefficient == 0) {
    return accumulator;
  }
  return Ops::Xor(accumulator, Ops::Product(input, factor, affine));
}

template <typename Ops>
size_t FFTVectors(Element* x0,
                  Element* x1,
                  Element* x2,
                  Element* x3,
                  size_t byte_count,
                  Element top,
                  Element low,
                  Element high,
                  bool affine,
                  const MultiplicationTables& tables) {
  const auto top_factor = Ops::Factor(top, affine, tables);
  const auto low_factor = Ops::Factor(low, affine, tables);
  const auto high_factor = Ops::Factor(high, affine, tables);
  size_t processed = 0;
  for (; processed + Ops::kWidth <= byte_count; processed += Ops::kWidth) {
    auto a = Ops::Load(x0 + processed);
    auto b = Ops::Load(x1 + processed);
    auto c = Ops::Load(x2 + processed);
    auto d = Ops::Load(x3 + processed);
    if (top != 0) {
      a = Ops::Xor(a, Ops::Product(c, top_factor, affine));
      b = Ops::Xor(b, Ops::Product(d, top_factor, affine));
    }
    c = Ops::Xor(c, a);
    d = Ops::Xor(d, b);
    if (low != 0) {
      a = Ops::Xor(a, Ops::Product(b, low_factor, affine));
    }
    b = Ops::Xor(b, a);
    if (high != 0) {
      c = Ops::Xor(c, Ops::Product(d, high_factor, affine));
    }
    d = Ops::Xor(d, c);
    Ops::Store(x0 + processed, a);
    Ops::Store(x1 + processed, b);
    Ops::Store(x2 + processed, c);
    Ops::Store(x3 + processed, d);
  }
  return processed;
}

template <typename Ops>
size_t IFFTVectors(Element* x0,
                   Element* x1,
                   Element* x2,
                   Element* x3,
                   size_t byte_count,
                   Element top,
                   Element low,
                   Element high,
                   bool affine,
                   const MultiplicationTables& tables) {
  const auto top_factor = Ops::Factor(top, affine, tables);
  const auto low_factor = Ops::Factor(low, affine, tables);
  const auto high_factor = Ops::Factor(high, affine, tables);
  size_t processed = 0;
  for (; processed + Ops::kWidth <= byte_count; processed += Ops::kWidth) {
    auto a = Ops::Load(x0 + processed);
    auto b = Ops::Load(x1 + processed);
    auto c = Ops::Load(x2 + processed);
    auto d = Ops::Load(x3 + processed);
    b = Ops::Xor(b, a);
    if (low != 0) {
      a = Ops::Xor(a, Ops::Product(b, low_factor, affine));
    }
    d = Ops::Xor(d, c);
    if (high != 0) {
      c = Ops::Xor(c, Ops::Product(d, high_factor, affine));
    }
    c = Ops::Xor(c, a);
    d = Ops::Xor(d, b);
    if (top != 0) {
      a = Ops::Xor(a, Ops::Product(c, top_factor, affine));
      b = Ops::Xor(b, Ops::Product(d, top_factor, affine));
    }
    Ops::Store(x0 + processed, a);
    Ops::Store(x1 + processed, b);
    Ops::Store(x2 + processed, c);
    Ops::Store(x3 + processed, d);
  }
  return processed;
}

template <typename Ops>
size_t IFFTRadix2CopyVectors(const Element* x,
                             const Element* y,
                             Element* output_x,
                             Element* output_y,
                             size_t byte_count,
                             Element coefficient,
                             bool affine,
                             const MultiplicationTables& tables) {
  const auto factor = Ops::Factor(coefficient, affine, tables);
  size_t processed = 0;
  for (; processed + Ops::kWidth <= byte_count; processed += Ops::kWidth) {
    auto a = Ops::Load(x + processed);
    auto b = Ops::Xor(Ops::Load(y + processed), a);
    a = AddProduct<Ops>(a, b, factor, coefficient, affine);
    Ops::Store(output_x + processed, a);
    Ops::Store(output_y + processed, b);
  }
  return processed;
}

template <typename Ops>
size_t IFFTRadix2XorVectors(const Element* x,
                            const Element* y,
                            Element* output_x,
                            Element* output_y,
                            size_t byte_count,
                            Element coefficient,
                            bool affine,
                            const MultiplicationTables& tables) {
  const auto factor = Ops::Factor(coefficient, affine, tables);
  size_t processed = 0;
  for (; processed + Ops::kWidth <= byte_count; processed += Ops::kWidth) {
    auto a = Ops::Load(x + processed);
    auto b = Ops::Load(y + processed);
    b = Ops::Xor(b, a);
    a = AddProduct<Ops>(a, b, factor, coefficient, affine);
    Ops::Store(output_x + processed,
               Ops::Xor(Ops::Load(output_x + processed), a));
    Ops::Store(output_y + processed,
               Ops::Xor(Ops::Load(output_y + processed), b));
  }
  return processed;
}

template <typename Ops>
size_t IFFTRadix4CopyVectors(const Element* x0,
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
                             bool affine,
                             const MultiplicationTables& tables) {
  const auto top_factor = Ops::Factor(top, affine, tables);
  const auto low_factor = Ops::Factor(low, affine, tables);
  const auto high_factor = Ops::Factor(high, affine, tables);
  size_t processed = 0;
  for (; processed + Ops::kWidth <= byte_count; processed += Ops::kWidth) {
    auto a = Ops::Load(x0 + processed);
    auto b = Ops::Xor(Ops::Load(x1 + processed), a);
    if (low != 0) {
      a = Ops::Xor(a, Ops::Product(b, low_factor, affine));
    }
    auto c = Ops::Load(x2 + processed);
    auto d = Ops::Xor(Ops::Load(x3 + processed), c);
    if (high != 0) {
      c = Ops::Xor(c, Ops::Product(d, high_factor, affine));
    }
    c = Ops::Xor(c, a);
    d = Ops::Xor(d, b);
    if (top != 0) {
      a = Ops::Xor(a, Ops::Product(c, top_factor, affine));
      b = Ops::Xor(b, Ops::Product(d, top_factor, affine));
    }
    Ops::Store(output0 + processed, a);
    Ops::Store(output1 + processed, b);
    Ops::Store(output2 + processed, c);
    Ops::Store(output3 + processed, d);
  }
  return processed;
}

template <typename Ops>
size_t IFFTRadix4XorVectors(const Element* x0,
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
                            bool affine,
                            const MultiplicationTables& tables) {
  const auto top_factor = Ops::Factor(top, affine, tables);
  const auto low_factor = Ops::Factor(low, affine, tables);
  const auto high_factor = Ops::Factor(high, affine, tables);
  size_t processed = 0;
  for (; processed + Ops::kWidth <= byte_count; processed += Ops::kWidth) {
    auto a = Ops::Load(x0 + processed);
    auto b = Ops::Load(x1 + processed);
    auto c = Ops::Load(x2 + processed);
    auto d = Ops::Load(x3 + processed);
    b = Ops::Xor(b, a);
    if (low != 0) {
      a = Ops::Xor(a, Ops::Product(b, low_factor, affine));
    }
    d = Ops::Xor(d, c);
    if (high != 0) {
      c = Ops::Xor(c, Ops::Product(d, high_factor, affine));
    }
    c = Ops::Xor(c, a);
    d = Ops::Xor(d, b);
    if (top != 0) {
      a = Ops::Xor(a, Ops::Product(c, top_factor, affine));
      b = Ops::Xor(b, Ops::Product(d, top_factor, affine));
    }
    Ops::Store(output0 + processed,
               Ops::Xor(Ops::Load(output0 + processed), a));
    Ops::Store(output1 + processed,
               Ops::Xor(Ops::Load(output1 + processed), b));
    Ops::Store(output2 + processed,
               Ops::Xor(Ops::Load(output2 + processed), c));
    Ops::Store(output3 + processed,
               Ops::Xor(Ops::Load(output3 + processed), d));
  }
  return processed;
}

#if defined(__GFNI__) && defined(__SSE2__)
struct GFNI128Ops {
  static constexpr size_t kWidth = 16;
  using Vector = __m128i;
  static Vector Factor(Element coefficient,
                       bool,
                       const MultiplicationTables& tables) {
    return _mm_set1_epi64x(static_cast<long long>(tables.affine[coefficient]));
  }
  static Vector Load(const Element* pointer) {
    return _mm_loadu_si128(reinterpret_cast<const Vector*>(pointer));
  }
  static void Store(Element* pointer, Vector value) {
    _mm_storeu_si128(reinterpret_cast<Vector*>(pointer), value);
  }
  static Vector Xor(Vector x, Vector y) { return _mm_xor_si128(x, y); }
  static Vector Product(Vector input, Vector factor, bool) {
    return _mm_gf2p8affine_epi64_epi8(input, factor, 0);
  }
};
#endif

#if defined(__GFNI__) && defined(__AVX2__)
struct GFNI256Ops {
  static constexpr size_t kWidth = 32;
  using Vector = __m256i;
  static Vector Factor(Element coefficient,
                       bool,
                       const MultiplicationTables& tables) {
    return _mm256_set1_epi64x(
        static_cast<long long>(tables.affine[coefficient]));
  }
  static Vector Load(const Element* pointer) {
    return _mm256_loadu_si256(reinterpret_cast<const Vector*>(pointer));
  }
  static void Store(Element* pointer, Vector value) {
    _mm256_storeu_si256(reinterpret_cast<Vector*>(pointer), value);
  }
  static Vector Xor(Vector x, Vector y) { return _mm256_xor_si256(x, y); }
  static Vector Product(Vector input, Vector factor, bool) {
    return _mm256_gf2p8affine_epi64_epi8(input, factor, 0);
  }
};
#endif

#if defined(__GFNI__) && defined(__AVX512F__) && defined(__AVX512BW__)
struct GFNI512Ops {
  static constexpr size_t kWidth = 64;
  using Vector = __m512i;
  static Vector Factor(Element coefficient,
                       bool,
                       const MultiplicationTables& tables) {
    return _mm512_set1_epi64(
        static_cast<long long>(tables.affine[coefficient]));
  }
  static Vector Load(const Element* pointer) {
    return _mm512_loadu_si512(pointer);
  }
  static void Store(Element* pointer, Vector value) {
    _mm512_storeu_si512(pointer, value);
  }
  static Vector Xor(Vector x, Vector y) { return _mm512_xor_si512(x, y); }
  static Vector Product(Vector input, Vector factor, bool) {
    return _mm512_gf2p8affine_epi64_epi8(input, factor, 0);
  }
};
#endif

}  // namespace

Status GFNIMultiplyAdd(Element* x,
                       const Element* y,
                       size_t byte_count,
                       Element coefficient,
                       Backend backend,
                       const MultiplicationTables& tables) {
  size_t processed = 0;
  switch (backend) {
    case Backend::gfni512_affine:
#if defined(__GFNI__) && defined(__AVX512F__) && defined(__AVX512BW__)
    {
      const __m512i matrix =
          _mm512_set1_epi64(static_cast<long long>(tables.affine[coefficient]));
      for (; processed + 64 <= byte_count; processed += 64) {
        const __m512i input = _mm512_loadu_si512(y + processed);
        const __m512i output = _mm512_loadu_si512(x + processed);
        _mm512_storeu_si512(
            x + processed,
            _mm512_xor_si512(output,
                             _mm512_gf2p8affine_epi64_epi8(input, matrix, 0)));
      }
      break;
    }
#else
      return Status::unsupported_backend;
#endif
    case Backend::gfni256_affine:
#if defined(__GFNI__) && defined(__AVX2__)
    {
      const __m256i matrix = _mm256_set1_epi64x(
          static_cast<long long>(tables.affine[coefficient]));
      for (; processed + 32 <= byte_count; processed += 32) {
        const __m256i input =
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(y + processed));
        const __m256i output =
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x + processed));
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(x + processed),
            _mm256_xor_si256(output,
                             _mm256_gf2p8affine_epi64_epi8(input, matrix, 0)));
      }
      break;
    }
#else
      return Status::unsupported_backend;
#endif
    case Backend::gfni128_affine:
#if defined(__GFNI__) && defined(__SSE2__)
    {
      const __m128i matrix =
          _mm_set1_epi64x(static_cast<long long>(tables.affine[coefficient]));
      for (; processed + 16 <= byte_count; processed += 16) {
        const __m128i input =
            _mm_loadu_si128(reinterpret_cast<const __m128i*>(y + processed));
        const __m128i output =
            _mm_loadu_si128(reinterpret_cast<const __m128i*>(x + processed));
        _mm_storeu_si128(reinterpret_cast<__m128i*>(x + processed),
                         _mm_xor_si128(output, _mm_gf2p8affine_epi64_epi8(
                                                   input, matrix, 0)));
      }
      break;
    }
#else
      return Status::unsupported_backend;
#endif
    default:
      return Status::invalid_argument;
  }

  ScalarTail(x, y, processed, byte_count, coefficient, tables);
  return Status::ok;
}

Status GFNIMultiply(Element* destination,
                    const Element* source,
                    size_t byte_count,
                    Element coefficient,
                    Backend backend,
                    const MultiplicationTables& tables) {
  size_t processed = 0;
  [[maybe_unused]] constexpr bool affine = true;
  switch (backend) {
    case Backend::gfni128_affine:
#if defined(__GFNI__) && defined(__SSE2__)
    {
      const auto factor = GFNI128Ops::Factor(coefficient, affine, tables);
      for (; processed + GFNI128Ops::kWidth <= byte_count;
           processed += GFNI128Ops::kWidth) {
        GFNI128Ops::Store(
            destination + processed,
            GFNI128Ops::Product(GFNI128Ops::Load(source + processed), factor,
                                affine));
      }
      break;
    }
#else
      return Status::unsupported_backend;
#endif
    case Backend::gfni256_affine:
#if defined(__GFNI__) && defined(__AVX2__)
    {
      const auto factor = GFNI256Ops::Factor(coefficient, affine, tables);
      for (; processed + GFNI256Ops::kWidth <= byte_count;
           processed += GFNI256Ops::kWidth) {
        GFNI256Ops::Store(
            destination + processed,
            GFNI256Ops::Product(GFNI256Ops::Load(source + processed), factor,
                                affine));
      }
      break;
    }
#else
      return Status::unsupported_backend;
#endif
    case Backend::gfni512_affine:
#if defined(__GFNI__) && defined(__AVX512F__) && defined(__AVX512BW__)
    {
      const auto factor = GFNI512Ops::Factor(coefficient, affine, tables);
      for (; processed + GFNI512Ops::kWidth <= byte_count;
           processed += GFNI512Ops::kWidth) {
        GFNI512Ops::Store(
            destination + processed,
            GFNI512Ops::Product(GFNI512Ops::Load(source + processed), factor,
                                affine));
      }
      break;
    }
#else
      return Status::unsupported_backend;
#endif
    default:
      return Status::invalid_argument;
  }
  ScalarMultiplyTail(destination, source, processed, byte_count, coefficient,
                     tables);
  return Status::ok;
}

Status GFNIIFFTRadix2Copy(const Element* x,
                          const Element* y,
                          Element* output_x,
                          Element* output_y,
                          size_t byte_count,
                          Element coefficient,
                          Backend backend,
                          const MultiplicationTables& tables) {
  size_t processed = 0;
  [[maybe_unused]] constexpr bool affine = true;
  switch (backend) {
    case Backend::gfni128_affine:
#if defined(__GFNI__) && defined(__SSE2__)
      processed = IFFTRadix2CopyVectors<GFNI128Ops>(
          x, y, output_x, output_y, byte_count, coefficient, affine, tables);
      break;
#else
      return Status::unsupported_backend;
#endif
    case Backend::gfni256_affine:
#if defined(__GFNI__) && defined(__AVX2__)
      processed = IFFTRadix2CopyVectors<GFNI256Ops>(
          x, y, output_x, output_y, byte_count, coefficient, affine, tables);
      break;
#else
      return Status::unsupported_backend;
#endif
    case Backend::gfni512_affine:
#if defined(__GFNI__) && defined(__AVX512F__) && defined(__AVX512BW__)
      processed = IFFTRadix2CopyVectors<GFNI512Ops>(
          x, y, output_x, output_y, byte_count, coefficient, affine, tables);
      break;
#else
      return Status::unsupported_backend;
#endif
    default:
      return Status::invalid_argument;
  }
  for (; processed < byte_count; ++processed) {
    const Element b = y[processed] ^ x[processed];
    const Element a = x[processed] ^ ScalarProduct(b, coefficient, tables);
    output_x[processed] = a;
    output_y[processed] = b;
  }
  return Status::ok;
}

Status GFNIIFFTRadix2Xor(const Element* x,
                         const Element* y,
                         Element* output_x,
                         Element* output_y,
                         size_t byte_count,
                         Element coefficient,
                         Backend backend,
                         const MultiplicationTables& tables) {
  size_t processed = 0;
  [[maybe_unused]] constexpr bool affine = true;
  switch (backend) {
    case Backend::gfni128_affine:
#if defined(__GFNI__) && defined(__SSE2__)
      processed = IFFTRadix2XorVectors<GFNI128Ops>(
          x, y, output_x, output_y, byte_count, coefficient, affine, tables);
      break;
#else
      return Status::unsupported_backend;
#endif
    case Backend::gfni256_affine:
#if defined(__GFNI__) && defined(__AVX2__)
      processed = IFFTRadix2XorVectors<GFNI256Ops>(
          x, y, output_x, output_y, byte_count, coefficient, affine, tables);
      break;
#else
      return Status::unsupported_backend;
#endif
    case Backend::gfni512_affine:
#if defined(__GFNI__) && defined(__AVX512F__) && defined(__AVX512BW__)
      processed = IFFTRadix2XorVectors<GFNI512Ops>(
          x, y, output_x, output_y, byte_count, coefficient, affine, tables);
      break;
#else
      return Status::unsupported_backend;
#endif
    default:
      return Status::invalid_argument;
  }
  for (; processed < byte_count; ++processed) {
    const Element b = y[processed] ^ x[processed];
    const Element a = x[processed] ^ ScalarProduct(b, coefficient, tables);
    output_x[processed] ^= a;
    output_y[processed] ^= b;
  }
  return Status::ok;
}

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
                          const MultiplicationTables& tables) {
  size_t processed = 0;
  [[maybe_unused]] constexpr bool affine = true;
  switch (backend) {
    case Backend::gfni128_affine:
#if defined(__GFNI__) && defined(__SSE2__)
      processed = IFFTRadix4CopyVectors<GFNI128Ops>(
          x0, x1, x2, x3, output0, output1, output2, output3, byte_count, top,
          low, high, affine, tables);
      break;
#else
      return Status::unsupported_backend;
#endif
    case Backend::gfni256_affine:
#if defined(__GFNI__) && defined(__AVX2__)
      processed = IFFTRadix4CopyVectors<GFNI256Ops>(
          x0, x1, x2, x3, output0, output1, output2, output3, byte_count, top,
          low, high, affine, tables);
      break;
#else
      return Status::unsupported_backend;
#endif
    case Backend::gfni512_affine:
#if defined(__GFNI__) && defined(__AVX512F__) && defined(__AVX512BW__)
      processed = IFFTRadix4CopyVectors<GFNI512Ops>(
          x0, x1, x2, x3, output0, output1, output2, output3, byte_count, top,
          low, high, affine, tables);
      break;
#else
      return Status::unsupported_backend;
#endif
    default:
      return Status::invalid_argument;
  }
  for (; processed < byte_count; ++processed) {
    Element a = x0[processed];
    Element b = x1[processed] ^ a;
    a ^= ScalarProduct(b, low, tables);
    Element c = x2[processed];
    Element d = x3[processed] ^ c;
    c ^= ScalarProduct(d, high, tables);
    c ^= a;
    a ^= ScalarProduct(c, top, tables);
    d ^= b;
    b ^= ScalarProduct(d, top, tables);
    output0[processed] = a;
    output1[processed] = b;
    output2[processed] = c;
    output3[processed] = d;
  }
  return Status::ok;
}

Status GFNIFFTRadix4(Element* x0,
                     Element* x1,
                     Element* x2,
                     Element* x3,
                     size_t byte_count,
                     Element top,
                     Element low,
                     Element high,
                     Backend backend,
                     const MultiplicationTables& tables) {
  size_t processed = 0;
  [[maybe_unused]] constexpr bool affine = true;
  switch (backend) {
    case Backend::gfni128_affine:
#if defined(__GFNI__) && defined(__SSE2__)
      processed = FFTVectors<GFNI128Ops>(x0, x1, x2, x3, byte_count, top, low,
                                         high, affine, tables);
      break;
#else
      return Status::unsupported_backend;
#endif
    case Backend::gfni256_affine:
#if defined(__GFNI__) && defined(__AVX2__)
      processed = FFTVectors<GFNI256Ops>(x0, x1, x2, x3, byte_count, top, low,
                                         high, affine, tables);
      break;
#else
      return Status::unsupported_backend;
#endif
    case Backend::gfni512_affine:
#if defined(__GFNI__) && defined(__AVX512F__) && defined(__AVX512BW__)
      processed = FFTVectors<GFNI512Ops>(x0, x1, x2, x3, byte_count, top, low,
                                         high, affine, tables);
      break;
#else
      return Status::unsupported_backend;
#endif
    default:
      return Status::invalid_argument;
  }
  FFTScalarTail(x0, x1, x2, x3, processed, byte_count, top, low, high, tables);
  return Status::ok;
}

Status GFNIIFFTRadix4(Element* x0,
                      Element* x1,
                      Element* x2,
                      Element* x3,
                      size_t byte_count,
                      Element top,
                      Element low,
                      Element high,
                      Backend backend,
                      const MultiplicationTables& tables) {
  size_t processed = 0;
  [[maybe_unused]] constexpr bool affine = true;
  switch (backend) {
    case Backend::gfni128_affine:
#if defined(__GFNI__) && defined(__SSE2__)
      processed = IFFTVectors<GFNI128Ops>(x0, x1, x2, x3, byte_count, top, low,
                                          high, affine, tables);
      break;
#else
      return Status::unsupported_backend;
#endif
    case Backend::gfni256_affine:
#if defined(__GFNI__) && defined(__AVX2__)
      processed = IFFTVectors<GFNI256Ops>(x0, x1, x2, x3, byte_count, top, low,
                                          high, affine, tables);
      break;
#else
      return Status::unsupported_backend;
#endif
    case Backend::gfni512_affine:
#if defined(__GFNI__) && defined(__AVX512F__) && defined(__AVX512BW__)
      processed = IFFTVectors<GFNI512Ops>(x0, x1, x2, x3, byte_count, top, low,
                                          high, affine, tables);
      break;
#else
      return Status::unsupported_backend;
#endif
    default:
      return Status::invalid_argument;
  }
  IFFTScalarTail(x0, x1, x2, x3, processed, byte_count, top, low, high, tables);
  return Status::ok;
}

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
                         const MultiplicationTables& tables) {
  size_t processed = 0;
  [[maybe_unused]] constexpr bool affine = true;
  switch (backend) {
    case Backend::gfni128_affine:
#if defined(__GFNI__) && defined(__SSE2__)
      processed = IFFTRadix4XorVectors<GFNI128Ops>(
          x0, x1, x2, x3, output0, output1, output2, output3, byte_count, top,
          low, high, affine, tables);
      break;
#else
      return Status::unsupported_backend;
#endif
    case Backend::gfni256_affine:
#if defined(__GFNI__) && defined(__AVX2__)
      processed = IFFTRadix4XorVectors<GFNI256Ops>(
          x0, x1, x2, x3, output0, output1, output2, output3, byte_count, top,
          low, high, affine, tables);
      break;
#else
      return Status::unsupported_backend;
#endif
    case Backend::gfni512_affine:
#if defined(__GFNI__) && defined(__AVX512F__) && defined(__AVX512BW__)
      processed = IFFTRadix4XorVectors<GFNI512Ops>(
          x0, x1, x2, x3, output0, output1, output2, output3, byte_count, top,
          low, high, affine, tables);
      break;
#else
      return Status::unsupported_backend;
#endif
    default:
      return Status::invalid_argument;
  }
  for (; processed < byte_count; ++processed) {
    Element a = x0[processed];
    Element b = x1[processed] ^ a;
    a ^= ScalarProduct(b, low, tables);
    Element c = x2[processed];
    Element d = x3[processed] ^ c;
    c ^= ScalarProduct(d, high, tables);
    c ^= a;
    a ^= ScalarProduct(c, top, tables);
    d ^= b;
    b ^= ScalarProduct(d, top, tables);
    output0[processed] ^= a;
    output1[processed] ^= b;
    output2[processed] ^= c;
    output3[processed] ^= d;
  }
  return Status::ok;
}

}  // namespace gf2p8::lch::detail
