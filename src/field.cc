#include "field.h"

#include <immintrin.h>

namespace gf2p8 {
namespace {

template <typename MultiplyFunction>
consteval LogarithmTables MakeLogarithmTables(Element primitive,
                                              MultiplyFunction multiply) {
  LogarithmTables tables;
  tables.logarithm.fill(255);
  tables.exponent[0] = 1;
  for (size_t i = 1; i < 255; ++i) {
    tables.exponent[i] = multiply(tables.exponent[i - 1], primitive);
  }
  tables.exponent[255] = 1;
  for (size_t i = 0; i < 255; ++i) {
    tables.logarithm[tables.exponent[i]] = static_cast<Element>(i);
  }
  return tables;
}

consteval uint64_t AffineMatrix(Element coefficient) {
  uint64_t columns = 0;
  for (size_t bit = 0; bit < 8; ++bit) {
    columns |= static_cast<uint64_t>(detail::MultiplyCantorDirect(
                   coefficient, static_cast<Element>(1U << bit)))
               << (8 * bit);
  }

  uint64_t matrix = 0;
  for (size_t row = 0; row < 8; ++row) {
    for (size_t column = 0; column < 8; ++column) {
      const uint64_t bit = (columns >> (8 * column + row)) & 1U;
      matrix |= bit << (8 * (7 - row) + column);
    }
  }
  return matrix;
}

consteval MultiplicationTables MakeTables() {
  MultiplicationTables tables;
  tables.standard = MakeLogarithmTables(2, [](Element a, Element b) {
    return detail::MultiplyStandardDirect(a, b);
  });
  tables.cantor = MakeLogarithmTables(
      StandardToCantor(2),
      [](Element a, Element b) { return detail::MultiplyCantorDirect(a, b); });

  for (size_t coefficient = 0; coefficient < 256; ++coefficient) {
    for (size_t nibble = 0; nibble < 16; ++nibble) {
      const Element low = detail::MultiplyCantorDirect(
          static_cast<Element>(coefficient), static_cast<Element>(nibble));
      const Element high = detail::MultiplyCantorDirect(
          static_cast<Element>(coefficient), static_cast<Element>(nibble << 4));
      tables.shuffle[coefficient][nibble] = low;
      tables.shuffle[coefficient][16 + nibble] = low;
      tables.shuffle[coefficient][32 + nibble] = high;
      tables.shuffle[coefficient][48 + nibble] = high;
    }
    tables.affine[coefficient] =
        AffineMatrix(static_cast<Element>(coefficient));
  }
  return tables;
}

alignas(32) constinit const MultiplicationTables kTables = MakeTables();

Element Product(const Element* row, Element value) {
  return row[value & 0x0f] ^ row[32 + (value >> 4)];
}

Element MultiplyWithLogs(Element a, Element b, const LogarithmTables& tables) {
  if (a == 0 || b == 0) {
    return 0;
  }
  const unsigned sum = tables.logarithm[a] + tables.logarithm[b];
  return tables.exponent[sum >= 255 ? sum - 255 : sum];
}

Element DivWithLogs(Element a, Element b, const LogarithmTables& tables) {
  if (a == 0) {
    return 0;
  }
  int log_diff = tables.logarithm[a] - tables.logarithm[b];
  if (log_diff < 0) {
    log_diff += 255;
  }
  return tables.exponent[log_diff];
}

Element InvWithLogs(Element a, const LogarithmTables& tables) {
  if (a == 0) {
    return 0;
  }
  return tables.exponent[255 - tables.logarithm[a]];
}

Element PowWithLogs(Element a, int n, const LogarithmTables& tables) {
  if (a == 0) {
    return n == 0 ? 1 : 0;
  }

  int64_t exponent = n;
  if (exponent < 0) {
    a = InvWithLogs(a, tables);
    exponent = -exponent;
  }
  const auto index =
      static_cast<size_t>((static_cast<uint64_t>(tables.logarithm[a]) *
                           static_cast<uint64_t>(exponent)) %
                          255);
  return tables.exponent[index];
}

}  // namespace

const MultiplicationTables& Tables() {
  return kTables;
}

Element MultiplyStandard(Element a, Element b) {
  return MultiplyWithLogs(a, b, Tables().standard);
}

Element MultiplyCantor(Element a, Element b) {
  return MultiplyWithLogs(a, b, Tables().cantor);
}

Element DivStandard(Element a, Element b) {
  return DivWithLogs(a, b, Tables().standard);
}

Element DivCantor(Element a, Element b) {
  return DivWithLogs(a, b, Tables().cantor);
}

Element InvStandard(Element a) {
  return InvWithLogs(a, Tables().standard);
}

Element InvCantor(Element a) {
  return InvWithLogs(a, Tables().cantor);
}

Element PowStandard(Element a, int n) {
  return PowWithLogs(a, n, Tables().standard);
}

Element PowCantor(Element a, int n) {
  return PowWithLogs(a, n, Tables().cantor);
}

void AddScaledRow(Element* destination,
                  const Element* source,
                  Element coefficient,
                  size_t length) {
  if (coefficient == 0) {
    return;
  }
  const Element* const row = Tables().shuffle[coefficient].data();
  size_t processed = 0;
#if defined(__AVX2__)
  const __m256i mask = _mm256_set1_epi8(0x0f);
  const __m256i table_lo =
      _mm256_load_si256(reinterpret_cast<const __m256i*>(row));
  const __m256i table_hi =
      _mm256_load_si256(reinterpret_cast<const __m256i*>(row + 32));
  while (processed + 32 <= length) {
    const __m256i source_vector = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(source + processed));
    const __m256i lo =
        _mm256_shuffle_epi8(table_lo, _mm256_and_si256(source_vector, mask));
    const __m256i hi = _mm256_shuffle_epi8(
        table_hi, _mm256_and_si256(_mm256_srli_epi64(source_vector, 4), mask));
    const __m256i destination_vector = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(destination + processed));
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(destination + processed),
        _mm256_xor_si256(destination_vector, _mm256_xor_si256(lo, hi)));
    processed += 32;
  }
#elif defined(__SSSE3__)
  const __m128i mask = _mm_set1_epi8(0x0f);
  const __m128i table_lo =
      _mm_load_si128(reinterpret_cast<const __m128i*>(row));
  const __m128i table_hi =
      _mm_load_si128(reinterpret_cast<const __m128i*>(row + 32));
  while (processed + 16 <= length) {
    const __m128i source_vector =
        _mm_loadu_si128(reinterpret_cast<const __m128i*>(source + processed));
    const __m128i lo =
        _mm_shuffle_epi8(table_lo, _mm_and_si128(source_vector, mask));
    const __m128i hi = _mm_shuffle_epi8(
        table_hi, _mm_and_si128(_mm_srli_epi64(source_vector, 4), mask));
    const __m128i destination_vector = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(destination + processed));
    _mm_storeu_si128(reinterpret_cast<__m128i*>(destination + processed),
                     _mm_xor_si128(destination_vector, _mm_xor_si128(lo, hi)));
    processed += 16;
  }
#endif
  for (; processed < length; ++processed) {
    destination[processed] ^= Product(row, source[processed]);
  }
}

}  // namespace gf2p8

namespace gf2p16 {

constexpr gf2p8::Element delta = 0xf0;

Element Zero() {
  return 0;
}

Element One() {
  return 1;
}

Element Add(Element a, Element b) {
  return a ^ b;
}

Element Sub(Element a, Element b) {
  return a ^ b;
}

Element Multiply(Element a, Element b) {
  // a = a_0 + a_1x, b = b_0 + b_1x
  // all four are from GF(256)
  gf2p8::Element a_0 = a & 255;
  gf2p8::Element a_1 = a >> 8;
  gf2p8::Element b_0 = b & 255;
  gf2p8::Element b_1 = b >> 8;

  auto t = gf2p8::MultiplyCantor(a_1, b_1);
  auto low_bits = gf2p8::Add(gf2p8::MultiplyCantor(a_0, b_0),
                             gf2p8::MultiplyCantor(t, delta));
  auto high_bits = gf2p8::Add(gf2p8::Add(gf2p8::MultiplyCantor(a_0, b_1),
                                         gf2p8::MultiplyCantor(a_1, b_0)),
                              t);
  return low_bits + (high_bits << 8);
}

Element Inv(Element a) {
  Element result = One();
  Element b = Multiply(a, a);
  for (size_t i = 1; i < 16; ++i) {
    result = Multiply(result, b);
    b = Multiply(b, b);
  }
  return result;
}

Element Div(Element a, Element b) {
  return Multiply(a, Inv(b));
}

Element Pow(Element a, size_t n) {
  Element result = One();
  while (n) {
    if (n & 1) {
      result = Multiply(result, a);
    }
    a = Multiply(a, a);
    n >>= 1;
  }
  return result;
}

Element InvIT(Element a) {
  Element a_r = a;
  for (size_t i = 0; i < 8; ++i) {
    a_r = Multiply(a_r, a_r);
  }
  gf2p8::Element a_r1 = Multiply(a_r, a);
  return Multiply(a_r, gf2p8::InvCantor(a_r1));
}

}  // namespace gf2p16
