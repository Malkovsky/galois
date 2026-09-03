#pragma once

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <utility>

/**
 * @brief GF(256) field operations.
 * @details Element is basis-neutral byte storage. Basis-specific operations
 * explicitly state whether bytes use standard polynomial coordinates or
 * Cantor coordinates over GF(2^8)/0x11D. Addition is byte XOR in either basis.
 * The Cantor representation is directly compatible with Leopard FF8.
 */
namespace gf2p8 {

using Element = uint8_t;

/** @cond */
namespace detail {

inline constexpr std::array<Element, 8> kCantorPolynomialBasis = {
    0x01, 0xd6, 0x98, 0x92, 0x56, 0xc8, 0x58, 0xe6};

consteval std::array<Element, 8> PolynomialToCantorRows() {
  std::array<uint16_t, 8> rows{};
  for (size_t polynomial_bit = 0; polynomial_bit < 8; ++polynomial_bit) {
    uint16_t row = static_cast<uint16_t>(1U << (8 + polynomial_bit));
    for (size_t coordinate_bit = 0; coordinate_bit < 8; ++coordinate_bit) {
      row |= static_cast<uint16_t>(
          ((kCantorPolynomialBasis[coordinate_bit] >> polynomial_bit) & 1U)
          << coordinate_bit);
    }
    rows[polynomial_bit] = row;
  }
  for (size_t column = 0; column < 8; ++column) {
    size_t pivot = column;
    while (((rows[pivot] >> column) & 1U) == 0) {
      ++pivot;
    }
    std::swap(rows[column], rows[pivot]);
    for (size_t row = 0; row < 8; ++row) {
      if (row != column && ((rows[row] >> column) & 1U) != 0) {
        rows[row] ^= rows[column];
      }
    }
  }
  std::array<Element, 8> inverse{};
  for (size_t row = 0; row < 8; ++row) {
    inverse[row] = static_cast<Element>(rows[row] >> 8);
  }
  return inverse;
}

inline constexpr auto kPolynomialToCantorRows = PolynomialToCantorRows();

constexpr Element CantorToStandardDirect(Element value) {
  Element result = 0;
  for (size_t bit = 0; bit < 8; ++bit) {
    if (((value >> bit) & 1U) != 0) {
      result ^= kCantorPolynomialBasis[bit];
    }
  }
  return result;
}

constexpr Element StandardToCantorDirect(Element value) {
  Element result = 0;
  for (size_t bit = 0; bit < 8; ++bit) {
    result |= static_cast<Element>((std::popcount(static_cast<unsigned>(
                                        value & kPolynomialToCantorRows[bit])) &
                                    1U)
                                   << bit);
  }
  return result;
}

constexpr Element MultiplyStandardDirect(Element a, Element b) {
  Element result = 0;
  while (a != 0) {
    result ^= static_cast<Element>(b * (a & 1U));
    a >>= 1;
    b = static_cast<Element>((b << 1) ^ (0x1dU * (b >> 7)));
  }
  return result;
}

constexpr Element MultiplyCantorDirect(Element a, Element b) {
  return StandardToCantorDirect(MultiplyStandardDirect(
      CantorToStandardDirect(a), CantorToStandardDirect(b)));
}

}  // namespace detail
/** @endcond */

/**
 * @brief Stores logarithm tables for one GF(256) coordinate basis.
 * @details Entries are indexed by bytes interpreted in the table's documented
 * coordinate basis.
 */
struct LogarithmTables {
  std::array<Element, 256> logarithm{};
  std::array<Element, 256> exponent{};
};

/**
 * @brief Stores the shared immutable GF(256) multiplication tables.
 * @details Standard and Cantor log tables support basis-explicit scalar
 * arithmetic. Each shuffle row contains lane-duplicated Cantor-coordinate
 * low-nibble products at offsets 0 and 16 and high-nibble products at offsets
 * 32 and 48. Affine rows encode fixed-coefficient Cantor multiplication for
 * GFNI.
 */
struct MultiplicationTables {
  LogarithmTables standard{};
  LogarithmTables cantor{};
  alignas(32) std::array<std::array<Element, 64>, 256> shuffle{};
  std::array<uint64_t, 256> affine{};
};

/**
 * @brief Returns the process-wide immutable multiplication tables.
 * @details Tables are generated at compile time and require no initialization.
 * @return The shared GF(256) table set.
 */
const MultiplicationTables& Tables();

/**
 * @brief Returns the GF(256) additive identity.
 * @details Both supported coordinate bases use the byte value zero.
 * @return The field zero element.
 */
constexpr Element Zero() {
  return 0;
}

/**
 * @brief Returns the GF(256) multiplicative identity.
 * @details Both supported coordinate bases use the byte value one.
 * @return The field unity element.
 */
constexpr Element One() {
  return 1;
}

/**
 * @brief Adds two GF(256) elements.
 * @details Field addition is bitwise XOR.
 * @param a First element.
 * @param b Second element.
 * @return The sum @p a + @p b.
 */
constexpr Element Add(Element a, Element b) {
  return a ^ b;
}

/**
 * @brief Subtracts two GF(256) elements.
 * @details In characteristic two, subtraction is the same as addition and is
 * implemented as bitwise XOR.
 * @param a First element.
 * @param b Second element.
 * @return The difference @p a - @p b.
 */
constexpr Element Sub(Element a, Element b) {
  return a ^ b;
}

/**
 * @brief Converts a standard-coordinate byte to Cantor coordinates.
 * @details The input bits are coefficients of `{1, x, ..., x^7}` modulo
 * `0x11D`. The output bits are coordinates in the Cantor basis
 * `{01, D6, 98, 92, 56, C8, 58, E6}`.
 * @param value Standard polynomial-coordinate byte.
 * @return Equivalent Cantor-coordinate byte.
 */
constexpr Element StandardToCantor(Element value) {
  return detail::StandardToCantorDirect(value);
}

/**
 * @brief Converts a Cantor-coordinate byte to standard coordinates.
 * @details The input bits are coordinates in the Cantor basis
 * `{01, D6, 98, 92, 56, C8, 58, E6}`. The output bits are coefficients of
 * `{1, x, ..., x^7}` modulo `0x11D`.
 * @param value Cantor-coordinate byte.
 * @return Equivalent standard polynomial-coordinate byte.
 */
constexpr Element CantorToStandard(Element value) {
  return detail::CantorToStandardDirect(value);
}

/**
 * @brief Multiplies two standard-coordinate GF(256) elements.
 * @details Both operands and the result use polynomial coefficients
 * `{1, x, ..., x^7}` modulo `0x11D`. Passing Cantor-coordinate bytes without
 * conversion is invalid. Uses the immutable standard-coordinate log tables.
 * @param a First element.
 * @param b Second element.
 * @return The product @p a * @p b.
 */
Element MultiplyStandard(Element a, Element b);

/**
 * @brief Multiplies two Cantor-coordinate GF(256) elements.
 * @details Both operands and the result use Cantor coordinates.
 * Passing standard-coordinate bytes without conversion is invalid. Uses the
 * immutable Cantor-coordinate log tables.
 * @param a First element.
 * @param b Second element.
 * @return The product @p a * @p b.
 */
Element MultiplyCantor(Element a, Element b);

/**
 * @brief Divides two standard-coordinate GF(256) elements.
 * @details Both operands and the result use polynomial coefficients modulo
 * `0x11D`. Computes multiplication by the inverse of @p b.
 * @param a Dividend.
 * @param b Nonzero divisor.
 * @return The quotient @p a / @p b.
 */
Element DivStandard(Element a, Element b);

/**
 * @brief Divides two Cantor-coordinate GF(256) elements.
 * @details Both operands and the result use Cantor coordinates.
 * Computes multiplication by the inverse of @p b.
 * @param a Dividend.
 * @param b Nonzero divisor.
 * @return The quotient @p a / @p b.
 */
Element DivCantor(Element a, Element b);

/**
 * @brief Computes a standard-coordinate GF(256) multiplicative inverse.
 * @details The input and result use polynomial coefficients modulo `0x11D`.
 * The input must be nonzero.
 * @param a Element to invert.
 * @return The inverse of @p a.
 */
Element InvStandard(Element a);

/**
 * @brief Computes a Cantor-coordinate GF(256) inverse.
 * @details The input and result use Cantor coordinates. The input must be
 * nonzero.
 * @param a Element to invert.
 * @return The inverse of @p a.
 */
Element InvCantor(Element a);

/**
 * @brief Raises a standard-coordinate GF(256) element to an integer power.
 * @details The base and result use polynomial coefficients modulo `0x11D`.
 * By convention, `0^0` returns one and every other power of zero returns zero.
 * @param a Base element.
 * @param n Exponent.
 * @return @p a raised to @p n.
 */
Element PowStandard(Element a, int n);

/**
 * @brief Raises a Cantor-coordinate GF(256) element to a power.
 * @details The base and result use Cantor coordinates.
 * By convention, `0^0` returns one and every other power of zero returns zero.
 * @param a Base element.
 * @param n Exponent.
 * @return @p a raised to @p n.
 */
Element PowCantor(Element a, int n);

/**
 * @brief Adds a scaled row using the best compiled shuffle-table kernel.
 * @details Computes `destination ^= coefficient * source` with exact scalar
 * tails. Source, coefficient, and destination use Cantor coordinates.
 * The operation uses the shared immutable table set.
 * @param destination Destination row.
 * @param source Source row.
 * @param coefficient Field scaling coefficient.
 * @param length Number of elements in each row.
 */
void AddScaledRow(Element* destination,
                  const Element* source,
                  Element coefficient,
                  size_t length);

}  // namespace gf2p8

/**
 * @brief GF(2^16) field operations.
 * @details Uses a quadratic extension over GF(2^8) with irreducible polynomial
 * `x^2 + x + delta`, where delta is polynomial element `0x20` and native
 * Cantor-coordinate byte `0xF0`.
 */
namespace gf2p16 {

using Element = uint16_t;

/**
 * @brief Returns the GF(2^16) additive identity.
 * @details The representation uses the value zero.
 * @return The field zero element.
 */
Element Zero();

/**
 * @brief Returns the GF(2^16) multiplicative identity.
 * @details The representation uses the value one.
 * @return The field unity element.
 */
Element One();

/**
 * @brief Adds two GF(2^16) elements.
 * @details Field addition is bitwise XOR.
 * @param a First element.
 * @param b Second element.
 * @return The sum @p a + @p b.
 */
Element Add(Element a, Element b);

/**
 * @brief Subtracts two GF(2^16) elements.
 * @details In characteristic two, subtraction is bitwise XOR.
 * @param a First element.
 * @param b Second element.
 * @return The difference @p a - @p b.
 */
Element Sub(Element a, Element b);

/**
 * @brief Multiplies two GF(2^16) elements.
 * @details Uses the quadratic-extension representation.
 * @param a First element.
 * @param b Second element.
 * @return The product @p a * @p b.
 */
Element Multiply(Element a, Element b);

/**
 * @brief Divides two GF(2^16) elements.
 * @details Computes multiplication by the inverse of @p b.
 * @param a Dividend.
 * @param b Nonzero divisor.
 * @return The quotient @p a / @p b.
 */
Element Div(Element a, Element b);

/**
 * @brief Computes a GF(2^16) multiplicative inverse.
 * @details The input must be nonzero.
 * @param a Element to invert.
 * @return The inverse of @p a.
 */
Element Inv(Element a);

/**
 * @brief Computes a GF(2^16) inverse with Itoh-Tsujii.
 * @details Uses the Itoh-Tsujii inversion algorithm; the input must be nonzero.
 * @param a Element to invert.
 * @return The inverse of @p a.
 */
Element InvIT(Element a);

/**
 * @brief Raises a GF(2^16) element to a power.
 * @details Exponentiation follows GF(2^16) multiplication semantics.
 * @param a Base element.
 * @param n Exponent.
 * @return @p a raised to @p n.
 */
Element Pow(Element a, size_t n);

}  // namespace gf2p16
