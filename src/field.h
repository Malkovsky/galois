#pragma once

#include <cstdint>

/**
 * GF(256) field implementation using polynomial representation
 * with the irreducible polynomial x^8 + x^4 + x^3 + x + 1 (0x11B)
 */
namespace gf_2_8 {

typedef uint8_t element_t;

/**
 * Field zero element, 0 for most implementations
 */
element_t Zero();

/**
 * Field unity element, 1 for most implementations
 */
element_t One();

/**
 * Initialize the GF(256) field tables (logarithm and exponentiation)
 * Must be called before using other field operations
 */
void Init();

/**
 * Add two elements in GF(256)
 * In GF(256), addition is XOR
 * @param a First element
 * @param b Second element
 * @return The sum a + b in GF(256)
 */
element_t Add(element_t a, element_t b);

/**
 * Subtract two elements in GF(256)
 * In GF(256), subtraction is the same as addition (XOR)
 * @param a First element
 * @param b Second element
 * @return The difference a - b in GF(256)
 */
element_t Sub(element_t a, element_t b);

/**
 * Multiply two elements in GF(256)
 * @param a First element
 * @param b Second element
 * @return The product a * b in GF(256)
 */
element_t Multiply(element_t a, element_t b);

/**
 * Multiply two elements in GF(256) using exp/log LUT
 * @param a First element
 * @param b Second element
 * @return The product a * b in GF(256)
 */
element_t MultiplyLUT(element_t a, element_t b);

/**
 * Divide two elements in GF(256)
 * @param a First element
 * @param b Second element (must be non-zero)
 * @return The quotient a / b in GF(256)
 */
element_t Div(element_t a, element_t b);

/**
 * Calculate the inverse of an element in GF(256)
 * @param a The element to invert (must be non-zero)
 * @return The inverse of a in GF(256)
 */
element_t Inv(element_t a);

/**
 * Exponentiate an element in GF(256)
 * @param a The base element
 * @param n The exponent
 * @return a^n in GF(256)
 */
element_t Pow(element_t a, int n);

} // namespace gf_2_8


/**
 * GF(2^16) field implementation using as an extension over GF(2^8)
 * with the irreducible polynomial x^2+x+32
 */
namespace gf_2_16 {

typedef uint16_t element_t;

/**
 * Field zero element, 0 for most implementations
 */
element_t Zero();

/**
 * Field unity element, 1 for most implementations
 */
element_t One();

/**
 * Add two elements in GF(2^16)
 * In GF(2^16), addition is XOR
 * @param a First element
 * @param b Second element
 * @return The sum a + b in GF(2^16)
 */
element_t Add(element_t a, element_t b);

/**
 * Subtract two elements in GF(2^16)
 * In GF(2^16), subtraction is the same as addition (XOR)
 * @param a First element
 * @param b Second element
 * @return The difference a - b in GF(2^16)
 */
element_t Sub(element_t a, element_t b);

/**
 * Multiply two elements in GF(2^16)
 * @param a First element
 * @param b Second element
 * @return The product a * b in GF(2^16)
 */
element_t Multiply(element_t a, element_t b);

/**
 * Divide two elements in GF(2^16)
 * @param a First element
 * @param b Second element (must be non-zero)
 * @return The quotient a / b in GF(2^16)
 */
element_t Div(element_t a, element_t b);

/**
 * Calculate the inverse of an element in GF(2^16)
 * @param a The element to invert (must be non-zero)
 * @return The inverse of a in GF(2^16)
 */
element_t Inv(element_t a);

/**
 * Calculate the inverse of an element in GF(2^16)
 * using Itoh-Tsujii algorithm
 * @param a The element to invert (must be non-zero)
 * @return The inverse of a in GF(2^16)
 */
element_t InvIT(element_t a);

/**
 * Exponentiate an element in GF(2^16)
 * @param a The base element
 * @param n The exponent
 * @return a^n in GF(2^16)
 */
element_t Pow(element_t a, size_t n);

} // namespace gf_2_16

