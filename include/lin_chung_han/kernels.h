#pragma once

#include <cstddef>
#include <cstdint>

#include "field.h"

/**
 * @brief LCH transform kernels over GF(256).
 * @details All field-valued bytes, coefficients, shards, and multiplication
 * tables in this namespace use Cantor coordinates. Callers must pass the shared
 * Cantor kernel tables returned by `gf2p8::Tables()`.
 */
namespace gf2p8::lch {

/**
 * @brief Reports LCH transform, kernel, and codec operation status.
 * @details Operations return explicit status values instead of throwing.
 */
enum class Status {
  ok,
  invalid_argument,
  unsupported_backend,
  insufficient_recovery_symbols,
};

/**
 * @brief Selects an LCH byte-kernel implementation.
 * @details Explicit backends report @ref Status::unsupported_backend when not
 * compiled. The tuned backend is selected from the shard byte count.
 */
enum class Backend {
  tuned,
  scalar,
  ssse3,
  avx2,
  gfni128_affine,
  gfni256_affine,
  gfni512_affine,
};

/**
 * @brief Selects radix-two or fused radix-four transform scheduling.
 * @details Radix four combines adjacent radix-two layers while preserving the
 * same transform semantics.
 */
enum class Radix { radix2, radix4 };

/**
 * @brief Tests whether a concrete backend is compiled.
 * @details The tuned pseudo-backend is not a concrete executable backend.
 * @param backend Backend to query.
 * @return True when the concrete backend can execute.
 */
bool BackendAvailable(Backend backend);

/**
 * @brief Selects the tuned backend for a shard byte count.
 * @details Selection uses compile-time ISA availability and measured static
 * byte-count thresholds.
 * @param byte_count Bytes in one shard.
 * @return Selected concrete backend.
 */
Backend SelectBackend(size_t byte_count);

/**
 * @brief Adds a scaled source byte vector to a destination.
 * @details Computes `destination ^= coefficient * source` in Cantor
 * coordinates, with an exact scalar tail for arbitrary byte counts.
 * @param destination Destination byte vector.
 * @param source Source byte vector.
 * @param byte_count Number of bytes to process.
 * @param coefficient Field scaling coefficient.
 * @param backend Kernel backend.
 * @param tables Shared Cantor kernel tables returned by `gf2p8::Tables()`.
 * @return Kernel status.
 */
Status AddScaled(Element* destination,
                 const Element* source,
                 size_t byte_count,
                 Element coefficient,
                 Backend backend,
                 const MultiplicationTables& tables);

/**
 * @brief Overwrites a destination with a scaled source byte vector.
 * @details Computes `destination = coefficient * source` in Cantor
 * coordinates, with an exact scalar tail for arbitrary byte counts.
 * @param destination Destination byte vector.
 * @param source Source byte vector.
 * @param byte_count Number of bytes to process.
 * @param coefficient Field scaling coefficient.
 * @param backend Kernel backend.
 * @param tables Shared Cantor kernel tables returned by `gf2p8::Tables()`.
 * @return Kernel status.
 */
Status Scale(Element* destination,
             const Element* source,
             size_t byte_count,
             Element coefficient,
             Backend backend,
             const MultiplicationTables& tables);

/**
 * @brief XORs one source byte vector into a destination.
 * @details Processes arbitrary byte counts with exact scalar tails.
 * @param destination Destination byte vector.
 * @param source Source byte vector.
 * @param byte_count Number of bytes to process.
 * @param backend Kernel backend.
 * @return Kernel status.
 */
Status Xor(Element* destination,
           const Element* source,
           size_t byte_count,
           Backend backend);

/**
 * @brief XORs four independent source vectors into four destinations.
 * @details Batches four streams to improve instruction-level parallelism while
 * preserving exact tails.
 * @param destination0 First destination byte vector.
 * @param source0 First source byte vector.
 * @param destination1 Second destination byte vector.
 * @param source1 Second source byte vector.
 * @param destination2 Third destination byte vector.
 * @param source2 Third source byte vector.
 * @param destination3 Fourth destination byte vector.
 * @param source3 Fourth source byte vector.
 * @param byte_count Number of bytes in each vector.
 * @param backend Kernel backend.
 * @return Kernel status.
 */
Status Xor4(Element* destination0,
            const Element* source0,
            Element* destination1,
            const Element* source1,
            Element* destination2,
            const Element* source2,
            Element* destination3,
            const Element* source3,
            size_t byte_count,
            Backend backend);

/**
 * @brief Applies one in-place radix-two FFT butterfly.
 * @details Uses Cantor-coordinate shards, skew coefficient, and shared
 * multiplication tables.
 * @param x First butterfly shard.
 * @param y Second butterfly shard.
 * @param byte_count Bytes in each shard.
 * @param coefficient LCH skew coefficient.
 * @param backend Kernel backend.
 * @param tables Shared Cantor kernel tables returned by `gf2p8::Tables()`.
 * @return Kernel status.
 */
Status FFTRadix2(Element* x,
                 Element* y,
                 size_t byte_count,
                 Element coefficient,
                 Backend backend,
                 const MultiplicationTables& tables);

/**
 * @brief Applies one in-place radix-two IFFT butterfly.
 * @details Reverses the normalized radix-two FFT butterfly in Cantor
 * coordinates.
 * @param x First butterfly shard.
 * @param y Second butterfly shard.
 * @param byte_count Bytes in each shard.
 * @param coefficient LCH skew coefficient.
 * @param backend Kernel backend.
 * @param tables Shared Cantor kernel tables returned by `gf2p8::Tables()`.
 * @return Kernel status.
 */
Status IFFTRadix2(Element* x,
                  Element* y,
                  size_t byte_count,
                  Element coefficient,
                  Backend backend,
                  const MultiplicationTables& tables);

/**
 * @brief Applies a radix-two IFFT butterfly into XOR accumulators.
 * @details Reads immutable Cantor-coordinate inputs and fuses output stores
 * with XOR accumulation.
 * @param x First immutable input shard.
 * @param y Second immutable input shard.
 * @param output_x First destination accumulator.
 * @param output_y Second destination accumulator.
 * @param byte_count Bytes in each shard.
 * @param coefficient LCH skew coefficient.
 * @param backend Kernel backend.
 * @param tables Shared Cantor kernel tables returned by `gf2p8::Tables()`.
 * @return Kernel status.
 */
Status IFFTRadix2Xor(const Element* x,
                     const Element* y,
                     Element* output_x,
                     Element* output_y,
                     size_t byte_count,
                     Element coefficient,
                     Backend backend,
                     const MultiplicationTables& tables);

/**
 * @brief Applies two fused in-place FFT layers to four shards.
 * @details Keeps all four Cantor-coordinate shards resident across the paired
 * radix-two layers.
 * @param x0 First butterfly shard.
 * @param x1 Second butterfly shard.
 * @param x2 Third butterfly shard.
 * @param x3 Fourth butterfly shard.
 * @param byte_count Bytes in each shard.
 * @param top Top-layer skew coefficient.
 * @param low Lower-group skew coefficient.
 * @param high Upper-group skew coefficient.
 * @param backend Kernel backend.
 * @param tables Shared Cantor kernel tables returned by `gf2p8::Tables()`.
 * @return Kernel status.
 */
Status FFTRadix4(Element* x0,
                 Element* x1,
                 Element* x2,
                 Element* x3,
                 size_t byte_count,
                 Element top,
                 Element low,
                 Element high,
                 Backend backend,
                 const MultiplicationTables& tables);

/**
 * @brief Applies two fused in-place IFFT layers to four shards.
 * @details Reverses the paired radix-four FFT operation while retaining four
 * Cantor-coordinate shards across both layers.
 * @param x0 First butterfly shard.
 * @param x1 Second butterfly shard.
 * @param x2 Third butterfly shard.
 * @param x3 Fourth butterfly shard.
 * @param byte_count Bytes in each shard.
 * @param top Top-layer skew coefficient.
 * @param low Lower-group skew coefficient.
 * @param high Upper-group skew coefficient.
 * @param backend Kernel backend.
 * @param tables Shared Cantor kernel tables returned by `gf2p8::Tables()`.
 * @return Kernel status.
 */
Status IFFTRadix4(Element* x0,
                  Element* x1,
                  Element* x2,
                  Element* x3,
                  size_t byte_count,
                  Element top,
                  Element low,
                  Element high,
                  Backend backend,
                  const MultiplicationTables& tables);

/**
 * @brief Applies fused radix-four IFFT layers into XOR accumulators.
 * @details Reads four immutable Cantor-coordinate shards and fuses final stores
 * with XOR into four disjoint destination accumulators.
 * @param x0 First immutable input shard.
 * @param x1 Second immutable input shard.
 * @param x2 Third immutable input shard.
 * @param x3 Fourth immutable input shard.
 * @param output0 First destination accumulator.
 * @param output1 Second destination accumulator.
 * @param output2 Third destination accumulator.
 * @param output3 Fourth destination accumulator.
 * @param byte_count Bytes in each shard.
 * @param top Top-layer skew coefficient.
 * @param low Lower-group skew coefficient.
 * @param high Upper-group skew coefficient.
 * @param backend Kernel backend.
 * @param tables Shared Cantor kernel tables returned by `gf2p8::Tables()`.
 * @return Kernel status.
 */
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
                     const MultiplicationTables& tables);

}  // namespace gf2p8::lch
