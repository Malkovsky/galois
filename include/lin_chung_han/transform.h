#pragma once

#include <array>
#include <cstddef>
#include <span>

#include "lin_chung_han/kernels.h"

namespace gf2p8::lch {

/**
 * @brief Selects the implementation and radix used by an LCH transform.
 * @details The tuned backend is resolved from the shard byte count before the
 * transform begins. Radix four fuses adjacent radix-two layers.
 */
struct TransformOptions {
  Backend backend = Backend::tuned;
  Radix radix = Radix::radix4;
};

/**
 * @brief Provides the immutable Cantor-coordinate LCH transform domain.
 * @details The process-wide domain uses native Cantor-coordinate evaluation
 * order and compile-time-generated skew and Walsh tables.
 */
class Context {
 public:
  static constexpr size_t kFieldSize = 256;
  static constexpr size_t kFieldBits = 8;

  /**
   * @brief Returns the process-wide immutable LCH context.
   * @details The shared context has native Cantor-coordinate semantics.
   * @return The canonical LCH transform context.
   */
  static const Context& Shared();

  /**
   * @brief Returns the native Cantor-coordinate basis vectors.
   * @details The basis determines the context's additive evaluation ordering.
   * @return Immutable array of basis vectors.
   */
  const std::array<Element, kFieldBits>& Basis() const;

  /**
   * @brief Returns an evaluation point by additive-space index.
   * @details The caller must provide an index smaller than @ref kFieldSize.
   * @param index Additive-space index.
   * @return Evaluation point in native Cantor coordinates.
   */
  Element EvaluationPoint(size_t index) const;

  /**
   * @brief Returns a normalized LCH skew coefficient.
   * @details The caller supplies a valid transform level and field index.
   * @param level Transform level.
   * @param index Evaluation-space index.
   * @return Skew coefficient for the requested butterfly group.
   */
  Element Skew(size_t level, size_t index) const;

  /**
   * @brief Returns the field logarithm of a value.
   * @details The zero entry follows the context's locator-table convention.
   * @param value Field element.
   * @return Logarithm-table entry.
   */
  Element Log(Element value) const;

  /**
   * @brief Returns a field exponent.
   * @details The exponent is reduced modulo 255.
   * @param exponent Integer exponent.
   * @return Corresponding nonzero field element.
   */
  Element Exp(size_t exponent) const;

  /**
   * @brief Returns a Walsh-transformed logarithm entry.
   * @details These entries support erasure-locator construction.
   * @param index Field index.
   * @return Walsh-domain logarithm value.
   */
  Element LogWalsh(size_t index) const;

  /**
   * @brief Returns immutable backend multiplication tables.
   * @details The tables include packed shuffle and GFNI affine forms.
   * @return The process-wide field multiplication tables.
   */
  const MultiplicationTables& Tables() const { return gf2p8::Tables(); }

 private:
  Context() = default;
};

/**
 * @brief Computes an in-place LCH FFT.
 * @details The shard count must be an aligned power of two from 1 through 256.
 * Arbitrary byte lengths are supported with exact scalar tails.
 * @param context Immutable transform context.
 * @param shards In-place evaluation shards.
 * @param byte_count Bytes in each shard.
 * @param evaluation_offset Aligned evaluation-space offset.
 * @param options Backend and radix selection.
 * @return Transform status.
 */
Status FFT(const Context& context,
           std::span<Element* const> shards,
           size_t byte_count,
           size_t evaluation_offset = 0,
           TransformOptions options = {});

/**
 * @brief Computes a prefix-truncated in-place LCH FFT.
 * @details Only the first @p output_count shards are guaranteed on return;
 * unrequested shard contents are unspecified.
 * @param context Immutable transform context.
 * @param shards In-place evaluation shards.
 * @param byte_count Bytes in each shard.
 * @param evaluation_offset Aligned evaluation-space offset.
 * @param output_count Number of guaranteed prefix outputs.
 * @param options Backend and radix selection.
 * @return Transform status.
 */
Status FFT(const Context& context,
           std::span<Element* const> shards,
           size_t byte_count,
           size_t evaluation_offset,
           size_t output_count,
           TransformOptions options = {});

/**
 * @brief Computes a sparse-output in-place LCH FFT.
 * @details Only outputs with a nonzero request-mask entry are guaranteed;
 * unrequested shard contents are unspecified.
 * @param context Immutable transform context.
 * @param shards In-place evaluation shards.
 * @param byte_count Bytes in each shard.
 * @param evaluation_offset Aligned evaluation-space offset.
 * @param requested_outputs Per-shard output request mask.
 * @param options Backend and radix selection.
 * @return Transform status.
 */
Status FFT(const Context& context,
           std::span<Element* const> shards,
           size_t byte_count,
           size_t evaluation_offset,
           std::span<const uint8_t> requested_outputs,
           TransformOptions options = {});

/**
 * @brief Computes an in-place LCH IFFT.
 * @details The shard count must be an aligned power of two from 1 through 256.
 * Arbitrary byte lengths are supported with exact scalar tails.
 * @param context Immutable transform context.
 * @param shards In-place evaluation shards.
 * @param byte_count Bytes in each shard.
 * @param evaluation_offset Aligned evaluation-space offset.
 * @param options Backend and radix selection.
 * @return Transform status.
 */
Status IFFT(const Context& context,
            std::span<Element* const> shards,
            size_t byte_count,
            size_t evaluation_offset = 0,
            TransformOptions options = {});

/**
 * @brief Computes an in-place LCH IFFT from a live input prefix.
 * @details The suffix `[input_count, shards.size())` must already contain zero
 * bytes. The transform may exploit the shorter live prefix.
 * @param context Immutable transform context.
 * @param shards In-place evaluation shards with a pre-zeroed suffix.
 * @param byte_count Bytes in each shard.
 * @param evaluation_offset Aligned evaluation-space offset.
 * @param input_count Number of live prefix shards.
 * @param options Backend and radix selection.
 * @return Transform status.
 */
Status IFFT(const Context& context,
            std::span<Element* const> shards,
            size_t byte_count,
            size_t evaluation_offset,
            size_t input_count,
            TransformOptions options = {});

/**
 * @brief Computes an LCH IFFT from immutable input into work shards.
 * @details @p input is a live prefix; remaining @p work shards are initialized
 * to zero. Input and work byte ranges must be disjoint.
 * @param context Immutable transform context.
 * @param input Immutable live-prefix source shards.
 * @param work Output and temporary transform shards.
 * @param byte_count Bytes in each shard.
 * @param evaluation_offset Aligned evaluation-space offset.
 * @param options Backend and radix selection.
 * @return Transform status.
 */
Status IFFT(const Context& context,
            std::span<const Element* const> input,
            std::span<Element* const> work,
            size_t byte_count,
            size_t evaluation_offset = 0,
            TransformOptions options = {});

/**
 * @brief Computes an immutable-input IFFT and XORs it into an accumulator.
 * @details Fuses the final transform layer into
 * `xor_accumulator ^= IFFT(input)`. Input, work, and accumulator byte ranges
 * must be pairwise disjoint; final work contents are unspecified.
 * @param context Immutable transform context.
 * @param input Immutable live-prefix source shards.
 * @param work Temporary transform shards.
 * @param xor_accumulator Destination accumulator shards.
 * @param byte_count Bytes in each shard.
 * @param evaluation_offset Aligned evaluation-space offset.
 * @param options Backend and radix selection.
 * @return Transform status.
 */
Status IFFT(const Context& context,
            std::span<const Element* const> input,
            std::span<Element* const> work,
            std::span<Element* const> xor_accumulator,
            size_t byte_count,
            size_t evaluation_offset = 0,
            TransformOptions options = {});

}  // namespace gf2p8::lch
