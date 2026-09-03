#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>

#include "lin_chung_han/transform.h"

namespace gf2p8::rs {

/**
 * @brief Recovers systematic Reed-Solomon data using LCH transforms.
 * @details Arbitrary supported dimensions are embedded into power-of-two
 * mother codes, then dispatched to low- or high-rate decoders derived from
 * XDRS. Rate-specific evaluation layouts remain private.
 */
class LCHDecoder {
 public:
  /**
   * @brief Constructs a decoder for the requested code dimensions.
   * @details Construction precomputes only coefficients required by the
   * selected decode strategy. Unsupported dimensions produce an invalid
   * decoder observable through @ref Valid.
   * @param data_count Number of systematic data shards K.
   * @param recovery_count Number of recovery shards R.
   */
  LCHDecoder(size_t data_count, size_t recovery_count);

  /** @brief Destroys the decoder and its precomputed coefficients. */
  ~LCHDecoder();

  /**
   * @brief Move-constructs a decoder.
   * @param other Decoder to move from.
   */
  LCHDecoder(LCHDecoder&& other) noexcept;

  /**
   * @brief Move-assigns a decoder.
   * @param other Decoder to move from.
   * @return This decoder.
   */
  LCHDecoder& operator=(LCHDecoder&& other) noexcept;
  LCHDecoder(const LCHDecoder&) = delete;
  LCHDecoder& operator=(const LCHDecoder&) = delete;

  /**
   * @brief Reports whether the requested dimensions are supported.
   * @return True when decoding may be used.
   */
  bool Valid() const;

  /**
   * @brief Returns the systematic data shard count K.
   * @return Data shard count.
   */
  size_t DataCount() const;

  /**
   * @brief Returns the recovery shard count R.
   * @return Recovery shard count.
   */
  size_t RecoveryCount() const;

  /**
   * @brief Returns required caller-owned workspace bytes.
   * @param byte_count Bytes in each shard.
   * @return Required bytes, zero for an invalid decoder, or `SIZE_MAX` on
   * overflow.
   */
  size_t WorkspaceSize(size_t byte_count) const;

  /**
   * @brief Recovers missing data shards.
   * @details When @p byte_count is nonzero, every data entry requires a
   * non-null pointer: present data is read and missing data is overwritten.
   * Missing data entries require a zero in @p data_present. Present recovery
   * entries require readable pointers; missing recovery entries may be null.
   * All shard and workspace byte ranges must be pairwise disjoint.
   * @param data Mutable systematic data shard pointers.
   * @param data_present Presence mask for data shards.
   * @param recovery Immutable recovery shard pointers.
   * @param recovery_present Presence mask for recovery shards.
   * @param byte_count Bytes in each shard.
   * @param workspace Caller-owned temporary byte storage.
   * @param backend Kernel backend.
   * @param radix Transform radix.
   * @return Decode status.
   */
  lch::Status Decode(std::span<Element* const> data,
                     std::span<const uint8_t> data_present,
                     std::span<const Element* const> recovery,
                     std::span<const uint8_t> recovery_present,
                     size_t byte_count,
                     std::span<Element> workspace,
                     lch::Backend backend = lch::Backend::tuned,
                     lch::Radix radix = lch::Radix::radix4) const;

 private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace gf2p8::rs
