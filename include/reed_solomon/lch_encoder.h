#pragma once

#include <cstddef>
#include <memory>
#include <span>

#include "lin_chung_han/transform.h"

namespace gf2p8::rs {

/**
 * @brief Encodes systematic Reed-Solomon codes using LCH transforms.
 * @details Arbitrary supported dimensions are embedded into power-of-two
 * mother codes. High-rate codes use folded recovery-subspace encoding;
 * low-rate codes use coefficient fan-out. These paths derive from the Leopard
 * and XDRS algorithms.
 */
class LCHEncoder {
 public:
  /**
   * @brief Constructs an encoder for the requested code dimensions.
   * @details Unsupported dimensions produce an invalid encoder observable
   * through @ref Valid.
   * @param data_count Number of systematic data shards K.
   * @param recovery_count Number of recovery shards R.
   */
  LCHEncoder(size_t data_count, size_t recovery_count);

  /** @brief Destroys the encoder. */
  ~LCHEncoder();

  /**
   * @brief Move-constructs an encoder.
   * @param other Encoder to move from.
   */
  LCHEncoder(LCHEncoder&& other) noexcept;

  /**
   * @brief Move-assigns an encoder.
   * @param other Encoder to move from.
   * @return This encoder.
   */
  LCHEncoder& operator=(LCHEncoder&& other) noexcept;
  LCHEncoder(const LCHEncoder&) = delete;
  LCHEncoder& operator=(const LCHEncoder&) = delete;

  /**
   * @brief Reports whether the requested dimensions are supported.
   * @return True when encoding may be used.
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
   * @return Required bytes, zero for an invalid encoder, or `SIZE_MAX` on
   * overflow.
   */
  size_t WorkspaceSize(size_t byte_count) const;

  /**
   * @brief Encodes data shards into recovery shards.
   * @details Data, recovery, and workspace byte ranges must be pairwise
   * disjoint. Internally selected low- and high-rate evaluation layouts do not
   * affect the separate public spans.
   * @param data Immutable systematic data shard pointers.
   * @param recovery Writable recovery shard pointers.
   * @param byte_count Bytes in each shard.
   * @param workspace Caller-owned temporary byte storage.
   * @param backend Kernel backend.
   * @param radix Transform radix.
   * @return Encode status.
   */
  lch::Status Encode(std::span<const Element* const> data,
                     std::span<Element* const> recovery,
                     size_t byte_count,
                     std::span<Element> workspace,
                     lch::Backend backend = lch::Backend::tuned,
                     lch::Radix radix = lch::Radix::radix4) const;

 private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace gf2p8::rs
