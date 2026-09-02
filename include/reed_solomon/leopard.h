#pragma once

#include <cstddef>
#include <memory>
#include <span>

#include "reed_solomon/lin_chung_han/transform.h"

namespace gf2p8::rs {

/**
 * @brief Implements shortened systematic Leopard Reed-Solomon coding.
 * @details The logical codeword layout is `[recovery][original]`. The codec
 * uses the shared native Cantor-coordinate LCH domain and supports
 * safe arbitrary shard byte lengths.
 */
class Leopard {
 public:
  /**
   * @brief Constructs a Leopard codec.
   * @details Unsupported dimensions produce an invalid codec observable
   * through @ref Valid.
   * @param original_count Number of systematic original shards.
   * @param recovery_count Number of recovery shards.
   */
  Leopard(size_t original_count, size_t recovery_count);

  /**
   * @brief Destroys the codec.
   * @details Releases the owned implementation.
   */
  ~Leopard();

  /**
   * @brief Move-constructs a codec.
   * @details Transfers ownership of the implementation. The moved-from object
   * may only be destroyed or assigned another value.
   * @param other Codec to move from.
   */
  Leopard(Leopard&& other) noexcept;

  /**
   * @brief Move-assigns a codec.
   * @details Transfers ownership of the implementation. The moved-from object
   * may only be destroyed or assigned another value.
   * @param other Codec to move from.
   * @return This codec.
   */
  Leopard& operator=(Leopard&& other) noexcept;
  Leopard(const Leopard&) = delete;
  Leopard& operator=(const Leopard&) = delete;

  /**
   * @brief Reports whether construction parameters are supported.
   * @details Valid Leopard dimensions satisfy `R > 0`, `R <= K`, and the FF8
   * transform-size limit.
   * @return True when encode and decode operations may be used.
   */
  bool Valid() const;

  /**
   * @brief Returns the number of original shards.
   * @details This is the systematic data dimension K.
   * @return Original shard count.
   */
  size_t OriginalCount() const;

  /**
   * @brief Returns the number of recovery shards.
   * @details This is the parity dimension R.
   * @return Recovery shard count.
   */
  size_t RecoveryCount() const;

  /**
   * @brief Returns required encode workspace bytes.
   * @details The workspace is caller-owned and depends on the rounded recovery
   * transform size and shard byte count.
   * @param byte_count Bytes in each shard.
   * @return Required workspace bytes, zero for an invalid codec, or `SIZE_MAX`
   * on overflow.
   */
  size_t EncodeWorkspaceSize(size_t byte_count) const;

  /**
   * @brief Returns required decode workspace bytes.
   * @details Workspace includes erasure-locator storage and transform shards.
   * @param byte_count Bytes in each shard.
   * @return Required workspace bytes, zero for an invalid codec, or `SIZE_MAX`
   * on overflow.
   */
  size_t DecodeWorkspaceSize(size_t byte_count) const;

  /**
   * @brief Encodes systematic original shards into recovery shards.
   * @details Workspace must contain at least @ref EncodeWorkspaceSize bytes.
   * Original, recovery, and workspace byte ranges must be pairwise disjoint.
   * @param original Immutable original shard pointers.
   * @param recovery Writable recovery shard pointers.
   * @param byte_count Bytes in each shard.
   * @param workspace Caller-owned temporary byte storage.
   * @param backend Kernel backend.
   * @param radix Transform radix.
   * @return Encode status.
   */
  lch::Status Encode(std::span<const Element* const> original,
                     std::span<Element* const> recovery,
                     size_t byte_count,
                     std::span<Element> workspace,
                     lch::Backend backend = lch::Backend::tuned,
                     lch::Radix radix = lch::Radix::radix4) const;

  /**
   * @brief Recovers missing original shards from a Leopard codeword.
   * @details The codeword layout is `[recovery][original]`. Missing original
   * entries require writable pointers and a zero entry in @p present. Codeword
   * and workspace byte ranges must be pairwise disjoint.
   * @param codeword Mutable recovery and original shard pointers.
   * @param present Presence mask for every codeword shard.
   * @param byte_count Bytes in each shard.
   * @param workspace Caller-owned temporary byte storage.
   * @param backend Kernel backend.
   * @param radix Transform radix.
   * @return Decode status.
   */
  lch::Status Decode(std::span<Element* const> codeword,
                     std::span<const uint8_t> present,
                     size_t byte_count,
                     std::span<Element> workspace,
                     lch::Backend backend = lch::Backend::tuned,
                     lch::Radix radix = lch::Radix::radix4) const;

 private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace gf2p8::rs
