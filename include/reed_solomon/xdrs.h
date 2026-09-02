#pragma once

#include <cstddef>
#include <memory>
#include <span>

#include "reed_solomon/lin_chung_han/transform.h"

namespace gf2p8::rs {

/**
 * @brief Selects the published low-rate or high-rate XDRS algorithm.
 * @details Low rate uses a K-point coefficient fan-out; high rate uses folded
 * recovery-sized LCH transforms.
 */
enum class XDRSRate { low, high };

/**
 * @brief Implements full-length systematic XDRS coding for N=256.
 * @details The codec uses the shared native Cantor-coordinate domain
 * while preserving the published rate-specific layouts and orchestration. It
 * is not byte-compatible with the released polynomial-coordinate artifact.
 */
class XDRS {
 public:
  /**
   * @brief Constructs an XDRS codec.
   * @details Low rate requires power-of-two K no greater than 128. High rate
   * requires power-of-two `256-K` no greater than 128.
   * @param data_count Number of systematic data shards K.
   * @param rate Published XDRS rate algorithm.
   */
  XDRS(size_t data_count, XDRSRate rate);

  /**
   * @brief Destroys the codec.
   * @details Releases the owned implementation and rate-specific coefficients.
   */
  ~XDRS();

  /**
   * @brief Move-constructs a codec.
   * @details Transfers ownership of the implementation. The moved-from object
   * may only be destroyed or assigned another value.
   * @param other Codec to move from.
   */
  XDRS(XDRS&& other) noexcept;

  /**
   * @brief Move-assigns a codec.
   * @details Transfers ownership of the implementation. The moved-from object
   * may only be destroyed or assigned another value.
   * @param other Codec to move from.
   * @return This codec.
   */
  XDRS& operator=(XDRS&& other) noexcept;
  XDRS(const XDRS&) = delete;
  XDRS& operator=(const XDRS&) = delete;

  /**
   * @brief Reports whether construction parameters are supported.
   * @details Invalid dimensions retain an inspectable but unusable codec.
   * @return True when encode and decode operations may be used.
   */
  bool Valid() const;

  /**
   * @brief Returns the systematic data shard count.
   * @details XDRS always has total codeword length 256.
   * @return Data shard count K.
   */
  size_t DataCount() const;

  /**
   * @brief Returns the recovery shard count.
   * @details The value is `256-K` for a valid codec.
   * @return Recovery shard count.
   */
  size_t RecoveryCount() const;

  /**
   * @brief Returns required encode workspace bytes.
   * @details Low-rate encode requires no workspace; high-rate encode uses one
   * recovery-sized scratch block.
   * @param byte_count Bytes in each shard.
   * @return Required workspace bytes, zero for an invalid codec, or `SIZE_MAX`
   * on overflow.
   */
  size_t EncodeWorkspaceSize(size_t byte_count) const;

  /**
   * @brief Returns required decode workspace bytes.
   * @details Workspace contains 256 locator bytes followed by 256 transform
   * scratch shards for either rate.
   * @param byte_count Bytes in each shard.
   * @return Required workspace bytes, zero for an invalid codec, or `SIZE_MAX`
   * on overflow.
   */
  size_t DecodeWorkspaceSize(size_t byte_count) const;

  /**
   * @brief Returns the selected XDRS rate algorithm.
   * @details The rate controls layout and transform orchestration.
   * @return Selected rate.
   */
  XDRSRate Rate() const;

  /**
   * @brief Encodes systematic data shards into recovery shards.
   * @details Low-rate layout is `[data][recovery]`; high-rate layout is
   * `[recovery][data]`. Encode writes only @p recovery. Data, recovery, and
   * workspace byte ranges must be pairwise disjoint.
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

  /**
   * @brief Recovers missing data shards from an XDRS codeword.
   * @details The codeword uses the selected rate's layout. Missing data entries
   * require writable pointers and a zero entry in @p present. Codeword and
   * workspace byte ranges must be pairwise disjoint.
   * @param codeword Mutable data and recovery shard pointers.
   * @param present Presence mask for all 256 codeword shards.
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
