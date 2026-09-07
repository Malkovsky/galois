#pragma once

#include <cstddef>
#include <cstdint>
#include <span>

#include "reed_solomon/lch_decoder.h"

namespace gf2p8::rs::detail::error_correction {

enum class CorrectionStatus {
  ok,
  invalid_argument,
  unsupported_dimensions,
  uncorrectable,
  reconstruction_failed,
};

struct CorrectionResult {
  CorrectionStatus status = CorrectionStatus::invalid_argument;
  size_t error_count = 0;
};

/**
 * @brief Corrects one scalar LCH Reed-Solomon codeword.
 *
 * @param decoder Decoder whose data and recovery dimensions define the code.
 * @param data Data symbols to repair transactionally on successful decoding.
 * @param recovery Immutable recovery symbols.
 * @param error_mask Output mask in public `[data][recovery]` position order.
 * @return The bounded-distance decoding status and detected error count.
 *
 * @note All spans must be pairwise disjoint. A mask that aliases an input is
 * rejected untouched; otherwise the mask is cleared before the remaining
 * contract is validated.
 */
CorrectionResult CorrectOne(const LCHDecoder& decoder,
                            std::span<Element> data,
                            std::span<const Element> recovery,
                            std::span<uint8_t> error_mask);

/**
 * @brief Corrects independent codewords stored across shard byte offsets.
 *
 * @param decoder Decoder whose data and recovery dimensions define the code.
 * @param data Mutable pointers to the K data shards.
 * @param recovery Immutable pointers to the R recovery shards.
 * @param byte_count Number of independent codewords, one per shard byte.
 * @param results Per-codeword bounded-distance outcomes, indexed by byte.
 * @param error_masks Position-major output masks in public `[data][recovery]`
 * order, indexed as `position * byte_count + byte`.
 * @return Call-level contract status. `CorrectionStatus::ok` means every entry
 * in `results` contains that codeword's decoding outcome.
 *
 * @note Every pointed-to shard range and both output spans must be pairwise
 * disjoint.
 */
CorrectionStatus CorrectBatch(const LCHDecoder& decoder,
                              std::span<Element* const> data,
                              std::span<const Element* const> recovery,
                              size_t byte_count,
                              std::span<CorrectionResult> results,
                              std::span<uint8_t> error_masks);

}  // namespace gf2p8::rs::detail::error_correction
