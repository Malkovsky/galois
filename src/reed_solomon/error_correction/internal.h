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

// Corrects one scalar codeword. All spans must be pairwise disjoint. A mask
// that aliases input is rejected untouched; otherwise the mask is cleared
// before validating the remaining contract.
CorrectionResult CorrectOne(const LCHDecoder& decoder,
                            std::span<Element> data,
                            std::span<const Element> recovery,
                            std::span<uint8_t> error_mask);

// Corrects byte_count independent codewords stored across shard byte offsets.
// Results are indexed by byte offset. The mask is position-major in public
// [data][recovery] order: error_masks[position * byte_count + byte]. A returned
// ok status means the call contract was accepted; each result reports that
// codeword's bounded-distance outcome. All pointed-to shard ranges and output
// spans must be pairwise disjoint.
CorrectionStatus CorrectBatch(const LCHDecoder& decoder,
                              std::span<Element* const> data,
                              std::span<const Element* const> recovery,
                              size_t byte_count,
                              std::span<CorrectionResult> results,
                              std::span<uint8_t> error_masks);

}  // namespace gf2p8::rs::detail::error_correction
