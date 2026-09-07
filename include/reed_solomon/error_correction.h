#pragma once

#include "reed_solomon/lch_decoder.h"

namespace gf2p8::rs {

/** @brief Bounded-distance correction outcome. */
enum class CorrectionStatus {
  ok,
  invalid_argument,
  unsupported_dimensions,
  uncorrectable,
  reconstruction_failed,
};

/** @brief Status and number of repaired symbols (zero on failure). */
struct CorrectionResult {
  CorrectionStatus status = CorrectionStatus::invalid_argument;
  size_t error_count = 0;
};

/**
 * @brief Transactionally repairs data AND parity of one Cantor RS codeword.
 * @param decoder Decoder defining K data and R recovery symbols.
 * @param codeword Exactly K+R mutable symbols in [data][recovery] order.
 * @return Status and corrected symbol count; zero count on an already valid word.
 * @details Requires unshortened power-of-two N=K+R <=256 and power-of-two
 * R<=K. Corrects up to floor(R/2) unknown symbol errors, with no erasures.
 * On any non-ok status every input byte is unchanged. Success means the result
 * is a codeword within that radius, NOT necessarily the transmitted codeword.
 */
CorrectionResult CorrectCodeword(const LCHDecoder& decoder,
                                std::span<Element> codeword);

}  // namespace gf2p8::rs
