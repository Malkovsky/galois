#pragma once

#include "reed_solomon/strong_weak_rs_product_code.h"

namespace gf2p8::rs::detail {

/** @brief Private differential-test and benchmark access to initial-pass
 * choices. */
struct ProductCorrectionAccess {
  /** @brief Exposes the full R=4 candidate for independent differential tests.
   */
  static CorrectionResult WeakCandidateR4(const StrongWeakRSProductCode& code,
                                          std::span<const Element> row,
                                          std::array<size_t, 2>& positions,
                                          std::array<Element, 2>& magnitudes) {
    return code.WeakCandidateR4(row, positions, magnitudes);
  }
  /**
   * @brief Direct-R2 experiment: bits 0/1/2 select sparse masks, vector
   * reduction, and padded strong batch lanes. Zero retains the R2 baseline.
   * @param code Product dimensions and component decoders.
   * @param block Mutable row-major block.
   * @param options Directional pass limit and independent weak gates.
   * @param optimizations Private per-call bitset, never shared mutable state.
   * @return The normal product outcome and exact work counts.
   */
  static ProductCorrectionResult Experiment(const StrongWeakRSProductCode& code,
                                            std::span<Element> block,
                                            ProductDecodeOptions options,
                                            unsigned optimizations) {
    const unsigned batches = lch::BackendAvailable(lch::Backend::avx2) ? 2 : 0;
    return code.CorrectImpl(block, options, batches, true, true, optimizations);
  }
  /**
   * @brief Private compact weak-row candidate for independent mother tests.
   * @param code Valid product code whose weak dimensions match row.
   * @param row Compact input, unchanged by this operation.
   * @param position Output public position for a one-error candidate only.
   * @param magnitude Output XOR delta for a one-error candidate only.
   * @return Clean, one-error, or uncorrectable outcome.
   */
  static CorrectionResult WeakCandidate(const StrongWeakRSProductCode& code,
                                        std::span<const Element> row,
                                        size_t& position,
                                        Element& magnitude) {
    return code.WeakCandidate(row, position, magnitude);
  }
  /**
   * @brief Runs the same scheduler with zero, one, or two initial batched
   * passes.
   * @param code Product dimensions and component decoders.
   * @param block Mutable row-major block.
   * @param cap Directional pass limit.
   * @param batch_passes Initial passes to batch (0: reference single path).
   * @param tracked_validation Skip known-clean lines; false retains the full
   * scan.
   * @return The normal product outcome and work counts.
   */
  static ProductCorrectionResult Correct(const StrongWeakRSProductCode& code,
                                         std::span<Element> block,
                                         size_t cap,
                                         unsigned batch_passes,
                                         bool tracked_validation = true) {
    return Correct(code, block, ProductDecodeOptions{cap}, batch_passes,
                   tracked_validation);
  }

  /**
   * @brief Runs a private initial-pass/validation variant with per-call gates.
   * @param code Product dimensions and component decoders.
   * @param block Mutable row-major block.
   * @param options Directional pass cap and independent weak gates.
   * @param batch_passes Initial passes to batch (0: reference single path).
   * @param tracked_validation Whether to skip known-clean lines at exit.
   * @return The normal product outcome and work counts.
   */
  static ProductCorrectionResult Correct(const StrongWeakRSProductCode& code,
                                         std::span<Element> block,
                                         ProductDecodeOptions options,
                                         unsigned batch_passes,
                                         bool tracked_validation = true) {
    return code.CorrectImpl(block, options, batch_passes, tracked_validation);
  }
};

}  // namespace gf2p8::rs::detail
