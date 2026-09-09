#pragma once

#include "reed_solomon/strong_weak_rs_product_code.h"

namespace gf2p8::rs::detail {

/** @brief Private differential-test and benchmark access to initial-pass
 * choices. */
struct ProductCorrectionAccess {
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
