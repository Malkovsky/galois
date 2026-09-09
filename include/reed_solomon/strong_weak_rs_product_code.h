#pragma once

#include "reed_solomon/error_correction.h"
#include "reed_solomon/lch_encoder.h"

namespace gf2p8::rs {

namespace detail {
struct ProductCorrectionAccess;
}

/** @brief Reason iterative product correction stopped, independent of validity.
 */
enum class ProductTermination { invalid_argument, no_change, pass_limit };

/** @brief Per-call product correction limits and independent weak acceptance
 * gates. */
struct ProductDecodeOptions {
  /** @brief Cap counting each direction separately; must be at least two. */
  size_t max_directional_passes = 16;
  /** @brief Reject weak repairs into columns protected by strong BDD success.
   */
  bool use_anchors = true;
  /** @brief Require weak repair deltas to have at most two set bits. */
  bool use_binary_image = true;
};

/** @brief Product correction outcome; validity does not prove original content.
 */
struct ProductCorrectionResult {
  ProductTermination termination = ProductTermination::invalid_argument;
  bool all_zero_syndromes = false;
  size_t directional_passes = 0;
  size_t strong_lines_visited = 0;
  size_t weak_lines_visited = 0;
  /** @brief Accepted symbol writes, counting repeated repairs separately. */
  size_t changed_symbols = 0;
  /** @brief Accepted bit toggles, including repeated committed repairs. */
  size_t changed_bits = 0;
  /** @brief Accepted strong-direction symbol writes and bit toggles. */
  size_t strong_changed_symbols = 0, strong_changed_bits = 0;
  /** @brief Accepted weak-direction symbol writes and bit toggles. */
  size_t weak_changed_symbols = 0, weak_changed_bits = 0;
};

/**
 * @brief Systematic Cantor RS product with strong columns and weak rows.
 * @details Row-major block has Nstrong rows and Nweak columns. The top-left
 * Kstrong by Kweak rectangle holds data; every column and every row, including
 * parity regions, is a component codeword. Errors only: no erasures or
 * backtracking. Strong N and R must be powers of two, N<=256 and 2<=R<=K.
 * Weak R=2, K>=2, N<=256 may be shortened from nextPow2(N): omitted data
 * [K,nextPow2(N)-2) are known zeros. Public rows remain compact [data][parity].
 * Mother-code candidates changing any omitted zero are rejected in full.
 */
class StrongWeakRSProductCode {
 public:
  /**
   * @brief Constructs a product code; unsupported dimensions make Valid false.
   * @param strong_n Number of rows.
   * @param strong_k Number of systematic rows.
   * @param weak_n Number of columns.
   * @param weak_k Number of systematic columns.
   */
  StrongWeakRSProductCode(size_t strong_n = 256,
                          size_t strong_k = 224,
                          size_t weak_n = 256,
                          size_t weak_k = 254);

  /** @brief Reports supported dimensions. @return Whether operations are valid.
   */
  bool Valid() const;
  /** @brief Returns required row-major block size. @return Bytes, or zero if
   * invalid. */
  size_t BlockSize() const;

  /**
   * @brief Fills all product parity, preserving the systematic rectangle.
   * @param block Exactly BlockSize() symbols, with data already in place.
   * @param backend Owned encoder backend, scalar available for reference
   * checks.
   * @return Encoding status; on failure block is unchanged.
   */
  lch::Status Encode(std::span<Element> block,
                     lch::Backend backend = lch::Backend::tuned) const;

  /**
   * @brief Applies alternating strong BDD and gated weak single-error passes.
   * @param block Exactly BlockSize() mutable symbols.
   * @param max_directional_passes Cap counting each direction separately; >=2.
   * @return Termination, pass/work counts, and final all-component validity.
   * @details Always starts with a full strong pass then a full weak pass, even
   * if strong changes nothing. From then on visits only lines intersecting
   * accepted byte changes in the preceding pass. Stops on a no-change pass
   * (including the initial weak pass), or the cap; no-change takes precedence.
   * Strong success, including zero syndromes, protects that column; failure
   * leaves its bytes unchanged and unprotects it. Unvisited protection
   * persists. Weak candidates commit only for one changed symbol with
   * popcount(old XOR new)<=2 in an unprotected column; rejected rows are
   * unchanged. Protection changes alone never activate lines. Accepted repairs
   * are retained at exit; the entire iterative operation is not transactional.
   * Invalid arguments leave the block untouched. Validity is checked separately
   * at exit and does not count as a directional pass or guarantee the original
   * message.
   */
  ProductCorrectionResult Correct(std::span<Element> block,
                                  size_t max_directional_passes = 16) const;

  /**
   * @brief Corrects with independent per-call weak acceptance gates.
   * @param block Exactly BlockSize() mutable symbols.
   * @param options Pass cap and optional anchor and binary-image gates.
   * @return Termination, pass/work counts, and final all-component validity.
   * @details Uses the same scheduling and stopping rules as the cap overload.
   * Weak BDD always requires exactly one changed symbol. Disabling anchors
   * bypasses only target-column protection; disabling binary image bypasses
   * only the two-bit delta limit. Every accepted write still invalidates
   * intersecting cached validity and activates the next direction.
   */
  ProductCorrectionResult Correct(std::span<Element> block,
                                  ProductDecodeOptions options) const;

 private:
  friend struct detail::ProductCorrectionAccess;
  ProductCorrectionResult CorrectImpl(std::span<Element> block,
                                      ProductDecodeOptions options,
                                      unsigned batch_passes,
                                      bool tracked_validation = true) const;
  size_t strong_n_, strong_k_, weak_n_, weak_k_;
  bool valid_;
  LCHEncoder strong_encoder_, weak_encoder_;
  LCHDecoder strong_decoder_, weak_decoder_;
};

}  // namespace gf2p8::rs
