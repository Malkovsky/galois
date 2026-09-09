#include "reed_solomon/strong_weak_rs_product_code.h"

#include <algorithm>
#include <array>
#include <bit>
#include <vector>

#include "reed_solomon/error_correction/internal.h"

namespace gf2p8::rs {
namespace {

bool Aligned(size_t n, size_t k) {
  return n <= 256 && std::has_single_bit(n) && k < n && n - k >= 2 &&
         n - k <= k && std::has_single_bit(n - k);
}

}  // namespace

StrongWeakRSProductCode::StrongWeakRSProductCode(size_t strong_n,
                                                 size_t strong_k,
                                                 size_t weak_n,
                                                 size_t weak_k)
    : strong_n_(strong_n),
      strong_k_(strong_k),
      weak_n_(weak_n),
      weak_k_(weak_k),
      valid_(Aligned(strong_n, strong_k) && Aligned(weak_n, weak_k) &&
             weak_n - weak_k == 2),
      strong_encoder_(valid_ ? strong_k : 0, valid_ ? strong_n - strong_k : 0),
      weak_encoder_(valid_ ? weak_k : 0, valid_ ? weak_n - weak_k : 0),
      strong_decoder_(valid_ ? strong_k : 0, valid_ ? strong_n - strong_k : 0),
      weak_decoder_(valid_ ? weak_k : 0, valid_ ? weak_n - weak_k : 0) {}

bool StrongWeakRSProductCode::Valid() const {
  return valid_ && strong_encoder_.Valid() && weak_encoder_.Valid() &&
         strong_decoder_.Valid() && weak_decoder_.Valid();
}

size_t StrongWeakRSProductCode::BlockSize() const {
  return Valid() ? strong_n_ * weak_n_ : 0;
}

lch::Status StrongWeakRSProductCode::Encode(std::span<Element> block,
                                            lch::Backend backend) const {
  if (!Valid() || block.size() != BlockSize()) {
    return lch::Status::invalid_argument;
  }
  std::vector<Element> candidate(block.begin(), block.end());
  std::array<const Element*, 256> data{};
  std::array<Element*, 256> recovery{};
  std::vector<Element> workspace(std::max(
      strong_encoder_.WorkspaceSize(weak_n_), weak_encoder_.WorkspaceSize(1)));
  for (size_t row = 0; row < strong_k_; ++row) {
    for (size_t col = 0; col < weak_k_; ++col) {
      data[col] = &candidate[row * weak_n_ + col];
    }
    for (size_t col = weak_k_; col < weak_n_; ++col) {
      recovery[col - weak_k_] = &candidate[row * weak_n_ + col];
    }
    const auto status = weak_encoder_.Encode(std::span(data).first(weak_k_),
                                             std::span(recovery).first(2), 1,
                                             workspace, backend);
    if (status != lch::Status::ok) {
      return status;
    }
  }
  for (size_t row = 0; row < strong_k_; ++row) {
    data[row] = &candidate[row * weak_n_];
  }
  for (size_t row = strong_k_; row < strong_n_; ++row) {
    recovery[row - strong_k_] = &candidate[row * weak_n_];
  }
  const auto status =
      strong_encoder_.Encode(std::span(data).first(strong_k_),
                             std::span(recovery).first(strong_n_ - strong_k_),
                             weak_n_, workspace, backend);
  if (status == lch::Status::ok) {
    std::copy(candidate.begin(), candidate.end(), block.begin());
  }
  return status;
}

ProductCorrectionResult StrongWeakRSProductCode::Correct(
    std::span<Element> block,
    size_t max_directional_passes) const {
  return Correct(block, ProductDecodeOptions{max_directional_passes});
}

ProductCorrectionResult StrongWeakRSProductCode::Correct(
    std::span<Element> block,
    ProductDecodeOptions options) const {
  // Avoid packing overhead when no complete SIMD batch can run.
  const unsigned batches = lch::BackendAvailable(lch::Backend::avx2) &&
                                   std::max(strong_n_, weak_n_) >= 32
                               ? 2
                               : 0;
  return CorrectImpl(block, options, batches);
}

ProductCorrectionResult StrongWeakRSProductCode::CorrectImpl(
    std::span<Element> block,
    ProductDecodeOptions options,
    unsigned batch_passes,
    bool tracked_validation) const {
  const size_t max_directional_passes = options.max_directional_passes;
  ProductCorrectionResult result;
  if (!Valid() || block.size() != BlockSize() || max_directional_passes < 2) {
    return result;
  }
  std::array<bool, 256> protected_columns{};
  // False means unknown, never proven nonzero: intersecting edits can cancel.
  std::array<bool, 256> clean_columns{}, clean_rows{};
  std::array<bool, 256> active{};
  std::array<Element, 256> candidate{};
  std::vector<Element> packed(batch_passes != 0 ? block.size() : 0);
  std::vector<uint8_t> masks(packed.size());
  std::array<CorrectionResult, 256> outcomes{};
  std::array<Element*, 256> shards{};
  for (size_t pass = 0; pass < max_directional_passes; ++pass) {
    const bool strong = pass % 2 == 0;
    const size_t lines = strong ? weak_n_ : strong_n_;
    const size_t length = strong ? strong_n_ : weak_n_;
    std::array<bool, 256> next{};
    size_t changes = 0;
    size_t bit_changes = 0;
    const bool batched = pass < std::min(batch_passes, 2u);
    if (batched) {
      // Columns already have position-major layout; rows need transposition.
      // Keep tentative repairs private until the existing weak gates accept.
      if (strong) {
        std::copy(block.begin(), block.end(), packed.begin());
      }
      for (size_t pos = 0; pos < length; ++pos) {
        shards[pos] = packed.data() + pos * lines;
        if (!strong) {
          for (size_t line = 0; line < lines; ++line) {
            shards[pos][line] = block[line * weak_n_ + pos];
          }
        }
      }
      const auto status = detail::error_correction::CorrectCodewordBatch(
          strong ? strong_decoder_ : weak_decoder_,
          std::span(shards).first(length), lines,
          std::span(outcomes).first(lines), masks);
      if (status != CorrectionStatus::ok) {
        return result;
      }
    }
    if (batched && strong) {
      result.strong_lines_visited += lines;
      for (size_t line = 0; line < lines; ++line) {
        clean_columns[line] = protected_columns[line] =
            outcomes[line].status == CorrectionStatus::ok;
      }
      // Batch failure is transactional and masks describe only verified edits.
      // Walk position-major output rather than gathering every column again.
      for (size_t pos = 0; pos < length; ++pos) {
        for (size_t line = 0; line < lines; ++line) {
          const size_t index = pos * lines + line;
          if (masks[index]) {
            bit_changes += std::popcount(
                static_cast<unsigned>(block[index] ^ packed[index]));
            block[index] = packed[index];
            next[pos] = true;
            clean_rows[pos] = false;
            ++changes;
          }
        }
      }
    } else {
      for (size_t line = 0; line < lines; ++line) {
        if (pass >= 2 && !active[line]) {
          continue;
        }
        if (strong) {
          ++result.strong_lines_visited;
        } else {
          ++result.weak_lines_visited;
        }
        const auto index = [&](size_t pos) {
          return strong ? pos * weak_n_ + line : line * weak_n_ + pos;
        };
        if (!batched) {
          for (size_t pos = 0; pos < length; ++pos) {
            candidate[pos] = block[index(pos)];
          }
        }
        const auto correction =
            batched ? outcomes[line]
                    : CorrectCodeword(strong ? strong_decoder_ : weak_decoder_,
                                      std::span(candidate).first(length));
        auto& clean = strong ? clean_columns[line] : clean_rows[line];
        clean = correction.status == CorrectionStatus::ok &&
                correction.error_count == 0;
        if (strong) {
          protected_columns[line] = correction.status == CorrectionStatus::ok;
        }
        if (correction.status != CorrectionStatus::ok) {
          continue;
        }
        if (correction.error_count == 0) {
          continue;
        }
        if (batched) {
          for (size_t pos = 0; pos < length; ++pos) {
            candidate[pos] = shards[pos][line];
          }
        }
        if (!strong) {
          if (correction.error_count != 1) {
            continue;
          }
          bool accept = true;
          for (size_t pos = 0; pos < length; ++pos) {
            const unsigned delta = block[index(pos)] ^ candidate[pos];
            if (delta != 0 &&
                ((options.use_anchors && protected_columns[pos]) ||
                 (options.use_binary_image && std::popcount(delta) > 2))) {
              accept = false;
            }
          }
          if (!accept) {
            continue;
          }
        }
        for (size_t pos = 0; pos < length; ++pos) {
          if (block[index(pos)] != candidate[pos]) {
            bit_changes += std::popcount(
                static_cast<unsigned>(block[index(pos)] ^ candidate[pos]));
            block[index(pos)] = candidate[pos];
            next[pos] = true;
            (strong ? clean_rows[pos] : clean_columns[pos]) = false;
            ++changes;
          }
        }
        // Successful BDD verifies its candidate internally; only committed
        // candidates establish validity. Rejected weak candidates do not.
        clean = true;
      }
    }
    ++result.directional_passes;
    result.changed_symbols += changes;
    result.changed_bits += bit_changes;
    (strong ? result.strong_changed_symbols : result.weak_changed_symbols) +=
        changes;
    (strong ? result.strong_changed_bits : result.weak_changed_bits) +=
        bit_changes;
    active = next;
    if (pass >= 1 && changes == 0) {
      result.termination = ProductTermination::no_change;
      break;
    }
    result.termination = ProductTermination::pass_limit;
  }
  // Check the actual final block, not stale protection or tentative candidates.
  result.all_zero_syndromes = true;
  for (bool strong : {true, false}) {
    const size_t lines = strong ? weak_n_ : strong_n_;
    const size_t length = strong ? strong_n_ : weak_n_;
    for (size_t line = 0; line < lines; ++line) {
      if (tracked_validation &&
          (strong ? clean_columns[line] : clean_rows[line])) {
        continue;
      }
      for (size_t pos = 0; pos < length; ++pos) {
        candidate[pos] =
            block[strong ? pos * weak_n_ + line : line * weak_n_ + pos];
      }
      const auto check =
          CorrectCodeword(strong ? strong_decoder_ : weak_decoder_,
                          std::span(candidate).first(length));
      if (check.status != CorrectionStatus::ok || check.error_count != 0) {
        result.all_zero_syndromes = false;
      }
    }
  }
  return result;
}

}  // namespace gf2p8::rs
