#include "reed_solomon/strong_weak_rs_product_code.h"

#include <algorithm>
#include <array>
#include <bit>
#include <vector>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

#include "reed_solomon/error_correction/internal.h"

namespace gf2p8::rs {
namespace {

bool Aligned(size_t n, size_t k) {
  return n <= 256 && std::has_single_bit(n) && k < n && n - k >= 2 &&
         n - k <= k && std::has_single_bit(n - k);
}

std::array<Element, 2> WeakDataSyndromes(std::span<const Element> data,
                                         bool vector_reduce = true) {
  Element s0 = 0, s1 = 0;
  const auto& tables = Tables().shuffle;
  size_t j = 0;
#if defined(__AVX2__)
  if (vector_reduce && data.size() >= 32 &&
      lch::BackendAvailable(lch::Backend::avx2)) {
    // Distribute the dot product over the eight coordinate bits of each
    // position's weight. Each group needs only one fixed-factor product.
    static constexpr auto masks = [] {
      std::array<std::array<std::array<uint8_t, 32>, 8>, 8> result{};
      for (size_t chunk = 0; chunk < 8; ++chunk) {
        for (size_t bit = 0; bit < 8; ++bit) {
          for (size_t lane = 0; lane < 32; ++lane) {
            result[chunk][bit][lane] =
                (((chunk * 32 + lane + 2) ^ 1) & (1u << bit)) ? 255 : 0;
          }
        }
      }
      return result;
    }();
    __m256i groups[8]{};
    auto sum = _mm256_setzero_si256();
    for (; j + 32 <= data.size(); j += 32) {
      const auto values =
          _mm256_loadu_si256(reinterpret_cast<const __m256i*>(data.data() + j));
      sum = _mm256_xor_si256(sum, values);
      for (size_t bit = 0; bit < 8; ++bit) {
        const auto mask = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(masks[j / 32][bit].data()));
        groups[bit] =
            _mm256_xor_si256(groups[bit], _mm256_and_si256(values, mask));
      }
    }
    const auto reduce = [](__m256i value) {
      auto half = _mm_xor_si128(_mm256_castsi256_si128(value),
                                _mm256_extracti128_si256(value, 1));
      half = _mm_xor_si128(half, _mm_srli_si128(half, 8));
      half = _mm_xor_si128(half, _mm_srli_si128(half, 4));
      half = _mm_xor_si128(half, _mm_srli_si128(half, 2));
      half = _mm_xor_si128(half, _mm_srli_si128(half, 1));
      return static_cast<Element>(_mm_cvtsi128_si32(half));
    };
    s1 = reduce(sum);
    for (size_t bit = 0; bit < 8; ++bit) {
      const auto value = reduce(groups[bit]);
      const auto& products = tables[1u << bit];
      s0 ^= products[value & 15] ^ products[32 + (value >> 4)];
    }
  }
#else
  (void)vector_reduce;
#endif
  for (; j < data.size(); ++j) {
    const Element value = data[j];
    const auto& products = tables[(j + 2) ^ 1];
    s0 ^= products[value & 15] ^ products[32 + (value >> 4)];
    s1 ^= value;
  }
  return {s0, s1};
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
      valid_(Aligned(strong_n, strong_k) && weak_n <= 256 && weak_k >= 2 &&
             weak_k < weak_n &&
             (weak_n - weak_k == 2 || (weak_n == 256 && weak_k == 252))),
      strong_encoder_(valid_ ? strong_k : 0, valid_ ? strong_n - strong_k : 0),
      weak_encoder_(valid_ ? weak_k : 0, valid_ ? weak_n - weak_k : 0),
      strong_decoder_(valid_ ? strong_k : 0, valid_ ? strong_n - strong_k : 0),
      weak_decoder_(valid_ ? std::bit_ceil(weak_n) - (weak_n - weak_k) : 0,
                    valid_ ? weak_n - weak_k : 0) {}

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
  std::vector<Element> workspace(strong_encoder_.WorkspaceSize(weak_n_));
  std::vector<Element> weak_workspace(
      weak_n_ - weak_k_ == 4 ? weak_encoder_.WorkspaceSize(1) : 0);
  for (size_t row = 0; row < strong_k_; ++row) {
    if (weak_n_ - weak_k_ == 4) {
      for (size_t j = 0; j < weak_k_; ++j) {
        data[j] = &candidate[row * weak_n_ + j];
      }
      for (size_t j = 0; j < 4; ++j) {
        recovery[j] = &candidate[row * weak_n_ + weak_k_ + j];
      }
      const auto status = weak_encoder_.Encode(std::span(data).first(weak_k_),
                                               std::span(recovery).first(4), 1,
                                               weak_workspace, backend);
      if (status != lch::Status::ok) {
        return status;
      }
      continue;
    }
    const auto [s0, s1] =
        WeakDataSyndromes(std::span(candidate).subspan(row * weak_n_, weak_k_));
    candidate[row * weak_n_ + weak_k_] = s0;
    candidate[row * weak_n_ + weak_k_ + 1] = s0 ^ s1;
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
  return CorrectImpl(block, options, batches, true, true);
}

CorrectionResult StrongWeakRSProductCode::WeakCandidateR4(
    std::span<const Element> row,
    std::array<size_t, 2>& positions,
    std::array<Element, 2>& magnitudes,
    bool locate) const {
  // Full-field evaluation code: parity at native 0..3, data at 4..255.
  // These are ordinary power moments, not the novel-basis syndrome entries.
  const auto& tables = Tables().shuffle;
  const auto mul = [&](Element a, Element b) -> Element {
    return tables[a][b & 15] ^ tables[a][32 + (b >> 4)];
  };
  std::array<Element, 4> s{};
  for (size_t pos = 0; pos < 256; ++pos) {
    const Element x = static_cast<Element>((pos + 4) % 256);
    Element value = row[pos];
    for (size_t j = 0; j < 4; ++j) {
      s[j] ^= value;
      if (j != 3) {
        value = mul(value, x);
      }
    }
  }
  if (s == std::array<Element, 4>{}) {
    return {CorrectionStatus::ok, 0};
  }
  const CorrectionResult failure{CorrectionStatus::uncorrectable, 0};
  if (!locate) {
    return failure;
  }
  const Element determinant = mul(s[1], s[1]) ^ mul(s[0], s[2]);
  std::array<Element, 2> roots{};
  size_t count = 1;
  if (determinant == 0) {
    if (s[0] == 0) {
      return failure;
    }
    roots[0] = mul(s[1], InvCantor(s[0]));
    magnitudes[0] = s[0];
  } else {
    const Element inverse = InvCantor(determinant);
    const Element a = mul(mul(s[1], s[2]) ^ mul(s[0], s[3]), inverse);
    const Element b = mul(mul(s[1], s[3]) ^ mul(s[2], s[2]), inverse);
    if (a == 0) {
      return failure;  // Repeated roots cannot describe two errors.
    }
    // -1 distinguishes insoluble values from the valid representative zero.
    static const auto artin_schreier = [] {
      std::array<int16_t, 256> result{};
      result.fill(-1);
      for (unsigned y = 0; y < 256; ++y) {
        result[MultiplyCantor(y, y) ^ y] = static_cast<int16_t>(y);
      }
      return result;
    }();
    const Element inv_a = InvCantor(a);
    const int y = artin_schreier[mul(b, mul(inv_a, inv_a))];
    if (y < 0) {
      return failure;
    }
    roots[0] = mul(a, static_cast<Element>(y));
    roots[1] = roots[0] ^ a;
    magnitudes[0] = mul(s[1] ^ mul(s[0], roots[1]), inv_a);
    magnitudes[1] = s[0] ^ magnitudes[0];
    count = 2;
  }
  for (size_t i = 0; i < count; ++i) {
    if (magnitudes[i] == 0) {
      return failure;
    }
    positions[i] = roots[i] < 4 ? 252 + roots[i] : roots[i] - 4;
    Element value = magnitudes[i];
    for (size_t j = 0; j < 4; ++j) {
      s[j] ^= value;
      value = mul(value, roots[i]);
    }
  }
  if (s != std::array<Element, 4>{}) {
    return failure;
  }
  return {CorrectionStatus::ok, count};
}

CorrectionResult StrongWeakRSProductCode::WeakCandidate(
    std::span<const Element> row,
    size_t& position,
    Element& magnitude,
    bool locate,
    bool vector_reduce) const {
  // Native parity positions are 0,1; public data position j is native j+2.
  // At K=R=2 the low-rate family is the same degree-one evaluation code:
  // translating every evaluation point by 2 leaves its parity checks unchanged.
  auto [s0, s1] = WeakDataSyndromes(row.first(weak_k_), vector_reduce);
  s0 ^= row[weak_k_];
  s1 ^= row[weak_k_] ^ row[weak_k_ + 1];
  if (s0 == 0 && s1 == 0) {
    return {CorrectionStatus::ok, 0};
  }
  if (!locate || s1 == 0) {
    return {CorrectionStatus::uncorrectable, 0};
  }
  const size_t native = MultiplyCantor(s0, InvCantor(s1)) ^ 1;
  if (native >= std::bit_ceil(weak_n_) ||
      (native >= 2 && native - 2 >= weak_k_)) {
    return {CorrectionStatus::uncorrectable, 0};
  }
  // XORing s1 here cancels both checks exactly, verifying the entire candidate.
  position = native < 2 ? weak_k_ + native : native - 2;
  magnitude = s1;
  return {CorrectionStatus::ok, 1};
}

ProductCorrectionResult StrongWeakRSProductCode::CorrectImpl(
    std::span<Element> block,
    ProductDecodeOptions options,
    unsigned batch_passes,
    bool tracked_validation,
    bool direct_weak,
    unsigned optimizations) const {
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
#if defined(__AVX2__)
  const bool sparse_masks =
      (optimizations & 1) && lch::BackendAvailable(lch::Backend::avx2);
#endif
  const size_t mother_n = std::bit_ceil(weak_n_);
  const size_t mother_k = mother_n - (weak_n_ - weak_k_);
  const auto weak_position = [&](size_t pos) {
    return pos < weak_k_ ? pos : mother_k + pos - weak_k_;
  };
  // Pad independent batch lanes, never code positions. Value-initialized
  // synthetic columns stay zero; only real lanes publish results or repairs.
  const size_t strong_stride =
      direct_weak && (optimizations & 4) && batch_passes != 0
          ? (weak_n_ + 31) / 32 * 32
          : weak_n_;
  std::vector<Element> packed(
      batch_passes != 0 ? strong_n_ * (direct_weak ? strong_stride : mother_n)
                        : 0);
  std::vector<uint8_t> masks(packed.size());
  std::array<CorrectionResult, 256> outcomes{};
  std::array<Element*, 256> shards{};
  for (size_t pass = 0; pass < max_directional_passes; ++pass) {
    const bool strong = pass % 2 == 0;
    const size_t lines = strong ? weak_n_ : strong_n_;
    const size_t length = strong ? strong_n_ : weak_n_;
    const size_t decoder_length = strong ? length : mother_n;
    const size_t stride = strong ? strong_stride : lines;
    std::array<bool, 256> next{};
    size_t changes = 0;
    size_t bit_changes = 0;
    const bool batched =
        pass < std::min(batch_passes, 2u) && (strong || !direct_weak);
    if (batched) {
      // Columns already have position-major layout; rows need transposition.
      // Keep tentative repairs private until the existing weak gates accept.
      if (strong) {
        if (stride == lines) {
          std::copy(block.begin(), block.end(), packed.begin());
        } else {
          for (size_t pos = 0; pos < length; ++pos) {
            std::copy_n(block.data() + pos * lines, lines,
                        packed.data() + pos * stride);
          }
        }
      }
      for (size_t pos = 0; pos < decoder_length; ++pos) {
        shards[pos] = packed.data() + pos * stride;
        if (!strong) {
          for (size_t line = 0; line < lines; ++line) {
            shards[pos][line] =
                pos >= weak_k_ && pos < mother_k
                    ? Element{0}
                    : block[line * weak_n_ +
                            (pos < weak_k_ ? pos : weak_k_ + pos - mother_k)];
          }
        }
      }
      const auto status = detail::error_correction::CorrectCodewordBatch(
          strong ? strong_decoder_ : weak_decoder_,
          std::span(shards).first(decoder_length), stride,
          std::span(outcomes).first(stride),
          std::span(masks).first(decoder_length * stride));
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
        const auto commit = [&](size_t line) {
          const size_t index = pos * lines + line;
          const size_t staged = pos * stride + line;
          bit_changes += std::popcount(
              static_cast<unsigned>(block[index] ^ packed[staged]));
          block[index] = packed[staged];
          next[pos] = true;
          clean_rows[pos] = false;
          ++changes;
        };
        size_t line = 0;
#if defined(__AVX2__)
        if (sparse_masks) {
          for (; line + 32 <= lines; line += 32) {
            const auto mask =
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                    masks.data() + pos * stride + line));
            auto bits = ~static_cast<uint32_t>(_mm256_movemask_epi8(
                _mm256_cmpeq_epi8(mask, _mm256_setzero_si256())));
            while (bits != 0) {
              commit(line + std::countr_zero(bits));
              bits &= bits - 1;
            }
          }
        }
#endif
        for (; line < lines; ++line) {
          if (masks[pos * stride + line]) {
            commit(line);
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
        if (!strong && direct_weak) {
          std::array<size_t, 2> positions{};
          std::array<Element, 2> magnitudes{};
          const auto correction =
              weak_n_ - weak_k_ == 4
                  ? WeakCandidateR4(block.subspan(line * weak_n_, weak_n_),
                                    positions, magnitudes)
                  : WeakCandidate(block.subspan(line * weak_n_, weak_n_),
                                  positions[0], magnitudes[0], true,
                                  optimizations & 2);
          auto& clean = clean_rows[line];
          clean = correction.status == CorrectionStatus::ok &&
                  correction.error_count == 0;
          if (correction.status != CorrectionStatus::ok || clean) {
            continue;
          }
          bool accept = true;
          for (size_t i = 0; i < correction.error_count; ++i) {
            if ((options.use_anchors && protected_columns[positions[i]]) ||
                (options.use_binary_image &&
                 std::popcount(static_cast<unsigned>(magnitudes[i])) > 2)) {
              accept = false;
            }
          }
          if (!accept) {
            continue;
          }
          for (size_t i = 0; i < correction.error_count; ++i) {
            block[line * weak_n_ + positions[i]] ^= magnitudes[i];
            bit_changes += std::popcount(static_cast<unsigned>(magnitudes[i]));
            ++changes;
            next[positions[i]] = true;
            clean_columns[positions[i]] = false;
          }
          clean = true;
          continue;
        }
        const auto index = [&](size_t pos) {
          return strong ? pos * weak_n_ + line : line * weak_n_ + pos;
        };
        if (!batched) {
          candidate.fill(0);
          for (size_t pos = 0; pos < length; ++pos) {
            candidate[strong ? pos : weak_position(pos)] = block[index(pos)];
          }
        }
        const auto correction =
            batched
                ? outcomes[line]
                : CorrectCodeword(strong ? strong_decoder_ : weak_decoder_,
                                  std::span(candidate).first(decoder_length));
        if (batched) {
          for (size_t pos = 0; pos < decoder_length; ++pos) {
            candidate[pos] = shards[pos][line];
          }
        }
        // A mother-code repair is not a shortened-code candidate if it
        // changes any known-zero data. Reject before validity or gate updates.
        const bool shortened_valid =
            strong || std::all_of(candidate.begin() + weak_k_,
                                  candidate.begin() + mother_k,
                                  [](Element value) { return value == 0; });
        auto& clean = strong ? clean_columns[line] : clean_rows[line];
        clean = shortened_valid && correction.status == CorrectionStatus::ok &&
                correction.error_count == 0;
        if (strong) {
          protected_columns[line] = correction.status == CorrectionStatus::ok;
        }
        if (!shortened_valid || correction.status != CorrectionStatus::ok) {
          continue;
        }
        if (correction.error_count == 0) {
          continue;
        }
        if (!strong) {
          if (correction.error_count > (weak_n_ - weak_k_) / 2) {
            continue;
          }
          bool accept = true;
          for (size_t pos = 0; pos < length; ++pos) {
            const unsigned delta =
                block[index(pos)] ^ candidate[weak_position(pos)];
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
          const auto value = candidate[strong ? pos : weak_position(pos)];
          if (block[index(pos)] != value) {
            bit_changes +=
                std::popcount(static_cast<unsigned>(block[index(pos)] ^ value));
            block[index(pos)] = value;
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
      if (!strong && direct_weak) {
        std::array<size_t, 2> positions{};
        std::array<Element, 2> magnitudes{};
        const auto check =
            weak_n_ - weak_k_ == 4
                ? WeakCandidateR4(block.subspan(line * weak_n_, weak_n_),
                                  positions, magnitudes, false)
                : WeakCandidate(block.subspan(line * weak_n_, weak_n_),
                                positions[0], magnitudes[0], false,
                                optimizations & 2);
        if (check.status != CorrectionStatus::ok || check.error_count != 0) {
          result.all_zero_syndromes = false;
        }
        continue;
      }
      candidate.fill(0);
      for (size_t pos = 0; pos < length; ++pos) {
        candidate[strong ? pos : weak_position(pos)] =
            block[strong ? pos * weak_n_ + line : line * weak_n_ + pos];
      }
      const auto check = CorrectCodeword(
          strong ? strong_decoder_ : weak_decoder_,
          std::span(candidate).first(strong ? length : mother_n));
      if (check.status != CorrectionStatus::ok || check.error_count != 0) {
        result.all_zero_syndromes = false;
      }
    }
  }
  return result;
}

}  // namespace gf2p8::rs
