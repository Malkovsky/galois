#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <limits>
#include <random>
#include <utility>
#include <vector>

#include "gtest/gtest.h"
#include "reed_solomon/error_correction/internal.h"
#include "reed_solomon/product_code_internal.h"
#include "reed_solomon/strong_weak_rs_product_code.h"

namespace {

using gf2p8::Element;
using gf2p8::lch::Backend;
using gf2p8::lch::Status;
using namespace gf2p8::rs;

std::vector<Element> Codeword(size_t n, size_t k, uint32_t seed) {
  LCHEncoder encoder(k, n - k);
  std::mt19937 random(seed);
  std::vector<Element> word(n);
  for (size_t i = 0; i < k; ++i) {
    word[i] = static_cast<Element>(random());
  }
  std::vector<const Element*> data(k);
  std::vector<Element*> recovery(n - k);
  for (size_t i = 0; i < k; ++i) {
    data[i] = &word[i];
  }
  for (size_t i = k; i < n; ++i) {
    recovery[i - k] = &word[i];
  }
  std::vector<Element> workspace(encoder.WorkspaceSize(1));
  EXPECT_EQ(encoder.Encode(data, recovery, 1, workspace, Backend::scalar),
            Status::ok);
  return word;
}

TEST(WholeCodeword, RepairsEveryPositionAndFullRadiusIncludingParity) {
  for (const auto [n, k] : {std::pair{4u, 2u},
                            {8u, 4u},
                            {16u, 12u},
                            {256u, 224u},
                            {256u, 128u},
                            {256u, 254u}}) {
    const auto original = Codeword(n, k, n + k);
    LCHDecoder decoder(k, n - k);
    for (size_t pos = 0; pos < n; ++pos) {
      auto word = original;
      word[pos] ^= 0xff;
      const auto result = CorrectCodeword(decoder, word);
      ASSERT_EQ(result.status, CorrectionStatus::ok) << n << ':' << pos;
      EXPECT_EQ(result.error_count, 1u);
      ASSERT_EQ(word, original);
    }
    std::mt19937 random(n);
    for (size_t trial = 0; trial < 20; ++trial) {
      auto word = original;
      std::vector<size_t> positions(n);
      for (size_t i = 0; i < n; ++i) {
        positions[i] = i;
      }
      if (trial != 0) {
        std::shuffle(positions.begin(), positions.end(), random);
      } else {
        std::reverse(positions.begin(), positions.end());
      }
      for (size_t i = 0; i < (n - k) / 2; ++i) {
        word[positions[i]] ^= static_cast<Element>(1 + random() % 255);
      }
      const auto result = CorrectCodeword(decoder, word);
      ASSERT_EQ(result.status, CorrectionStatus::ok);
      EXPECT_EQ(result.error_count, (n - k) / 2);
      EXPECT_EQ(word, original);
    }
    auto clean = original;
    EXPECT_EQ(CorrectCodeword(decoder, clean).error_count, 0u);
    EXPECT_EQ(clean, original);
  }
}

TEST(WholeCodeword, TransactionalFailureAndInvalidDimensions) {
  LCHDecoder decoder(6, 2);
  auto word = Codeword(8, 6, 42);
  // Equal magnitudes cancel the leading syndrome: no one-error candidate.
  word[0] ^= 7;
  word[7] ^= 7;
  const auto before = word;
  const auto result = CorrectCodeword(decoder, word);
  EXPECT_EQ(result.status, CorrectionStatus::uncorrectable);
  EXPECT_EQ(result.error_count, 0u);
  EXPECT_EQ(word, before);
  EXPECT_EQ(CorrectCodeword(decoder, std::span(word).first(7)).status,
            CorrectionStatus::invalid_argument);
  EXPECT_EQ(CorrectCodeword(LCHDecoder(0, 0), word).status,
            CorrectionStatus::invalid_argument);
  EXPECT_EQ(CorrectCodeword(LCHDecoder(5, 3), word).status,
            CorrectionStatus::unsupported_dimensions);
  EXPECT_EQ(word, before);
  std::vector<Element> shortened(7, 1);
  const auto saved = shortened;
  EXPECT_EQ(CorrectCodeword(LCHDecoder(5, 2), shortened).status,
            CorrectionStatus::unsupported_dimensions);
  EXPECT_EQ(shortened, saved);
}

TEST(ProductCode, DimensionsAndInvalidCalls) {
  EXPECT_TRUE(StrongWeakRSProductCode().Valid());
  EXPECT_EQ(StrongWeakRSProductCode().BlockSize(), 65536u);
  for (const auto [n, k] : {std::pair{0u, 0u},
                            {8u, 8u},
                            {8u, 9u},
                            {7u, 5u},
                            {8u, 5u},
                            {8u, 2u},
                            {512u, 480u},
                            {8u, 7u}}) {
    StrongWeakRSProductCode code(n, k, 8, 6);
    EXPECT_FALSE(code.Valid());
    EXPECT_EQ(code.BlockSize(), 0u);
    std::vector<Element> block(32, 17);
    const auto before = block;
    EXPECT_EQ(code.Encode(block), Status::invalid_argument);
    EXPECT_EQ(code.Correct(block).termination,
              ProductTermination::invalid_argument);
    EXPECT_EQ(block, before);
  }
  EXPECT_FALSE(StrongWeakRSProductCode(8, 4, 8, 4).Valid());
  EXPECT_FALSE(
      StrongWeakRSProductCode(std::numeric_limits<size_t>::max(), 1).Valid());
  StrongWeakRSProductCode code(4, 2, 8, 6);
  std::vector<Element> block(32, 17);
  const auto before = block;
  for (size_t cap : {0u, 1u}) {
    const auto result = code.Correct(block, cap);
    EXPECT_EQ(result.termination, ProductTermination::invalid_argument);
    EXPECT_EQ(result.directional_passes, 0u);
    EXPECT_FALSE(result.all_zero_syndromes);
    EXPECT_EQ(block, before);
  }
  EXPECT_EQ(code.Correct(std::span(block).first(31)).termination,
            ProductTermination::invalid_argument);
  EXPECT_EQ(code.Encode(std::span(block).first(31)), Status::invalid_argument);
  EXPECT_EQ(block, before);
}

TEST(ProductCode, SystematicEncodingScalarAgreementAndAllComponentValidity) {
  for (const auto [ns, ks, nw, kw] : {std::array<size_t, 4>{4, 2, 8, 6},
                                      {16, 12, 16, 14},
                                      {256, 224, 256, 254}}) {
    StrongWeakRSProductCode code(ns, ks, nw, kw);
    std::mt19937 random(901);
    std::vector<Element> block(code.BlockSize());
    for (auto& value : block) {
      value = static_cast<Element>(random());
    }
    const auto input = block;
    auto scalar = block;
    ASSERT_EQ(code.Encode(block), Status::ok);
    ASSERT_EQ(code.Encode(scalar, Backend::scalar), Status::ok);
    EXPECT_EQ(block, scalar);
    for (size_t row = 0; row < ks; ++row) {
      for (size_t col = 0; col < kw; ++col) {
        ASSERT_EQ(block[row * nw + col], input[row * nw + col]);
      }
    }
    const auto result = code.Correct(block);
    EXPECT_TRUE(result.all_zero_syndromes);
    EXPECT_EQ(result.directional_passes, 2u);
    EXPECT_EQ(result.strong_lines_visited, nw);
    EXPECT_EQ(result.weak_lines_visited, ns);
    EXPECT_EQ(result.changed_symbols, 0u);
    EXPECT_EQ(result.termination, ProductTermination::no_change);
    EXPECT_EQ(block, scalar);
    // Arbitrary symbol damage in each of the four product regions.
    for (auto [row, col] : {std::pair{size_t{0}, size_t{0}},
                            {ks, size_t{0}},
                            {size_t{0}, kw},
                            {ks, kw}}) {
      block[row * nw + col] ^= 0xff;
      const auto repaired = code.Correct(block);
      EXPECT_TRUE(repaired.all_zero_syndromes);
      EXPECT_EQ(repaired.changed_symbols, 1u);
      EXPECT_EQ(repaired.changed_bits, 8u);
      EXPECT_EQ(repaired.strong_changed_bits, 8u);
      EXPECT_EQ(repaired.weak_changed_bits, 0u);
      EXPECT_EQ(block, scalar);
    }
  }
}

TEST(ProductCode, InitialWeakPassSelectiveActivationAndCap) {
  StrongWeakRSProductCode code(4, 2, 8, 6);
  for (size_t col : {0u, 6u, 7u}) {
    for (Element magnitude : {Element{1}, Element{3}}) {
      std::vector<Element> damaged(32);
      // Strong cannot repair two equal errors; weak can repair both rows,
      // including strong-parity rows and the parity/parity corner.
      damaged[2 * 8 + col] = magnitude;
      damaged[3 * 8 + col] = magnitude;
      auto capped = damaged;
      const auto cap = code.Correct(capped, 2);
      EXPECT_EQ(cap.termination, ProductTermination::pass_limit);
      EXPECT_EQ(cap.directional_passes, 2u);
      EXPECT_TRUE(cap.all_zero_syndromes);
      EXPECT_EQ(capped, std::vector<Element>(32));
      const auto result = code.Correct(damaged);
      EXPECT_EQ(result.directional_passes, 3u);
      EXPECT_EQ(result.strong_lines_visited, 9u);
      EXPECT_EQ(result.weak_lines_visited, 4u);
      EXPECT_EQ(result.changed_symbols, 2u);
      EXPECT_EQ(result.termination, ProductTermination::no_change);
      EXPECT_TRUE(result.all_zero_syndromes);
      EXPECT_EQ(damaged, capped);
    }
  }
}

TEST(ProductCode, DefaultStrongRadiusRepairsAllColumnsIncludingParity) {
  StrongWeakRSProductCode code;
  std::vector<Element> block(code.BlockSize());
  std::mt19937 random(1983);
  for (auto& value : block) {
    value = static_cast<Element>(random());
  }
  ASSERT_EQ(code.Encode(block, Backend::scalar), Status::ok);
  const auto original = block;
  for (size_t col = 0; col < 256; ++col) {
    for (size_t error = 0; error < 16; ++error) {
      const size_t row = (col + error * 17) % 256;
      block[row * 256 + col] ^= static_cast<Element>(1 + random() % 255);
    }
  }
  const auto result = code.Correct(block);
  EXPECT_EQ(result.directional_passes, 2u);
  EXPECT_EQ(result.changed_symbols, 4096u);
  EXPECT_TRUE(result.all_zero_syndromes);
  EXPECT_EQ(block, original);
}

TEST(ProductCode, WeakBitGateRejectsTransactionallyWithoutActivation) {
  StrongWeakRSProductCode code(4, 2, 8, 6);
  for (Element magnitude : {Element{7}, Element{255}}) {
    std::vector<Element> block(32);
    block[6] = magnitude;
    block[3 * 8 + 6] = magnitude;
    const auto before = block;
    const auto result = code.Correct(block);
    EXPECT_FALSE(result.all_zero_syndromes);
    EXPECT_EQ(result.termination, ProductTermination::no_change);
    EXPECT_EQ(result.directional_passes, 2u);
    EXPECT_EQ(result.changed_symbols, 0u);
    EXPECT_EQ(block, before);
  }
}

TEST(ProductCode, ChangesPropagateThroughFourDirectionalPasses) {
  StrongWeakRSProductCode code(4, 2, 8, 6);
  std::vector<Element> block(32);
  block[0] = block[8] = block[9] = block[17] = 1;
  // Both columns initially fail, and the middle row fails. Weak fixes the
  // outer two rows; strong then fixes the middle row in both active columns.
  auto capped = block;
  const auto limit = code.Correct(capped, 2);
  EXPECT_EQ(limit.termination, ProductTermination::pass_limit);
  EXPECT_FALSE(limit.all_zero_syndromes);
  EXPECT_EQ(limit.changed_symbols, 2u);
  EXPECT_EQ(capped[8], 1);
  EXPECT_EQ(capped[9], 1);
  const auto result = code.Correct(block);
  EXPECT_TRUE(result.all_zero_syndromes);
  EXPECT_EQ(result.termination, ProductTermination::no_change);
  EXPECT_EQ(result.directional_passes, 4u);
  EXPECT_EQ(result.strong_lines_visited, 10u);
  EXPECT_EQ(result.weak_lines_visited, 5u);
  EXPECT_EQ(result.changed_symbols, 4u);
  EXPECT_EQ(block, std::vector<Element>(32));
}

TEST(ProductCode, SuccessfulStrongRepairProtectsItsColumn) {
  StrongWeakRSProductCode code(4, 2, 8, 6);
  std::vector<Element> block(32);
  for (size_t row = 0; row < 4; ++row) {
    block[row * 8] = 1;
  }
  const auto protected_word = block;
  block[0] ^= 2;
  const auto result = code.Correct(block);
  EXPECT_FALSE(result.all_zero_syndromes);
  EXPECT_EQ(result.changed_symbols, 1u);
  EXPECT_EQ(result.directional_passes, 2u);
  EXPECT_EQ(block, protected_word);
}

TEST(ProductCode, UnvisitedColumnRetainsProtectionOnLaterWeakPass) {
  StrongWeakRSProductCode code(8, 4, 8, 6);
  // A valid strong column with systematic symbols [0,1,0,0]. Obtain its
  // parity using the owned scalar encoder, independently of product Encode.
  LCHEncoder encoder(4, 4);
  std::array<Element, 8> protected_column{0, 1, 0, 0};
  std::array<const Element*, 4> data{};
  std::array<Element*, 4> recovery{};
  for (size_t i = 0; i < 4; ++i) {
    data[i] = &protected_column[i];
    recovery[i] = &protected_column[4 + i];
  }
  std::vector<Element> workspace(encoder.WorkspaceSize(1));
  ASSERT_EQ(encoder.Encode(data, recovery, 1, workspace, Backend::scalar),
            Status::ok);
  std::vector<Element> block(64);
  for (size_t row = 0; row < 8; ++row) {
    block[row * 8 + 1] = protected_column[row];
  }
  const auto expected = block;
  for (size_t row = 0; row < 4; ++row) {
    block[row * 8] = 1;
  }
  const auto result = code.Correct(block);
  EXPECT_EQ(result.directional_passes, 4u);
  EXPECT_EQ(result.strong_lines_visited, 9u);
  EXPECT_EQ(result.weak_lines_visited, 9u);
  EXPECT_EQ(result.changed_symbols, 4u);
  EXPECT_FALSE(result.all_zero_syndromes);
  EXPECT_EQ(block, expected);
}

TEST(WholeCodeword, OverRadiusOutcomesAreTransactionalOrVerifiedCodewords) {
  std::mt19937 random(8181);
  for (const auto [n, k] : {std::pair{8u, 6u}, {16u, 12u}, {32u, 16u}}) {
    LCHDecoder decoder(k, n - k);
    LCHEncoder encoder(k, n - k);
    for (size_t trial = 0; trial < 100; ++trial) {
      auto word = Codeword(n, k, random());
      for (size_t pos = 0; pos <= (n - k) / 2; ++pos) {
        word[n - 1 - pos] ^= static_cast<Element>(1 + random() % 255);
      }
      const auto before = word;
      const auto result = CorrectCodeword(decoder, word);
      if (result.status != CorrectionStatus::ok) {
        EXPECT_EQ(word, before);
        EXPECT_EQ(result.error_count, 0u);
        continue;
      }
      size_t distance = 0;
      for (size_t pos = 0; pos < n; ++pos) {
        distance += word[pos] != before[pos];
      }
      EXPECT_EQ(distance, result.error_count);
      EXPECT_LE(distance, (n - k) / 2);
      std::vector<const Element*> data(k);
      std::vector<Element> parity(n - k);
      std::vector<Element*> recovery(n - k);
      for (size_t i = 0; i < k; ++i) {
        data[i] = &word[i];
      }
      for (size_t i = 0; i < n - k; ++i) {
        recovery[i] = &parity[i];
      }
      std::vector<Element> workspace(encoder.WorkspaceSize(1));
      ASSERT_EQ(encoder.Encode(data, recovery, 1, workspace, Backend::scalar),
                Status::ok);
      EXPECT_TRUE(std::equal(parity.begin(), parity.end(), word.begin() + k));
    }
  }
}

TEST(ProductCode, ProtectedColumnRejectsEvenLowWeightWeakCandidate) {
  StrongWeakRSProductCode code(4, 2, 8, 6);
  // Constant columns are valid RS words. Every row proposes a one-bit repair
  // at column 0, but its zero-syndrome strong protection must reject them.
  std::vector<Element> block(32);
  for (size_t row = 0; row < 4; ++row) {
    block[row * 8] = 1;
  }
  const auto before = block;
  const auto result = code.Correct(block);
  EXPECT_FALSE(result.all_zero_syndromes);
  EXPECT_EQ(result.changed_symbols, 0u);
  EXPECT_EQ(result.directional_passes, 2u);
  EXPECT_EQ(block, before);
}

TEST(WholeCodewordBatch, DifferentialDataParityFailuresAndTails) {
  using detail::error_correction::CorrectCodewordBatch;
  std::mt19937 random(0xba7c224);
  for (const auto [n, k] : {std::pair{4u, 2u},
                            {8u, 4u},
                            {16u, 12u},
                            {256u, 128u},
                            {256u, 224u},
                            {256u, 254u}}) {
    LCHDecoder decoder(k, n - k);
    for (size_t lanes : {1u, 31u, 32u, 33u, 65u, 256u}) {
      std::vector<Element> packed(n * lanes);
      auto expected = packed;
      std::vector<CorrectionResult> results(lanes), reference(lanes);
      std::vector<uint8_t> masks(packed.size(), 0xff);
      std::vector<Element*> shards(n);
      for (size_t pos = 0; pos < n; ++pos) {
        shards[pos] = &packed[pos * lanes];
      }
      for (size_t lane = 0; lane < lanes; ++lane) {
        auto word = Codeword(n, k, random());
        const size_t radius = (n - k) / 2;
        const size_t errors = lane % 6 == 0   ? 0
                              : lane % 6 == 1 ? 1
                              : lane % 6 == 2 ? radius
                              : lane % 6 == 3 ? radius + 1
                              : lane % 6 == 4 ? n
                                              : radius;
        std::vector<size_t> positions(n);
        for (size_t i = 0; i < n; ++i) {
          positions[i] = i;
        }
        std::shuffle(positions.begin(), positions.end(), random);
        // Include parity-only full-radius candidates in every vector chunk.
        if (lane % 6 == 5) {
          std::sort(positions.rbegin(), positions.rend());
        }
        for (size_t i = 0; i < errors; ++i) {
          word[positions[i]] ^= static_cast<Element>(1 + random() % 255);
        }
        for (size_t pos = 0; pos < n; ++pos) {
          packed[pos * lanes + lane] = word[pos];
        }
        const auto before = word;
        reference[lane] = CorrectCodeword(decoder, word);
        if (reference[lane].status != CorrectionStatus::ok) {
          EXPECT_EQ(word, before);
        }
        for (size_t pos = 0; pos < n; ++pos) {
          expected[pos * lanes + lane] = word[pos];
        }
      }
      const auto before = packed;
      ASSERT_EQ(CorrectCodewordBatch(decoder, shards, lanes, results, masks),
                CorrectionStatus::ok);
      ASSERT_EQ(packed, expected) << n << ':' << lanes;
      for (size_t lane = 0; lane < lanes; ++lane) {
        EXPECT_EQ(results[lane].status, reference[lane].status) << lane;
        EXPECT_EQ(results[lane].error_count, reference[lane].error_count)
            << lane;
        for (size_t pos = 0; pos < n; ++pos) {
          const auto index = pos * lanes + lane;
          EXPECT_EQ(masks[index], before[index] != packed[index]);
        }
      }
    }
  }
}

TEST(WholeCodewordBatch, InvalidRangesAreUntouched) {
  using detail::error_correction::CorrectCodewordBatch;
  LCHDecoder decoder(6, 2);
  std::vector<Element> packed(8 * 33, 7);
  std::array<Element*, 8> shards{};
  for (size_t i = 0; i < 8; ++i) {
    shards[i] = packed.data() + i * 33;
  }
  std::vector<CorrectionResult> results(33);
  std::vector<uint8_t> masks(packed.size(), 42);
  const auto before = packed;
  shards[7] = shards[0];
  EXPECT_EQ(CorrectCodewordBatch(decoder, shards, 33, results, masks),
            CorrectionStatus::invalid_argument);
  shards[7] = nullptr;
  EXPECT_EQ(CorrectCodewordBatch(decoder, shards, 33, results, masks),
            CorrectionStatus::invalid_argument);
  shards[7] = packed.data() + 7 * 33;
  EXPECT_EQ(CorrectCodewordBatch(decoder, shards, 33, results, packed),
            CorrectionStatus::invalid_argument);
  EXPECT_EQ(CorrectCodewordBatch(decoder, shards, 32, results, masks),
            CorrectionStatus::invalid_argument);
  EXPECT_EQ(CorrectCodewordBatch(LCHDecoder(5, 3), shards, 33, results, masks),
            CorrectionStatus::unsupported_dimensions);
  EXPECT_EQ(CorrectCodewordBatch(decoder, shards, 0, {}, {}),
            CorrectionStatus::ok);
  EXPECT_EQ(packed, before);
  EXPECT_EQ(masks, std::vector<uint8_t>(packed.size(), 42));
  for (const auto& result : results) {
    EXPECT_EQ(result.status, CorrectionStatus::invalid_argument);
  }
}

// Independent parity oracle: no decoder outcomes or scheduler bookkeeping.
bool AllComponentsValid(const std::vector<Element>& block,
                        size_t ns,
                        size_t ks,
                        size_t nw,
                        size_t kw) {
  bool valid = true;
  for (bool strong : {true, false}) {
    const size_t n = strong ? ns : nw, k = strong ? ks : kw;
    LCHEncoder encoder(k, n - k);
    std::vector<Element> parity(n - k);
    std::vector<const Element*> data(k);
    std::vector<Element*> recovery(n - k);
    std::vector<Element> workspace(encoder.WorkspaceSize(1));
    for (size_t i = 0; i < n - k; ++i) {
      recovery[i] = &parity[i];
    }
    for (size_t line = 0; line < (strong ? nw : ns); ++line) {
      const auto index = [&](size_t pos) {
        return strong ? pos * nw + line : line * nw + pos;
      };
      for (size_t i = 0; i < k; ++i) {
        data[i] = &block[index(i)];
      }
      EXPECT_EQ(encoder.Encode(data, recovery, 1, workspace, Backend::scalar),
                Status::ok);
      for (size_t i = k; i < n; ++i) {
        valid &= parity[i - k] == block[index(i)];
      }
    }
  }
  return valid;
}

TEST(ProductCode, InitialBatchChoicesMatchSingleOutputsAndAllCounters) {
  std::mt19937 random(0x5b5c0224);
  for (const auto dims : {std::array<size_t, 4>{4, 2, 8, 6},
                          {32, 16, 64, 62},
                          {256, 224, 256, 254}}) {
    const auto [ns, ks, nw, kw] = dims;
    StrongWeakRSProductCode code(ns, ks, nw, kw);
    for (size_t trial = 0; trial < 16; ++trial) {
      std::vector<Element> input(code.BlockSize());
      for (auto& value : input) {
        value = static_cast<Element>(random());
      }
      ASSERT_EQ(code.Encode(input), Status::ok);
      if (trial < 8) {
        for (auto& value : input) {
          for (unsigned bit = 0; bit < 8; ++bit) {
            if (random() % 200 == 0) {
              value ^= static_cast<Element>(1u << bit);
            }
          }
        }
      } else {
        // Valid strong columns propose weak repairs into protected columns;
        // overloaded columns exercise rejection, bit gates, and activation.
        std::fill(input.begin(), input.end(), Element{0});
        for (size_t row = 0; row < ns; ++row) {
          input[row * nw] = 1;
        }
        for (size_t row = 0; row <= (ns - ks) / 2; ++row) {
          input[row * nw + 1] = trial % 2 ? 7 : 3;
          if (row % 2) {
            input[row * nw + 2] = 1;
          }
        }
        if (trial % 3 == 0) {
          input[0] ^= 2;
        }
      }
      for (size_t cap : {2u, 3u, 4u, 5u, 6u, 16u}) {
        for (bool anchors : {false, true}) {
          for (bool binary : {false, true}) {
            const ProductDecodeOptions options{cap, anchors, binary};
            auto reference = input;
            const auto expected = detail::ProductCorrectionAccess::Correct(
                code, reference, options, 0, false);
            EXPECT_EQ(expected.all_zero_syndromes,
                      AllComponentsValid(reference, ns, ks, nw, kw));
            for (unsigned batches : {0u, 1u, 2u}) {
              auto actual = input;
              const auto result = detail::ProductCorrectionAccess::Correct(
                  code, actual, options, batches);
              ASSERT_EQ(actual, reference)
                  << ns << ':' << trial << ':' << cap << ':' << batches;
              EXPECT_EQ(result.termination, expected.termination);
              EXPECT_EQ(result.all_zero_syndromes, expected.all_zero_syndromes);
              EXPECT_EQ(result.directional_passes, expected.directional_passes);
              EXPECT_EQ(result.strong_lines_visited,
                        expected.strong_lines_visited);
              EXPECT_EQ(result.weak_lines_visited, expected.weak_lines_visited);
              EXPECT_EQ(result.changed_symbols, expected.changed_symbols);
              EXPECT_EQ(result.changed_bits, expected.changed_bits);
              EXPECT_EQ(result.strong_changed_bits,
                        expected.strong_changed_bits);
              EXPECT_EQ(result.weak_changed_bits, expected.weak_changed_bits);
              EXPECT_EQ(result.strong_changed_symbols,
                        expected.strong_changed_symbols);
              EXPECT_EQ(result.weak_changed_symbols,
                        expected.weak_changed_symbols);
            }
          }
        }
      }
    }
  }
}

TEST(ProductCode,
     TrackedValidityMatchesIndependentParityAcrossCapsAndCancellations) {
  std::mt19937 random(0xc1ea0224);
  for (const auto dims :
       {std::array<size_t, 4>{4, 2, 8, 6}, {8, 4, 4, 2}, {32, 28, 32, 30}}) {
    const auto [ns, ks, nw, kw] = dims;
    StrongWeakRSProductCode code(ns, ks, nw, kw);
    for (size_t trial = 0; trial < 128; ++trial) {
      std::vector<Element> input(code.BlockSize());
      // Include undetected valid words, equal-magnitude cancellations,
      // parity damage, dense failures, and both weak rejection gates.
      if (trial % 4 == 0) {
        for (auto& value : input) {
          value = static_cast<Element>(random());
        }
        ASSERT_EQ(code.Encode(input, Backend::scalar), Status::ok);
      }
      const size_t errors = trial % (ns * 2);
      for (size_t i = 0; i < errors; ++i) {
        input[random() % input.size()] ^=
            trial % 3 == 0   ? Element{1}
            : trial % 3 == 1 ? Element{7}
                             : static_cast<Element>(1 + random() % 255);
      }
      for (size_t cap : {2u, 3u, 4u, 5u, 6u, 7u, 8u, 16u}) {
        for (bool anchors : {false, true}) {
          for (bool binary : {false, true}) {
            const ProductDecodeOptions options{cap, anchors, binary};
            auto reference = input;
            const auto expected = detail::ProductCorrectionAccess::Correct(
                code, reference, options, 0, false);
            const bool valid = AllComponentsValid(reference, ns, ks, nw, kw);
            EXPECT_EQ(expected.all_zero_syndromes, valid);
            for (unsigned batches : {0u, 1u, 2u, 3u}) {
              auto actual = input;
              const auto result =
                  batches == 3 ? code.Correct(actual, options)
                               : detail::ProductCorrectionAccess::Correct(
                                     code, actual, options, batches);
              ASSERT_EQ(actual, reference)
                  << ns << ':' << trial << ':' << cap << ':' << batches;
              EXPECT_EQ(result.all_zero_syndromes, valid);
              EXPECT_EQ(result.termination, expected.termination);
              EXPECT_EQ(result.directional_passes, expected.directional_passes);
              EXPECT_EQ(result.strong_lines_visited,
                        expected.strong_lines_visited);
              EXPECT_EQ(result.weak_lines_visited, expected.weak_lines_visited);
              EXPECT_EQ(result.changed_symbols, expected.changed_symbols);
              EXPECT_EQ(result.changed_bits, expected.changed_bits);
              EXPECT_EQ(result.strong_changed_bits,
                        expected.strong_changed_bits);
              EXPECT_EQ(result.weak_changed_bits, expected.weak_changed_bits);
              EXPECT_EQ(result.strong_changed_symbols,
                        expected.strong_changed_symbols);
              EXPECT_EQ(result.weak_changed_symbols,
                        expected.weak_changed_symbols);
            }
          }
        }
      }
    }
  }
}

TEST(ProductCode, IndependentGatesProtectParityAndKeepSingleSymbolBDD) {
  for (size_t n : {4u, 32u}) {
    StrongWeakRSProductCode code(n, n - 2, n, n - 2);
    for (bool anchors : {false, true}) {
      for (bool binary : {false, true}) {
        for (unsigned batches : {0u, 1u, 2u, 3u}) {
          for (size_t cap : {2u, 16u}) {
            const ProductDecodeOptions options{cap, anchors, binary};
            const auto correct = [&](std::vector<Element>& block) {
              return batches == 3 ? code.Correct(block, options)
                                  : detail::ProductCorrectionAccess::Correct(
                                        code, block, options, batches);
            };
            for (Element magnitude : {Element{1}, Element{7}}) {
              // A clean strong parity column: each weak row proposes one
              // repair.
              std::vector<Element> block(n * n);
              for (size_t row = 0; row < n; ++row) {
                block[row * n + n - 1] = magnitude;
              }
              const auto before = block;
              const bool accept = !anchors && (!binary || magnitude == 1);
              const auto result = correct(block);
              EXPECT_EQ(block, accept ? std::vector<Element>(n * n) : before);
              EXPECT_EQ(result.changed_symbols, accept ? n : 0u);
              EXPECT_EQ(result.changed_bits,
                        accept ? n * std::popcount(magnitude) : 0u);
              EXPECT_EQ(result.weak_changed_bits, result.changed_bits);
              EXPECT_EQ(result.strong_changed_bits, 0u);
              EXPECT_EQ(result.all_zero_syndromes, accept);
              EXPECT_EQ(result.directional_passes, accept && cap > 2 ? 3u : 2u);
              EXPECT_EQ(result.termination,
                        accept && cap == 2 ? ProductTermination::pass_limit
                                           : ProductTermination::no_change);

              // Strong failure leaves this parity column unprotected. Only the
              // binary-image option may reject these weak repairs.
              block.assign(n * n, 0);
              block[(n - 2) * n + n - 1] = magnitude;
              block[(n - 1) * n + n - 1] = magnitude;
              const auto damaged = block;
              const bool repair = !binary || magnitude == 1;
              const auto unprotected = correct(block);
              EXPECT_EQ(block, repair ? std::vector<Element>(n * n) : damaged);
              EXPECT_EQ(unprotected.changed_symbols, repair ? 2u : 0u);
              EXPECT_EQ(unprotected.changed_bits,
                        repair ? 2u * std::popcount(magnitude) : 0u);
              EXPECT_EQ(unprotected.all_zero_syndromes, repair);
              EXPECT_EQ(unprotected.strong_lines_visited,
                        n + (repair && cap > 2 ? 1 : 0));
            }

            // Equal errors in a 2x2 corner defeat both one-symbol component
            // decoders even when both acceptance gates are disabled.
            std::vector<Element> block(n * n);
            for (size_t row : {n - 2, n - 1}) {
              for (size_t col : {n - 2, n - 1}) {
                block[row * n + col] = 1;
              }
            }
            const auto before = block;
            const auto result = correct(block);
            EXPECT_EQ(block, before);
            EXPECT_EQ(result.changed_symbols, 0u);
            EXPECT_FALSE(result.all_zero_syndromes);
          }
        }
      }
    }
  }
}

TEST(ProductCode, AnchorFreeWritesInvalidatePreviouslyCleanColumnsAtCap) {
  for (size_t n : {4u, 32u}) {
    StrongWeakRSProductCode code(n, n - 2, n, n - 2);
    for (bool binary : {false, true}) {
      for (unsigned batches : {0u, 1u, 2u, 3u}) {
        std::vector<Element> block(n * n);
        for (size_t row = 0; row < n; ++row) {
          block[row * n] = 1;
        }
        block[1] = block[n + 1] = 1;
        // Column 0 starts clean. Weak repairs only rows 2..N-1, making
        // column 0 invalid. Cached strong success must not hide those writes.
        const ProductDecodeOptions options{2, false, binary};
        const auto result = batches == 3
                                ? code.Correct(block, options)
                                : detail::ProductCorrectionAccess::Correct(
                                      code, block, options, batches);
        EXPECT_EQ(result.changed_symbols, n - 2);
        EXPECT_EQ(result.termination, ProductTermination::pass_limit);
        EXPECT_FALSE(result.all_zero_syndromes);
        EXPECT_FALSE(AllComponentsValid(block, n, n - 2, n, n - 2));
      }
    }
  }
}

TEST(ProductCode, CountsRepeatedCommittedWritesByIndependentPassDifferences) {
  StrongWeakRSProductCode code(4, 2, 8, 6);
  LCHDecoder strong(2, 2);
  std::mt19937 random(0xacc37);
  bool repeated = false;
  for (size_t trial = 0; trial < 512; ++trial) {
    std::vector<Element> input(32);
    for (size_t i = 0; i < 12; ++i) {
      input[random() % input.size()] ^= static_cast<Element>(1 + random() % 7);
    }
    auto previous = input;
    // Independently reproduce just the initial strong pass. Each following
    // cap extends the same deterministic prefix by one directional pass.
    for (size_t col = 0; col < 8; ++col) {
      std::array<Element, 4> column{};
      for (size_t row = 0; row < 4; ++row) {
        column[row] = input[row * 8 + col];
      }
      CorrectCodeword(strong, column);
      for (size_t row = 0; row < 4; ++row) {
        previous[row * 8 + col] = column[row];
      }
    }
    size_t bits = 0, symbols = 0, strong_bits = 0, strong_symbols = 0;
    std::array<size_t, 32> writes{};
    for (size_t pos = 0; pos < 32; ++pos) {
      bits += std::popcount(static_cast<unsigned>(input[pos] ^ previous[pos]));
      symbols += input[pos] != previous[pos];
      writes[pos] += input[pos] != previous[pos];
    }
    strong_bits = bits;
    strong_symbols = symbols;
    for (size_t cap = 2; cap <= 8; ++cap) {
      auto actual = input;
      const auto result =
          code.Correct(actual, ProductDecodeOptions{cap, false, false});
      for (size_t pos = 0; pos < 32; ++pos) {
        const size_t delta_bits =
            std::popcount(static_cast<unsigned>(previous[pos] ^ actual[pos]));
        const size_t delta_symbols = previous[pos] != actual[pos];
        bits += delta_bits;
        symbols += delta_symbols;
        if (cap % 2 == 1) {
          strong_bits += delta_bits;
          strong_symbols += delta_symbols;
        }
        writes[pos] += delta_symbols;
        repeated |= writes[pos] > 1;
      }
      EXPECT_EQ(result.changed_bits, bits);
      EXPECT_EQ(result.changed_symbols, symbols);
      EXPECT_EQ(result.strong_changed_bits, strong_bits);
      EXPECT_EQ(result.strong_changed_symbols, strong_symbols);
      EXPECT_EQ(result.weak_changed_bits, bits - strong_bits);
      EXPECT_EQ(result.weak_changed_symbols, symbols - strong_symbols);
      previous = actual;
    }
  }
  EXPECT_TRUE(repeated);
}

TEST(ProductCode, OptionsDefaultsAndInvalidCapsArePerCall) {
  StrongWeakRSProductCode code(4, 2, 8, 6);
  std::vector<Element> input(32);
  input[6] = input[30] = 7;
  for (bool anchors : {false, true}) {
    for (bool binary : {false, true}) {
      for (size_t cap : {0u, 1u}) {
        auto block = input;
        const auto result =
            code.Correct(block, ProductDecodeOptions{cap, anchors, binary});
        EXPECT_EQ(result.termination, ProductTermination::invalid_argument);
        EXPECT_EQ(result.directional_passes, 0u);
        EXPECT_FALSE(result.all_zero_syndromes);
        EXPECT_EQ(block, input);
      }
      auto block = input;
      code.Correct(block, ProductDecodeOptions{16, anchors, binary});
      auto defaults = input;
      auto legacy = input;
      const auto result = code.Correct(defaults, ProductDecodeOptions{});
      const auto old = code.Correct(legacy);
      EXPECT_EQ(defaults, input);
      EXPECT_EQ(defaults, legacy);
      EXPECT_EQ(result.termination, old.termination);
      EXPECT_EQ(result.all_zero_syndromes, old.all_zero_syndromes);
      EXPECT_EQ(result.directional_passes, old.directional_passes);
      EXPECT_EQ(result.changed_symbols, old.changed_symbols);
    }
  }
}

}  // namespace
