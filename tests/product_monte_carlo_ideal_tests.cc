#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <bit>
#include <vector>

#include "product_monte_carlo_trials.h"

namespace {
struct TrialResult {
  std::array<uint64_t, 22> metrics{};
  std::array<uint64_t, 3> supplemental{};
  std::vector<uint8_t> residual;
};

TrialResult Decode(const std::vector<uint8_t>& noise,
                   size_t n1 = 4,
                   size_t k1 = 2,
                   size_t n2 = 4,
                   size_t k2 = 2,
                   bool postprocessing = false,
                   size_t passes = 64,
                   bool random = false,
                   bool gates = true) {
  std::vector<uint32_t> positions;
  for (size_t i = 0; i < noise.size(); ++i) {
    for (size_t bit = 0; bit < 8; ++bit) {
      if (noise[i] & (1u << bit)) {
        positions.push_back(static_cast<uint32_t>(8 * i + bit));
      }
    }
  }
  TrialResult result;
  result.residual.resize(noise.size());
  uint32_t empty = 0;
  EXPECT_EQ(
      product_trial_decoder(42, 0, 0, positions.size(), passes, gates, gates,
                            result.metrics.data(), 2,
                            positions.empty() ? &empty : positions.data(),
                            random, result.residual.data(), n1, k1, n2, k2,
                            postprocessing, 1, result.supplemental.data()),
      0);
  return result;
}
}  // namespace

TEST(ProductMonteCarloIdeal, ExhaustiveSmallSupportsMatchIndependentPeeling) {
  // Independent full sweeps, without the implementation's active-line cache.
  // Exhaust every support, including parity cells, and both final-repair modes.
  for (unsigned mask = 0; mask < (1u << 16); ++mask) {
    std::vector<uint8_t> noise(16);
    for (size_t i = 0; i < noise.size(); ++i) {
      noise[i] = (mask >> i) & 1;
    }
    auto expected = noise;
    while (true) {
      const auto before = expected;
      for (size_t c = 0; c < 4; ++c) {
        size_t degree = 0;
        for (size_t r = 0; r < 4; ++r) {
          degree += expected[4 * r + c] != 0;
        }
        if (degree <= 1) {
          for (size_t r = 0; r < 4; ++r) {
            expected[4 * r + c] = 0;
          }
        }
      }
      for (size_t r = 0; r < 4; ++r) {
        if (std::count(expected.begin() + 4 * r, expected.begin() + 4 * (r + 1),
                       1) <= 1) {
          std::fill(expected.begin() + 4 * r, expected.begin() + 4 * (r + 1),
                    0);
        }
      }
      if (expected == before) {
        break;
      }
    }
    size_t columns = 0;
    for (size_t c = 0; c < 4; ++c) {
      bool wrong = false;
      for (size_t r = 0; r < 4; ++r) {
        wrong |= expected[4 * r + c] != 0;
      }
      columns += wrong;
    }
    for (bool postprocessing : {false, true}) {
      auto oracle = expected;
      const bool repaired = postprocessing && columns > 0 && columns <= 2;
      if (repaired) {
        std::fill(oracle.begin(), oracle.end(), 0);
      }
      const auto result = Decode(noise, 4, 2, 4, 2, postprocessing);
      ASSERT_EQ(result.residual, oracle) << mask << " pp=" << postprocessing;
      const size_t remaining = std::count(oracle.begin(), oracle.end(), 1);
      ASSERT_EQ(result.metrics[4], remaining) << mask;
      ASSERT_EQ(result.metrics[10], remaining == 0) << mask;
      ASSERT_EQ(result.metrics[13], std::popcount(mask) - remaining) << mask;
      ASSERT_EQ(result.metrics[14], result.metrics[13]) << mask;
      ASSERT_EQ(result.metrics[11], 0u);
      ASSERT_EQ(result.metrics[21], 0u);
      ASSERT_EQ(result.supplemental[0], repaired);
      ASSERT_EQ(result.supplemental[1], 0u);
      ASSERT_EQ(result.supplemental[2], 0u);
    }
  }
}

TEST(ProductMonteCarloIdeal, CountsBytesAndRepairsAllTheirBitsAtRadius) {
  std::vector<uint8_t> noise(8 * 5);
  // Strong radius two: 16 wrong bits in two bytes still correct completely.
  noise[0] = noise[5] = 255;
  for (bool random : {false, true}) {
    const auto result = Decode(noise, 8, 4, 5, 3, false, 64, random);
    EXPECT_EQ(result.metrics[4], 0u);
    EXPECT_EQ(result.metrics[15], 16u);
    EXPECT_EQ(result.metrics[16], 2u);
  }
  // Weak radius one repairs an eight-bit error in a failed strong column.
  noise[10] = 255;
  const auto weak = Decode(noise, 8, 4, 5, 3);
  EXPECT_EQ(weak.metrics[4], 0u);
  EXPECT_EQ(weak.metrics[15], 0u);
  EXPECT_EQ(weak.metrics[17], 24u);
  EXPECT_EQ(weak.metrics[18], 3u);
  // One beyond each radius stalls. Exactly two faulty columns are erased.
  noise[1] = noise[6] = noise[11] = 255;
  const auto stalled = Decode(noise, 8, 4, 5, 3);
  EXPECT_EQ(stalled.residual, noise);
  EXPECT_EQ(stalled.metrics[4], 48u);
  const auto repaired = Decode(noise, 8, 4, 5, 3, true);
  EXPECT_EQ(repaired.metrics[4], 0u);
  EXPECT_EQ(repaired.supplemental[0], 1u);
  // Gates are irrelevant, including for a transmitted random codeword.
  EXPECT_EQ(Decode(noise, 8, 4, 5, 3, true, 64, true, false).metrics,
            repaired.metrics);
}

TEST(ProductMonteCarloIdeal, PassCapPreservesErrorsThatLaterPassesRepair) {
  std::vector<uint8_t> noise(16);
  for (size_t i : {0, 1, 5, 6, 10}) {
    noise[i] = 1;
  }
  const auto capped = Decode(noise, 4, 2, 4, 2, false, 2);
  EXPECT_EQ(capped.metrics[4], 2u);
  EXPECT_EQ(capped.metrics[12], 2u);
  EXPECT_EQ(capped.metrics[21], 1u);
  const auto converged = Decode(noise);
  EXPECT_EQ(converged.metrics[4], 0u);
  EXPECT_EQ(converged.metrics[12], 4u);
  EXPECT_EQ(converged.metrics[21], 0u);
  const auto post = Decode(noise, 4, 2, 4, 2, true, 2);
  EXPECT_EQ(post.metrics[4], 0u);
  EXPECT_EQ(post.metrics[12], 2u);
  EXPECT_EQ(post.metrics[21], 1u);
  EXPECT_EQ(post.supplemental[0], 1u);
}
