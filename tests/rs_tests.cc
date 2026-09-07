#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

#include "gtest/gtest.h"
#if defined(GF256_ENABLE_GFNI512_RADIX8_EXPERIMENT)
#include "lin_chung_han/experiment/gfni512_radix8.h"
#include "reed_solomon/experiment/gfni512_radix8.h"
#endif
#include "reed_solomon/code_parameters.h"
#if defined(__AVX2__)
#include "reed_solomon/error_correction/avx2_internal.h"
#endif
#include "reed_solomon/error_correction/internal.h"
#include "reed_solomon/lch_decoder.h"
#include "reed_solomon/lch_encoder.h"

namespace {

using gf2p8::Element;
using gf2p8::lch::Backend;
using gf2p8::lch::Radix;
using gf2p8::lch::Status;
using gf2p8::rs::LCHDecoder;
using gf2p8::rs::LCHEncoder;
using gf2p8::rs::detail::error_correction::CorrectBatch;
using gf2p8::rs::detail::error_correction::CorrectionStatus;
using gf2p8::rs::detail::error_correction::CorrectOne;

std::vector<Element*> MutablePointers(
    std::vector<std::vector<Element>>& shards) {
  std::vector<Element*> result;
  result.reserve(shards.size());
  for (auto& shard : shards) {
    result.push_back(shard.data());
  }
  return result;
}

std::vector<const Element*> ConstPointers(
    const std::vector<std::vector<Element>>& shards) {
  std::vector<const Element*> result;
  result.reserve(shards.size());
  for (const auto& shard : shards) {
    result.push_back(shard.data());
  }
  return result;
}

std::vector<std::vector<Element>> RandomShards(size_t count,
                                               size_t bytes,
                                               uint32_t seed) {
  std::mt19937 random(seed);
  std::vector<std::vector<Element>> result(count, std::vector<Element>(bytes));
  for (auto& shard : result) {
    std::generate(shard.begin(), shard.end(),
                  [&random] { return static_cast<Element>(random()); });
  }
  return result;
}

std::vector<Backend> AvailableBackends() {
  constexpr Backend candidates[] = {
      Backend::scalar,         Backend::ssse3,          Backend::avx2,
      Backend::gfni128_affine, Backend::gfni256_affine, Backend::gfni512_affine,
  };
  std::vector<Backend> result;
  for (const Backend backend : candidates) {
    if (gf2p8::lch::BackendAvailable(backend)) {
      result.push_back(backend);
    }
  }
  return result;
}

size_t NextPowerOfTwo(size_t value) {
  size_t result = 1;
  while (result < value) {
    result *= 2;
  }
  return result;
}

std::vector<std::vector<Element>> Encode(
    const LCHEncoder& encoder,
    const std::vector<std::vector<Element>>& data,
    size_t bytes,
    Backend backend = Backend::scalar,
    Radix radix = Radix::radix2) {
  std::vector<std::vector<Element>> recovery(encoder.RecoveryCount(),
                                             std::vector<Element>(bytes));
  const auto data_pointers = ConstPointers(data);
  auto recovery_pointers = MutablePointers(recovery);
  std::vector<Element> workspace(encoder.WorkspaceSize(bytes));
  EXPECT_EQ(encoder.Encode(data_pointers, recovery_pointers, bytes, workspace,
                           backend, radix),
            Status::ok);
  return recovery;
}

std::vector<Element> EncodeOne(const LCHEncoder& encoder,
                               std::span<const Element> data) {
  std::vector<std::vector<Element>> data_shards(data.size(),
                                                std::vector<Element>(1));
  for (size_t i = 0; i < data.size(); ++i) {
    data_shards[i][0] = data[i];
  }
  const auto recovery_shards = Encode(encoder, data_shards, 1);
  std::vector<Element> recovery(recovery_shards.size());
  for (size_t i = 0; i < recovery.size(); ++i) {
    recovery[i] = recovery_shards[i][0];
  }
  return recovery;
}

void ExpectCorrects(size_t data_count,
                    size_t recovery_count,
                    std::span<const Element> expected_data,
                    std::span<const std::pair<size_t, Element>> errors) {
  const LCHEncoder encoder(data_count, recovery_count);
  const LCHDecoder decoder(data_count, recovery_count);
  ASSERT_TRUE(encoder.Valid());
  ASSERT_TRUE(decoder.Valid());

  const std::vector<Element> expected_recovery =
      EncodeOne(encoder, expected_data);
  std::vector<Element> data(expected_data.begin(), expected_data.end());
  std::vector<Element> recovery = expected_recovery;
  std::vector<uint8_t> expected_mask(data_count + recovery_count, 0);
  for (const auto [position, magnitude] : errors) {
    ASSERT_LT(position, expected_mask.size());
    ASSERT_NE(magnitude, 0);
    ASSERT_EQ(expected_mask[position], 0);
    expected_mask[position] = 1;
    if (position < data_count) {
      data[position] ^= magnitude;
    } else {
      recovery[position - data_count] ^= magnitude;
    }
  }
  const std::vector<Element> corrupted_recovery = recovery;

  std::vector<uint8_t> actual_mask(expected_mask.size(), 0xa5);
  const auto result = CorrectOne(decoder, data, recovery, actual_mask);
  ASSERT_EQ(result.status, CorrectionStatus::ok);
  EXPECT_EQ(result.error_count, errors.size());
  EXPECT_EQ(actual_mask, expected_mask);
  EXPECT_TRUE(std::equal(data.begin(), data.end(), expected_data.begin(),
                         expected_data.end()));
  EXPECT_EQ(recovery, corrupted_recovery);

  const std::vector<Element> rebuilt_recovery = EncodeOne(encoder, data);
  for (size_t i = 0; i < recovery_count; ++i) {
    if (expected_mask[data_count + i] == 0) {
      EXPECT_EQ(rebuilt_recovery[i], recovery[i]);
    }
  }
}

std::vector<std::vector<Element>> InterpolationOracle(
    const std::vector<std::vector<Element>>& data,
    size_t recovery_count) {
  const size_t k = data.size();
  const bool low_rate = recovery_count >= k;
  const size_t transform_size = NextPowerOfTwo(low_rate ? k : recovery_count);
  const size_t mother_size =
      NextPowerOfTwo(transform_size + (low_rate ? recovery_count : k));
  const size_t systematic_count =
      low_rate ? transform_size : mother_size - transform_size;
  const size_t systematic_offset = low_rate ? 0 : transform_size;
  const size_t recovery_offset = low_rate ? transform_size : 0;
  const size_t bytes = data.empty() ? 0 : data.front().size();
  std::vector<std::vector<Element>> recovery(recovery_count,
                                             std::vector<Element>(bytes));

  for (size_t output = 0; output < recovery_count; ++output) {
    const Element x = static_cast<Element>(recovery_offset + output);
    for (size_t source = 0; source < k; ++source) {
      const Element source_x = static_cast<Element>(systematic_offset + source);
      Element numerator = 1;
      Element denominator = 1;
      for (size_t other = 0; other < systematic_count; ++other) {
        if (other == source) {
          continue;
        }
        const Element other_x = static_cast<Element>(systematic_offset + other);
        numerator = gf2p8::MultiplyCantor(numerator, x ^ other_x);
        denominator = gf2p8::MultiplyCantor(denominator, source_x ^ other_x);
      }
      const Element coefficient = gf2p8::DivCantor(numerator, denominator);
      for (size_t byte = 0; byte < bytes; ++byte) {
        recovery[output][byte] ^=
            gf2p8::MultiplyCantor(coefficient, data[source][byte]);
      }
    }
  }
  return recovery;
}

void Recover(const LCHDecoder& decoder,
             const std::vector<std::vector<Element>>& expected_data,
             const std::vector<std::vector<Element>>& recovery,
             std::span<const size_t> canonical_losses,
             size_t bytes,
             Backend backend,
             Radix radix) {
  auto data = expected_data;
  auto data_pointers = MutablePointers(data);
  auto recovery_pointers = ConstPointers(recovery);
  std::vector<uint8_t> data_present(data.size(), 1);
  std::vector<uint8_t> recovery_present(recovery.size(), 1);
  for (const size_t loss : canonical_losses) {
    if (loss < recovery.size()) {
      recovery_present[loss] = 0;
      recovery_pointers[loss] = nullptr;
    } else {
      const size_t data_index = loss - recovery.size();
      data_present[data_index] = 0;
      std::fill(data[data_index].begin(), data[data_index].end(), 0xa5);
    }
  }
  std::vector<Element> workspace(decoder.WorkspaceSize(bytes));
  ASSERT_EQ(decoder.Decode(data_pointers, data_present, recovery_pointers,
                           recovery_present, bytes, workspace, backend, radix),
            Status::ok);
  EXPECT_EQ(data, expected_data);
}

TEST(LCHCode, ValidatesDimensionsAndNormalizesFamilies) {
  EXPECT_FALSE(LCHEncoder(0, 1).Valid());
  EXPECT_FALSE(LCHDecoder(8, 0).Valid());
  EXPECT_FALSE(LCHEncoder(5, 249).Valid());
  EXPECT_FALSE(LCHDecoder(130, 126).Valid());
  EXPECT_FALSE(LCHEncoder(std::numeric_limits<size_t>::max(), 1).Valid());

  for (const auto [k, r] : {std::pair<size_t, size_t>{1, 1},
                            {5, 3},
                            {5, 5},
                            {5, 6},
                            {37, 100},
                            {128, 128},
                            {192, 64},
                            {255, 1}}) {
    EXPECT_TRUE(LCHEncoder(k, r).Valid()) << k << '/' << r;
    EXPECT_TRUE(LCHDecoder(k, r).Valid()) << k << '/' << r;
  }

  const auto high_shortened = gf2p8::rs::detail::MakeCodeParameters(5, 3);
  EXPECT_EQ(high_shortened.family, gf2p8::rs::detail::CodeFamily::high_rate);
  EXPECT_EQ(high_shortened.transform_size, 4);
  EXPECT_EQ(high_shortened.mother_size, 16);
  EXPECT_EQ(gf2p8::rs::detail::MakeCodeParameters(5, 5).family,
            gf2p8::rs::detail::CodeFamily::low_rate);
  const auto low_shortened = gf2p8::rs::detail::MakeCodeParameters(5, 6);
  EXPECT_EQ(low_shortened.family, gf2p8::rs::detail::CodeFamily::low_rate);
  EXPECT_EQ(low_shortened.transform_size, 8);
  EXPECT_EQ(low_shortened.mother_size, 16);
  EXPECT_EQ(gf2p8::rs::detail::MakeCodeParameters(192, 64).family,
            gf2p8::rs::detail::CodeFamily::high_rate);
  const auto high_padded = gf2p8::rs::detail::MakeCodeParameters(129, 64);
  EXPECT_EQ(high_padded.family, gf2p8::rs::detail::CodeFamily::high_rate);
  EXPECT_EQ(high_padded.transform_size, 64);
  EXPECT_EQ(high_padded.mother_size, 256);
}

TEST(LCHCode, MatchesIndependentInterpolationOracle) {
  for (const auto [k, r] : {std::pair<size_t, size_t>{3, 3},
                            {3, 5},
                            {5, 3},
                            {6, 5},
                            {5, 6},
                            {9, 17},
                            {37, 11}}) {
    constexpr size_t kBytes = 7;
    const auto data = RandomShards(k, kBytes, 1000 + k * 257 + r);
    const LCHEncoder encoder(k, r);
    ASSERT_TRUE(encoder.Valid());
    EXPECT_EQ(Encode(encoder, data, kBytes), InterpolationOracle(data, r))
        << k << '/' << r;
  }
}

TEST(LCHCode, MatchesPinnedLeopardCodewordByteForByte) {
  constexpr size_t kDataCount = 5;
  constexpr size_t kRecoveryCount = 3;
  constexpr size_t kBytes = 64;
  // Generated with catid/leopard@6e5725e via leo_encode().
  constexpr std::array<std::array<Element, kBytes>, kRecoveryCount>
      kUpstreamRecovery{{
          {0x8d, 0x0f, 0x26, 0xcd, 0xf4, 0xb2, 0x2e, 0x17, 0xdc, 0x23, 0x17,
           0x93, 0xbf, 0xcb, 0xc4, 0x86, 0x04, 0x21, 0xc2, 0x87, 0x52, 0x2d,
           0x11, 0xd7, 0xa7, 0x8d, 0x9c, 0xb6, 0xc5, 0xc7, 0x37, 0x0f, 0x2a,
           0xc5, 0x88, 0xb5, 0x70, 0x12, 0xd1, 0xac, 0x52, 0xf3, 0xb9, 0xcc,
           0xc9, 0x34, 0x3b, 0x21, 0xce, 0x8f, 0xba, 0xd8, 0xf2, 0xd2, 0xaa,
           0x59, 0x6d, 0x23, 0xc3, 0xc0, 0x3a, 0x38, 0x90, 0xc5},
          {0x3c, 0xde, 0x26, 0xa0, 0x93, 0xcf, 0xe0, 0xbc, 0xcc, 0x0f, 0xd1,
           0xf5, 0xaf, 0x4f, 0xaa, 0x3e, 0xd3, 0x21, 0xaa, 0x7b, 0x7e, 0xee,
           0xb2, 0xce, 0xa3, 0x61, 0xff, 0xaa, 0x49, 0xa4, 0xad, 0xd1, 0x2c,
           0xad, 0x71, 0xaf, 0x48, 0xbc, 0xc0, 0xa1, 0x82, 0x8e, 0xa0, 0x4c,
           0xa2, 0xa3, 0x8a, 0x2e, 0xa0, 0x76, 0xa5, 0xc2, 0x0d, 0xce, 0xaf,
           0x80, 0xd0, 0x10, 0x46, 0xa7, 0xa5, 0x84, 0xbd, 0xa2},
          {0x2d, 0x05, 0xb7, 0xe5, 0x28, 0x29, 0x04, 0x9d, 0x43, 0xcd, 0x9c,
           0x5a, 0x09, 0xd2, 0xaa, 0x2d, 0x00, 0xbd, 0xe4, 0x79, 0xf6, 0x0c,
           0x9b, 0x43, 0xbe, 0x78, 0x5b, 0x0d, 0xd9, 0xa2, 0xff, 0x00, 0xb8,
           0xee, 0x78, 0x6f, 0x92, 0x93, 0x45, 0xbe, 0xca, 0xf0, 0x0c, 0xdd,
           0xa9, 0xf7, 0xeb, 0xb8, 0xeb, 0x72, 0x6e, 0xfe, 0x4c, 0x4d, 0xb8,
           0xca, 0x55, 0xe8, 0xdc, 0xad, 0xfc, 0xe3, 0x6a, 0xeb},
      }};

  std::vector<std::vector<Element>> data(kDataCount,
                                         std::vector<Element>(kBytes));
  for (size_t i = 0; i < kDataCount; ++i) {
    for (size_t j = 0; j < kBytes; ++j) {
      data[i][j] = static_cast<Element>((i * 53 + j * 17 + 7) & 255);
    }
  }
  const LCHEncoder encoder(kDataCount, kRecoveryCount);
  const LCHDecoder decoder(kDataCount, kRecoveryCount);
  std::vector<std::vector<Element>> pinned_recovery;
  pinned_recovery.reserve(kRecoveryCount);
  for (const auto& shard : kUpstreamRecovery) {
    pinned_recovery.emplace_back(shard.begin(), shard.end());
  }
  for (const Backend backend : AvailableBackends()) {
    for (const Radix radix : {Radix::radix2, Radix::radix4}) {
      const auto recovery = Encode(encoder, data, kBytes, backend, radix);
      EXPECT_EQ(recovery, pinned_recovery);
      const std::array<size_t, 2> losses{0, kRecoveryCount + 2};
      Recover(decoder, data, pinned_recovery, losses, kBytes, backend, radix);
      const std::array<size_t, 3> maximum_losses{
          0, kRecoveryCount, kRecoveryCount + kDataCount - 1};
      Recover(decoder, data, pinned_recovery, maximum_losses, kBytes, backend,
              radix);
    }
  }
}

TEST(LCHCode, RecoversAcrossFamiliesBackendsAndTails) {
  struct Parameters {
    size_t data_count;
    size_t recovery_count;
    size_t bytes;
  };
  constexpr Parameters parameters[] = {
      {1, 1, 0},    {1, 3, 17},    {1, 255, 1},    {5, 3, 31},    {6, 5, 63},
      {5, 5, 32},   {5, 6, 33},    {5, 9, 65},     {9, 17, 17},   {37, 11, 33},
      {65, 66, 33}, {127, 128, 1}, {128, 128, 17}, {129, 64, 33}, {192, 64, 17},
      {248, 8, 65}, {255, 1, 33},
  };

  for (const auto [k, r, bytes] : parameters) {
    const LCHEncoder encoder(k, r);
    const LCHDecoder decoder(k, r);
    ASSERT_TRUE(encoder.Valid()) << k << '/' << r;
    ASSERT_TRUE(decoder.Valid()) << k << '/' << r;
    const auto data = RandomShards(k, bytes, 2000 + k * 257 + r + bytes);
    const auto recovery = Encode(encoder, data, bytes);

    std::vector<size_t> losses(k + r);
    std::iota(losses.begin(), losses.end(), 0);
    std::mt19937 random(static_cast<uint32_t>(3000 + k * 257 + r));
    std::shuffle(losses.begin(), losses.end(), random);
    const auto missing_data =
        std::find(losses.begin(), losses.end(), r + k / 2);
    std::iter_swap(losses.begin(), missing_data);
    losses.resize(r);

    for (const Backend backend : AvailableBackends()) {
      for (const Radix radix : {Radix::radix2, Radix::radix4}) {
        const auto backend_recovery =
            Encode(encoder, data, bytes, backend, radix);
        EXPECT_EQ(backend_recovery, recovery);
        Recover(decoder, data, backend_recovery, losses, bytes, backend, radix);
      }
    }
  }
}

TEST(LCHCode, ExhaustivelyRecoversSmallShortenedCodes) {
  for (const auto [k, r] : {std::pair<size_t, size_t>{2, 1},
                            {3, 2},
                            {4, 3},
                            {5, 3},
                            {6, 5},
                            {3, 3},
                            {3, 5},
                            {5, 6}}) {
    constexpr size_t kBytes = 3;
    const LCHDecoder decoder(k, r);
    const auto data = RandomShards(k, kBytes, 3500 + k * 257 + r);
    const auto recovery = InterpolationOracle(data, r);
    const size_t logical_count = k + r;
    for (uint32_t mask = 1; mask < (uint32_t{1} << logical_count); ++mask) {
      if (std::popcount(mask) > r || (mask >> r) == 0) {
        continue;
      }
      std::vector<size_t> losses;
      for (size_t i = 0; i < logical_count; ++i) {
        if ((mask & (uint32_t{1} << i)) != 0) {
          losses.push_back(i);
        }
      }
      Recover(decoder, data, recovery, losses, kBytes, Backend::scalar,
              Radix::radix2);
    }
  }
}

TEST(LCHCode, RecoversEverySupportedDimension) {
  constexpr size_t kBytes = 1;
  size_t valid_dimension_count = 0;
  for (size_t k = 1; k < gf2p8::lch::Context::kFieldSize; ++k) {
    for (size_t r = 1; r < gf2p8::lch::Context::kFieldSize; ++r) {
      const LCHEncoder encoder(k, r);
      const LCHDecoder decoder(k, r);
      if (!encoder.Valid()) {
        EXPECT_FALSE(decoder.Valid()) << k << '/' << r;
        continue;
      }
      ASSERT_TRUE(decoder.Valid()) << k << '/' << r;
      ++valid_dimension_count;

      std::vector<std::vector<Element>> data(k, std::vector<Element>(kBytes));
      for (size_t i = 0; i < k; ++i) {
        data[i][0] = static_cast<Element>((k * 17 + r * 29 + i * 53) & 255);
      }
      const auto recovery = Encode(encoder, data, kBytes);

      std::vector<size_t> data_heavy_losses;
      const size_t data_losses = std::min(k, r);
      for (size_t i = 0; i < data_losses; ++i) {
        data_heavy_losses.push_back(r + i);
      }
      for (size_t i = 0; data_heavy_losses.size() < r; ++i) {
        data_heavy_losses.push_back(i);
      }
      Recover(decoder, data, recovery, data_heavy_losses, kBytes,
              Backend::scalar, Radix::radix2);

      std::vector<size_t> shuffled_losses(k + r);
      std::iota(shuffled_losses.begin(), shuffled_losses.end(), 0);
      std::mt19937 random(static_cast<uint32_t>(k * 65537 + r));
      std::shuffle(shuffled_losses.begin(), shuffled_losses.end(), random);
      const auto required_data =
          std::find(shuffled_losses.begin(), shuffled_losses.end(), r + k - 1);
      std::iter_swap(shuffled_losses.begin(), required_data);
      shuffled_losses.resize(r);
      Recover(decoder, data, recovery, shuffled_losses, kBytes, Backend::scalar,
              Radix::radix2);
    }
  }
  EXPECT_EQ(valid_dimension_count, 27306);
}

TEST(LCHCode, HighRateZeroByteDecodeUsesLogicalPadding) {
  constexpr size_t k = 5;
  constexpr size_t r = 3;
  const LCHDecoder decoder(k, r);
  EXPECT_EQ(decoder.WorkspaceSize(0), 256);

  std::array<Element*, k> data{};
  std::array<const Element*, r> recovery{};
  std::array<uint8_t, k> data_present{};
  std::array<uint8_t, r> recovery_present{};
  data_present.fill(1);
  recovery_present.fill(1);
  data_present[0] = 0;
  data_present[k - 1] = 0;
  recovery_present[0] = 0;
  std::array<Element, 256> workspace{};

  EXPECT_EQ(decoder.Decode(data, data_present, recovery, recovery_present, 0,
                           std::span<Element>(workspace).first(255)),
            Status::invalid_argument);
  for (const Backend backend : AvailableBackends()) {
    for (const Radix radix : {Radix::radix2, Radix::radix4}) {
      EXPECT_EQ(decoder.Decode(data, data_present, recovery, recovery_present,
                               0, workspace, backend, radix),
                Status::ok);
    }
  }
}

TEST(LCHCode, EnforcesWorkspaceAndSeparateSpanContracts) {
  constexpr size_t kBytes = 17;
  const LCHEncoder folded_encoder(5, 3);
  const LCHEncoder low_encoder(5, 6);
  const LCHDecoder decoder(5, 6);
  EXPECT_EQ(folded_encoder.WorkspaceSize(kBytes), 5 * kBytes);
  EXPECT_EQ(low_encoder.WorkspaceSize(kBytes), 2 * kBytes);
  EXPECT_EQ(decoder.WorkspaceSize(kBytes), 256 + 16 * kBytes);
  EXPECT_EQ(low_encoder.WorkspaceSize(std::numeric_limits<size_t>::max()),
            std::numeric_limits<size_t>::max());
  EXPECT_EQ(decoder.WorkspaceSize(std::numeric_limits<size_t>::max()),
            std::numeric_limits<size_t>::max());

  const auto data = RandomShards(5, kBytes, 4000);
  std::vector<std::vector<Element>> recovery(6, std::vector<Element>(kBytes));
  const auto data_input = ConstPointers(data);
  auto recovery_output = MutablePointers(recovery);
  std::vector<Element> encode_workspace(low_encoder.WorkspaceSize(kBytes));
  ASSERT_EQ(
      low_encoder.Encode(data_input, recovery_output, kBytes, encode_workspace),
      Status::ok);
  EXPECT_EQ(low_encoder.Encode(data_input, recovery_output, kBytes,
                               std::span<Element>(encode_workspace)
                                   .first(encode_workspace.size() - 1)),
            Status::invalid_argument);

  auto mutable_data = data;
  auto data_output = MutablePointers(mutable_data);
  auto recovery_input = ConstPointers(recovery);
  std::vector<uint8_t> data_present(5, 1);
  std::vector<uint8_t> recovery_present(6, 1);
  data_present[2] = 0;
  std::vector<Element> decode_workspace(decoder.WorkspaceSize(kBytes));
  EXPECT_EQ(decoder.Decode(data_output, data_present, recovery_input,
                           recovery_present, kBytes,
                           std::span<Element>(decode_workspace)
                               .first(decode_workspace.size() - 1)),
            Status::invalid_argument);
  EXPECT_EQ(decoder.Decode(data_output, data_present,
                           std::span<const Element* const>(recovery_input)
                               .first(recovery_input.size() - 1),
                           recovery_present, kBytes, decode_workspace),
            Status::invalid_argument);
  EXPECT_EQ(decoder.Decode(
                data_output, std::span<const uint8_t>(data_present).first(4),
                recovery_input, recovery_present, kBytes, decode_workspace),
            Status::invalid_argument);
  EXPECT_EQ(decoder.Decode(data_output, data_present, recovery_input,
                           std::span<const uint8_t>(recovery_present).first(5),
                           kBytes, decode_workspace),
            Status::invalid_argument);

  recovery_present[0] = 0;
  recovery_input[0] = nullptr;
  ASSERT_EQ(decoder.Decode(data_output, data_present, recovery_input,
                           recovery_present, kBytes, decode_workspace),
            Status::ok);
  EXPECT_EQ(mutable_data, data);
}

TEST(LCHCode, MovedFromObjectsAreSafelyInvalid) {
  LCHEncoder encoder(5, 6);
  LCHEncoder moved_encoder(std::move(encoder));
  EXPECT_TRUE(moved_encoder.Valid());
  EXPECT_FALSE(encoder.Valid());
  EXPECT_EQ(encoder.DataCount(), 0);
  EXPECT_EQ(encoder.RecoveryCount(), 0);
  EXPECT_EQ(encoder.WorkspaceSize(17), 0);

  LCHDecoder decoder(5, 6);
  LCHDecoder moved_decoder(std::move(decoder));
  EXPECT_TRUE(moved_decoder.Valid());
  EXPECT_FALSE(decoder.Valid());
  EXPECT_EQ(decoder.DataCount(), 0);
  EXPECT_EQ(decoder.RecoveryCount(), 0);
  EXPECT_EQ(decoder.WorkspaceSize(17), 0);
}

TEST(LCHCode, ReportsInsufficientSymbolsAndRejectsNullRequiredPointers) {
  constexpr size_t kBytes = 17;
  const LCHEncoder encoder(8, 4);
  const LCHDecoder decoder(8, 4);
  const auto expected_data = RandomShards(8, kBytes, 5000);
  const auto recovery = Encode(encoder, expected_data, kBytes);
  auto data = expected_data;
  auto data_pointers = MutablePointers(data);
  auto recovery_pointers = ConstPointers(recovery);
  std::vector<uint8_t> data_present(8, 1);
  std::vector<uint8_t> recovery_present(4, 0);
  data_present[0] = 0;
  std::vector<Element> workspace(decoder.WorkspaceSize(kBytes));
  EXPECT_EQ(decoder.Decode(data_pointers, data_present, recovery_pointers,
                           recovery_present, kBytes, workspace),
            Status::insufficient_recovery_symbols);

  std::fill(recovery_present.begin(), recovery_present.end(), 1);
  data_pointers[0] = nullptr;
  EXPECT_EQ(decoder.Decode(data_pointers, data_present, recovery_pointers,
                           recovery_present, kBytes, workspace),
            Status::invalid_argument);
  data_pointers[0] = data[0].data();
  data_pointers[1] = nullptr;
  EXPECT_EQ(decoder.Decode(data_pointers, data_present, recovery_pointers,
                           recovery_present, kBytes, workspace),
            Status::invalid_argument);
  data_pointers[1] = data[1].data();
  recovery_pointers[0] = nullptr;
  EXPECT_EQ(decoder.Decode(data_pointers, data_present, recovery_pointers,
                           recovery_present, kBytes, workspace),
            Status::invalid_argument);
}

TEST(LCHCode, ZeroByteOperationsAllowNullShardPointers) {
  const LCHEncoder encoder(5, 6);
  const LCHDecoder decoder(5, 6);
  std::array<const Element*, 5> data_input{};
  std::array<Element*, 6> recovery_output{};
  std::vector<Element> encode_workspace(encoder.WorkspaceSize(0));
  EXPECT_EQ(encoder.Encode(data_input, recovery_output, 0, encode_workspace),
            Status::ok);

  std::array<Element*, 5> data_output{};
  std::array<const Element*, 6> recovery_input{};
  std::array<uint8_t, 5> data_present{};
  std::array<uint8_t, 6> recovery_present{};
  data_present.fill(1);
  recovery_present.fill(1);
  std::vector<Element> decode_workspace(decoder.WorkspaceSize(0));
  EXPECT_EQ(decoder.Decode(data_output, data_present, recovery_input,
                           recovery_present, 0, decode_workspace),
            Status::ok);
}

TEST(LCHErrorCorrection, ValidatesArgumentsAndSupportedDimensions) {
  std::array<Element, 6> data = {1, 2, 3, 4, 5, 6};
  std::array<Element, 3> recovery = {7, 8, 9};
  std::array<uint8_t, 9> mask{};
  mask.fill(0xa5);

  const LCHDecoder invalid_decoder(0, 1);
  auto result =
      CorrectOne(invalid_decoder, std::span(data).first(0),
                 std::span(recovery).first(1), std::span(mask).first(1));
  EXPECT_EQ(result.status, CorrectionStatus::invalid_argument);
  EXPECT_EQ(mask[0], 0);

  const LCHDecoder supported_decoder(6, 2);
  const auto original_data = data;
  mask.fill(0xa5);
  result = CorrectOne(supported_decoder, std::span(data).first(5),
                      std::span(recovery).first(2), std::span(mask).first(8));
  EXPECT_EQ(result.status, CorrectionStatus::invalid_argument);
  EXPECT_EQ(data, original_data);
  EXPECT_TRUE(std::all_of(mask.begin(), mask.begin() + 8,
                          [](uint8_t value) { return value == 0; }));

  mask.fill(0xa5);
  result = CorrectOne(supported_decoder, data, std::span(recovery).first(1),
                      std::span(mask).first(8));
  EXPECT_EQ(result.status, CorrectionStatus::invalid_argument);
  EXPECT_TRUE(std::all_of(mask.begin(), mask.begin() + 8,
                          [](uint8_t value) { return value == 0; }));

  mask.fill(0xa5);
  result = CorrectOne(supported_decoder, data, std::span(recovery).first(2),
                      std::span(mask).first(7));
  EXPECT_EQ(result.status, CorrectionStatus::invalid_argument);
  EXPECT_TRUE(std::all_of(mask.begin(), mask.begin() + 7,
                          [](uint8_t value) { return value == 0; }));

  const LCHDecoder unsupported_decoder(5, 3);
  mask.fill(0xa5);
  result = CorrectOne(unsupported_decoder, std::span(data).first(5), recovery,
                      std::span(mask).first(8));
  EXPECT_EQ(result.status, CorrectionStatus::unsupported_dimensions);
  EXPECT_EQ(data, original_data);
  EXPECT_TRUE(std::all_of(mask.begin(), mask.begin() + 8,
                          [](uint8_t value) { return value == 0; }));

  const LCHDecoder shortened_decoder(4, 2);
  mask.fill(0xa5);
  result = CorrectOne(shortened_decoder, std::span(data).first(4),
                      std::span(recovery).first(2), std::span(mask).first(6));
  EXPECT_EQ(result.status, CorrectionStatus::unsupported_dimensions);
  EXPECT_EQ(data, original_data);
  EXPECT_TRUE(std::all_of(mask.begin(), mask.begin() + 6,
                          [](uint8_t value) { return value == 0; }));

  const LCHDecoder overlap_decoder(2, 2);
  std::array<Element, 4> overlapping_storage = {1, 2, 3, 4};
  const auto original_storage = overlapping_storage;
  std::array<Element, 2> separate_recovery = {5, 6};
  result = CorrectOne(overlap_decoder, std::span(overlapping_storage).first(2),
                      separate_recovery, overlapping_storage);
  EXPECT_EQ(result.status, CorrectionStatus::invalid_argument);
  EXPECT_EQ(overlapping_storage, original_storage);

  std::array<Element, 2> separate_data = {7, 8};
  std::array<Element, 4> overlapping_recovery = {9, 10, 11, 12};
  const auto original_overlapping_recovery = overlapping_recovery;
  result = CorrectOne(overlap_decoder, separate_data,
                      std::span<const Element>(overlapping_recovery).first(2),
                      overlapping_recovery);
  EXPECT_EQ(result.status, CorrectionStatus::invalid_argument);
  EXPECT_EQ(overlapping_recovery, original_overlapping_recovery);

  std::array<uint8_t, 4> separate_mask{};
  separate_mask.fill(0xa5);
  result = CorrectOne(overlap_decoder, std::span(overlapping_storage).first(2),
                      std::span<const Element>(overlapping_storage).first(2),
                      separate_mask);
  EXPECT_EQ(result.status, CorrectionStatus::invalid_argument);
  EXPECT_EQ(overlapping_storage, original_storage);
  EXPECT_TRUE(std::all_of(separate_mask.begin(), separate_mask.end(),
                          [](uint8_t value) { return value == 0; }));
}

TEST(LCHErrorCorrection, AcceptsUncorruptedFullCodes) {
  for (const auto [data_count, recovery_count] :
       {std::pair<size_t, size_t>{1, 1}, {2, 2}, {6, 2}, {12, 4}, {8, 8}}) {
    const LCHEncoder encoder(data_count, recovery_count);
    const LCHDecoder decoder(data_count, recovery_count);
    std::vector<Element> data(data_count);
    for (size_t i = 0; i < data_count; ++i) {
      data[i] = static_cast<Element>(17 * i + data_count);
    }
    const std::vector<Element> recovery = EncodeOne(encoder, data);
    const std::vector<Element> expected_data = data;
    std::vector<uint8_t> mask(data_count + recovery_count, 0xa5);
    const auto result = CorrectOne(decoder, data, recovery, mask);
    EXPECT_EQ(result.status, CorrectionStatus::ok)
        << data_count << '/' << recovery_count;
    EXPECT_EQ(result.error_count, 0) << data_count << '/' << recovery_count;
    EXPECT_EQ(data, expected_data) << data_count << '/' << recovery_count;
    EXPECT_TRUE(std::all_of(mask.begin(), mask.end(),
                            [](uint8_t value) { return value == 0; }))
        << data_count << '/' << recovery_count;
  }
}

TEST(LCHErrorCorrection, RadiusZeroRejectsNonCodewords) {
  const LCHEncoder encoder(3, 1);
  const LCHDecoder decoder(3, 1);
  const std::array<Element, 3> expected_data = {0x12, 0x34, 0x56};
  std::vector<Element> data(expected_data.begin(), expected_data.end());
  std::vector<Element> recovery = EncodeOne(encoder, data);
  recovery[0] ^= 0x7b;
  const auto corrupted_recovery = recovery;
  std::array<uint8_t, 4> mask{};
  mask.fill(0xa5);

  const auto result = CorrectOne(decoder, data, recovery, mask);
  EXPECT_EQ(result.status, CorrectionStatus::uncorrectable);
  EXPECT_EQ(result.error_count, 0);
  EXPECT_TRUE(std::equal(data.begin(), data.end(), expected_data.begin(),
                         expected_data.end()));
  EXPECT_EQ(recovery, corrupted_recovery);
  EXPECT_TRUE(std::all_of(mask.begin(), mask.end(),
                          [](uint8_t value) { return value == 0; }));
}

TEST(LCHErrorCorrection, ExhaustivelyCorrectsOneErrorForTwoPlusTwo) {
  constexpr std::array<Element, 2> data = {0x31, 0xa7};
  for (size_t position = 0; position < 4; ++position) {
    for (unsigned magnitude = 1; magnitude < 256; ++magnitude) {
      const std::array<std::pair<size_t, Element>, 1> errors = {
          std::pair<size_t, Element>{position,
                                     static_cast<Element>(magnitude)}};
      ExpectCorrects(2, 2, data, errors);
    }
  }
}

TEST(LCHErrorCorrection, ExhaustsLocationSubsetsThroughRadiusForFourPlusFour) {
  constexpr std::array<Element, 4> data = {0x09, 0x53, 0xa1, 0xfe};
  for (size_t first = 0; first < 8; ++first) {
    const std::array<std::pair<size_t, Element>, 1> one_error = {
        std::pair<size_t, Element>{
            first, static_cast<Element>((37 * first + 11) % 255 + 1)}};
    ExpectCorrects(4, 4, data, one_error);

    for (size_t second = first + 1; second < 8; ++second) {
      const std::array<std::pair<size_t, Element>, 2> two_errors = {
          std::pair<size_t, Element>{
              first, static_cast<Element>((37 * first + 11) % 255 + 1)},
          std::pair<size_t, Element>{
              second, static_cast<Element>((53 * second + 19) % 255 + 1)},
      };
      ExpectCorrects(4, 4, data, two_errors);
    }
  }
}

TEST(LCHErrorCorrection, CorrectsDeterministicLargerFullCodes) {
  constexpr std::array<std::pair<size_t, size_t>, 7> dimensions = {
      std::pair<size_t, size_t>{6, 2},
      {12, 4},
      {24, 8},
      {48, 16},
      {96, 32},
      {192, 64},
      {128, 128},
  };
  std::mt19937 random(0x5eed1234);

  for (const auto [data_count, recovery_count] : dimensions) {
    std::vector<Element> data(data_count);
    std::generate(data.begin(), data.end(),
                  [&random] { return static_cast<Element>(random()); });
    const size_t radius = recovery_count / 2;

    for (size_t trial = 0; trial < 6; ++trial) {
      const size_t error_count =
          trial < 3 ? radius : 1 + static_cast<size_t>(random()) % radius;
      std::vector<size_t> positions;
      if (trial % 3 == 0) {
        positions.resize(data_count);
        std::iota(positions.begin(), positions.end(), size_t{0});
      } else if (trial % 3 == 1) {
        positions.resize(recovery_count);
        std::iota(positions.begin(), positions.end(), data_count);
      } else {
        std::vector<size_t> data_positions(data_count);
        std::vector<size_t> recovery_positions(recovery_count);
        std::iota(data_positions.begin(), data_positions.end(), size_t{0});
        std::iota(recovery_positions.begin(), recovery_positions.end(),
                  data_count);
        std::shuffle(data_positions.begin(), data_positions.end(), random);
        std::shuffle(recovery_positions.begin(), recovery_positions.end(),
                     random);
        positions.push_back(data_positions.front());
        if (error_count > 1) {
          positions.push_back(recovery_positions.front());
        }
        data_positions.erase(data_positions.begin());
        if (error_count > 1) {
          recovery_positions.erase(recovery_positions.begin());
        }
        positions.insert(positions.end(), data_positions.begin(),
                         data_positions.end());
        positions.insert(positions.end(), recovery_positions.begin(),
                         recovery_positions.end());
        std::shuffle(positions.begin() + std::min(error_count, size_t{2}),
                     positions.end(), random);
      }
      if (trial % 3 != 2) {
        std::shuffle(positions.begin(), positions.end(), random);
      }

      std::vector<std::pair<size_t, Element>> errors;
      errors.reserve(error_count);
      for (size_t i = 0; i < error_count; ++i) {
        Element magnitude = 0;
        while (magnitude == 0) {
          magnitude = static_cast<Element>(random());
        }
        errors.emplace_back(positions[i], magnitude);
      }
      ExpectCorrects(data_count, recovery_count, data, errors);
    }
  }
}

TEST(LCHErrorCorrection, CorrectsRandomMultiErrorLocationsAndMagnitudes) {
  constexpr std::array<std::pair<size_t, size_t>, 5> dimensions = {
      std::pair<size_t, size_t>{4, 4},
      {12, 4},
      {224, 32},
      {192, 64},
      {128, 128},
  };
  std::mt19937 random(0x7a6b5c4dU);

  for (const auto [data_count, recovery_count] : dimensions) {
    for (size_t trial = 0; trial < 24; ++trial) {
      SCOPED_TRACE(::testing::Message()
                   << "K=" << data_count << " R=" << recovery_count
                   << " trial=" << trial);
      std::vector<Element> data(data_count);
      for (Element& value : data) {
        value = static_cast<Element>(random());
      }

      const size_t radius = recovery_count / 2;
      const size_t error_count = 2 + random() % (radius - 1);
      const size_t first_position = random() % data_count;
      std::vector<size_t> remaining_positions;
      remaining_positions.reserve(data_count + recovery_count - 1);
      for (size_t position = 0; position < data_count + recovery_count;
           ++position) {
        if (position != first_position) {
          remaining_positions.push_back(position);
        }
      }
      std::shuffle(remaining_positions.begin(), remaining_positions.end(),
                   random);

      std::vector<std::pair<size_t, Element>> errors;
      errors.reserve(error_count);
      errors.emplace_back(first_position,
                          static_cast<Element>(1 + random() % 255));
      for (size_t i = 1; i < error_count; ++i) {
        errors.emplace_back(remaining_positions[i - 1],
                            static_cast<Element>(1 + random() % 255));
      }
      ExpectCorrects(data_count, recovery_count, data, errors);
    }
  }
}

#if defined(__AVX2__)
TEST(LCHErrorCorrection, AVX2VariableProductsMatchCantorField) {
  namespace avx2 = gf2p8::rs::detail::error_correction::avx2;
  alignas(32) std::array<Element, 32> first{};
  alignas(32) std::array<Element, 32> second{};
  alignas(32) std::array<Element, 32> products{};
  const auto& tables = gf2p8::Tables();

  for (size_t first_value = 0; first_value < 256; ++first_value) {
    first.fill(static_cast<Element>(first_value));
    const __m256i first_vector =
        _mm256_load_si256(reinterpret_cast<const __m256i*>(first.data()));
    for (size_t second_base = 0; second_base < 256; second_base += 32) {
      for (size_t lane = 0; lane < second.size(); ++lane) {
        second[lane] = static_cast<Element>(second_base + lane);
      }
      const __m256i product = avx2::MultiplyVariable(
          first_vector,
          _mm256_load_si256(reinterpret_cast<const __m256i*>(second.data())),
          tables);
      _mm256_store_si256(reinterpret_cast<__m256i*>(products.data()), product);
      for (size_t lane = 0; lane < products.size(); ++lane) {
        ASSERT_EQ(products[lane],
                  gf2p8::MultiplyCantor(first[lane], second[lane]))
            << "first=" << first_value
            << " second=" << static_cast<unsigned>(second[lane]);
      }
    }
  }
}
#endif

TEST(LCHErrorCorrection, BatchCorrectsDivergentIndependentCodewords) {
  constexpr size_t kBytes = 65;
  constexpr std::array<std::pair<size_t, size_t>, 5> dimensions = {
      std::pair<size_t, size_t>{6, 2},
      {12, 4},
      {224, 32},
      {192, 64},
      {128, 128},
  };

  for (const auto [data_count, recovery_count] : dimensions) {
    SCOPED_TRACE(::testing::Message()
                 << "K=" << data_count << " R=" << recovery_count);
    const LCHEncoder encoder(data_count, recovery_count);
    const LCHDecoder decoder(data_count, recovery_count);
    const auto expected_data = RandomShards(
        data_count, kBytes,
        static_cast<uint32_t>(0xb47c0000U + data_count + recovery_count));
    const auto expected_recovery =
        Encode(encoder, expected_data, kBytes, Backend::scalar, Radix::radix2);
    auto data = expected_data;
    auto recovery = expected_recovery;
    std::vector<uint8_t> expected_mask((data_count + recovery_count) * kBytes,
                                       0);
    std::vector<size_t> expected_counts(kBytes, 0);

    const size_t radius = recovery_count / 2;
    const size_t codeword_size = data_count + recovery_count;
    for (size_t byte = 0; byte < kBytes; ++byte) {
      size_t error_count = 0;
      if (byte % 4 == 1) {
        error_count = 1;
      } else if (byte % 4 == 2) {
        error_count = 1;
      } else if (byte % 4 == 3) {
        error_count = radius;
      }
      expected_counts[byte] = error_count;

      std::vector<size_t> positions(codeword_size);
      std::iota(positions.begin(), positions.end(), size_t{0});
      std::mt19937 random(static_cast<uint32_t>(
          0x51d30000U ^ (data_count << 12U) ^ (recovery_count << 4U) ^ byte));
      std::shuffle(positions.begin(), positions.end(), random);
      if (byte % 4 == 1) {
        positions[0] = byte % data_count;
      } else if (byte % 4 == 2) {
        positions[0] = data_count + byte % recovery_count;
      }

      for (size_t i = 0; i < error_count; ++i) {
        const size_t position = positions[i];
        ASSERT_EQ(expected_mask[position * kBytes + byte], 0);
        expected_mask[position * kBytes + byte] = 1;
        const Element magnitude =
            static_cast<Element>(1 + ((73 * byte + 41 * i + 19) % 255));
        if (position < data_count) {
          data[position][byte] ^= magnitude;
        } else {
          recovery[position - data_count][byte] ^= magnitude;
        }
      }
    }
    const auto corrupted_recovery = recovery;

    auto data_pointers = MutablePointers(data);
    const auto recovery_pointers = ConstPointers(recovery);
    std::vector<gf2p8::rs::detail::error_correction::CorrectionResult> results(
        kBytes);
    std::vector<uint8_t> actual_mask(expected_mask.size(), 0xa5);
    ASSERT_EQ(CorrectBatch(decoder, data_pointers, recovery_pointers, kBytes,
                           results, actual_mask),
              CorrectionStatus::ok);
    for (size_t byte = 0; byte < kBytes; ++byte) {
      EXPECT_EQ(results[byte].status, CorrectionStatus::ok) << byte;
      EXPECT_EQ(results[byte].error_count, expected_counts[byte]) << byte;
    }
    EXPECT_EQ(data, expected_data);
    EXPECT_EQ(recovery, corrupted_recovery);
    EXPECT_EQ(actual_mask, expected_mask);

    const auto rebuilt_recovery =
        Encode(encoder, data, kBytes, Backend::scalar, Radix::radix2);
    for (size_t shard = 0; shard < recovery_count; ++shard) {
      for (size_t byte = 0; byte < kBytes; ++byte) {
        if (expected_mask[(data_count + shard) * kBytes + byte] == 0) {
          EXPECT_EQ(rebuilt_recovery[shard][byte], recovery[shard][byte]);
        }
      }
    }
  }
}

TEST(LCHErrorCorrection, BatchValidatesContractWithoutMutation) {
  constexpr size_t kBytes = 4;
  const LCHEncoder encoder(6, 2);
  const LCHDecoder decoder(6, 2);
  const auto expected_data = RandomShards(6, kBytes, 0x93810000U);
  const auto recovery = Encode(encoder, expected_data, kBytes);
  auto data = expected_data;
  auto data_pointers = MutablePointers(data);
  auto recovery_pointers = ConstPointers(recovery);
  std::array<gf2p8::rs::detail::error_correction::CorrectionResult, kBytes>
      results{};
  std::array<uint8_t, 8 * kBytes> masks{};
  masks.fill(0xa5);

  EXPECT_EQ(CorrectBatch(decoder, std::span(data_pointers).first(5),
                         recovery_pointers, kBytes, results, masks),
            CorrectionStatus::invalid_argument);
  EXPECT_EQ(data, expected_data);
  EXPECT_TRUE(std::all_of(masks.begin(), masks.end(),
                          [](uint8_t value) { return value == 0xa5; }));

  Element* saved = data_pointers[1];
  data_pointers[1] = data_pointers[0];
  EXPECT_EQ(CorrectBatch(decoder, data_pointers, recovery_pointers, kBytes,
                         results, masks),
            CorrectionStatus::invalid_argument);
  data_pointers[1] = saved;
  EXPECT_EQ(data, expected_data);
  EXPECT_TRUE(std::all_of(masks.begin(), masks.end(),
                          [](uint8_t value) { return value == 0xa5; }));

  const LCHDecoder unsupported(4, 2);
  EXPECT_EQ(CorrectBatch(unsupported, std::span(data_pointers).first(4),
                         recovery_pointers, kBytes, results,
                         std::span(masks).first(6 * kBytes)),
            CorrectionStatus::unsupported_dimensions);
  EXPECT_EQ(data, expected_data);
}

TEST(LCHErrorCorrection, BatchIsTransactionalPerCodeword) {
  constexpr size_t kBytes = 32;
  const LCHEncoder encoder(4, 4);
  const LCHDecoder decoder(4, 4);
  const auto expected_data = RandomShards(4, kBytes, 0x71a50000U);
  auto data = expected_data;
  auto recovery = Encode(encoder, expected_data, kBytes);
  const auto expected_recovery = recovery;

  data[0][0] ^= 0x5b;
  data[0][1] ^= 0x01;
  data[1][1] ^= 0x02;
  data[2][1] ^= 0x03;
  const auto corrupted_data = data;
  const auto corrupted_recovery = recovery;

  auto data_pointers = MutablePointers(data);
  const auto recovery_pointers = ConstPointers(recovery);
  std::array<gf2p8::rs::detail::error_correction::CorrectionResult, kBytes>
      results{};
  std::array<uint8_t, 8 * kBytes> masks{};
  masks.fill(0xa5);
  ASSERT_EQ(CorrectBatch(decoder, data_pointers, recovery_pointers, kBytes,
                         results, masks),
            CorrectionStatus::ok);

  EXPECT_EQ(results[0].status, CorrectionStatus::ok);
  EXPECT_EQ(results[0].error_count, 1);
  EXPECT_EQ(results[1].status, CorrectionStatus::uncorrectable);
  EXPECT_EQ(results[1].error_count, 0);
  for (size_t byte = 2; byte < kBytes; ++byte) {
    EXPECT_EQ(results[byte].status, CorrectionStatus::ok);
    EXPECT_EQ(results[byte].error_count, 0);
  }
  for (size_t shard = 0; shard < data.size(); ++shard) {
    EXPECT_EQ(data[shard][0], expected_data[shard][0]);
    EXPECT_EQ(data[shard][1], corrupted_data[shard][1]);
    for (size_t byte = 2; byte < kBytes; ++byte) {
      EXPECT_EQ(data[shard][byte], expected_data[shard][byte]);
    }
  }
  EXPECT_EQ(recovery, corrupted_recovery);
  EXPECT_EQ(masks[0 * kBytes + 0], 1);
  for (size_t position = 0; position < 8; ++position) {
    EXPECT_EQ(masks[position * kBytes + 1], 0);
  }
  EXPECT_EQ(expected_recovery, corrupted_recovery);
}

TEST(LCHErrorCorrection, RejectsSelectedOverRadiusPatternWithoutMutation) {
  const LCHEncoder encoder(4, 4);
  const LCHDecoder decoder(4, 4);
  constexpr std::array<Element, 4> expected_data = {0x22, 0x47, 0x91, 0xd3};
  std::vector<Element> data(expected_data.begin(), expected_data.end());
  std::vector<Element> recovery = EncodeOne(encoder, data);
  data[0] ^= 0x01;
  data[1] ^= 0x02;
  data[2] ^= 0x03;
  const std::vector<Element> corrupted_data = data;
  const std::vector<Element> corrupted_recovery = recovery;
  std::array<uint8_t, 8> mask{};
  mask.fill(0xa5);

  const auto result = CorrectOne(decoder, data, recovery, mask);
  EXPECT_EQ(result.status, CorrectionStatus::uncorrectable);
  EXPECT_EQ(result.error_count, 0);
  EXPECT_EQ(data, corrupted_data);
  EXPECT_EQ(recovery, corrupted_recovery);
  EXPECT_TRUE(std::all_of(mask.begin(), mask.end(),
                          [](uint8_t value) { return value == 0; }));
}

#if defined(GF256_ENABLE_GFNI512_RADIX8_EXPERIMENT)
TEST(LCHRadix8Experiment, FoldedEncoderMatchesProduction) {
  const auto* radix8 = gf2p8::lch::detail::experiment::radix8::ResolveKernels();
  if (radix8 == nullptr) {
    GTEST_SKIP() << "GFNI512 radix-8 was not compiled";
  }
  constexpr size_t k = 224;
  constexpr size_t r = 32;
  constexpr size_t bytes = 65;
  const LCHEncoder encoder(k, r);
  const auto data = RandomShards(k, bytes, 6000);
  const auto expected =
      Encode(encoder, data, bytes, Backend::gfni512_affine, Radix::radix4);
  std::vector<std::vector<Element>> actual(r, std::vector<Element>(bytes));
  const auto input = ConstPointers(data);
  auto output = MutablePointers(actual);
  std::vector<Element> workspace(encoder.WorkspaceSize(bytes));
  const auto& base =
      *gf2p8::lch::detail::ResolveKernels(Backend::gfni512_affine, bytes);
  ASSERT_EQ(gf2p8::rs::detail::experiment::radix8::EncodeLCH(
                gf2p8::lch::Context::Shared(), input, output, bytes, r,
                workspace, base, *radix8),
            Status::ok);
  EXPECT_EQ(actual, expected);
}
#endif

}  // namespace
