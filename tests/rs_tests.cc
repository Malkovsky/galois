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
#include "reed_solomon/lch_decoder.h"
#include "reed_solomon/lch_encoder.h"

namespace {

using gf2p8::Element;
using gf2p8::lch::Backend;
using gf2p8::lch::Radix;
using gf2p8::lch::Status;
using gf2p8::rs::LCHDecoder;
using gf2p8::rs::LCHEncoder;

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
