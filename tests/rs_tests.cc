#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

#include "gtest/gtest.h"
#include "reed_solomon/leopard.h"
#include "reed_solomon/xdrs.h"

namespace {

using gf2p8::Element;
using gf2p8::lch::Backend;
using gf2p8::lch::Radix;
using gf2p8::lch::Status;
using gf2p8::rs::Leopard;
using gf2p8::rs::XDRS;
using gf2p8::rs::XDRSRate;

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
  const Backend candidates[] = {
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

TEST(LeopardRS, ValidatesFf8Parameters) {
  EXPECT_FALSE(Leopard(0, 1).Valid());
  EXPECT_FALSE(Leopard(8, 0).Valid());
  EXPECT_FALSE(Leopard(8, 9).Valid());
  EXPECT_FALSE(Leopard(130, 128).Valid());
  EXPECT_FALSE(Leopard(std::numeric_limits<size_t>::max(),
                       std::numeric_limits<size_t>::max())
                   .Valid());
  EXPECT_TRUE(Leopard(128, 128).Valid());
  EXPECT_TRUE(Leopard(129, 64).Valid());
}

TEST(LeopardRS, MatchesPinnedUpstreamCodewordByteForByte) {
  constexpr size_t kOriginalCount = 5;
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

  Leopard code(kOriginalCount, kRecoveryCount);
  ASSERT_TRUE(code.Valid());
  std::vector<std::vector<Element>> original(kOriginalCount,
                                             std::vector<Element>(kBytes));
  for (size_t i = 0; i < kOriginalCount; ++i) {
    for (size_t j = 0; j < kBytes; ++j) {
      original[i][j] = static_cast<Element>((i * 53 + j * 17 + 7) & 255);
    }
  }
  auto original_pointers = ConstPointers(original);
  std::vector<Element> encode_workspace(code.EncodeWorkspaceSize(kBytes));
  std::vector<Element> decode_workspace(code.DecodeWorkspaceSize(kBytes));
  EXPECT_EQ(encode_workspace.size(), (2 * 4 - kRecoveryCount) * kBytes);

  for (const Backend backend : AvailableBackends()) {
    for (const Radix radix : {Radix::radix2, Radix::radix4}) {
      std::vector<std::vector<Element>> recovery(kRecoveryCount,
                                                 std::vector<Element>(kBytes));
      auto recovery_pointers = MutablePointers(recovery);
      ASSERT_EQ(code.Encode(original_pointers, recovery_pointers, kBytes,
                            encode_workspace, backend, radix),
                Status::ok);
      for (size_t i = 0; i < kRecoveryCount; ++i) {
        EXPECT_EQ(recovery[i],
                  std::vector<Element>(kUpstreamRecovery[i].begin(),
                                       kUpstreamRecovery[i].end()));
      }

      std::vector<std::vector<Element>> codeword = recovery;
      codeword.insert(codeword.end(), original.begin(), original.end());
      std::vector<uint8_t> present(codeword.size(), 1);
      present[0] = 0;
      present[kRecoveryCount + 2] = 0;
      std::fill(codeword[kRecoveryCount + 2].begin(),
                codeword[kRecoveryCount + 2].end(), 0);
      auto codeword_pointers = MutablePointers(codeword);
      ASSERT_EQ(code.Decode(codeword_pointers, present, kBytes,
                            decode_workspace, backend, radix),
                Status::ok);
      EXPECT_EQ(codeword[kRecoveryCount + 2], original[2]);
    }
  }
}

TEST(LeopardRS, NoLossOneLossShortenedAndMaximumLoss) {
  struct Parameters {
    size_t original;
    size_t recovery;
    size_t bytes;
  };
  const Parameters parameters[] = {
      {1, 1, 0},  {1, 1, 17}, {8, 1, 17},   {5, 3, 15},    {5, 3, 16},
      {5, 3, 17}, {5, 3, 31}, {5, 3, 32},   {5, 3, 33},    {5, 3, 63},
      {5, 3, 64}, {5, 3, 65}, {37, 11, 33}, {128, 128, 1},
  };

  for (const auto [original_count, recovery_count, bytes] : parameters) {
    Leopard code(original_count, recovery_count);
    ASSERT_TRUE(code.Valid());
    const auto original = RandomShards(original_count, bytes, 100 + bytes);
    auto recovery = RandomShards(recovery_count, bytes, 200 + bytes);
    const auto original_pointers = ConstPointers(original);
    auto recovery_pointers = MutablePointers(recovery);
    std::vector<Element> encode_workspace(code.EncodeWorkspaceSize(bytes));
    ASSERT_EQ(code.Encode(original_pointers, recovery_pointers, bytes,
                          encode_workspace, Backend::scalar),
              Status::ok);

    std::vector<std::vector<Element>> codeword = recovery;
    codeword.insert(codeword.end(), original.begin(), original.end());
    auto codeword_pointers = MutablePointers(codeword);
    std::vector<uint8_t> present(codeword.size(), 1);
    std::vector<Element> workspace(code.DecodeWorkspaceSize(bytes));

    ASSERT_EQ(code.Decode(codeword_pointers, present, bytes, workspace,
                          Backend::scalar),
              Status::ok);

    for (size_t lost = 0; lost < recovery_count; ++lost) {
      present[recovery_count + lost] = 0;
      std::fill(codeword[recovery_count + lost].begin(),
                codeword[recovery_count + lost].end(), 0);
    }
    ASSERT_EQ(code.Decode(codeword_pointers, present, bytes, workspace,
                          Backend::tuned),
              Status::ok);
    for (size_t i = 0; i < original_count; ++i) {
      EXPECT_EQ(codeword[recovery_count + i], original[i]);
    }

    for (uint32_t trial = 0; trial < 3; ++trial) {
      codeword = recovery;
      codeword.insert(codeword.end(), original.begin(), original.end());
      codeword_pointers = MutablePointers(codeword);
      std::fill(present.begin(), present.end(), 1);
      std::vector<size_t> losses(codeword.size());
      std::iota(losses.begin(), losses.end(), 0);
      std::mt19937 random(static_cast<uint32_t>(900 + bytes + trial));
      std::shuffle(losses.begin(), losses.end(), random);
      const auto original_position = std::find(
          losses.begin(), losses.end(), recovery_count + original_count / 2);
      std::iter_swap(losses.begin(), original_position);
      losses.resize(recovery_count);
      for (const size_t loss : losses) {
        present[loss] = 0;
        std::fill(codeword[loss].begin(), codeword[loss].end(), 0);
      }
      ASSERT_EQ(code.Decode(codeword_pointers, present, bytes, workspace,
                            Backend::tuned),
                Status::ok);
      for (size_t i = 0; i < original_count; ++i) {
        EXPECT_EQ(codeword[recovery_count + i], original[i]);
      }
    }
  }
}

TEST(LeopardRS, ReportsInsufficientAndInvalidInputs) {
  Leopard code(8, 4);
  ASSERT_TRUE(code.Valid());
  auto original = RandomShards(8, 17, 300);
  auto recovery = RandomShards(4, 17, 301);
  const auto original_pointers = ConstPointers(original);
  auto recovery_pointers = MutablePointers(recovery);
  std::vector<Element> encode_workspace(code.EncodeWorkspaceSize(17));
  ASSERT_EQ(
      code.Encode(original_pointers, recovery_pointers, 17, encode_workspace),
      Status::ok);
  EXPECT_EQ(code.Encode(original_pointers, recovery_pointers, 17,
                        std::span<Element>(encode_workspace)
                            .first(encode_workspace.size() - 1)),
            Status::invalid_argument);

  std::vector<std::vector<Element>> codeword = recovery;
  codeword.insert(codeword.end(), original.begin(), original.end());
  auto pointers = MutablePointers(codeword);
  std::vector<uint8_t> present(codeword.size(), 1);
  for (size_t i = 0; i < 5; ++i) {
    present[i] = 0;
  }
  std::vector<Element> workspace(code.DecodeWorkspaceSize(17));
  EXPECT_EQ(code.Decode(pointers, present, 17, workspace),
            Status::insufficient_recovery_symbols);
  EXPECT_EQ(
      code.Decode(pointers, present, 17,
                  std::span<Element>(workspace).first(workspace.size() - 1)),
      Status::invalid_argument);
}

void TestXDRSRecovery(size_t data_count, XDRSRate rate, size_t bytes) {
  XDRS code(data_count, rate);
  ASSERT_TRUE(code.Valid());
  const size_t recovery_count = code.RecoveryCount();
  const auto data = RandomShards(data_count, bytes, 400 + data_count);
  auto recovery = RandomShards(recovery_count, bytes, 500 + data_count);
  const auto data_pointers = ConstPointers(data);
  auto recovery_pointers = MutablePointers(recovery);
  std::vector<Element> encode_workspace(code.EncodeWorkspaceSize(bytes));
  ASSERT_EQ(code.Encode(data_pointers, recovery_pointers, bytes,
                        encode_workspace, Backend::scalar),
            Status::ok);

  const size_t data_offset = rate == XDRSRate::low ? 0 : recovery_count;
  const auto decode_losses = [&](std::span<const size_t> losses) {
    std::vector<std::vector<Element>> codeword;
    if (rate == XDRSRate::low) {
      codeword = data;
      codeword.insert(codeword.end(), recovery.begin(), recovery.end());
    } else {
      codeword = recovery;
      codeword.insert(codeword.end(), data.begin(), data.end());
    }
    auto pointers = MutablePointers(codeword);
    std::vector<uint8_t> present(256, 1);
    for (const size_t loss : losses) {
      present[loss] = 0;
      std::fill(codeword[loss].begin(), codeword[loss].end(), 0);
    }
    std::vector<Element> workspace(code.DecodeWorkspaceSize(bytes));
    ASSERT_EQ(code.Decode(pointers, present, bytes, workspace, Backend::tuned),
              Status::ok);
    for (size_t i = 0; i < data_count; ++i) {
      EXPECT_EQ(codeword[data_offset + i], data[i]);
    }
  };

  decode_losses(std::span<const size_t>{});
  const std::array<size_t, 1> one_loss{data_offset};
  decode_losses(one_loss);

  for (uint32_t trial = 0; trial < 3; ++trial) {
    std::vector<size_t> maximum_losses(256);
    std::iota(maximum_losses.begin(), maximum_losses.end(), 0);
    std::mt19937 random(static_cast<uint32_t>(700 + data_count + 31 * trial));
    std::shuffle(maximum_losses.begin(), maximum_losses.end(), random);
    const auto data_position =
        std::find(maximum_losses.begin(), maximum_losses.end(), data_offset);
    std::iter_swap(maximum_losses.begin(), data_position);
    maximum_losses.resize(recovery_count);
    decode_losses(maximum_losses);
  }
}

TEST(XDRSRS, EnforcesPublishedInitialParameterSet) {
  EXPECT_TRUE(XDRS(8, XDRSRate::low).Valid());
  EXPECT_TRUE(XDRS(128, XDRSRate::low).Valid());
  EXPECT_FALSE(XDRS(12, XDRSRate::low).Valid());
  EXPECT_FALSE(XDRS(192, XDRSRate::low).Valid());

  EXPECT_TRUE(XDRS(248, XDRSRate::high).Valid());
  EXPECT_TRUE(XDRS(128, XDRSRate::high).Valid());
  EXPECT_FALSE(XDRS(239, XDRSRate::high).Valid());
  EXPECT_FALSE(XDRS(64, XDRSRate::high).Valid());
}

TEST(XDRSRS, RecoversLowAndHighRateData) {
  for (const size_t data_count : {8, 16, 32, 64, 128}) {
    TestXDRSRecovery(data_count, XDRSRate::low, 33);
  }
  for (const size_t data_count : {128, 192, 224, 240, 248}) {
    TestXDRSRecovery(data_count, XDRSRate::high, 17);
  }
}

TEST(XDRSRS, EncodeBackendsMatchAcrossCopyFusionBoundaries) {
  struct Parameters {
    size_t data_count;
    XDRSRate rate;
  };
  constexpr Parameters parameters[] = {
      {8, XDRSRate::low},
      {128, XDRSRate::high},
      {248, XDRSRate::high},
      {255, XDRSRate::high},
  };

  for (const auto [data_count, rate] : parameters) {
    XDRS code(data_count, rate);
    ASSERT_TRUE(code.Valid());
    for (const size_t bytes :
         {size_t{0}, size_t{31}, size_t{32}, size_t{33}, size_t{65}}) {
      const auto data =
          RandomShards(data_count, bytes, 1200 + data_count + 17 * bytes);
      const auto data_pointers = ConstPointers(data);
      std::vector<std::vector<Element>> expected(code.RecoveryCount(),
                                                 std::vector<Element>(bytes));
      auto expected_pointers = MutablePointers(expected);
      std::vector<Element> expected_workspace(code.EncodeWorkspaceSize(bytes));
      ASSERT_EQ(code.Encode(data_pointers, expected_pointers, bytes,
                            expected_workspace, Backend::scalar, Radix::radix2),
                Status::ok);
      if (data_count == 255) {
        for (size_t byte = 0; byte < bytes; ++byte) {
          Element parity = 0;
          for (const auto& shard : data) {
            parity ^= shard[byte];
          }
          EXPECT_EQ(expected[0][byte], parity);
        }
      }

      for (const Backend backend : AvailableBackends()) {
        for (const Radix radix : {Radix::radix2, Radix::radix4}) {
          auto recovery = RandomShards(code.RecoveryCount(), bytes,
                                       1300 + data_count + bytes);
          auto recovery_pointers = MutablePointers(recovery);
          std::vector<Element> workspace(code.EncodeWorkspaceSize(bytes));
          ASSERT_EQ(code.Encode(data_pointers, recovery_pointers, bytes,
                                workspace, backend, radix),
                    Status::ok);
          EXPECT_EQ(recovery, expected)
              << "K=" << data_count << " bytes=" << bytes;
        }
      }
    }
  }
}

TEST(XDRSRS, HalfRateLowAndHighEncodersAgree) {
  constexpr size_t kDataCount = 128;
  const XDRS low(kDataCount, XDRSRate::low);
  const XDRS high(kDataCount, XDRSRate::high);
  ASSERT_TRUE(low.Valid());
  ASSERT_TRUE(high.Valid());

  for (const size_t bytes :
       {size_t{0}, size_t{31}, size_t{32}, size_t{33}, size_t{65}}) {
    const auto data = RandomShards(kDataCount, bytes, 1500 + bytes);
    const auto data_pointers = ConstPointers(data);
    for (const Backend backend : AvailableBackends()) {
      for (const Radix radix : {Radix::radix2, Radix::radix4}) {
        SCOPED_TRACE(testing::Message()
                     << "bytes=" << bytes
                     << " backend=" << static_cast<int>(backend)
                     << " radix=" << static_cast<int>(radix));
        std::vector<std::vector<Element>> low_recovery(
            low.RecoveryCount(), std::vector<Element>(bytes));
        std::vector<std::vector<Element>> high_recovery(
            high.RecoveryCount(), std::vector<Element>(bytes));
        auto low_pointers = MutablePointers(low_recovery);
        auto high_pointers = MutablePointers(high_recovery);
        std::vector<Element> low_workspace(low.EncodeWorkspaceSize(bytes));
        std::vector<Element> high_workspace(high.EncodeWorkspaceSize(bytes));
        ASSERT_EQ(low.Encode(data_pointers, low_pointers, bytes, low_workspace,
                             backend, radix),
                  Status::ok);
        ASSERT_EQ(high.Encode(data_pointers, high_pointers, bytes,
                              high_workspace, backend, radix),
                  Status::ok);
        EXPECT_EQ(low_recovery, high_recovery);
      }
    }
  }
}

TEST(XDRSRS, HighRateMatchesLeopardCodeFamily) {
  for (const size_t data_count : {128, 192, 224, 240, 248}) {
    const size_t recovery_count = 256 - data_count;
    const Leopard leopard(data_count, recovery_count);
    const XDRS xdrs(data_count, XDRSRate::high);
    ASSERT_TRUE(leopard.Valid());
    ASSERT_TRUE(xdrs.Valid());
    for (const size_t bytes : {size_t{0}, size_t{31}, size_t{32}, size_t{33}}) {
      const auto data =
          RandomShards(data_count, bytes, 1900 + data_count + bytes);
      const auto data_pointers = ConstPointers(data);
      for (const Backend backend : AvailableBackends()) {
        for (const Radix radix : {Radix::radix2, Radix::radix4}) {
          std::vector<std::vector<Element>> leopard_recovery(
              recovery_count, std::vector<Element>(bytes));
          std::vector<std::vector<Element>> xdrs_recovery(
              recovery_count, std::vector<Element>(bytes));
          auto leopard_pointers = MutablePointers(leopard_recovery);
          auto xdrs_pointers = MutablePointers(xdrs_recovery);
          std::vector<Element> leopard_workspace(
              leopard.EncodeWorkspaceSize(bytes));
          std::vector<Element> xdrs_workspace(xdrs.EncodeWorkspaceSize(bytes));
          ASSERT_EQ(leopard.Encode(data_pointers, leopard_pointers, bytes,
                                   leopard_workspace, backend, radix),
                    Status::ok);
          ASSERT_EQ(xdrs.Encode(data_pointers, xdrs_pointers, bytes,
                                xdrs_workspace, backend, radix),
                    Status::ok);
          EXPECT_EQ(xdrs_recovery, leopard_recovery)
              << "K=" << data_count << " bytes=" << bytes;
        }
      }
    }
  }
}

TEST(XDRSRS, EncodeWorkspaceContractSurvivesLCHUnification) {
  constexpr size_t kBytes = 17;
  const XDRS low(8, XDRSRate::low);
  const XDRS high(248, XDRSRate::high);
  const XDRS single_recovery(255, XDRSRate::high);
  ASSERT_TRUE(low.Valid());
  ASSERT_TRUE(high.Valid());
  ASSERT_TRUE(single_recovery.Valid());
  EXPECT_EQ(low.EncodeWorkspaceSize(kBytes), 0);
  EXPECT_EQ(high.EncodeWorkspaceSize(kBytes), 8 * kBytes);
  EXPECT_EQ(single_recovery.EncodeWorkspaceSize(kBytes), kBytes);

  const auto data = RandomShards(high.DataCount(), kBytes, 1700);
  const auto data_pointers = ConstPointers(data);
  std::vector<std::vector<Element>> recovery(high.RecoveryCount(),
                                             std::vector<Element>(kBytes));
  auto recovery_pointers = MutablePointers(recovery);
  std::vector<Element> workspace(high.EncodeWorkspaceSize(kBytes));
  ASSERT_EQ(high.Encode(data_pointers, recovery_pointers, kBytes, workspace),
            Status::ok);
  EXPECT_EQ(
      high.Encode(data_pointers, recovery_pointers, kBytes,
                  std::span<Element>(workspace).first(workspace.size() - 1)),
      Status::invalid_argument);
}

TEST(XDRSRS, SingleRecoveryDecodeAndWorkspaceContract) {
  constexpr size_t kDataCount = 255;
  constexpr size_t kBytes = 33;
  const XDRS code(kDataCount, XDRSRate::high);
  ASSERT_TRUE(code.Valid());
  ASSERT_EQ(code.RecoveryCount(), 1);

  const auto data = RandomShards(kDataCount, kBytes, 2200);
  const auto data_pointers = ConstPointers(data);
  std::vector<std::vector<Element>> recovery(1, std::vector<Element>(kBytes));
  auto recovery_pointers = MutablePointers(recovery);
  std::vector<Element> encode_workspace(code.EncodeWorkspaceSize(kBytes));
  ASSERT_EQ(code.Encode(data_pointers, recovery_pointers, kBytes,
                        encode_workspace, Backend::scalar, Radix::radix2),
            Status::ok);

  const size_t workspace_size = code.DecodeWorkspaceSize(kBytes);
  ASSERT_GT(workspace_size, 0);
  for (const Backend backend : AvailableBackends()) {
    for (const Radix radix : {Radix::radix2, Radix::radix4}) {
      auto codeword = recovery;
      codeword.insert(codeword.end(), data.begin(), data.end());
      auto pointers = MutablePointers(codeword);
      std::vector<uint8_t> present(256, 1);
      constexpr size_t kMissingData = 137;
      present[1 + kMissingData] = 0;
      std::fill(codeword[1 + kMissingData].begin(),
                codeword[1 + kMissingData].end(), 0xa5);
      std::vector<Element> workspace(workspace_size);

      EXPECT_EQ(
          code.Decode(pointers, present, kBytes,
                      std::span<Element>(workspace).first(workspace.size() - 1),
                      backend, radix),
          Status::invalid_argument);
      ASSERT_EQ(
          code.Decode(pointers, present, kBytes, workspace, backend, radix),
          Status::ok);
      EXPECT_EQ(codeword[1 + kMissingData], data[kMissingData]);
    }
  }
}

TEST(XDRSRS, DecodeBackendsMatchAcrossSparseHighRateBoundaries) {
  constexpr size_t kDataCount = 248;
  constexpr size_t kRecoveryCount = 8;
  const XDRS code(kDataCount, XDRSRate::high);
  ASSERT_TRUE(code.Valid());
  const std::vector<std::vector<size_t>> loss_patterns{
      {kRecoveryCount},
      {kRecoveryCount, kRecoveryCount + 8, kRecoveryCount + 16,
       kRecoveryCount + 24, kRecoveryCount + 32, kRecoveryCount + 40,
       kRecoveryCount + 48, kRecoveryCount + 56},
      {0, kRecoveryCount, kRecoveryCount + 1, kRecoveryCount + 8,
       kRecoveryCount + 24, kRecoveryCount + 56, kRecoveryCount + 120,
       kRecoveryCount + 232},
  };

  for (const size_t bytes :
       {size_t{0}, size_t{31}, size_t{32}, size_t{33}, size_t{65}}) {
    const auto data = RandomShards(kDataCount, bytes, 1400 + bytes);
    const auto data_pointers = ConstPointers(data);
    std::vector<std::vector<Element>> recovery(kRecoveryCount,
                                               std::vector<Element>(bytes));
    auto recovery_pointers = MutablePointers(recovery);
    std::vector<Element> encode_workspace(code.EncodeWorkspaceSize(bytes));
    ASSERT_EQ(code.Encode(data_pointers, recovery_pointers, bytes,
                          encode_workspace, Backend::scalar, Radix::radix2),
              Status::ok);

    for (const Backend backend : AvailableBackends()) {
      for (const Radix radix : {Radix::radix2, Radix::radix4}) {
        for (const auto& losses : loss_patterns) {
          SCOPED_TRACE(testing::Message()
                       << "bytes=" << bytes
                       << " backend=" << static_cast<int>(backend)
                       << " radix=" << static_cast<int>(radix)
                       << " losses=" << losses.size());
          auto codeword = recovery;
          codeword.insert(codeword.end(), data.begin(), data.end());
          auto pointers = MutablePointers(codeword);
          std::vector<uint8_t> present(256, 1);
          for (const size_t loss : losses) {
            present[loss] = 0;
            std::fill(codeword[loss].begin(), codeword[loss].end(), 0xa5);
          }
          std::vector<Element> workspace(code.DecodeWorkspaceSize(bytes));
          ASSERT_EQ(
              code.Decode(pointers, present, bytes, workspace, backend, radix),
              Status::ok);
          for (size_t i = 0; i < kDataCount; ++i) {
            EXPECT_EQ(codeword[kRecoveryCount + i], data[i]);
          }
        }
      }
    }
  }
}

TEST(XDRSRS, DecodeBackendsMatchAcrossLowRateDerivativeBoundaries) {
  for (const size_t data_count : {size_t{8}, size_t{128}}) {
    const XDRS code(data_count, XDRSRate::low);
    ASSERT_TRUE(code.Valid());
    for (const size_t bytes :
         {size_t{0}, size_t{31}, size_t{32}, size_t{33}, size_t{65}}) {
      EXPECT_EQ(code.DecodeWorkspaceSize(bytes), 256 + 256 * bytes);
      const auto data =
          RandomShards(data_count, bytes, 2100 + data_count + bytes);
      const auto data_pointers = ConstPointers(data);
      std::vector<std::vector<Element>> recovery(code.RecoveryCount(),
                                                 std::vector<Element>(bytes));
      auto recovery_pointers = MutablePointers(recovery);
      std::vector<Element> encode_workspace(code.EncodeWorkspaceSize(bytes));
      ASSERT_EQ(code.Encode(data_pointers, recovery_pointers, bytes,
                            encode_workspace, Backend::scalar, Radix::radix2),
                Status::ok);

      std::vector<size_t> losses;
      losses.reserve(code.RecoveryCount());
      const size_t missing_data = data_count / 2;
      for (size_t i = 0; i < missing_data; ++i) {
        losses.push_back(i);
      }
      for (size_t i = 0; losses.size() < code.RecoveryCount(); ++i) {
        losses.push_back(data_count + i);
      }

      for (const Backend backend : AvailableBackends()) {
        for (const Radix radix : {Radix::radix2, Radix::radix4}) {
          auto codeword = data;
          codeword.insert(codeword.end(), recovery.begin(), recovery.end());
          auto pointers = MutablePointers(codeword);
          std::vector<uint8_t> present(256, 1);
          for (const size_t loss : losses) {
            present[loss] = 0;
            std::fill(codeword[loss].begin(), codeword[loss].end(), 0xa5);
          }
          std::vector<Element> workspace(code.DecodeWorkspaceSize(bytes));
          ASSERT_EQ(
              code.Decode(pointers, present, bytes, workspace, backend, radix),
              Status::ok);
          for (size_t i = 0; i < data_count; ++i) {
            EXPECT_EQ(codeword[i], data[i]);
          }
        }
      }
    }
  }
}

}  // namespace
