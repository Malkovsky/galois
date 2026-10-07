#include "product_monte_carlo_data.h"

#include <gtest/gtest.h>

#include "product_monte_carlo_sha256.h"

TEST(ProductMonteCarloData,
     PostprocessingMetadataCountersAndTransactionalOverflow) {
  mc::Settings settings;
  EXPECT_FALSE(mc::Settings::FromJson(settings.ToJson()).postprocessing);
  EXPECT_FALSE(settings.ToJson().contains("postprocessing"));
  settings.postprocessing = true;
  EXPECT_TRUE(mc::Settings::FromJson(settings.ToJson()).postprocessing);
  auto invalid = settings.ToJson();
  invalid["postprocessing"] = 1;
  EXPECT_THROW(mc::Settings::FromJson(invalid), std::runtime_error);
  mc::Stats increment{1, 2, 0, 0, 0, 0, {1, 2, 2}};
  EXPECT_EQ(mc::Stats::FromJson(increment.ToJson(settings), 1, 0, 16, settings)
                .postprocessing,
            increment.postprocessing);
  mc::Aggregate aggregate("pp", 2, settings);
  aggregate.Add(0, increment);
  aggregate.Add(1, increment);
  EXPECT_EQ(aggregate.overall.postprocessing,
            (std::array<uint64_t, 3>{2, 4, 4}));
  EXPECT_EQ(aggregate.by_k.at(0).postprocessing, increment.postprocessing);
  for (size_t i = 0; i < 3; ++i) {
    mc::Stats huge;
    huge.postprocessing[i] = UINT64_MAX;
    const auto before = aggregate.Summary();
    EXPECT_THROW(aggregate.Add(0, huge), std::runtime_error);
    EXPECT_EQ(aggregate.Summary(), before);
    // Also exercise failure in per-k staging after overall succeeds.
    auto copy = aggregate;
    copy.by_k.at(0).postprocessing[i] = UINT64_MAX;
    const auto saved = copy.Summary();
    EXPECT_THROW(copy.Add(0, increment), std::runtime_error);
    EXPECT_EQ(copy.Summary(), saved);
  }
  auto old = mc::Stats{}.ToJson();
  for (auto name : mc::kPostprocessingMetrics) {
    old.erase(name);
  }
  EXPECT_NO_THROW(mc::Stats::FromJson(old, 0, 0, 16));
  EXPECT_THROW(mc::Stats::FromJson(old, 0, 0, 16, settings),
               std::runtime_error);
  old[mc::kPostprocessingMetrics[0]] = 0;
  EXPECT_THROW(mc::Stats::FromJson(old, 0, 0, 16), std::runtime_error);
}

TEST(ProductMonteCarloData, Sha256KnownVectorsAndPadding) {
  const auto check = [](std::string_view input, std::string_view expected) {
    const auto digest = mc::Sha256(
        {reinterpret_cast<const uint8_t*>(input.data()), input.size()});
    constexpr char digits[] = "0123456789abcdef";
    std::string hex;
    for (auto byte : digest) {
      hex += digits[byte >> 4];
      hex += digits[byte & 15];
    }
    EXPECT_EQ(hex, expected) << "input length " << input.size();
  };
  check({}, "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
  check("abc",
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
  check("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq",
        "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
  check(
      "abcdefghbcdefghicdefghijdefghijkefghijklfghijklmghijklmnhijklmno"
      "ijklmnopjklmnopqklmnopqrlmnopqrsmnopqrstnopqrstu",
      "cf5b16a778af8380036ce59e7b0492370b249b11e8f07a51afac45037afee9d1");
  check(std::string(1000000, 'a'),
        "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
  // Python hashlib.sha256(bytes(range(n))).hexdigest(), including NUL bytes.
  const std::pair<size_t, std::string_view> boundaries[] = {
      {55, "463eb28e72f82e0a96c0a4cc53690c571281131f672aa229e0d45ae59b598b59"},
      {56, "da2ae4d6b36748f2a318f23e7ab1dfdf45acdc9d049bd80e59de82a60895f562"},
      {63, "29af2686fd53374a36b0846694cc342177e428d1647515f078784d69cdb9e488"},
      {64, "fdeab9acf3710362bd2658cdc9a29e8f9c757fcf9811603a8c447cd1d9151108"},
      {65, "4bfd2c8b6f1eec7a2afeb48b934ee4b2694182027e6d0fc075074f2fabb31781"}};
  for (const auto& [length, expected] : boundaries) {
    std::string input(length, '\0');
    for (size_t i = 0; i < length; ++i) {
      input[i] = static_cast<char>(i);
    }
    check(input, expected);
  }
}

TEST(ProductMonteCarloData, DimensionsMetadataAndDynamicOverflowBounds) {
  mc::Settings settings;
  const auto legacy = settings.ToJson();
  EXPECT_FALSE(legacy.contains("n1"));
  EXPECT_EQ(mc::Settings::FromJson(legacy).ToJson(), legacy);
  settings.n2 = 175;
  settings.k2 = 173;
  settings.Validate();
  EXPECT_EQ(settings.FullBits(), 358400u);
  EXPECT_EQ(settings.InfoBits(), 310016u);
  EXPECT_EQ(mc::Settings::FromJson(settings.ToJson()).ToJson(),
            settings.ToJson());
  auto partial = settings.ToJson();
  partial.erase("k2");
  EXPECT_THROW(mc::Settings::FromJson(partial), std::runtime_error);
  const uint64_t count = UINT64_MAX / settings.FullBits();
  mc::Stats stats{count, 2 * count, count, 0, count, 0};
  const auto json = stats.ToJson(settings);
  EXPECT_NO_THROW(mc::Stats::FromJson(json, count, 1, 16, settings));
  EXPECT_THROW(mc::Stats::FromJson(json, count, 1, 16), std::runtime_error);
  mc::Aggregate aggregate("shortened", 2, settings);
  aggregate.Add(1, stats);
  const auto before = aggregate.Summary();
  EXPECT_EQ(before.at("code parameters"),
            (mc::Json{{"n1", 256}, {"k1", 224}, {"n2", 175}, {"k2", 173}}));
  EXPECT_THROW(aggregate.Add(1, mc::Stats{1}), std::runtime_error);
  EXPECT_EQ(aggregate.Summary(), before);
  EXPECT_THROW(stats.Add(mc::Stats{1}, settings.FullBits()),
               std::runtime_error);
  EXPECT_EQ(stats.ToJson(settings), json);
  settings.n1 = 4;
  settings.k1 = 2;
  settings.n2 = 5;
  settings.k2 = 3;
  EXPECT_THROW(settings.Validate(), std::runtime_error);
  settings.lo = 0;
  settings.hi = settings.FullBits();
  EXPECT_NO_THROW(settings.Validate());
}

TEST(ProductMonteCarloData, ExactBoundedIntegers) {
  const uint64_t n = UINT64_MAX / 524288;
  mc::Stats s;
  s.blocks = n;
  s.iterations = n * 2;
  s.info_raw = n;
  s.full_raw = n;
  s.info_post = 0;
  s.full_post = 0;
  auto json = s.ToJson();
  auto text = mc::Dump(json);
  EXPECT_NE(text.find(std::to_string(n)), std::string::npos);
  auto parsed = mc::Parse(text);
  auto recovered = mc::Stats::FromJson(parsed, n, 1, 16);
  EXPECT_THROW(recovered.Add(s), std::runtime_error);
  EXPECT_EQ(recovered.ToJson(), s.ToJson());
  EXPECT_EQ(recovered.blocks, n);
  EXPECT_EQ(
      mc::Natural(recovered.ToJson().at("information bits").at("total bits")),
      n * 455168);
  EXPECT_EQ(mc::Dump(mc::Parse(mc::Dump(recovered.ToJson(), true))),
            mc::Dump(recovered.ToJson()));
  s.iterations = UINT64_MAX;
  EXPECT_NO_THROW(
      mc::Stats::FromJson(mc::Parse(mc::Dump(s.ToJson())), n, 1, 1000000));
  EXPECT_THROW(mc::Stats::FromJson(s.ToJson(), n, 1, 16), std::runtime_error);
}

TEST(ProductMonteCarloData, CheckedArithmeticAndTransactionalUpdates) {
  EXPECT_EQ(mc::CheckedAdd(UINT64_MAX - 1, 1), UINT64_MAX);
  EXPECT_EQ(mc::CheckedMultiply(UINT64_MAX, 1), UINT64_MAX);
  EXPECT_EQ(mc::CheckedMultiply(UINT64_MAX, 0), 0);
  EXPECT_THROW(mc::CheckedAdd(UINT64_MAX, 1), std::runtime_error);
  EXPECT_THROW(mc::CheckedMultiply(UINT64_MAX, 2), std::runtime_error);
  EXPECT_EQ(mc::Decimal("18446744073709551615"), UINT64_MAX);
  EXPECT_EQ(mc::Decimal("00001"), 1);
  for (auto text : {"18446744073709551616", "-1", "", "1.0", "+1", "1 "}) {
    EXPECT_THROW(mc::Decimal(text), std::runtime_error);
  }
  for (auto member :
       {&mc::Stats::iterations, &mc::Stats::info_raw, &mc::Stats::info_post,
        &mc::Stats::full_raw, &mc::Stats::full_post}) {
    mc::Stats s, delta;
    s.*member = UINT64_MAX;
    delta.blocks = 1;
    delta.*member = 1;
    const auto before = s.ToJson();
    EXPECT_THROW(s.Add(delta), std::runtime_error);
    EXPECT_EQ(s.ToJson(), before);
  }
  mc::Stats s;
  s.blocks = UINT64_MAX;
  EXPECT_THROW(s.Add(mc::Stats{1}), std::runtime_error);
  EXPECT_EQ(s.blocks, UINT64_MAX);
  EXPECT_THROW(s.ToJson(), std::runtime_error);
  s.blocks = UINT64_MAX / 524288;
  const auto before = s.ToJson();
  EXPECT_THROW(s.Add(std::array<uint64_t, 22>{}), std::runtime_error);
  EXPECT_EQ(s.ToJson(), before);
  auto j = before;
  j["completed blocks"] = s.blocks + 1;
  j["total iterations"] = (s.blocks + 1) * 2;
  EXPECT_THROW(mc::Stats::FromJson(j, s.blocks + 1, 0, 16), std::runtime_error);
}

TEST(ProductMonteCarloData, LegacyBoundedMomentsAndTransactionality) {
  auto old = mc::LegacyStats();
  old[mc::kMetrics.back()]["squared sum"] = UINT64_MAX;
  auto delta = mc::LegacyStats();
  delta[mc::kMetrics.back()]["squared sum"] = uint64_t{1};
  const auto before = old;
  EXPECT_THROW(mc::AddLegacy(old, delta), std::runtime_error);
  EXPECT_EQ(old, before);
  std::array<uint64_t, 22> m{};
  m.back() = UINT64_MAX;
  EXPECT_THROW(mc::AddLegacyTrial(old, m), std::runtime_error);
  EXPECT_EQ(old, before);
  // A theoretical bound may exceed u64 even though every stored moment fits.
  old = mc::LegacyStats();
  m = {};
  m[12] = 2;
  mc::AddLegacyTrial(old, m);
  EXPECT_NO_THROW(mc::ValidateLegacy(old, 1, 0, 1000000));
  old = mc::LegacyStats();
  old[mc::kMetrics[13]]["sum"] = UINT64_MAX;
  old[mc::kMetrics[13]]["squared sum"] = UINT64_MAX;
  old[mc::kMetrics[15]] = old[mc::kMetrics[13]];
  EXPECT_NO_THROW(mc::ValidateLegacy(old, UINT64_MAX, 0, 1000000));
}

TEST(ProductMonteCarloData, StrictJsonAndNaturalNumbers) {
  for (const auto* text :
       {"{\"a\":1,\"a\":2}", "{\"a\":{\"b\":0,\"b\":1}}", "1.0", "1e3", "NaN",
        "/*comment*/1", "[1,]", "1 2", "Infinity", "1e9999",
        "{\"a\":0,\"\\u0061\":1}", "{\"squared sum\":18446744073709551616}"}) {
    EXPECT_THROW(mc::Parse(text), std::exception) << text;
  }
  for (const auto* text : {"true", "\"123\"", "-1", "-184467440737095516160"}) {
    EXPECT_THROW(mc::Natural(mc::Parse(text)), std::exception) << text;
  }
  EXPECT_EQ(mc::U64(mc::Parse("18446744073709551615")), UINT64_MAX);
  EXPECT_EQ(mc::Natural(mc::Json(1)), 1);
  EXPECT_EQ(mc::Dump(mc::Number(UINT64_MAX)), "18446744073709551615");
  EXPECT_THROW(mc::Natural(mc::Json(1.0)), std::exception);
  EXPECT_THROW(mc::U64(mc::Parse("18446744073709551616")), std::exception);
  EXPECT_EQ(mc::Dump(mc::Parse("{\"z\":1,\"a\":\"x\"}")),
            "{\"a\":\"x\",\"z\":1}");
  EXPECT_EQ(
      mc::Dump(mc::Parse(
          R"({"z":true,"a":"\u00e9\u000f/\ud83d\ude00","n":18446744073709551615})")),
      R"({"a":"\u00e9\u000f/\ud83d\ude00","n":18446744073709551615,"z":true})");
  EXPECT_NO_THROW(mc::Parse(std::string(32, '[') + "0" + std::string(32, ']')));
  EXPECT_THROW(mc::Parse(std::string(33, '[') + "0" + std::string(33, ']')),
               std::exception);
  EXPECT_NO_THROW(mc::Parse(R"([{"a":1},{"a":2,"b":[{"a":3}]}])"));
}

TEST(ProductMonteCarloData, AggregatePreparesOnlyOneStratum) {
  for (unsigned schema : {1, 2}) {
    mc::Aggregate aggregate("test", schema);
    const auto old = mc::LegacyStats();
    for (uint64_t k = 0; k < 1024; ++k) {
      aggregate.Add(k, mc::Stats{1}, old);
    }
    const auto* untouched = &aggregate.by_k.at(0);
    const auto before = aggregate.Summary();
    for (uint64_t k : {512, 1024}) {
      {
        auto abandoned = aggregate.PrepareAdd(k, mc::Stats{1}, old);
        EXPECT_EQ(abandoned.row.key(), k);
        EXPECT_EQ(aggregate.Summary(), before);
      }
      EXPECT_EQ(aggregate.Summary(), before);
    }
    aggregate.Commit(aggregate.PrepareAdd(512, mc::Stats{1}, old));
    aggregate.Add(1024, mc::Stats{1}, old);
    EXPECT_EQ(aggregate.by_k.size(), 1025);
    EXPECT_EQ(aggregate.overall.blocks, 1026);
    EXPECT_EQ(aggregate.by_k.at(512).blocks, 2);
    EXPECT_EQ(aggregate.by_k.at(1024).blocks, 1);
    EXPECT_EQ(&aggregate.by_k.at(0), untouched);
    if (schema == 1) {
      EXPECT_EQ(aggregate.legacy_k.size(), 1025);
    }
  }
}

TEST(ProductMonteCarloData, AggregateCounterOverflowIsTransactional) {
  for (auto member :
       {&mc::Stats::blocks, &mc::Stats::iterations, &mc::Stats::info_raw,
        &mc::Stats::info_post, &mc::Stats::full_raw, &mc::Stats::full_post}) {
    for (bool overall : {false, true}) {
      mc::Aggregate aggregate("test", 2);
      aggregate.Add(7, mc::Stats{});
      const uint64_t maximum =
          member == &mc::Stats::blocks ? UINT64_MAX / 524288 : UINT64_MAX;
      if (overall) {
        aggregate.overall.*member = maximum;
        ASSERT_EQ(aggregate.overall.*member, maximum);
      } else {
        aggregate.by_k.at(7).*member = maximum;
        ASSERT_EQ(aggregate.by_k.at(7).*member, maximum);
      }
      const auto before = aggregate.Summary();
      mc::Stats delta;
      delta.*member = 1;
      EXPECT_THROW(aggregate.Add(7, delta), std::runtime_error);
      EXPECT_EQ(aggregate.Summary(), before);
      if (overall) {
        EXPECT_THROW(aggregate.Add(8, delta), std::runtime_error);
        EXPECT_EQ(aggregate.Summary(), before);
        EXPECT_FALSE(aggregate.by_k.contains(8));
      }
    }
  }
}

TEST(ProductMonteCarloData, AggregateLegacyOverflowIsTransactional) {
  for (const auto* field : {"sum", "squared sum"}) {
    for (bool overall : {false, true}) {
      mc::Aggregate aggregate("test", 1);
      auto delta = mc::LegacyStats();
      aggregate.Add(7, mc::Stats{1}, delta);
      (overall ? aggregate.legacy
               : aggregate.legacy_k.at(7))[mc::kMetrics.back()][field] =
          UINT64_MAX;
      delta[mc::kMetrics.back()][field] = uint64_t{1};
      const auto before = aggregate.Summary();
      EXPECT_THROW(aggregate.Add(7, mc::Stats{1}, delta), std::runtime_error);
      EXPECT_EQ(aggregate.Summary(), before);
      EXPECT_EQ(aggregate.overall.blocks, 1);
      if (overall) {
        EXPECT_THROW(aggregate.Add(8, mc::Stats{1}, delta), std::runtime_error);
        EXPECT_EQ(aggregate.Summary(), before);
        EXPECT_FALSE(aggregate.by_k.contains(8));
        EXPECT_FALSE(aggregate.legacy_k.contains(8));
      }
    }
  }
}
