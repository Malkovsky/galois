#include "product_monte_carlo_data.h"

#include <gtest/gtest.h>

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
