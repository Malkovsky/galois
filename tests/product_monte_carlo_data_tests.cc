#include "product_monte_carlo_data.h"

#include <gtest/gtest.h>

TEST(ProductMonteCarloData, ExactUnboundedIntegers) {
  const auto n = mc::Big::from_string(
      "184467440737095516160000000000000000000000000000001");
  mc::Stats s;
  s.blocks = n;
  s.iterations = n * 2;
  s.info_raw = n;
  s.full_raw = n;
  s.info_post = 0;
  s.full_post = 0;
  auto json = s.ToJson();
  auto text = mc::Dump(json);
  EXPECT_NE(text.find(n.to_string()), std::string::npos);
  auto parsed = mc::Parse(text);
  auto recovered = mc::Stats::FromJson(parsed, n, 1, 16);
  recovered.Add(s);
  EXPECT_EQ(recovered.blocks, n * 2);
  EXPECT_EQ(
      mc::Natural(recovered.ToJson().at("information bits").at("total bits")),
      n * 910336);
  EXPECT_EQ(mc::Dump(mc::Parse(mc::Dump(recovered.ToJson(), true))),
            mc::Dump(recovered.ToJson()));
}

TEST(ProductMonteCarloData, StrictJsonAndNaturalNumbers) {
  for (const auto* text :
       {"{\"a\":1,\"a\":2}", "{\"a\":{\"b\":0,\"b\":1}}", "1.0", "1e3", "NaN",
        "/*comment*/1", "[1,]", "1 2"}) {
    EXPECT_THROW(mc::Parse(text), std::exception) << text;
  }
  for (const auto* text : {"true", "\"123\"", "-1", "-184467440737095516160"}) {
    EXPECT_THROW(mc::Natural(mc::Parse(text)), std::exception) << text;
  }
  EXPECT_EQ(mc::U64(mc::Parse("18446744073709551615")), UINT64_MAX);
  EXPECT_THROW(mc::U64(mc::Parse("18446744073709551616")), std::exception);
  EXPECT_EQ(mc::Dump(mc::Parse("{\"z\":1,\"a\":\"x\"}")),
            "{\"a\":\"x\",\"z\":1}");
}
