#include "field.h"

#include <random>

#include "gtest/gtest.h"

namespace {

TEST(GF_2_8, Distributivity) {
  for (uint16_t x = 0; x < 256; ++x) {
    for (uint16_t y = 0; y < 256; ++y) {
      for (uint16_t z = 0; z < 256; ++z) {
        ASSERT_EQ(gf_2_8::Multiply(z, gf_2_8::Add(x, y)),
                  gf_2_8::Add(gf_2_8::Multiply(z, x), gf_2_8::Multiply(z, y)));
      }
    }
  }
}

TEST(GF_2_8, LUT) {
  gf_2_8::Init();
  for (uint16_t x = 0; x < 256; ++x) {
    for (uint16_t y = 0; y < 256; ++y) {
      ASSERT_EQ(gf_2_8::Multiply(x, y), gf_2_8::MultiplyLUT(x, y));
    }
  }
}

TEST(GF_2_8, GFNI) {
  gf_2_8::InitGFNI();
  for (uint16_t x = 0; x < 256; ++x) {
    for (uint16_t y = 0; y < 256; ++y) {
      ASSERT_EQ(gf_2_8::Multiply(x, y), gf_2_8::MultiplyGFNI(x, y));
    }
  }
}

TEST(GF_2_8, Inverse) {
  gf_2_8::Init();
  for (uint16_t x = 1; x < 256; ++x) {
    ASSERT_EQ(gf_2_8::Multiply(x, gf_2_8::Inv(x)), gf_2_8::One());
  }
}

TEST(GF_2_8, Irreducibly) {
  for (uint16_t x = 0; x < 256; ++x) {
    ASSERT_NE(gf_2_8::Add(gf_2_8::Add(gf_2_8::Multiply(x, x), x), 0x20),
              gf_2_8::Zero());
  }
}

TEST(GF_2_16, Distributivity) {
  std::mt19937 rng(42);
  for (size_t i = 0; i < 100000; ++i) {
    gf_2_16::element_t x = rng();
    gf_2_16::element_t y = rng();
    gf_2_16::element_t z = rng();
    ASSERT_EQ(gf_2_16::Multiply(z, gf_2_16::Add(x, y)),
              gf_2_16::Add(gf_2_16::Multiply(z, x), gf_2_16::Multiply(z, y)));
  }
}

TEST(GF_2_16, Inverse) {
  gf_2_8::Init();
  for (uint16_t x = 1; x != 0; ++x) {
    ASSERT_EQ(gf_2_16::Multiply(x, gf_2_16::Inv(x)), gf_2_16::One());
  }
}

TEST(GF_2_16, Pow_2_8) {
  gf_2_8::Init();
  for (uint16_t x = 1; x != 0; ++x) {
    ASSERT_EQ(gf_2_16::Pow(x, 257) >> 8, 0);
  }  
}

TEST(GF_2_16, InverseIT) {
  gf_2_8::Init();
  for (uint16_t x = 1; x != 0; ++x) {
    ASSERT_EQ(gf_2_16::Inv(x), gf_2_16::InvIT(x));
  }
}



} // namespace