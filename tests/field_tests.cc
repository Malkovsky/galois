#include "field.h"
#include "gf256/gf256.h"

#include <algorithm>
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

// TEST(GF_2_8, GFNI) {
//   gf_2_8::InitGFNI();
//   for (uint16_t x = 0; x < 256; ++x) {
//     for (uint16_t y = 0; y < 256; ++y) {
//       ASSERT_EQ(gf_2_8::Multiply(x, y), gf_2_8::MultiplyGFNI(x, y));
//     }
//   }
// }

TEST(GF_2_8, RowMulAdd) {
  gf_2_8::InitGFNI();
  gf_2_8::Init();
  gf256_init();
  std::mt19937 rng(42);
  constexpr size_t length = 1000;
  std::array<gf_2_8::element_t, length> data, x, y, result, ref;
  for (size_t i = 0; i < length; ++i) {
    for (size_t j = 0; j < length; ++j) {
      data[j] = rng();
      y[j] = rng();
    }
    gf_2_8::element_t z = 1;

    std::copy(data.begin(), data.end(), x.begin());
    gf_2_8::AddScaledRowBase(x.data(), y.data(), z, length);
    std::copy(x.begin(), x.end(), ref.begin());

    std::copy(data.begin(), data.end(), x.begin());
    gf_2_8::AddScaledRowGFNIGeneral(x.data(), y.data(), z, length);
    ASSERT_EQ(std::equal(ref.begin(), ref.end(), x.begin()), true);

    std::copy(data.begin(), data.end(), x.begin());
    gf_2_8::AddScaledRowGFNIDedicated(x.data(), y.data(), z, length);
    ASSERT_EQ(std::equal(ref.begin(), ref.end(), x.begin()), true);

    std::copy(data.begin(), data.end(), x.begin());
    gf_2_8::AddScaledRowSIMD(x.data(), y.data(), z, length);
    ASSERT_EQ(std::equal(ref.begin(), ref.end(), x.begin()), true);
  }
}

TEST(GF_2_8, Inverse) {
  gf_2_8::Init();
  for (uint16_t x = 1; x < 256; ++x) {
    ASSERT_EQ(gf_2_8::Multiply(x, gf_2_8::Inv(x)), gf_2_8::One());
  }
}

TEST(GF_2_8, MatMul) {
  gf_2_8::Init();
  gf_2_8::InitGFNI();
  std::mt19937 rng(42);

  std::vector<gf_2_8::element_t> left;
  std::vector<gf_2_8::element_t> right;
  std::vector<gf_2_8::element_t> result;
  std::vector<gf_2_8::element_t> ref;

  for (size_t n = 5; n < 10; ++n) {
    for (size_t m = 7; m < 12; ++m) {
      for (size_t l = 11; l < 17; ++l) {
        left.resize(n * m);
        right.resize(m * l);
        result.resize(n * l);
        ref.resize(n * l);
        for (auto& x: left) {
          x = rng();
        }
        for (auto& x: right) {
          x = rng();
        }
        
        gf_2_8::MatMul(left.data(), right.data(), n, m, l,
                       gf_2_8::AddScaledRowBase, result.data());
        std::fill(ref.begin(), ref.end(), 0);
        for (size_t i = 0; i < n; ++i) {
          for (size_t j = 0; j < l; ++j) {
            for (size_t k = 0; k < m; ++k) {
              ref[i * l + j] ^=
                  gf_2_8::Multiply(left[i * m + k], right[k * l + j]);
            }
          }
        }
        ASSERT_EQ(std::equal(ref.begin(), ref.end(), result.begin()), true);
      }
    }
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