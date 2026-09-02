#include "field.h"

#include <algorithm>
#include <array>
#include <limits>
#include <random>
#include <vector>

#include "gtest/gtest.h"
#include "matrix.h"

namespace {

gf2p8::Element PowStandardPositiveReference(gf2p8::Element value,
                                            uint64_t exponent) {
  gf2p8::Element result = 1;
  while (exponent != 0) {
    if ((exponent & 1U) != 0) {
      result = gf2p8::detail::MultiplyStandardDirect(result, value);
    }
    value = gf2p8::detail::MultiplyStandardDirect(value, value);
    exponent >>= 1;
  }
  return result;
}

gf2p8::Element PowStandardReference(gf2p8::Element value, int exponent) {
  if (value == 0) {
    return exponent == 0 ? 1 : 0;
  }
  int64_t magnitude = exponent;
  if (magnitude < 0) {
    value = PowStandardPositiveReference(value, 254);
    magnitude = -magnitude;
  }
  return PowStandardPositiveReference(value, static_cast<uint64_t>(magnitude));
}

TEST(GF2P8, Distributivity) {
  for (uint16_t x = 0; x < 256; ++x) {
    for (uint16_t y = 0; y < 256; ++y) {
      for (uint16_t z = 0; z < 256; ++z) {
        ASSERT_EQ(gf2p8::MultiplyStandard(z, gf2p8::Add(x, y)),
                  gf2p8::Add(gf2p8::MultiplyStandard(z, x),
                             gf2p8::MultiplyStandard(z, y)));
        ASSERT_EQ(gf2p8::MultiplyCantor(z, gf2p8::Add(x, y)),
                  gf2p8::Add(gf2p8::MultiplyCantor(z, x),
                             gf2p8::MultiplyCantor(z, y)));
      }
    }
  }
}

TEST(GF2P8, ExplicitBasisMultiplicationAndConversion) {
  for (uint16_t x = 0; x < 256; ++x) {
    for (uint16_t y = 0; y < 256; ++y) {
      const auto standard_x = static_cast<gf2p8::Element>(x);
      const auto standard_y = static_cast<gf2p8::Element>(y);
      const auto cantor_x = gf2p8::StandardToCantor(standard_x);
      const auto cantor_y = gf2p8::StandardToCantor(standard_y);
      ASSERT_EQ(gf2p8::MultiplyStandard(standard_x, standard_y),
                gf2p8::detail::MultiplyStandardDirect(standard_x, standard_y));
      ASSERT_EQ(gf2p8::MultiplyCantor(cantor_x, cantor_y),
                gf2p8::detail::MultiplyCantorDirect(cantor_x, cantor_y));
      ASSERT_EQ(gf2p8::StandardToCantor(
                    gf2p8::MultiplyStandard(standard_x, standard_y)),
                gf2p8::MultiplyCantor(cantor_x, cantor_y));
    }
  }
}

TEST(GF2P8, NativeCantorCoordinateMapsAndTables) {
  constexpr std::array<gf2p8::Element, 8> kCantorPolynomialBasis{
      0x01, 0xd6, 0x98, 0x92, 0x56, 0xc8, 0x58, 0xe6};
  EXPECT_EQ(gf2p8::detail::kCantorPolynomialBasis, kCantorPolynomialBasis);
  for (uint16_t value = 0; value < 256; ++value) {
    const auto element = static_cast<gf2p8::Element>(value);
    EXPECT_EQ(gf2p8::StandardToCantor(gf2p8::CantorToStandard(element)),
              element);
    EXPECT_EQ(gf2p8::CantorToStandard(gf2p8::StandardToCantor(element)),
              element);
  }

  const auto& tables = gf2p8::Tables();
  for (uint16_t coefficient = 0; coefficient < 256; ++coefficient) {
    const auto c = static_cast<gf2p8::Element>(coefficient);
    for (uint16_t value = 0; value < 256; ++value) {
      const auto x = static_cast<gf2p8::Element>(value);
      const auto& row = tables.shuffle[c];
      EXPECT_EQ(static_cast<gf2p8::Element>(row[x & 0x0f] ^ row[32 + (x >> 4)]),
                gf2p8::MultiplyCantor(c, x));

      gf2p8::Element affine_product = 0;
      for (size_t output_bit = 0; output_bit < 8; ++output_bit) {
        const auto mask = static_cast<gf2p8::Element>(tables.affine[c] >>
                                                      (8 * (7 - output_bit)));
        affine_product |= static_cast<gf2p8::Element>(
            (std::popcount(static_cast<unsigned>(mask & x)) & 1U)
            << output_bit);
      }
      EXPECT_EQ(affine_product, gf2p8::MultiplyCantor(c, x));
    }
  }
}

TEST(GF2P8, RowMulAdd) {
  std::mt19937 rng(42);
  constexpr size_t length = 1000;
  std::array<gf2p8::Element, length> data, x, y, result, ref;
  for (size_t i = 0; i < length; ++i) {
    for (size_t j = 0; j < length; ++j) {
      data[j] = rng();
      y[j] = rng();
    }
    gf2p8::Element z = rng();

    for (size_t j = 0; j < length; ++j) {
      ref[j] = data[j] ^ gf2p8::MultiplyCantor(y[j], z);
    }

    std::copy(data.begin(), data.end(), x.begin());
    gf2p8::AddScaledRow(x.data(), y.data(), z, length);
    ASSERT_EQ(std::equal(ref.begin(), ref.end(), x.begin()), true);
  }
}

TEST(GF2P8, Inverse) {
  for (uint16_t x = 1; x < 256; ++x) {
    const auto standard = static_cast<gf2p8::Element>(x);
    const auto cantor = gf2p8::StandardToCantor(standard);
    ASSERT_EQ(gf2p8::MultiplyStandard(standard, gf2p8::InvStandard(standard)),
              gf2p8::One());
    ASSERT_EQ(gf2p8::MultiplyCantor(cantor, gf2p8::InvCantor(cantor)),
              gf2p8::One());
    ASSERT_EQ(gf2p8::StandardToCantor(gf2p8::InvStandard(standard)),
              gf2p8::InvCantor(cantor));
  }
}

TEST(GF2P8, DivisionAndPowersPreserveCoordinateMapping) {
  constexpr std::array<int, 10> kExponents{
      std::numeric_limits<int>::min(), -257, -1, 0, 1, 2, 17, 254, 255, 256};
  for (uint16_t a = 0; a < 256; ++a) {
    const auto standard_a = static_cast<gf2p8::Element>(a);
    const auto cantor_a = gf2p8::StandardToCantor(standard_a);
    for (const int exponent : kExponents) {
      const auto expected = PowStandardReference(standard_a, exponent);
      EXPECT_EQ(gf2p8::PowStandard(standard_a, exponent), expected);
      EXPECT_EQ(gf2p8::CantorToStandard(gf2p8::PowCantor(cantor_a, exponent)),
                expected);
    }
    for (uint16_t b = 1; b < 256; ++b) {
      const auto standard_b = static_cast<gf2p8::Element>(b);
      const auto cantor_b = gf2p8::StandardToCantor(standard_b);
      const auto standard_quotient = gf2p8::DivStandard(standard_a, standard_b);
      const auto cantor_quotient = gf2p8::DivCantor(cantor_a, cantor_b);
      EXPECT_EQ(
          gf2p8::detail::MultiplyStandardDirect(standard_quotient, standard_b),
          standard_a);
      EXPECT_EQ(gf2p8::detail::MultiplyCantorDirect(cantor_quotient, cantor_b),
                cantor_a);
      EXPECT_EQ(gf2p8::StandardToCantor(standard_quotient), cantor_quotient);
    }
  }
}

TEST(GF2P8, MatMul) {
  std::mt19937 rng(42);

  std::vector<gf2p8::Element> left;
  std::vector<gf2p8::Element> right;
  std::vector<gf2p8::Element> result;
  std::vector<gf2p8::Element> ref;

  for (size_t n = 5; n < 10; ++n) {
    for (size_t m = 7; m < 12; ++m) {
      for (size_t l = 11; l < 17; ++l) {
        left.resize(n * m);
        right.resize(m * l);
        result.resize(n * l);
        ref.resize(n * l);
        for (auto& x : left) {
          x = rng();
        }
        for (auto& x : right) {
          x = rng();
        }

        gf2p8::MatMul(left.data(), right.data(), n, m, l, gf2p8::AddScaledRow,
                      result.data());
        std::fill(ref.begin(), ref.end(), 0);
        for (size_t i = 0; i < n; ++i) {
          for (size_t j = 0; j < l; ++j) {
            for (size_t k = 0; k < m; ++k) {
              ref[i * l + j] ^=
                  gf2p8::MultiplyCantor(left[i * m + k], right[k * l + j]);
            }
          }
        }
        ASSERT_EQ(std::equal(ref.begin(), ref.end(), result.begin()), true);
      }
    }
  }
}

TEST(GF2P8, EmptyMatricesAcceptNullStorage) {
  gf2p8::MatMul(nullptr, nullptr, 0, 0, 0, gf2p8::AddScaledRow, nullptr);
  gf2p8::MatMulBlockedLUT(nullptr, nullptr, 0, 0, 0, nullptr);
  gf2p8::MatMulBlockedGFNI(nullptr, nullptr, 0, 0, 0, nullptr);

  std::array<gf2p8::Element, 1> storage{};
  gf2p8::MatMul(storage.data(), storage.data(), 0, 1, 1, gf2p8::AddScaledRow,
                nullptr);
  gf2p8::MatMulBlockedLUT(storage.data(), storage.data(), 1, 1, 0, nullptr);
  gf2p8::MatMulBlockedGFNI(storage.data(), storage.data(), 1, 1, 0, nullptr);
}

TEST(GF2P8, MatMulBlockedGFNI) {
  std::mt19937 rng(42);

  struct Shape {
    size_t rows;
    size_t inner;
    size_t cols;
  };

  const Shape shapes[] = {
      {1, 0, 3},  {1, 1, 1},    {3, 5, 7},    {4, 6, 64},
      {5, 9, 65}, {7, 11, 130}, {9, 13, 512}, {3, 7, 2053},
  };

  std::vector<gf2p8::Element> left;
  std::vector<gf2p8::Element> right;
  std::vector<gf2p8::Element> result;
  std::vector<gf2p8::Element> ref;

  for (const auto shape : shapes) {
    left.resize(shape.rows * shape.inner);
    right.resize(shape.inner * shape.cols);
    result.resize(shape.rows * shape.cols);
    ref.resize(shape.rows * shape.cols);

    for (auto& x : left) {
      x = rng();
    }
    for (auto& x : right) {
      x = rng();
    }

    gf2p8::MatMulBlockedGFNI(left.data(), right.data(), shape.rows, shape.inner,
                             shape.cols, result.data());

    std::fill(ref.begin(), ref.end(), 0);
    for (size_t i = 0; i < shape.rows; ++i) {
      for (size_t j = 0; j < shape.cols; ++j) {
        for (size_t k = 0; k < shape.inner; ++k) {
          ref[i * shape.cols + j] ^= gf2p8::MultiplyCantor(
              left[i * shape.inner + k], right[k * shape.cols + j]);
        }
      }
    }

    ASSERT_EQ(std::equal(ref.begin(), ref.end(), result.begin()), true)
        << "shape=" << shape.rows << "x" << shape.inner << "x" << shape.cols;
  }
}

TEST(GF2P8, MatMulBlockedLUT) {
  std::mt19937 rng(43);

  struct Shape {
    size_t rows;
    size_t inner;
    size_t cols;
  };

  const Shape shapes[] = {
      {1, 0, 3},  {1, 1, 1},    {3, 5, 7},    {4, 6, 64},
      {5, 9, 65}, {7, 11, 130}, {9, 13, 512}, {3, 7, 2053},
  };

  std::vector<gf2p8::Element> left;
  std::vector<gf2p8::Element> right;
  std::vector<gf2p8::Element> result;
  std::vector<gf2p8::Element> ref;

  for (const auto shape : shapes) {
    left.resize(shape.rows * shape.inner);
    right.resize(shape.inner * shape.cols);
    result.resize(shape.rows * shape.cols);
    ref.resize(shape.rows * shape.cols);

    for (auto& x : left) {
      x = rng();
    }
    for (auto& x : right) {
      x = rng();
    }

    gf2p8::MatMulBlockedLUT(left.data(), right.data(), shape.rows, shape.inner,
                            shape.cols, result.data());

    std::fill(ref.begin(), ref.end(), 0);
    for (size_t i = 0; i < shape.rows; ++i) {
      for (size_t j = 0; j < shape.cols; ++j) {
        for (size_t k = 0; k < shape.inner; ++k) {
          ref[i * shape.cols + j] ^= gf2p8::MultiplyCantor(
              left[i * shape.inner + k], right[k * shape.cols + j]);
        }
      }
    }

    ASSERT_EQ(std::equal(ref.begin(), ref.end(), result.begin()), true)
        << "shape=" << shape.rows << "x" << shape.inner << "x" << shape.cols;
  }
}

TEST(GF2P8, Irreducibly) {
  for (uint16_t x = 0; x < 256; ++x) {
    ASSERT_NE(gf2p8::Add(gf2p8::Add(gf2p8::MultiplyCantor(x, x), x), 0xf0),
              gf2p8::Zero());
  }
}

TEST(GF2P16, Distributivity) {
  std::mt19937 rng(42);
  for (size_t i = 0; i < 100000; ++i) {
    gf2p16::Element x = rng();
    gf2p16::Element y = rng();
    gf2p16::Element z = rng();
    ASSERT_EQ(gf2p16::Multiply(z, gf2p16::Add(x, y)),
              gf2p16::Add(gf2p16::Multiply(z, x), gf2p16::Multiply(z, y)));
  }
}

TEST(GF2P16, Inverse) {
  for (uint16_t x = 1; x != 0; ++x) {
    ASSERT_EQ(gf2p16::Multiply(x, gf2p16::Inv(x)), gf2p16::One());
  }
}

TEST(GF2P16, Pow_2_8) {
  for (uint16_t x = 1; x != 0; ++x) {
    ASSERT_EQ(gf2p16::Pow(x, 257) >> 8, 0);
  }
}

TEST(GF2P16, InverseIT) {
  for (uint16_t x = 1; x != 0; ++x) {
    ASSERT_EQ(gf2p16::Inv(x), gf2p16::InvIT(x));
  }
}

}  // namespace
