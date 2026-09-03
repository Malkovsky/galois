#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

#include "gtest/gtest.h"
#if defined(GF256_ENABLE_GFNI512_RADIX8_EXPERIMENT)
#include "lin_chung_han/experiment/gfni512_radix8.h"
#endif
#include "lin_chung_han/transform.h"
#include "lin_chung_han/transform_internal.h"

namespace {

using gf2p8::Element;
using gf2p8::lch::Backend;
using gf2p8::lch::Context;
using gf2p8::lch::FFT;
using gf2p8::lch::IFFT;
using gf2p8::lch::Radix;
using gf2p8::lch::Status;
using gf2p8::lch::TransformOptions;
#if defined(GF256_ENABLE_GFNI512_RADIX8_EXPERIMENT)
using Radix8Kernels = gf2p8::lch::detail::experiment::radix8::Kernels;

std::vector<Element*> Pointers(std::vector<std::vector<Element>>& shards);
std::vector<const Element*> ConstPointers(
    const std::vector<std::vector<Element>>& shards);

void ScalarIFFTButterfly(Element& x, Element& y, Element coefficient) {
  y ^= x;
  x ^= gf2p8::MultiplyCantor(y, coefficient);
}

void ScalarIFFTRadix8(Element* x0,
                      Element* x1,
                      Element* x2,
                      Element* x3,
                      Element* x4,
                      Element* x5,
                      Element* x6,
                      Element* x7,
                      size_t byte_count,
                      Element top,
                      Element middle0,
                      Element middle1,
                      Element low0,
                      Element low1,
                      Element low2,
                      Element low3,
                      const gf2p8::MultiplicationTables&) {
  for (size_t i = 0; i < byte_count; ++i) {
    ScalarIFFTButterfly(x0[i], x1[i], low0);
    ScalarIFFTButterfly(x2[i], x3[i], low1);
    ScalarIFFTButterfly(x4[i], x5[i], low2);
    ScalarIFFTButterfly(x6[i], x7[i], low3);
    ScalarIFFTButterfly(x0[i], x2[i], middle0);
    ScalarIFFTButterfly(x1[i], x3[i], middle0);
    ScalarIFFTButterfly(x4[i], x6[i], middle1);
    ScalarIFFTButterfly(x5[i], x7[i], middle1);
    ScalarIFFTButterfly(x0[i], x4[i], top);
    ScalarIFFTButterfly(x1[i], x5[i], top);
    ScalarIFFTButterfly(x2[i], x6[i], top);
    ScalarIFFTButterfly(x3[i], x7[i], top);
  }
}

void ScalarIFFTRadix8Xor(const Element* x0,
                         const Element* x1,
                         const Element* x2,
                         const Element* x3,
                         const Element* x4,
                         const Element* x5,
                         const Element* x6,
                         const Element* x7,
                         Element* output0,
                         Element* output1,
                         Element* output2,
                         Element* output3,
                         Element* output4,
                         Element* output5,
                         Element* output6,
                         Element* output7,
                         size_t byte_count,
                         Element top,
                         Element middle0,
                         Element middle1,
                         Element low0,
                         Element low1,
                         Element low2,
                         Element low3,
                         const gf2p8::MultiplicationTables& tables) {
  for (size_t i = 0; i < byte_count; ++i) {
    std::array<Element, 8> values{x0[i], x1[i], x2[i], x3[i],
                                  x4[i], x5[i], x6[i], x7[i]};
    ScalarIFFTRadix8(&values[0], &values[1], &values[2], &values[3], &values[4],
                     &values[5], &values[6], &values[7], 1, top, middle0,
                     middle1, low0, low1, low2, low3, tables);
    output0[i] ^= values[0];
    output1[i] ^= values[1];
    output2[i] ^= values[2];
    output3[i] ^= values[3];
    output4[i] ^= values[4];
    output5[i] ^= values[5];
    output6[i] ^= values[6];
    output7[i] ^= values[7];
  }
}

const Radix8Kernels kScalarRadix8Kernels{ScalarIFFTRadix8, ScalarIFFTRadix8Xor};

struct Radix8Coefficients {
  Element top;
  Element middle0;
  Element middle1;
  Element low0;
  Element low1;
  Element low2;
  Element low3;
};

Radix8Coefficients Coefficients(const Context& context,
                                size_t level,
                                size_t offset,
                                size_t block,
                                size_t distance) {
  return {
      context.Skew(level + 2, offset ^ block),
      context.Skew(level + 1, offset ^ block),
      context.Skew(level + 1, offset ^ (block + 4 * distance)),
      context.Skew(level, offset ^ block),
      context.Skew(level, offset ^ (block + 2 * distance)),
      context.Skew(level, offset ^ (block + 4 * distance)),
      context.Skew(level, offset ^ (block + 6 * distance)),
  };
}

void RunRadix8Control(std::vector<std::vector<Element>>& shards,
                      size_t byte_count,
                      const Radix8Coefficients& coefficients,
                      const gf2p8::lch::detail::ResolvedKernels& base,
                      const gf2p8::MultiplicationTables& tables) {
  auto pointers = Pointers(shards);
  base.ifft_radix4(pointers[0], pointers[1], pointers[2], pointers[3],
                   byte_count, coefficients.middle0, coefficients.low0,
                   coefficients.low1, tables);
  base.ifft_radix4(pointers[4], pointers[5], pointers[6], pointers[7],
                   byte_count, coefficients.middle1, coefficients.low2,
                   coefficients.low3, tables);
  for (size_t i = 0; i < 4; ++i) {
    base.ifft_radix2(pointers[i], pointers[i + 4], byte_count, coefficients.top,
                     tables);
  }
}

void RunRadix8Leaf(const Radix8Kernels& kernels,
                   std::vector<std::vector<Element>>& shards,
                   size_t byte_count,
                   const Radix8Coefficients& coefficients,
                   const gf2p8::MultiplicationTables& tables) {
  auto pointers = Pointers(shards);
  kernels.ifft_radix8(
      pointers[0], pointers[1], pointers[2], pointers[3], pointers[4],
      pointers[5], pointers[6], pointers[7], byte_count, coefficients.top,
      coefficients.middle0, coefficients.middle1, coefficients.low0,
      coefficients.low1, coefficients.low2, coefficients.low3, tables);
}

void RunRadix8XorLeaf(const Radix8Kernels& kernels,
                      const std::vector<std::vector<Element>>& input,
                      std::vector<std::vector<Element>>& output,
                      size_t byte_count,
                      const Radix8Coefficients& coefficients,
                      const gf2p8::MultiplicationTables& tables) {
  const auto input_pointers = ConstPointers(input);
  auto output_pointers = Pointers(output);
  kernels.ifft_radix8_xor(
      input_pointers[0], input_pointers[1], input_pointers[2],
      input_pointers[3], input_pointers[4], input_pointers[5],
      input_pointers[6], input_pointers[7], output_pointers[0],
      output_pointers[1], output_pointers[2], output_pointers[3],
      output_pointers[4], output_pointers[5], output_pointers[6],
      output_pointers[7], byte_count, coefficients.top, coefficients.middle0,
      coefficients.middle1, coefficients.low0, coefficients.low1,
      coefficients.low2, coefficients.low3, tables);
}
#endif

Element Multiply11d(Element a, Element b) {
  Element result = 0;
  while (a != 0) {
    result ^= static_cast<Element>(b * (a & 1));
    a >>= 1;
    b = static_cast<Element>((b << 1) ^ (0x1d * (b >> 7)));
  }
  return result;
}

Element CantorToPolynomial(Element value) {
  constexpr std::array<Element, 8> kCantorBasis{0x01, 0xd6, 0x98, 0x92,
                                                0x56, 0xc8, 0x58, 0xe6};
  Element result = 0;
  for (size_t bit = 0; bit < 8; ++bit) {
    if ((value & (1U << bit)) != 0) {
      result ^= kCantorBasis[bit];
    }
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

std::vector<Element*> Pointers(std::vector<std::vector<Element>>& shards) {
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

std::vector<size_t> BoundaryCounts(size_t size) {
  std::vector<size_t> result{0, 1, size / 2, size - 1, size};
  std::sort(result.begin(), result.end());
  result.erase(std::unique(result.begin(), result.end()), result.end());
  return result;
}

TEST(LCHMath, NativeCantorCoordinatesMatchPolynomialReference) {
  for (uint16_t a = 0; a < 256; ++a) {
    for (uint16_t b = 0; b < 256; ++b) {
      const auto aa = static_cast<Element>(a);
      const auto bb = static_cast<Element>(b);
      EXPECT_EQ(CantorToPolynomial(aa ^ bb),
                CantorToPolynomial(aa) ^ CantorToPolynomial(bb));
      EXPECT_EQ(CantorToPolynomial(gf2p8::MultiplyCantor(aa, bb)),
                Multiply11d(CantorToPolynomial(aa), CantorToPolynomial(bb)));
    }
  }
}

TEST(LCHMath, CantorBasisAndSubspacePolynomials) {
  const Context& context = Context::Shared();
  std::array<bool, 256> seen{};
  for (size_t i = 0; i < 256; ++i) {
    seen[context.EvaluationPoint(i)] = true;
  }
  EXPECT_TRUE(
      std::all_of(seen.begin(), seen.end(), [](bool value) { return value; }));

  for (size_t i = 0; i < Context::kFieldBits; ++i) {
    EXPECT_EQ(context.Basis()[i], static_cast<Element>(1U << i));
  }
  for (size_t i = 1; i < Context::kFieldBits; ++i) {
    const Element basis = context.Basis()[i];
    EXPECT_EQ(gf2p8::MultiplyCantor(basis, basis) ^ basis,
              context.Basis()[i - 1]);
  }

  for (size_t level = 0; level < Context::kFieldBits; ++level) {
    EXPECT_EQ(context.Skew(level, size_t{1} << level), 1);
    for (size_t i = 0; i < (size_t{1} << level); ++i) {
      EXPECT_EQ(context.Skew(level, i), 0);
    }
    for (size_t a = 0; a < 256; ++a) {
      for (size_t b = 0; b < 256; ++b) {
        EXPECT_EQ(context.Skew(level, a ^ b),
                  static_cast<Element>(context.Skew(level, a) ^
                                       context.Skew(level, b)));
      }
    }
  }
}

TEST(LCHTransform, MatchesIndependentNovelBasisEvaluation) {
  const Context& context = Context::Shared();
  std::mt19937 random(42);

  for (size_t shard_count = 1; shard_count <= 256; shard_count *= 2) {
    const size_t offset = 256 - shard_count;
    std::vector<std::vector<Element>> shards(shard_count,
                                             std::vector<Element>(1));
    std::vector<Element> coefficients(shard_count);
    for (size_t i = 0; i < shard_count; ++i) {
      coefficients[i] = static_cast<Element>(random());
      shards[i][0] = coefficients[i];
    }
    auto pointers = Pointers(shards);

    ASSERT_EQ(gf2p8::lch::FFT(context, pointers, 1, offset,
                              {Backend::scalar, Radix::radix2}),
              Status::ok);

    for (size_t evaluation = 0; evaluation < shard_count; ++evaluation) {
      const size_t point_index = offset ^ evaluation;
      Element expected = 0;
      for (size_t coefficient = 0; coefficient < shard_count; ++coefficient) {
        Element novel_basis_value = 1;
        for (size_t bit = 0; (size_t{1} << bit) < shard_count; ++bit) {
          if ((coefficient & (size_t{1} << bit)) != 0) {
            novel_basis_value = gf2p8::MultiplyCantor(
                novel_basis_value, context.Skew(bit, point_index));
          }
        }
        expected ^=
            gf2p8::MultiplyCantor(coefficients[coefficient], novel_basis_value);
      }
      EXPECT_EQ(shards[evaluation][0], expected)
          << "size=" << shard_count << " evaluation=" << evaluation;
    }
  }
}

TEST(LCHTransform, EveryBackendMatchesScalarAndRoundTrips) {
  std::mt19937 random(43);
  const size_t lengths[] = {0, 1, 15, 16, 17, 31, 32, 33, 63, 64, 65, 257};

  const Context& context = Context::Shared();
  for (size_t shard_count = 1; shard_count <= 256; shard_count *= 2) {
    for (size_t offset = 0; offset + shard_count <= 256;
         offset += shard_count) {
      for (const size_t length : lengths) {
        std::vector<std::vector<Element>> input(shard_count,
                                                std::vector<Element>(length));
        for (auto& shard : input) {
          std::generate(shard.begin(), shard.end(),
                        [&random] { return static_cast<Element>(random()); });
        }

        auto reference = input;
        auto reference_pointers = Pointers(reference);
        ASSERT_EQ(gf2p8::lch::FFT(context, reference_pointers, length, offset,
                                  {Backend::scalar, Radix::radix2}),
                  Status::ok);

        for (const Backend backend : AvailableBackends()) {
          for (const Radix radix : {Radix::radix2, Radix::radix4}) {
            auto actual = input;
            auto pointers = Pointers(actual);
            ASSERT_EQ(gf2p8::lch::FFT(context, pointers, length, offset,
                                      {backend, radix}),
                      Status::ok);
            EXPECT_EQ(actual, reference)
                << "size=" << shard_count << " length=" << length;
            ASSERT_EQ(gf2p8::lch::IFFT(context, pointers, length, offset,
                                       {backend, radix}),
                      Status::ok);
            EXPECT_EQ(actual, input)
                << "size=" << shard_count << " length=" << length;
          }
        }
      }
    }
  }
}

TEST(LCHTransform, AdvancedOverloadsMatchFullTransforms) {
  std::mt19937 random(44);
  const Context& context = Context::Shared();
  for (const size_t shard_count : {size_t{1}, size_t{2}, size_t{4}, size_t{8},
                                   size_t{16}, size_t{64}, size_t{256}}) {
    for (const size_t offset : {size_t{0}, 256 - shard_count}) {
      for (const size_t length : {size_t{0}, size_t{1}, size_t{17}, size_t{32},
                                  size_t{33}, size_t{65}}) {
        std::vector<std::vector<Element>> input(shard_count,
                                                std::vector<Element>(length));
        for (auto& shard : input) {
          std::generate(shard.begin(), shard.end(),
                        [&random] { return static_cast<Element>(random()); });
        }

        auto full_fft = input;
        auto full_fft_pointers = Pointers(full_fft);
        ASSERT_EQ(FFT(context, full_fft_pointers, length, offset,
                      {Backend::scalar, Radix::radix2}),
                  Status::ok);

        for (const Backend backend : AvailableBackends()) {
          for (const Radix radix : {Radix::radix2, Radix::radix4}) {
            const TransformOptions options{backend, radix};
            for (const size_t output_count : BoundaryCounts(shard_count)) {
              auto actual = input;
              auto pointers = Pointers(actual);
              ASSERT_EQ(
                  FFT(context, pointers, length, offset, output_count, options),
                  Status::ok);
              for (size_t i = 0; i < output_count; ++i) {
                EXPECT_EQ(actual[i], full_fft[i]);
              }
            }

            std::array<std::vector<uint8_t>, 4> masks{
                std::vector<uint8_t>(shard_count, 0),
                std::vector<uint8_t>(shard_count, 0),
                std::vector<uint8_t>(shard_count, 0),
                std::vector<uint8_t>(shard_count, 1)};
            masks[1][shard_count - 1] = 1;
            for (size_t i = 0; i < shard_count; ++i) {
              masks[2][i] = (i % 3) == 1;
            }
            for (const auto& mask : masks) {
              auto actual = input;
              auto pointers = Pointers(actual);
              ASSERT_EQ(FFT(context, pointers, length, offset, mask, options),
                        Status::ok);
              for (size_t i = 0; i < shard_count; ++i) {
                if (mask[i] != 0) {
                  EXPECT_EQ(actual[i], full_fft[i]);
                }
              }
            }

            for (const size_t input_count : BoundaryCounts(shard_count)) {
              auto padded = input;
              for (size_t i = input_count; i < shard_count; ++i) {
                std::fill(padded[i].begin(), padded[i].end(), 0);
              }
              auto expected = padded;
              auto expected_pointers = Pointers(expected);
              ASSERT_EQ(IFFT(context, expected_pointers, length, offset,
                             {Backend::scalar, Radix::radix2}),
                        Status::ok);

              auto in_place = padded;
              auto in_place_pointers = Pointers(in_place);
              ASSERT_EQ(IFFT(context, in_place_pointers, length, offset,
                             input_count, options),
                        Status::ok);
              EXPECT_EQ(in_place, expected);

              auto source_pointers = ConstPointers(padded);
              auto work = input;
              auto work_pointers = Pointers(work);
              ASSERT_EQ(IFFT(context,
                             std::span<const Element* const>(source_pointers)
                                 .first(input_count),
                             work_pointers, length, offset, options),
                        Status::ok);
              EXPECT_EQ(work, expected);

              auto accumulator = input;
              const auto initial_accumulator = accumulator;
              auto accumulator_pointers = Pointers(accumulator);
              work = input;
              work_pointers = Pointers(work);
              ASSERT_EQ(IFFT(context,
                             std::span<const Element* const>(source_pointers)
                                 .first(input_count),
                             work_pointers, accumulator_pointers, length,
                             offset, options),
                        Status::ok);
              for (size_t i = 0; i < shard_count; ++i) {
                for (size_t j = 0; j < length; ++j) {
                  EXPECT_EQ(accumulator[i][j],
                            static_cast<Element>(initial_accumulator[i][j] ^
                                                 expected[i][j]));
                }
              }
            }
          }
        }
      }
    }
  }
}

TEST(LCHKernels, ScaleAndXor4MatchScalar) {
  const Context& context = Context::Shared();
  std::mt19937 random(45);
  for (const Backend backend : AvailableBackends()) {
    for (const size_t length :
         {size_t{0}, size_t{1}, size_t{15}, size_t{16}, size_t{17}, size_t{31},
          size_t{32}, size_t{33}, size_t{63}, size_t{64}, size_t{65}}) {
      std::vector<Element> source(length);
      std::generate(source.begin(), source.end(),
                    [&random] { return static_cast<Element>(random()); });
      for (const Element coefficient :
           {Element{0}, Element{1}, Element{0x80}}) {
        std::vector<Element> actual(length, 0xa5);
        ASSERT_EQ(gf2p8::lch::Scale(actual.data(), source.data(), length,
                                    coefficient, backend, context.Tables()),
                  Status::ok);
        for (size_t i = 0; i < length; ++i) {
          EXPECT_EQ(actual[i], gf2p8::MultiplyCantor(source[i], coefficient));
        }
      }

      std::array<std::vector<Element>, 4> destination;
      std::array<std::vector<Element>, 4> xor_source;
      std::array<std::vector<Element>, 4> expected;
      std::array<Element*, 4> destination_pointers{};
      std::array<const Element*, 4> source_pointers{};
      for (size_t stream = 0; stream < 4; ++stream) {
        destination[stream].resize(length);
        xor_source[stream].resize(length);
        std::generate(destination[stream].begin(), destination[stream].end(),
                      [&random] { return static_cast<Element>(random()); });
        std::generate(xor_source[stream].begin(), xor_source[stream].end(),
                      [&random] { return static_cast<Element>(random()); });
        expected[stream] = destination[stream];
        for (size_t i = 0; i < length; ++i) {
          expected[stream][i] ^= xor_source[stream][i];
        }
        destination_pointers[stream] = destination[stream].data();
        source_pointers[stream] = xor_source[stream].data();
      }
      ASSERT_EQ(gf2p8::lch::Xor4(destination_pointers[0], source_pointers[0],
                                 destination_pointers[1], source_pointers[1],
                                 destination_pointers[2], source_pointers[2],
                                 destination_pointers[3], source_pointers[3],
                                 length, backend),
                Status::ok);
      EXPECT_EQ(destination, expected);
    }
  }
}

TEST(LCHKernels, ResolvedTablesCoverEveryBackend) {
  EXPECT_EQ(gf2p8::lch::detail::ResolveKernels(Backend::tuned), nullptr);
  for (const Backend backend :
       {Backend::scalar, Backend::ssse3, Backend::avx2, Backend::gfni128_affine,
        Backend::gfni256_affine, Backend::gfni512_affine}) {
    const auto* kernels = gf2p8::lch::detail::ResolveKernels(backend);
    EXPECT_EQ(kernels != nullptr, gf2p8::lch::BackendAvailable(backend));
    if (kernels != nullptr) {
      const bool has_copy_first = backend == Backend::avx2 ||
                                  backend == Backend::gfni128_affine ||
                                  backend == Backend::gfni256_affine ||
                                  backend == Backend::gfni512_affine;
      EXPECT_NE(kernels->add_scaled, nullptr);
      EXPECT_NE(kernels->scale, nullptr);
      EXPECT_NE(kernels->xor_one, nullptr);
      EXPECT_NE(kernels->xor_four, nullptr);
      EXPECT_EQ(kernels->ifft_radix2_copy != nullptr, has_copy_first);
      EXPECT_EQ(kernels->ifft_radix4_copy != nullptr, has_copy_first);
      EXPECT_NE(kernels->fft_radix2, nullptr);
      EXPECT_NE(kernels->ifft_radix2, nullptr);
      EXPECT_NE(kernels->ifft_radix2_xor, nullptr);
      EXPECT_NE(kernels->fft_radix4, nullptr);
      EXPECT_NE(kernels->ifft_radix4, nullptr);
      EXPECT_NE(kernels->ifft_radix4_xor, nullptr);
    }
  }
  const auto* avx2_exact =
      gf2p8::lch::detail::ResolveKernels(Backend::avx2, 32);
  const auto* avx2_tailed =
      gf2p8::lch::detail::ResolveKernels(Backend::avx2, 33);
  if (gf2p8::lch::BackendAvailable(Backend::avx2)) {
    EXPECT_NE(avx2_exact, avx2_tailed);
  } else {
    EXPECT_EQ(avx2_exact, nullptr);
    EXPECT_EQ(avx2_tailed, nullptr);
  }
}

TEST(LCHKernels, TunedBackendUsesFinalByteThresholds) {
  const Backend at_16 = gf2p8::lch::BackendAvailable(Backend::ssse3)
                            ? Backend::ssse3
                            : Backend::scalar;
  const Backend at_32 =
      gf2p8::lch::BackendAvailable(Backend::avx2) ? Backend::avx2 : at_16;
  const Backend at_128 = gf2p8::lch::BackendAvailable(Backend::gfni256_affine)
                             ? Backend::gfni256_affine
                             : at_32;

  EXPECT_EQ(gf2p8::lch::SelectBackend(0), Backend::scalar);
  EXPECT_EQ(gf2p8::lch::SelectBackend(15), Backend::scalar);
  EXPECT_EQ(gf2p8::lch::SelectBackend(16), at_16);
  EXPECT_EQ(gf2p8::lch::SelectBackend(31), at_16);
  EXPECT_EQ(gf2p8::lch::SelectBackend(32), at_32);
  EXPECT_EQ(gf2p8::lch::SelectBackend(127), at_32);
  EXPECT_EQ(gf2p8::lch::SelectBackend(128), at_128);
  EXPECT_EQ(gf2p8::lch::SelectBackend(std::numeric_limits<size_t>::max()),
            at_128);
}

#if defined(GF256_ENABLE_GFNI512_RADIX8_EXPERIMENT)
TEST(LCHRadix8Experiment, ResolverMatchesGFNI512Gate) {
  const auto* kernels =
      gf2p8::lch::detail::experiment::radix8::ResolveKernels();
  EXPECT_EQ(kernels != nullptr,
            gf2p8::lch::BackendAvailable(Backend::gfni512_affine));
  if (kernels != nullptr) {
    EXPECT_NE(kernels->ifft_radix8, nullptr);
    EXPECT_NE(kernels->ifft_radix8_xor, nullptr);
  }
}

TEST(LCHRadix8Experiment, LeafMatchesExactRadix4AndRadix2Decomposition) {
  const Context& context = Context::Shared();
  const auto& base =
      *gf2p8::lch::detail::ResolveKernels(Backend::scalar, size_t{1});
  std::vector<const Radix8Kernels*> candidates{&kScalarRadix8Kernels};
  if (const auto* handwritten =
          gf2p8::lch::detail::experiment::radix8::ResolveKernels()) {
    candidates.push_back(handwritten);
  }
  std::mt19937 random(46);
  for (const size_t length :
       {size_t{0}, size_t{1}, size_t{63}, size_t{64}, size_t{65}, size_t{127},
        size_t{128}, size_t{129}}) {
    for (const size_t offset : {size_t{0}, size_t{8}}) {
      const Radix8Coefficients coefficients =
          Coefficients(context, 0, offset, 0, 1);
      if (offset == 0) {
        EXPECT_EQ(coefficients.low0, 0);
        EXPECT_EQ(coefficients.middle0, 0);
        EXPECT_EQ(coefficients.top, 0);
      }
      std::vector<std::vector<Element>> input(8, std::vector<Element>(length));
      for (auto& shard : input) {
        std::generate(shard.begin(), shard.end(),
                      [&random] { return static_cast<Element>(random()); });
      }
      auto expected = input;
      RunRadix8Control(expected, length, coefficients, base, context.Tables());

      for (const Radix8Kernels* candidate : candidates) {
        auto actual = input;
        RunRadix8Leaf(*candidate, actual, length, coefficients,
                      context.Tables());
        EXPECT_EQ(actual, expected)
            << "length=" << length << " offset=" << offset;

        std::vector<std::vector<Element>> accumulator(
            8, std::vector<Element>(length));
        for (auto& shard : accumulator) {
          std::generate(shard.begin(), shard.end(),
                        [&random] { return static_cast<Element>(random()); });
        }
        auto expected_accumulator = accumulator;
        for (size_t shard = 0; shard < 8; ++shard) {
          for (size_t byte = 0; byte < length; ++byte) {
            expected_accumulator[shard][byte] ^= expected[shard][byte];
          }
        }
        const auto immutable_input = input;
        RunRadix8XorLeaf(*candidate, input, accumulator, length, coefficients,
                         context.Tables());
        EXPECT_EQ(accumulator, expected_accumulator)
            << "length=" << length << " offset=" << offset;
        EXPECT_EQ(input, immutable_input);
      }
    }
  }
}

TEST(LCHRadix8Experiment, FullSchedulerMatchesCurrentRadix4) {
  const Context& context = Context::Shared();
  const auto& base =
      *gf2p8::lch::detail::ResolveKernels(Backend::scalar, size_t{65});
  std::vector<const Radix8Kernels*> candidates{&kScalarRadix8Kernels};
  const auto* handwritten =
      gf2p8::lch::detail::experiment::radix8::ResolveKernels();
  if (handwritten != nullptr) {
    candidates.push_back(handwritten);
  }
  std::mt19937 random(47);
  for (const size_t shard_count : {size_t{8}, size_t{16}, size_t{32},
                                   size_t{64}, size_t{128}, size_t{256}}) {
    for (size_t offset = 0; offset + shard_count <= Context::kFieldSize;
         offset += shard_count) {
      std::vector<std::vector<Element>> source(shard_count,
                                               std::vector<Element>(65));
      for (auto& shard : source) {
        std::generate(shard.begin(), shard.end(),
                      [&random] { return static_cast<Element>(random()); });
      }
      for (const size_t input_count : BoundaryCounts(shard_count)) {
        auto padded = source;
        for (size_t i = input_count; i < shard_count; ++i) {
          std::fill(padded[i].begin(), padded[i].end(), 0);
        }
        auto expected = padded;
        auto expected_pointers = Pointers(expected);
        ASSERT_EQ(gf2p8::lch::detail::IFFTResolved(context, expected_pointers,
                                                   65, offset, input_count,
                                                   base, Radix::radix4),
                  Status::ok);

        for (const Radix8Kernels* candidate : candidates) {
          const auto& candidate_base =
              candidate == handwritten
                  ? *gf2p8::lch::detail::ResolveKernels(Backend::gfni512_affine,
                                                        size_t{65})
                  : base;
          auto in_place = padded;
          auto in_place_pointers = Pointers(in_place);
          ASSERT_EQ(gf2p8::lch::detail::experiment::radix8::IFFTResolved(
                        context, in_place_pointers, 65, offset, input_count,
                        candidate_base, *candidate),
                    Status::ok);
          EXPECT_EQ(in_place, expected)
              << "N=" << shard_count << " offset=" << offset
              << " input=" << input_count;

          const auto immutable_source = source;
          const auto source_pointers = ConstPointers(source);
          auto work = source;
          auto work_pointers = Pointers(work);
          ASSERT_EQ(gf2p8::lch::detail::experiment::radix8::IFFTResolved(
                        context,
                        std::span<const Element* const>(source_pointers)
                            .first(input_count),
                        work_pointers, 65, offset, candidate_base, *candidate),
                    Status::ok);
          EXPECT_EQ(work, expected)
              << "N=" << shard_count << " offset=" << offset
              << " input=" << input_count;
          EXPECT_EQ(source, immutable_source);

          std::vector<std::vector<Element>> accumulator(
              shard_count, std::vector<Element>(65));
          for (auto& shard : accumulator) {
            std::generate(shard.begin(), shard.end(),
                          [&random] { return static_cast<Element>(random()); });
          }
          auto expected_accumulator = accumulator;
          for (size_t shard = 0; shard < shard_count; ++shard) {
            for (size_t byte = 0; byte < 65; ++byte) {
              expected_accumulator[shard][byte] ^= expected[shard][byte];
            }
          }
          work = source;
          work_pointers = Pointers(work);
          auto accumulator_pointers = Pointers(accumulator);
          ASSERT_EQ(gf2p8::lch::detail::experiment::radix8::IFFTResolved(
                        context,
                        std::span<const Element* const>(source_pointers)
                            .first(input_count),
                        work_pointers, accumulator_pointers, 65, offset,
                        candidate_base, *candidate),
                    Status::ok);
          EXPECT_EQ(accumulator, expected_accumulator)
              << "N=" << shard_count << " offset=" << offset
              << " input=" << input_count;
          EXPECT_EQ(source, immutable_source);
        }
      }
    }
  }
}

TEST(LCHRadix8Experiment, ZeroBytesTouchNoPointers) {
  const Context& context = Context::Shared();
  const auto& base =
      *gf2p8::lch::detail::ResolveKernels(Backend::scalar, size_t{0});
  std::array<Element*, 8> null_work{};
  std::array<Element*, 8> null_accumulator{};
  std::array<const Element*, 8> null_input{};
  EXPECT_EQ(gf2p8::lch::detail::experiment::radix8::IFFTResolved(
                context, std::span<Element* const>(null_work), 0, 0, 8, base,
                kScalarRadix8Kernels),
            Status::ok);
  EXPECT_EQ(gf2p8::lch::detail::experiment::radix8::IFFTResolved(
                context, std::span<const Element* const>(null_input),
                std::span<Element* const>(null_work), 0, 0, base,
                kScalarRadix8Kernels),
            Status::ok);
  EXPECT_EQ(gf2p8::lch::detail::experiment::radix8::IFFTResolved(
                context, std::span<const Element* const>(null_input),
                std::span<Element* const>(null_work),
                std::span<Element* const>(null_accumulator), 0, 0, base,
                kScalarRadix8Kernels),
            Status::ok);
}
#endif

TEST(LCHKernels, CopyFirstPropagatesResolvedKernelFailure) {
  auto kernels =
      *gf2p8::lch::detail::ResolveKernels(Backend::scalar, size_t{1});
  kernels.ifft_radix2_copy =
      +[](const Element*, const Element*, Element*, Element*, size_t, Element,
          const gf2p8::MultiplicationTables&) {
        return Status::unsupported_backend;
      };
  kernels.ifft_radix4_copy =
      +[](const Element*, const Element*, const Element*, const Element*,
          Element*, Element*, Element*, Element*, size_t, Element, Element,
          Element, const gf2p8::MultiplicationTables&) {
        return Status::unsupported_backend;
      };

  std::array<std::array<Element, 1>, 4> input_storage{};
  std::array<std::array<Element, 1>, 4> work_storage{};
  std::array<const Element*, 4> input{};
  std::array<Element*, 4> work{};
  for (size_t i = 0; i < input.size(); ++i) {
    input[i] = input_storage[i].data();
    work[i] = work_storage[i].data();
  }

  const Context& context = Context::Shared();
  EXPECT_EQ(gf2p8::lch::detail::IFFTResolved(
                context, std::span<const Element* const>(input).first(2),
                std::span<Element* const>(work).first(2), 1, 0, kernels,
                Radix::radix2),
            Status::unsupported_backend);
  EXPECT_EQ(gf2p8::lch::detail::IFFTResolved(
                context, std::span<const Element* const>(input),
                std::span<Element* const>(work), 1, 0, kernels, Radix::radix4),
            Status::unsupported_backend);
}

TEST(LCHSchedule, MatchesLeopardPairedAndSparseClosures) {
  const auto full = gf2p8::lch::testing::FFTSchedule(256, 256, Radix::radix4);
  EXPECT_EQ(full.radix4_groups, 85);
  EXPECT_EQ(full.radix4_kernel_calls, 256);
  EXPECT_EQ(full.radix2_groups, 0);
  EXPECT_EQ(full.shard_kernel_visits, 1024);

  std::array<uint8_t, 64> requested64{};
  requested64[16] = 1;
  const auto sparse64 =
      gf2p8::lch::testing::FFTSchedule(64, requested64, Radix::radix4);
  EXPECT_EQ(sparse64.radix4_groups, 3);
  EXPECT_EQ(sparse64.radix4_kernel_calls, 21);
  EXPECT_EQ(sparse64.shard_kernel_visits, 84);

  std::array<uint8_t, 256> requested256{};
  requested256[64] = 1;
  const auto sparse256 =
      gf2p8::lch::testing::FFTSchedule(256, requested256, Radix::radix4);
  EXPECT_EQ(sparse256.radix4_groups, 4);
  EXPECT_EQ(sparse256.radix4_kernel_calls, 85);
  EXPECT_EQ(sparse256.shard_kernel_visits, 340);

  const auto partial_ifft =
      gf2p8::lch::testing::IFFTSchedule(64, 1, Radix::radix4, true);
  EXPECT_EQ(partial_ifft.radix4_groups, 3);
  EXPECT_EQ(partial_ifft.radix4_kernel_calls, 21);
  EXPECT_EQ(partial_ifft.fused_accumulation_groups, 1);
  EXPECT_EQ(partial_ifft.fused_accumulation_kernel_calls, 16);
  EXPECT_EQ(partial_ifft.temporary_shard_stores_avoided, 64);

  const auto full_ifft =
      gf2p8::lch::testing::IFFTSchedule(16, 16, Radix::radix4, true);
  EXPECT_EQ(full_ifft.radix4_groups, 5);
  EXPECT_EQ(full_ifft.radix4_kernel_calls, 8);
  EXPECT_EQ(full_ifft.fused_accumulation_groups, 1);
  EXPECT_EQ(full_ifft.fused_accumulation_kernel_calls, 4);
  EXPECT_EQ(full_ifft.temporary_shard_stores_avoided, 16);

  const auto radix2 = gf2p8::lch::testing::FFTSchedule(256, 256, Radix::radix2);
  EXPECT_EQ(radix2.radix2_groups, 255);
  EXPECT_EQ(radix2.radix2_kernel_calls, 1024);
  EXPECT_EQ(radix2.shard_kernel_visits, 2048);

  struct LeopardShape {
    size_t k;
    size_t r;
    size_t m;
    size_t n;
    size_t source_blocks;
    size_t encode_shard_visits;
    size_t fused_kernel_calls;
    size_t decode_shard_visits;
  };
  constexpr LeopardShape shapes[] = {
      {32, 16, 16, 64, 2, 96, 4, 244},
      {129, 64, 64, 256, 3, 660, 32, 1256},
      {128, 128, 128, 256, 1, 1024, 0, 1364},
  };
  for (const auto& shape : shapes) {
    size_t encode_visits = 0;
    size_t fused_calls = 0;
    size_t source_offset = 0;
    size_t source_blocks = 0;
    while (source_offset < shape.k) {
      const size_t live = std::min(shape.m, shape.k - source_offset);
      const bool fused = source_blocks != 0;
      const auto schedule = gf2p8::lch::testing::IFFTSchedule(
          shape.m, live, Radix::radix4, fused);
      encode_visits += schedule.shard_kernel_visits;
      fused_calls += schedule.fused_accumulation_kernel_calls;
      source_offset += live;
      ++source_blocks;
    }
    const auto parity =
        gf2p8::lch::testing::FFTSchedule(shape.m, shape.r, Radix::radix4);
    encode_visits += parity.shard_kernel_visits;

    const auto decode_ifft = gf2p8::lch::testing::IFFTSchedule(
        shape.n, shape.m + shape.k, Radix::radix4, false);
    std::array<uint8_t, 256> requested{};
    requested[shape.m + shape.k / 2] = 1;
    const auto sparse = gf2p8::lch::testing::FFTSchedule(
        shape.n, std::span<const uint8_t>(requested).first(shape.n),
        Radix::radix4);

    EXPECT_EQ(source_blocks, shape.source_blocks) << shape.k << '/' << shape.r;
    EXPECT_EQ(encode_visits, shape.encode_shard_visits)
        << shape.k << '/' << shape.r;
    EXPECT_EQ(fused_calls, shape.fused_kernel_calls)
        << shape.k << '/' << shape.r;
    EXPECT_EQ(decode_ifft.shard_kernel_visits + sparse.shard_kernel_visits,
              shape.decode_shard_visits)
        << shape.k << '/' << shape.r;
  }
}

TEST(LCHTransform, ValidatesContract) {
  const Context& context = Context::Shared();
  std::array<Element, 1> byte{};
  std::array<Element*, 1> one{byte.data()};
  std::array<Element*, 2> two{byte.data(), nullptr};
  std::array<Element*, 3> three{byte.data(), byte.data(), byte.data()};

  EXPECT_EQ(gf2p8::lch::FFT(context, {}, 0), Status::invalid_argument);
  EXPECT_EQ(gf2p8::lch::FFT(context, three, 1), Status::invalid_argument);
  EXPECT_EQ(gf2p8::lch::FFT(context, two, 1), Status::invalid_argument);
  EXPECT_EQ(gf2p8::lch::FFT(context, one, 0, 255), Status::ok);
  EXPECT_EQ(gf2p8::lch::FFT(context, one, 1, 0,
                            {Backend::scalar, static_cast<Radix>(255)}),
            Status::invalid_argument);

  const Backend unavailable[] = {
      Backend::ssse3,          Backend::avx2,           Backend::gfni128_affine,
      Backend::gfni256_affine, Backend::gfni512_affine,
  };
  for (const Backend backend : unavailable) {
    if (!gf2p8::lch::BackendAvailable(backend)) {
      EXPECT_EQ(gf2p8::lch::FFT(context, one, 1, 0, {backend, Radix::radix2}),
                Status::unsupported_backend);
      EXPECT_EQ(gf2p8::lch::Xor(byte.data(), byte.data(), 1, backend),
                Status::unsupported_backend);
      EXPECT_EQ(gf2p8::lch::FFTRadix2(byte.data(), byte.data(), 1, 1, backend,
                                      context.Tables()),
                Status::unsupported_backend);
      EXPECT_EQ(gf2p8::lch::IFFTRadix2(byte.data(), byte.data(), 1, 1, backend,
                                       context.Tables()),
                Status::unsupported_backend);
      EXPECT_EQ(gf2p8::lch::FFTRadix4(byte.data(), byte.data(), byte.data(),
                                      byte.data(), 1, 1, 1, 1, backend,
                                      context.Tables()),
                Status::unsupported_backend);
      EXPECT_EQ(gf2p8::lch::IFFTRadix4(byte.data(), byte.data(), byte.data(),
                                       byte.data(), 1, 1, 1, 1, backend,
                                       context.Tables()),
                Status::unsupported_backend);
    }
  }
}

}  // namespace
