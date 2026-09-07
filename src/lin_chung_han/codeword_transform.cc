#if defined(__i386__) || defined(__x86_64__) || defined(_M_IX86) || \
    defined(_M_X64)
#include <immintrin.h>
#endif

#include <array>
#include <bit>
#include <cstddef>
#include <span>

#include "field.h"
#include "lin_chung_han/codeword_transform_internal.h"
#include "lin_chung_han/kernels_internal.h"

namespace gf2p8::lch::detail {
namespace {

constexpr Element AESMultiply(Element a, Element b) {
  Element result = 0;
  while (a != 0) {
    result ^= static_cast<Element>(b * (a & 1U));
    a >>= 1;
    b = static_cast<Element>((b << 1) ^ (0x1bU * (b >> 7)));
  }
  return result;
}

constexpr Element Evaluate11d(Element value) {
  Element result = 1;
  for (int bit = 7; bit >= 0; --bit) {
    result = AESMultiply(result, value);
    if (((0x11dU >> bit) & 1U) != 0) {
      result ^= 1;
    }
  }
  return result;
}

static_assert(Evaluate11d(0x03) == 0);

constexpr Element Standard11dToAES(Element value) {
  Element result = 0;
  Element power = 1;
  for (size_t bit = 0; bit < 8; ++bit) {
    if (((value >> bit) & 1U) != 0) {
      result ^= power;
    }
    power = AESMultiply(power, 0x03);
  }
  return result;
}

constexpr Element CantorToAES(Element value) {
  return Standard11dToAES(gf2p8::detail::CantorToStandardDirect(value));
}

template <typename Map>
consteval uint64_t AffineMatrix(Map map) {
  uint64_t matrix = 0;
  for (size_t output_bit = 0; output_bit < 8; ++output_bit) {
    for (size_t input_bit = 0; input_bit < 8; ++input_bit) {
      const uint64_t bit =
          (map(static_cast<Element>(1U << input_bit)) >> output_bit) & 1U;
      matrix |= bit << (8 * (7 - output_bit) + input_bit);
    }
  }
  return matrix;
}

struct IsomorphismTables {
  std::array<Element, 256> cantor_to_aes{};
  std::array<Element, 256> aes_to_cantor{};
  uint64_t cantor_to_aes_matrix = 0;
  uint64_t aes_to_cantor_matrix = 0;
};

consteval IsomorphismTables MakeIsomorphismTables() {
  IsomorphismTables tables;
  for (size_t value = 0; value < 256; ++value) {
    tables.cantor_to_aes[value] = CantorToAES(static_cast<Element>(value));
  }
  for (size_t value = 0; value < 256; ++value) {
    tables.aes_to_cantor[tables.cantor_to_aes[value]] =
        static_cast<Element>(value);
  }
  tables.cantor_to_aes_matrix =
      AffineMatrix([](Element value) { return CantorToAES(value); });
  tables.aes_to_cantor_matrix = AffineMatrix(
      [&tables](Element value) { return tables.aes_to_cantor[value]; });
  return tables;
}

inline constexpr IsomorphismTables kIsomorphism = MakeIsomorphismTables();
static_assert(kIsomorphism.cantor_to_aes_matrix == 0xefd0aaca822e5caeULL);
static_assert(kIsomorphism.aes_to_cantor_matrix == 0xffb08466dc727ea0ULL);

consteval bool IsomorphismIsValid() {
  for (size_t value = 0; value < 256; ++value) {
    if (kIsomorphism.aes_to_cantor[kIsomorphism.cantor_to_aes[value]] !=
        value) {
      return false;
    }
  }
  // Multiplication is bilinear, so checking all basis-vector pairs proves the
  // generated linear map is a field isomorphism for every byte pair.
  for (size_t a_bit = 0; a_bit < 8; ++a_bit) {
    const Element a = static_cast<Element>(1U << a_bit);
    for (size_t b_bit = 0; b_bit < 8; ++b_bit) {
      const Element b = static_cast<Element>(1U << b_bit);
      const Element cantor_product = gf2p8::detail::MultiplyCantorDirect(a, b);
      if (AESMultiply(kIsomorphism.cantor_to_aes[a],
                      kIsomorphism.cantor_to_aes[b]) !=
          kIsomorphism.cantor_to_aes[cantor_product]) {
        return false;
      }
    }
  }
  return true;
}

static_assert(IsomorphismIsValid());

Status Validate(std::span<Element> values,
                size_t block_size,
                size_t evaluation_offset,
                Backend backend) {
  if (values.empty() || !std::has_single_bit(values.size()) ||
      !std::has_single_bit(block_size) || values.size() % block_size != 0 ||
      values.size() > Context::kFieldSize ||
      evaluation_offset >= Context::kFieldSize ||
      evaluation_offset + values.size() > Context::kFieldSize ||
      (evaluation_offset & (block_size - 1)) != 0) {
    return Status::invalid_argument;
  }
  if (backend != Backend::scalar && backend != Backend::tuned &&
      backend != Backend::avx2 && backend != Backend::gfni256_affine) {
    return Status::unsupported_backend;
  }
  if ((backend == Backend::avx2 || backend == Backend::gfni256_affine) &&
      !BackendAvailable(backend)) {
    return Status::unsupported_backend;
  }
  return Status::ok;
}

size_t Log2(size_t value) {
  return std::bit_width(value) - 1;
}

Element Product(Element value,
                Element coefficient,
                const MultiplicationTables& tables) {
  if (coefficient == 0) {
    return 0;
  }
  const auto& row = tables.shuffle[coefficient];
  return row[value & 0x0f] ^ row[32 + (value >> 4)];
}

void FFTScalar(const Context& context,
               std::span<Element> values,
               size_t block_size,
               size_t evaluation_offset) {
  const MultiplicationTables& tables = context.Tables();
  for (size_t half = block_size / 2; half != 0; half /= 2) {
    const size_t group_size = 2 * half;
    const size_t level = Log2(half);
    for (size_t outer = 0; outer < values.size(); outer += block_size) {
      const size_t transform_offset = evaluation_offset + outer;
      for (size_t group = 0; group < block_size; group += group_size) {
        const Element coefficient =
            context.Skew(level, transform_offset ^ group);
        for (size_t i = 0; i < half; ++i) {
          Element& x = values[outer + group + i];
          Element& y = values[outer + group + half + i];
          x ^= Product(y, coefficient, tables);
          y ^= x;
        }
      }
    }
  }
}

void IFFTScalar(const Context& context,
                std::span<Element> values,
                size_t block_size,
                size_t evaluation_offset) {
  const MultiplicationTables& tables = context.Tables();
  for (size_t half = 1; half < block_size; half *= 2) {
    const size_t group_size = 2 * half;
    const size_t level = Log2(half);
    for (size_t outer = 0; outer < values.size(); outer += block_size) {
      const size_t transform_offset = evaluation_offset + outer;
      for (size_t group = 0; group < block_size; group += group_size) {
        const Element coefficient =
            context.Skew(level, transform_offset ^ group);
        for (size_t i = 0; i < half; ++i) {
          Element& x = values[outer + group + i];
          Element& y = values[outer + group + half + i];
          y ^= x;
          x ^= Product(y, coefficient, tables);
        }
      }
    }
  }
}

#if defined(__AVX2__)
template <bool Inverse>
void TransformAVX2(const Context& context,
                   std::span<Element> values,
                   size_t block_size,
                   size_t evaluation_offset) {
  const ResolvedKernels& kernels = *ResolveKernels(Backend::avx2, size_t{32});
  const ResolvedKernels* short_kernels =
      ResolveKernels(Backend::ssse3, size_t{16});
  const MultiplicationTables& tables = context.Tables();
  for (size_t half = Inverse ? 1 : block_size / 2;
       Inverse ? half < block_size : half != 0;
       half = Inverse ? half * 2 : half / 2) {
    const size_t group_size = 2 * half;
    const size_t level = Log2(half);
    for (size_t outer = 0; outer < values.size(); outer += block_size) {
      const size_t transform_offset = evaluation_offset + outer;
      for (size_t group = 0; group < block_size; group += group_size) {
        const Element coefficient =
            context.Skew(level, transform_offset ^ group);
        Element* x = values.data() + outer + group;
        Element* y = x + half;
        const ResolvedKernels* stage_kernels = half >= 32   ? &kernels
                                               : half == 16 ? short_kernels
                                                            : nullptr;
        if (stage_kernels != nullptr) {
          if constexpr (Inverse) {
            stage_kernels->ifft_radix2(x, y, half, coefficient, tables);
          } else {
            stage_kernels->fft_radix2(x, y, half, coefficient, tables);
          }
          continue;
        }
        for (size_t i = 0; i < half; ++i) {
          if constexpr (Inverse) {
            y[i] ^= x[i];
            x[i] ^= Product(y[i], coefficient, tables);
          } else {
            x[i] ^= Product(y[i], coefficient, tables);
            y[i] ^= x[i];
          }
        }
      }
    }
  }
}
#endif

#if defined(__GFNI__) && defined(__AVX2__)

void ConvertBasis(std::span<Element> values,
                  uint64_t matrix_value,
                  const std::array<Element, 256>& scalar_table) {
  const __m256i matrix =
      _mm256_set1_epi64x(static_cast<long long>(matrix_value));
  size_t i = 0;
  for (; i + 32 <= values.size(); i += 32) {
    const __m256i input =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(values.data() + i));
    const __m256i output = _mm256_gf2p8affine_epi64_epi8(input, matrix, 0);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(values.data() + i), output);
  }
  for (; i < values.size(); ++i) {
    values[i] = scalar_table[values[i]];
  }
}

struct alignas(32) SmallStageMasks {
  std::array<std::array<Element, 32>, 4> low{};
  std::array<std::array<Element, 32>, 4> high{};
  std::array<std::array<Element, 32>, 4> upper{};
#if defined(GF256_ENABLE_CODEWORD_CANTOR_AFFINE_EXPERIMENT)
  std::array<std::array<std::array<Element, 32>, 4>, 4> affine_pass{};
#endif
};

consteval SmallStageMasks MakeSmallStageMasks() {
  SmallStageMasks masks;
  for (size_t level = 0; level < 4; ++level) {
    const size_t half = size_t{1} << level;
    const size_t group_size = 2 * half;
    for (size_t index = 0; index < 32; ++index) {
      const size_t lane_index = index & 15U;
      const size_t group = lane_index & ~(group_size - 1);
      const size_t low = group + (lane_index & (half - 1));
      masks.low[level][index] = static_cast<Element>(low);
      masks.high[level][index] = static_cast<Element>(low + half);
      masks.upper[level][index] = (lane_index & half) != 0 ? 0xff : 0;

#if defined(GF256_ENABLE_CODEWORD_CANTOR_AFFINE_EXPERIMENT)
      const size_t pass_count = half < 4 ? 4 / half : 1;
      for (size_t pass = 0; pass < pass_count; ++pass) {
        bool selected = false;
        if (half == 8) {
          selected = (lane_index & 15U) < 8;
        } else {
          const size_t qword_index = index & 7U;
          const size_t first = pass * group_size;
          selected = qword_index >= first && qword_index < first + half;
        }
        masks.affine_pass[level][pass][index] = selected ? 0xff : 0;
      }
#endif
    }
  }
  return masks;
}

inline constexpr SmallStageMasks kSmallStageMasks = MakeSmallStageMasks();

struct alignas(32) AESSkewTables {
  std::array<std::array<Element, Context::kFieldSize>, Context::kFieldBits>
      factors{};
};

const AESSkewTables& GetAESSkewTables() {
  static const AESSkewTables tables = [] {
    AESSkewTables result;
    const Context& context = Context::Shared();
    for (size_t level = 0; level < Context::kFieldBits; ++level) {
      const size_t group_size = size_t{2} << level;
      for (size_t position = 0; position < Context::kFieldSize; ++position) {
        const size_t group = position & ~(group_size - 1);
        result.factors[level][position] =
            kIsomorphism.cantor_to_aes[context.Skew(level, group)];
      }
    }
    return result;
  }();
  return tables;
}

#if defined(GF256_ENABLE_CODEWORD_CANTOR_AFFINE_EXPERIMENT)
struct alignas(32) CantorAffineSkewTables {
  using PackedMatrices = std::array<uint64_t, 4>;

  std::array<std::array<uint64_t, Context::kFieldSize>, Context::kFieldBits>
      factors{};
  std::array<std::array<std::array<PackedMatrices, 8>, 4>, 4> small_aligned{};
};

const CantorAffineSkewTables& GetCantorAffineSkewTables() {
  static const CantorAffineSkewTables tables = [] {
    CantorAffineSkewTables result;
    const Context& context = Context::Shared();
    const MultiplicationTables& multiplication = context.Tables();
    for (size_t level = 0; level < Context::kFieldBits; ++level) {
      const size_t group_size = size_t{2} << level;
      for (size_t position = 0; position < Context::kFieldSize; ++position) {
        const size_t group = position & ~(group_size - 1);
        result.factors[level][position] =
            multiplication.affine[context.Skew(level, group)];
      }
    }
    for (size_t level = 0; level < 4; ++level) {
      const size_t half = size_t{1} << level;
      const size_t pass_count = half < 4 ? 4 / half : 1;
      for (size_t pass = 0; pass < pass_count; ++pass) {
        for (size_t chunk = 0; chunk < 8; ++chunk) {
          const size_t base = 32 * chunk;
          auto& packed = result.small_aligned[level][pass][chunk];
          if (half == 8) {
            packed[0] = result.factors[level][base];
            packed[2] = result.factors[level][base + 16];
          } else {
            const size_t group_size = 2 * half;
            for (size_t qword = 0; qword < packed.size(); ++qword) {
              packed[qword] =
                  result.factors[level][base + 8 * qword + pass * group_size];
            }
          }
        }
      }
    }
    return result;
  }();
  return tables;
}
#endif

template <bool Inverse>
void Butterfly(__m256i& x, __m256i& y, __m256i coefficient) {
  if constexpr (Inverse) {
    y = _mm256_xor_si256(y, x);
    x = _mm256_xor_si256(x, _mm256_gf2p8mul_epi8(y, coefficient));
  } else {
    x = _mm256_xor_si256(x, _mm256_gf2p8mul_epi8(y, coefficient));
    y = _mm256_xor_si256(y, x);
  }
}

template <bool Inverse>
void TransformGFNIImpl(std::span<Element> values,
                       size_t block_size,
                       size_t evaluation_offset) {
  if (block_size == 1 || values.size() < 32) {
    if constexpr (Inverse) {
      IFFTScalar(Context::Shared(), values, block_size, evaluation_offset);
    } else {
      FFTScalar(Context::Shared(), values, block_size, evaluation_offset);
    }
    return;
  }

  ConvertBasis(values, kIsomorphism.cantor_to_aes_matrix,
               kIsomorphism.cantor_to_aes);
  const AESSkewTables& skew = GetAESSkewTables();

  for (size_t half = Inverse ? 1 : block_size / 2;
       Inverse ? half < block_size : half != 0;
       half = Inverse ? half * 2 : half / 2) {
    const size_t group_size = 2 * half;
    const size_t level = Log2(half);
    if (half <= 8) {
      const __m256i low_mask = _mm256_load_si256(
          reinterpret_cast<const __m256i*>(kSmallStageMasks.low[level].data()));
      const __m256i high_mask =
          _mm256_load_si256(reinterpret_cast<const __m256i*>(
              kSmallStageMasks.high[level].data()));
      const __m256i upper_mask =
          _mm256_load_si256(reinterpret_cast<const __m256i*>(
              kSmallStageMasks.upper[level].data()));
      for (size_t base = 0; base < values.size(); base += 32) {
        const __m256i input = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(values.data() + base));
        __m256i x = _mm256_shuffle_epi8(input, low_mask);
        __m256i y = _mm256_shuffle_epi8(input, high_mask);
        const __m256i coefficient =
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                skew.factors[level].data() + evaluation_offset + base));
        Butterfly<Inverse>(x, y, coefficient);
        const __m256i output = _mm256_blendv_epi8(x, y, upper_mask);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(values.data() + base),
                            output);
      }
      continue;
    }

    if (half == 16) {
      for (size_t base = 0; base < values.size(); base += 32) {
        const __m256i input = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(values.data() + base));
        __m256i x = _mm256_permute2x128_si256(input, input, 0x00);
        __m256i y = _mm256_permute2x128_si256(input, input, 0x11);
        const __m256i coefficient =
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                skew.factors[level].data() + evaluation_offset + base));
        Butterfly<Inverse>(x, y, coefficient);
        const __m256i output = _mm256_permute2x128_si256(x, y, 0x20);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(values.data() + base),
                            output);
      }
      continue;
    }

    for (size_t group = 0; group < values.size(); group += group_size) {
      const __m256i coefficient = _mm256_set1_epi8(
          static_cast<char>(skew.factors[level][evaluation_offset + group]));
      for (size_t i = 0; i < half; i += 32) {
        __m256i x = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(values.data() + group + i));
        __m256i y = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(values.data() + group + half + i));
        Butterfly<Inverse>(x, y, coefficient);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(values.data() + group + i), x);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(values.data() + group + half + i), y);
      }
    }
  }

  ConvertBasis(values, kIsomorphism.aes_to_cantor_matrix,
               kIsomorphism.aes_to_cantor);
}

void TransformGFNI(std::span<Element> values,
                   size_t block_size,
                   size_t evaluation_offset,
                   bool inverse) {
  if (inverse) {
    TransformGFNIImpl<true>(values, block_size, evaluation_offset);
  } else {
    TransformGFNIImpl<false>(values, block_size, evaluation_offset);
  }
}

#if defined(GF256_ENABLE_CODEWORD_CANTOR_AFFINE_EXPERIMENT)
template <size_t Half>
__m256i ShiftRightBytes(__m256i value) {
  return _mm256_srli_si256(value, Half);
}

template <size_t Half>
__m256i ShiftLeftBytes(__m256i value) {
  return _mm256_slli_si256(value, Half);
}

template <size_t Half>
__m256i CantorAffineFactors(const CantorAffineSkewTables& skew,
                            size_t level,
                            size_t absolute_base,
                            size_t pass) {
  if ((absolute_base & 31U) == 0) {
    return _mm256_load_si256(reinterpret_cast<const __m256i*>(
        skew.small_aligned[level][pass][absolute_base / 32].data()));
  }

  std::array<uint64_t, 4> matrix{};
  if constexpr (Half == 8) {
    matrix[0] = skew.factors[level][absolute_base];
    matrix[2] = skew.factors[level][absolute_base + 16];
  } else {
    constexpr size_t kGroupSize = 2 * Half;
    for (size_t qword = 0; qword < matrix.size(); ++qword) {
      matrix[qword] =
          skew.factors[level][absolute_base + 8 * qword + pass * kGroupSize];
    }
  }
  return _mm256_set_epi64x(
      static_cast<long long>(matrix[3]), static_cast<long long>(matrix[2]),
      static_cast<long long>(matrix[1]), static_cast<long long>(matrix[0]));
}

template <bool Inverse, size_t Half>
__m256i CantorAffineSmallStage(__m256i values,
                               const CantorAffineSkewTables& skew,
                               size_t level,
                               size_t absolute_base) {
  const __m256i upper_mask = _mm256_load_si256(
      reinterpret_cast<const __m256i*>(kSmallStageMasks.upper[level].data()));

  if constexpr (Inverse) {
    values = _mm256_xor_si256(
        values, _mm256_and_si256(ShiftLeftBytes<Half>(values), upper_mask));
  }

  const __m256i y_to_x = ShiftRightBytes<Half>(values);
  constexpr size_t kPassCount = Half < 4 ? 4 / Half : 1;
  for (size_t pass = 0; pass < kPassCount; ++pass) {
    const __m256i pass_mask =
        _mm256_load_si256(reinterpret_cast<const __m256i*>(
            kSmallStageMasks.affine_pass[level][pass].data()));
    const __m256i matrix =
        CantorAffineFactors<Half>(skew, level, absolute_base, pass);
    const __m256i product = _mm256_gf2p8affine_epi64_epi8(y_to_x, matrix, 0);
    values = _mm256_xor_si256(values, _mm256_and_si256(product, pass_mask));
  }

  if constexpr (!Inverse) {
    values = _mm256_xor_si256(
        values, _mm256_and_si256(ShiftLeftBytes<Half>(values), upper_mask));
  }
  return values;
}

template <bool Inverse>
void CantorAffineButterfly(__m256i& x, __m256i& y, __m256i matrix) {
  if constexpr (Inverse) {
    y = _mm256_xor_si256(y, x);
    x = _mm256_xor_si256(x, _mm256_gf2p8affine_epi64_epi8(y, matrix, 0));
  } else {
    x = _mm256_xor_si256(x, _mm256_gf2p8affine_epi64_epi8(y, matrix, 0));
    y = _mm256_xor_si256(y, x);
  }
}

template <bool Inverse>
void TransformCantorAffineImpl(std::span<Element> values,
                               size_t block_size,
                               size_t evaluation_offset) {
  if (block_size == 1 || values.size() < 32) {
    if constexpr (Inverse) {
      IFFTScalar(Context::Shared(), values, block_size, evaluation_offset);
    } else {
      FFTScalar(Context::Shared(), values, block_size, evaluation_offset);
    }
    return;
  }

  const CantorAffineSkewTables& skew = GetCantorAffineSkewTables();
  for (size_t half = Inverse ? 1 : block_size / 2;
       Inverse ? half < block_size : half != 0;
       half = Inverse ? half * 2 : half / 2) {
    const size_t group_size = 2 * half;
    const size_t level = Log2(half);
    if (half <= 8) {
      for (size_t base = 0; base < values.size(); base += 32) {
        __m256i packed = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(values.data() + base));
        switch (half) {
          case 1:
            packed = CantorAffineSmallStage<Inverse, 1>(
                packed, skew, level, evaluation_offset + base);
            break;
          case 2:
            packed = CantorAffineSmallStage<Inverse, 2>(
                packed, skew, level, evaluation_offset + base);
            break;
          case 4:
            packed = CantorAffineSmallStage<Inverse, 4>(
                packed, skew, level, evaluation_offset + base);
            break;
          case 8:
            packed = CantorAffineSmallStage<Inverse, 8>(
                packed, skew, level, evaluation_offset + base);
            break;
        }
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(values.data() + base),
                            packed);
      }
      continue;
    }

    if (half == 16) {
      for (size_t base = 0; base < values.size(); base += 32) {
        const __m256i input = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(values.data() + base));
        __m256i x = _mm256_permute2x128_si256(input, input, 0x00);
        __m256i y = _mm256_permute2x128_si256(input, input, 0x11);
        const __m256i matrix = _mm256_set1_epi64x(static_cast<long long>(
            skew.factors[level][evaluation_offset + base]));
        CantorAffineButterfly<Inverse>(x, y, matrix);
        const __m256i output = _mm256_permute2x128_si256(x, y, 0x20);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(values.data() + base),
                            output);
      }
      continue;
    }

    for (size_t group = 0; group < values.size(); group += group_size) {
      const __m256i matrix = _mm256_set1_epi64x(static_cast<long long>(
          skew.factors[level][evaluation_offset + group]));
      for (size_t i = 0; i < half; i += 32) {
        __m256i x = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(values.data() + group + i));
        __m256i y = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(values.data() + group + half + i));
        CantorAffineButterfly<Inverse>(x, y, matrix);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(values.data() + group + i), x);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(values.data() + group + half + i), y);
      }
    }
  }
}

void TransformCantorAffine(std::span<Element> values,
                           size_t block_size,
                           size_t evaluation_offset,
                           bool inverse) {
  if (inverse) {
    TransformCantorAffineImpl<true>(values, block_size, evaluation_offset);
  } else {
    TransformCantorAffineImpl<false>(values, block_size, evaluation_offset);
  }
}
#endif

#endif

Backend ResolveBackend(Backend backend) {
  if (backend != Backend::tuned) {
    return backend;
  }
  if (BackendAvailable(Backend::gfni256_affine)) {
    return Backend::gfni256_affine;
  }
  return BackendAvailable(Backend::avx2) ? Backend::avx2 : Backend::scalar;
}

Status Run(const Context& context,
           std::span<Element> values,
           size_t block_size,
           size_t evaluation_offset,
           Backend backend,
           bool inverse) {
  backend = ResolveBackend(backend);
  if (backend == Backend::scalar) {
    if (inverse) {
      IFFTScalar(context, values, block_size, evaluation_offset);
    } else {
      FFTScalar(context, values, block_size, evaluation_offset);
    }
    return Status::ok;
  }
#if defined(__AVX2__)
  if (backend == Backend::avx2) {
    if (inverse) {
      TransformAVX2<true>(context, values, block_size, evaluation_offset);
    } else {
      TransformAVX2<false>(context, values, block_size, evaluation_offset);
    }
    return Status::ok;
  }
#endif
#if defined(__GFNI__) && defined(__AVX2__)
  if (backend == Backend::gfni256_affine) {
    TransformGFNI(values, block_size, evaluation_offset, inverse);
    return Status::ok;
  }
#endif
  return Status::unsupported_backend;
}

}  // namespace

#if defined(__GFNI__) && defined(__AVX2__)
const std::array<Element, Context::kFieldSize>& CantorToAESMap() {
  return kIsomorphism.cantor_to_aes;
}

const std::array<Element, Context::kFieldSize>& AESToCantorMap() {
  return kIsomorphism.aes_to_cantor;
}

void ConvertCantorToAES(std::span<Element> values) {
  ConvertBasis(values, kIsomorphism.cantor_to_aes_matrix,
               kIsomorphism.cantor_to_aes);
}

void ConvertAESToCantor(std::span<Element> values) {
  ConvertBasis(values, kIsomorphism.aes_to_cantor_matrix,
               kIsomorphism.aes_to_cantor);
}
#endif

Status FFTCodewordBlocks(const Context& context,
                         std::span<Element> values,
                         size_t block_size,
                         size_t evaluation_offset,
                         Backend backend) {
  const Status status =
      Validate(values, block_size, evaluation_offset, backend);
  if (status != Status::ok) {
    return status;
  }
  return Run(context, values, block_size, evaluation_offset, backend, false);
}

Status IFFTCodewordBlocks(const Context& context,
                          std::span<Element> values,
                          size_t block_size,
                          size_t evaluation_offset,
                          Backend backend) {
  const Status status =
      Validate(values, block_size, evaluation_offset, backend);
  if (status != Status::ok) {
    return status;
  }
  return Run(context, values, block_size, evaluation_offset, backend, true);
}

#if defined(GF256_ENABLE_CODEWORD_CANTOR_AFFINE_EXPERIMENT)
Status FFTCodewordBlocksCantorAffine(const Context& context,
                                     std::span<Element> values,
                                     size_t block_size,
                                     size_t evaluation_offset) {
  (void)context;
  const Status status =
      Validate(values, block_size, evaluation_offset, Backend::gfni256_affine);
  if (status != Status::ok) {
    return status;
  }
#if defined(__GFNI__) && defined(__AVX2__)
  TransformCantorAffine(values, block_size, evaluation_offset, false);
  return Status::ok;
#else
  return Status::unsupported_backend;
#endif
}

Status IFFTCodewordBlocksCantorAffine(const Context& context,
                                      std::span<Element> values,
                                      size_t block_size,
                                      size_t evaluation_offset) {
  (void)context;
  const Status status =
      Validate(values, block_size, evaluation_offset, Backend::gfni256_affine);
  if (status != Status::ok) {
    return status;
  }
#if defined(__GFNI__) && defined(__AVX2__)
  TransformCantorAffine(values, block_size, evaluation_offset, true);
  return Status::ok;
#else
  return Status::unsupported_backend;
#endif
}
#endif

}  // namespace gf2p8::lch::detail
