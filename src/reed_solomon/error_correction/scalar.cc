#if defined(__i386__) || defined(__x86_64__) || defined(_M_IX86) || \
    defined(_M_X64)
#include <immintrin.h>
#endif

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>

#include "field.h"
#include "lin_chung_han/codeword_transform_internal.h"
#include "lin_chung_han/transform.h"
#include "reed_solomon/code_parameters.h"
#include "reed_solomon/error_correction/avx2_internal.h"
#include "reed_solomon/error_correction/internal.h"

namespace gf2p8::rs::detail::error_correction {
namespace {

constexpr size_t kFieldSize = lch::Context::kFieldSize;

using Values = std::array<Element, kFieldSize>;

CorrectionResult Result(CorrectionStatus status, size_t error_count = 0) {
  return {.status = status, .error_count = error_count};
}

bool RangesOverlap(std::span<const Element> first,
                   std::span<const Element> second) {
  if (first.empty() || second.empty()) {
    return false;
  }
  const uintptr_t first_begin = reinterpret_cast<uintptr_t>(first.data());
  const uintptr_t first_end = first_begin + first.size_bytes();
  const uintptr_t second_begin = reinterpret_cast<uintptr_t>(second.data());
  const uintptr_t second_end = second_begin + second.size_bytes();
  return first_begin < second_end && second_begin < first_end;
}

bool Transform(Values& values,
               size_t count,
               size_t evaluation_offset,
               bool inverse) {
  const std::span<Element> block(values.data(), count);
  const lch::Status status =
      inverse ? lch::detail::IFFTCodewordBlocks(lch::Context::Shared(), block,
                                                count, evaluation_offset,
                                                lch::Backend::tuned)
              : lch::detail::FFTCodewordBlocks(lch::Context::Shared(), block,
                                               count, evaluation_offset,
                                               lch::Backend::tuned);
  return status == lch::Status::ok;
}

Element Product(Element value,
                Element coefficient,
                const MultiplicationTables& tables) {
  const auto& row = tables.shuffle[coefficient];
  return row[value & 0x0f] ^ row[32 + (value >> 4)];
}

bool PolynomialDegree(const Values& coefficients,
                      size_t maximum_degree,
                      size_t& degree) {
  bool is_nonzero = false;
  degree = 0;
  for (size_t i = 0; i <= maximum_degree; ++i) {
    if (coefficients[i] != 0) {
      is_nonzero = true;
      degree = i;
    }
  }
  return is_nonzero;
}

void DifferentiateNovel(const Values& coefficients,
                        size_t degree,
                        Values& derivative) {
  std::fill_n(derivative.begin(), degree, Element{0});
  for (size_t source = 1; source <= degree; ++source) {
    for (size_t bits = source; bits != 0; bits &= bits - 1) {
      const size_t basis_term = std::countr_zero(bits);
      derivative[source ^ (size_t{1} << basis_term)] ^= coefficients[source];
    }
  }
}

#if defined(__GFNI__) && defined(__AVX2__)
template <bool UseFirst>
void UpdateSamplesGFNI(Element* first,
                       Element* second,
                       size_t begin,
                       size_t end,
                       size_t constraint,
                       Element first_discrepancy,
                       Element second_discrepancy,
                       const MultiplicationTables& tables) {
  const auto& cantor_to_aes = lch::detail::CantorToAESMap();
  const auto& aes_to_cantor = lch::detail::AESToCantorMap();
  const __m256i first_factor256 =
      _mm256_set1_epi8(static_cast<char>(first_discrepancy));
  const __m256i second_factor256 =
      _mm256_set1_epi8(static_cast<char>(second_discrepancy));
  const __m256i constraint256 =
      _mm256_set1_epi8(static_cast<char>(cantor_to_aes[constraint]));

  size_t i = begin;
  for (; i + 32 <= end; i += 32) {
    const __m256i old_first =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(first + i));
    const __m256i old_second =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(second + i));
    const __m256i points = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(cantor_to_aes.data() + i));
    const __m256i differences = _mm256_xor_si256(points, constraint256);
    const __m256i new_first =
        _mm256_xor_si256(_mm256_gf2p8mul_epi8(old_first, second_factor256),
                         _mm256_gf2p8mul_epi8(old_second, first_factor256));
    const __m256i source = UseFirst ? old_first : old_second;
    const __m256i new_second = _mm256_gf2p8mul_epi8(source, differences);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(first + i), new_first);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(second + i), new_second);
  }

  const __m128i first_factor128 =
      _mm_set1_epi8(static_cast<char>(first_discrepancy));
  const __m128i second_factor128 =
      _mm_set1_epi8(static_cast<char>(second_discrepancy));
  const __m128i constraint128 =
      _mm_set1_epi8(static_cast<char>(cantor_to_aes[constraint]));
  for (; i + 16 <= end; i += 16) {
    const __m128i old_first =
        _mm_loadu_si128(reinterpret_cast<const __m128i*>(first + i));
    const __m128i old_second =
        _mm_loadu_si128(reinterpret_cast<const __m128i*>(second + i));
    const __m128i points = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(cantor_to_aes.data() + i));
    const __m128i differences = _mm_xor_si128(points, constraint128);
    const __m128i new_first =
        _mm_xor_si128(_mm_gf2p8mul_epi8(old_first, second_factor128),
                      _mm_gf2p8mul_epi8(old_second, first_factor128));
    const __m128i source = UseFirst ? old_first : old_second;
    const __m128i new_second = _mm_gf2p8mul_epi8(source, differences);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(first + i), new_first);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(second + i), new_second);
  }

  const Element first_discrepancy_cantor = aes_to_cantor[first_discrepancy];
  const Element second_discrepancy_cantor = aes_to_cantor[second_discrepancy];
  for (; i < end; ++i) {
    const Element old_first_cantor = aes_to_cantor[first[i]];
    const Element old_second_cantor = aes_to_cantor[second[i]];
    const Element new_first_cantor =
        Product(old_first_cantor, second_discrepancy_cantor, tables) ^
        Product(old_second_cantor, first_discrepancy_cantor, tables);
    const Element source_cantor =
        UseFirst ? old_first_cantor : old_second_cantor;
    const Element new_second_cantor =
        Product(source_cantor, static_cast<Element>(i ^ constraint), tables);
    first[i] = cantor_to_aes[new_first_cantor];
    second[i] = cantor_to_aes[new_second_cantor];
  }
}
#endif

#if defined(__AVX2__)
template <bool UseFirst>
void UpdateSamplesAVX2(Element* first,
                       Element* second,
                       size_t begin,
                       size_t end,
                       size_t constraint,
                       Element first_discrepancy,
                       Element second_discrepancy,
                       const MultiplicationTables& tables) {
  static constexpr std::array<Element, kFieldSize> kEvaluationPoints = [] {
    std::array<Element, kFieldSize> points{};
    for (size_t i = 0; i < points.size(); ++i) {
      points[i] = static_cast<Element>(i);
    }
    return points;
  }();
  const __m256i constraint_vector =
      _mm256_set1_epi8(static_cast<char>(constraint));

  size_t i = begin;
  for (; i + 32 <= end; i += 32) {
    const __m256i old_first =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(first + i));
    const __m256i old_second =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(second + i));
    const __m256i new_first = _mm256_xor_si256(
        avx2::MultiplyFixed(old_first, second_discrepancy, tables),
        avx2::MultiplyFixed(old_second, first_discrepancy, tables));
    const __m256i source = UseFirst ? old_first : old_second;
    const __m256i point_differences = _mm256_xor_si256(
        _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(kEvaluationPoints.data() + i)),
        constraint_vector);
    const __m256i new_second =
        avx2::MultiplyVariable(source, point_differences, tables);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(first + i), new_first);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(second + i), new_second);
  }

  for (; i < end; ++i) {
    const Element old_first = first[i];
    const Element old_second = second[i];
    first[i] = Product(old_first, second_discrepancy, tables) ^
               Product(old_second, first_discrepancy, tables);
    second[i] = Product(UseFirst ? old_first : old_second,
                        static_cast<Element>(i ^ constraint), tables);
  }
}
#endif

size_t PublicPosition(CodeFamily family,
                      size_t data_count,
                      size_t recovery_count,
                      size_t native_position) {
  if (family == CodeFamily::low_rate) {
    return native_position;
  }
  return native_position < recovery_count ? data_count + native_position
                                          : native_position - recovery_count;
}

CorrectionStatus RecoverWithEvaluator(
    CodeFamily family,
    std::span<Element> data,
    std::span<const Element> recovery,
    std::span<Element> mutable_recovery,
    size_t recovery_count,
    std::span<const uint8_t> root_positions,
    const Values& locator_samples,
    const Values& locator_coefficients,
    size_t locator_degree,
    const Values& syndrome_samples,
    const MultiplicationTables& tables) {
  const size_t data_count = data.size();
  const size_t codeword_size = data_count + recovery_count;
  const size_t correction_radius = recovery_count / 2;

  // Tang-Han Algorithm 4 reconstructs the evaluator from z(a)=u(a)lambda(a)
  // at the syndrome points.  The common scale of z and lambda cancels in the
  // Forney ratios below.
  Values evaluator_samples;
  for (size_t i = 0; i <= correction_radius; ++i) {
    evaluator_samples[i] =
        Product(locator_samples[i], syndrome_samples[i], tables);
  }
  Values transform_values;
  std::copy_n(evaluator_samples.begin(), correction_radius,
              transform_values.begin());
  if (!Transform(transform_values, correction_radius, 0, true)) {
    return CorrectionStatus::reconstruction_failed;
  }
  Values evaluator_coefficients;
  std::copy_n(transform_values.begin(), correction_radius,
              evaluator_coefficients.begin());
  if (!Transform(transform_values, correction_radius, correction_radius,
                 false)) {
    return CorrectionStatus::reconstruction_failed;
  }
  evaluator_coefficients[correction_radius] =
      evaluator_samples[correction_radius] ^ transform_values[0];

  size_t evaluator_degree = 0;
  const bool evaluator_is_nonzero = PolynomialDegree(
      evaluator_coefficients, correction_radius, evaluator_degree);
  if (evaluator_is_nonzero && evaluator_degree >= locator_degree) {
    return CorrectionStatus::uncorrectable;
  }

  Values evaluator_values;
  std::fill_n(evaluator_values.begin(), codeword_size, Element{0});
  if (evaluator_is_nonzero) {
    for (size_t block = 0; block < codeword_size; block += recovery_count) {
      std::copy_n(evaluator_coefficients.begin(), evaluator_degree + 1,
                  evaluator_values.begin() + block);
    }
  }
  if (lch::detail::FFTCodewordBlocks(
          lch::Context::Shared(),
          std::span(evaluator_values).first(codeword_size), recovery_count, 0,
          lch::Backend::tuned) != lch::Status::ok) {
    return CorrectionStatus::reconstruction_failed;
  }

  std::fill_n(transform_values.begin(), recovery_count, Element{0});
  DifferentiateNovel(evaluator_coefficients, evaluator_degree,
                     transform_values);
  if (!Transform(transform_values, recovery_count, 0, false)) {
    return CorrectionStatus::reconstruction_failed;
  }

  Values corrections;
  const Element locator_leading = locator_coefficients[locator_degree];
  const size_t syndrome_level = std::countr_zero(recovery_count);
  const lch::Context& context = lch::Context::Shared();
  for (size_t root_index = 0; root_index < root_positions.size();
       ++root_index) {
    const size_t position = root_positions[root_index];
    Element locator_derivative = locator_leading;
    for (size_t other_index = 0; other_index < root_positions.size();
         ++other_index) {
      if (other_index != root_index) {
        locator_derivative = Product(
            locator_derivative,
            static_cast<Element>(position ^ root_positions[other_index]),
            tables);
      }
    }
    if (locator_derivative == 0) {
      return CorrectionStatus::uncorrectable;
    }

    Element correction = 0;
    if (position < recovery_count) {
      correction = syndrome_samples[position] ^
                   DivCantor(transform_values[position], locator_derivative);
    } else {
      const Element denominator = Product(
          context.Skew(syndrome_level, position), locator_derivative, tables);
      if (denominator == 0) {
        return CorrectionStatus::uncorrectable;
      }
      correction = DivCantor(evaluator_values[position], denominator);
    }
    if (correction == 0) {
      return CorrectionStatus::uncorrectable;
    }
    corrections[position] = correction;
  }

  // Verify the candidate correction transactionally before mutating caller
  // data.  This also rejects algebraically valid-looking over-radius results.
  if (family == CodeFamily::low_rate) {
    std::copy(data.begin(), data.end(), transform_values.begin());
    std::copy(recovery.begin(), recovery.end(),
              transform_values.begin() + data_count);
  } else {
    std::copy(recovery.begin(), recovery.end(), transform_values.begin());
    std::copy(data.begin(), data.end(),
              transform_values.begin() + recovery_count);
  }
  for (const uint8_t position : root_positions) {
    transform_values[position] ^= corrections[position];
  }
  if (lch::detail::IFFTCodewordBlocks(
          context, std::span(transform_values).first(codeword_size),
          recovery_count, 0, lch::Backend::tuned) != lch::Status::ok) {
    return CorrectionStatus::reconstruction_failed;
  }
  for (size_t i = 0; i < recovery_count; ++i) {
    Element syndrome = 0;
    for (size_t block = 0; block < codeword_size; block += recovery_count) {
      syndrome ^= transform_values[block + i];
    }
    if (syndrome != 0) {
      return CorrectionStatus::uncorrectable;
    }
  }

  for (const uint8_t position : root_positions) {
    const size_t data_position =
        PublicPosition(family, data_count, recovery_count, position);
    if (data_position < data_count) {
      data[data_position] ^= corrections[position];
    } else if (!mutable_recovery.empty()) {
      mutable_recovery[data_position - data_count] ^= corrections[position];
    }
  }
  return CorrectionStatus::ok;
}

}  // namespace

static CorrectionResult CorrectOneImpl(const LCHDecoder& decoder,
                                       std::span<Element> data,
                                       std::span<const Element> recovery,
                                       std::span<uint8_t> error_mask,
                                       std::span<Element> mutable_recovery) {
  if (RangesOverlap(error_mask, data) || RangesOverlap(error_mask, recovery)) {
    return Result(CorrectionStatus::invalid_argument);
  }
  std::fill(error_mask.begin(), error_mask.end(), uint8_t{0});

  if (RangesOverlap(data, recovery)) {
    return Result(CorrectionStatus::invalid_argument);
  }

  if (!decoder.Valid()) {
    return Result(CorrectionStatus::invalid_argument);
  }

  const size_t data_count = decoder.DataCount();
  const size_t recovery_count = decoder.RecoveryCount();
  const size_t codeword_size = data_count + recovery_count;
  if (data.size() != data_count || recovery.size() != recovery_count ||
      error_mask.size() != codeword_size) {
    return Result(CorrectionStatus::invalid_argument);
  }

  const CodeParameters parameters =
      MakeCodeParameters(data_count, recovery_count);
  if (!parameters.valid || !std::has_single_bit(codeword_size) ||
      !std::has_single_bit(recovery_count) ||
      parameters.transform_size != recovery_count ||
      parameters.mother_size != codeword_size) {
    return Result(CorrectionStatus::unsupported_dimensions);
  }

  Values syndrome;
  std::fill_n(syndrome.begin(), recovery_count, Element{0});
  Values transform_values;
  if (parameters.family == CodeFamily::low_rate) {
    std::copy(data.begin(), data.end(), transform_values.begin());
    std::copy(recovery.begin(), recovery.end(),
              transform_values.begin() + data_count);
  } else {
    std::copy(recovery.begin(), recovery.end(), transform_values.begin());
    std::copy(data.begin(), data.end(),
              transform_values.begin() + recovery_count);
  }
  if (lch::detail::IFFTCodewordBlocks(
          lch::Context::Shared(),
          std::span(transform_values).first(codeword_size), recovery_count, 0,
          lch::Backend::tuned) != lch::Status::ok) {
    return Result(CorrectionStatus::reconstruction_failed);
  }
  for (size_t block = 0; block < codeword_size; block += recovery_count) {
    for (size_t i = 0; i < recovery_count; ++i) {
      syndrome[i] ^= transform_values[block + i];
    }
  }

  const bool syndrome_is_zero =
      std::all_of(syndrome.begin(), syndrome.begin() + recovery_count,
                  [](Element value) { return value == 0; });
  if (syndrome_is_zero) {
    return Result(CorrectionStatus::ok);
  }
  if (recovery_count == 1) {
    return Result(CorrectionStatus::uncorrectable);
  }

  std::copy_n(syndrome.begin(), recovery_count, transform_values.begin());
  if (!Transform(transform_values, recovery_count, 0, false)) {
    return Result(CorrectionStatus::reconstruction_failed);
  }
  Values syndrome_samples;
  std::copy_n(transform_values.begin(), recovery_count,
              syndrome_samples.begin());

  Values d;
  Values g;
  std::copy_n(transform_values.begin(), recovery_count, d.begin());
  std::fill_n(g.begin(), recovery_count, Element{1});

  const size_t correction_radius = recovery_count / 2;
  Values locator_first;
  Values locator_second;
  std::fill_n(locator_first.begin(), correction_radius + 1, Element{1});
  std::fill_n(locator_second.begin(), correction_radius + 1, Element{0});

  size_t first_rank = 0;
  size_t second_rank = 1;
  const MultiplicationTables& tables = lch::Context::Shared().Tables();
#if defined(__GFNI__) && defined(__AVX2__)
  const bool use_gfni_fdma =
      recovery_count >= 64 &&
      lch::BackendAvailable(lch::Backend::gfni256_affine);
  if (use_gfni_fdma) {
    lch::detail::ConvertCantorToAES(std::span(d).first(recovery_count));
    lch::detail::ConvertCantorToAES(std::span(g).first(recovery_count));
    lch::detail::ConvertCantorToAES(
        std::span(locator_first).first(correction_radius + 1));
    lch::detail::ConvertCantorToAES(
        std::span(locator_second).first(correction_radius + 1));
  }
#endif
#if defined(__AVX2__)
#if defined(__GFNI__)
  const bool use_avx2_fdma = !use_gfni_fdma && recovery_count >= 32;
#else
  const bool use_avx2_fdma = recovery_count >= 32;
#endif
#endif
  for (size_t constraint = 0; constraint < recovery_count; ++constraint) {
    const Element first_discrepancy = d[constraint];
    const Element second_discrepancy = g[constraint];
    const bool first_update =
        second_discrepancy == 0 ||
        (first_discrepancy != 0 && first_rank < second_rank);

#if defined(__GFNI__) && defined(__AVX2__)
    if (use_gfni_fdma) {
      if (first_update) {
        UpdateSamplesGFNI<true>(d.data(), g.data(), constraint + 1,
                                recovery_count, constraint, first_discrepancy,
                                second_discrepancy, tables);
        UpdateSamplesGFNI<true>(locator_first.data(), locator_second.data(), 0,
                                correction_radius + 1, constraint,
                                first_discrepancy, second_discrepancy, tables);
      } else {
        UpdateSamplesGFNI<false>(d.data(), g.data(), constraint + 1,
                                 recovery_count, constraint, first_discrepancy,
                                 second_discrepancy, tables);
        UpdateSamplesGFNI<false>(locator_first.data(), locator_second.data(), 0,
                                 correction_radius + 1, constraint,
                                 first_discrepancy, second_discrepancy, tables);
      }
    } else
#endif
#if defined(__AVX2__)
        if (use_avx2_fdma) {
      if (first_update) {
        UpdateSamplesAVX2<true>(d.data(), g.data(), constraint + 1,
                                recovery_count, constraint, first_discrepancy,
                                second_discrepancy, tables);
        UpdateSamplesAVX2<true>(locator_first.data(), locator_second.data(), 0,
                                correction_radius + 1, constraint,
                                first_discrepancy, second_discrepancy, tables);
      } else {
        UpdateSamplesAVX2<false>(d.data(), g.data(), constraint + 1,
                                 recovery_count, constraint, first_discrepancy,
                                 second_discrepancy, tables);
        UpdateSamplesAVX2<false>(locator_first.data(), locator_second.data(), 0,
                                 correction_radius + 1, constraint,
                                 first_discrepancy, second_discrepancy, tables);
      }
    } else
#endif
    {
      for (size_t i = constraint + 1; i < recovery_count; ++i) {
        const Element old_first = d[i];
        const Element old_second = g[i];
        const Element point_difference = static_cast<Element>(i ^ constraint);
        d[i] = Product(old_first, second_discrepancy, tables) ^
               Product(old_second, first_discrepancy, tables);
        g[i] = Product(first_update ? old_first : old_second, point_difference,
                       tables);
      }

      for (size_t i = 0; i <= correction_radius; ++i) {
        const Element old_first = locator_first[i];
        const Element old_second = locator_second[i];
        const Element point_difference = static_cast<Element>(i ^ constraint);
        locator_first[i] = Product(old_first, second_discrepancy, tables) ^
                           Product(old_second, first_discrepancy, tables);
        locator_second[i] = Product(first_update ? old_first : old_second,
                                    point_difference, tables);
      }
    }

    if (first_update) {
      const size_t old_first_rank = first_rank;
      first_rank = second_rank;
      second_rank = old_first_rank + 2;
    } else {
      second_rank += 2;
    }
  }

  const bool use_first = first_rank < second_rank;
  const size_t locator_rank = use_first ? first_rank : second_rank;
  if ((locator_rank & 1U) != 0) {
    return Result(CorrectionStatus::uncorrectable);
  }
  const size_t locator_degree = locator_rank / 2;
  if (locator_degree == 0 || locator_degree > correction_radius) {
    return Result(CorrectionStatus::uncorrectable);
  }
#if defined(__GFNI__) && defined(__AVX2__)
  if (use_gfni_fdma) {
    Values& selected_locator = use_first ? locator_first : locator_second;
    lch::detail::ConvertAESToCantor(
        std::span(selected_locator).first(correction_radius + 1));
  }
#endif
  const Values& locator_samples = use_first ? locator_first : locator_second;

  std::copy_n(locator_samples.begin(), correction_radius,
              transform_values.begin());
  if (!Transform(transform_values, correction_radius, 0, true)) {
    return Result(CorrectionStatus::reconstruction_failed);
  }
  Values locator_coefficients;
  std::copy_n(transform_values.begin(), correction_radius,
              locator_coefficients.begin());
  if (!Transform(transform_values, correction_radius, correction_radius,
                 false)) {
    return Result(CorrectionStatus::reconstruction_failed);
  }
  locator_coefficients[correction_radius] =
      locator_samples[correction_radius] ^ transform_values[0];

  size_t reconstructed_degree = 0;
  if (!PolynomialDegree(locator_coefficients, correction_radius,
                        reconstructed_degree) ||
      reconstructed_degree != locator_degree) {
    return Result(CorrectionStatus::uncorrectable);
  }

  std::array<uint8_t, kFieldSize> root_positions;
  size_t root_count = 0;
  for (size_t block = 0; block < codeword_size; block += recovery_count) {
    std::fill_n(transform_values.begin() + block, recovery_count, Element{0});
    std::copy_n(locator_coefficients.begin(), locator_degree + 1,
                transform_values.begin() + block);
  }
  if (lch::detail::FFTCodewordBlocks(
          lch::Context::Shared(),
          std::span(transform_values).first(codeword_size), recovery_count, 0,
          lch::Backend::tuned) != lch::Status::ok) {
    return Result(CorrectionStatus::reconstruction_failed);
  }
  for (size_t block = 0; block < codeword_size; block += recovery_count) {
    for (size_t i = 0; i < recovery_count; ++i) {
      if (transform_values[block + i] == 0) {
        const size_t position = block + i;
        root_positions[root_count] = static_cast<uint8_t>(position);
        ++root_count;
      }
    }
  }
  if (root_count != locator_degree) {
    return Result(CorrectionStatus::uncorrectable);
  }

  bool has_data_error = false;
  for (size_t i = 0; i < root_count; ++i) {
    const size_t position = PublicPosition(parameters.family, data_count,
                                           recovery_count, root_positions[i]);
    has_data_error |= position < data_count;
  }

  CorrectionStatus recovery_status = CorrectionStatus::ok;
  // Whole-codeword mode verifies every candidate, including parity-only roots.
  // Retain the existing data-only fast path for CorrectOne/CorrectBatch callers.
  if (has_data_error && root_count == 1 && mutable_recovery.empty()) {
    // Every aligned R-point native Cantor IFFT has unit leading Lagrange
    // coefficient. Therefore the highest syndrome coefficient is the error
    // magnitude when exactly one error is present.
    const Element magnitude = syndrome[recovery_count - 1];
    const size_t data_index = PublicPosition(parameters.family, data_count,
                                             recovery_count, root_positions[0]);
    data[data_index] ^= magnitude;
  } else if (has_data_error || !mutable_recovery.empty()) {
    recovery_status = RecoverWithEvaluator(
        parameters.family, data, recovery, mutable_recovery, recovery_count,
        std::span(root_positions).first(root_count), locator_samples,
        locator_coefficients, locator_degree, syndrome_samples, tables);
  }
  if (recovery_status != CorrectionStatus::ok) {
    return Result(recovery_status);
  }
  for (size_t i = 0; i < root_count; ++i) {
    const size_t position = PublicPosition(parameters.family, data_count,
                                           recovery_count, root_positions[i]);
    error_mask[position] = 1;
  }
  return Result(CorrectionStatus::ok, root_count);
}

CorrectionResult CorrectOne(const LCHDecoder& decoder,
                            std::span<Element> data,
                            std::span<const Element> recovery,
                            std::span<uint8_t> error_mask) {
  return CorrectOneImpl(decoder, data, recovery, error_mask, {});
}

}  // namespace gf2p8::rs::detail::error_correction

namespace gf2p8::rs {

CorrectionResult CorrectCodeword(const LCHDecoder& decoder,
                                std::span<Element> codeword) {
  if (!decoder.Valid() ||
      codeword.size() != decoder.DataCount() + decoder.RecoveryCount()) {
    return {.status = CorrectionStatus::invalid_argument};
  }
  std::array<uint8_t, 256> mask{};
  auto recovery = codeword.subspan(decoder.DataCount());
  return detail::error_correction::CorrectOneImpl(
      decoder, codeword.first(decoder.DataCount()), recovery,
      std::span(mask).first(codeword.size()), recovery);
}

}  // namespace gf2p8::rs
