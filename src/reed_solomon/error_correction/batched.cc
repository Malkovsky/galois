#if defined(__i386__) || defined(__x86_64__) || defined(_M_IX86) || \
    defined(_M_X64)
#include <immintrin.h>
#endif

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <span>

#include "lin_chung_han/codeword_transform_internal.h"
#include "lin_chung_han/kernels_internal.h"
#include "lin_chung_han/transform.h"
#include "reed_solomon/code_parameters.h"
#include "reed_solomon/error_correction/internal.h"

namespace gf2p8::rs::detail::error_correction {
namespace {

constexpr size_t kFieldSize = lch::Context::kFieldSize;

using Values = std::array<Element, kFieldSize>;

struct ByteRange {
  uintptr_t begin;
  uintptr_t end;
};

bool AddRange(const void* pointer,
              size_t byte_count,
              std::span<ByteRange> ranges,
              size_t& range_count) {
  if (byte_count == 0) {
    return true;
  }
  if (pointer == nullptr) {
    return false;
  }
  const uintptr_t begin = reinterpret_cast<uintptr_t>(pointer);
  if (begin > std::numeric_limits<uintptr_t>::max() - byte_count) {
    return false;
  }
  ranges[range_count++] = {.begin = begin, .end = begin + byte_count};
  return true;
}

bool RangesAreDisjoint(std::span<ByteRange> ranges) {
  std::sort(ranges.begin(), ranges.end(),
            [](const ByteRange& first, const ByteRange& second) {
              return first.begin < second.begin;
            });
  for (size_t i = 1; i < ranges.size(); ++i) {
    if (ranges[i - 1].end > ranges[i].begin) {
      return false;
    }
  }
  return true;
}

bool SupportedDimensions(const LCHDecoder& decoder,
                         CodeParameters& parameters) {
  if (!decoder.Valid()) {
    return false;
  }
  const size_t data_count = decoder.DataCount();
  const size_t recovery_count = decoder.RecoveryCount();
  const size_t codeword_size = data_count + recovery_count;
  parameters = MakeCodeParameters(data_count, recovery_count);
  return parameters.valid && std::has_single_bit(codeword_size) &&
         std::has_single_bit(recovery_count) &&
         parameters.transform_size == recovery_count &&
         parameters.mother_size == codeword_size;
}

#if defined(__GFNI__) && defined(__AVX2__)
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

const Element* NativeSource(std::span<Element* const> data,
                            std::span<const Element* const> recovery,
                            CodeFamily family,
                            size_t native_position) {
  if (family == CodeFamily::low_rate) {
    return native_position < data.size()
               ? data[native_position]
               : recovery[native_position - data.size()];
  }
  return native_position < recovery.size()
             ? recovery[native_position]
             : data[native_position - recovery.size()];
}
#endif

CorrectionStatus CorrectColumnsScalar(const LCHDecoder& decoder,
                                      std::span<Element* const> data,
                                      std::span<const Element* const> recovery,
                                      size_t byte_count,
                                      size_t first_column,
                                      std::span<CorrectionResult> results,
                                      std::span<uint8_t> error_masks) {
  const size_t data_count = data.size();
  const size_t recovery_count = recovery.size();
  const size_t codeword_size = data_count + recovery_count;
  Values data_values{};
  Values recovery_values{};
  std::array<uint8_t, kFieldSize> mask{};

  for (size_t column = first_column; column < byte_count; ++column) {
    for (size_t i = 0; i < data_count; ++i) {
      data_values[i] = data[i][column];
    }
    for (size_t i = 0; i < recovery_count; ++i) {
      recovery_values[i] = recovery[i][column];
    }
    const CorrectionResult result = CorrectOne(
        decoder, std::span(data_values).first(data_count),
        std::span<const Element>(recovery_values).first(recovery_count),
        std::span(mask).first(codeword_size));
    results[column] = result;
    if (result.status == CorrectionStatus::ok) {
      for (size_t i = 0; i < data_count; ++i) {
        data[i][column] = data_values[i];
      }
    }
    for (size_t position = 0; position < codeword_size; ++position) {
      error_masks[position * byte_count + column] = mask[position];
    }
  }
  return CorrectionStatus::ok;
}

#if defined(__GFNI__) && defined(__AVX2__)

constexpr size_t kBatchLanes = 32;
constexpr size_t kMaximumLocatorSamples = kFieldSize / 4 + 1;

using Rows = std::array<Element, kFieldSize * kBatchLanes>;
using SyndromeRows = std::array<Element, (kFieldSize / 2) * kBatchLanes>;
using SampleRows = std::array<Element, kMaximumLocatorSamples * kBatchLanes>;

template <typename Storage>
Element* Row(Storage& storage, size_t row) {
  return storage.data() + row * kBatchLanes;
}

template <typename Storage>
const Element* Row(const Storage& storage, size_t row) {
  return storage.data() + row * kBatchLanes;
}

void FFT32Rows(Element* rows,
               size_t row_count,
               size_t evaluation_offset,
               const lch::detail::ResolvedKernels& kernels) {
  const lch::Context& context = lch::Context::Shared();
  size_t group_size = row_count;
  for (; group_size >= 4; group_size /= 4) {
    const size_t distance = group_size / 4;
    const size_t top_level = std::countr_zero(group_size) - 1;
    const size_t low_level = top_level - 1;
    for (size_t block = 0; block < row_count; block += group_size) {
      const Element top = context.Skew(top_level, evaluation_offset ^ block);
      const Element low = context.Skew(low_level, evaluation_offset ^ block);
      const Element high =
          context.Skew(low_level, evaluation_offset ^ (block + 2 * distance));
      for (size_t i = 0; i < distance; ++i) {
        kernels.fft_radix4(rows + (block + i) * kBatchLanes,
                           rows + (block + distance + i) * kBatchLanes,
                           rows + (block + 2 * distance + i) * kBatchLanes,
                           rows + (block + 3 * distance + i) * kBatchLanes,
                           kBatchLanes, top, low, high, context.Tables());
      }
    }
  }
  if (group_size == 2) {
    for (size_t block = 0; block < row_count; block += 2) {
      kernels.fft_radix2(rows + block * kBatchLanes,
                         rows + (block + 1) * kBatchLanes, kBatchLanes,
                         context.Skew(0, evaluation_offset ^ block),
                         context.Tables());
    }
  }
}

void IFFT32Rows(Element* rows,
                size_t row_count,
                size_t evaluation_offset,
                const lch::detail::ResolvedKernels& kernels) {
  const lch::Context& context = lch::Context::Shared();
  size_t distance = 1;
  for (; 4 * distance <= row_count; distance *= 4) {
    const size_t group_size = 4 * distance;
    const size_t low_level = std::countr_zero(distance);
    const size_t top_level = low_level + 1;
    for (size_t block = 0; block < row_count; block += group_size) {
      const Element top = context.Skew(top_level, evaluation_offset ^ block);
      const Element low = context.Skew(low_level, evaluation_offset ^ block);
      const Element high =
          context.Skew(low_level, evaluation_offset ^ (block + 2 * distance));
      for (size_t i = 0; i < distance; ++i) {
        kernels.ifft_radix4(rows + (block + i) * kBatchLanes,
                            rows + (block + distance + i) * kBatchLanes,
                            rows + (block + 2 * distance + i) * kBatchLanes,
                            rows + (block + 3 * distance + i) * kBatchLanes,
                            kBatchLanes, top, low, high, context.Tables());
      }
    }
  }
  if (distance < row_count) {
    const size_t half = distance;
    const Element coefficient =
        context.Skew(std::countr_zero(half), evaluation_offset);
    for (size_t i = 0; i < half; ++i) {
      kernels.ifft_radix2(rows + i * kBatchLanes,
                          rows + (half + i) * kBatchLanes, kBatchLanes,
                          coefficient, context.Tables());
    }
  }
}

void FFT32Blocks(Rows& rows,
                 size_t codeword_size,
                 size_t block_size,
                 const lch::detail::ResolvedKernels& kernels) {
  for (size_t block = 0; block < codeword_size; block += block_size) {
    FFT32Rows(Row(rows, block), block_size, block, kernels);
  }
}

void IFFT32Blocks(Rows& rows,
                  size_t codeword_size,
                  size_t block_size,
                  const lch::detail::ResolvedKernels& kernels) {
  for (size_t block = 0; block < codeword_size; block += block_size) {
    IFFT32Rows(Row(rows, block), block_size, block, kernels);
  }
}

__m256i LaneMask(uint32_t bits) {
  alignas(32) std::array<Element, kBatchLanes> bytes;
  for (size_t lane = 0; lane < kBatchLanes; ++lane) {
    bytes[lane] = static_cast<Element>(0U - ((bits >> lane) & 1U));
  }
  return _mm256_load_si256(reinterpret_cast<const __m256i*>(bytes.data()));
}

uint32_t ZeroMask(__m256i values) {
  return static_cast<uint32_t>(
      _mm256_movemask_epi8(_mm256_cmpeq_epi8(values, _mm256_setzero_si256())));
}

void UpdateBatchRows(Element* first,
                     Element* second,
                     size_t begin,
                     size_t end,
                     size_t constraint,
                     __m256i first_discrepancy,
                     __m256i second_discrepancy,
                     __m256i first_update_mask) {
  const auto& cantor_to_aes = lch::detail::CantorToAESMap();
  for (size_t i = begin; i < end; ++i) {
    __m256i old_first = _mm256_load_si256(
        reinterpret_cast<const __m256i*>(first + i * kBatchLanes));
    __m256i old_second = _mm256_load_si256(
        reinterpret_cast<const __m256i*>(second + i * kBatchLanes));
    const __m256i new_first =
        _mm256_xor_si256(_mm256_gf2p8mul_epi8(old_first, second_discrepancy),
                         _mm256_gf2p8mul_epi8(old_second, first_discrepancy));
    const __m256i source =
        _mm256_blendv_epi8(old_second, old_first, first_update_mask);
    const __m256i point_difference = _mm256_set1_epi8(
        static_cast<char>(cantor_to_aes[static_cast<Element>(i ^ constraint)]));
    const __m256i new_second = _mm256_gf2p8mul_epi8(source, point_difference);
    _mm256_store_si256(reinterpret_cast<__m256i*>(first + i * kBatchLanes),
                       new_first);
    _mm256_store_si256(reinterpret_cast<__m256i*>(second + i * kBatchLanes),
                       new_second);
  }
}

void DifferentiateRows(const Element* coefficients,
                       size_t maximum_degree,
                       Element* derivative,
                       size_t output_count) {
  std::fill_n(derivative, output_count * kBatchLanes, Element{0});
  for (size_t source = 1; source <= maximum_degree; ++source) {
    const __m256i coefficient = _mm256_load_si256(
        reinterpret_cast<const __m256i*>(coefficients + source * kBatchLanes));
    for (size_t bits = source; bits != 0; bits &= bits - 1) {
      const size_t basis_term = std::countr_zero(bits);
      Element* destination =
          derivative + (source ^ (size_t{1} << basis_term)) * kBatchLanes;
      const __m256i old =
          _mm256_load_si256(reinterpret_cast<const __m256i*>(destination));
      _mm256_store_si256(reinterpret_cast<__m256i*>(destination),
                         _mm256_xor_si256(old, coefficient));
    }
  }
}

__m256i InvertAES(__m256i value) {
  const __m256i power2 = _mm256_gf2p8mul_epi8(value, value);
  const __m256i power3 = _mm256_gf2p8mul_epi8(power2, value);
  const __m256i power6 = _mm256_gf2p8mul_epi8(power3, power3);
  const __m256i power12 = _mm256_gf2p8mul_epi8(power6, power6);
  const __m256i power15 = _mm256_gf2p8mul_epi8(power12, power3);
  const __m256i power30 = _mm256_gf2p8mul_epi8(power15, power15);
  const __m256i power60 = _mm256_gf2p8mul_epi8(power30, power30);
  const __m256i power120 = _mm256_gf2p8mul_epi8(power60, power60);
  const __m256i power240 = _mm256_gf2p8mul_epi8(power120, power120);
  const __m256i power252 = _mm256_gf2p8mul_epi8(power240, power12);
  return _mm256_gf2p8mul_epi8(power252, power2);
}

uint32_t DivideRows(const Element* numerator,
                    const Element* denominator,
                    Element denominator_scale,
                    Element* quotient,
                    uint32_t active_lanes) {
  alignas(32) std::array<Element, kBatchLanes> numerator_aes;
  alignas(32) std::array<Element, kBatchLanes> denominator_aes;
  std::memcpy(numerator_aes.data(), numerator, kBatchLanes);
  std::memcpy(denominator_aes.data(), denominator, kBatchLanes);
  lch::detail::ConvertCantorToAES(numerator_aes);
  lch::detail::ConvertCantorToAES(denominator_aes);

  const __m256i numerator_vector =
      _mm256_load_si256(reinterpret_cast<const __m256i*>(numerator_aes.data()));
  __m256i denominator_vector = _mm256_load_si256(
      reinterpret_cast<const __m256i*>(denominator_aes.data()));
  if (denominator_scale != 1) {
    const Element scale_aes = lch::detail::CantorToAESMap()[denominator_scale];
    denominator_vector = _mm256_gf2p8mul_epi8(
        denominator_vector, _mm256_set1_epi8(static_cast<char>(scale_aes)));
  }
  const uint32_t zero_denominators = ZeroMask(denominator_vector);
  const __m256i quotient_aes =
      _mm256_gf2p8mul_epi8(numerator_vector, InvertAES(denominator_vector));
  _mm256_store_si256(reinterpret_cast<__m256i*>(numerator_aes.data()),
                     quotient_aes);
  lch::detail::ConvertAESToCantor(numerator_aes);
  std::memcpy(quotient, numerator_aes.data(), kBatchLanes);
  return active_lanes & zero_denominators;
}

void CopyNativeCodeword(std::span<Element* const> data,
                        std::span<const Element* const> recovery,
                        CodeFamily family,
                        size_t column,
                        Rows& destination) {
  const size_t codeword_size = data.size() + recovery.size();
  for (size_t position = 0; position < codeword_size; ++position) {
    std::memcpy(Row(destination, position),
                NativeSource(data, recovery, family, position) + column,
                kBatchLanes);
  }
}

uint32_t SyndromeZeroMask(const Rows& transformed,
                          size_t codeword_size,
                          size_t recovery_count) {
  __m256i nonzero = _mm256_setzero_si256();
  for (size_t i = 0; i < recovery_count; ++i) {
    __m256i value = _mm256_load_si256(
        reinterpret_cast<const __m256i*>(Row(transformed, i)));
    for (size_t block = recovery_count; block < codeword_size;
         block += recovery_count) {
      value = _mm256_xor_si256(
          value, _mm256_load_si256(reinterpret_cast<const __m256i*>(
                     Row(transformed, block + i))));
    }
    nonzero = _mm256_or_si256(nonzero, value);
  }
  return ZeroMask(nonzero);
}

void FoldSyndrome(Rows& transformed,
                  size_t codeword_size,
                  size_t recovery_count,
                  const lch::detail::ResolvedKernels& kernels) {
  for (size_t block = recovery_count; block < codeword_size;
       block += recovery_count) {
    for (size_t i = 0; i < recovery_count; ++i) {
      kernels.xor_one(Row(transformed, i), Row(transformed, block + i),
                      kBatchLanes);
    }
  }
}

void PublishChunkResults(std::span<CorrectionResult> results,
                         std::span<uint8_t> error_masks,
                         size_t byte_count,
                         size_t column,
                         CodeFamily family,
                         size_t data_count,
                         size_t recovery_count,
                         uint32_t clean_lanes,
                         uint32_t corrected_lanes,
                         const std::array<uint32_t, kFieldSize>& root_masks,
                         const std::array<uint16_t, kBatchLanes>& root_counts) {
  const size_t codeword_size = data_count + recovery_count;
  for (size_t lane = 0; lane < kBatchLanes; ++lane) {
    const uint32_t bit = uint32_t{1} << lane;
    if ((clean_lanes & bit) != 0) {
      results[column + lane] = {.status = CorrectionStatus::ok,
                                .error_count = 0};
    } else if ((corrected_lanes & bit) != 0) {
      results[column + lane] = {
          .status = CorrectionStatus::ok,
          .error_count = static_cast<size_t>(root_counts[lane])};
    } else {
      results[column + lane] = {.status = CorrectionStatus::uncorrectable,
                                .error_count = 0};
    }
  }
  for (size_t native_position = 0; native_position < codeword_size;
       ++native_position) {
    const size_t public_position =
        PublicPosition(family, data_count, recovery_count, native_position);
    const uint32_t roots = root_masks[native_position] & corrected_lanes;
    for (size_t lane = 0; lane < kBatchLanes; ++lane) {
      error_masks[public_position * byte_count + column + lane] =
          static_cast<uint8_t>((roots >> lane) & 1U);
    }
  }
}

void CorrectChunk32(std::span<Element* const> data,
                    std::span<const Element* const> recovery,
                    size_t byte_count,
                    size_t column,
                    const CodeParameters& parameters,
                    const lch::detail::ResolvedKernels& kernels,
                    std::span<CorrectionResult> results,
                    std::span<uint8_t> error_masks) {
  const size_t data_count = data.size();
  const size_t recovery_count = recovery.size();
  const size_t codeword_size = data_count + recovery_count;
  const size_t correction_radius = recovery_count / 2;
  const lch::Context& context = lch::Context::Shared();

  alignas(32) Rows work{};
  alignas(32) SyndromeRows syndrome_or_scratch{};
  alignas(32) Rows derivative{};
  alignas(32) SampleRows locator_first{};
  alignas(32) SampleRows locator_second{};
  std::array<uint32_t, kFieldSize> root_masks{};
  std::array<uint16_t, kBatchLanes> first_ranks{};
  std::array<uint16_t, kBatchLanes> second_ranks{};
  std::array<uint16_t, kBatchLanes> locator_degrees{};
  std::array<uint16_t, kBatchLanes> root_counts{};
  second_ranks.fill(1);

  CopyNativeCodeword(data, recovery, parameters.family, column, work);
  IFFT32Blocks(work, codeword_size, recovery_count, kernels);
  FoldSyndrome(work, codeword_size, recovery_count, kernels);

  __m256i syndrome_nonzero = _mm256_setzero_si256();
  for (size_t i = 0; i < recovery_count; ++i) {
    syndrome_nonzero = _mm256_or_si256(
        syndrome_nonzero,
        _mm256_load_si256(reinterpret_cast<const __m256i*>(Row(work, i))));
  }
  const uint32_t clean_lanes = ZeroMask(syndrome_nonzero);
  if (clean_lanes == std::numeric_limits<uint32_t>::max()) {
    PublishChunkResults(results, error_masks, byte_count, column,
                        parameters.family, data_count, recovery_count,
                        clean_lanes, 0, root_masks, root_counts);
    return;
  }
  if (recovery_count == 1) {
    PublishChunkResults(results, error_masks, byte_count, column,
                        parameters.family, data_count, recovery_count,
                        clean_lanes, 0, root_masks, root_counts);
    return;
  }

  FFT32Rows(work.data(), recovery_count, 0, kernels);
  std::memcpy(syndrome_or_scratch.data(), work.data(),
              recovery_count * kBatchLanes);
  std::fill_n(derivative.data(), recovery_count * kBatchLanes, Element{1});
  std::fill_n(locator_first.data(), (correction_radius + 1) * kBatchLanes,
              Element{1});
  std::fill_n(locator_second.data(), (correction_radius + 1) * kBatchLanes,
              Element{0});

  lch::detail::ConvertCantorToAES(
      std::span(syndrome_or_scratch).first(recovery_count * kBatchLanes));
  lch::detail::ConvertCantorToAES(
      std::span(derivative).first(recovery_count * kBatchLanes));
  lch::detail::ConvertCantorToAES(
      std::span(locator_first).first((correction_radius + 1) * kBatchLanes));
  lch::detail::ConvertCantorToAES(
      std::span(locator_second).first((correction_radius + 1) * kBatchLanes));

  for (size_t constraint = 0; constraint < recovery_count; ++constraint) {
    const __m256i first_discrepancy = _mm256_load_si256(
        reinterpret_cast<const __m256i*>(Row(syndrome_or_scratch, constraint)));
    const __m256i second_discrepancy = _mm256_load_si256(
        reinterpret_cast<const __m256i*>(Row(derivative, constraint)));
    alignas(32) std::array<Element, kBatchLanes> first_values;
    alignas(32) std::array<Element, kBatchLanes> second_values;
    alignas(32) std::array<Element, kBatchLanes> update_bytes;
    _mm256_store_si256(reinterpret_cast<__m256i*>(first_values.data()),
                       first_discrepancy);
    _mm256_store_si256(reinterpret_cast<__m256i*>(second_values.data()),
                       second_discrepancy);
    for (size_t lane = 0; lane < kBatchLanes; ++lane) {
      const bool first_update =
          second_values[lane] == 0 ||
          (first_values[lane] != 0 && first_ranks[lane] < second_ranks[lane]);
      update_bytes[lane] = first_update ? 0xff : 0;
      if (first_update) {
        const uint16_t old_first_rank = first_ranks[lane];
        first_ranks[lane] = second_ranks[lane];
        second_ranks[lane] = static_cast<uint16_t>(old_first_rank + 2);
      } else {
        second_ranks[lane] = static_cast<uint16_t>(second_ranks[lane] + 2);
      }
    }
    const __m256i update_mask = _mm256_load_si256(
        reinterpret_cast<const __m256i*>(update_bytes.data()));
    UpdateBatchRows(syndrome_or_scratch.data(), derivative.data(),
                    constraint + 1, recovery_count, constraint,
                    first_discrepancy, second_discrepancy, update_mask);
    UpdateBatchRows(locator_first.data(), locator_second.data(), 0,
                    correction_radius + 1, constraint, first_discrepancy,
                    second_discrepancy, update_mask);
  }

  // FDMA mutated its discrepancy rows, but work still holds the original
  // syndrome evaluations needed to construct z and inside-subspace magnitudes.
  std::memcpy(syndrome_or_scratch.data(), work.data(),
              recovery_count * kBatchLanes);

  alignas(32) std::array<Element, kBatchLanes> selection_bytes;
  uint32_t candidate_lanes = ~clean_lanes;
  for (size_t lane = 0; lane < kBatchLanes; ++lane) {
    const bool use_first = first_ranks[lane] < second_ranks[lane];
    selection_bytes[lane] = use_first ? 0xff : 0;
    const size_t locator_rank =
        use_first ? first_ranks[lane] : second_ranks[lane];
    const size_t locator_degree = locator_rank / 2;
    locator_degrees[lane] = static_cast<uint16_t>(locator_degree);
    if ((locator_rank & 1U) != 0 || locator_degree == 0 ||
        locator_degree > correction_radius) {
      candidate_lanes &= ~(uint32_t{1} << lane);
    }
  }
  const __m256i selection_mask = _mm256_load_si256(
      reinterpret_cast<const __m256i*>(selection_bytes.data()));
  for (size_t i = 0; i <= correction_radius; ++i) {
    const __m256i first = _mm256_load_si256(
        reinterpret_cast<const __m256i*>(Row(locator_first, i)));
    const __m256i second = _mm256_load_si256(
        reinterpret_cast<const __m256i*>(Row(locator_second, i)));
    _mm256_store_si256(reinterpret_cast<__m256i*>(Row(locator_first, i)),
                       _mm256_blendv_epi8(second, first, selection_mask));
  }
  lch::detail::ConvertAESToCantor(
      std::span(locator_first).first((correction_radius + 1) * kBatchLanes));
  std::memcpy(locator_second.data(), locator_first.data(),
              (correction_radius + 1) * kBatchLanes);

  IFFT32Rows(locator_first.data(), correction_radius, 0, kernels);
  std::memcpy(derivative.data(), locator_first.data(),
              correction_radius * kBatchLanes);
  FFT32Rows(derivative.data(), correction_radius, correction_radius, kernels);
  const __m256i top_locator = _mm256_xor_si256(
      _mm256_load_si256(reinterpret_cast<const __m256i*>(
          Row(locator_second, correction_radius))),
      _mm256_load_si256(reinterpret_cast<const __m256i*>(Row(derivative, 0))));
  _mm256_store_si256(
      reinterpret_cast<__m256i*>(Row(locator_first, correction_radius)),
      top_locator);

  for (size_t lane = 0; lane < kBatchLanes; ++lane) {
    size_t degree = 0;
    bool nonzero = false;
    for (size_t coefficient = 0; coefficient <= correction_radius;
         ++coefficient) {
      if (Row(locator_first, coefficient)[lane] != 0) {
        nonzero = true;
        degree = coefficient;
      }
    }
    if (!nonzero || degree != locator_degrees[lane]) {
      candidate_lanes &= ~(uint32_t{1} << lane);
    }
  }

  std::fill_n(work.data(), codeword_size * kBatchLanes, Element{0});
  for (size_t block = 0; block < codeword_size; block += recovery_count) {
    std::memcpy(Row(work, block), locator_first.data(),
                (correction_radius + 1) * kBatchLanes);
  }
  FFT32Blocks(work, codeword_size, recovery_count, kernels);
  for (size_t position = 0; position < codeword_size; ++position) {
    const __m256i values = _mm256_load_si256(
        reinterpret_cast<const __m256i*>(Row(work, position)));
    root_masks[position] = ZeroMask(values) & candidate_lanes;
    uint32_t roots = root_masks[position];
    while (roots != 0) {
      const size_t lane = std::countr_zero(roots);
      ++root_counts[lane];
      roots &= roots - 1;
    }
  }
  for (size_t lane = 0; lane < kBatchLanes; ++lane) {
    if (root_counts[lane] != locator_degrees[lane]) {
      candidate_lanes &= ~(uint32_t{1} << lane);
    }
  }

  uint32_t data_error_lanes = 0;
  for (size_t native_position = 0; native_position < codeword_size;
       ++native_position) {
    if (PublicPosition(parameters.family, data_count, recovery_count,
                       native_position) < data_count) {
      data_error_lanes |= root_masks[native_position];
    }
  }
  const uint32_t location_only_lanes = candidate_lanes & ~data_error_lanes;
  candidate_lanes &= data_error_lanes;
  if (candidate_lanes == 0) {
    PublishChunkResults(results, error_masks, byte_count, column,
                        parameters.family, data_count, recovery_count,
                        clean_lanes, location_only_lanes, root_masks,
                        root_counts);
    return;
  }

  // Reuse the second locator row for evaluator coefficients. First construct
  // z(a)=u(a)lambda(a) in the AES-isomorphic basis.
  std::memcpy(derivative.data(), syndrome_or_scratch.data(),
              (correction_radius + 1) * kBatchLanes);
  lch::detail::ConvertCantorToAES(
      std::span(locator_second).first((correction_radius + 1) * kBatchLanes));
  lch::detail::ConvertCantorToAES(
      std::span(derivative).first((correction_radius + 1) * kBatchLanes));
  for (size_t i = 0; i <= correction_radius; ++i) {
    const __m256i locator = _mm256_load_si256(
        reinterpret_cast<const __m256i*>(Row(locator_second, i)));
    const __m256i syndrome =
        _mm256_load_si256(reinterpret_cast<const __m256i*>(Row(derivative, i)));
    _mm256_store_si256(reinterpret_cast<__m256i*>(Row(locator_second, i)),
                       _mm256_gf2p8mul_epi8(locator, syndrome));
  }
  lch::detail::ConvertAESToCantor(
      std::span(locator_second).first((correction_radius + 1) * kBatchLanes));

  IFFT32Rows(locator_second.data(), correction_radius, 0, kernels);
  std::memcpy(derivative.data(), locator_second.data(),
              correction_radius * kBatchLanes);
  FFT32Rows(derivative.data(), correction_radius, correction_radius, kernels);
  const __m256i top_evaluator = _mm256_xor_si256(
      _mm256_load_si256(reinterpret_cast<const __m256i*>(
          Row(locator_second, correction_radius))),
      _mm256_load_si256(reinterpret_cast<const __m256i*>(Row(derivative, 0))));
  _mm256_store_si256(
      reinterpret_cast<__m256i*>(Row(locator_second, correction_radius)),
      top_evaluator);
  for (size_t lane = 0; lane < kBatchLanes; ++lane) {
    size_t degree = 0;
    bool nonzero = false;
    for (size_t coefficient = 0; coefficient <= correction_radius;
         ++coefficient) {
      if (Row(locator_second, coefficient)[lane] != 0) {
        nonzero = true;
        degree = coefficient;
      }
    }
    if (nonzero && degree >= locator_degrees[lane]) {
      candidate_lanes &= ~(uint32_t{1} << lane);
    }
  }

  DifferentiateRows(locator_first.data(), correction_radius, derivative.data(),
                    recovery_count);
  for (size_t block = recovery_count; block < codeword_size;
       block += recovery_count) {
    std::memcpy(Row(derivative, block), derivative.data(),
                recovery_count * kBatchLanes);
  }
  FFT32Blocks(derivative, codeword_size, recovery_count, kernels);

  std::fill_n(work.data(), codeword_size * kBatchLanes, Element{0});
  for (size_t block = 0; block < codeword_size; block += recovery_count) {
    std::memcpy(Row(work, block), locator_second.data(),
                (correction_radius + 1) * kBatchLanes);
  }
  FFT32Blocks(work, codeword_size, recovery_count, kernels);

  const size_t syndrome_level = std::countr_zero(recovery_count);
  for (size_t position = recovery_count; position < codeword_size; ++position) {
    const uint32_t active = root_masks[position] & candidate_lanes;
    if (active == 0) {
      continue;
    }
    candidate_lanes &= ~DivideRows(
        Row(work, position), Row(derivative, position),
        context.Skew(syndrome_level, position), Row(work, position), active);
    const __m256i correction = _mm256_load_si256(
        reinterpret_cast<const __m256i*>(Row(work, position)));
    candidate_lanes &= ~(active & ZeroMask(correction));
  }

  // Inside the syndrome subspace use u(a)+z'(a)/lambda'(a).
  DifferentiateRows(locator_second.data(), correction_radius, work.data(),
                    recovery_count);
  FFT32Rows(work.data(), recovery_count, 0, kernels);
  for (size_t position = 0; position < recovery_count; ++position) {
    const uint32_t active = root_masks[position] & candidate_lanes;
    if (active == 0) {
      continue;
    }
    candidate_lanes &=
        ~DivideRows(Row(work, position), Row(derivative, position), 1,
                    Row(work, position), active);
    const __m256i ratio = _mm256_load_si256(
        reinterpret_cast<const __m256i*>(Row(work, position)));
    const __m256i syndrome = _mm256_load_si256(
        reinterpret_cast<const __m256i*>(Row(syndrome_or_scratch, position)));
    const __m256i correction = _mm256_xor_si256(ratio, syndrome);
    _mm256_store_si256(reinterpret_cast<__m256i*>(Row(work, position)),
                       correction);
    candidate_lanes &= ~(active & ZeroMask(correction));
  }

  // Verify all candidate codewords after correcting a private gathered copy.
  CopyNativeCodeword(data, recovery, parameters.family, column, derivative);
  for (size_t position = 0; position < codeword_size; ++position) {
    const uint32_t active = root_masks[position] & candidate_lanes;
    if (active == 0) {
      continue;
    }
    const __m256i received = _mm256_load_si256(
        reinterpret_cast<const __m256i*>(Row(derivative, position)));
    const __m256i correction = _mm256_and_si256(
        _mm256_load_si256(
            reinterpret_cast<const __m256i*>(Row(work, position))),
        LaneMask(active));
    _mm256_store_si256(reinterpret_cast<__m256i*>(Row(derivative, position)),
                       _mm256_xor_si256(received, correction));
  }
  IFFT32Blocks(derivative, codeword_size, recovery_count, kernels);
  const uint32_t verified_lanes =
      SyndromeZeroMask(derivative, codeword_size, recovery_count);
  candidate_lanes &= verified_lanes;

  for (size_t native_position = 0; native_position < codeword_size;
       ++native_position) {
    const size_t public_position = PublicPosition(
        parameters.family, data_count, recovery_count, native_position);
    if (public_position >= data_count) {
      continue;
    }
    const uint32_t active = root_masks[native_position] & candidate_lanes;
    uint32_t lanes = active;
    while (lanes != 0) {
      const size_t lane = std::countr_zero(lanes);
      data[public_position][column + lane] ^= Row(work, native_position)[lane];
      lanes &= lanes - 1;
    }
  }
  PublishChunkResults(results, error_masks, byte_count, column,
                      parameters.family, data_count, recovery_count,
                      clean_lanes, candidate_lanes | location_only_lanes,
                      root_masks, root_counts);
}

#endif

}  // namespace

CorrectionStatus CorrectBatch(const LCHDecoder& decoder,
                              std::span<Element* const> data,
                              std::span<const Element* const> recovery,
                              size_t byte_count,
                              std::span<CorrectionResult> results,
                              std::span<uint8_t> error_masks) {
  if (!decoder.Valid()) {
    return CorrectionStatus::invalid_argument;
  }
  const size_t data_count = decoder.DataCount();
  const size_t recovery_count = decoder.RecoveryCount();
  if (data_count > kFieldSize || recovery_count > kFieldSize - data_count ||
      byte_count >
          std::numeric_limits<size_t>::max() / sizeof(CorrectionResult) ||
      (byte_count != 0 &&
       data_count + recovery_count >
           std::numeric_limits<size_t>::max() / byte_count)) {
    return CorrectionStatus::invalid_argument;
  }
  const size_t codeword_size = data_count + recovery_count;
  if (data.size() != data_count || recovery.size() != recovery_count ||
      results.size() != byte_count ||
      error_masks.size() != codeword_size * byte_count) {
    return CorrectionStatus::invalid_argument;
  }

  CodeParameters parameters;
  if (!SupportedDimensions(decoder, parameters)) {
    return CorrectionStatus::unsupported_dimensions;
  }
  if (byte_count == 0) {
    return CorrectionStatus::ok;
  }

  std::array<ByteRange, kFieldSize + 2> ranges;
  size_t range_count = 0;
  for (Element* shard : data) {
    if (!AddRange(shard, byte_count, ranges, range_count)) {
      return CorrectionStatus::invalid_argument;
    }
  }
  for (const Element* shard : recovery) {
    if (!AddRange(shard, byte_count, ranges, range_count)) {
      return CorrectionStatus::invalid_argument;
    }
  }
  if (!AddRange(results.data(), results.size_bytes(), ranges, range_count) ||
      !AddRange(error_masks.data(), error_masks.size_bytes(), ranges,
                range_count) ||
      !RangesAreDisjoint(std::span(ranges).first(range_count))) {
    return CorrectionStatus::invalid_argument;
  }

  std::fill(error_masks.begin(), error_masks.end(), uint8_t{0});
  std::fill(results.begin(), results.end(), CorrectionResult{});

  size_t column = 0;
#if defined(__GFNI__) && defined(__AVX2__)
  if (lch::BackendAvailable(lch::Backend::gfni256_affine)) {
    const lch::Backend transform_backend = lch::SelectBackend(kBatchLanes);
    const lch::detail::ResolvedKernels* kernels =
        lch::detail::ResolveKernels(transform_backend, kBatchLanes);
    if (kernels != nullptr) {
      for (; column + kBatchLanes <= byte_count; column += kBatchLanes) {
        CorrectChunk32(data, recovery, byte_count, column, parameters, *kernels,
                       results, error_masks);
      }
    }
  }
#endif
  return CorrectColumnsScalar(decoder, data, recovery, byte_count, column,
                              results, error_masks);
}

}  // namespace gf2p8::rs::detail::error_correction
