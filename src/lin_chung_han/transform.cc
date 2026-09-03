#include <algorithm>

#include "lin_chung_han/transform_internal.h"

namespace gf2p8::lch {
namespace {

bool IsPowerOfTwo(size_t value) {
  return value != 0 && (value & (value - 1)) == 0;
}

constexpr Element ProductOverSubspace(Element value,
                                      const std::array<Element, 256>& points,
                                      size_t dimension) {
  Element result = 1;
  const size_t count = size_t{1} << dimension;
  for (size_t i = 0; i < count; ++i) {
    result = gf2p8::detail::MultiplyCantorDirect(result, value ^ points[i]);
  }
  return result;
}

constexpr Element Power(Element value, unsigned exponent) {
  Element result = 1;
  while (exponent != 0) {
    if ((exponent & 1U) != 0) {
      result = gf2p8::detail::MultiplyCantorDirect(result, value);
    }
    value = gf2p8::detail::MultiplyCantorDirect(value, value);
    exponent >>= 1;
  }
  return result;
}

constexpr void Walsh(std::array<Element, 256>& values) {
  for (size_t width = 1; width < values.size(); width *= 2) {
    for (size_t block = 0; block < values.size(); block += 2 * width) {
      for (size_t i = 0; i < width; ++i) {
        const unsigned a = values[block + i];
        const unsigned b = values[block + width + i];
        values[block + i] = static_cast<Element>((a + b) % 255);
        values[block + width + i] = static_cast<Element>((a + 255 - b) % 255);
      }
    }
  }
}

struct DomainTables {
  std::array<Element, Context::kFieldBits> basis{};
  std::array<Element, Context::kFieldSize> points{};
  std::array<std::array<Element, Context::kFieldSize>, Context::kFieldBits>
      skew{};
  std::array<Element, Context::kFieldSize> log_walsh{};
};

consteval DomainTables MakeDomainTables() {
  DomainTables tables;
  for (size_t bit = 0; bit < Context::kFieldBits; ++bit) {
    tables.basis[bit] = static_cast<Element>(1U << bit);
  }
  for (size_t index = 0; index < Context::kFieldSize; ++index) {
    tables.points[index] = static_cast<Element>(index);
  }

  std::array<Element, Context::kFieldSize> exponent{};
  std::array<Element, Context::kFieldSize> logarithm{};
  logarithm.fill(255);
  const Element primitive = gf2p8::StandardToCantor(2);
  exponent[0] = 1;
  for (size_t i = 1; i < 255; ++i) {
    exponent[i] =
        gf2p8::detail::MultiplyCantorDirect(exponent[i - 1], primitive);
  }
  for (size_t i = 0; i < 255; ++i) {
    logarithm[exponent[i]] = static_cast<Element>(i);
  }
  for (size_t i = 0; i < Context::kFieldSize; ++i) {
    tables.log_walsh[i] = i == 0 ? 0 : logarithm[tables.points[i]];
  }
  Walsh(tables.log_walsh);

  for (size_t level = 0; level < Context::kFieldBits; ++level) {
    const Element denominator =
        ProductOverSubspace(tables.basis[level], tables.points, level);
    const Element denominator_inverse = Power(denominator, 254);
    for (size_t index = 0; index < Context::kFieldSize; ++index) {
      tables.skew[level][index] = gf2p8::detail::MultiplyCantorDirect(
          ProductOverSubspace(tables.points[index], tables.points, level),
          denominator_inverse);
    }
  }
  return tables;
}

constinit const DomainTables kDomainTables = MakeDomainTables();

Status Validate(std::span<Element* const> shards,
                size_t byte_count,
                size_t offset,
                TransformOptions options) {
  if (!IsPowerOfTwo(shards.size()) || shards.size() > Context::kFieldSize ||
      offset >= Context::kFieldSize ||
      offset + shards.size() > Context::kFieldSize ||
      (offset & (shards.size() - 1)) != 0) {
    return Status::invalid_argument;
  }
  if (byte_count != 0 &&
      std::any_of(shards.begin(), shards.end(),
                  [](const Element* shard) { return shard == nullptr; })) {
    return Status::invalid_argument;
  }
  if (options.backend != Backend::tuned && !BackendAvailable(options.backend)) {
    return Status::unsupported_backend;
  }
  if (options.radix != Radix::radix2 && options.radix != Radix::radix4) {
    return Status::invalid_argument;
  }
  return Status::ok;
}

Status ValidateInput(std::span<const Element* const> input,
                     std::span<Element* const> work,
                     size_t byte_count,
                     size_t offset,
                     TransformOptions options) {
  const Status status = Validate(work, byte_count, offset, options);
  if (status != Status::ok || input.size() > work.size()) {
    return status == Status::ok ? Status::invalid_argument : status;
  }
  if (byte_count != 0 &&
      std::any_of(input.begin(), input.end(),
                  [](const Element* shard) { return shard == nullptr; })) {
    return Status::invalid_argument;
  }
  return Status::ok;
}

size_t Log2(size_t value) {
  size_t result = 0;
  while (value > 1) {
    value >>= 1;
    ++result;
  }
  return result;
}

void RecordRadix2(testing::ScheduleStats* stats,
                  size_t kernel_calls,
                  bool fused = false) {
  if (stats == nullptr) {
    return;
  }
  ++stats->radix2_groups;
  stats->radix2_kernel_calls += kernel_calls;
  stats->shard_kernel_visits += 2 * kernel_calls;
  if (fused) {
    ++stats->fused_accumulation_groups;
    stats->fused_accumulation_kernel_calls += kernel_calls;
    stats->temporary_shard_stores_avoided += 2 * kernel_calls;
  }
}

void RecordRadix4(testing::ScheduleStats* stats,
                  size_t kernel_calls,
                  bool fused = false) {
  if (stats == nullptr) {
    return;
  }
  ++stats->radix4_groups;
  stats->radix4_kernel_calls += kernel_calls;
  stats->shard_kernel_visits += 4 * kernel_calls;
  if (fused) {
    ++stats->fused_accumulation_groups;
    stats->fused_accumulation_kernel_calls += kernel_calls;
    stats->temporary_shard_stores_avoided += 4 * kernel_calls;
  }
}

template <bool CollectStats = false,
          typename GroupNeeded,
          typename... StatsArgs>
Status RunFFTImpl(const Context& context,
                  std::span<Element* const> shards,
                  size_t byte_count,
                  size_t offset,
                  const detail::ResolvedKernels& kernels,
                  Radix radix,
                  GroupNeeded group_needed,
                  StatsArgs... stats_args) {
  static_assert(sizeof...(StatsArgs) == (CollectStats ? 1 : 0));
  testing::ScheduleStats* stats = nullptr;
  if constexpr (CollectStats) {
    stats = (stats_args, ...);
  }
  if (shards.size() == 1 || byte_count == 0) {
    return Status::ok;
  }

  if (radix == Radix::radix2) {
    for (size_t group_size = shards.size(); group_size >= 2; group_size /= 2) {
      const size_t half = group_size / 2;
      const size_t level = Log2(half);
      for (size_t block = 0; block < shards.size(); block += group_size) {
        if (!group_needed(block, group_size)) {
          continue;
        }
        if constexpr (CollectStats) {
          RecordRadix2(stats, half);
        }
        const Element coefficient = context.Skew(level, offset ^ block);
        for (size_t i = 0; i < half; ++i) {
          kernels.fft_radix2(shards[block + i], shards[block + half + i],
                             byte_count, coefficient, context.Tables());
        }
      }
    }
    return Status::ok;
  }

  size_t group_size = shards.size();
  for (; group_size >= 4; group_size /= 4) {
    const size_t distance = group_size / 4;
    const size_t top_level = Log2(group_size) - 1;
    const size_t low_level = top_level - 1;
    for (size_t block = 0; block < shards.size(); block += group_size) {
      if (!group_needed(block, group_size)) {
        continue;
      }
      if constexpr (CollectStats) {
        RecordRadix4(stats, distance);
      }
      const Element top = context.Skew(top_level, offset ^ block);
      const Element low = context.Skew(low_level, offset ^ block);
      const Element high =
          context.Skew(low_level, offset ^ (block + 2 * distance));
      for (size_t i = 0; i < distance; ++i) {
        kernels.fft_radix4(shards[block + i], shards[block + distance + i],
                           shards[block + 2 * distance + i],
                           shards[block + 3 * distance + i], byte_count, top,
                           low, high, context.Tables());
      }
    }
  }

  if (group_size == 2) {
    for (size_t block = 0; block < shards.size(); block += 2) {
      if (!group_needed(block, 2)) {
        continue;
      }
      if constexpr (CollectStats) {
        RecordRadix2(stats, 1);
      }
      kernels.fft_radix2(shards[block], shards[block + 1], byte_count,
                         context.Skew(0, offset ^ block), context.Tables());
    }
  }
  return Status::ok;
}

template <typename GroupNeeded>
Status RunFFT(const Context& context,
              std::span<Element* const> shards,
              size_t byte_count,
              size_t offset,
              const detail::ResolvedKernels& kernels,
              Radix radix,
              GroupNeeded group_needed) {
  return RunFFTImpl(context, shards, byte_count, offset, kernels, radix,
                    group_needed);
}

template <bool CollectStats = false, typename... StatsArgs>
Status RunIFFTImpl(const Context& context,
                   std::span<Element* const> shards,
                   size_t byte_count,
                   size_t offset,
                   size_t input_count,
                   const detail::ResolvedKernels& kernels,
                   Radix radix,
                   size_t initial_distance,
                   std::span<Element* const> xor_accumulator = {},
                   StatsArgs... stats_args) {
  static_assert(sizeof...(StatsArgs) == (CollectStats ? 1 : 0));
  testing::ScheduleStats* stats = nullptr;
  if constexpr (CollectStats) {
    stats = (stats_args, ...);
  }
  if (input_count == 0 || byte_count == 0) {
    return Status::ok;
  }
  if (shards.size() == 1) {
    if constexpr (CollectStats) {
      if (!xor_accumulator.empty()) {
        stats->fused_accumulation_groups = 1;
        stats->fused_accumulation_kernel_calls = 1;
        stats->temporary_shard_stores_avoided = 1;
        stats->shard_kernel_visits = 1;
      }
    }
    if (!xor_accumulator.empty()) {
      kernels.xor_one(xor_accumulator[0], shards[0], byte_count);
    }
    return Status::ok;
  }

  if (radix == Radix::radix2) {
    for (size_t half = initial_distance; half < shards.size(); half *= 2) {
      const size_t group_size = 2 * half;
      const size_t level = Log2(half);
      const bool final_layer = group_size == shards.size();
      for (size_t block = 0; block < input_count; block += group_size) {
        const bool fused = final_layer && !xor_accumulator.empty();
        if constexpr (CollectStats) {
          RecordRadix2(stats, half, fused);
        }
        const Element coefficient = context.Skew(level, offset ^ block);
        for (size_t i = 0; i < half; ++i) {
          if (fused) {
            kernels.ifft_radix2_xor(shards[block + i], shards[block + half + i],
                                    xor_accumulator[block + i],
                                    xor_accumulator[block + half + i],
                                    byte_count, coefficient, context.Tables());
          } else {
            kernels.ifft_radix2(shards[block + i], shards[block + half + i],
                                byte_count, coefficient, context.Tables());
          }
        }
      }
    }
    return Status::ok;
  }

  size_t distance = initial_distance;
  for (; 4 * distance <= shards.size(); distance *= 4) {
    const size_t group_size = 4 * distance;
    const size_t low_level = Log2(distance);
    const size_t top_level = low_level + 1;
    const bool final_layer = group_size == shards.size();
    for (size_t block = 0; block < input_count; block += group_size) {
      const bool fused = final_layer && !xor_accumulator.empty();
      if constexpr (CollectStats) {
        RecordRadix4(stats, distance, fused);
      }
      const Element top = context.Skew(top_level, offset ^ block);
      const Element low = context.Skew(low_level, offset ^ block);
      const Element high =
          context.Skew(low_level, offset ^ (block + 2 * distance));
      for (size_t i = 0; i < distance; ++i) {
        if (fused) {
          kernels.ifft_radix4_xor(
              shards[block + i], shards[block + distance + i],
              shards[block + 2 * distance + i],
              shards[block + 3 * distance + i], xor_accumulator[block + i],
              xor_accumulator[block + distance + i],
              xor_accumulator[block + 2 * distance + i],
              xor_accumulator[block + 3 * distance + i], byte_count, top, low,
              high, context.Tables());
        } else {
          kernels.ifft_radix4(shards[block + i], shards[block + distance + i],
                              shards[block + 2 * distance + i],
                              shards[block + 3 * distance + i], byte_count, top,
                              low, high, context.Tables());
        }
      }
    }
  }

  if (distance < shards.size()) {
    const size_t half = distance;
    if constexpr (CollectStats) {
      RecordRadix2(stats, half, !xor_accumulator.empty());
    }
    const Element coefficient = context.Skew(Log2(half), offset);
    for (size_t i = 0; i < half; ++i) {
      if (xor_accumulator.empty()) {
        kernels.ifft_radix2(shards[i], shards[half + i], byte_count,
                            coefficient, context.Tables());
      } else {
        kernels.ifft_radix2_xor(shards[i], shards[half + i], xor_accumulator[i],
                                xor_accumulator[half + i], byte_count,
                                coefficient, context.Tables());
      }
    }
  }
  return Status::ok;
}

Status RunIFFT(const Context& context,
               std::span<Element* const> shards,
               size_t byte_count,
               size_t offset,
               size_t input_count,
               const detail::ResolvedKernels& kernels,
               Radix radix,
               std::span<Element* const> xor_accumulator = {},
               size_t initial_distance = 1) {
  return RunIFFTImpl(context, shards, byte_count, offset, input_count, kernels,
                     radix, initial_distance, xor_accumulator);
}

Status RunIFFTCopyFirst(const Context& context,
                        std::span<const Element* const> input,
                        std::span<Element* const> work,
                        size_t byte_count,
                        size_t offset,
                        const detail::ResolvedKernels& kernels,
                        Radix radix,
                        std::span<Element* const> xor_accumulator = {}) {
  if (byte_count == 0) {
    return Status::ok;
  }
  if (work.size() == 1) {
    if (xor_accumulator.empty()) {
      std::copy_n(input[0], byte_count, work[0]);
    } else {
      kernels.xor_one(xor_accumulator[0], input[0], byte_count);
    }
    return Status::ok;
  }

  if (radix == Radix::radix2 || work.size() == 2) {
    const bool fused = work.size() == 2 && !xor_accumulator.empty();
    for (size_t block = 0; block < work.size(); block += 2) {
      const Element coefficient = context.Skew(0, offset ^ block);
      if (fused) {
        kernels.ifft_radix2_xor(input[block], input[block + 1],
                                xor_accumulator[block],
                                xor_accumulator[block + 1], byte_count,
                                coefficient, context.Tables());
      } else {
        const Status status = kernels.ifft_radix2_copy(
            input[block], input[block + 1], work[block], work[block + 1],
            byte_count, coefficient, context.Tables());
        if (status != Status::ok) {
          return status;
        }
      }
    }
    if (work.size() == 2) {
      return Status::ok;
    }
    return RunIFFT(context, work, byte_count, offset, work.size(), kernels,
                   radix, xor_accumulator, 2);
  }

  const bool fused = work.size() == 4 && !xor_accumulator.empty();
  for (size_t block = 0; block < work.size(); block += 4) {
    const Element top = context.Skew(1, offset ^ block);
    const Element low = context.Skew(0, offset ^ block);
    const Element high = context.Skew(0, offset ^ (block + 2));
    if (fused) {
      kernels.ifft_radix4_xor(
          input[block], input[block + 1], input[block + 2], input[block + 3],
          xor_accumulator[block], xor_accumulator[block + 1],
          xor_accumulator[block + 2], xor_accumulator[block + 3], byte_count,
          top, low, high, context.Tables());
    } else {
      const Status status = kernels.ifft_radix4_copy(
          input[block], input[block + 1], input[block + 2], input[block + 3],
          work[block], work[block + 1], work[block + 2], work[block + 3],
          byte_count, top, low, high, context.Tables());
      if (status != Status::ok) {
        return status;
      }
    }
  }
  if (work.size() == 4) {
    return Status::ok;
  }
  return RunIFFT(context, work, byte_count, offset, work.size(), kernels, radix,
                 xor_accumulator, 4);
}

Backend ResolveBackend(Backend backend, size_t byte_count) {
  return backend == Backend::tuned ? SelectBackend(byte_count) : backend;
}

const detail::ResolvedKernels& ResolveKernelSet(Backend backend,
                                                size_t byte_count) {
  return *detail::ResolveKernels(ResolveBackend(backend, byte_count),
                                 byte_count);
}

}  // namespace

const Context& Context::Shared() {
  static const Context context;
  return context;
}

const std::array<Element, Context::kFieldBits>& Context::Basis() const {
  return kDomainTables.basis;
}

Element Context::EvaluationPoint(size_t index) const {
  return kDomainTables.points[index];
}

Element Context::Skew(size_t level, size_t index) const {
  return kDomainTables.skew[level][index];
}

Element Context::Log(Element value) const {
  return Tables().cantor.logarithm[value];
}

Element Context::Exp(size_t exponent) const {
  return Tables().cantor.exponent[exponent % 255];
}

Element Context::LogWalsh(size_t index) const {
  return kDomainTables.log_walsh[index];
}

Status detail::FFTResolved(const Context& context,
                           std::span<Element* const> shards,
                           size_t byte_count,
                           size_t evaluation_offset,
                           size_t output_count,
                           const ResolvedKernels& kernels,
                           Radix radix) {
  return RunFFT(
      context, shards, byte_count, evaluation_offset, kernels, radix,
      [output_count](size_t block, size_t) { return block < output_count; });
}

Status detail::FFTResolved(const Context& context,
                           std::span<Element* const> shards,
                           size_t byte_count,
                           size_t evaluation_offset,
                           std::span<const uint8_t> requested_outputs,
                           const ResolvedKernels& kernels,
                           Radix radix) {
  return RunFFT(context, shards, byte_count, evaluation_offset, kernels, radix,
                [requested_outputs](size_t block, size_t group_size) {
                  return std::any_of(
                      requested_outputs.begin() + block,
                      requested_outputs.begin() + block + group_size,
                      [](uint8_t requested) { return requested != 0; });
                });
}

Status detail::IFFTResolved(const Context& context,
                            std::span<Element* const> shards,
                            size_t byte_count,
                            size_t evaluation_offset,
                            size_t input_count,
                            const ResolvedKernels& kernels,
                            Radix radix) {
  return RunIFFT(context, shards, byte_count, evaluation_offset, input_count,
                 kernels, radix);
}

Status detail::IFFTResolved(const Context& context,
                            std::span<const Element* const> input,
                            std::span<Element* const> work,
                            size_t byte_count,
                            size_t evaluation_offset,
                            const ResolvedKernels& kernels,
                            Radix radix) {
  if (input.size() == work.size() && kernels.ifft_radix2_copy != nullptr &&
      kernels.ifft_radix4_copy != nullptr) {
    return RunIFFTCopyFirst(context, input, work, byte_count, evaluation_offset,
                            kernels, radix);
  }
  for (size_t i = 0; i < input.size(); ++i) {
    std::copy_n(input[i], byte_count, work[i]);
  }
  for (size_t i = input.size(); i < work.size(); ++i) {
    std::fill_n(work[i], byte_count, 0);
  }
  return RunIFFT(context, work, byte_count, evaluation_offset, input.size(),
                 kernels, radix);
}

Status detail::IFFTResolved(const Context& context,
                            std::span<const Element* const> input,
                            std::span<Element* const> work,
                            std::span<Element* const> xor_accumulator,
                            size_t byte_count,
                            size_t evaluation_offset,
                            const ResolvedKernels& kernels,
                            Radix radix) {
  if (input.size() == work.size() && kernels.ifft_radix2_copy != nullptr &&
      kernels.ifft_radix4_copy != nullptr) {
    return RunIFFTCopyFirst(context, input, work, byte_count, evaluation_offset,
                            kernels, radix, xor_accumulator);
  }
  for (size_t i = 0; i < input.size(); ++i) {
    std::copy_n(input[i], byte_count, work[i]);
  }
  for (size_t i = input.size(); i < work.size(); ++i) {
    std::fill_n(work[i], byte_count, 0);
  }
  if (input.empty()) {
    return Status::ok;
  }
  if (work.size() == 1) {
    kernels.xor_one(xor_accumulator[0], work[0], byte_count);
    return Status::ok;
  }
  return RunIFFT(context, work, byte_count, evaluation_offset, input.size(),
                 kernels, radix, xor_accumulator);
}

Status FFT(const Context& context,
           std::span<Element* const> shards,
           size_t byte_count,
           size_t evaluation_offset,
           TransformOptions options) {
  return FFT(context, shards, byte_count, evaluation_offset, shards.size(),
             options);
}

Status FFT(const Context& context,
           std::span<Element* const> shards,
           size_t byte_count,
           size_t evaluation_offset,
           size_t output_count,
           TransformOptions options) {
  const Status validation =
      Validate(shards, byte_count, evaluation_offset, options);
  if (validation != Status::ok || output_count > shards.size()) {
    return validation == Status::ok ? Status::invalid_argument : validation;
  }
  return detail::FFTResolved(
      context, shards, byte_count, evaluation_offset, output_count,
      ResolveKernelSet(options.backend, byte_count), options.radix);
}

Status FFT(const Context& context,
           std::span<Element* const> shards,
           size_t byte_count,
           size_t evaluation_offset,
           std::span<const uint8_t> requested_outputs,
           TransformOptions options) {
  const Status validation =
      Validate(shards, byte_count, evaluation_offset, options);
  if (validation != Status::ok || requested_outputs.size() != shards.size()) {
    return validation == Status::ok ? Status::invalid_argument : validation;
  }
  return detail::FFTResolved(
      context, shards, byte_count, evaluation_offset, requested_outputs,
      ResolveKernelSet(options.backend, byte_count), options.radix);
}

Status IFFT(const Context& context,
            std::span<Element* const> shards,
            size_t byte_count,
            size_t evaluation_offset,
            TransformOptions options) {
  return IFFT(context, shards, byte_count, evaluation_offset, shards.size(),
              options);
}

Status IFFT(const Context& context,
            std::span<Element* const> shards,
            size_t byte_count,
            size_t evaluation_offset,
            size_t input_count,
            TransformOptions options) {
  const Status validation =
      Validate(shards, byte_count, evaluation_offset, options);
  if (validation != Status::ok || input_count > shards.size()) {
    return validation == Status::ok ? Status::invalid_argument : validation;
  }
  return detail::IFFTResolved(
      context, shards, byte_count, evaluation_offset, input_count,
      ResolveKernelSet(options.backend, byte_count), options.radix);
}

Status IFFT(const Context& context,
            std::span<const Element* const> input,
            std::span<Element* const> work,
            size_t byte_count,
            size_t evaluation_offset,
            TransformOptions options) {
  const Status validation =
      ValidateInput(input, work, byte_count, evaluation_offset, options);
  if (validation != Status::ok) {
    return validation;
  }
  return detail::IFFTResolved(
      context, input, work, byte_count, evaluation_offset,
      ResolveKernelSet(options.backend, byte_count), options.radix);
}

Status IFFT(const Context& context,
            std::span<const Element* const> input,
            std::span<Element* const> work,
            std::span<Element* const> xor_accumulator,
            size_t byte_count,
            size_t evaluation_offset,
            TransformOptions options) {
  const Status validation =
      ValidateInput(input, work, byte_count, evaluation_offset, options);
  if (validation != Status::ok || xor_accumulator.size() != work.size() ||
      (byte_count != 0 &&
       std::any_of(xor_accumulator.begin(), xor_accumulator.end(),
                   [](const Element* shard) { return shard == nullptr; }))) {
    return validation == Status::ok ? Status::invalid_argument : validation;
  }
  return detail::IFFTResolved(
      context, input, work, xor_accumulator, byte_count, evaluation_offset,
      ResolveKernelSet(options.backend, byte_count), options.radix);
}

namespace testing {
namespace {

std::span<Element* const> ScheduleShards(
    size_t shard_count,
    std::array<Element, Context::kFieldSize>& storage,
    std::array<Element*, Context::kFieldSize>& pointers) {
  for (size_t i = 0; i < shard_count; ++i) {
    pointers[i] = &storage[i];
  }
  return std::span<Element* const>(pointers).first(shard_count);
}

}  // namespace

ScheduleStats FFTSchedule(size_t shard_count,
                          size_t output_count,
                          Radix radix) {
  if (output_count > shard_count) {
    return {};
  }
  ScheduleStats stats;
  if (!IsPowerOfTwo(shard_count) || shard_count > Context::kFieldSize ||
      (radix != Radix::radix2 && radix != Radix::radix4)) {
    return stats;
  }
  const Context& context = Context::Shared();
  std::array<Element, Context::kFieldSize> storage{};
  std::array<Element*, Context::kFieldSize> pointers{};
  const auto shards = ScheduleShards(shard_count, storage, pointers);
  const auto& kernels = *detail::ResolveKernels(Backend::scalar);
  static_cast<void>(RunFFTImpl<true>(
      context, shards, 1, 0, kernels, radix,
      [output_count](size_t block, size_t) { return block < output_count; },
      &stats));
  return stats;
}

ScheduleStats FFTSchedule(size_t shard_count,
                          std::span<const uint8_t> requested_outputs,
                          Radix radix) {
  if (requested_outputs.size() != shard_count) {
    return {};
  }
  ScheduleStats stats;
  if (!IsPowerOfTwo(shard_count) || shard_count > Context::kFieldSize ||
      (radix != Radix::radix2 && radix != Radix::radix4)) {
    return stats;
  }
  const Context& context = Context::Shared();
  std::array<Element, Context::kFieldSize> storage{};
  std::array<Element*, Context::kFieldSize> pointers{};
  const auto shards = ScheduleShards(shard_count, storage, pointers);
  const auto& kernels = *detail::ResolveKernels(Backend::scalar);
  static_cast<void>(RunFFTImpl<true>(
      context, shards, 1, 0, kernels, radix,
      [requested_outputs](size_t block, size_t group_size) {
        return std::any_of(requested_outputs.begin() + block,
                           requested_outputs.begin() + block + group_size,
                           [](uint8_t requested) { return requested != 0; });
      },
      &stats));
  return stats;
}

ScheduleStats IFFTSchedule(size_t shard_count,
                           size_t input_count,
                           Radix radix,
                           bool fused_accumulation) {
  ScheduleStats stats;
  if (!IsPowerOfTwo(shard_count) || shard_count > Context::kFieldSize ||
      input_count > shard_count || input_count == 0) {
    return stats;
  }
  if (radix != Radix::radix2 && radix != Radix::radix4) {
    return stats;
  }
  const Context& context = Context::Shared();
  std::array<Element, Context::kFieldSize> storage{};
  std::array<Element*, Context::kFieldSize> pointers{};
  std::array<Element, Context::kFieldSize> accumulator_storage{};
  std::array<Element*, Context::kFieldSize> accumulator_pointers{};
  const auto shards = ScheduleShards(shard_count, storage, pointers);
  const auto accumulator =
      ScheduleShards(shard_count, accumulator_storage, accumulator_pointers);
  const auto& kernels = *detail::ResolveKernels(Backend::scalar);
  static_cast<void>(RunIFFTImpl<true>(
      context, shards, 1, 0, input_count, kernels, radix, 1,
      fused_accumulation ? accumulator : std::span<Element* const>{}, &stats));
  return stats;
}

}  // namespace testing

}  // namespace gf2p8::lch
