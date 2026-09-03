#pragma once

#include <cstddef>
#include <span>

#include "lin_chung_han/kernels_internal.h"
#include "lin_chung_han/transform.h"

namespace gf2p8::lch::detail {

Status FFTResolved(const Context& context,
                   std::span<Element* const> shards,
                   size_t byte_count,
                   size_t evaluation_offset,
                   size_t output_count,
                   const ResolvedKernels& kernels,
                   Radix radix);
Status FFTResolved(const Context& context,
                   std::span<Element* const> shards,
                   size_t byte_count,
                   size_t evaluation_offset,
                   std::span<const uint8_t> requested_outputs,
                   const ResolvedKernels& kernels,
                   Radix radix);
Status IFFTResolved(const Context& context,
                    std::span<Element* const> shards,
                    size_t byte_count,
                    size_t evaluation_offset,
                    size_t input_count,
                    const ResolvedKernels& kernels,
                    Radix radix);
Status IFFTResolved(const Context& context,
                    std::span<const Element* const> input,
                    std::span<Element* const> work,
                    size_t byte_count,
                    size_t evaluation_offset,
                    const ResolvedKernels& kernels,
                    Radix radix);
Status IFFTResolved(const Context& context,
                    std::span<const Element* const> input,
                    std::span<Element* const> work,
                    std::span<Element* const> xor_accumulator,
                    size_t byte_count,
                    size_t evaluation_offset,
                    const ResolvedKernels& kernels,
                    Radix radix);

}  // namespace gf2p8::lch::detail

namespace gf2p8::lch::testing {

struct ScheduleStats {
  size_t radix2_groups = 0;
  size_t radix4_groups = 0;
  size_t radix2_kernel_calls = 0;
  size_t radix4_kernel_calls = 0;
  size_t shard_kernel_visits = 0;
  size_t fused_accumulation_groups = 0;
  size_t fused_accumulation_kernel_calls = 0;
  size_t temporary_shard_stores_avoided = 0;
};

ScheduleStats FFTSchedule(size_t shard_count, size_t output_count, Radix radix);
ScheduleStats FFTSchedule(size_t shard_count,
                          std::span<const uint8_t> requested_outputs,
                          Radix radix);
ScheduleStats IFFTSchedule(size_t shard_count,
                           size_t input_count,
                           Radix radix,
                           bool fused_accumulation);

}  // namespace gf2p8::lch::testing
