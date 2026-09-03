#include "reed_solomon/experiment/gfni512_radix8.h"

#include <algorithm>
#include <array>

namespace gf2p8::rs::detail::experiment::radix8 {

lch::Status EncodeLCH(const lch::Context& context,
                      std::span<const Element* const> data,
                      std::span<Element* const> recovery,
                      size_t byte_count,
                      size_t transform_size,
                      std::span<Element> workspace,
                      const lch::detail::ResolvedKernels& base,
                      const lch::detail::experiment::radix8::Kernels& radix8) {
  if (byte_count == 0) {
    return lch::Status::ok;
  }
  if (transform_size == 1) {
    std::copy_n(data[0], byte_count, recovery[0]);
    for (size_t i = 1; i < data.size(); ++i) {
      base.xor_one(recovery[0], data[i], byte_count);
    }
    return lch::Status::ok;
  }

  std::array<Element*, lch::Context::kFieldSize> accumulator_pointers;
  std::array<Element*, lch::Context::kFieldSize> scratch_pointers;
  Element* const workspace_base = workspace.data();
  auto accumulator = recovery;
  if (recovery.size() != transform_size) {
    for (size_t i = 0; i < recovery.size(); ++i) {
      accumulator_pointers[i] = recovery[i];
    }
    for (size_t i = recovery.size(); i < transform_size; ++i) {
      accumulator_pointers[i] =
          workspace_base + (i - recovery.size()) * byte_count;
    }
    accumulator =
        std::span<Element* const>(accumulator_pointers).first(transform_size);
  }
  Element* const scratch_base =
      workspace_base + (transform_size - recovery.size()) * byte_count;
  for (size_t i = 0; i < transform_size; ++i) {
    scratch_pointers[i] = scratch_base + i * byte_count;
  }
  auto scratch =
      std::span<Element* const>(scratch_pointers).first(transform_size);

  size_t source_offset = 0;
  size_t live_count = std::min(transform_size, data.size());
  lch::Status status = lch::detail::experiment::radix8::IFFTResolved(
      context, data.first(live_count), accumulator, byte_count, transform_size,
      base, radix8);
  if (status != lch::Status::ok) {
    return status;
  }
  source_offset += live_count;
  while (source_offset < data.size()) {
    live_count = std::min(transform_size, data.size() - source_offset);
    status = lch::detail::experiment::radix8::IFFTResolved(
        context, data.subspan(source_offset, live_count), scratch, accumulator,
        byte_count, transform_size + source_offset, base, radix8);
    if (status != lch::Status::ok) {
      return status;
    }
    source_offset += live_count;
  }
  return lch::detail::FFTResolved(context, accumulator, byte_count, 0,
                                  recovery.size(), base, lch::Radix::radix4);
}

}  // namespace gf2p8::rs::detail::experiment::radix8
