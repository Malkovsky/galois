#include "reed_solomon/common.h"

#include <algorithm>
#include <array>
#include <limits>

namespace gf2p8::rs::detail {

size_t ShardWorkspaceSize(size_t shard_count,
                          size_t byte_count,
                          size_t prefix_bytes) {
  if (shard_count == 0) {
    return prefix_bytes;
  }
  if (byte_count >
      (std::numeric_limits<size_t>::max() - prefix_bytes) / shard_count) {
    return std::numeric_limits<size_t>::max();
  }
  return prefix_bytes + shard_count * byte_count;
}

lch::Status ValidateBackend(lch::Backend backend) {
  if (backend != lch::Backend::tuned && !lch::BackendAvailable(backend)) {
    return lch::Status::unsupported_backend;
  }
  return lch::Status::ok;
}

lch::Backend ResolveBackend(lch::Backend backend, size_t byte_count) {
  return backend == lch::Backend::tuned ? lch::SelectBackend(byte_count)
                                        : backend;
}

lch::Status EncodeLCH(const lch::Context& context,
                      std::span<const Element* const> data,
                      std::span<Element* const> recovery,
                      size_t byte_count,
                      size_t transform_size,
                      std::span<Element> workspace,
                      const lch::detail::ResolvedKernels& kernels,
                      lch::Radix radix) {
  if (byte_count == 0) {
    return lch::Status::ok;
  }
  if (transform_size == 1) {
    std::copy_n(data[0], byte_count, recovery[0]);
    for (size_t i = 1; i < data.size(); ++i) {
      kernels.xor_one(recovery[0], data[i], byte_count);
    }
    return lch::Status::ok;
  }

  // Recovery storage backs the requested accumulator prefix. Workspace backs
  // only the padded prefix and one temporary IFFT block.
  std::array<Element*, lch::Context::kFieldSize> accumulator_pointers;
  std::array<Element*, lch::Context::kFieldSize> scratch_pointers;
  Element* const base = workspace.data();
  auto accumulator = recovery;
  if (recovery.size() != transform_size) {
    for (size_t i = 0; i < recovery.size(); ++i) {
      accumulator_pointers[i] = recovery[i];
    }
    for (size_t i = recovery.size(); i < transform_size; ++i) {
      accumulator_pointers[i] = base + (i - recovery.size()) * byte_count;
    }
    accumulator =
        std::span<Element* const>(accumulator_pointers).first(transform_size);
  }
  Element* const scratch_base =
      base + (transform_size - recovery.size()) * byte_count;
  for (size_t i = 0; i < transform_size; ++i) {
    scratch_pointers[i] = scratch_base + i * byte_count;
  }
  auto scratch =
      std::span<Element* const>(scratch_pointers).first(transform_size);

  size_t source_offset = 0;
  size_t live_count = std::min(transform_size, data.size());
  lch::Status status =
      lch::detail::IFFTResolved(context, data.first(live_count), accumulator,
                                byte_count, transform_size, kernels, radix);
  if (status != lch::Status::ok) {
    return status;
  }
  source_offset += live_count;
  while (source_offset < data.size()) {
    live_count = std::min(transform_size, data.size() - source_offset);
    status = lch::detail::IFFTResolved(
        context, data.subspan(source_offset, live_count), scratch, accumulator,
        byte_count, transform_size + source_offset, kernels, radix);
    if (status != lch::Status::ok) {
      return status;
    }
    source_offset += live_count;
  }
  return lch::detail::FFTResolved(context, accumulator, byte_count, 0,
                                  recovery.size(), kernels, radix);
}

lch::Status MultiplyCopy(Element* destination,
                         const Element* source,
                         size_t byte_count,
                         Element coefficient,
                         const lch::detail::ResolvedKernels& kernels,
                         const lch::Context& context) {
  kernels.scale(destination, source, byte_count, coefficient, context.Tables());
  return lch::Status::ok;
}

lch::Status FormalDerivative(std::span<Element* const> coefficients,
                             size_t byte_count,
                             const lch::detail::ResolvedKernels& kernels) {
  // Novel-basis formal derivative from Leopard-RS, expressed using the
  // repository-owned XOR primitive.
  for (size_t i = 1; i < coefficients.size(); ++i) {
    const size_t width = ((i ^ (i - 1)) + 1) / 2;
    size_t j = i - width;
    for (; j + 4 <= i; j += 4) {
      kernels.xor_four(coefficients[j], coefficients[j + width],
                       coefficients[j + 1], coefficients[j + width + 1],
                       coefficients[j + 2], coefficients[j + width + 2],
                       coefficients[j + 3], coefficients[j + width + 3],
                       byte_count);
    }
    for (; j < i; ++j) {
      kernels.xor_one(coefficients[j], coefficients[j + width], byte_count);
    }
  }
  return lch::Status::ok;
}

namespace {

Element AddMod255(Element a, Element b) {
  unsigned result = static_cast<unsigned>(a) + b;
  if (result >= 255) {
    result -= 255;
  }
  return static_cast<Element>(result);
}

Element SubMod255(Element a, Element b) {
  return static_cast<Element>(a >= b ? a - b : a + 255 - b);
}

void Walsh(std::span<Element, lch::Context::kFieldSize> data,
           size_t live_count) {
  size_t distance = 1;
  for (; 4 * distance <= data.size(); distance *= 4) {
    const size_t group_size = 4 * distance;
    for (size_t block = 0; block < live_count; block += group_size) {
      for (size_t i = 0; i < distance; ++i) {
        Element a = data[block + i];
        Element b = data[block + distance + i];
        Element c = data[block + 2 * distance + i];
        Element d = data[block + 3 * distance + i];
        const Element ab_sum = AddMod255(a, b);
        const Element ab_difference = SubMod255(a, b);
        const Element cd_sum = AddMod255(c, d);
        const Element cd_difference = SubMod255(c, d);
        a = AddMod255(ab_sum, cd_sum);
        b = AddMod255(ab_difference, cd_difference);
        c = SubMod255(ab_sum, cd_sum);
        d = SubMod255(ab_difference, cd_difference);
        data[block + i] = a;
        data[block + distance + i] = b;
        data[block + 2 * distance + i] = c;
        data[block + 3 * distance + i] = d;
      }
    }
  }
  if (distance < data.size()) {
    for (size_t i = 0; i < distance; ++i) {
      const Element a = data[i];
      const Element b = data[distance + i];
      data[i] = AddMod255(a, b);
      data[distance + i] = SubMod255(a, b);
    }
  }
}

}  // namespace

void ErrorLocatorLogs(
    const lch::Context& context,
    const std::array<uint8_t, lch::Context::kFieldSize>& erased,
    size_t live_count,
    std::span<Element, lch::Context::kFieldSize> values) {
  for (size_t i = 0; i < values.size(); ++i) {
    values[i] = erased[i] != 0 ? 1 : 0;
  }
  Walsh(values, live_count);
  for (size_t i = 0; i < values.size(); ++i) {
    values[i] = static_cast<Element>(
        (static_cast<unsigned>(values[i]) * context.LogWalsh(i)) % 255);
  }
  Walsh(values, values.size());
}

void ErrorLocator(const lch::Context& context,
                  const std::array<uint8_t, lch::Context::kFieldSize>& erased,
                  std::span<Element, lch::Context::kFieldSize> values) {
  ErrorLocatorLogs(context, erased, values.size(), values);
  for (Element& value : values) {
    value = context.Exp(value);
  }
}

Element Inverse(const lch::Context& context, Element value) {
  return value == 0 ? 0 : context.Exp(255 - context.Log(value));
}

}  // namespace gf2p8::rs::detail
