#include "reed_solomon/lch_encoder.h"

#include <algorithm>
#include <array>

#include "lin_chung_han/transform_internal.h"
#include "reed_solomon/code_parameters.h"

namespace gf2p8::rs {
namespace {

bool HasNull(std::span<const Element* const> pointers) {
  return std::any_of(pointers.begin(), pointers.end(),
                     [](const Element* pointer) { return pointer == nullptr; });
}

bool HasNull(std::span<Element* const> pointers) {
  return std::any_of(pointers.begin(), pointers.end(),
                     [](const Element* pointer) { return pointer == nullptr; });
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

lch::Status EncodeFolded(const detail::CodeParameters& parameters,
                         std::span<const Element* const> data,
                         std::span<Element* const> recovery,
                         size_t byte_count,
                         std::span<Element> workspace,
                         const lch::detail::ResolvedKernels& kernels,
                         lch::Radix radix) {
  const size_t transform_size = parameters.transform_size;
  if (transform_size == 1) {
    std::copy_n(data[0], byte_count, recovery[0]);
    for (size_t i = 1; i < data.size(); ++i) {
      kernels.xor_one(recovery[0], data[i], byte_count);
    }
    return lch::Status::ok;
  }

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

  const lch::Context& context = lch::Context::Shared();
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

lch::Status EncodeCoefficientFanout(const detail::CodeParameters& parameters,
                                    std::span<const Element* const> data,
                                    std::span<Element* const> recovery,
                                    size_t byte_count,
                                    std::span<Element> workspace,
                                    const lch::detail::ResolvedKernels& kernels,
                                    lch::Radix radix) {
  const size_t transform_size = parameters.transform_size;
  std::array<Element*, lch::Context::kFieldSize> block_pointers;
  Element* const padding = workspace.data();
  const auto make_block = [&](size_t recovery_offset, size_t live_count) {
    if (live_count == transform_size) {
      return std::span<Element* const>(recovery).subspan(recovery_offset,
                                                         transform_size);
    }
    for (size_t i = 0; i < live_count; ++i) {
      block_pointers[i] = recovery[recovery_offset + i];
    }
    for (size_t i = live_count; i < transform_size; ++i) {
      block_pointers[i] = padding + (i - live_count) * byte_count;
    }
    return std::span<Element* const>(block_pointers).first(transform_size);
  };

  const size_t first_count = std::min(transform_size, recovery.size());
  auto coefficients = make_block(0, first_count);
  const lch::Context& context = lch::Context::Shared();
  lch::Status status = lch::detail::IFFTResolved(context, data, coefficients,
                                                 byte_count, 0, kernels, radix);
  if (status != lch::Status::ok) {
    return status;
  }

  for (size_t output_offset = transform_size; output_offset < recovery.size();
       output_offset += transform_size) {
    const size_t live_count =
        std::min(transform_size, recovery.size() - output_offset);
    auto output = make_block(output_offset, live_count);
    for (size_t i = 0; i < transform_size; ++i) {
      std::copy_n(coefficients[i], byte_count, output[i]);
    }
    status = lch::detail::FFTResolved(context, output, byte_count,
                                      transform_size + output_offset,
                                      live_count, kernels, radix);
    if (status != lch::Status::ok) {
      return status;
    }
  }
  return lch::detail::FFTResolved(context, coefficients, byte_count,
                                  transform_size, first_count, kernels, radix);
}

}  // namespace

class LCHEncoder::Impl {
 public:
  explicit Impl(size_t data_count, size_t recovery_count)
      : parameters(detail::MakeCodeParameters(data_count, recovery_count)) {}

  detail::CodeParameters parameters;
};

LCHEncoder::LCHEncoder(size_t data_count, size_t recovery_count)
    : impl_(std::make_unique<Impl>(data_count, recovery_count)) {}

LCHEncoder::~LCHEncoder() = default;
LCHEncoder::LCHEncoder(LCHEncoder&&) noexcept = default;
LCHEncoder& LCHEncoder::operator=(LCHEncoder&&) noexcept = default;

bool LCHEncoder::Valid() const {
  return impl_ != nullptr && impl_->parameters.valid;
}

size_t LCHEncoder::DataCount() const {
  return impl_ == nullptr ? 0 : impl_->parameters.data_count;
}

size_t LCHEncoder::RecoveryCount() const {
  return impl_ == nullptr ? 0 : impl_->parameters.recovery_count;
}

size_t LCHEncoder::WorkspaceSize(size_t byte_count) const {
  if (!Valid()) {
    return 0;
  }
  const size_t transform_size = impl_->parameters.transform_size;
  if (impl_->parameters.family == detail::CodeFamily::high_rate) {
    if (transform_size == 1) {
      return 0;
    }
    return detail::CheckedWorkspaceSize(2 * transform_size - RecoveryCount(),
                                        byte_count);
  }
  const size_t remainder = RecoveryCount() % transform_size;
  const size_t padding_count = remainder == 0 ? 0 : transform_size - remainder;
  return detail::CheckedWorkspaceSize(padding_count, byte_count);
}

lch::Status LCHEncoder::Encode(std::span<const Element* const> data,
                               std::span<Element* const> recovery,
                               size_t byte_count,
                               std::span<Element> workspace,
                               lch::Backend backend,
                               lch::Radix radix) const {
  if (!Valid() || data.size() != DataCount() ||
      recovery.size() != RecoveryCount() ||
      workspace.size() < WorkspaceSize(byte_count) ||
      (radix != lch::Radix::radix2 && radix != lch::Radix::radix4) ||
      (byte_count != 0 && (HasNull(data) || HasNull(recovery)))) {
    return lch::Status::invalid_argument;
  }
  const lch::Status backend_status = ValidateBackend(backend);
  if (backend_status != lch::Status::ok) {
    return backend_status;
  }
  backend = ResolveBackend(backend, byte_count);
  const lch::detail::ResolvedKernels& kernels =
      *lch::detail::ResolveKernels(backend, byte_count);
  if (byte_count == 0) {
    return lch::Status::ok;
  }

  if (impl_->parameters.family == detail::CodeFamily::low_rate) {
    return EncodeCoefficientFanout(impl_->parameters, data, recovery,
                                   byte_count, workspace, kernels, radix);
  }
  return EncodeFolded(impl_->parameters, data, recovery, byte_count, workspace,
                      kernels, radix);
}

}  // namespace gf2p8::rs
