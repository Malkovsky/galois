#include "reed_solomon/lch_decoder.h"

#include <algorithm>
#include <vector>

#include "lin_chung_han/transform_internal.h"
#include "reed_solomon/code_parameters.h"
#include "reed_solomon/decoder/internal.h"

namespace gf2p8::rs {
namespace {

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

}  // namespace

class LCHDecoder::Impl {
 public:
  explicit Impl(size_t data_count, size_t recovery_count)
      : parameters(detail::MakeCodeParameters(data_count, recovery_count)) {
    if (!parameters.valid) {
      return;
    }
    if (parameters.family == detail::CodeFamily::low_rate) {
      coefficients = detail::decoder::MakeLowRateCoefficients(parameters);
    } else {
      coefficients = detail::decoder::MakeHighRateCoefficients(parameters);
    }
  }

  detail::CodeParameters parameters;
  std::vector<Element> coefficients;
};

LCHDecoder::LCHDecoder(size_t data_count, size_t recovery_count)
    : impl_(std::make_unique<Impl>(data_count, recovery_count)) {}

LCHDecoder::~LCHDecoder() = default;
LCHDecoder::LCHDecoder(LCHDecoder&&) noexcept = default;
LCHDecoder& LCHDecoder::operator=(LCHDecoder&&) noexcept = default;

bool LCHDecoder::Valid() const {
  return impl_ != nullptr && impl_->parameters.valid;
}

size_t LCHDecoder::DataCount() const {
  return impl_ == nullptr ? 0 : impl_->parameters.data_count;
}

size_t LCHDecoder::RecoveryCount() const {
  return impl_ == nullptr ? 0 : impl_->parameters.recovery_count;
}

size_t LCHDecoder::WorkspaceSize(size_t byte_count) const {
  if (!Valid()) {
    return 0;
  }
  return detail::CheckedWorkspaceSize(impl_->parameters.mother_size, byte_count,
                                      lch::Context::kFieldSize);
}

lch::Status LCHDecoder::Decode(std::span<Element* const> data,
                               std::span<const uint8_t> data_present,
                               std::span<const Element* const> recovery,
                               std::span<const uint8_t> recovery_present,
                               size_t byte_count,
                               std::span<Element> workspace,
                               lch::Backend backend,
                               lch::Radix radix) const {
  if (!Valid() || data.size() != DataCount() ||
      data_present.size() != DataCount() ||
      recovery.size() != RecoveryCount() ||
      recovery_present.size() != RecoveryCount() ||
      workspace.size() < WorkspaceSize(byte_count) ||
      (radix != lch::Radix::radix2 && radix != lch::Radix::radix4)) {
    return lch::Status::invalid_argument;
  }
  if (byte_count != 0) {
    if (std::any_of(data.begin(), data.end(), [](const Element* pointer) {
          return pointer == nullptr;
        })) {
      return lch::Status::invalid_argument;
    }
    for (size_t i = 0; i < recovery.size(); ++i) {
      if (recovery_present[i] != 0 && recovery[i] == nullptr) {
        return lch::Status::invalid_argument;
      }
    }
  }
  const lch::Status backend_status = ValidateBackend(backend);
  if (backend_status != lch::Status::ok) {
    return backend_status;
  }

  size_t missing_count = 0;
  size_t missing_data_count = 0;
  for (const uint8_t present : data_present) {
    missing_data_count += present == 0;
  }
  missing_count = missing_data_count;
  for (const uint8_t present : recovery_present) {
    missing_count += present == 0;
  }
  if (missing_count > RecoveryCount()) {
    return lch::Status::insufficient_recovery_symbols;
  }
  if (missing_data_count == 0) {
    return lch::Status::ok;
  }

  backend = ResolveBackend(backend, byte_count);
  const lch::detail::ResolvedKernels& kernels =
      *lch::detail::ResolveKernels(backend, byte_count);
  const detail::decoder::Inputs inputs{
      data, data_present, recovery, recovery_present, byte_count, workspace};
  if (impl_->parameters.family == detail::CodeFamily::low_rate) {
    return detail::decoder::DecodeLowRate(
        impl_->parameters, impl_->coefficients, inputs, kernels, radix);
  }
  return detail::decoder::DecodeHighRate(impl_->parameters, impl_->coefficients,
                                         inputs, kernels, radix);
}

}  // namespace gf2p8::rs
