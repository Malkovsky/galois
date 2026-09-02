#include "reed_solomon/xdrs.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <vector>

#include "reed_solomon/common.h"

namespace gf2p8::rs {
namespace {

bool IsPowerOfTwo(size_t value) {
  return value != 0 && (value & (value - 1)) == 0;
}

bool HasNull(std::span<const Element* const> pointers) {
  return std::any_of(pointers.begin(), pointers.end(),
                     [](const Element* pointer) { return pointer == nullptr; });
}

bool HasNull(std::span<Element* const> pointers) {
  return std::any_of(pointers.begin(), pointers.end(),
                     [](const Element* pointer) { return pointer == nullptr; });
}

lch::Status FormalDerivative(std::span<Element* const> coefficients,
                             size_t byte_count,
                             const lch::detail::ResolvedKernels& kernels) {
  // The artifact's B pre/post factors are unity in native Cantor coordinates;
  // its clear-and-XOR derivative can therefore run directly in place.
  for (size_t i = 1; i < coefficients.size(); ++i) {
    std::fill_n(coefficients[i - 1], byte_count, 0);
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
  std::fill_n(coefficients.back(), byte_count, 0);
  return lch::Status::ok;
}

}  // namespace

class XDRS::Impl {
 public:
  Impl(size_t data_count, XDRSRate rate)
      : data_count(data_count),
        recovery_count(data_count < 256 ? 256 - data_count : 0),
        rate(rate) {
    valid = data_count != 0 && data_count < 256 &&
            ((rate == XDRSRate::low && data_count <= 128 &&
              IsPowerOfTwo(data_count)) ||
             (rate == XDRSRate::high && recovery_count <= 128 &&
              IsPowerOfTwo(recovery_count)));
    if (!valid) {
      return;
    }

    const lch::Context& context = lch::Context::Shared();

    if (rate == XDRSRate::low) {
      coefficient.resize(256 / data_count);
      Element numerator = 1;
      for (size_t i = 1; i < data_count; ++i) {
        numerator = MultiplyCantor(numerator, context.EvaluationPoint(i));
      }
      for (size_t block = 1; block < 256 / data_count; ++block) {
        Element denominator = 1;
        for (size_t i = 0; i < data_count; ++i) {
          denominator = MultiplyCantor(
              denominator, context.EvaluationPoint(block * data_count) ^
                               context.EvaluationPoint(i));
        }
        coefficient[block] =
            MultiplyCantor(numerator, detail::Inverse(context, denominator));
      }
      return;
    }

    coefficient.resize(data_count);
    Element common = 1;
    for (size_t block = recovery_count; block < 256; block *= 2) {
      for (size_t i = 0; i < block; ++i) {
        common = MultiplyCantor(common, context.EvaluationPoint(block) ^
                                            context.EvaluationPoint(i));
      }
    }
    for (size_t position = recovery_count; position < 256; ++position) {
      Element denominator = common;
      for (size_t i = 0; i < recovery_count; ++i) {
        denominator =
            MultiplyCantor(denominator, context.EvaluationPoint(position) ^
                                            context.EvaluationPoint(i));
      }
      coefficient[position - recovery_count] =
          detail::Inverse(context, denominator);
    }
  }

  size_t data_count;
  size_t recovery_count;
  XDRSRate rate;
  bool valid = false;
  std::vector<Element> coefficient;
};

XDRS::XDRS(size_t data_count, XDRSRate rate)
    : impl_(std::make_unique<Impl>(data_count, rate)) {}

XDRS::~XDRS() = default;
XDRS::XDRS(XDRS&&) noexcept = default;
XDRS& XDRS::operator=(XDRS&&) noexcept = default;

bool XDRS::Valid() const {
  return impl_->valid;
}
size_t XDRS::DataCount() const {
  return impl_->data_count;
}
size_t XDRS::RecoveryCount() const {
  return impl_->recovery_count;
}
XDRSRate XDRS::Rate() const {
  return impl_->rate;
}

size_t XDRS::EncodeWorkspaceSize(size_t byte_count) const {
  if (!Valid() || Rate() == XDRSRate::low) {
    return 0;
  }
  return detail::ShardWorkspaceSize(RecoveryCount(), byte_count);
}

size_t XDRS::DecodeWorkspaceSize(size_t byte_count) const {
  if (!Valid()) {
    return 0;
  }
  return detail::ShardWorkspaceSize(256, byte_count, 256);
}

lch::Status XDRS::Encode(std::span<const Element* const> data,
                         std::span<Element* const> recovery,
                         size_t byte_count,
                         std::span<Element> workspace,
                         lch::Backend backend,
                         lch::Radix radix) const {
  if (!Valid() || data.size() != DataCount() ||
      recovery.size() != RecoveryCount() ||
      workspace.size() < EncodeWorkspaceSize(byte_count) ||
      (radix != lch::Radix::radix2 && radix != lch::Radix::radix4) ||
      (byte_count != 0 && (HasNull(data) || HasNull(recovery)))) {
    return lch::Status::invalid_argument;
  }
  const lch::Status backend_status = detail::ValidateBackend(backend);
  if (backend_status != lch::Status::ok) {
    return backend_status;
  }
  backend = detail::ResolveBackend(backend, byte_count);
  const lch::detail::ResolvedKernels& kernels =
      *lch::detail::ResolveKernels(backend, byte_count);
  const lch::Context& context = lch::Context::Shared();

  // Published XDRS low-rate path: one K-point IFFT followed by independent
  // K-point FFT evaluations for all parity blocks.
  if (Rate() == XDRSRate::low) {
    const size_t k = DataCount();
    lch::Status status = lch::detail::IFFTResolved(
        context, data, std::span<Element* const>(recovery).first(k), byte_count,
        0, kernels, radix);
    if (status != lch::Status::ok) {
      return status;
    }
    for (size_t offset = 2 * k; offset < 256; offset += k) {
      auto block = std::span<Element* const>(recovery).subspan(offset - k, k);
      for (size_t i = 0; i < k; ++i) {
        std::copy_n(recovery[i], byte_count, block[i]);
      }
      status = lch::detail::FFTResolved(context, block, byte_count, offset,
                                        block.size(), kernels, radix);
      if (status != lch::Status::ok) {
        return status;
      }
    }
    return lch::detail::FFTResolved(
        context, std::span<Element* const>(recovery).first(k), byte_count, k, k,
        kernels, radix);
  }

  return detail::EncodeLCH(context, data, recovery, byte_count, RecoveryCount(),
                           workspace, kernels, radix);
}

lch::Status XDRS::Decode(std::span<Element* const> codeword,
                         std::span<const uint8_t> present,
                         size_t byte_count,
                         std::span<Element> workspace,
                         lch::Backend backend,
                         lch::Radix radix) const {
  if (!Valid() || codeword.size() != 256 || present.size() != 256 ||
      workspace.size() < DecodeWorkspaceSize(byte_count) ||
      (radix != lch::Radix::radix2 && radix != lch::Radix::radix4)) {
    return lch::Status::invalid_argument;
  }
  const lch::Status backend_status = detail::ValidateBackend(backend);
  if (backend_status != lch::Status::ok) {
    return backend_status;
  }

  const size_t data_offset = Rate() == XDRSRate::low ? 0 : RecoveryCount();
  Element zero_byte = 0;
  size_t missing_count = 0;
  size_t missing_data_count = 0;
  Element* zero_source = byte_count == 0 ? &zero_byte : nullptr;
  std::array<uint8_t, 256> erased{};
  for (size_t i = 0; i < 256; ++i) {
    erased[i] = present[i] == 0;
    missing_count += erased[i];
    if (byte_count != 0 && present[i] != 0 && codeword[i] == nullptr) {
      return lch::Status::invalid_argument;
    }
  }
  for (size_t i = 0; i < DataCount(); ++i) {
    const size_t position = data_offset + i;
    if (erased[position] != 0) {
      if (byte_count != 0 && codeword[position] == nullptr) {
        return lch::Status::invalid_argument;
      }
      ++missing_data_count;
      if (zero_source == nullptr) {
        zero_source = codeword[position];
      }
    }
  }
  if (missing_count > RecoveryCount()) {
    return lch::Status::insufficient_recovery_symbols;
  }
  if (missing_data_count == 0) {
    return lch::Status::ok;
  }
  std::fill_n(zero_source, byte_count, 0);

  backend = detail::ResolveBackend(backend, byte_count);
  const lch::detail::ResolvedKernels& kernels =
      *lch::detail::ResolveKernels(backend, byte_count);
  const lch::Context& context = lch::Context::Shared();
  auto locator = std::span<Element, 256>(workspace.data(), 256);
  detail::ErrorLocator(context, erased, locator);

  std::array<Element*, 256> scratch{};
  Element* const scratch_base =
      byte_count == 0 ? &zero_byte : workspace.data() + 256;
  for (size_t i = 0; i < scratch.size(); ++i) {
    scratch[i] = scratch_base + i * byte_count;
  }

  if (Rate() == XDRSRate::low) {
    for (size_t i = 0; i < 256; ++i) {
      if (erased[i] != 0) {
        std::fill_n(scratch[i], byte_count, 0);
      } else {
        const lch::Status status = detail::MultiplyCopy(
            scratch[i], codeword[i], byte_count, locator[i], kernels, context);
        if (status != lch::Status::ok) {
          return status;
        }
      }
    }

    const size_t k = DataCount();
    for (size_t offset = 0; offset < 256; offset += k) {
      const lch::Status status = lch::detail::IFFTResolved(
          context, std::span<Element* const>(scratch).subspan(offset, k),
          byte_count, offset, k, kernels, radix);
      if (status != lch::Status::ok) {
        return status;
      }
    }
    lch::Status status = FormalDerivative(
        std::span<Element* const>(scratch).first(k), byte_count, kernels);
    if (status != lch::Status::ok) {
      return status;
    }
    for (size_t block = 1; block < 256 / k; ++block) {
      for (size_t i = 0; i < k; ++i) {
        kernels.add_scaled(scratch[i], scratch[block * k + i], byte_count,
                           impl_->coefficient[block], context.Tables());
      }
    }
    status = lch::detail::FFTResolved(
        context, std::span<Element* const>(scratch).first(k), byte_count, 0, k,
        kernels, radix);
    if (status != lch::Status::ok) {
      return status;
    }
    for (size_t i = 0; i < k; ++i) {
      if (erased[i] != 0) {
        status = detail::MultiplyCopy(codeword[i], scratch[i], byte_count,
                                      detail::Inverse(context, locator[i]),
                                      kernels, context);
        if (status != lch::Status::ok) {
          return status;
        }
      }
    }
    return lch::Status::ok;
  }

  // Published high-rate path computes a T-point companion polynomial. Fold
  // each input block directly into the accumulator instead of copying the
  // full codeword and making a separate XOR pass.
  const size_t t = RecoveryCount();
  std::array<const Element*, 128> input{};
  auto accumulator = std::span<Element* const>(scratch).first(t);
  lch::Status status = lch::Status::ok;
  for (size_t offset = 0; offset < 256; offset += t) {
    for (size_t i = 0; i < t; ++i) {
      const size_t position = offset + i;
      input[i] = erased[position] != 0 ? zero_source : codeword[position];
    }
    const auto source = std::span<const Element* const>(input).first(t);
    if (offset == 0) {
      status = lch::detail::IFFTResolved(context, source, accumulator,
                                         byte_count, 0, kernels, radix);
    } else {
      auto block = std::span<Element* const>(scratch).subspan(offset, t);
      status = lch::detail::IFFTResolved(context, source, block, accumulator,
                                         byte_count, offset, kernels, radix);
    }
    if (status != lch::Status::ok) {
      return status;
    }
  }
  status = lch::detail::FFTResolved(context, accumulator, byte_count, 0, t,
                                    kernels, radix);
  if (status != lch::Status::ok) {
    return status;
  }

  for (size_t i = 0; i < t; ++i) {
    if (erased[i] != 0) {
      std::fill_n(scratch[t + i], byte_count, 0);
    } else {
      status = detail::MultiplyCopy(scratch[t + i], scratch[i], byte_count,
                                    locator[i], kernels, context);
      if (status != lch::Status::ok) {
        return status;
      }
    }
  }
  auto coefficients = std::span<Element* const>(scratch).subspan(t, t);
  status = lch::detail::IFFTResolved(context, coefficients, byte_count, 0, t,
                                     kernels, radix);
  if (status != lch::Status::ok) {
    return status;
  }

  std::array<uint8_t, 128> requested{};
  for (size_t offset = t; offset < 256; offset += t) {
    bool has_requested_output = false;
    for (size_t i = 0; i < t; ++i) {
      requested[i] = erased[offset + i];
      has_requested_output |= requested[i] != 0;
    }
    if (!has_requested_output) {
      continue;
    }
    auto output = offset == t
                      ? std::span<Element* const>(scratch).first(t)
                      : std::span<Element* const>(scratch).subspan(offset, t);
    for (size_t i = 0; i < t; ++i) {
      std::copy_n(coefficients[i], byte_count, output[i]);
    }
    status = lch::detail::FFTResolved(
        context, output, byte_count, offset,
        std::span<const uint8_t>(requested).first(t), kernels, radix);
    if (status != lch::Status::ok) {
      return status;
    }
    for (size_t i = 0; i < t; ++i) {
      const size_t position = offset + i;
      if (erased[position] == 0) {
        continue;
      }
      const Element factor =
          MultiplyCantor(impl_->coefficient[position - t],
                         detail::Inverse(context, locator[position]));
      status = detail::MultiplyCopy(codeword[position], output[i], byte_count,
                                    factor, kernels, context);
      if (status != lch::Status::ok) {
        return status;
      }
    }
  }
  return lch::Status::ok;
}

}  // namespace gf2p8::rs
