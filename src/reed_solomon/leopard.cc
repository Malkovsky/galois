#include "reed_solomon/leopard.h"

#include <algorithm>
#include <array>

#include "reed_solomon/common.h"

namespace gf2p8::rs {
namespace {

size_t NextPowerOfTwo(size_t value) {
  if (value == 0 || value > lch::Context::kFieldSize) {
    return 0;
  }
  size_t result = 1;
  while (result < value) {
    result *= 2;
  }
  return result;
}

size_t WorkCount(size_t transform_size, size_t original_count) {
  if (transform_size == 0 ||
      original_count > lch::Context::kFieldSize - transform_size) {
    return 0;
  }
  return NextPowerOfTwo(transform_size + original_count);
}

bool HasNull(std::span<const Element* const> pointers) {
  return std::any_of(pointers.begin(), pointers.end(),
                     [](const Element* pointer) { return pointer == nullptr; });
}

bool HasNull(std::span<Element* const> pointers) {
  return std::any_of(pointers.begin(), pointers.end(),
                     [](const Element* pointer) { return pointer == nullptr; });
}

}  // namespace

class Leopard::Impl {
 public:
  Impl(size_t original_count, size_t recovery_count)
      : original_count(original_count),
        recovery_count(recovery_count),
        transform_size(NextPowerOfTwo(recovery_count)),
        work_count(WorkCount(transform_size, original_count)) {
    valid = original_count != 0 && recovery_count != 0 &&
            recovery_count <= original_count && transform_size != 0 &&
            work_count != 0;
  }

  size_t original_count;
  size_t recovery_count;
  size_t transform_size;
  size_t work_count;
  bool valid = false;
};

Leopard::Leopard(size_t original_count, size_t recovery_count)
    : impl_(std::make_unique<Impl>(original_count, recovery_count)) {}

Leopard::~Leopard() = default;
Leopard::Leopard(Leopard&&) noexcept = default;
Leopard& Leopard::operator=(Leopard&&) noexcept = default;

bool Leopard::Valid() const {
  return impl_->valid;
}
size_t Leopard::OriginalCount() const {
  return impl_->original_count;
}
size_t Leopard::RecoveryCount() const {
  return impl_->recovery_count;
}

size_t Leopard::EncodeWorkspaceSize(size_t byte_count) const {
  if (!Valid()) {
    return 0;
  }
  if (impl_->transform_size == 1) {
    return 0;
  }
  return detail::ShardWorkspaceSize(
      2 * impl_->transform_size - impl_->recovery_count, byte_count);
}

size_t Leopard::DecodeWorkspaceSize(size_t byte_count) const {
  if (!Valid()) {
    return 0;
  }
  return detail::ShardWorkspaceSize(impl_->work_count, byte_count, 256);
}

lch::Status Leopard::Encode(std::span<const Element* const> original,
                            std::span<Element* const> recovery,
                            size_t byte_count,
                            std::span<Element> workspace,
                            lch::Backend backend,
                            lch::Radix radix) const {
  if (!Valid() || original.size() != OriginalCount() ||
      recovery.size() != RecoveryCount() ||
      workspace.size() < EncodeWorkspaceSize(byte_count) ||
      (radix != lch::Radix::radix2 && radix != lch::Radix::radix4) ||
      (byte_count != 0 && (HasNull(original) || HasNull(recovery)))) {
    return lch::Status::invalid_argument;
  }
  const lch::Status backend_status = detail::ValidateBackend(backend);
  if (backend_status != lch::Status::ok) {
    return backend_status;
  }
  backend = detail::ResolveBackend(backend, byte_count);
  const lch::detail::ResolvedKernels& kernels =
      *lch::detail::ResolveKernels(backend, byte_count);

  return detail::EncodeLCH(lch::Context::Shared(), original, recovery,
                           byte_count, impl_->transform_size, workspace,
                           kernels, radix);
}

lch::Status Leopard::Decode(std::span<Element* const> codeword,
                            std::span<const uint8_t> present,
                            size_t byte_count,
                            std::span<Element> workspace,
                            lch::Backend backend,
                            lch::Radix radix) const {
  const size_t logical_count = RecoveryCount() + OriginalCount();
  if (!Valid() || codeword.size() != logical_count ||
      present.size() != logical_count ||
      workspace.size() < DecodeWorkspaceSize(byte_count) ||
      (radix != lch::Radix::radix2 && radix != lch::Radix::radix4)) {
    return lch::Status::invalid_argument;
  }
  const lch::Status backend_status = detail::ValidateBackend(backend);
  if (backend_status != lch::Status::ok) {
    return backend_status;
  }

  size_t missing_count = 0;
  size_t missing_original_count = 0;
  for (size_t i = 0; i < logical_count; ++i) {
    missing_count += present[i] == 0;
    if (byte_count != 0 && present[i] != 0 && codeword[i] == nullptr) {
      return lch::Status::invalid_argument;
    }
  }
  for (size_t i = 0; i < OriginalCount(); ++i) {
    const size_t logical_index = RecoveryCount() + i;
    missing_original_count += present[logical_index] == 0;
    if (byte_count != 0 && present[logical_index] == 0 &&
        codeword[logical_index] == nullptr) {
      return lch::Status::invalid_argument;
    }
  }
  if (missing_count > RecoveryCount()) {
    return lch::Status::insufficient_recovery_symbols;
  }
  if (missing_original_count == 0) {
    return lch::Status::ok;
  }

  backend = detail::ResolveBackend(backend, byte_count);
  const lch::detail::ResolvedKernels& kernels =
      *lch::detail::ResolveKernels(backend, byte_count);
  const size_t m = impl_->transform_size;
  const size_t n = impl_->work_count;
  const lch::Context& context = lch::Context::Shared();
  std::array<uint8_t, 256> erased{};
  for (size_t i = 0; i < RecoveryCount(); ++i) {
    erased[i] = present[i] == 0;
  }
  for (size_t i = RecoveryCount(); i < m; ++i) {
    erased[i] = 1;
  }
  for (size_t i = 0; i < OriginalCount(); ++i) {
    erased[m + i] = present[RecoveryCount() + i] == 0;
  }

  auto locator = std::span<Element, 256>(workspace.data(), 256);
  detail::ErrorLocatorLogs(context, erased, m + OriginalCount(), locator);
  std::array<Element*, 256> work{};
  Element zero_byte = 0;
  Element* const base = byte_count == 0 ? &zero_byte : workspace.data() + 256;
  for (size_t i = 0; i < n; ++i) {
    work[i] = base + i * byte_count;
  }

  // Each work shard is written exactly once: present symbols are scale-copied
  // directly, while erased, padded, and trailing positions are cleared.
  for (size_t position = 0; position < n; ++position) {
    const Element* source = nullptr;
    if (position < RecoveryCount() && present[position] != 0) {
      source = codeword[position];
    } else if (position >= m && position < m + OriginalCount()) {
      const size_t logical_index = RecoveryCount() + position - m;
      if (present[logical_index] != 0) {
        source = codeword[logical_index];
      }
    }
    if (source == nullptr) {
      std::fill_n(work[position], byte_count, 0);
      continue;
    }
    kernels.scale(work[position], source, byte_count,
                  context.Exp(locator[position]), context.Tables());
  }

  auto transform = std::span<Element* const>(work).first(n);
  lch::Status status = lch::detail::IFFTResolved(
      context, transform, byte_count, 0, m + OriginalCount(), kernels, radix);
  if (status != lch::Status::ok) {
    return status;
  }
  status = detail::FormalDerivative(transform, byte_count, kernels);
  if (status != lch::Status::ok) {
    return status;
  }
  std::array<uint8_t, 256> requested{};
  for (size_t i = 0; i < OriginalCount(); ++i) {
    if (present[RecoveryCount() + i] == 0) {
      requested[m + i] = 1;
    }
  }
  status = lch::detail::FFTResolved(
      context, transform, byte_count, 0,
      std::span<const uint8_t>(requested).first(n), kernels, radix);
  if (status != lch::Status::ok) {
    return status;
  }

  for (size_t i = 0; i < OriginalCount(); ++i) {
    const size_t logical_index = RecoveryCount() + i;
    if (present[logical_index] == 0) {
      status = detail::MultiplyCopy(
          codeword[logical_index], work[m + i], byte_count,
          context.Exp(255 - locator[m + i]), kernels, context);
      if (status != lch::Status::ok) {
        return status;
      }
    }
  }
  return lch::Status::ok;
}

}  // namespace gf2p8::rs
