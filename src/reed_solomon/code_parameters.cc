#include "reed_solomon/code_parameters.h"

#include <limits>

#include "lin_chung_han/transform.h"

namespace gf2p8::rs::detail {
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

}  // namespace

CodeParameters MakeCodeParameters(size_t data_count, size_t recovery_count) {
  CodeParameters result;
  result.data_count = data_count;
  result.recovery_count = recovery_count;
  if (data_count == 0 || recovery_count == 0) {
    return result;
  }

  if (recovery_count >= data_count) {
    result.transform_size = NextPowerOfTwo(data_count);
    if (result.transform_size == 0 ||
        recovery_count > lch::Context::kFieldSize - result.transform_size) {
      return result;
    }
    result.mother_size = NextPowerOfTwo(result.transform_size + recovery_count);
    result.family = CodeFamily::low_rate;
    result.valid = result.mother_size != 0;
    return result;
  }

  result.transform_size = NextPowerOfTwo(recovery_count);
  if (result.transform_size == 0 ||
      data_count > lch::Context::kFieldSize - result.transform_size) {
    return result;
  }
  result.mother_size = NextPowerOfTwo(result.transform_size + data_count);
  result.family = CodeFamily::high_rate;
  result.valid = result.mother_size != 0;
  return result;
}

size_t CheckedWorkspaceSize(size_t shard_count,
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

}  // namespace gf2p8::rs::detail
