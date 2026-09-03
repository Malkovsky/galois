#pragma once

#include <cstddef>

namespace gf2p8::rs::detail {

enum class CodeFamily { low_rate, high_rate };

struct CodeParameters {
  size_t data_count = 0;
  size_t recovery_count = 0;
  size_t transform_size = 0;
  size_t mother_size = 0;
  CodeFamily family = CodeFamily::high_rate;
  bool valid = false;
};

CodeParameters MakeCodeParameters(size_t data_count, size_t recovery_count);
size_t CheckedWorkspaceSize(size_t shard_count,
                            size_t byte_count,
                            size_t prefix_bytes = 0);

}  // namespace gf2p8::rs::detail
