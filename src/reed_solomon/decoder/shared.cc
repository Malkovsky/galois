#include <algorithm>

#include "reed_solomon/decoder/internal.h"

namespace gf2p8::rs::detail::decoder {

lch::Status MultiplyCopy(Element* destination,
                         const Element* source,
                         size_t byte_count,
                         Element coefficient,
                         const lch::detail::ResolvedKernels& kernels,
                         const lch::Context& context) {
  kernels.scale(destination, source, byte_count, coefficient, context.Tables());
  return lch::Status::ok;
}

lch::Status LowRateFormalDerivative(
    std::span<Element* const> coefficients,
    size_t byte_count,
    const lch::detail::ResolvedKernels& kernels) {
  // XDRS Algorithm 2 uses this clear-and-XOR derivative. Its published B
  // pre/post factors are identity in the shared Cantor-coordinate domain.
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

void ErrorLocatorLogs(const lch::Context& context,
                      const Erased& erased,
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

}  // namespace

void ErrorLocator(const lch::Context& context,
                  const Erased& erased,
                  size_t live_count,
                  std::span<Element, lch::Context::kFieldSize> values) {
  ErrorLocatorLogs(context, erased, live_count, values);
  for (Element& value : values) {
    value = context.Exp(value);
  }
}

Element Inverse(const lch::Context& context, Element value) {
  return value == 0 ? 0 : context.Exp(255 - context.Log(value));
}

}  // namespace gf2p8::rs::detail::decoder
