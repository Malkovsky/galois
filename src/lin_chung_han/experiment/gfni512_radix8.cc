#include "lin_chung_han/experiment/gfni512_radix8.h"

#include <algorithm>

#if defined(__i386__) || defined(__x86_64__) || defined(_M_IX86) || \
    defined(_M_X64)
#include <immintrin.h>
#endif

namespace gf2p8::lch::detail::experiment::radix8 {
namespace {

size_t Log2(size_t value) {
  size_t result = 0;
  while (value > 1) {
    value >>= 1;
    ++result;
  }
  return result;
}

#if defined(__GFNI__) && defined(__AVX512F__) && defined(__AVX512BW__)

Element Product(Element value,
                Element coefficient,
                const MultiplicationTables& tables) {
  const auto& row = tables.shuffle[coefficient];
  return row[value & 0x0f] ^ row[32 + (value >> 4)];
}

void ScalarIFFT(Element& x,
                Element& y,
                Element coefficient,
                const MultiplicationTables& tables) {
  y ^= x;
  if (coefficient != 0) {
    x ^= Product(y, coefficient, tables);
  }
}

void IFFT(Element* x0,
          Element* x1,
          Element* x2,
          Element* x3,
          Element* x4,
          Element* x5,
          Element* x6,
          Element* x7,
          size_t byte_count,
          Element top,
          Element middle0,
          Element middle1,
          Element low0,
          Element low1,
          Element low2,
          Element low3,
          const MultiplicationTables& tables) {
  const __m512i top_matrix =
      _mm512_set1_epi64(static_cast<long long>(tables.affine[top]));
  const __m512i middle0_matrix =
      _mm512_set1_epi64(static_cast<long long>(tables.affine[middle0]));
  const __m512i middle1_matrix =
      _mm512_set1_epi64(static_cast<long long>(tables.affine[middle1]));
  const __m512i low0_matrix =
      _mm512_set1_epi64(static_cast<long long>(tables.affine[low0]));
  const __m512i low1_matrix =
      _mm512_set1_epi64(static_cast<long long>(tables.affine[low1]));
  const __m512i low2_matrix =
      _mm512_set1_epi64(static_cast<long long>(tables.affine[low2]));
  const __m512i low3_matrix =
      _mm512_set1_epi64(static_cast<long long>(tables.affine[low3]));

  size_t processed = 0;
  for (; processed + 64 <= byte_count; processed += 64) {
    __m512i a = _mm512_loadu_si512(x0 + processed);
    __m512i b = _mm512_loadu_si512(x1 + processed);
    __m512i c = _mm512_loadu_si512(x2 + processed);
    __m512i d = _mm512_loadu_si512(x3 + processed);
    __m512i e = _mm512_loadu_si512(x4 + processed);
    __m512i f = _mm512_loadu_si512(x5 + processed);
    __m512i g = _mm512_loadu_si512(x6 + processed);
    __m512i h = _mm512_loadu_si512(x7 + processed);

    b = _mm512_xor_si512(b, a);
    if (low0 != 0) {
      a = _mm512_xor_si512(a, _mm512_gf2p8affine_epi64_epi8(b, low0_matrix, 0));
    }
    d = _mm512_xor_si512(d, c);
    if (low1 != 0) {
      c = _mm512_xor_si512(c, _mm512_gf2p8affine_epi64_epi8(d, low1_matrix, 0));
    }
    f = _mm512_xor_si512(f, e);
    if (low2 != 0) {
      e = _mm512_xor_si512(e, _mm512_gf2p8affine_epi64_epi8(f, low2_matrix, 0));
    }
    h = _mm512_xor_si512(h, g);
    if (low3 != 0) {
      g = _mm512_xor_si512(g, _mm512_gf2p8affine_epi64_epi8(h, low3_matrix, 0));
    }

    c = _mm512_xor_si512(c, a);
    if (middle0 != 0) {
      a = _mm512_xor_si512(a,
                           _mm512_gf2p8affine_epi64_epi8(c, middle0_matrix, 0));
    }
    d = _mm512_xor_si512(d, b);
    if (middle0 != 0) {
      b = _mm512_xor_si512(b,
                           _mm512_gf2p8affine_epi64_epi8(d, middle0_matrix, 0));
    }
    g = _mm512_xor_si512(g, e);
    if (middle1 != 0) {
      e = _mm512_xor_si512(e,
                           _mm512_gf2p8affine_epi64_epi8(g, middle1_matrix, 0));
    }
    h = _mm512_xor_si512(h, f);
    if (middle1 != 0) {
      f = _mm512_xor_si512(f,
                           _mm512_gf2p8affine_epi64_epi8(h, middle1_matrix, 0));
    }

    e = _mm512_xor_si512(e, a);
    if (top != 0) {
      a = _mm512_xor_si512(a, _mm512_gf2p8affine_epi64_epi8(e, top_matrix, 0));
    }
    f = _mm512_xor_si512(f, b);
    if (top != 0) {
      b = _mm512_xor_si512(b, _mm512_gf2p8affine_epi64_epi8(f, top_matrix, 0));
    }
    g = _mm512_xor_si512(g, c);
    if (top != 0) {
      c = _mm512_xor_si512(c, _mm512_gf2p8affine_epi64_epi8(g, top_matrix, 0));
    }
    h = _mm512_xor_si512(h, d);
    if (top != 0) {
      d = _mm512_xor_si512(d, _mm512_gf2p8affine_epi64_epi8(h, top_matrix, 0));
    }

    _mm512_storeu_si512(x0 + processed, a);
    _mm512_storeu_si512(x1 + processed, b);
    _mm512_storeu_si512(x2 + processed, c);
    _mm512_storeu_si512(x3 + processed, d);
    _mm512_storeu_si512(x4 + processed, e);
    _mm512_storeu_si512(x5 + processed, f);
    _mm512_storeu_si512(x6 + processed, g);
    _mm512_storeu_si512(x7 + processed, h);
  }

  for (size_t i = processed; i < byte_count; ++i) {
    ScalarIFFT(x0[i], x1[i], low0, tables);
    ScalarIFFT(x2[i], x3[i], low1, tables);
    ScalarIFFT(x4[i], x5[i], low2, tables);
    ScalarIFFT(x6[i], x7[i], low3, tables);
    ScalarIFFT(x0[i], x2[i], middle0, tables);
    ScalarIFFT(x1[i], x3[i], middle0, tables);
    ScalarIFFT(x4[i], x6[i], middle1, tables);
    ScalarIFFT(x5[i], x7[i], middle1, tables);
    ScalarIFFT(x0[i], x4[i], top, tables);
    ScalarIFFT(x1[i], x5[i], top, tables);
    ScalarIFFT(x2[i], x6[i], top, tables);
    ScalarIFFT(x3[i], x7[i], top, tables);
  }
}

void IFFTXor(const Element* x0,
             const Element* x1,
             const Element* x2,
             const Element* x3,
             const Element* x4,
             const Element* x5,
             const Element* x6,
             const Element* x7,
             Element* output0,
             Element* output1,
             Element* output2,
             Element* output3,
             Element* output4,
             Element* output5,
             Element* output6,
             Element* output7,
             size_t byte_count,
             Element top,
             Element middle0,
             Element middle1,
             Element low0,
             Element low1,
             Element low2,
             Element low3,
             const MultiplicationTables& tables) {
  const __m512i top_matrix =
      _mm512_set1_epi64(static_cast<long long>(tables.affine[top]));
  const __m512i middle0_matrix =
      _mm512_set1_epi64(static_cast<long long>(tables.affine[middle0]));
  const __m512i middle1_matrix =
      _mm512_set1_epi64(static_cast<long long>(tables.affine[middle1]));
  const __m512i low0_matrix =
      _mm512_set1_epi64(static_cast<long long>(tables.affine[low0]));
  const __m512i low1_matrix =
      _mm512_set1_epi64(static_cast<long long>(tables.affine[low1]));
  const __m512i low2_matrix =
      _mm512_set1_epi64(static_cast<long long>(tables.affine[low2]));
  const __m512i low3_matrix =
      _mm512_set1_epi64(static_cast<long long>(tables.affine[low3]));

  size_t processed = 0;
  for (; processed + 64 <= byte_count; processed += 64) {
    __m512i a = _mm512_loadu_si512(x0 + processed);
    __m512i b = _mm512_loadu_si512(x1 + processed);
    __m512i c = _mm512_loadu_si512(x2 + processed);
    __m512i d = _mm512_loadu_si512(x3 + processed);
    __m512i e = _mm512_loadu_si512(x4 + processed);
    __m512i f = _mm512_loadu_si512(x5 + processed);
    __m512i g = _mm512_loadu_si512(x6 + processed);
    __m512i h = _mm512_loadu_si512(x7 + processed);

    b = _mm512_xor_si512(b, a);
    if (low0 != 0) {
      a = _mm512_xor_si512(a, _mm512_gf2p8affine_epi64_epi8(b, low0_matrix, 0));
    }
    d = _mm512_xor_si512(d, c);
    if (low1 != 0) {
      c = _mm512_xor_si512(c, _mm512_gf2p8affine_epi64_epi8(d, low1_matrix, 0));
    }
    f = _mm512_xor_si512(f, e);
    if (low2 != 0) {
      e = _mm512_xor_si512(e, _mm512_gf2p8affine_epi64_epi8(f, low2_matrix, 0));
    }
    h = _mm512_xor_si512(h, g);
    if (low3 != 0) {
      g = _mm512_xor_si512(g, _mm512_gf2p8affine_epi64_epi8(h, low3_matrix, 0));
    }

    c = _mm512_xor_si512(c, a);
    if (middle0 != 0) {
      a = _mm512_xor_si512(a,
                           _mm512_gf2p8affine_epi64_epi8(c, middle0_matrix, 0));
    }
    d = _mm512_xor_si512(d, b);
    if (middle0 != 0) {
      b = _mm512_xor_si512(b,
                           _mm512_gf2p8affine_epi64_epi8(d, middle0_matrix, 0));
    }
    g = _mm512_xor_si512(g, e);
    if (middle1 != 0) {
      e = _mm512_xor_si512(e,
                           _mm512_gf2p8affine_epi64_epi8(g, middle1_matrix, 0));
    }
    h = _mm512_xor_si512(h, f);
    if (middle1 != 0) {
      f = _mm512_xor_si512(f,
                           _mm512_gf2p8affine_epi64_epi8(h, middle1_matrix, 0));
    }

    e = _mm512_xor_si512(e, a);
    if (top != 0) {
      a = _mm512_xor_si512(a, _mm512_gf2p8affine_epi64_epi8(e, top_matrix, 0));
    }
    f = _mm512_xor_si512(f, b);
    if (top != 0) {
      b = _mm512_xor_si512(b, _mm512_gf2p8affine_epi64_epi8(f, top_matrix, 0));
    }
    g = _mm512_xor_si512(g, c);
    if (top != 0) {
      c = _mm512_xor_si512(c, _mm512_gf2p8affine_epi64_epi8(g, top_matrix, 0));
    }
    h = _mm512_xor_si512(h, d);
    if (top != 0) {
      d = _mm512_xor_si512(d, _mm512_gf2p8affine_epi64_epi8(h, top_matrix, 0));
    }

    _mm512_storeu_si512(
        output0 + processed,
        _mm512_xor_si512(_mm512_loadu_si512(output0 + processed), a));
    _mm512_storeu_si512(
        output1 + processed,
        _mm512_xor_si512(_mm512_loadu_si512(output1 + processed), b));
    _mm512_storeu_si512(
        output2 + processed,
        _mm512_xor_si512(_mm512_loadu_si512(output2 + processed), c));
    _mm512_storeu_si512(
        output3 + processed,
        _mm512_xor_si512(_mm512_loadu_si512(output3 + processed), d));
    _mm512_storeu_si512(
        output4 + processed,
        _mm512_xor_si512(_mm512_loadu_si512(output4 + processed), e));
    _mm512_storeu_si512(
        output5 + processed,
        _mm512_xor_si512(_mm512_loadu_si512(output5 + processed), f));
    _mm512_storeu_si512(
        output6 + processed,
        _mm512_xor_si512(_mm512_loadu_si512(output6 + processed), g));
    _mm512_storeu_si512(
        output7 + processed,
        _mm512_xor_si512(_mm512_loadu_si512(output7 + processed), h));
  }

  for (size_t i = processed; i < byte_count; ++i) {
    Element a = x0[i];
    Element b = x1[i];
    Element c = x2[i];
    Element d = x3[i];
    Element e = x4[i];
    Element f = x5[i];
    Element g = x6[i];
    Element h = x7[i];
    ScalarIFFT(a, b, low0, tables);
    ScalarIFFT(c, d, low1, tables);
    ScalarIFFT(e, f, low2, tables);
    ScalarIFFT(g, h, low3, tables);
    ScalarIFFT(a, c, middle0, tables);
    ScalarIFFT(b, d, middle0, tables);
    ScalarIFFT(e, g, middle1, tables);
    ScalarIFFT(f, h, middle1, tables);
    ScalarIFFT(a, e, top, tables);
    ScalarIFFT(b, f, top, tables);
    ScalarIFFT(c, g, top, tables);
    ScalarIFFT(d, h, top, tables);
    output0[i] ^= a;
    output1[i] ^= b;
    output2[i] ^= c;
    output3[i] ^= d;
    output4[i] ^= e;
    output5[i] ^= f;
    output6[i] ^= g;
    output7[i] ^= h;
  }
}

const Kernels kKernels{IFFT, IFFTXor};

#endif

Status RunIFFT(const Context& context,
               std::span<Element* const> shards,
               size_t byte_count,
               size_t offset,
               size_t input_count,
               const detail::ResolvedKernels& base,
               const Kernels& radix8,
               std::span<Element* const> xor_accumulator = {},
               size_t initial_distance = 1) {
  if (input_count == 0 || byte_count == 0) {
    return Status::ok;
  }
  if (shards.size() == 1) {
    if (!xor_accumulator.empty()) {
      base.xor_one(xor_accumulator[0], shards[0], byte_count);
    }
    return Status::ok;
  }

  size_t distance = initial_distance;
  for (; 8 * distance <= shards.size(); distance *= 8) {
    const size_t group_size = 8 * distance;
    const size_t low_level = Log2(distance);
    const size_t middle_level = low_level + 1;
    const size_t top_level = low_level + 2;
    const bool fused = group_size == shards.size() && !xor_accumulator.empty();
    for (size_t block = 0; block < input_count; block += group_size) {
      const Element low0 = context.Skew(low_level, offset ^ block);
      const Element low1 =
          context.Skew(low_level, offset ^ (block + 2 * distance));
      const Element low2 =
          context.Skew(low_level, offset ^ (block + 4 * distance));
      const Element low3 =
          context.Skew(low_level, offset ^ (block + 6 * distance));
      const Element middle0 = context.Skew(middle_level, offset ^ block);
      const Element middle1 =
          context.Skew(middle_level, offset ^ (block + 4 * distance));
      const Element top = context.Skew(top_level, offset ^ block);
      for (size_t i = 0; i < distance; ++i) {
        if (fused) {
          radix8.ifft_radix8_xor(
              shards[block + i], shards[block + distance + i],
              shards[block + 2 * distance + i],
              shards[block + 3 * distance + i],
              shards[block + 4 * distance + i],
              shards[block + 5 * distance + i],
              shards[block + 6 * distance + i],
              shards[block + 7 * distance + i], xor_accumulator[block + i],
              xor_accumulator[block + distance + i],
              xor_accumulator[block + 2 * distance + i],
              xor_accumulator[block + 3 * distance + i],
              xor_accumulator[block + 4 * distance + i],
              xor_accumulator[block + 5 * distance + i],
              xor_accumulator[block + 6 * distance + i],
              xor_accumulator[block + 7 * distance + i], byte_count, top,
              middle0, middle1, low0, low1, low2, low3, context.Tables());
        } else {
          radix8.ifft_radix8(shards[block + i], shards[block + distance + i],
                             shards[block + 2 * distance + i],
                             shards[block + 3 * distance + i],
                             shards[block + 4 * distance + i],
                             shards[block + 5 * distance + i],
                             shards[block + 6 * distance + i],
                             shards[block + 7 * distance + i], byte_count, top,
                             middle0, middle1, low0, low1, low2, low3,
                             context.Tables());
        }
      }
    }
  }

  if (4 * distance <= shards.size()) {
    const bool fused = !xor_accumulator.empty();
    const Element top = context.Skew(Log2(distance) + 1, offset);
    const Element low = context.Skew(Log2(distance), offset);
    const Element high = context.Skew(Log2(distance), offset ^ (2 * distance));
    for (size_t i = 0; i < distance; ++i) {
      if (fused) {
        base.ifft_radix4_xor(shards[i], shards[distance + i],
                             shards[2 * distance + i], shards[3 * distance + i],
                             xor_accumulator[i], xor_accumulator[distance + i],
                             xor_accumulator[2 * distance + i],
                             xor_accumulator[3 * distance + i], byte_count, top,
                             low, high, context.Tables());
      } else {
        base.ifft_radix4(shards[i], shards[distance + i],
                         shards[2 * distance + i], shards[3 * distance + i],
                         byte_count, top, low, high, context.Tables());
      }
    }
    distance *= 4;
  }

  if (2 * distance <= shards.size()) {
    const Element coefficient = context.Skew(Log2(distance), offset);
    for (size_t i = 0; i < distance; ++i) {
      if (xor_accumulator.empty()) {
        base.ifft_radix2(shards[i], shards[distance + i], byte_count,
                         coefficient, context.Tables());
      } else {
        base.ifft_radix2_xor(shards[i], shards[distance + i],
                             xor_accumulator[i], xor_accumulator[distance + i],
                             byte_count, coefficient, context.Tables());
      }
    }
  }
  return Status::ok;
}

Status RunCopyFirst(const Context& context,
                    std::span<const Element* const> input,
                    std::span<Element* const> work,
                    size_t byte_count,
                    size_t offset,
                    const detail::ResolvedKernels& base,
                    const Kernels& radix8,
                    std::span<Element* const> xor_accumulator = {}) {
  if (byte_count == 0) {
    return Status::ok;
  }
  if (work.size() == 1) {
    if (xor_accumulator.empty()) {
      std::copy_n(input[0], byte_count, work[0]);
    } else {
      base.xor_one(xor_accumulator[0], input[0], byte_count);
    }
    return Status::ok;
  }
  if (work.size() == 2) {
    const Element coefficient = context.Skew(0, offset);
    if (!xor_accumulator.empty()) {
      base.ifft_radix2_xor(input[0], input[1], xor_accumulator[0],
                           xor_accumulator[1], byte_count, coefficient,
                           context.Tables());
      return Status::ok;
    }
    return base.ifft_radix2_copy(input[0], input[1], work[0], work[1],
                                 byte_count, coefficient, context.Tables());
  }

  const bool fused = work.size() == 4 && !xor_accumulator.empty();
  for (size_t block = 0; block < work.size(); block += 4) {
    const Element top = context.Skew(1, offset ^ block);
    const Element low = context.Skew(0, offset ^ block);
    const Element high = context.Skew(0, offset ^ (block + 2));
    if (fused) {
      base.ifft_radix4_xor(
          input[block], input[block + 1], input[block + 2], input[block + 3],
          xor_accumulator[block], xor_accumulator[block + 1],
          xor_accumulator[block + 2], xor_accumulator[block + 3], byte_count,
          top, low, high, context.Tables());
    } else {
      const Status status = base.ifft_radix4_copy(
          input[block], input[block + 1], input[block + 2], input[block + 3],
          work[block], work[block + 1], work[block + 2], work[block + 3],
          byte_count, top, low, high, context.Tables());
      if (status != Status::ok) {
        return status;
      }
    }
  }
  if (work.size() == 4) {
    return Status::ok;
  }
  return RunIFFT(context, work, byte_count, offset, work.size(), base, radix8,
                 xor_accumulator, 4);
}

}  // namespace

const Kernels* ResolveKernels() {
#if defined(__GFNI__) && defined(__AVX512F__) && defined(__AVX512BW__)
  return &kKernels;
#else
  return nullptr;
#endif
}

Status IFFTResolved(const Context& context,
                    std::span<Element* const> shards,
                    size_t byte_count,
                    size_t evaluation_offset,
                    size_t input_count,
                    const detail::ResolvedKernels& base,
                    const Kernels& radix8) {
  return RunIFFT(context, shards, byte_count, evaluation_offset, input_count,
                 base, radix8);
}

Status IFFTResolved(const Context& context,
                    std::span<const Element* const> input,
                    std::span<Element* const> work,
                    size_t byte_count,
                    size_t evaluation_offset,
                    const detail::ResolvedKernels& base,
                    const Kernels& radix8) {
  if (input.size() == work.size() && base.ifft_radix2_copy != nullptr &&
      base.ifft_radix4_copy != nullptr) {
    return RunCopyFirst(context, input, work, byte_count, evaluation_offset,
                        base, radix8);
  }
  if (byte_count == 0) {
    return Status::ok;
  }
  for (size_t i = 0; i < input.size(); ++i) {
    std::copy_n(input[i], byte_count, work[i]);
  }
  for (size_t i = input.size(); i < work.size(); ++i) {
    std::fill_n(work[i], byte_count, 0);
  }
  return RunIFFT(context, work, byte_count, evaluation_offset, input.size(),
                 base, radix8);
}

Status IFFTResolved(const Context& context,
                    std::span<const Element* const> input,
                    std::span<Element* const> work,
                    std::span<Element* const> xor_accumulator,
                    size_t byte_count,
                    size_t evaluation_offset,
                    const detail::ResolvedKernels& base,
                    const Kernels& radix8) {
  if (input.size() == work.size() && base.ifft_radix2_copy != nullptr &&
      base.ifft_radix4_copy != nullptr) {
    return RunCopyFirst(context, input, work, byte_count, evaluation_offset,
                        base, radix8, xor_accumulator);
  }
  if (byte_count == 0) {
    return Status::ok;
  }
  for (size_t i = 0; i < input.size(); ++i) {
    std::copy_n(input[i], byte_count, work[i]);
  }
  for (size_t i = input.size(); i < work.size(); ++i) {
    std::fill_n(work[i], byte_count, 0);
  }
  return RunIFFT(context, work, byte_count, evaluation_offset, input.size(),
                 base, radix8, xor_accumulator);
}

}  // namespace gf2p8::lch::detail::experiment::radix8
