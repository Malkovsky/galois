#pragma once

#include "reed_solomon/lin_chung_han/kernels.h"

namespace gf2p8::lch::detail {

struct ResolvedKernels {
  using ScaleFunction = void (*)(Element*,
                                 const Element*,
                                 size_t,
                                 Element,
                                 const MultiplicationTables&);
  using XorFunction = void (*)(Element*, const Element*, size_t);
  using Xor4Function = void (*)(Element*,
                                const Element*,
                                Element*,
                                const Element*,
                                Element*,
                                const Element*,
                                Element*,
                                const Element*,
                                size_t);
  using Radix2Function = void (*)(Element*,
                                  Element*,
                                  size_t,
                                  Element,
                                  const MultiplicationTables&);
  using Radix2XorFunction = void (*)(const Element*,
                                     const Element*,
                                     Element*,
                                     Element*,
                                     size_t,
                                     Element,
                                     const MultiplicationTables&);
  using Radix2CopyFunction = Status (*)(const Element*,
                                        const Element*,
                                        Element*,
                                        Element*,
                                        size_t,
                                        Element,
                                        const MultiplicationTables&);
  using Radix4Function = void (*)(Element*,
                                  Element*,
                                  Element*,
                                  Element*,
                                  size_t,
                                  Element,
                                  Element,
                                  Element,
                                  const MultiplicationTables&);
  using Radix4XorFunction = void (*)(const Element*,
                                     const Element*,
                                     const Element*,
                                     const Element*,
                                     Element*,
                                     Element*,
                                     Element*,
                                     Element*,
                                     size_t,
                                     Element,
                                     Element,
                                     Element,
                                     const MultiplicationTables&);
  using Radix4CopyFunction = Status (*)(const Element*,
                                        const Element*,
                                        const Element*,
                                        const Element*,
                                        Element*,
                                        Element*,
                                        Element*,
                                        Element*,
                                        size_t,
                                        Element,
                                        Element,
                                        Element,
                                        const MultiplicationTables&);

  ScaleFunction add_scaled;
  ScaleFunction scale;
  XorFunction xor_one;
  Xor4Function xor_four;
  Radix2CopyFunction ifft_radix2_copy;
  Radix4CopyFunction ifft_radix4_copy;
  Radix2Function fft_radix2;
  Radix2Function ifft_radix2;
  Radix2XorFunction ifft_radix2_xor;
  Radix4Function fft_radix4;
  Radix4Function ifft_radix4;
  Radix4XorFunction ifft_radix4_xor;
};

const ResolvedKernels* ResolveKernels(Backend backend);
const ResolvedKernels* ResolveKernels(Backend backend, size_t byte_count);

Status ScaleUnchecked(Element* destination,
                      const Element* source,
                      size_t byte_count,
                      Element coefficient,
                      Backend backend,
                      const MultiplicationTables& tables);
void XorUnchecked(Element* destination,
                  const Element* source,
                  size_t byte_count,
                  Backend backend);
void Xor4Unchecked(Element* destination0,
                   const Element* source0,
                   Element* destination1,
                   const Element* source1,
                   Element* destination2,
                   const Element* source2,
                   Element* destination3,
                   const Element* source3,
                   size_t byte_count,
                   Backend backend);
void FFTRadix2Unchecked(Element* x,
                        Element* y,
                        size_t byte_count,
                        Element coefficient,
                        Backend backend,
                        const MultiplicationTables& tables);
void IFFTRadix2Unchecked(Element* x,
                         Element* y,
                         size_t byte_count,
                         Element coefficient,
                         Backend backend,
                         const MultiplicationTables& tables);
void IFFTRadix2XorUnchecked(const Element* x,
                            const Element* y,
                            Element* output_x,
                            Element* output_y,
                            size_t byte_count,
                            Element coefficient,
                            Backend backend,
                            const MultiplicationTables& tables);
void FFTRadix4Unchecked(Element* x0,
                        Element* x1,
                        Element* x2,
                        Element* x3,
                        size_t byte_count,
                        Element top,
                        Element low,
                        Element high,
                        Backend backend,
                        const MultiplicationTables& tables);
void IFFTRadix4Unchecked(Element* x0,
                         Element* x1,
                         Element* x2,
                         Element* x3,
                         size_t byte_count,
                         Element top,
                         Element low,
                         Element high,
                         Backend backend,
                         const MultiplicationTables& tables);
void IFFTRadix4XorUnchecked(const Element* x0,
                            const Element* x1,
                            const Element* x2,
                            const Element* x3,
                            Element* output0,
                            Element* output1,
                            Element* output2,
                            Element* output3,
                            size_t byte_count,
                            Element top,
                            Element low,
                            Element high,
                            Backend backend,
                            const MultiplicationTables& tables);

}  // namespace gf2p8::lch::detail
