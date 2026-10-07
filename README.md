# Galois

Practical library for some common algorithms over $\mathbb{F}_{2^8}$.

## Lin-Chung-Han variant of Reed-Solomon codes

The [Lin-Chung-Han (LCH) transform](https://arxiv.org/abs/1404.3458) is an
additive FFT-like transform over fields $\mathbb{F}_{2^m}$. For an $RS(n, k)$
code, it enables practical encoding and erasure decoding in
$\mathcal{O}(n\log\min(k, n-k))$ time. A canonical implementation is
[Leopard-RS](https://github.com/catid/leopard). This library provides:

- [XDRS](https://github.com/fastecc/xdrs)-derived low- and high-rate LCH decoders.
- Arbitrary positive data/recovery dimensions whose normalized power-of-two
  mother code fits 256 symbols.
- Fine-tuned AVX2 and [GFNI](https://builders.intel.com/docs/networkbuilders/galois-field-new-instructions-gfni-technology-guide-1-1639042826.pdf) kernels.

High-rate codewords use Leopard-compatible encoding.

### Benchmarks

Performance comparison of $RS(256, k)$ implementations:

![RS backend comparison](docs/images/rs_backend_comparison.svg)

[Open the benchmark report on GitHub Pages](https://malkovsky.github.io/galois/).

### Reed-Solomon product codes

We provide a special decoder for Reed-Solomon product codes over, error correction decoder for LCH aprroach is an implementation of a $\mathcal{O}((n-k)^2)$ variant of the modular approach by [Tang and Han](https://arxiv.org/pdf/2207.11079). See [Product code sampling experiments](manuals/rs-product-test.md) for details.

## Build

```bash
cmake --preset release
cmake --build --preset release
ctest --preset release
```

## Manuals

[Product code sampling experiments](manuals/rs-product-test.md)
