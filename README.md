# Galois

Practical library for some common algorithms over $\mathbb{F}_{2^8}$.

## Lin-Chung-Han variant of Reed-Solomon codes

[Lin-Chung-Han](https://arxiv.org/abs/1404.3458) transform is the basis for several practical algorithms
in fields $\mathbb{F}_{2^m}$, most notably Reed-Solomon codes with the canonical implementation in [Leopard-RS](https://github.com/catid/leopard). We provide

- XDRS-derived low- and high-rate LCH decoders.
- Arbitrary positive data/recovery dimensions whose normalized power-of-two
  mother code fits 256 symbols.
- Fine-tuned AVX2 and [GFNI](https://builders.intel.com/docs/networkbuilders/galois-field-new-instructions-gfni-technology-guide-1-1639042826.pdf) kernels.

High-rate codewords use Leopard-compatible encoding.

### Benchmarks

Comparison of $RS(256, k)$ code with common libraries

![RS backend comparison](docs/images/rs_backend_comparison.svg)

[Open the benchmark report on GitHub Pages](https://malkovsky.github.io/galois/).


## Build

```bash
cmake --preset release
cmake --build --preset release
ctest --preset release
```
