# Galois
Implementation of Galois fields $GF(2^8)$ and $GF(2^{16})$.

## $GF(2^8)$

### Base arithmetics

An element is represented by a basis-neutral `uint8_t`.
* Addition/subtraction is bitwise xor.
* `MultiplyStandard` interprets bytes as polynomial coefficients modulo
  $x^8+x^4+x^3+x^2+1$ (`0x11d`).
* `MultiplyCantor` interprets bytes as Cantor coordinates over the same field.
* Both multiplication forms use basis-specific logarithm tables:
```math
a\times b=\exp[\log[a]+\log[b]]
```
where $\exp$, $\log$ are constructed from primitive element $\alpha=x$ as
```math
\exp[i]=\alpha^i, \\
\log[\alpha^i]=i.
```
Each basis uses $2\times 256$ bytes. Division, inversion, and exponentiation
likewise have explicit `Standard` and `Cantor` variants.

### Vector operations

The row operation $\mathbf{x} += c\mathbf{y}$ uses one compile-time-generated,
immutable packed shuffle table with scalar, SSSE3, and AVX2 kernels.

GFNI uses generated fixed-coefficient `GF2P8AFFINEQB` matrices. Direct
`GF2P8MULB` is fixed to the AES polynomial `0x11b`, so it matches neither the
standard `0x11d` representation nor Cantor coordinates.

![Matrix multiplication benchmarks](https://malkovsky.github.io/galois/images/benchmarks.svg)

## $GF(2^{16})$

Implementation of $GF(2^{16})$ is extension over $GF(2^8)$ via polynomial $x^2+x+\delta$ where $\delta=x^5$ in polynomial coordinates (`0x20`) and `0xf0` in native Cantor coordinates.
* Multiplication is via fomula
```math
\begin{array}{rl} (a_0+a_1x)(b_0+b_1x)&=a_0b_0+(a_0b_1+a_1b_0)x+a_1b_1x^2\\&=a_0b_0+(a_0b_1+a_1b_0)x+a_1x_1(x+\delta) \\&=a_0b_0+a_1b_1\delta+(a_0b_1+a_1b_0+a_1b_1)x \end{array}
```
* Inverse is via powering and Itoh–Tsujii algorithm.
