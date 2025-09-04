# Galois
Implementation of Galois fields $GF(2^8)$ and $GF(2^{16})$.

## $GF(2^8)$

### Base arithmetics

Implementation of $GF(2^8)$ is pretty standard
* Element is represented by a `uint8_t`.
* Addition/subtraction is bitwise xor.
* Base multiplication is polynomial multiplication module $x^8+x^4+x^3+x+1$ which is `0x11b`. The choice is most common and is compatible with Intel GFNI.
* LUT multiplication is
```math
a\times b=\exp[\log[a]+\log[b]]
```
where $\exp$, $\log$ are constructed using primitive element $\alpha=x+1$ (which is $3$) as
```math
\exp[i]=\alpha^i, \\
\log[\alpha^i]=i.
```
$2\times 256$ bytes in total.
* Inversion via $\log$ table.

### Vector operations

Currently only implement operation $\mathbf{x} += c\mathbf{y}$ with $\mathbf{x}, \mathbf{y}\in \mathbb{F}^k, c\in \mathbb{F}$ in several variants for benchmraking

- Baseline: multiplication via binary multiplication tables
- SIMD: Intel implementation using high/low tables (taken from [catid/gf256](https://github.com/catid/gf256/))
- GFNIAffine/GFNIMul: multiplication via `GF2P8AFFINEQB`/`GF2P8MULB` instructions from [GFNI](https://builders.intel.com/docs/networkbuilders/galois-field-new-instructions-gfni-technology-guide-1-1639042826.pdf)

## $GF(2^{16})$

Implementation of $GF(2^{16})$ is extension over $GF(2^8)$ via polynomial $x^2+x+\delta$ where $\delta=x^5$ in $GF(2^8)$ (or `32`).
* Multiplication is via fomula
```math
\begin{array}{rl} (a_0+a_1x)(b_0+b_1x)&=a_0b_0+(a_0b_1+a_1b_0)x+a_1b_1x^2\\&=a_0b_0+(a_0b_1+a_1b_0)x+a_1x_1(x+\delta) \\&=a_0b_0+a_1b_1\delta+(a_0b_1+a_1b_0+a_1b_1)x \end{array}
```
* Inverse is via powering and Itoh–Tsujii algorithm.


