# Galois
Implementation of Galois fields $GF(2^8)$ and $GF(2^{16})$.

## Details

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
* Addtionally implemented GFNI-like multiplication with precomputed $8\times 8$ matrices packed into `unit64_t` totalling in $8\times 256$ byte. The multiplication itself is manual matrix-vector multiplication.
* Inversion via $\log$ table.

Implementation of $GF(2^{16})$ is extension over $GF(2^8)$ via polynomial $x^2+x+\delta$ where $\delta=x^5$ in $GF(2^8)$ (or `32`).
* Multiplication is via fomula
```math
\begin{array}{rl} (a_0+a_1x)(b_0+b_1x)&=a_0b_0+(a_0b_1+a_1b_0)x+a_1b_1x^2\\&=a_0b_0+(a_0b_1+a_1b_0)x+a_1x_1(x+\delta) \\&=a_0b_0+a_1b_1\delta+(a_0b_1+a_1b_0+a_1b_1)x \end{array}
```
* Inverse is via powering and Itoh–Tsujii algorithm.

## Currently not implemented
* No actual GFNI. `GF2P8MULB` instruction is by far the fastest way of multiplication in $GF(2^8)$ if available, expecially batched.
* No SIMD optimizations. This is a clear miss as most applications operate on sequences, as far as I know used multiplication implementations are poorly optimized by the compiler. A special mention should be made to [catid/gf256](https://github.com/catid/gf256/blob/master/gf256.cpp#L471) with heavy manual SIMD optimizations and dedicated SIMD-based multiplication algorithm.
* No benchmarks. Will add them probably at some point. 

