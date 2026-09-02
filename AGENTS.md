# Repository Agent Notes

This repository is a small C++ GF arithmetic and matrix benchmark project.

## Shared Agent Reference

- Shared C++ workflows live in `agentic/cpp`.
- Do not put repository-specific paths, benchmark names, or commands inside
  `agentic/cpp`; keep that submodule generic.
- Put repository-local examples under
  `agentic/local/cpp/skills/<skill-name>/EXAMPLES.md`.
- When using a shared skill, read both the shared `SKILL.md` and the local
  `EXAMPLES.md` if present.

## Source Layout

- `include/field.h` and `src/field.cc` own field operations and row FMA
  primitives.
- `include/matrix.h` and `src/matrix.cc` own matrix multiplication APIs and
  kernels.
- Public RS APIs live under `include/reed_solomon`; implementations and private
  orchestration live under the mirrored `src/reed_solomon` tree.
- LCH public APIs live under `include/reed_solomon/lin_chung_han`; private
  dispatch and schedule declarations stay under `src/reed_solomon/lin_chung_han`.
- `benchmarks/matrix_multiplication.cc` owns Google Benchmark coverage.
- `tests/field_tests.cc` owns current correctness tests.
- `third_party/gf256` is a reference implementation; do not edit it unless the
  task explicitly targets that submodule.
- `third_party/gf256` emits SSSE3 intrinsics unconditionally; keep it out of
  scalar-only core targets or compile it separately with SSSE3 enabled.
- `third_party/xdrs` has global mutable setup, AVX2-only kernels with unsafe
  tails, and no `Algorithm4` definition; isolate it behind benchmark adapters.
- XDRS standalone builds need a forced `<cstring>` include; high-rate encode
  needs `2 * (256 - K)` parity buffers because its second half is scratch.
- `gf2p8::Element` is basis-neutral byte storage; scalar field APIs explicitly
  distinguish standard polynomial coordinates from Cantor coordinates.
- Owned LCH, matrix, Leopard, and XDRS code uses one immutable
  Cantor-coordinate domain.
  Native XDRS remains polynomial-coordinate and is not codeword-compatible.
- Native Cantor coordinates make XDRS derivative `B` scales identity, but its
  low-rate clear-and-XOR derivative still differs from Leopard's derivative.
- Field shuffle, affine, log, and exponent tables are compile-time-generated
  once in `src/field.cc`; no production field initialization is required.
- Leopard forms `FFTSkew - 1` sentinel pointers and triggers UBSan; sanitize the
  safe XDRS adapter rows and verify Leopard adapters with output-checked smoke.
- Leopard's speed depends on radix-4 fusion at every layer pair, truncated
  transforms, fused encoder IFFT/XOR, and sparse decode FFT; leaf-only fusion
  and full transforms lose most of the kernel-level GFNI gain.
- Leopard decode should scale-copy each work shard once. Clearing all work and
  then zeroing again inside multiply-copy adds a major full-shard memory pass.
- Iterative LCH coefficients are indexed by group start: low uses
  `Skew(level, offset ^ block)`, high uses `block + 2 * distance`, and top uses
  the next level. Prefix/sparse schedules must execute whole dependency groups.
- AVX2 LCH parity needs lane-duplicated 32-byte factor rows in `gf2p8::Tables`;
  per-kernel broadcasts and multiplying zero skew factors cost several percent.
- Validated RS loops should use resolved-backend scale/XOR internals; repeating
  public backend checks across every shard or derivative batch is measurable.
- Resolve ISA and byte-count kernel tables once per transform/RS operation;
  repeating enum dispatch inside every butterfly measurably raises branches.
- AVX2 rows divisible by 32 should use exact-vector radix entries; routing them
  through scalar-tail-capable entries costs roughly 2-4% retired instructions.
- Full-prefix immutable AVX2/GFNI IFFT should fuse its source copy into the
  first radix layer; separate copy and fold passes dominate RS encode.
- High-rate XDRS decode should fuse present/zero source-block IFFTs into the
  companion accumulator and sparse-evaluate only blocks with missing data.
- Leopard and high-rate XDRS encode share `rs::detail::EncodeLCH`; low-rate XDRS
  remains a separate coefficient fan-out algorithm.
- Concise owned/native RS rows share `benchmarks/rs_benchmark_cases.h`; XDRS
  uses the full log grid while Leopard uses only its `R<=K` subset.
- Pinned XDRS decode timing shuffles exactly `N-K` erasures across the full
  codeword; top-level `DecodeMax` must use the shared equivalent workload.
- The XDRS paper plots decode input throughput (`K * bytes / time`), not the
  artifact's separate recovered-output metric (`K * (N-K) / N`).
- The XDRS paper's Leopard comparison harness is absent from the released XDRS
  artifact; do not treat its plotted Leopard series as source-reproducible.

## Matrix/GF(256) Guidance

- Keep exact GF(256) semantics: field addition is XOR, not integer addition.
- For lookup-table SIMD paths, use the low/high-nibble `pshufb` decomposition:
  `product = table_lo[a][b & 0x0f] ^ table_hi[a][b >> 4]`.
- Native Cantor-coordinate GFNI uses generated fixed-factor affine matrices;
  direct `GF2P8MULB` uses AES polynomial `0x11B`, not standard `0x11D` or
  Cantor coordinates.
- GFNI paths must have correct non-GFNI fallbacks.
- Remove or isolate slow experimental variants before committing production
  kernels.
- Matrix changes should be validated against scalar/reference multiplication
  across empty, small, non-multiple, and wide shapes.

## Common Checks

```bash
cmake --build /tmp/gf256-linux --target gf_unittests benchmarks
/tmp/gf256-linux/gf_unittests
```

Compile fallback checks when touching SIMD gates:

```bash
g++ -std=c++20 -I include -I src -mno-gfni -mno-avx512f -mno-avx512bw -c src/field.cc -o /tmp/gf256-field-nogfni.o
g++ -std=c++20 -I include -I src -mno-gfni -mno-avx512f -mno-avx512bw -c src/matrix.cc -o /tmp/gf256-matrix-nogfni.o
g++ -std=c++20 -I include -I src -mno-avx2 -mssse3 -c src/field.cc -o /tmp/gf256-field-ssse3.o
g++ -std=c++20 -I include -I src -mno-avx2 -mssse3 -c src/matrix.cc -o /tmp/gf256-matrix-ssse3.o
g++ -std=c++20 -I include -I src -mno-avx2 -mno-ssse3 -c src/field.cc -o /tmp/gf256-field-scalar.o
g++ -std=c++20 -I include -I src -mno-avx2 -mno-ssse3 -c src/matrix.cc -o /tmp/gf256-matrix-scalar.o
```
