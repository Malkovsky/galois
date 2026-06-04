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

- `src/field.h` and `src/field.cc` own field operations and row FMA primitives.
- `src/matrix.h` and `src/matrix.cc` own matrix multiplication APIs and kernels.
- `benchmarks/matrix_multiplication.cc` owns Google Benchmark coverage.
- `tests/field_tests.cc` owns current correctness tests.
- `third_party/gf256` is a reference implementation; do not edit it unless the
  task explicitly targets that submodule.

## Matrix/GF(256) Guidance

- Keep exact GF(256) semantics: field addition is XOR, not integer addition.
- For lookup-table SIMD paths, use the low/high-nibble `pshufb` decomposition:
  `product = table_lo[a][b & 0x0f] ^ table_hi[a][b >> 4]`.
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
g++ -std=c++20 -I src -I third_party -mno-gfni -mno-avx512f -mno-avx512bw -c src/field.cc -o /tmp/gf256-field-nogfni.o
g++ -std=c++20 -I src -I third_party -mno-gfni -mno-avx512f -mno-avx512bw -c src/matrix.cc -o /tmp/gf256-matrix-nogfni.o
g++ -std=c++20 -I src -I third_party -mno-avx2 -mssse3 -c src/field.cc -o /tmp/gf256-field-ssse3.o
g++ -std=c++20 -I src -I third_party -mno-avx2 -mssse3 -c src/matrix.cc -o /tmp/gf256-matrix-ssse3.o
g++ -std=c++20 -I src -I third_party -mno-avx2 -mno-ssse3 -c src/field.cc -o /tmp/gf256-field-scalar.o
g++ -std=c++20 -I src -I third_party -mno-avx2 -mno-ssse3 -c src/matrix.cc -o /tmp/gf256-matrix-scalar.o
```
