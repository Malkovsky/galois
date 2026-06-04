# GF256 Benchmark Examples

This repository has one Google Benchmark binary named `benchmarks`.

## Build

```bash
cmake -B /tmp/gf256-linux -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/gf256-linux --target benchmarks gf_unittests
```

## List Matrix Benchmarks

```bash
/tmp/gf256-linux/benchmarks --benchmark_list_tests=true
```

## Focused Matrix Timing

Use a narrow filter while iterating on matrix kernels:

```bash
BENCH_CPU=${BENCH_CPU:-0}
taskset -c "${BENCH_CPU}" /tmp/gf256-linux/benchmarks \
  --benchmark_filter='^(BlockedLowHighSIMD|LowHighSIMDTables|BlockedGFNI)/n:(512|1024|2048)$' \
  --benchmark_repetitions=5 \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_out=/tmp/gf256_matrix.json \
  --benchmark_out_format=json
```

If `taskset` fails on the host, rerun without it and report that the timing is
not CPU-pinned.

## Useful Benchmark Names

- `BinaryTable`: baseline scalar table row-FMA matrix multiply.
- `LowHighSIMDTables`: row-FMA matrix multiply using low/high-nibble SIMD LUTs.
- `BlockedLowHighSIMD`: blocked non-GFNI matrix kernel.
- `GFNIAffine`: matrix multiply using the general GFNI row FMA.
- `GFNIMul`: matrix multiply using the dedicated GFNI row FMA.
- `BlockedGFNI`: blocked GFNI matrix kernel with non-GFNI fallback.
