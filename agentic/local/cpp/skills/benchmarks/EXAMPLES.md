# GF256 Benchmark Examples

This repository has a matrix binary named `benchmarks`, a concise high-level RS
binary named `rs_benchmarks`, and an exhaustive tuning binary named
`rs_verbose_benchmarks`.

## Build

```bash
cmake -B /tmp/gf256-linux -DCMAKE_BUILD_TYPE=Release \
  -DGF256_ENABLE_NATIVE_ISA=ON
cmake --build /tmp/gf256-linux --target benchmarks rs_benchmarks rs_verbose_benchmarks gf_unittests lch_rs_unittests
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

## Focused LCH Timing

For ordinary encode/decode throughput, run the concise high-level suite:

```bash
BENCH_CPU=${BENCH_CPU:-0}
taskset -c "${BENCH_CPU}" /tmp/gf256-linux/rs_benchmarks \
  --benchmark_min_warmup_time=0.1 \
  --benchmark_repetitions=15 \
  --benchmark_min_time=0.1s \
  --benchmark_enable_random_interleaving=true \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true
```

With native references enabled, this registers twelve top-level families:
native/AVX2/GFNI256 encode/decode for both Leopard and XDRS. XDRS uses the
full logarithmic grid from `8/248` through `248/8`; Leopard uses its valid
`R<=K` subset from `128/128` through `248/8`. `DecodeMax` uses exactly `R=N-K` deterministic shuffled erasures
across data and recovery symbols. Decode throughput uses the paper's input
normalization, `K*bytes` per codeword, rather than its separate recovered-output
metric.

Use the exhaustive binary only when comparing kernels or tuning thresholds:

```bash
BENCH_CPU=${BENCH_CPU:-0}
taskset -c "${BENCH_CPU}" /tmp/gf256-linux/rs_verbose_benchmarks \
  --benchmark_filter='^LCH/(FFT|IFFT)/(Scalar|SSSE3|AVX2|GFNI(128|256|512)(Mul|Affine))/(Radix2|Radix4)/transform:(2|4|8|16|32|64|128|256)/bytes:(1|8|15|16|32|64|128|256|512|1024|4096|65536)$' \
  --benchmark_min_warmup_time=0.1 \
  --benchmark_repetitions=15 \
  --benchmark_min_time=0.1s \
  --benchmark_enable_random_interleaving=true \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_out=/tmp/gf256_lch.json \
  --benchmark_out_format=json
```

Build with `-DGF256_BUILD_REFERENCE_BENCHMARKS=ON` to add
`RS/{Leopard,XDRS}/Native/{Encode,DecodeMax}` rows. This requires x86
GNU/Clang because the pinned XDRS artifact is AVX2-only.

The repository preset configures and builds both concise and verbose binaries
with native ISA and both reference backends:

```bash
cmake --preset benchmarks
cmake --build --preset benchmarks
```

Compare every polished backend/radix row with the pinned native adapters:

```bash
taskset -c "${BENCH_CPU}" /tmp/gf256-linux/rs_verbose_benchmarks \
  --benchmark_filter='^(Leopard|XDRS)/(Polished|Native)/' \
  --benchmark_min_warmup_time=0.1 \
  --benchmark_repetitions=15 \
  --benchmark_min_time=0.1s \
  --benchmark_enable_random_interleaving=true \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_out=/tmp/gf256_rs_comparison.json \
  --benchmark_out_format=json
```

Use fixed iteration counts for Leopard perf-counter ratios so owned and native
execute equal work. Split core and cache events to avoid multiplexing:

```bash
PERF_BIN=${PERF_BIN:-perf}
BENCH_ITERATIONS=${BENCH_ITERATIONS:-10000}
ROW='Leopard/Polished/DecodeMax/AVX2/Radix4/K:128/R:128/bytes:1024'
for EVENTS in \
  cycles,instructions,branches,branch-misses \
  L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses; do
  "${PERF_BIN}" stat -x, -r 5 -e "${EVENTS}" -- \
    taskset -c "${BENCH_CPU}" /tmp/gf256-linux/rs_verbose_benchmarks \
      --benchmark_filter="^${ROW}$" \
      --benchmark_min_time="${BENCH_ITERATIONS}x" \
      --benchmark_min_warmup_time=0.1
done
```

Treat unsupported LLC events as a blocked cache-counter gate instead of
substituting generic cache events.
