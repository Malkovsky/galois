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

## LLVM Kernel-Screen Preflight

Use the LLVM 18 tools for a compute-bound LCH kernel screen and record the host
CPU before selecting the model:

```bash
command -v clang++-18 llvm-mca-18 perf taskset
clang++-18 --version
llvm-mca-18 --version
lscpu
```

Match the benchmark build's `-O3 -DNDEBUG -march=native` policy. Temporarily add
`LLVM-MCA-BEGIN` and `LLVM-MCA-END` inline-assembly comments around only the hot
loop, then generate and analyze the marked assembly:

```bash
MCA_CPU=${MCA_CPU:-native}
clang++-18 -std=c++20 -O3 -DNDEBUG -march=native \
  -I include -I src -S src/lin_chung_han/kernels.cc \
  -o /tmp/gf256_lch_kernels.s
llvm-mca-18 -mcpu="${MCA_CPU}" --iterations=100 \
  --bottleneck-analysis --resource-pressure --timeline \
  /tmp/gf256_lch_kernels.s
```

Inspect the assembly before modeling and reject new hot-loop spills or excessive
code growth. Treat throughput, dependency, and port-pressure output only as a
cheap static screen. Skip this step for memory-bound rows, complex control flow,
and complete encode/decode operations.

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

- `SharedShuffleRows`: row-FMA matrix multiply using shared packed tables.
- `BlockedLowHighSIMD`: blocked non-GFNI matrix kernel.
- `BlockedGFNIAffine`: blocked GFNI affine kernel with non-GFNI fallback.

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

With native references enabled, the C++ binary registers twelve top-level
families: owned LCH AVX2/GFNI256 encode/decode plus native Leopard, XDRS,
ISA-L, and Jerasure encode/decode. The standalone klauspost runner adds two
more families. Owned LCH, XDRS, ISA-L, Jerasure, and klauspost use the full
logarithmic grid from `8/248` through `248/8`; native Leopard uses its valid
`R<=K` subset from `128/128` through `248/8`. `DecodeMax` uses exactly `R=N-K`
deterministic shuffled erasures across data and recovery symbols. Decode
throughput uses the paper's input normalization, `K*bytes` per codeword, rather
than its separate recovered-output metric. ISA-L and Jerasure top-level decode
include matrix construction/inversion and data repair.

Use the exhaustive binary only when comparing kernels or tuning thresholds:

```bash
BENCH_CPU=${BENCH_CPU:-0}
taskset -c "${BENCH_CPU}" /tmp/gf256-linux/rs_verbose_benchmarks \
  --benchmark_filter='^LCH/(FFT|IFFT)/(Scalar|SSSE3|AVX2|GFNI(128|256|512)Affine)/(Radix2|Radix4)/transform:(2|4|8|16|32|64|128|256)/bytes:(1|8|15|16|32|64|128|256|512|1024|4096|65536)$' \
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
`RS/{Leopard,XDRS,ISA-L,Jerasure}/Native/{Encode,DecodeMax}` rows. This requires
x86 GNU/Clang and NASM because the pinned XDRS artifact and ISA-L kernels are
ISA-specific. Jerasure/GF-Complete are direct static source builds.

The repository preset configures and builds both concise and verbose binaries
with native ISA and all four C/C++ reference backends:

```bash
cmake --preset benchmarks
cmake --build --preset benchmarks
```

Build the standalone pinned klauspost runner when Go 1.24+ is available:

```bash
cmake --preset benchmarks -DGF256_GO_EXECUTABLE="$(command -v go)"
cmake --build build/benchmarks-preset --target klauspost_benchmark
build/benchmarks-preset/klauspost_rs_benchmarks --list
```

Compare every owned backend/radix row with the pinned native adapters:

```bash
taskset -c "${BENCH_CPU}" /tmp/gf256-linux/rs_verbose_benchmarks \
  --benchmark_filter='^LCH/Owned/|^(Leopard|XDRS)/Native/|^ISA-L/Native/' \
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
execute equal work. Use explicit event groups and split core and cache events to
prevent multiplexing:

```bash
PERF_BIN=${PERF_BIN:-perf}
BENCH_CPU=${BENCH_CPU:-0}
BENCH_ITERATIONS=${BENCH_ITERATIONS:-10000}
ROW='LCH/Owned/DecodeMax/AVX2/Radix4/K:128/R:128/bytes:1024'
for EVENTS in \
  '{cycles,instructions,branches,branch-misses}' \
  '{L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses}'; do
  "${PERF_BIN}" stat -x, -r 5 -e "${EVENTS}" -- \
    taskset -c "${BENCH_CPU}" /tmp/gf256-linux/rs_verbose_benchmarks \
      --benchmark_filter="^${ROW}$" \
      --benchmark_min_time="${BENCH_ITERATIONS}x" \
      --benchmark_min_warmup_time=0.1
done
```

Reject scaled/multiplexed or `<not counted>` output; split an oversized group and
rerun. Treat unsupported LLC events as a blocked cache-counter gate instead of
substituting generic cache events. Treat these counters and `llvm-mca` as
diagnostic screens; keep pinned `rs_verbose_benchmarks` timing and the concise
integrated `rs_benchmarks` suite authoritative.
