# GF256 Optimization Experiment Examples

Use this overlay for GF(256) matrix and row-FMA experiments.

## Matrix Kernel Contract

- Public matrix APIs live in `include/matrix.h`.
- Matrix implementations live in `src/matrix.cc`.
- `src/field.*` should stay focused on field arithmetic and row FMA primitives.
- Keep production callers on `MatMul`, `MatMulBlockedLUT`, or
  `MatMulBlockedGFNI`; avoid leaving slow experimental variants in the public
  API.

## Correctness Loop

Before trusting timings, run:

```bash
cmake --build /tmp/gf256-linux --target gf_unittests benchmarks
/tmp/gf256-linux/gf_unittests
```

For SIMD gate changes, also compile representative fallback modes:

```bash
g++ -std=c++20 -I include -I src -mno-gfni -mno-avx512f -mno-avx512bw -c src/field.cc -o /tmp/gf256-field-nogfni.o
g++ -std=c++20 -I include -I src -mno-gfni -mno-avx512f -mno-avx512bw -c src/matrix.cc -o /tmp/gf256-matrix-nogfni.o
g++ -std=c++20 -I include -I src -mno-avx2 -mssse3 -c src/field.cc -o /tmp/gf256-field-ssse3.o
g++ -std=c++20 -I include -I src -mno-avx2 -mssse3 -c src/matrix.cc -o /tmp/gf256-matrix-ssse3.o
g++ -std=c++20 -I include -I src -mno-avx2 -mno-ssse3 -c src/field.cc -o /tmp/gf256-field-scalar.o
g++ -std=c++20 -I include -I src -mno-avx2 -mno-ssse3 -c src/matrix.cc -o /tmp/gf256-matrix-scalar.o
```

## Benchmark Loop

Compare candidate matrix kernels against the current `BlockedLowHighSIMD`,
`LowHighSIMDTables`, and `BlockedGFNI` rows:

```bash
BENCH_CPU=${BENCH_CPU:-0}
taskset -c "${BENCH_CPU}" /tmp/gf256-linux/benchmarks \
  --benchmark_filter='^(BlockedLowHighSIMD|LowHighSIMDTables|BlockedGFNI)/n:(512|1024|2048)$' \
  --benchmark_repetitions=5 \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_out=/tmp/gf256_matrix_candidate.json \
  --benchmark_out_format=json
```

Keep notes with the command, CPU pinning status, and JSON path when a result
drives an implementation choice.

## Compute-Bound LCH Kernel Screen

Preflight the pinned LLVM 18 tools and confirm the host CPU model before
screening an AVX2 or GFNI kernel:

```bash
command -v clang++-18 llvm-mca-18 perf taskset
clang++-18 --version
llvm-mca-18 --version
lscpu
llvm-mca-18 -mcpu=help
```

Temporarily bracket only the candidate's steady-state loop with
`asm volatile("# LLVM-MCA-BEGIN lch_kernel")` and
`asm volatile("# LLVM-MCA-END lch_kernel")`. Remove these markers after the
experiment. Generate assembly with the same native ISA policy as the benchmark
build, inspect it for spills and code growth, then model the marked loop:

```bash
MCA_CPU=${MCA_CPU:-native}
clang++-18 -std=c++20 -O3 -DNDEBUG -march=native \
  -I include -I src -S src/lin_chung_han/kernels.cc \
  -o /tmp/gf256_lch_kernels.s
llvm-mca-18 -mcpu="${MCA_CPU}" --iterations=100 \
  --bottleneck-analysis --resource-pressure --timeline \
  /tmp/gf256_lch_kernels.s
```

Use this screen for compact compute-bound radix kernels. Skip it for full
encode/decode orchestration, sparse data-dependent schedules, or rows dominated
by shard memory traffic. Advance a candidate only after correctness, no
unexplained spills or bloat, fixed-equal-work `perf stat`, and pinned
`rs_verbose_benchmarks` timing.
