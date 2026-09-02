# LCH and Reed-Solomon Experiments

`gf2p8::Element` is basis-neutral byte storage. LCH, matrices, and owned RS
codecs use Cantor coordinates over GF(2^8)/`0x11D`; the scalar field API also
provides explicit standard polynomial-coordinate operations. One process-wide,
compile-time-generated field table set and one immutable LCH domain are shared
by transforms, matrices, Leopard, and XDRS. The released native XDRS artifact
uses polynomial-coordinate bytes and a different evaluation ordering, so owned
XDRS intentionally preserves its algorithms and layouts rather than artifact
codeword compatibility.

## Algorithms

- `gf2p8::rs::Leopard` implements shortened FF8 encode/decode with
  `R > 0`, `R <= K`, and `next_power_of_two(R) + K <= 256`.
- `gf2p8::rs::XDRS` implements the published full-length `N=256` paths.
  Low rate requires power-of-two `K <= 128`; high rate requires power-of-two
  `N-K <= 128`.
- Leopard and high-rate XDRS share the internal folded-IFFT `EncodeLCH`
  orchestration. Low-rate XDRS retains its distinct coefficient fan-out path.
- XDRS low-rate decode uses the blockwise `O(N log K)` construction. High-rate
  decode uses the companion-polynomial `O(N log (N-K))` construction, fuses
  source-block IFFTs into its accumulator, and sparse-evaluates only data
  blocks containing missing outputs.
- Encode and decode perform no allocation. Callers size byte workspaces with
  `EncodeWorkspaceSize(byte_count)` and `DecodeWorkspaceSize(byte_count)`.
- Leopard original, recovery, and workspace byte ranges must be pairwise
  disjoint. Decode likewise requires every codeword shard to be disjoint from
  its workspace; missing original shards still provide writable output buffers.
- XDRS encode likewise requires pairwise-disjoint data, recovery, and workspace
  byte ranges.

The tuned transform backend is selected statically from pinned Release
measurements: scalar below 16 bytes, SSSE3 for 16-31 bytes, AVX2 for 32-127
bytes, and 256-bit GFNI affine multiplication from 128 bytes when each ISA was
compiled. Missing compiled backends fall through to the next portable choice.
Explicit backend/radix variants remain available for transforms and RS
operations. There is no runtime CPU probing or autotuning.

Each public transform resolves its ISA and byte-count kernel table once. RS
encode/decode resolves once for the complete operation and reuses that table
across every transform and row primitive. AVX2/GFNI immutable-input IFFT fuses
a full-prefix source copy into its first radix layer; divisible-by-32 AVX2 rows
use exact-vector entries while arbitrary lengths retain scalar tails. Both
Leopard and XDRS use the copy-first path and fuse later block IFFTs into their
XOR accumulators.

## Transform Contracts

- Ordinary `FFT` and `IFFT` operate in place on a power-of-two shard span of
  size 1-256 at an aligned evaluation offset.
- Prefix `FFT(..., output_count, ...)` guarantees only the first
  `output_count` shards. Sparse `FFT(..., requested_outputs, ...)` guarantees
  only mask-selected outputs; all unrequested work is unspecified.
- Prefix-live in-place `IFFT(..., input_count, ...)` requires the suffix after
  `input_count` to be zero before the call.
- Immutable-input `IFFT(input, work, ...)` treats `input` as the live prefix,
  initializes the remaining work shards to zero, and leaves the transformed
  result in `work`.
- The three-span overload fuses the last layer into
  `xor_accumulator ^= IFFT(input)`. Its temporary `work` contents are
  unspecified afterward.
- Advanced input, work, and accumulator byte ranges must be pairwise disjoint.
  All overloads preserve arbitrary byte lengths, including zero, with exact
  scalar tails.

## Build and Benchmark

Portable tests and benchmarks:

```bash
cmake -S . -B /tmp/gf256-lch -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/gf256-lch --target gf_unittests lch_rs_unittests rs_benchmarks rs_verbose_benchmarks
ctest --test-dir /tmp/gf256-lch --output-on-failure
```

Run the concise high-level encode/decode throughput suite:

```bash
cmake --build /tmp/gf256-lch --target rs_benchmarks
BENCH_CPU=${BENCH_CPU:-0}
taskset -c "${BENCH_CPU}" /tmp/gf256-lch/rs_benchmarks \
  --benchmark_min_warmup_time=0.1 \
  --benchmark_repetitions=15 \
  --benchmark_min_time=0.1s \
  --benchmark_enable_random_interleaving=true \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true
```

With `GF256_BUILD_REFERENCE_BENCHMARKS=ON`, the concise binary has twelve
top-level families: native, owned AVX2, and owned GFNI256 encode/decode for
both Leopard and XDRS. XDRS uses its full logarithmic grid,
`K=8,16,32,64,128,192,224,240,248`; Leopard uses only the valid `R<=K` range,
`K=128,192,224,240,248`. `DecodeMax` uses exactly `R=N-K` deterministic shuffled erasures
across data and recovery symbols, matching the pinned XDRS workload shape.
Decode throughput uses the paper's input normalization, `K*bytes` per
codeword. The released XDRS benchmark also reports a separate recovered-output
rate, `K*R/N*bytes`, but that is not the metric plotted here.

The exhaustive `rs_verbose_benchmarks` binary remains available for backend,
vector-width, and radix tuning.

Compile host ISA variants and optional native references on x86 GNU/Clang:

```bash
cmake --preset benchmarks
cmake --build --preset benchmarks
BENCH_CPU=${BENCH_CPU:-0}
taskset -c "${BENCH_CPU}" build/benchmarks-preset/rs_verbose_benchmarks \
  --benchmark_filter='^LCH/(FFT|IFFT)/(Scalar|SSSE3|AVX2|GFNI(128|256|512)Affine)/(Radix2|Radix4)/transform:(2|4|8|16|32|64|128|256)/bytes:(1|8|15|16|32|64|128|256|512|1024|4096|65536)$' \
  --benchmark_min_warmup_time=0.1 \
  --benchmark_repetitions=15 \
  --benchmark_min_time=0.1s \
  --benchmark_enable_random_interleaving=true \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_out=/tmp/gf256_lch_tuning.json \
  --benchmark_out_format=json

taskset -c "${BENCH_CPU}" build/benchmarks-preset/rs_verbose_benchmarks \
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

The pinned native Leopard source constructs `FFTSkew - 1` sentinel pointers
and therefore emits a known UBSan diagnostic. Keep that reference unmodified,
do not claim its rows are sanitizer-clean, and rely on the adapter's untimed
decode/output verification. Owned rows and safe native XDRS rows remain
sanitizer-testable.

For the Leopard AVX2 parity gate, narrow the timing filter to the explicit
owned AVX2 radix-4 and native rows while retaining the settings above:

```bash
taskset -c "${BENCH_CPU}" build/benchmarks-preset/rs_verbose_benchmarks \
  --benchmark_filter='^Leopard/(Polished/(Encode|DecodeMax)/AVX2/Radix4|Native/(Encode|DecodeMax))/K:(32|129|128)/R:(16|64|128)/bytes:1024$' \
  --benchmark_min_warmup_time=0.1 \
  --benchmark_repetitions=15 \
  --benchmark_min_time=0.1s \
  --benchmark_enable_random_interleaving=true \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_out=/tmp/gf256_leopard_avx2_parity.json \
  --benchmark_out_format=json
```

Use equal fixed iteration counts when comparing hardware counters. Run each
owned/native row separately and repeat once for each event group; splitting the
groups avoids counter multiplexing. Instructions divided by cycles gives IPC.
If a distribution's `perf` wrapper does not match its running kernel, set
`PERF_BIN` to the installed real binary under `/usr/lib/linux-tools-*`.

```bash
PERF_BIN=${PERF_BIN:-perf}
BENCH_ITERATIONS=${BENCH_ITERATIONS:-10000}
ROW='Leopard/Polished/DecodeMax/AVX2/Radix4/K:128/R:128/bytes:1024'
for EVENTS in \
  cycles,instructions,branches,branch-misses \
  L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses; do
  "${PERF_BIN}" stat -x, -r 5 -e "${EVENTS}" -- \
    taskset -c "${BENCH_CPU}" build/benchmarks-preset/rs_verbose_benchmarks \
      --benchmark_filter="^${ROW}$" \
      --benchmark_min_time="${BENCH_ITERATIONS}x" \
      --benchmark_min_warmup_time=0.1
done
# Repeat with ROW='Leopard/Native/DecodeMax/K:128/R:128/bytes:1024'.
```

If the LLC events report `<not supported>`, record the cache-counter portion of
the acceptance gate as blocked rather than substituting generic cache events.

Reference adapters are benchmark-only. XDRS is compiled with AVX2 and an
adapter-level forced `<cstring>` include; it is never linked into production
or ordinary test targets.

## Attribution

The transform organization and Leopard-style orchestration derive from
Christopher A. Taylor's Leopard-RS, pinned under `third_party/leopard` at
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198` (BSD-3-Clause). The low-/high-rate
algorithm structure derives from Chao Chen's XDRS publication artifact, pinned
under `third_party/xdrs` at `ae05a779e7f44be13c3d34e14d15b08b4bc02404`
(Apache-2.0). The upstream source and license notices remain unmodified in
those submodules.
