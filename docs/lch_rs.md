# LCH and Reed-Solomon Experiments

`gf2p8::Element` is basis-neutral byte storage. LCH, matrices, and owned RS
codecs use Cantor coordinates over GF(2^8)/`0x11D`; the scalar field API also
provides explicit standard polynomial-coordinate operations. One process-wide,
compile-time-generated field table set and one immutable LCH domain are shared
by transforms, matrices, and the owned RS encoder/decoder. The released native
XDRS artifact uses polynomial-coordinate bytes and a different evaluation
ordering, so it remains a benchmark reference rather than a compatible wire
format.

## Algorithms

- `gf2p8::rs::LCHEncoder` and `gf2p8::rs::LCHDecoder` expose separate data and
  recovery spans, so rate-specific transform layouts are private.
- For `R<K`, encoding and decoding use a generalized high-rate construction.
  Let `M=next_power_of_two(R)` and `N=next_power_of_two(M+K)`; supported
  dimensions satisfy `M+K<=256`. Positions `[R,M)` are punctured recovery and
  `[M+K,N)` are shortened known-zero data.
- For `R>=K`, encoding and decoding use a generalized low-rate construction
  derived from XDRS Algorithm 2. Let `D=next_power_of_two(K)`; supported
  dimensions satisfy `D+R<=256`. Positions `[K,D)` are shortened known zeros,
  and unused recovery positions in the rounded mother code are punctured.
- The high-rate companion decoder derives from XDRS Algorithm 3 and uses the
  logical `K` and `R` to truncate source blocks and sparse output evaluation.
- Low-rate decode is `O(N log D)`. High-rate decode fuses source-block IFFTs
  into its accumulator and sparse-evaluates only blocks containing missing
  outputs.
- Encode and decode perform no allocation. Callers size byte workspaces with
  each object's `WorkspaceSize(byte_count)`.
- Data, recovery, and workspace byte ranges must be pairwise disjoint. Missing
  data shards provide writable output pointers; missing recovery pointers may
  be null because decode repairs data only.

The tuned transform backend is selected statically from pinned Release
measurements: scalar below 16 bytes, SSSE3 for 16-31 bytes, AVX2 for 32-127
bytes, and 256-bit GFNI affine multiplication from 128 bytes when each ISA was
compiled. These thresholds were reconfirmed after the Cantor affine migration
against AVX2 and 128/256/512-bit GFNI implementations. Missing compiled
backends fall through to the next portable choice. Explicit backend/radix
variants remain available for transforms and RS operations. There is no
runtime CPU probing or autotuning.

Each public transform resolves its ISA and byte-count kernel table once. RS
encode/decode resolves once for the complete operation and reuses that table
across every transform and row primitive. AVX2/GFNI immutable-input IFFT fuses
a full-prefix source copy into its first radix layer; divisible-by-32 AVX2 rows
use exact-vector entries while arbitrary lengths retain scalar tails. Folded
encoding uses the copy-first path and fuses later block IFFTs into its XOR
accumulator.

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

Portable tests and benchmark compilation:

```bash
cmake -S . -B /tmp/gf256-lch -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/gf256-lch --target gf_unittests lch_rs_unittests rs_benchmarks rs_verbose_benchmarks
ctest --test-dir /tmp/gf256-lch --output-on-failure
```

Explicit AVX2/GFNI benchmark rows require a native-ISA build:

```bash
cmake -S . -B /tmp/gf256-lch-native \
  -DCMAKE_BUILD_TYPE=Release \
  -DGF256_ENABLE_NATIVE_ISA=ON
cmake --build /tmp/gf256-lch-native --target rs_benchmarks
```

Run the concise high-level encode/decode throughput suite:

```bash
BENCH_CPU=${BENCH_CPU:-0}
taskset -c "${BENCH_CPU}" /tmp/gf256-lch-native/rs_benchmarks \
  --benchmark_min_warmup_time=0.1 \
  --benchmark_repetitions=15 \
  --benchmark_min_time=0.1s \
  --benchmark_enable_random_interleaving=true \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true
```

With `GF256_BUILD_REFERENCE_BENCHMARKS=ON`, the concise binary has ten
top-level families: owned LCH AVX2/GFNI256 encode/decode, plus native
Leopard, XDRS, and ISA-L encode/decode. Owned LCH, XDRS, and ISA-L use the full
logarithmic grid, `K=8,16,32,64,128,192,224,240,248`; native Leopard uses only
the valid `R<=K` range, `K=128,192,224,240,248`. `DecodeMax` uses exactly
`R=N-K` deterministic shuffled erasures across data and recovery symbols,
matching the pinned XDRS workload shape.
Decode throughput uses the paper's input normalization, `K*bytes` per
codeword. The released XDRS benchmark also reports a separate recovered-output
rate, `K*R/N*bytes`, but that is not the metric plotted here.

The exhaustive `rs_verbose_benchmarks` binary remains available for backend,
vector-width, and radix tuning.

### Retained GFNI512 radix-8 experiment

The handwritten GFNI512 radix-8 IFFT is retained as a private experiment. It
fuses three IFFT layers over eight shards with seven hoisted affine matrices;
it never uses `GF2P8MULB`, because multiplication remains in native Cantor
coordinates. Unaligned 64-byte loads and stores are followed by exact scalar
tails, including zero-length and non-multiple-of-64 shards. The immutable-input
variant changes only its XOR outputs. Full immutable inputs still use the
existing radix-4 copy-first leaves before the experimental scheduler resumes at
distance four.

Production backend selection, `ResolvedKernels`, public `Radix`, public codec
APIs, and final radix-4 FFT scheduling are unchanged. Direct integrated
comparisons against the production GFNI256 path did not establish a stable
encoder win, despite strong isolated IFFT results, so this iteration retains
GFNI256 and radix-4 in production. An Intel measurement is still required
before reconsidering promotion. Ordinary builds contain no experiment objects,
tests, or rows. Enable the isolated private library, tests, and verbose rows
with the default-OFF
`GF256_ENABLE_GFNI512_RADIX8_EXPERIMENT` option or the dedicated preset:

```bash
cmake --preset gfni512-radix8-experiment
cmake --build --preset gfni512-radix8-experiment
ctest --preset gfni512-radix8-experiment
```

The exact benchmark families are:

```text
LCH/IFFTRadix8Leaf/GFNI512Affine/Radix4Control
LCH/IFFTRadix8Leaf/GFNI512Affine/Radix8Experiment
LCH/IFFT/GFNI512Affine/Radix4Control
LCH/IFFT/GFNI512Affine/Radix8Experiment
LCH/Experiment/Encode/GFNI512Affine/Radix4Control
LCH/Experiment/Encode/GFNI512Affine/Radix8
```

On a GFNI+AVX512F+AVX512BW Intel host, rerun the retained comparison pinned to
CPU 0 and preserve JSON evidence:

```bash
BENCH_CPU=${BENCH_CPU:-0}
taskset -c "${BENCH_CPU}" \
  build/gfni512-radix8-experiment-preset/rs_verbose_benchmarks \
  --benchmark_filter='^(LCH/(IFFTRadix8Leaf|IFFT)/GFNI512Affine/(Radix4Control|Radix8Experiment)|LCH/Experiment/Encode/GFNI512Affine/(Radix4Control|Radix8))/' \
  --benchmark_min_warmup_time=0.1 \
  --benchmark_repetitions=15 \
  --benchmark_min_time=0.1s \
  --benchmark_enable_random_interleaving=true \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_out=/tmp/gf256_gfni512_radix8_intel.json \
  --benchmark_out_format=json

for ROW in \
  'LCH/IFFTRadix8Leaf/GFNI512Affine/Radix4Control/bytes:1024' \
  'LCH/IFFTRadix8Leaf/GFNI512Affine/Radix8Experiment/bytes:1024'; do
  perf stat -r 5 -e '{cycles,instructions}' -- \
    taskset -c "${BENCH_CPU}" \
    build/gfni512-radix8-experiment-preset/rs_verbose_benchmarks \
      --benchmark_filter="^${ROW}$" --benchmark_min_time=10000x
done
```

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

The pinned native Leopard source constructs `FFTSkew - 1` sentinel pointers
and therefore emits a known UBSan diagnostic. Keep that reference unmodified,
do not claim its rows are sanitizer-clean, and rely on the adapter's untimed
decode/output verification. Owned rows and safe native XDRS rows remain
sanitizer-testable.

For the native-Leopard AVX2 parity gate, narrow the timing filter to the
explicit owned LCH AVX2 radix-4 and native rows while retaining the settings
above:

```bash
taskset -c "${BENCH_CPU}" build/benchmarks-preset/rs_verbose_benchmarks \
  --benchmark_filter='^(LCH/Owned/(Encode|DecodeMax)/AVX2/Radix4|Leopard/Native/(Encode|DecodeMax))/K:(32|129|128)/R:(16|64|128)/bytes:1024$' \
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
ROW='LCH/Owned/DecodeMax/AVX2/Radix4/K:128/R:128/bytes:1024'
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
adapter-level forced `<cstring>` include. ISA-L requires NASM and is built as a
static reference library. Neither is linked into production or ordinary test
targets. ISA-L uses its standard-coordinate systematic Cauchy code rather than
the owned LCH code family. Its top-level `DecodeMax` row includes survivor
selection, decode-matrix construction and inversion, table initialization, and
repair; `ISA-L/Native/RepairPrepared` is the verbose reconstruction-only row.
When owned AVX2 is compiled, verbose
`Kernel/{Scale,AddScaled}/{Owned,ISA-L}/AVX2` rows compare ISA-L's 32-byte
Cantor-table ABI against the project's already-resolved AVX2 entries.

## Attribution

The transform organization and Leopard-style orchestration derive from
Christopher A. Taylor's Leopard-RS, pinned under `third_party/leopard` at
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198` (BSD-3-Clause). The low-/high-rate
algorithm structure derives from Chao Chen's XDRS publication artifact, pinned
under `third_party/xdrs` at `ae05a779e7f44be13c3d34e14d15b08b4bc02404`
(Apache-2.0). The upstream source and license notices remain unmodified in
those submodules. Intel ISA-L is pinned under `third_party/isa-l` at
`7c3479e0a9dac17f448603ec1ad64c7c625f530c` (BSD-3-Clause).
