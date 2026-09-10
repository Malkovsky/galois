# Strong-Weak Reed-Solomon Product Code Acceleration Report

## Applied R=4 Weak Code (2026-09-10)

**Implemented and retained:** `StrongWeakRSProductCode(256,224,256,252)` now
encodes and corrects with a direct two-error weak decoder. This section supersedes
the earlier statements that R=4 was deferred. Strong kernels and behavior are
unchanged. Default construction remains 256/224 x 256/254, and all existing R=2
shortenings, including 175/173, retain their coordinate mapping and behavior.
R=4 support is deliberately restricted to full RS(256,252); smaller or shortened
R=4 codes, including the low-rate K=R boundary, are rejected rather than assumed
equivalent. Other already-supported strong dimensions can use this full weak code.

### Exact Native Parity Checks

All arithmetic below is in native Cantor GF(256); addition is XOR. Public weak
rows are `[252 data][4 parity]`. Native evaluation points for public data index
`j` are `j+4`, and parity index `252+i` has point `i`. Thus public index `p`
maps to native byte `(p+4) mod 256`. This is an index permutation, not field
addition by integer 4.

The full native LCH code is evaluation of polynomials of degree at most 251.
Changing from its monic novel basis to ordinary monomials preserves that
polynomial space. For the full field, the vanishing polynomial is
\(P(X)=X^{256}+X\), with \(P'(X)=1\). Consequently all dual evaluation weights
are one. Equivalently, \(\sum_{x\in GF(256)}x^m=0\) for \(0\leq m\leq254\),
including \(m=0\), since 256 is zero in characteristic two. Therefore these
four ordinary power moments vanish on every encoded row:

\[
 S_j=\sum_x c(x)x^j,\qquad j=0,1,2,3.
\]

The four check rows have rank four: a nonzero polynomial of degree at most
three cannot vanish at all 256 distinct points. Their kernel thus has dimension
252 and is **exactly** the native code, not merely a necessary validity test.
These moments are not asserted to be the individual novel-basis syndrome
entries used internally by generic FDMA. They are an independent complete
parity-check system for the same code. Scalar LCH encoding and the untouched
generic `CorrectCodeword(LCHDecoder(252,4), ...)` supply independent checks.

### Locator, Magnitudes, and Degeneracies

For errors at distinct points \(x_1,x_2\) with nonzero magnitudes \(e_1,e_2\),
\(S_j=e_1x_1^j+e_2x_2^j\). Write the locator as
\(L(X)=X^2+aX+b\), where \(a=x_1+x_2\) and \(b=x_1x_2\). Its recurrence gives

\[
 \begin{pmatrix}S_1&S_0\\S_2&S_1\end{pmatrix}
 \begin{pmatrix}a\\b\end{pmatrix}
 =\begin{pmatrix}S_2\\S_3\end{pmatrix},\qquad
 D=S_1^2+S_0S_2=e_1e_2(x_1+x_2)^2.
\]

Thus a genuine two-error pattern always has nonzero determinant, even when
equal magnitudes make \(S_0=0\). No division by \(S_0\) occurs on this branch:

\[
 a=(S_1S_2+S_0S_3)/D,\qquad
 b=(S_1S_3+S_2^2)/D.
\]

For \(a\ne0\), substitute \(X=ay\) and solve
\(y^2+y=b/a^2\). The Artin-Schreier map has kernel \(\{0,1\}\) and its
128-element image consists exactly of trace-zero field elements. A private
512-byte immutable table stores a representative for each solvable value and
`-1` for insoluble values, keeping solvable zero distinct from failure. It is
generated once, thread-safely, using native `MultiplyCantor`; there is no global
mutable field setup or large multiplication table. The roots are \(ay\) and
\(ay+a\), and their magnitudes are

\[
 e_1=(S_1+S_0x_2)/a,\qquad e_2=S_0+e_1.
\]

The implementation handles every branch explicitly:

- All four moments zero: already a codeword, not proof of original content.
- Nonzero syndrome with `D=0`: require `S0!=0`, propose `e=S0`, `x=S1/S0`,
  and verify all four moments. Rank-one-looking but inconsistent moments fail.
- `D!=0` with `a=0`: reject the repeated-root locator; square roots cannot
  produce two distinct error positions.
- Insoluble quadratic or zero proposed magnitude: reject.
- For either candidate size, XOR its contribution out of all four moments and
  require zero before reporting success. No input bytes are modified here.

All byte-valued roots belong to this full mother code; there are no omitted
positions. That fact is specific to the retained full-length R=4 scope. The
existing R=2 virtual-position rejection remains unchanged. Distance five makes
a candidate within radius two unique, but over-radius received words can still
lie within radius two of a different codeword. The direct and generic paths
preserve that BDD behavior; neither promises detection of all larger errors.

### Integration and Coverage

`WeakCandidateR4` in `src/reed_solomon/strong_weak_rs_product_code.cc` performs
three native table products per received symbol to accumulate four moments.
It performs no weak transposition or full candidate copy, and final weak
validation requests moments only, without locator solving. R=4 encoding uses
the existing `LCHEncoder(252,4)` on each information row before the unchanged
strong encode; no closed-form encoder speedup is claimed.

The scheduler checks **every** proposed byte against the optional binary-image
gate (`popcount(delta)<=2`) and its protected column before committing **any**
byte in that row. Rejection is transactional for the entire candidate. Accepted
writes update exact bit/symbol counters, invalidate intersecting column validity,
and activate those columns. Generic reference scheduling now accepts one or two
weak repairs for R=4 while R=2 still accepts only one.

Tests cover 65,280 single errors; all 32,640 position pairs with magnitudes
`(1,1)` and `(3,128)`; every normalized Artin-Schreier input, both soluble and
insoluble; explicit inconsistent rank-one and repeated-root cases; and 4,096
random one-through-nine-injection comparisons against generic BDD. Exhaustive
position tests use independently scalar-encoded nonzero rows as their oracle.
All Artin-Schreier and randomized cases compare generic status, count, and full
output. Product scheduler differential coverage includes strong 4/2 and 256/224
with weak 256/252, all gate combinations, caps 2/3/4/5/6/16, generic single and
batch paths, public direct correction, and every retained optimization bitset.
It compares final bytes, independent component validity, and every result field.
Dedicated tests reject both repairs if only one delta exceeds two bits or only
one target is protected, including data/parity targets.

Native CLI dimension validation/help, atomic summary metadata, report/replay,
and plotting accept the new dimensions without changing schema or old defaults.
R=2 and R=4 reports never pool into one code group. Random-message native trials
already route nondefault dimensions through `code->Encode`, so the legacy
default-only R=2 encoder helper requires no change. Extended tests check this
route against all-zero trials, all gate combinations and saved-position replay.
Floyd remains thread-count reproducible; Fisher-Yates intentionally retains
worker-local permutations and is tested by saved-position replay instead.

```cpp
gf2p8::rs::StrongWeakRSProductCode code(256, 224, 256, 252);
// code.Encode(block), then code.Correct(block, options), as before.
```

```bash
rs-product-monte-carlo --output /tmp/kilo/mc-r4 --n2 256 --k2 252 --seed 42 --batches 1 --batch-size 64 --threads 1
```

### Matched R=4 Measurements

AMD Ryzen 7 8845HS, WSL2, GCC 13.3.0, Google Benchmark 1.9.4, Release
`-O3 -DNDEBUG`. Native uses `-march=native`; AVX2-only disables GFNI and AVX-512.
Runs were CPU-0 pinned, sequential, three randomly interleaved repetitions of
at least 0.3s, with `CPP_JOBS=2`. Process inspection showed active editor indexing
and agent processes, but no Monte Carlo worker. Nothing was killed or re-pinned.
Host load and clocks are uncontrolled; these are bounded local measurements,
not Intel-host or multithread scaling evidence.

Both rows use the same 64-block BSC(0.005) corpus, seed `0x5b5c0224`, strong
256/224, weak 256/252, cap 16, both gates enabled. Corpus setup compares every
output byte and result field with generic scalar/full-validation scheduling.
Timing includes correction and final validity, excludes encode/reset/accounting,
and normalizes to **224*252 = 56,448 information bytes per block**.

| Profile | Generic batch CPU ms / 64 blocks | Direct CPU ms / 64 blocks | Generic MiB/s | Direct MiB/s | Throughput gain | Generic / direct CPU CV |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Native | 20.290 | 17.756 | 169.801 | 194.033 | 14.3% | 0.48% / 2.30% |
| AVX2-only | 23.501 | 20.852 | 146.600 | 165.230 | 12.7% | 3.38% / 0.34% |

All measured corpus residual bits and message/block/validity/valid-wrong/pass-limit
failures were zero. Both paths had mean 3.734375 passes, 262.218750 strong lines,
258.765625 weak lines, 2570.328125 byte writes and 2615.734375 bit toggles;
weak-only writes/toggles were 103.453125/105.453125. These are corpus results,
not error-floor estimates. R=2 numbers elsewhere describe a different code and
are not the denominator of these speedup claims. No new instruction-count,
assembly-model, or isolated arithmetic throughput claim is made.

Reproduction, with exported `CPP_JOBS=2 CPP_BENCH_CPU=0`,
`BENCHMARK_OUT_FORMAT=json` and distinct `BENCHMARK_OUT` paths:

```bash
/home/user/.config/kilo/scripts/cpp-bench native matrix '^LCH/Owned/StrongWeakRSProductCode/R4/(Direct|Generic)256252$' 3 0.3s
/home/user/.config/kilo/scripts/cpp-bench avx2 matrix '^LCH/Owned/StrongWeakRSProductCode/R4/(Direct|Generic)256252$' 3 0.3s
```

Artifacts: `/tmp/kilo/product-r4-native.json` and
`/tmp/kilo/product-r4-avx2.json`. `/check quick` passed portable and native
unfiltered CTest suites; `/check avx2` passed its unfiltered suite. Focused
ASan/UBSan `/check sanitize` with
`GTEST_FILTER='ProductCode.*:WholeCodewordBatch.*'` passed (204.40s selected RS
tests, 268.02s total including CLI/legacy/plot integration). This sanitizer run
excludes unrelated GTest suites and is not described as a full sanitizer pass.
`/check cli` passed all seven tests, including the unrelated fragmenter CLI.
`/check fallback` passed; because that wrapper omits the product source, these
additional strict compilations also passed:

```bash
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -I include -I src -mno-avx2 -mno-ssse3 -mno-gfni -mno-avx512f -mno-avx512bw -c src/reed_solomon/strong_weak_rs_product_code.cc -o /tmp/kilo/product-r4-scalar.o
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -I include -I src -mavx2 -mssse3 -mno-gfni -mno-avx512f -mno-avx512bw -c src/reed_solomon/strong_weak_rs_product_code.cc -o /tmp/kilo/product-r4-avx2.o
```

Formatting uses clang-format 18.1.3; `git diff --check` passed. No commits made.

## Post-Profile Experiments (2026-09-10)

All three requested experiments were implemented, differentially tested and
measured. Retained: sparse strong-mask traversal, coordinate-bit weak syndrome
reduction, and padding independent strong batch lanes. Public APIs, gates,
transactionality, stopping rules, counters, mother codes and R=2 policy are
unchanged. The earlier uncommitted direct-R2 implementation and the report below
are preserved. No commits were made.

### Retained Techniques and Limits

1. **Strong masks:** AVX2 compares 32 mask bytes against zero and extracts a
   bitset. Set-bit traversal visits only verified repairs; old/new XOR popcount,
   byte write counts, next-direction activation and clean-state invalidation
   still occur for every committed byte. Real partial tails use the exact scalar
   scan. Runtime availability is checked once per product correction, not once
   per position row. Staging remains transactional per lane; failed lanes are
   not committed. This is a modest, workload-dependent improvement, not removal
   of the entire old 19-23% profiled region.
2. **Weak reduction:** for each complete 32-byte group, XOR symbols into eight
   accumulators selected by the coordinate bits of weight `(j+2) XOR 1`.
   Linearity gives `S0 = XOR_b MultiplyCantor(1<<b, XOR_selected_symbols_b)`.
   `S1` is the unweighted XOR. Thus the vectorized prefix needs only eight
   fixed-factor scalar table products after horizontal reduction. A 2 KiB
   compile-time mask array selects positions; there is no mutable field setup,
   lane-varying fixed-factor shuffle, GFNI matrix misuse, or basis conversion.
   Remaining symbols retain the old scalar table loop. The helper serves
   encoding, correction and final weak validation. N=4 stays scalar with the
   original parity weights; all shortened lengths keep their original weights.
3. **Strong tail:** staging rounds only the batch lane stride to 32. For width
   175, five full chunks plus 15 scalar columns become six full chunks:
   **15 real tail lanes plus 17 synthetic zero codewords**, not 15 synthetic
   lanes. Only 175 real outcomes, masks and bytes are published. Combined packed
   and mask storage grows from 87.5 to 96 KiB at strong N=256; allocation and the
   strided copy are included in timing. Width 256 remains 128 KiB with a single
   contiguous copy. Generic mother-code differential paths are not lane-padded.

The private `ProductCorrectionAccess::Experiment` selects a per-call bitset:
0 = direct-R2 baseline, 1 = sparse masks, 2 = vector reduction, 4 = padded lanes;
3/6/7 combine the corresponding bits. Public correction enables all three.
No environment switch, shared mutable setting, workspace API or public benchmark
API was introduced. These small private controls and rows remain for reproducible
same-executable comparison. Both the vector reduction and sparse scan have
compile-time AVX2 guards and runtime availability checks. The existing generic
batch kernel, including its scalar fallback, is unchanged.

The first mask variant repeated backend availability checks per row; that was
removed before the confirmation runs. Alternative GFNI/AES dot products,
lane-varying multiply tables and a generic scratch-padded `CorrectChunk32` tail
were not implemented or measured. The latter is deferred because product-local
staging already addresses the requested 175 case without modifying generic BDD
scratch or validation. No losing arithmetic implementation is left in production.
No R=4, allocation-reuse or Monte Carlo algorithm changes were made.

### Matched Measurements

Same host and flags as the post-R2 profile below: Ryzen 7 8845HS, WSL2, GCC
13.3.0, Release `-O3 -DNDEBUG`, native `-march=native`; AVX2-only disables GFNI
and AVX-512. CPU 0 pinned, `CPP_JOBS=2`. Process inspection before measurement
found no active product-code Monte Carlo worker; `mc` was idle and two Kilo
processes were active. No process was stopped or had its affinity changed.
Windows host load and actual clocks remain uncontrolled.

Every row uses the existing deterministic 64-block BSC(0.005) corpus,
seed `0x5b5c0224`, strong 256/224, cap 16, both weak gates enabled. Throughput
counts **224*Kweak information bytes**, not encoded bytes. Benchmark setup
checks full output bytes and all result fields against generic scalar mother
BDD. Reset and accounting are paused; correction includes final validity.
Encoding setup is not timed. Results below are median CPU ms per 64-block sweep.

Before editing, exact public default rows (3 x 0.3s) measured 32.037 / 28.614 ms,
108.396 / 82.659 MiB/s, CPU CV 5.79% / 7.59% (256 / 175). Subsequent comparisons
use the private direct-R2 baseline in the **same executable**, not historical
generic-R2 timings or these drifting initial measurements.

| Variant | First native 256 ms | First native 175 ms | Confirm native 256 ms | Confirm native 175 ms | AVX2-only 256 ms | AVX2-only 175 ms |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0: Direct-R2 baseline | 27.898 | 26.840 | 31.480 | 27.849 | 33.967 | 37.731 |
| 1: Masks only | 27.025 | 23.022 | 30.158 | 26.278 | 34.207 | 34.174 |
| 2: Weak reduction only | 22.556 | 21.466 | 24.702 | 24.420 | 30.090 | 32.069 |
| 3: Masks + weak reduction | not run | not run | 24.359 | 21.621 | not run | not run |
| 4: Padded lanes only | 28.208 | 20.685 | 31.675 | 21.174 | 36.866 | 26.687 |
| 6: Weak reduction + padding | not run | not run | 25.389 | 18.457 | 32.455 | 22.922 |
| 7: All three, private row | 22.536 | 15.175 | 21.592 | 17.546 | not run | not run |
| All three, public default | not run | not run | 21.922 | 17.006 | 29.031 | 21.921 |

First native: 3 x 0.3s, CPU CV 1.30-7.03%. Confirmation: 5 x 0.5s, CPU CV
3.51-28.55%; baseline CV 19.74% / 5.56%, public default CV 22.61% / 18.25%.
AVX2-only: 3 x 0.3s, CPU CV 0.84-9.37%. All rows randomly interleaved by the
wrapper. The noisy confirmation is retained here, not silently discarded.
Padding does no extra lane work at 256; its differing control-row timings are
an important warning against treating these individual deltas as exact costs.

One bounded narrower native confirmation (3 x 0.5s) compared baseline, variant 6
and public default. Its final results, alongside the AVX2-only run above:

| Profile / weak code | Baseline ms | Final public ms | Baseline MiB/s | Final MiB/s | Throughput gain | Baseline / final CPU CV |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Native 256/254 | 30.099 | 22.911 | 115.376 | 151.571 | 31.4% | 1.40% / 2.56% |
| Native 175/173 | 26.663 | 17.358 | 88.708 | 136.265 | 53.6% | 3.45% / 5.76% |
| AVX2-only 256/254 | 33.967 | 29.031 | 102.235 | 119.617 | 17.0% | 1.35% / 7.90% |
| AVX2-only 175/173 | 37.731 | 21.921 | 62.686 | 107.899 | 72.1% | 4.03% / 3.20% |

Native variant 6 was 25.064 / 17.696 ms (CV 1.70% / 1.24%). Adding masks to
that combination reduced medians by 8.6% / 1.9%; the latter is smaller than
public-row variability. Mask-only AVX2 at 256 regressed 0.7%, within run noise,
while 175 improved 9.4%. Retention is based on the repeated shortened gains
and combined comparisons, not a claim of a universal standalone mask speedup.
Weak reduction and shortened lane padding repeatedly improved their target rows.
No independent new encoding throughput gain or small-width speedup is claimed.

All measured corpora had zero residual bits, message/block/validity failures,
valid-wrong blocks and pass-limit failures. Mean passes stayed
3.953125 / 3.906250; strong lines 262.218750 / 179.703125; weak lines
273.640625 / 268.156250; byte writes 2570.890625 / 1761.031250; bit toggles
2616.765625 / 1792.265625. These are corpus checks, not a no-miscorrection proof.

### Assembly and Fixed-Work Counters

Production-equivalent GCC assembly was inspected for native and AVX2-only.
The vector-prefix loop has one 32-byte input load, eight mask selections and
XOR accumulations, no helper calls and no stack spills. AVX2-only uses
`VPAND`/`VPXOR`; native GCC can fuse these into `VPTERNLOGQ` under its host ISA
flags. Source does not require that fusion. Horizontal reductions and table
products occur after the prefix loop. Mask traversal uses zero comparison,
`VPMOVMSKB` and set-bit iteration. No llvm-mca prediction was used for the
data-dependent product scheduler or as performance evidence.

A temporary `/bench` perf adapter replaced only benchmark minimum time with
`20x` and selected grouped `{cycles,instructions}`. Each case executed exactly
20 sweeps / 1280 timed blocks, with one diagnostic repetition. Installed perf
was `/usr/lib/linux-tools-6.8.0-139/perf`; events reported user-space counts and
100.00% running time, without multiplexing. Divide whole-process counts by 20:

| Weak width / variant | Process-amortized million cycles / sweep | Process-amortized million instructions / sweep |
| --- | ---: | ---: |
| 256 / baseline | 141.893 | 481.971 |
| 256 / combined | 123.257 | 392.306 |
| 175 / baseline | 113.525 | 384.678 |
| 175 / combined | 90.108 | 281.988 |

These include corpus creation, generic differential checks, reset, accounting,
startup and paused code. They are **not decoder-only counts** or isolated
hotspot savings. Diagnostic CPU times were 17.345 / 12.970 ms at 256 and
15.372 / 10.460 ms at 175, substantially different from earlier timing runs;
they are not substituted into the final table. No new cycle-sample attribution
was collected. The original hotspot percentages below describe the old binary.

### Validation and Reproduction

`/check native` and `/check avx2` passed during experimentation. Final
`CPP_JOBS=2 /check full` passed all portable, AVX2 and native CTest suites but
hit the 600-second command limit during sanitizer LCHRSUnittests, before reaching
CLI/fallback. It is **not recorded as a full pass**. One bounded sanitizer retry
with `GTEST_FILTER='ProductCode.*:WholeCodewordBatch.*'` passed (251.64s for
selected RS tests, 340.58s total including MC integration). This filter excludes
unrelated GTest suites. Separate unfiltered `/check cli` (7/7) and `/check
fallback` passed. Sanitizer reported no errors in the selected tests.

Extended product scheduler tests compare all eight optimization combinations
against generic scalar/full-validation output and every counter across gates,
caps 2/3/4/5/6/16, clean/noisy/overloaded inputs and widths 31/32/33/65/175/256.
Existing whole-codeword batch tests retain widths 1/31/32/33/65/175/256 and
transactional failure, parity repair and invalid-range checks; width 1 is not a
valid product weak code. Existing direct-R2 tests cover every weak length 4..256,
every position, selected exhaustive magnitudes, virtual locators, over-radius
inputs and independent scalar encoding parity. Repeated writes, stale-clean
invalidation, pending work at caps, parity gates, const-instance concurrency,
MC legacy fixtures, native CLI, data aggregation and plotting remain covered.

The wrapper omits the product source from strict fallback compilation, so both
additional commands passed:

```bash
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -I include -I src -mno-avx2 -mno-ssse3 -mno-gfni -mno-avx512f -mno-avx512bw -c src/reed_solomon/strong_weak_rs_product_code.cc -o /tmp/kilo/product-three-scalar.o
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -I include -I src -mavx2 -mssse3 -mno-gfni -mno-avx512f -mno-avx512bw -c src/reed_solomon/strong_weak_rs_product_code.cc -o /tmp/kilo/product-three-avx2.o
```

Commands below ran through `bash -c` with exported `CPP_JOBS=2`,
`CPP_BENCH_CPU=0`, `BENCHMARK_OUT_FORMAT=json` and a distinct `BENCHMARK_OUT`.
This matches existing command permissions without editing permission policy.

```bash
# Initial exact defaults, before edits:
/home/user/.config/kilo/scripts/cpp-bench native matrix '^LCH/Owned/StrongWeakRSProductCode/Correct/BSC005/Nstrong:256/Kstrong:224/Nweak:(256/Kweak:254|175/Kweak:173)$' 3 0.3s
# Initial isolated experiments:
/home/user/.config/kilo/scripts/cpp-bench native matrix '^LCH/Owned/StrongWeakRSProductCode/Experiment/(0|1|2|4|7)/(256|175)$' 3 0.3s
# Expanded native confirmation:
/home/user/.config/kilo/scripts/cpp-bench native matrix '^LCH/Owned/StrongWeakRSProductCode/(Experiment/(0|1|2|3|4|6|7)/(256|175)|Correct/BSC005/Nstrong:256/Kstrong:224/Nweak:(256/Kweak:254|175/Kweak:173))$' 5 0.5s
# AVX2 isolated and combined:
/home/user/.config/kilo/scripts/cpp-bench avx2 matrix '^LCH/Owned/StrongWeakRSProductCode/(Experiment/(0|1|2|4|6)/(256|175)|Correct/BSC005/Nstrong:256/Kstrong:224/Nweak:(256/Kweak:254|175/Kweak:173))$' 3 0.3s
# Bounded narrow native retry:
/home/user/.config/kilo/scripts/cpp-bench native matrix '^LCH/Owned/StrongWeakRSProductCode/(Experiment/(0|6)/(256|175)|Correct/BSC005/Nstrong:256/Kstrong:224/Nweak:(256/Kweak:254|175/Kweak:173))$' 3 0.5s
# Validation (each through bash -c, CPP_JOBS=2):
/home/user/.config/kilo/scripts/cpp-check full
export GTEST_FILTER='ProductCode.*:WholeCodewordBatch.*'
/home/user/.config/kilo/scripts/cpp-check sanitize
unset GTEST_FILTER
/home/user/.config/kilo/scripts/cpp-check cli
/home/user/.config/kilo/scripts/cpp-check fallback
```

JSON artifacts: `/tmp/kilo/product-three-{baseline,isolated-native,confirm-native,
final-native,final-avx2}.json`; fixed-work diagnostics:
`/tmp/kilo/product-three-perf-{256,175}-{0,7}.json`. Assembly:
`/tmp/kilo/product-three-{native,avx2}.s`. The temporary perf adapter strips the
wrapper's first four `stat` arguments, replaces `--benchmark_min_time=*` with
`--benchmark_min_time=20x`, and invokes installed perf with
`stat -x, -e '{cycles,instructions}' --` followed by the unchanged pinned command.
Set `CPP_PERF_BIN=/tmp/kilo/product-three-perf` and use exact
`'^LCH/Owned/StrongWeakRSProductCode/Experiment/0/256$'` (or 7 and/or 175),
`1 0.3s perf`. No temporary instrumentation was added to the library.

Formatting used clang-format 18.1.3; `git diff --check` passed. WSL emitted
brief clock-skew warnings during sanitizer linking; changed sources were visibly
compiled and the retry relinked affected outputs. Results remain local AMD/WSL
evidence, not Intel validation, multithread scaling or allocator-contention
evidence. Allocation remains in both Encode and Correct.

## Implementation Experiment Log (2026-09-10)

The original timings, profiling percentages, and experimental proof claims below
are historical, unverified inputs, not validation evidence for this change.
In particular, lane-varying weights cannot use a single fixed-factor shuffle or
GFNI affine matrix. Clean generic BDD already returns early on zero syndromes.

Milestones complete: baseline captured before implementation edits; R=2 direct
correction and encoding retained; reusable caller scratch measured and removed;
full validation passed. R=4 was not implemented or promoted. No child-agent tool
was available; no delegated R=4 derivation or reference experiment is claimed.
Public weak dimensions remain R=2, including 175/173. CLI/plot/metadata behavior
is unchanged.

Host: AMD Ryzen 7 8845HS, WSL2, CPU 0 pinned, GCC 13.3.0, Google Benchmark 1.9.4,
Release `-O3 -DNDEBUG`. Native flags: `-march=native` (AVX2/GFNI available).
AVX2-only flags: `-mavx2 -mssse3 -mno-gfni -mno-avx512f -mno-avx512bw`.
Weak reductions use scalar nibble-table lookups in both profiles; strong batched
BDD uses the existing backend dispatch. No new SIMD assumptions or field tables.

Baseline command: `CPP_JOBS=2 CPP_BENCH_CPU=0
/home/user/.config/kilo/scripts/cpp-bench native matrix
'^LCH/Owned/StrongWeakRSProductCode/Correct/BSC005/(Nstrong:256/Kstrong:224/Nweak:(256/Kweak:254|175/Kweak:173)|Single|StrongBatch|BothBatch)$'
5 0.2s`. This wrapper does not expose JSON output arguments. Terminal output is
the measurement record for this first run. Timing is CPU median per 64 blocks;
throughput counts information bytes (224*254 or 224*173), not full block bytes.

| Baseline row | CPU ms / 64 blocks | Information MiB/s | Mean passes |
| --- | ---: | ---: | ---: |
| Default 256/254 | 23.519 | 147.656 | 3.953125 |
| Default 175/173 | 21.738 | 108.808 | 3.906250 |
| Generic Single 256/254 | 71.608 | 48.495 | 3.953125 |
| Generic StrongBatch 256/254 | 28.302 | 122.699 | 3.953125 |
| Generic BothBatch 256/254 | 24.579 | 141.286 | 3.953125 |

All baseline corpus recovery/validity failures were zero. CPU timing CV ranged
from 5.0% to 12.9%; these are local noisy measurements, not universal estimates.

### Final Measurements

The final comparison interleaves the retained generic BDD rows with the direct
implementation on identical deterministic BSC(0.005) corpora: 64 blocks,
seed `0x5b5c0224`, cap 16, both gates enabled. Reset and result accounting are
outside timing. Every timed corpus is first compared against the independent
generic mother decoder, including every result field. Encoding setup is not
part of correction timing. All correction times below are per 64-block sweep.

| Profile / row | Generic CPU ms | Direct CPU ms | Generic MiB/s | Direct MiB/s |
| --- | ---: | ---: | ---: | ---: |
| Native 256/254 | 23.850 | 16.823 | 145.606 | 206.421 |
| Native 175/173 | 21.825 | 14.419 | 108.372 | 164.036 |
| AVX2-only 256/254 | 27.548 | 19.878 | 126.057 | 174.702 |
| AVX2-only 175/173 | 27.541 | 19.388 | 85.882 | 121.993 |

Native throughput gains in this final interleaved comparison are 41.8% and
51.4%. Against the initial default baselines they are 39.8% and 50.8%.
Final native correction CPU CV: direct 0.97%/0.81%, generic 5.23%/1.07%.
AVX2-only CPU CV: 0.24%-1.13%. These remain local measurements, not Intel-host
evidence or universal speedups.

| Native encoding | Baseline CPU us/block | Direct CPU us/block | Speed ratio |
| --- | ---: | ---: | ---: |
| 256/254 | 356.189 | 49.910 | 7.14 |
| 175/173 | 250.083 | 42.392 | 5.90 |

Encoding includes the candidate allocation/copy, weak parity, strong encode,
and commit. Baseline was measured before changing parity generation (5 x 0.3s);
final encoding uses 7 x 0.5s. Final encoding CPU CV was 0.67% and 8.70%.
No isolated weak-kernel timings, hardware-counter attribution, or llvm-mca
speed predictions are claimed for this orchestration change.
Production-equivalent GCC assembly (`-std=c++20 -O3 -DNDEBUG -march=native`)
was inspected in `/tmp/kilo/product-r2-native.s`: the reduction is inlined into
`WeakCandidate`, has two scalar table loads per symbol, and has no calls or
stack spills inside its loop. The table accessor is outside the loop and field
inversion is restricted to nonzero correction candidates, not final validation.

Final commands (run via the allowed `bash -c` wrapper, with `CPP_JOBS=2` and
`CPP_BENCH_CPU=0` exported):

```bash
/home/user/.config/kilo/scripts/cpp-bench native matrix '^LCH/Owned/StrongWeakRSProductCode/(Encode/.*|Correct/BSC005/(Nstrong:.*|BothBatch|Generic175))$' 7 0.5s
/home/user/.config/kilo/scripts/cpp-bench avx2 matrix '^LCH/Owned/StrongWeakRSProductCode/Correct/BSC005/(Nstrong:.*|BothBatch|Generic175)$' 5 0.3s
```

Google Benchmark accepts `BENCHMARK_OUT` and `BENCHMARK_OUT_FORMAT=json` via the
environment even though the wrapper has no output-file argument. Verified JSON
artifacts under `/tmp/kilo/`: `product-r2-stage1.json`,
`product-r2-workspace.json`, `product-encode-baseline.json`,
`product-r2-final-native.json`, and `product-r2-final-avx2.json`.
This report preserves the measurement summary independently of those local files.

All final correction corpora had zero message/block/validity failures and zero
residual bits. Default means: 3.953125 passes, 262.218750 strong lines,
273.640625 weak lines, 2570.890625 accepted byte writes, 2616.765625 bit toggles.
Shortened means: 3.906250 passes, 179.703125 strong lines, 268.156250 weak lines,
1761.031250 byte writes, 1792.265625 bit toggles. Generic and direct statistics
match; this is corpus evidence, not a guarantee against miscorrection.

### Scratch Experiment and Retained Scope

A caller-owned workspace overload was implemented and tested with reuse across
dimensions, invalid calls, reset, and separate concurrent callers. In a pinned
7 x 0.5s interleaved run, allocating/reusable throughput was 213.031/212.967
MiB/s at 256/254 and 166.289/162.981 MiB/s at 175/173. No worthwhile gain was
established, so the workspace API and its benchmark variants were removed.
There is no mutable shared scratch or hidden thread-local state. Ordinary
`code.Correct(block, options)` and `code.Encode(block, backend)` remain the APIs.

Strong batch staging is intentionally retained: old bytes are needed for exact
bit-change counters and only verified masked repairs commit. Direct shortened
correction allocates 87.5 KiB combined packed/mask storage rather than 128 KiB,
because only the actual strong-column batch remains; allocation reuse itself
is not claimed as a speedup. Weak correction no longer transposes, constructs a
padded mother row, or scans a full candidate to commit one byte. Weak exit checks
reduce two syndromes without inversion/BDD. Strong exit checks are unchanged.

The R=2 proof uses native parity positions 0/1 and data j+2. A repair with
magnitude S1 at native `(S0/S1) XOR 1` cancels both parity checks exactly;
out-of-mother and omitted-data positions are rejected. At N=4,K=2 the low-rate
family has the same degree-one evaluation code: translating points by 2
preserves it. Scalar LCH encoding and generic mother decoding test this boundary
independently of the new reduction helper. Over-radius inputs can still have a
valid BDD candidate; the direct decoder preserves that behavior, not an
impossible guarantee to detect all multiple errors.

### Validation

`CPP_JOBS=2 /home/user/.config/kilo/scripts/cpp-check full` passed portable,
AVX2-only, native, ASan/UBSan, CLI, and wrapper fallback checks. The wrapper omits
the product source from fallback compilation, so it was additionally compiled
with strict scalar and AVX2-no-GFNI flags, `-O2 -Wall -Wextra -Wpedantic -Werror`.
Formatting: clang-format 18.1.3. CMake emitted occasional sub-second clock-skew
warnings on WSL; changed sources were visibly compiled and tests passed.

New coverage: all weak lengths 4..256; every error position; all 255 magnitudes
at N=4,175,256; all 256 synthetic native locators at every length (including
virtual/out-of-mother rejections); randomized over-radius cases; complete
encoding parity against scalar LCH at every length. Existing gate/cap/scheduling
tests now include public direct correction and compare every counter against
the generic path. Shared const-instance concurrent calls use distinct blocks.
Legacy MC known-answer, native CLI, data aggregation, and plotting suites pass.

## Post-R2 Profiling (2026-09-10)

This section measures the retained direct contiguous R=2 implementation, not the
historical generic weak decoder below. No production or benchmark source was
edited for profiling. Base revision was `63f3aed`, with the existing uncommitted
R2/175 changes. Both rows use strong RS(256,224), the public default `Correct`,
and the same 64-block BSC(0.005) corpus, seed `0x5b5c0224`, cap 16, both gates on.

**Result:** identifiable field-algorithm work accounts for about **54.3% / 55.1%**
of resolved decoder-only cycle samples (weak lengths 256 / 175). This includes
the instructions, loads, and loop/schedule control of those algorithms, **not
just GF multiply/XOR instructions**. Another **16.6% / 19.3%** remains mixed
inside the BDD routines. The largest individually isolated non-field region is
the strong mask/transactional-writeback scan, **22.7% / 18.6%**. Neither "almost
everything is field arithmetic" nor the old "mostly orchestration" claim is
supported. This is software-region attribution, not proof of an ALU-bound or
memory-bandwidth-bound processor bottleneck.

### Environment and Timing

AMD Ryzen 7 8845HS, 8 cores/16 logical CPUs, CPU 0 affinity; WSL2 kernel
`6.6.87.2-microsoft-standard-WSL2`; GCC 13.3.0; Google Benchmark 1.9.4.
The wrapper's existing native Release build uses
`-march=native -O3 -DNDEBUG -std=gnu++20 -fPIC -march=native`, with symbol names
and unwind information but no requested source debug information. Native GFNI
FDMA and AVX2 32-byte transform leaves were observed. Google Benchmark reports
3792.64 MHz; actual operating clocks/governor were not controlled or measured.

Read-only process inspection found no running product-code Monte Carlo worker;
`mc` was idle, and two Kilo processes were active. No process was stopped or
modified. Builds used two jobs. Windows host load is not established by the
Linux process list. Timing drift was substantial, so these measurements must
not be used to claim a regression against the earlier implementation log.

| Native row | Initial median CPU ms / sweep | Initial MiB/s | Later median CPU ms / sweep | Later MiB/s | Later CPU CV |
| --- | ---: | ---: | ---: | ---: | ---: |
| 256/254 | 78.581 | 44.192 | 34.420 | 100.890 | 3.91% |
| 175/173 | 63.845 | 37.047 | 32.358 | 73.096 | 4.08% |

Each unprofiled run used three randomly interleaved repetitions, minimum 0.3s.
Initial CPU CV was 6.53% / 2.44%; initial sweeps per repetition were 5 / 7 and
later sweeps were 13 / 12. Throughput counts `224 * Kweak` information bytes per
block. Under profiling, the final repetitions ran 555 / 599 sweeps (35,520 /
38,336 blocks), with CPU times 30.650 / 19.950 ms per sweep. Those profiled times
are diagnostic only, not a speed comparison with the unprofiled runs.

All runs passed the benchmark's independent generic-decoder differential check,
including every result field. Both corpora had zero message, block, validity,
valid-wrong-block and pass-limit failures, and zero residual bits. Mean passes
were 3.953125 / 3.906250; strong lines 262.218750 / 179.703125; weak lines
273.640625 / 268.156250; accepted byte writes 2570.890625 / 1761.031250; accepted
bit toggles 2616.765625 / 1792.265625. Maximum passes were four.

### Sampling Denominator

The requested `/usr/lib/linux-tools-6.8.0-138/perf` did not exist. The installed
`/usr/lib/linux-tools-6.8.0-139/perf` (reports version 6.8.12) successfully opened
hardware cycles/instructions. The mismatched `/usr/bin/perf` wrapper was not
used. Sampling used `cycles:u`, 499 Hz and DWARF callchains with a 65,528-byte
stack capture because the batched scratch frame alone is about 31 KiB.

Whole-process profiling includes paused benchmark code. To exclude it, the
analysis selected only stacks whose benchmark ancestor is the **timed public
correction callsite**: `.constprop.0+0x29ef` for 256, `.constprop.1+0x272f` for
175. `objdump -d -C` confirms the corresponding call/return sites at
`0x12ddb/0x12de0` and `0x15e1b/0x15e20`, inside the 64-result sweep. These offsets
are specific to this ELF. Clean checks and generic/direct setup comparisons
have different return sites and are excluded even though they call `Correct`.

Calibration is removed by also rejecting samples at or before the last sampled
setup/differential correction call in the process: perf timestamps
52292.539350 / 52335.464989. Those setup stacks precede the final repetition's
timed loop; the callsite filter excludes any setup remaining after that sampled
boundary. No phase timing or source instrumentation was necessary.

Every selected sample contributes its event **period once, to its leaf only**.
Percentages are sums of periods, not sums of inclusive callgraph percentages
and not raw sample-count percentages. Known inline instruction regions are
split below; other mixed functions are left mixed.

| Capture accounting | 256/254 | 175/173 |
| --- | ---: | ---: |
| Whole-process samples | 10,809 | 9,059 |
| Final-repetition resolved decoder samples | 8,144 | 5,856 |
| Decoder event-period denominator | 44,428,472,344 | 39,834,648,551 |
| Whole-process event-period total | 62,048,114,875 | 62,105,304,175 |
| Before final setup boundary, % whole process | 20.12% | 18.71% |
| Selected decoder, % whole process | 71.60% | 64.14% |
| Post-boundary outside selected decoder, % whole process | 8.27% | 17.15% |
| Post-boundary GF stacks lacking benchmark ancestor, % whole process | 1.85% | 1.31% |

The last row is a **subset** of the outside-decoder row, not an additional cost.
After removing the pre-boundary portion, selected decoder work is 89.64% /
78.91%; the remaining 10.36% / 21.09% includes reset, result accounting, other
harness work and unresolved ancestry. It is not all measured harness cost.
Incomplete GF ancestry is 2.32% / 1.61% of that post-boundary total and can bias
the normalized decoder distribution. Selected unknown leaf weight is zero /
about 0.02%. The default capture reports one lost sample in `perf report` and
one lost chunk warning in `perf record`; shortened reports no loss. Sampling
skid, stack truncation and roughly 519 / 432 MB of trace writes are additional
limitations. Kernel cycles, including kernel allocation/page-fault work, are
not sampled by `cycles:u`. No confidence interval or exact phase split is claimed.

### Decoder-Only Self Attribution

| Disjoint leaf/instruction category | 256/254 | 175/173 |
| --- | ---: | ---: |
| Direct weak syndrome/table reduction loop | 20.14% | 15.50% |
| FFT/IFFT kernels and their scheduling/skew lookup | 20.81% | 17.95% |
| FDMA update helpers, magnitudes/division, basis conversion, field/XOR kernels | 13.34% | 21.64% |
| BDD self: inseparable arithmetic, locator/root logic, control and memory | 16.56% | 19.27% |
| Strong mask scan, bit counting and transactional writeback | 22.70% | 18.58% |
| Explicit native-codeword copying and mask/result publication | 2.50% | 1.91% |
| Remaining product control, gathers and final-validation bookkeeping | 1.97% | 2.47% |
| Weak candidate remainder: checks, locator mapping, call/return | 0.97% | 1.35% |
| Visible allocation/free/zeroing symbols | 0.19% | 0.14% |
| Other dispatch, argument checks, table access and unknown | 0.83% | 1.20% |

The first three rows give the 54.3% / 55.1% field-algorithm subtotal. It is not
a strict "arithmetic instructions" percentage: table loads, loop branches,
FFT scheduling, and local movement are necessary parts of those regions.
Conversely, inlined FDMA, differentiation and syndrome/root checking inside
`CorrectBatchImpl` are **not** assigned wholesale to orchestration or added to
that subtotal. `CorrectBatchImpl` self is 14.76% / 12.16%, and `CorrectOneImpl`
self is 1.80% / 7.11%; both remain in the mixed BDD row. Inlined zeroing remains
in its containing mixed category, so the allocator-symbol row is not an upper
bound on all buffer initialization costs.

`perf annotate` plus the matching assembly identifies `WeakCandidate+0x50..0x93`
as the inlined `WeakDataSyndromes` loop (`strong_weak_rs_product_code.cc:18-27`):
two scalar nibble-table accesses per symbol, XOR accumulation, and loop control.
It has no helper calls in the loop. The whole `WeakCandidate` self share is
21.11% / 16.84%; it must not be labeled mere candidate/gate overhead. This ELF
does not inline that reduction into `CorrectImpl`; the latter calls the helper.

`CorrectImpl+0x14a0..0x14f8` is the nested position-major mask/writeback loop
(`strong_weak_rs_product_code.cc:203-215`), including `cmpb`/conditional branch,
loads of old/new bytes, popcount, stores, and clean/active flags. In particular,
the hot `+0x14b5` branch follows the mask test; it is not GF arithmetic.
The rest of `CorrectImpl` is only 1.97% / 2.47%. These ranges are classified as
whole loops, not as exact instruction-latency or branch-misprediction costs.

Final validity is included, not subtracted: weak checks contribute to the weak
syndrome category and strong checks to BDD/transform categories; their caller
bookkeeping is mixed with the remaining product control. This run does not
provide a separate final-validity, gate, or allocation-zeroing phase percentage.

### Evidence-Based Priorities

1. The strong mask/writeback scan and direct weak syndrome loop are the largest
   isolated loops. Investigate both before assuming another generic field
   multiply optimization will dominate. Any writeback experiment must preserve
   transactional repairs, exact old/new bit counts and clean/active flags.
2. Strong arithmetic remains material. Shortening to 175 leaves five complete
   32-lane batches plus 15 scalar columns (`batched.cc:1053` and its scalar-tail
   path). The larger scalar `CorrectOneImpl` and `UpdateSamplesAVX2` shares are
   consistent with that execution path, not evidence of extra weak padding.
3. Allocation reuse is not the demonstrated single-core bottleneck. Visible
   allocator work is small, consistent with the earlier reuse experiment's lack
   of improvement. Neither experiment establishes multithread allocator
   contention; that needs its own concurrent measurement.

No optimization was implemented. AVX2-only profiling was omitted to keep this
study bounded, given the native results and unstable cross-run timings. No
llvm-mca throughput estimate substitutes for these samples.

### Reproduction and Artifacts

Commands ran through `bash -c` using the local `/bench` wrapper. Set
`CPP_JOBS=2 CPP_BENCH_CPU=0` and `BENCHMARK_OUT_FORMAT=json` in its environment.
Unprofiled initial/later outputs were `/tmp/kilo/product-profile-baseline.json`
and `/tmp/kilo/product-profile-confirm.json`:

```bash
/home/user/.config/kilo/scripts/cpp-bench native matrix '^LCH/Owned/StrongWeakRSProductCode/Correct/BSC005/Nstrong:256/Kstrong:224/Nweak:(256/Kweak:254|175/Kweak:173)$' 3 0.3s
```

The wrapper's `perf` mode validates exactly one matched row but ordinarily runs
`perf stat`. Its supported `CPP_PERF_BIN` override pointed to the temporary
`/tmp/kilo/product-perf-record` adapter, whose full body after the shebang was:

```bash
set -euo pipefail
# Remove cpp-bench's stat -x, -e EVENTS arguments, retain taskset and benchmark.
shift 4
exec /usr/lib/linux-tools-6.8.0-139/perf record -e cycles:u -F 499 --call-graph dwarf,65528 -o "${PRODUCT_PERF_DATA}" -- "$@"
```

For each dimension, set `PRODUCT_PERF_DATA=/tmp/kilo/product-profile-256.data`
or `product-profile-175.data`, and `BENCHMARK_OUT` to the corresponding `.json`:

```bash
/home/user/.config/kilo/scripts/cpp-bench native matrix '^LCH/Owned/StrongWeakRSProductCode/Correct/BSC005/Nstrong:256/Kstrong:224/Nweak:256/Kweak:254$' 1 8s perf
/home/user/.config/kilo/scripts/cpp-bench native matrix '^LCH/Owned/StrongWeakRSProductCode/Correct/BSC005/Nstrong:256/Kstrong:224/Nweak:175/Kweak:173$' 1 8s perf
/usr/lib/linux-tools-6.8.0-139/perf script -i /tmp/kilo/product-profile-256.data --demangle
/usr/lib/linux-tools-6.8.0-139/perf report -i /tmp/kilo/product-profile-256.data --stdio --no-children --percent-limit 1
```

Raw stacks, assembly, exact-symbol `perf annotate --stdio -s` output and the
period-weighted postprocessor remain under `/tmp/kilo/product-profile-*`:
`product-profile-analyze.py`, `product-profile-analysis.txt`,
`product-profile.asm`, and `product-profile-{weak,correct}.annotate`.
Annotation is whole-capture context only; the table is recomputed from filtered
stacks, not copied from annotate's whole-capture percentages. Full tests were
not rerun because no code was changed; both measured rows performed their
differential smoke validation. Existing source and unrelated dirty files were
preserved.

## Historical Analysis (Unverified Measurements)

The following original proposal is preserved as historical context. Its timing,
profiling, and projected speed claims are not results of the experiment above;
its descriptions of the "current" implementation predate the retained changes.

## 1. Executive Summary

This report analyzes the performance bottlenecks of the systematic Cantor Reed-Solomon product code implementation (`StrongWeakRSProductCode` in `src/reed_solomon/strong_weak_rs_product_code.cc` and `include/reed_solomon/strong_weak_rs_product_code.h`) and details closed-form mathematical accelerations for both decoding and encoding.

### Baseline Benchmark Measurements

Measurements taken on an AMD Ryzen 7 8845HS (x86_64, AVX2 + GFNI, Linux/WSL2) using the pinned benchmark harness (`benchmarks/strong_weak_rs_product_code_benchmarks.cc`):

| Configuration | Throughput | Mean Time / Block (64 KiB) | Directional Passes |
| :--- | :--- | :--- | :--- |
| **Default (`BothBatch`)** | **141.7 MiB/s** | **~383 $\mu$s** | 3.95 |
| **`StrongBatch` (Pass 0 batched)** | **110.8 MiB/s** | **~490 $\mu$s** | 3.95 |
| **`Single` (Scalar reference)** | **48.1 MiB/s** | **~1,128 $\mu$s** | 3.95 |

### Hardware Profiling & Hotspots

A sampling profile (`perf record`) of the default product correction benchmark revealed that CPU cycles are heavily concentrated in the orchestration layer rather than the low-level Galois field arithmetic:

- **30.82%** in `StrongWeakRSProductCode::CorrectImpl` (orchestration, transpositions, gathers, and gate checks).
- **11.14%** in `CorrectBatchImpl` (batched column/row BDD).
- **10.76%** in test benchmark harness and noise generation.
- **4.38%** in `memmove` (buffer copying and transpositions).
- **4.27%** in `CorrectOneImpl` (scalar BDD for active passes and validation).
- **4.07%** in `DivideRows` (evaluator polynomial division).
- Remaining cycles (~34%) across FFT/IFFT kernels, sample updates, and basis conversions.

---

## 2. Identified Bottlenecks

### 2.1 Heavyweight Generic BDD for Weak Single-Error Correction ($R = 2$)
The product code mandates that the weak component code has exactly $R = 2$ recovery symbols ($weak\_n - weak\_k == 2$). The bounded distance decoding (BDD) radius is:
$$t = \lfloor R / 2 \rfloor = 1$$
Every weak row can correct at most a single symbol error. Currently, weak correction executes the generic Cantor RS decoding pipeline:
1. In Pass 1 (batched): Transposes all 256 rows into a 64 KiB temporary buffer, executes 8 chunks of 32-lane SIMD BDD (evaluating locator polynomials across all 256 positions, performing polynomial root searches, evaluator reconstruction, and verification IFFTs), then performs strided gathers of 256 bytes per row before evaluating gates.
2. In Pass 3, scalar fallback, and exit validation: Calls generic `CorrectCodeword` / `CorrectOneImpl`, running scalar FDMA discrepancy loops, degree checks, root searches, and evaluator recovery on every visited row.

### 2.2 Dynamic Memory Allocations Per Block
In `CorrectImpl`:
```cpp
std::vector<Element> packed(batch_passes != 0 ? strong_n_ * mother_n : 0);
std::vector<uint8_t> masks(packed.size());
```
On every single call to `Correct` on a 64 KiB block, 128 KiB of heap memory (64 KiB `packed` + 64 KiB `masks`) is allocated, zeroed, and freed. Across thousands of Monte Carlo trials or benchmark sweeps, this causes memory allocator lock contention and cache churn.

### 2.3 Unconditional Strided Gathers in Pass 1
In Pass 1, after `CorrectCodewordBatch` finishes:
```cpp
if (batched) {
  for (size_t pos = 0; pos < decoder_length; ++pos) {
    candidate[pos] = shards[pos][line];
  }
}
```
This strided 256-byte gather is executed unconditionally for all 256 lines, even though >90% of weak rows have zero errors following the strong pass.

### 2.4 Unvectorized Scalar Weak Encoding
In `StrongWeakRSProductCode::Encode`:
```cpp
for (size_t row = 0; row < strong_k_; ++row) {
  ...
  weak_encoder_.Encode(..., 1, workspace, backend);
}
```
`weak_encoder_.Encode` is called 224 times with `byte_count = 1`. Each invocation incurs span construction, parameter checking, and scalar transform loops.

### 2.5 Redundant BDD Invocations During Exit Validation
`tracked_validation` checks unverified rows and columns at exit by calling `CorrectCodeword`. `CorrectCodeword` runs full BDD even when the purpose is solely to verify whether syndromes are all zero.

---

## 3. Mathematical Foundations: Closed-Form Weak Code ($R = 2$)

### 3.1 Cantor Basis Syndrome Properties
In the Lin-Chung-Han additive FFT framework over $GF(2^8)$, the skew factor at level 0 is:
$$\text{Skew}(0, 2b) = 2b$$
For an $R = 2$ code, the transform size is 2, and the syndrome vector $(S_0, S_1)$ is obtained by a single Radix-2 IFFT stage over pairs of elements. In native coordinates, the syndrome equations for a vector $c_{native}$ of length $N_{mother} \le 256$ reduce to:

$$S_1 = \bigoplus_{i=0}^{N_{mother}-1} c_{native}[i]$$

$$S_0 = \bigoplus_{i=0}^{N_{mother}-1} c_{native}[i] \cdot (i \oplus 1)$$

### 3.2 High-Rate Native Coordinate Mapping
For high-rate codes ($K \ge R$, such as weak code $K = 254, R = 2$):
- Native positions $0$ and $1$ are recovery symbols (public indices $K$ and $K+1$).
- Native positions $2 \dots N_{mother}-1$ are data symbols (public indices $0 \dots K-1$, with omitted zero padding in $[K, N_{mother}-2)$ for shortened codes).

Weights $w[i] = i \oplus 1$ for native positions:
- Native 0 (recovery 0): weight $= 0 \oplus 1 = 1$.
- Native 1 (recovery 1): weight $= 1 \oplus 1 = 0$.
- Native $k+2$ (data $k$): weight $= (k + 2) \oplus 1$.

### 3.3 Closed-Form Single-Error Locator & Magnitude
If an unknown error of magnitude $e \in GF(256) \setminus \{0\}$ occurs at native position $p$:
$$S_1 = e$$
$$S_0 = e \cdot (p \oplus 1)$$

Consequently:
1. **Zero Errors**: $S_0 = 0$ and $S_1 = 0 \implies$ codeword is valid.
2. **Uncorrectable ($\ge 2$ errors)**: $S_1 = 0$ and $S_0 \neq 0 \implies$ cannot be a single error; immediately uncorrectable.
3. **Single Error Candidate**: When $S_1 \neq 0$:
   $$e = S_1$$
   $$Q = S_0 / S_1 = \text{MultiplyCantor}(S_0, \text{InvertCantor}(S_1))$$
   $$p_{native} = Q \oplus 1$$

Converting to public coordinates:
- If $p_{native} < 2$: error is in recovery symbol $p_{public} = K + p_{native}$.
- If $p_{native} \ge 2$: error is in data symbol $p_{mother} = p_{native} - 2$.
  - If $p_{mother} \ge K$ and $p_{mother} < N_{mother}-2$: error falls into the shortened zero-padding subspace; rejected as invalid.
  - Otherwise: $p_{public} = p_{mother}$.

### 3.4 Experimental Proof of Equivalence
The closed-form locator and magnitude formula was verified against `LCHDecoder(254, 2)`:
- Tested against **all 65,280 single-error configurations** (256 positions $\times$ 255 non-zero magnitudes): **100% exact match**.
- Tested against shortened codes ($N = 175, K = 173$): **100% exact match**.
- Tested against random multi-error patterns and over-radius cases: produces the exact same BDD candidates and uncorrectable outcomes as `CorrectCodeword`.

### 3.5 Closed-Form Weak Encoding ($R = 2$)
By enforcing $S_0 = 0$ and $S_1 = 0$ on clean codewords:
$$S_0 = rec_0 \cdot 1 \oplus rec_1 \cdot 0 \oplus \bigoplus_{k=0}^{K-1} data[k] \cdot ((k+2) \oplus 1) = 0$$
$$S_1 = rec_0 \oplus rec_1 \oplus \bigoplus_{k=0}^{K-1} data[k] = 0$$

Direct parity calculation requires no matrix inversion or transform workspace:
$$rec_0 = \bigoplus_{k=0}^{K-1} data[k] \cdot ((k+2) \oplus 1)$$
$$rec_1 = rec_0 \oplus \bigoplus_{k=0}^{K-1} data[k]$$

Verified against `LCHEncoder(254, 2)`: **100% exact bitwise match**.

---

## 4. Acceleration Architecture & Opportunities

```
+-----------------------------------------------------------------------------+
|                               CURRENT PIPELINE                              |
|  Pass 0: Strong Batch (packed copy + CorrectCodewordBatch)                  |
|  Pass 1: Weak Batch (transpose 64KiB + batched BDD + strided gather 256x256)|
|  Pass 2: Strong Scalar (visited active columns)                             |
|  Pass 3: Weak Scalar (generic CorrectCodeword with polynomial BDD)          |
|  Exit:   Tracked validation via full CorrectCodeword                        |
+-----------------------------------------------------------------------------+
                                      |
                                      v
+-----------------------------------------------------------------------------+
|                             ACCELERATED PIPELINE                            |
|  Pass 0: In-Place Strong Batch (direct block pointers, zero malloc)         |
|  Pass 1: Direct Contiguous Weak Pass (closed-form S0/S1, zero transpose)    |
|  Pass 2: Strong Scalar (visited active columns)                             |
|  Pass 3: Direct Contiguous Weak Pass (active rows only, <1us)               |
|  Exit:   Syndrome-only validation (S0==0 && S1==0 for weak, S==0 for strong)|
+-----------------------------------------------------------------------------+
```

### Opportunity 1: In-Place Contiguous Weak Row Correction
- **Mechanism**: Weak rows are already contiguous in row-major memory `block[row * weak_n + col]`. Evaluating $S_1$ is an 8-instruction horizontal XOR reduction using AVX2 `_mm256_xor_si256`. Evaluating $S_0$ is a vector product against constant weight vector $w[i] = (i_{native} \oplus 1)$ using AVX2 `vpshufb` or GFNI affine.
- **Speedup**: Microbenchmarks show that evaluating all 256 weak rows takes ~210 $\mu$s in pure scalar and under 5 $\mu$s with AVX2 SIMD reduction, compared to ~100-150 $\mu$s for the transposed batched BDD. In Pass 3, only active lines (~17 rows) are scanned, taking under 350 ns.

### Opportunity 2: Zero-Allocation Execution
- **Mechanism**: Strong columns are already in position-major format in row-major memory: `shards[pos] = block.data() + pos * weak_n`. `shards` can point directly to `block`.
- A small reusable thread-local or stack-allocated scratch space replaces the 128 KiB heap allocations (`packed` and `masks`) entirely.

### Opportunity 3: Closed-Form Weak Encoder
- **Mechanism**: Replace 224 individual `weak_encoder_.Encode(..., 1)` calls with direct computation of $rec_0$ and $rec_1$ across each contiguous row in `block`.
- **Speedup**: Eliminates row-by-row function overhead and temporary workspace allocations, reducing weak encoding time from ~270 $\mu$s to <20 $\mu$s per block.

### Opportunity 4: Fast Syndrome-Only Exit Validation
- **Mechanism**: In `tracked_validation`, unverified lines only need syndrome checks rather than error correction.
  - Weak rows: check if $S_1 == 0$ and $S_0 == 0$ (~10 ns per row).
  - Strong columns: compute the 32 syndrome values and check if all are zero, skipping locator polynomials and root searches when clean.

---

## 5. Architectural Recommendations & Compatibility

### Implementation Placement Options
1. **Localized within `StrongWeakRSProductCode`**:
   - Implement the direct $R = 2$ weak syndrome/locator formulas as private helper methods in `strong_weak_rs_product_code.cc`.
   - Leaves public and core RS APIs unchanged.
   - Immediately accelerates product code correction and encoding.
2. **Generalizing into `CorrectOneImpl`**:
   - Add an $R = 2$ fast path branch inside `src/reed_solomon/error_correction/scalar.cc`.
   - All external callers of `CorrectCodeword` with $R = 2$ benefit from closed-form acceleration.

### Contract & Differential Testing Invariants
- `tests/product_code_tests.cc` (`InitialBatchChoicesMatchSingleOutputsAndAllCounters`) enforces that `batches = 0, 1, 2` produce identical metrics, line counts, and termination reasons.
- Implementing the closed form in the scalar weak path guarantees that both single and batched modes produce bit-for-bit identical results and preserve all counter semantics while delivering a substantial speedup.
