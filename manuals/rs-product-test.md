# Product code experiments

## Quick start

We provide two CLI tools for RS product code experiments: `rs-product-test` and
`rs-product-monte-carlo`. Build them with the `experimental` preset:

```bash
cmake --preset experimental
cmake --build --preset experimental
```

On macOS, install Xcode Command Line Tools and CMake, then use the same preset
commands. Apple Silicon uses the existing scalar field and transform kernels;
Intel Macs can use the available x86 SIMD paths. CMake supplies Clang's required
constexpr step budget automatically. No intrinsics shim or manual compiler
flags are needed. The C++ Monte-Carlo runner does not require Python.

A short end-to-end check of sampling, report files, and saved-flip replay is:

```bash
python3 -B tests/product_monte_carlo_smoke_test.py \
  build/experimental-preset/rs-product-monte-carlo
```

Plotting the reports requires Python with NumPy and Matplotlib. The plotting
script uses a headless backend, so it does not require a GUI session.

## Testing stall patterns and miscorrections

`rs-product-test sample` plants a verified error scenario, adds exactly `k` uniformly
selected distinct bit flips, and saves complete error samples independently of
testing:

For $RS(256,224)\times RS(175,173)$, each block contains
$256\times175\times8=358{,}400$ bits. A bit-error rate of
$5\times10^{-3}$ corresponds to exactly $1,792$ random flips per block:

```bash
build/experimental-preset/rs-product-test sample stall \
  --strong-n 256 --strong-k 224 --weak-n 175 --weak-k 173 \
  --random-bits 1792 --samples 100 --seed 42 > stall.json
```

Replace `stall` with another
scenario below, or `none` for background noise alone. For scenarios with
protected miscorrection columns, the same count is distributed over fewer
eligible bits, giving a slightly higher background rate in those columns that
remain eligible.

Predefined error patterns are:

- `none` meaning no special error patterns
- `stall` meaning $\left(\frac{n_1-k_1}{2}+1\right)\times \left(\frac{n_2-k_2}{2}+1\right)$ submatrix full of errors formed by random rows/columns
- `miscorrection` meaning a single miscorrection of strong component
- `mixed` meaning a single miscorrection plus errors on another strong component aligned with the miscorrection witness positions
- `miscorrection-stall` meaning miscorrection at one strong component and a stall pattern elsewhere
- `two-miscorrections` meaning two miscorrections of strong component

Use `--strong-n`, `--strong-k`, `--weak-n`, and `--weak-k` to change
dimensions, by default $RS(256, 224)\times RS(175, 173)$ is used. Planted scenarios require weak redundancy 2 and strong redundancy
at least 4. Selected columns include parity columns and are randomized per sample.

`stall.json` can then be inspected. The tool uses an all-zero transmitted
codeword; each sample's `error_hex` contains the erroneous block as a hex string
of row-major Cantor-coordinate bytes, including parity positions.

Use `rs-product-test test` to decode the erroneous blocks with postprocessing:

```bash
build/experimental-preset/rs-product-test test stall.json > results.json
```

Start with `--random-bits 0` to inspect the planted mechanisms, then sweep the
noise count and seed. Random flips exclude planted strong-miscorrection columns
(one column for `miscorrection`, `mixed`, and `miscorrection-stall`; two for
`two-miscorrections`). Their initial strong-decoder miscorrections are therefore
preserved.

In `result.json` we're interested in outcomes:

- `recovered` means that the the sample was correctly recovered
- `detected-failure` means that we couldn't recover the sample and detected failure
- `undetected-failure` means that we couldn't recover the sample but detected it as a recovery

## Monte-Carlo simulation for post-decoding BER curves

`rs-product-monte-carlo` combines
random sampling and decoding for larger experiments. It supports batches,
multiple worker threads, progress reporting, and saved statistics. It does not
emit the editable sample JSON used by `rs-product-test`.

The following finite run tests 10,000 random blocks with exactly 1,000 flipped
bits each. Dimensions are explicit to match the examples above: the Monte Carlo
runner otherwise defaults to weak $RS(256,254)$, rather than $RS(175,173)$.

```bash
build/experimental-preset/rs-product-monte-carlo \
  --output mc-random-42 --n1 256 --k1 224 --n2 175 --k2 173 \
  --minimum-flipped-bits 1000 --maximum-flipped-bits 1000 \
  --batch-size 1000 --batches 10 --threads 4 --seed 42 \
  --postprocessing --sampler fisher-yates
```

The output directory must not already exist. Postprocessing is disabled by
default. Set different minimum and
maximum flip counts to sample one uniformly chosen count per batch. Every block
in that batch receives exactly that many random flips. Omitting `--batches`
starts an unlimited run; interrupting with Ctrl-C drains in-flight work and
saves a summary. Runs cannot be resumed.

### Comparing decoder options

These three runs use the same code dimensions, seed and sampling settings,
varying only the decoder options. Each tests 100,000 blocks, with a uniformly
selected count of 1,500–2,100 flips per batch (input BER approximately
$4.19\times10^{-3}$–$5.86\times10^{-3}$). Use fresh output directories for
each run. Set both flip-count bounds to `1792` for a single input-BER point
at $5\times10^{-3}$ instead of a range.

All three mechanisms enabled: postprocessing, anchors and binary-image gating;
up to 16 half-iterations:

```bash
build/experimental-preset/rs-product-monte-carlo \
  --output mc-all-on-42 --n1 256 --k1 224 --n2 175 --k2 173 \
  --minimum-flipped-bits 1500 --maximum-flipped-bits 2100 \
  --batch-size 1000 --batches 100 --threads 4 --seed 42 \
  --sampler fisher-yates --max-directional-passes 16 \
  --postprocessing --anchors --binary-image
```

All three mechanisms disabled; up to 16 half-iterations. Omitting
`--postprocessing` disables it, while the other two mechanisms require explicit
`--no-` flags:

```bash
build/experimental-preset/rs-product-monte-carlo \
  --output mc-all-off-42 --n1 256 --k1 224 --n2 175 --k2 173 \
  --minimum-flipped-bits 1500 --maximum-flipped-bits 2100 \
  --batch-size 1000 --batches 100 --threads 4 --seed 42 \
  --sampler fisher-yates --max-directional-passes 16 \
  --no-anchors --no-binary-image
```

All three mechanisms disabled, capped at four half-iterations:

```bash
build/experimental-preset/rs-product-monte-carlo \
  --output mc-all-off-four-passes-42 --n1 256 --k1 224 --n2 175 --k2 173 \
  --minimum-flipped-bits 1500 --maximum-flipped-bits 2100 \
  --batch-size 1000 --batches 100 --threads 4 --seed 42 \
  --sampler fisher-yates --max-directional-passes 4 \
  --no-anchors --no-binary-image
```

One directional pass is one half-iteration. Four passes allow strong, weak,
strong, then weak decoding (two full iterations). The decoder may stop early
when a pass makes no changes after the initial strong pass; the current CLI
cannot force exactly four passes. Postprocessing, when enabled, is a separate
final stage and does not consume a directional pass. These switches control
decoding rules, not SIMD acceleration.

### Plotting Monte-Carlo reports and the error-floor reference

To simulate ideal bounded-distance decoding without RS arithmetic or
miscorrections, add `--ideal-decoder` to `rs-product-monte-carlo`:

```bash
build/experimental-preset/rs-product-monte-carlo --output mc-ideal \
  --n1 64 --k1 56 --n2 175 --k2 173 --ideal-decoder \
  --minimum-flipped-bits 300 --maximum-flipped-bits 450 \
  --batches 100 --batch-size 1000 --threads 4 --seed 42
```

The oracle counts a byte as erroneous if any of its bits differ from the
transmitted byte. A column with at most `floor((n1-k1)/2)` wrong bytes or a row
with at most `floor((n2-k2)/2)` wrong bytes is restored completely. It starts
with columns, then rows, and repeats until a pass makes no changes or
`--max-directional-passes` is reached, using the same pass counting as the RS
decoder. Anchors and binary-image gates do not affect the ideal decoder.

Combine `--ideal-decoder --postprocessing` to run one final ideal erasure
repair: if residual errors occupy at most `n2-k2` columns, restore them across
all rows. This oracle knows the erroneous columns; it counts one repaired stall
and zero miscorrections. It does not add directional passes or restart peeling.
In ideal mode the clean-block counter means zero residual errors, since no
syndromes are computed. Corrected bits and bytes still use their actual counts.

Metadata records `"ideal decoder": true`; report validation and saved-flip replay
retain that choice. The plotter keeps ideal and RS runs in separate curves,
even when their seeds and code dimensions match. Sampling and the conversion
from fixed-weight results to BSC contributions are the same in both modes.

```bash
python3 -B scripts/plot_product_monte_carlo.py \
  mc-all-on-42 mc-all-off-42 mc-all-off-four-passes-42 \
  --mode both --metric information --postprocessing-error-floor --output product-ber.svg
```

The plot includes one minimum-stall BSC error-floor estimate per distinct code
size, using the report dimensions. With input bit error probability $p$, byte
error probability $q=1-(1-p)^8$, and $t_i=\lfloor(n_i-k_i)/2\rfloor$, it is

$$
B_0(p)=\binom{n_1-1}{t_1}\binom{n_2-1}{t_2}
       p\,q^{(t_1+1)(t_2+1)-1}.
$$

This counts the smallest $(t_1+1)\times(t_2+1)$ bad-byte rectangles and uses
the exact conditional mean $8p/q$ wrong bits per erroneous byte. See
[Condo et al., Section III](https://arxiv.org/abs/1611.04834) for stall-pattern
counting. The reference assumes converged bounded-distance decoding without
miscorrections or postprocessing. It is an asymptotic estimate, not a finite-$p$
bound or a prediction for every decoder option. Both information and full-block
BER use the same reference by symmetry; their Monte-Carlo curves retain their
own denominators. For $RS(256,224)\times RS(175,173)$ at $p=0.005$, the estimate
is approximately $3.41\times10^{-22}$.

The default input BER range is 0.0045–0.008. Use `--p-min` and `--p-max` to
change both the axis and BSC/reference grid, and `--y-min`/`--y-max` to change
the residual BER axis. For the smaller code's floor region, for example, use
`--p-min 0.0035 --p-max 0.0055 --y-min 1e-10 --y-max 1e-2`.

The reference appears in all plotting modes. In conditional mode its horizontal
coordinate is the independent-bit channel probability $p$, while the measured
points use their fixed-weight ratio $k/(8n_1n_2)$. It uses 200 reference points
in conditional mode and the requested `--points` grid in BER/both modes.
`--no-error-floor` hides the plain BDD reference. The accompanying CSV identifies
reference rows with `record_type=error_floor` and stores `log10_error_floor`,
independently of sample coverage; `error_floor_model=bdd` identifies this model.
Observed zeros remain omitted from the plot;
they are never replaced by the theoretical estimate. The `--export` JSON
continues to contain only conditional Monte-Carlo values.

The independent `--postprocessing-error-floor` switch adds a second reference
for the failed-column erasure repair in `StrongWeakRSProductCode::CorrectImpl`.
The stage erases up to $R_\text{weak}=n_2-k_2$ syndrome-nonzero columns across
all weak rows. In the detected-column ideal model, it removes all converged
stalls confined to that many columns, including the original minimum rectangles.
The first remaining family uses $c=R_\text{weak}+1$ columns with exactly
$a=t_1+1$ bad bytes in each, and at least $b=t_2+1$ bad bytes in every affected
row. Its weight is $w=ac$. If $M_\text{PP}$ counts these supports, the estimate is

$$
B_\text{PP}(p)=\frac{w M_\text{PP}}{n_1n_2}\,p\,q^{w-1}.
$$

The script counts all leading supports, rather than just full rectangles.
For weak redundancy two, each active row contains two or three bad bytes, giving

$$
M_\text{PP}=\binom{n_2}{3}
\sum_{r=a}^{\min(n_1,\lfloor3a/2\rfloor)}
\binom{n_1}{r}\frac{r!}{(3a-2r)!\,((r-a)!)^3}.
$$

For weak redundancy four, the script counts five-column supports with row degree
at least three by grouping row degrees in the first three columns and counting
valid completions in the last two. For $RS(256,224)\times RS(175,173)$, the
leading surviving weight increases from 34 to 51 bad bytes; at $p=0.005$ this
second estimate is approximately $2.69\times10^{-24}$.

The new reference retains the converged, no-miscorrection idealization. It does
not model syndrome-clean wrong columns, value-dependent consensus, or finite
pass caps. Both references are separate leading-term estimates and may cross
outside the regime where those terms dominate; neither is a finite-$p$ bound.
The new model uses `error_floor_model=postprocessing` in the CSV and a distinct
dash-dot curve. Reference selection is independent of the decoder flags in
the reports:

| Plot references | Options |
| --- | --- |
| Plain BDD only (default) | No additional options |
| Both | `--postprocessing-error-floor` |
| Postprocessing only | `--no-error-floor --postprocessing-error-floor` |
| Neither | `--no-error-floor` |
