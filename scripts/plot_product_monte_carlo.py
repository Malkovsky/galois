#!/usr/bin/env python3
"""Plot sampled-stratum BER contributions, not extrapolated decoder BER."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

N = 524288
DENOMINATORS = {"information": 455168, "full": N}
METRICS = {"information": "residual information bits", "full": "residual full block bits"}
CODE = "RS256,224 x RS256,254 Cantor systematic row major"
RANDOM = "splitmix64 domain seeds; mt19937_64; rejection modulo; Floyd complement v1"
RANDOM_FY = "splitmix64 domain seeds; mt19937_64; rejection modulo; persistent Fisher-Yates complement v1; replay saved flips"
NEG_INF = -math.inf


def require(condition, message):
    if not condition:
        raise ValueError(message)


def natural(value, maximum=None):
    return type(value) is int and value >= 0 and (maximum is None or value <= maximum)


def strict_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f"duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_number(value):
    raise ValueError(f"noninteger JSON number: {value}")


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=strict_object,
                      parse_float=reject_number, parse_constant=reject_number)


def identity(metadata):
    text = json.dumps(metadata, sort_keys=True, ensure_ascii=True, separators=(",", ":"))
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def minimal_rows(summary, settings):
    """Validate additive schema-2 counters and project the two BER numerators."""
    rows = summary["by flipped bit count"]
    require(type(rows) is list, "invalid per-k rows")
    pooled, totals = {}, None
    for row in rows + [summary["overall"]]:
        overall = row is summary["overall"]
        require(type(row) is dict and set(row) == ({"statistics"} if overall else
                {"flipped bit count", "statistics"}), "invalid row fields")
        stats = row["statistics"]
        require(type(stats) is dict and set(stats) == {"completed blocks", "total iterations",
            "information bits", "full-codeword bits"}, "invalid minimal statistics")
        trials, iterations = stats["completed blocks"], stats["total iterations"]
        require(natural(trials) and (overall or trials > 0), "invalid completed blocks")
        require(natural(iterations) and trials * 2 <= iterations <= trials * settings[
            "maximum directional passes"], "invalid total iterations")
        flat = {"completed blocks": trials, "total iterations": iterations}
        residuals = {}
        for metric, name in (("information", "information bits"), ("full", "full-codeword bits")):
            bits = stats[name]
            require(type(bits) is dict and set(bits) == {"total bits", "raw corrupted bits",
                "post decoding corrupted bits"}, "invalid bit fields")
            total = trials * DENOMINATORS[metric]
            require(bits["total bits"] == total and natural(bits["total bits"]), "invalid total bits")
            for field in ("raw corrupted bits", "post decoding corrupted bits"):
                require(natural(bits[field], total), "invalid corrupted bits")
            flat.update({name + ":" + field: value for field, value in bits.items()})
            residuals[metric] = bits["post decoding corrupted bits"]
        for field in ("raw corrupted bits", "post decoding corrupted bits"):
            difference = stats["full-codeword bits"][field] - stats["information bits"][field]
            require(0 <= difference <= trials * (N - DENOMINATORS["information"]),
                    "inconsistent full/information bits")
        if overall:
            require(flat == totals if rows else all(v == 0 for v in flat.values()),
                    "overall/per-k reconciliation failed")
            continue
        k = row["flipped bit count"]
        require(natural(k, N) and settings["minimum flipped bits"] <= k <= settings[
            "maximum flipped bits"], "invalid flipped bit count")
        require(k not in pooled, f"duplicate k: {k}")
        require(stats["full-codeword bits"]["raw corrupted bits"] == trials*k,
                "initial count is not exact k")
        if totals is None:
            totals = dict.fromkeys(flat, 0)
        for key, value in flat.items():
            totals[key] += value
        pooled[k] = {"trials": trials, **residuals}
    return pooled


def load_report(path):
    """Validate a summary snapshot without loading the native library or journal."""
    path = Path(path)
    if path.is_dir():
        path /= "summary.json"
    try:
        metadata = read_json(path.parent / "metadata.json")
        summary = read_json(path)
        require(type(metadata) is dict and set(metadata) - {"codeword"} == {
            "schema revision", "created at", "settings", "code", "random algorithm"},
            "incompatible metadata fields")
        codeword = metadata.get("codeword", "random")
        require(codeword in ("zero", "random"), "incompatible codeword convention")
        require(type(metadata["schema revision"]) is int and metadata["schema revision"] in (1, 2)
                and metadata["code"] == CODE and metadata["random algorithm"] in (RANDOM, RANDOM_FY),
                "incompatible schema/code/random algorithm")
        require(isinstance(metadata["created at"], str), "invalid created at")
        settings = metadata["settings"]
        integer_settings = {"root seed", "batch size", "batches", "threads",
                            "minimum flipped bits", "maximum flipped bits",
                            "maximum directional passes", "checkpoint trials",
                            "report seconds", "fsync seconds"}
        require(type(settings) is dict and set(settings) == integer_settings | {
            "anchors", "binary image"}, "incompatible settings")
        for name in integer_settings:
            require(natural(settings[name], (1 << 64) - 1), f"invalid {name}")
        for name in ("anchors", "binary image"):
            require(type(settings[name]) is bool, f"invalid {name}")
        require(0 <= settings["minimum flipped bits"] <= settings["maximum flipped bits"] <= N,
                "invalid sampled k range")
        for name, low, high in (("batch size", 1, (1 << 64) - 1), ("threads", 1, 1024),
                                ("maximum directional passes", 2, 1000000),
                                ("checkpoint trials", 1, 4096), ("report seconds", 1, 86400),
                                ("fsync seconds", 1, 86400)):
            require(low <= settings[name] <= high, f"invalid {name}")
        require(type(summary) is dict and set(summary) == {
            "schema revision", "run identity", "overall", "by flipped bit count"},
            "incompatible summary fields")
        require(type(summary["schema revision"]) is int and summary["schema revision"] == metadata["schema revision"]
                and summary["run identity"] == identity(metadata), "identity/schema mismatch")
        if summary["schema revision"] == 2:
            pooled = minimal_rows(summary, settings)
            config = (settings["maximum directional passes"], settings["anchors"], settings["binary image"])
            return summary["run identity"], settings["root seed"], config, pooled, codeword
        rows = summary["by flipped bit count"]
        require(type(rows) is list, "invalid per-k rows")
        pooled = {}
        totals = {}
        count = 0
        for row in rows + [summary["overall"]]:
            overall = row is summary["overall"]
            require(type(row) is dict and set(row) == ({"trial count", "statistics"} if overall
                    else {"flipped bit count", "trial count", "statistics"}), "invalid row fields")
            trials = row["trial count"]
            require(natural(trials) and (overall or trials > 0), "invalid trial count")
            stats = row["statistics"]
            require(type(stats) is dict and all(m in stats for m in METRICS.values()),
                    "missing residual statistics")
            for name, moment in stats.items():
                require(type(moment) is dict and set(moment) == {"sum", "squared sum"},
                        f"invalid moment fields: {name}")
                total, square = moment["sum"], moment["squared sum"]
                require(natural(total) and natural(square) and total <= square
                        and total * total <= trials * square, f"invalid moments: {name}")
                require(trials > 0 or total == square == 0, "nonzero empty statistics")
                if name in METRICS.values():
                    d = DENOMINATORS[next(m for m in METRICS if METRICS[m] == name)]
                    require(total <= trials * d and square <= d * total,
                            f"residual moments exceed bit count: {name}")
            require(stats[METRICS["information"]]["sum"] <= stats[METRICS["full"]]["sum"],
                    "information residual exceeds full residual")
            if overall:
                require(trials == count and (stats == totals if rows else all(
                    item == {"sum": 0, "squared sum": 0} for item in stats.values())),
                    "overall/per-k reconciliation failed")
                continue
            k = row["flipped bit count"]
            require(natural(k, N) and settings["minimum flipped bits"] <= k <= settings[
                "maximum flipped bits"], "invalid flipped bit count")
            require(k not in pooled, f"duplicate k: {k}")
            if "initial full block corrupted bits" in stats:
                require(stats["initial full block corrupted bits"] == {
                    "sum": trials * k, "squared sum": trials * k * k}, "initial count is not exact k")
            if not totals:
                totals = {name: {"sum": 0, "squared sum": 0} for name in stats}
            require(set(stats) == set(totals), "inconsistent metric sets")
            for name in stats:
                for field in ("sum", "squared sum"):
                    totals[name][field] += stats[name][field]
            count += trials
            pooled[k] = {"trials": trials, **{m: stats[name]["sum"] for m, name in METRICS.items()}}
        config = (settings["maximum directional passes"], settings["anchors"], settings["binary image"])
        return summary["run identity"], settings["root seed"], config, pooled, codeword
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{path}: {error}") from error


def pool_reports(paths, allow_mixed_codewords=False):
    groups, identities, seeds, conventions = {}, set(), set(), {}
    for path in paths:
        digest, seed, config, rows, codeword = load_report(path)
        require(digest not in identities, f"duplicate run identity: {path}")
        require((config, seed) not in seeds,
                f"repeated seed {seed} within configuration {config}: {path}; "
                "sample overlap cannot be excluded, even with different k ranges")
        identities.add(digest)
        seeds.add((config, seed))
        previous = conventions.setdefault(config, codeword)
        require(allow_mixed_codewords or previous == codeword,
                "mixed zero/random codewords require --allow-mixed-codewords")
        group = groups.setdefault(config, {})
        for k, row in rows.items():
            target = group.setdefault(k, dict.fromkeys(row, 0))
            for name, value in row.items():
                target[name] += value
    return groups


def logsumexp(values):
    values = list(values)
    maximum = max(values, default=NEG_INF)
    if maximum == NEG_INF:
        return maximum
    return maximum + math.log(math.fsum(math.exp(v - maximum) for v in values))


def log_binomial(n, k, p):
    require(natural(n) and natural(k, n) and 0 < p < 1, "invalid binomial arguments")
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
            + k * math.log(p) + (n - k) * math.log1p(-p))


def log_interval(n, low, high, p):
    """Sum a binomial interval from its mode, bounding discarded relative tails."""
    if low > high:
        return NEG_INF
    mode = min(high, max(low, int((n + 1) * p)))
    terms = [1.0]
    for direction, end in ((-1, low), (1, high)):
        k, term = mode, 1.0
        while k != end:
            ratio = (k / (n - k + 1) * (1 - p) / p if direction < 0
                     else (n - k) / (k + 1) * p / (1 - p))
            # Ratios decrease away from the mode. The geometric bound covers
            # the entire omitted tail, not merely the next term.
            if ratio < 1 and term * ratio / (1 - ratio) < 1e-16:
                break
            term *= ratio
            terms.append(term)
            k += direction
    return log_binomial(n, mode, p) + math.log(math.fsum(terms))


def log1mexp(value):
    if value == NEG_INF:
        return 0.0
    if value == 0:
        return NEG_INF
    return (math.log1p(-math.exp(value)) if value < -math.log(2)
            else math.log(-math.expm1(value)))


def evaluate(rows, p, n=N, denominators=None):
    """Return natural-log contributions and masses; never renormalize weights."""
    denominators = DENOMINATORS if denominators is None else denominators
    weights = {k: log_binomial(n, k, p) for k in rows}
    covered = logsumexp(weights.values())
    gaps, start = [], 0
    for k in sorted(rows):
        if start < k:
            gaps.append((start, k - 1))
        start = k + 1
    if start <= n:
        gaps.append((start, n))
    missing = logsumexp(log_interval(n, low, high, p) for low, high in gaps)
    # Derive only the larger probability by subtraction. Tiny missing tails
    # must not be inferred from a rounded covered mass close to one.
    if covered <= missing:
        missing = log1mexp(min(0.0, covered))
    else:
        covered = log1mexp(min(0.0, missing))
    contributions = {metric: logsumexp(
        weights[k] + math.log(row[metric]) - math.log(row["trials"]) - math.log(d)
        for k, row in rows.items() if row[metric] > 0)
        for metric, d in denominators.items()}
    return contributions, covered, missing


def config_label(config):
    passes, anchors, binary = config
    return f"passes={passes}, anchors={'on' if anchors else 'off'}, binary-image={'on' if binary else 'off'}"


def plot_value(value):
    result = math.exp(value)
    return result if result > 0 else math.nan


def main(argv=None):
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Weights are Binomial(N,p), without renormalization. Missing k are unknown; "
        "zero observed errors are not certainty. Arbitrarily low plotted values are not reliability "
        "evidence. No extrapolation or MSE fit: a fitting model has not been specified. "
        "Matching decoder configurations pool per-k sums/counts; different flags/caps stay separate. "
        "Repeated seeds within a configuration are rejected. Metadata has no source revision.")
    parser.add_argument("inputs", nargs="+", type=Path, help="run directories or summary.json with sibling metadata.json")
    parser.add_argument("--output", type=Path, default=Path("product_monte_carlo.svg"), help="SVG output")
    parser.add_argument("--csv", type=Path, help="CSV output (default: output path with .csv suffix)")
    parser.add_argument("--points", type=int, default=200)
    parser.add_argument("--metric", choices=("information", "full", "both"), default="information")
    parser.add_argument("--allow-mixed-codewords", action="store_true",
                        help="pool zero/random inputs using linear-code BDD translation equivariance; "
                        "absent legacy convention means random; repeated seeds remain forbidden")
    args = parser.parse_args(argv)
    require(args.points >= 2, "--points must be at least 2")
    require(args.output.suffix.lower() == ".svg", "--output must be an .svg path")
    csv_path = args.csv or args.output.with_suffix(".csv")
    protected = {p.resolve() for path in args.inputs for p in (
        (path / "summary.json" if path.is_dir() else path),
        (path if path.is_dir() else path.parent) / "metadata.json")}
    require(args.output.resolve() != csv_path.resolve() and not protected.intersection(
        {args.output.resolve(), csv_path.resolve()}), "output paths must be distinct from each other and inputs")
    groups = pool_reports(args.inputs, args.allow_mixed_codewords)
    print("Warning: metadata lacks source revision; decoder implementation compatibility cannot be verified. "
          "Missing strata are unknown; zero observed errors do not establish zero BER. "
          "No extrapolation or MSE fit is performed.", file=sys.stderr)
    metrics = list(METRICS) if args.metric == "both" else [args.metric]
    ps = [0.008 + (0.0045 - 0.008) * i / (args.points - 1) for i in range(args.points)]
    results = []
    for config, rows in groups.items():
        values = [evaluate(rows, p) for p in ps]
        label = config_label(config)
        for metric in metrics:
            zeros = sum(row[metric] == 0 for row in rows.values())
            print(f"{label}; {metric}: {len(rows)}/{N + 1} sampled strata, "
                  f"{zeros} zero-observed strata; log10 missing mass range "
                  f"[{min(v[2] for v in values) / math.log(10):.6g}, "
                  f"{max(v[2] for v in values) / math.log(10):.6g}]", file=sys.stderr)
            results.append((label, metric, zeros, len(rows), values))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise ValueError("plotting requires matplotlib; numerical core uses only the standard library") from error
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["configuration", "metric", "p", "log10_contribution", "log10_covered_mass",
                         "log10_missing_mass", "sampled_strata", "zero_observed_strata"])
        for label, metric, zeros, sampled, values in results:
            for p, (contributions, covered, missing) in zip(ps, values):
                writer.writerow([label, metric, p, contributions[metric] / math.log(10),
                                 covered / math.log(10), missing / math.log(10), sampled, zeros])
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "svg.hashsalt": "gf256-product-monte-carlo"})
    figure, axis = plt.subplots(figsize=(11, 7))
    for label, metric, _, _, values in results:
        axis.plot(ps, [plot_value(v[0][metric]) for v in values], linewidth=2,
                  label=f"{metric}: {label}")
    axis.set(xlim=(0.008, 0.0045), ylim=(1e-30, 1e-1), yscale="log",
             xlabel="Channel bit-flip probability p",
             ylabel="Sampled-stratum BER contribution",
             title="Product-code sampled-stratum BER contribution")
    axis.grid(True, which="major", color="#d1d5db", alpha=0.75)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(fontsize=8)
    figure.text(0.5, 0.02, "Missing strata unknown; zero observed errors are not certainty. No extrapolation.",
                ha="center", fontsize=9)
    figure.tight_layout(rect=(0, 0.04, 1, 1))
    figure.savefig(args.output, facecolor="white", metadata={"Creator": __file__, "Date": None})
    plt.close(figure)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        sys.exit(1)
