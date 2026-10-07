#!/usr/bin/env python3
"""Pool reports and plot fixed-weight conditional BER and/or BSC contributions."""

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
SNAPSHOTS = "atomic summary v1"
NEG_INF = -math.inf
DEFAULT_DIMENSIONS = (256, 224, 256, 254)
POSTPROCESSING_METRICS = ("stall patterns corrected", "strong miscorrections detected",
                        "strong miscorrections corrected")


def dimensions(settings):
    require(type(settings) is dict, "incompatible settings")
    names = ("n1", "k1", "n2", "k2")
    require(not any(name in settings for name in names) or all(name in settings for name in names),
            "incomplete dimensions")
    dims = tuple(settings.get(name, default) for name, default in zip(names, DEFAULT_DIMENSIONS))
    require(all(natural(value, 256) for value in dims), "invalid dimensions")
    n1, k1, n2, k2 = dims
    power2 = lambda n: n > 0 and n & (n - 1) == 0
    require(power2(n1) and 2 <= n1 - k1 <= k1 and power2(n1 - k1)
            and 2 <= k2 < n2 <= 256 and (n2 - k2 == 2 or (n2, k2) == (256, 252)), "unsupported dimensions")
    return dims


def denominators(dims):
    n1, k1, n2, k2 = dims
    return {"information": 8 * k1 * k2, "full": 8 * n1 * n2}


def configuration(settings):
    flags = (settings["maximum directional passes"], settings["anchors"], settings["binary image"],
             settings.get("postprocessing", False))
    dims = dimensions(settings)
    return flags if dims == DEFAULT_DIMENSIONS else flags + dims


def config_dimensions(config):
    return config[4:] if len(config) == 8 else DEFAULT_DIMENSIONS


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
    ds = denominators(dimensions(settings))
    n = ds["full"]
    rows = summary["by flipped bit count"]
    require(type(rows) is list, "invalid per-k rows")
    pooled, totals = {}, None
    for row in rows + [summary["overall"]]:
        overall = row is summary["overall"]
        require(type(row) is dict and set(row) == ({"statistics"} if overall else
                {"flipped bit count", "statistics"}), "invalid row fields")
        stats = row["statistics"]
        fields = {"completed blocks", "total iterations", "information bits", "full-codeword bits"}
        extra = set(POSTPROCESSING_METRICS)
        require(type(stats) is dict and (set(stats) == fields | extra or
                (set(stats) == fields and not settings.get("postprocessing", False))),
                "invalid minimal statistics")
        trials, iterations = stats["completed blocks"], stats["total iterations"]
        require(natural(trials) and (overall or trials > 0), "invalid completed blocks")
        require(natural(iterations) and trials * 2 <= iterations <= trials * settings[
            "maximum directional passes"], "invalid total iterations")
        flat = {"completed blocks": trials, "total iterations": iterations}
        for name in POSTPROCESSING_METRICS:
            value = stats.get(name, 0)
            require(natural(value, (1 << 64) - 1), "invalid postprocessing counter")
            require(settings.get("postprocessing", False) or value == 0,
                    "postprocessing counters with postprocessing disabled")
            flat[name] = value
        require(flat[POSTPROCESSING_METRICS[0]] <= trials and
                flat[POSTPROCESSING_METRICS[2]] <= flat[POSTPROCESSING_METRICS[1]] <= trials * dimensions(settings)[2],
                "inconsistent postprocessing counters")
        residuals = {}
        for metric, name in (("information", "information bits"), ("full", "full-codeword bits")):
            bits = stats[name]
            require(type(bits) is dict and set(bits) == {"total bits", "raw corrupted bits",
                "post decoding corrupted bits"}, "invalid bit fields")
            total = trials * ds[metric]
            require(bits["total bits"] == total and natural(bits["total bits"]), "invalid total bits")
            for field in ("raw corrupted bits", "post decoding corrupted bits"):
                require(natural(bits[field], total), "invalid corrupted bits")
            flat.update({name + ":" + field: value for field, value in bits.items()})
            residuals[metric] = bits["post decoding corrupted bits"]
        for field in ("raw corrupted bits", "post decoding corrupted bits"):
            difference = stats["full-codeword bits"][field] - stats["information bits"][field]
            require(0 <= difference <= trials * (n - ds["information"]),
                    "inconsistent full/information bits")
        if overall:
            require(flat == totals if rows else all(v == 0 for v in flat.values()),
                    "overall/per-k reconciliation failed")
            continue
        k = row["flipped bit count"]
        require(natural(k, n) and settings["minimum flipped bits"] <= k <= settings[
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
        require(type(metadata) is dict and set(metadata) - {"codeword", "storage"} == {
            "schema revision", "created at", "settings", "code", "random algorithm"},
            "incompatible metadata fields")
        snapshot = "storage" in metadata
        require(not snapshot or (metadata["storage"] == SNAPSHOTS and metadata["schema revision"] == 2),
                "incompatible snapshot storage/schema")
        codeword = metadata.get("codeword", "random")
        require(codeword in ("zero", "random"), "incompatible codeword convention")
        require(type(metadata["schema revision"]) is int and metadata["schema revision"] in (1, 2)
                and metadata["random algorithm"] in (RANDOM, RANDOM_FY),
                "incompatible schema/code/random algorithm")
        require(isinstance(metadata["created at"], str), "invalid created at")
        settings = metadata["settings"]
        dims = dimensions(settings)
        n1, k1, n2, k2 = dims
        ds = denominators(dims)
        n = ds["full"]
        require(metadata["code"] == f"RS{n1},{k1} x RS{n2},{k2} Cantor systematic row major",
                "incompatible schema/code/random algorithm: code/dimensions/coordinates mismatch")
        integer_settings = {"root seed", "batch size", "batches", "threads",
                            "minimum flipped bits", "maximum flipped bits",
                            "maximum directional passes", "checkpoint trials",
                             "report seconds", "fsync seconds"}
        if "n1" in settings:
            integer_settings |= {"n1", "k1", "n2", "k2"}
        require(type(settings) is dict and set(settings) - {"postprocessing"} == integer_settings | {
            "anchors", "binary image"}, "incompatible settings")
        for name in integer_settings:
            require(natural(settings[name], (1 << 64) - 1), f"invalid {name}")
        for name in ("anchors", "binary image"):
            require(type(settings[name]) is bool, f"invalid {name}")
        require(type(settings.get("postprocessing", False)) is bool, "invalid postprocessing")
        require(not settings.get("postprocessing", False) or snapshot,
                "postprocessing requires schema-2 snapshot storage")
        require(0 <= settings["minimum flipped bits"] <= settings["maximum flipped bits"] <= n,
                "invalid sampled k range")
        for name, low, high in (("batch size", 1, (1 << 64) - 1), ("threads", 1, 1024),
                                ("maximum directional passes", 2, 1000000),
                                ("checkpoint trials", 1, 4096), ("report seconds", 1, 86400),
                                ("fsync seconds", 1, 86400)):
            require(low <= settings[name] <= high, f"invalid {name}")
        require(type(summary) is dict and set(summary) - {"code parameters", "checkpoint"} == {
            "schema revision", "run identity", "overall", "by flipped bit count"},
            "incompatible summary fields")
        recorded_snapshot = snapshot and metadata["random algorithm"] == RANDOM_FY
        require(("checkpoint" in summary) == recorded_snapshot, "invalid snapshot checkpoint")
        if recorded_snapshot:
            checkpoint = summary["checkpoint"]
            require(type(checkpoint) is dict and set(checkpoint) == {"flip end"}
                    and natural(checkpoint["flip end"], (1 << 64) - 1)
                    and checkpoint["flip end"] >= 40, "invalid committed flip boundary")
        if "code parameters" in summary:
            parameters = summary["code parameters"]
            require(type(parameters) is dict and set(parameters) == {"n1", "k1", "n2", "k2"}
                    and dimensions(parameters) == dims, "summary code parameters mismatch")
        require(type(summary["schema revision"]) is int and summary["schema revision"] == metadata["schema revision"]
                and summary["run identity"] == identity(metadata), "identity/schema mismatch")
        if summary["schema revision"] == 2:
            pooled = minimal_rows(summary, settings)
            config = configuration(settings)
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
                    d = ds[next(m for m in METRICS if METRICS[m] == name)]
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
            require(natural(k, n) and settings["minimum flipped bits"] <= k <= settings[
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
        config = configuration(settings)
        return summary["run identity"], settings["root seed"], config, pooled, codeword
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{path}: {error}") from error


def discover_reports(paths):
    """Expand containers in sorted order, stopping at runs and rejecting overlap."""
    reports, sources = [], {}

    def visit(path, source):
        if path.is_dir():
            metadata, summary = path / "metadata.json", path / "summary.json"
            if metadata.exists() or summary.exists():
                require(metadata.is_file() and summary.is_file(),
                        f"{path}: incomplete run; expected both metadata.json and summary.json "
                        "as files; finish the report or move the partial run outside the input tree")
                visit(summary, source)
            else:
                for child in sorted(path.iterdir()):
                    if not child.is_symlink() and child.is_dir():
                        visit(child, source)
            return
        require(path.is_file(), f"input does not exist or is not a report file: {path}")
        resolved = path.resolve()
        require(resolved not in sources,
                f"duplicate report source (duplicate run identity): {path}; "
                f"overlapping inputs {sources.get(resolved)} and {source}; supply each run only once")
        sources[resolved] = source
        reports.append(path)

    for source in paths:
        before = len(reports)
        visit(Path(source), source)
        require(len(reports) > before,
                f"{source}: no runs found; expected metadata.json + summary.json pairs "
                "in this directory or its descendants (subdirectory symlinks are not followed)")
    require(reports, "no runs found: supply run directories, containers, or summary files")
    return reports


def pool_reports(paths, allow_mixed_codewords=False):
    groups, identities, seeds, conventions = {}, set(), set(), {}
    for path in discover_reports(paths):
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
    passes, anchors, binary, postprocessing = config[:4]
    n1, k1, n2, k2 = config_dimensions(config)
    return (f"RS{n1},{k1} x RS{n2},{k2}; passes={passes}, "
            f"anchors={'on' if anchors else 'off'}, binary-image={'on' if binary else 'off'}, "
            f"postprocessing={'on' if postprocessing else 'off'}")


def plot_value(value):
    result = math.exp(value)
    return result if result > 0 else math.nan


def conditional_values(rows, metric, n=N, ds=None):
    """Yield provenance and normalized conditional BER in ascending k."""
    ds = DENOMINATORS if ds is None else ds
    for k, row in sorted(rows.items()):
        total, trials = row[metric], row["trials"]
        mean = total / trials
        log_mean = math.log(total) - math.log(trials) if total else NEG_INF
        log_ber = log_mean - math.log(ds[metric])
        yield {"k": k, "residual_bits_sum": total, "completed_blocks": trials,
               "mean": mean, "log10_mean": log_mean / math.log(10), "raw_ber": k / n,
               "conditional_ber": math.exp(log_ber),
               "log10_conditional_ber": log_ber / math.log(10)}


def export_conditional(groups, metrics, path):
    """Export pooled conditional BER points, independent of plotted mode/limits."""
    records = []
    for config, rows in groups.items():
        ds = denominators(config_dimensions(config))
        for metric in metrics:
            points = []
            for value in conditional_values(rows, metric, ds["full"], ds):
                log_ber = value["log10_conditional_ber"]
                points.append({"raw_ber": value["raw_ber"],
                               "residual_ber": value["conditional_ber"],
                               "log10_residual_ber": log_ber if math.isfinite(log_ber) else None})
            records.append({"configuration": config_label(config), "metric": metric,
                            "points": points})
    with path.open("w", encoding="utf-8") as stream:
        json.dump(records, stream, indent=2, allow_nan=False)
        stream.write("\n")


def plot_results(groups, metrics, mode, ps, output, csv_path, plt, thin=0.5):
    require(math.isfinite(thin) and thin > 0, "--thin must be finite and greater than zero")
    figure, axis = plt.subplots(figsize=(11, 7))
    positive = False
    total_underflows = 0
    sampled = set()
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=[
            "configuration", "metric", "p", "log10_contribution", "log10_covered_mass",
            "log10_missing_mass", "sampled_strata", "zero_observed_strata", "record_type",
            "k", "residual_bits_sum", "completed_blocks", "mean", "log10_mean",
            "raw_ber", "conditional_ber", "log10_conditional_ber"])
        writer.writeheader()
        for index, (config, rows) in enumerate(groups.items()):
            ds = denominators(config_dimensions(config))
            n = ds["full"]
            label = config_label(config)
            color = colors[index % len(colors)]
            sampled.update(rows)
            weighted = [evaluate(rows, p, n, ds) for p in ps] if mode != "conditional" else []
            for metric in metrics:
                zeros = sum(row[metric] == 0 for row in rows.values())
                base = {"configuration": label, "metric": metric,
                        "sampled_strata": len(rows), "zero_observed_strata": zeros}
                if mode != "ber":
                    values = list(conditional_values(rows, metric, n, ds))
                    for value in values:
                        writer.writerow({**base, "record_type": "conditional", **value})
                    ys = [plot_value(v["log10_conditional_ber"] * math.log(10)) for v in values]
                    visible = [(v["raw_ber"], y) for v, y in zip(values, ys) if math.isfinite(y)]
                    underflows = sum(v["residual_bits_sum"] > 0 and not math.isfinite(y)
                                     for v, y in zip(values, ys))
                    total_underflows += underflows
                    positive |= bool(visible)
                    clipped = sum(not (0.0045 <= x <= 0.008 and 1e-30 <= y <= 1e-1)
                                  for x, y in visible)
                    note = f"{zeros}/{len(values)} zero-observed strata omitted"
                    print(f"{label}; {metric}: {note}; {underflows} positive BERs underflowed "
                          f"(see CSV logs); {clipped} points outside axis limits", file=sys.stderr)
                    axis.scatter([x for x, _ in visible], [y for _, y in visible], s=18 * thin**2,
                                 linewidths=thin,
                                 color=color, marker="o" if metric == "information" else "x",
                                 label=f"Fixed-weight conditional BER; {metric}: {label}\n{note}")
                if mode != "conditional":
                    print(f"{label}; {metric}: {len(rows)}/{n + 1} sampled strata, "
                          f"{zeros} zero-observed strata; log10 missing mass range "
                          f"[{min(v[2] for v in weighted) / math.log(10):.6g}, "
                          f"{max(v[2] for v in weighted) / math.log(10):.6g}]", file=sys.stderr)
                    for p, (contributions, covered, missing) in zip(ps, weighted):
                        writer.writerow({**base, "record_type": "ber", "p": p,
                                         "log10_contribution": contributions[metric] / math.log(10),
                                         "log10_covered_mass": covered / math.log(10),
                                         "log10_missing_mass": missing / math.log(10)})
                    axis.plot(ps, [plot_value(v[0][metric]) for v in weighted], linewidth=2 * thin,
                              color=color, linestyle="-" if metric == "information" else "--",
                              label=f"Sampled-stratum BSC contribution; {metric}: {label}")
    axis.set(xscale="linear", yscale="log", xlim=(0.008, 0.0045), ylim=(1e-30, 1e-1),
             xlabel="Raw BER (k / transmitted bits for conditional; p for BSC)",
             ylabel="Residual BER", title="Product-code BER")
    if mode == "conditional" and not positive:
        message = ("All observed residual sums are zero; no positive BERs to plot."
                    if sampled else "No completed blocks in these report snapshots.")
        if total_underflows:
            message = "No representable positive BERs to plot; see CSV logs."
        axis.text(0.5, 0.5, message + "\nNo artificial floor is used.",
                  transform=axis.transAxes, ha="center", va="center")
    axis.grid(True, which="major", color="#d1d5db", alpha=0.75)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(fontsize=8)
    figure.text(0.5, 0.02,
                "Fixed-weight conditional BER is not a full BSC expectation. Missing strata unknown.\n"
                "Zeros omitted, not floored; zero observed errors are not certainty. "
                "No weight renormalization or extrapolation.", ha="center", fontsize=9)
    figure.tight_layout(rect=(0, 0.06, 1, 1))
    figure.savefig(output, facecolor="white", metadata={"Creator": __file__, "Date": None})
    plt.close(figure)


def main(argv=None):
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="BER mode uses Binomial(N,p) weights, without renormalization. Conditional mode uses "
        "raw BER k/N and residual BER sum(residual bits)/(sum(completed blocks)*D), "
        "N=8*n1*n2; D=8*k1*k2 for information or N for full, per code. Both mode overlays conditional scatter and "
        "sampled-stratum BSC contribution lines; conditional BER is not a full BSC expectation. "
        "All modes use descending linear x 0.008 to 0.0045 and log y 1e-30 to 1e-1; "
        "out-of-range points are clipped, not discarded from CSV. "
        "Zeros are omitted and counted, never floored. The unified CSV uses record_type conditional/ber, "
        "retains BER columns and conditional exact sums/counts and mean/log10_mean, and adds "
        "raw_ber, conditional_ber/log10_conditional_ber. Inapplicable fields are blank. "
        "Missing k are unknown; "
        "zero observed errors are not certainty. Arbitrarily low plotted values are not reliability "
        "evidence. No extrapolation or MSE fit: a fitting model has not been specified. "
        "Matching decoder configurations pool per-k sums/counts; different dimensions/flags/caps stay separate. "
        "Repeated seeds within a configuration and overlapping inputs are rejected. "
        "Cantor systematic row-major codes only: strong N,R powers of two, N<=256, 2<=R<=K; "
        "weak N<=256, K>=2, R=2, including shortening, or RS(256,252). Absent dimensions mean 256,224,256,254. "
        "Zero/random conventions must match unless --allow-mixed-codewords is explicit. "
        "Metadata has no source revision. "
        "Discovery stops at run directories, ignores unrelated files, and does not follow subdirectory "
        "symlinks. Partial or malformed runs are errors, not silently skipped. "
        "Example: python3 -B scripts/plot_product_monte_carlo.py experiments --mode conditional "
        "--metric both --output merged.svg (also writes merged.csv). "
        "Direct inputs remain supported: run-a run-b/summary.json --output merged.svg.")
    parser.add_argument("inputs", nargs="+", type=Path,
                        help="run directories, summary files with sibling metadata.json, or containers "
                        "recursively searched in sorted order for metadata.json + summary.json pairs")
    parser.add_argument("--output", type=Path, default=Path("product_monte_carlo.svg"), help="SVG output")
    parser.add_argument("--csv", type=Path, help="unified CSV output, both row types in both mode "
                        "(default: output path with .csv suffix; no sidecar)")
    parser.add_argument("--export", type=Path, metavar="JSON",
                        help="also export only pooled conditional raw/residual BER values to indented JSON, "
                        "with null log10 for zero estimates, "
                        "including zeros and out-of-axis points, regardless of --mode")
    parser.add_argument("--mode", choices=("conditional", "ber", "both"), default="ber",
                        help="conditional: fixed-weight BER scatter; ber: BSC contribution lines; both: overlay")
    parser.add_argument("--points", type=int, default=200, help="BER p-grid size; ignored in conditional mode")
    parser.add_argument("--metric", choices=("information", "full", "both"), default="information")
    parser.add_argument("--thin", type=float, default=0.5,
                        help="positive size multiplier for curve widths, marker diameters and marker strokes; "
                        "smaller is thinner, 1 uses the previous curve width and marker area")
    parser.add_argument("--allow-mixed-codewords", action="store_true",
                        help="pool zero/random inputs using linear-code BDD translation equivariance; "
                        "absent legacy convention means random; repeated seeds remain forbidden")
    args = parser.parse_args(argv)
    require(math.isfinite(args.thin) and args.thin > 0,
            "--thin must be finite and greater than zero")
    require(args.mode == "conditional" or args.points >= 2, "--points must be at least 2")
    require(args.output.suffix.lower() == ".svg", "--output must be an .svg path")
    csv_path = args.csv or args.output.with_suffix(".csv")
    reports = discover_reports(args.inputs)
    protected = {p.resolve() for path in reports for p in (path, path.parent / "metadata.json")}
    outputs = [args.output.resolve(), csv_path.resolve()]
    if args.export is not None:
        outputs.append(args.export.resolve())
    require(len(set(outputs)) == len(outputs) and not protected.intersection(outputs),
            "output paths must be distinct from each other and inputs")
    groups = pool_reports(reports, args.allow_mixed_codewords)
    print("Warning: metadata lacks source revision; decoder implementation compatibility cannot be verified. "
          "Missing strata are unknown; zero observed errors do not establish zero BER. "
          "No extrapolation or MSE fit is performed. "
          f"Validated {len(reports)} report(s); dimensions and Cantor coordinates checked.", file=sys.stderr)
    metrics = list(METRICS) if args.metric == "both" else [args.metric]
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise ValueError("plotting requires matplotlib; numerical core uses only the standard library") from error
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "svg.hashsalt": "gf256-product-monte-carlo"})
    ps = ([0.008 + (0.0045 - 0.008) * i / (args.points - 1) for i in range(args.points)]
          if args.mode != "conditional" else [])
    plot_results(groups, metrics, args.mode, ps, args.output, csv_path, plt, thin=args.thin)
    if args.export is not None:
        export_conditional(groups, metrics, args.export)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        sys.exit(1)
