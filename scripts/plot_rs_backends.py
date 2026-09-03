#!/usr/bin/env python3
"""Generate the LCH and external RS backend comparison report graphs.

The combined report accepts matched C++ and Go JSON inputs and labels the owned
implementations LCH+AVX2 and LCH+GFNI. The GFNI-width report retains its separate
matched input so system-load windows are not mixed.
"""

import argparse
import json
import re
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

GIB = 1024**3
K_VALUES = (8, 16, 32, 64, 128, 192, 224, 240, 248)
K_POSITIONS = {value: index for index, value in enumerate(K_VALUES)}
LEOPARD_K_VALUES = (128, 192, 224, 240, 248)
COLORS = {
    "Native": "#5f6368",
    "AVX2": "#2563eb",
    "GFNI256": "#7c3aed",
    "GFNI512": "#db2777",
    "LCH+AVX2": "#2563eb",
    "LCH+GFNI": "#7c3aed",
    "XDRS": "#5f6368",
    "ISA-L": "#ea580c",
    "Jerasure": "#16a34a",
    "Klauspost": "#dc2626",
    "Leopard": "#111827",
}
MARKERS = {
    "Native": "o",
    "AVX2": "s",
    "GFNI256": "D",
    "GFNI512": "^",
    "LCH+AVX2": "s",
    "LCH+GFNI": "D",
    "XDRS": "o",
    "ISA-L": "v",
    "Jerasure": "P",
    "Klauspost": "X",
    "Leopard": "h",
}
COMPARISON_BACKENDS = (
    "LCH+AVX2",
    "LCH+GFNI",
    "XDRS",
    "ISA-L",
    "Jerasure",
    "Klauspost",
)

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.titlesize": 15,
        "axes.labelsize": 12,
        "legend.fontsize": 11,
        "svg.hashsalt": "gf256-rs-backend-report",
    }
)


def load_json(path):
    with path.open("r", encoding="utf-8") as source:
        return json.load(source)


def load_raw_medians(path, pattern, make_key):
    values = defaultdict(list)
    regex = re.compile(pattern)
    for row in load_json(path)["benchmarks"]:
        if "aggregate_name" in row:
            continue
        match = regex.match(row["name"])
        if match:
            values[make_key(match)].append(row["bytes_per_second"] / GIB)
    return {key: statistics.median(samples) for key, samples in values.items()}


def load_lch_and_xdrs(path):
    return load_raw_medians(
        path,
        r"^RS/(?:LCH|XDRS)/(Owned/(AVX2|GFNI256)|Native)/"
        r"(Encode|DecodeMax)/K:(\d+)/",
        lambda match: (
            match.group(3),
            int(match.group(4)),
            "Native" if match.group(1) == "Native" else match.group(2),
        ),
    )


def load_comparison(paths):
    values = defaultdict(list)
    aggregate_medians = {}
    regex = re.compile(
        r"^RS/(?:(?:LCH|XDRS)/Owned/(AVX2|GFNI256)|"
        r"(XDRS|ISA-L|Jerasure|Klauspost)/Native)/"
        r"(Encode|DecodeMax)/K:(\d+)/"
    )
    for path in paths:
        for row in load_json(path)["benchmarks"]:
            run_name = row.get("run_name", row["name"])
            match = regex.match(run_name)
            if not match:
                continue
            if match.group(1) == "AVX2":
                backend = "LCH+AVX2"
            elif match.group(1) == "GFNI256":
                backend = "LCH+GFNI"
            else:
                backend = match.group(2)
            key = (match.group(3), int(match.group(4)), backend)
            throughput = row["bytes_per_second"] / GIB
            if "aggregate_name" not in row:
                values[key].append(throughput)
            elif row.get("aggregate_name") == "median":
                aggregate_medians[key] = throughput

    result = dict(aggregate_medians)
    result.update(
        {key: statistics.median(samples) for key, samples in values.items()}
    )
    return result


def validate_comparison_contexts(paths):
    contexts = [load_json(path).get("context", {}) for path in paths]
    hosts = {context.get("host_name") for context in contexts}
    if None in hosts or len(hosts) != 1:
        raise ValueError(f"comparison inputs use different hosts: {hosts}")
    pinning = {context.get("cpu_pinning") for context in contexts}
    if None in pinning or len(pinning) != 1:
        raise ValueError(
            f"comparison inputs use different CPU-pinning states: {pinning}"
        )


def load_gfni_widths(path):
    return load_raw_medians(
        path,
        r"^(?:LCH/Owned|XDRS/Polished/(?:Low|High))/(Encode|DecodeMax)/"
        r"GFNI(256|512)Affine/Radix4/K:(\d+)/",
        lambda match: (
            match.group(1),
            int(match.group(3)),
            "GFNI" + match.group(2),
        ),
    )


def load_leopard(path):
    values = defaultdict(list)
    aggregate_medians = {}
    regex = re.compile(r"^RS/Leopard/Native/(Encode|DecodeMax)/K:(\d+)/")
    for row in load_json(path)["benchmarks"]:
        match = regex.match(row.get("run_name", row["name"]))
        if not match:
            continue
        key = (match.group(1), int(match.group(2)), "Leopard")
        throughput = row["bytes_per_second"] / GIB
        if "aggregate_name" not in row:
            values[key].append(throughput)
        elif row.get("aggregate_name") == "median":
            aggregate_medians[key] = throughput
    result = dict(aggregate_medians)
    result.update(
        {key: statistics.median(samples) for key, samples in values.items()}
    )
    return result


def require_points(data, operations, backends, k_values, description):
    missing = [
        (operation, k, backend)
        for operation in operations
        for backend in backends
        for k in k_values
        if (operation, k, backend) not in data
    ]
    if missing:
        raise ValueError(f"{description} is missing benchmark points: {missing}")


def configure_axis(axis, title, ticks, limits):
    axis.set_title(title, loc="left", fontweight="semibold")
    axis.set_yscale("log", base=2)
    axis.set_ylim(*limits)
    axis.set_yticks(ticks)
    axis.set_yticklabels([str(value) for value in ticks])
    axis.set_ylabel("Input throughput (GiB/s)")
    axis.set_xticks(range(len(K_VALUES)))
    axis.set_xticklabels(K_VALUES)
    axis.set_xlim(-0.2, len(K_VALUES) - 0.8)
    axis.grid(
        True,
        which="major",
        color="#d1d5db",
        linewidth=0.8,
        alpha=0.75,
    )
    axis.grid(
        True,
        which="minor",
        axis="y",
        color="#e5e7eb",
        linewidth=0.5,
        alpha=0.5,
    )
    axis.set_axisbelow(True)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


def plot_line(axis, data, backend, operation, k_values, label, dashed=False):
    return axis.plot(
        [K_POSITIONS[k] for k in k_values],
        [data[(operation, k, backend)] for k in k_values],
        color=COLORS[backend],
        linestyle="--" if dashed else "-",
        linewidth=2.4,
        marker=MARKERS[backend],
        markersize=7,
        markerfacecolor="white" if dashed else COLORS[backend],
        markeredgecolor=COLORS[backend],
        markeredgewidth=1.6,
        label=label,
    )[0]


def make_figure(title):
    figure, axes = plt.subplots(2, 1, figsize=(10.5, 12), sharex=True)
    figure.suptitle(
        title,
        x=0.08,
        y=0.965,
        ha="left",
        fontsize=22,
        fontweight="bold",
    )
    configure_axis(
        axes[0],
        "Encode throughput",
        (0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32),
        (0.05, 32),
    )
    configure_axis(
        axes[1],
        "DecodeMax input throughput",
        (0.015625, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16),
        (0.012, 16),
    )
    axes[1].set_xlabel("Dimension K")
    return figure, axes


def save_figure(figure, output_dir, stem):
    metadata = {"Creator": "scripts/plot_rs_backends.py", "Date": None}
    figure.savefig(
        output_dir / f"{stem}.svg",
        facecolor="white",
        metadata=metadata,
    )
    svg_path = output_dir / f"{stem}.svg"
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_path.read_text().splitlines())
        + "\n",
        encoding="utf-8",
    )
    plt.close(figure)


def generate_combined_report(data, leopard, backends, output_dir):
    figure, axes = make_figure("RS(256, K) backend comparison")
    for backend in backends:
        plot_line(
            axes[0],
            data,
            backend,
            "Encode",
            K_VALUES,
            backend,
        )

    handles = []
    labels = []
    for backend in backends:
        line = plot_line(
            axes[1],
            data,
            backend,
            "DecodeMax",
            K_VALUES,
            backend,
        )
        handles.append(line)
        labels.append(line.get_label())
    leopard_line = plot_line(
        axes[1],
        leopard,
        "Leopard",
        "DecodeMax",
        LEOPARD_K_VALUES,
        "Leopard",
        dashed=True,
    )
    handles.append(leopard_line)
    labels.append(leopard_line.get_label())

    figure.legend(
        handles,
        labels,
        loc="lower center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 0.025),
    )
    figure.subplots_adjust(
        top=0.91, bottom=0.13, left=0.11, right=0.97, hspace=0.28
    )
    save_figure(figure, output_dir, "rs_backend_comparison")


def generate_gfni_report(data, output_dir):
    figure, axes = make_figure("LCH GFNI vector-width comparison")
    for operation, axis in (("Encode", axes[0]), ("DecodeMax", axes[1])):
        plot_line(
            axis,
            data,
            "GFNI256",
            operation,
            K_VALUES,
            "GFNI256 affine",
        )
        plot_line(
            axis,
            data,
            "GFNI512",
            operation,
            K_VALUES,
            "GFNI512 affine",
        )
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="lower center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 0.035),
    )
    figure.subplots_adjust(
        top=0.91, bottom=0.11, left=0.11, right=0.97, hspace=0.28
    )
    save_figure(figure, output_dir, "xdrs_gfni_width_comparison")


def parse_arguments():
    project_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--comparison-results",
        type=Path,
        nargs="+",
        help="matched C++ and Go JSON files for the full backend report",
    )
    parser.add_argument(
        "--combined-only",
        action="store_true",
        help="generate only rs_backend_comparison from comparison inputs",
    )
    parser.add_argument(
        "--xdrs-results",
        type=Path,
        default=project_root
        / "build/benchmarks-preset/xdrs_night_interleaved.json",
    )
    parser.add_argument(
        "--gfni-results",
        type=Path,
        default=project_root
        / "build/benchmarks-preset/xdrs_gfni256_vs_gfni512_night.json",
    )
    parser.add_argument(
        "--leopard-results",
        type=Path,
        default=project_root
        / "build/benchmarks-preset/rs_unified_backends.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=project_root / "docs/images",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    xdrs = None
    if args.comparison_results:
        validate_comparison_contexts(args.comparison_results)
        comparison = load_comparison(args.comparison_results)
        leopard = {}
        for path in args.comparison_results:
            leopard.update(load_leopard(path))
        comparison_backends = COMPARISON_BACKENDS
    else:
        xdrs = load_lch_and_xdrs(args.xdrs_results)
        comparison = {
            (operation, k, "XDRS" if backend == "Native" else
             "LCH+AVX2" if backend == "AVX2" else "LCH+GFNI"): value
            for (operation, k, backend), value in xdrs.items()
        }
        leopard = load_leopard(args.leopard_results)
        comparison_backends = ("LCH+AVX2", "LCH+GFNI", "XDRS")

    require_points(
        leopard,
        ("DecodeMax",),
        ("Leopard",),
        LEOPARD_K_VALUES,
        "native Leopard results",
    )
    require_points(
        comparison,
        ("Encode", "DecodeMax"),
        comparison_backends,
        K_VALUES,
        "combined RS results",
    )

    generate_combined_report(
        comparison, leopard, comparison_backends, args.output_dir
    )
    if not args.combined_only:
        if xdrs is None:
            xdrs = load_lch_and_xdrs(args.xdrs_results)
        gfni = load_gfni_widths(args.gfni_results)
        require_points(
            xdrs,
            ("Encode", "DecodeMax"),
            ("Native", "AVX2", "GFNI256"),
            K_VALUES,
            "XDRS results",
        )
        require_points(
            gfni,
            ("Encode", "DecodeMax"),
            ("GFNI256", "GFNI512"),
            K_VALUES,
            "GFNI width results",
        )
        generate_gfni_report(gfni, args.output_dir)


if __name__ == "__main__":
    main()
