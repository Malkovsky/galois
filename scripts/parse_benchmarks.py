import argparse
import json
import matplotlib.pyplot as plt
import re
from collections import defaultdict


def plot_google_benchmark(json_file, out_file, metric="cpu_time", title=None):
    """
    Plot Google Benchmark results as log-scale line plots.

    Parameters
    ----------
    json_file : str
        Path to Google Benchmark JSON output.
    out_file : str
        Path to dump plot graph.
    metric : str
        Which metric to plot ("real_time" or "cpu_time").
    title : str or None
        Plot title. If None, defaults to metric.
    """
    with open(json_file, "r") as f:
        data = json.load(f)

    benchmarks = data.get("benchmarks", [])
    if not benchmarks:
        raise ValueError("No 'benchmarks' found in JSON file")

    groups = defaultdict(list)

    # Regex to split name into "base" and "param"
    name_re = re.compile(r"([^/]+)/?(n:\d+)?")

    for b in benchmarks:
        if metric not in b:
            continue
        match = name_re.match(b["name"])
        if not match:
            continue
        base, param = match.groups()
        param = int(param[2:]) if param else None
        groups[base].append((param, b[metric]))

    plt.figure(figsize=(10, 6))
    plt.yscale("log")

    for base, items in groups.items():
        # Sort by parameter if available
        items.sort(key=lambda x: (x[0] is None, x[0]))
        xs = [p if p is not None else i for i, (p, _) in enumerate(items)]
        ys = [v for _, v in items]
        plt.plot(xs, ys, marker="o", label=base)

    plt.xscale("log")
    plt.xlabel("Size")
    plt.ylabel("CPU time, ns")
    plt.title(title if title else f"Google Benchmark: {metric}")
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(out_file)


def main():
    parser = argparse.ArgumentParser(
        description="Plot Google Benchmark JSON output with log-scale y-axis"
    )
    parser.add_argument(
        "filename",
        help="Path to Google Benchmark JSON output file (generated with --benchmark_format=json)",
    )
    parser.add_argument(
        "--metric",
        "-m",
        default="real_time",
        choices=["real_time", "cpu_time"],
        help="Metric to plot (default: real_time)",
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Output file name to save graph",
    )

    parser.add_argument("--title", "-t", default=None, help="Optional plot title")
    args = parser.parse_args()

    # Call the plotting function defined earlier
    plot_google_benchmark(
        args.filename, args.output, metric=args.metric, title=args.title
    )


if __name__ == "__main__":
    main()
