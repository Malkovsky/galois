#!/usr/bin/env python3
"""Collect a chunk-interleaved C++/Go RS benchmark comparison."""

import argparse
import json
import os
import random
import shutil
import statistics
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


CPP_FILTER = (
    r"^RS/(LCH/Owned/(AVX2|GFNI256)|XDRS/Native|ISA-L/Native|"
    r"Jerasure/Native|Leopard/Native)/(Encode|DecodeMax)/"
)
COMPARISON_CASES = (
    (8, 248),
    (16, 240),
    (32, 224),
    (64, 192),
    (128, 128),
    (192, 64),
    (224, 32),
    (240, 16),
    (248, 8),
)
LEOPARD_CASES = COMPARISON_CASES[4:]


def expected_run_names():
    result = set()
    full_grid_families = (
        "RS/LCH/Owned/AVX2",
        "RS/LCH/Owned/GFNI256",
        "RS/XDRS/Native",
        "RS/ISA-L/Native",
        "RS/Jerasure/Native",
        "RS/Klauspost/Native",
    )
    for family in full_grid_families:
        for operation in ("Encode", "DecodeMax"):
            for data, recovery in COMPARISON_CASES:
                result.add(
                    f"{family}/{operation}/K:{data}/R:{recovery}/bytes:1024"
                )
    for operation in ("Encode", "DecodeMax"):
        for data, recovery in LEOPARD_CASES:
            result.add(
                f"RS/Leopard/Native/{operation}/K:{data}/R:{recovery}/bytes:1024"
            )
    return result


def positive_int(value):
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return result


def nonnegative_float(value):
    result = float(value)
    if result < 0:
        raise argparse.ArgumentTypeError("value must be non-negative")
    return result


def integer(value):
    return int(value, 0)


def parse_arguments():
    project_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cpp-binary",
        type=Path,
        default=project_root / "build/benchmarks-preset/rs_benchmarks",
    )
    parser.add_argument(
        "--go-binary",
        type=Path,
        default=project_root
        / "build/benchmarks-preset/klauspost_rs_benchmarks",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=project_root / "build/benchmarks-preset/rs_long.json",
    )
    parser.add_argument(
        "--cpu", type=int, default=int(os.environ.get("BENCH_CPU", "0"))
    )
    parser.add_argument("--repetitions", type=positive_int, default=100)
    parser.add_argument("--chunk-repetitions", type=positive_int, default=10)
    parser.add_argument("--warmup-seconds", type=nonnegative_float, default=1.0)
    parser.add_argument("--min-time-seconds", type=nonnegative_float, default=0.2)
    parser.add_argument("--seed", type=integer, default=0x6D2B79F5)
    parser.add_argument(
        "--dry-run", action="store_true", help="print the randomized chunk plan"
    )
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="run one verified iteration per row instead of collecting timings",
    )
    return parser.parse_args()


def repetition_chunks(total, chunk_size):
    chunks = []
    remaining = total
    while remaining:
        current = min(remaining, chunk_size)
        chunks.append(current)
        remaining -= current
    return chunks


def cpp_command(args, repetitions, output):
    command = [
        "taskset",
        "-c",
        str(args.cpu),
        str(args.cpp_binary),
        f"--benchmark_filter={CPP_FILTER}",
        f"--benchmark_min_warmup_time={0 if args.smoke else args.warmup_seconds}",
        f"--benchmark_repetitions={repetitions}",
        "--benchmark_enable_random_interleaving=true",
        "--benchmark_display_aggregates_only=true",
        f"--benchmark_out={output}",
        "--benchmark_out_format=json",
    ]
    command.append(
        "--benchmark_min_time=1x"
        if args.smoke
        else f"--benchmark_min_time={args.min_time_seconds}s"
    )
    return command


def go_command(args, repetitions, output, chunk_index):
    command = [
        "taskset",
        "-c",
        str(args.cpu),
        str(args.go_binary),
        f"--warmup={0 if args.smoke else args.warmup_seconds}s",
        f"--repetitions={repetitions}",
        f"--seed={args.seed + chunk_index}",
        f"--output={output}",
    ]
    if args.smoke:
        command.append("--iterations=1")
    else:
        command.append(f"--min-time={args.min_time_seconds}s")
    return command


def load_chunk(path, kind, expected_repetitions):
    with path.open(encoding="utf-8") as source:
        document = json.load(source)
    context = document.get("context", {})
    if context.get("cpu_pinning") != "pinned":
        raise RuntimeError(f"{path} was not pinned to one CPU")

    rows = []
    counts = defaultdict(int)
    for row in document.get("benchmarks", []):
        if row.get("aggregate_name") is not None:
            continue
        if row.get("error_occurred"):
            raise RuntimeError(f"{path} contains a skipped row: {row}")
        name = row.get("run_name", row.get("name", ""))
        if kind == "go" and not name.startswith("RS/Klauspost/Native/"):
            continue
        if kind == "cpp" and not name.startswith("RS/"):
            continue
        counts[name] += 1
        rows.append(dict(row))

    if not rows or any(count != expected_repetitions for count in counts.values()):
        raise RuntimeError(f"{path} has incomplete raw repetitions: {dict(counts)}")
    return context, rows


def aggregate(values, name):
    if name == "mean":
        return statistics.fmean(values)
    if name == "median":
        return statistics.median(values)
    if name == "stddev":
        return statistics.stdev(values) if len(values) > 1 else 0.0
    if name == "cv":
        mean = statistics.fmean(values)
        return 0.0 if mean == 0 else aggregate(values, "stddev") / mean
    raise ValueError(f"unknown aggregate {name}")


def aggregate_rows(run_name, rows):
    result = []
    for name in ("mean", "median", "stddev", "cv"):
        row = {
            "name": f"{run_name}_{name}",
            "family_index": rows[0]["family_index"],
            "per_family_instance_index": rows[0]["per_family_instance_index"],
            "run_name": run_name,
            "run_type": "aggregate",
            "repetitions": len(rows),
            "threads": 1,
            "aggregate_name": name,
            "aggregate_unit": "percentage" if name == "cv" else "time",
            "iterations": len(rows),
            "real_time": aggregate([row["real_time"] for row in rows], name),
            "cpu_time": aggregate([row["cpu_time"] for row in rows], name),
            "time_unit": "ns",
            "bytes_per_second": aggregate(
                [row["bytes_per_second"] for row in rows], name
            ),
        }
        result.append(row)
    return result


def merge_chunks(chunks, args):
    contexts = {}
    rows_by_name = defaultdict(list)
    for kind, path, repetitions in chunks:
        context, rows = load_chunk(path, kind, repetitions)
        contexts.setdefault(kind, context)
        if contexts[kind].get("host_name") != context.get("host_name"):
            raise RuntimeError(f"{kind} chunks use different hosts")
        for row in rows:
            rows_by_name[row.get("run_name", row["name"])].append(row)

    if set(contexts) != {"cpp", "go"}:
        raise RuntimeError("both C++ and Go chunks are required")
    if contexts["cpp"].get("host_name") != contexts["go"].get("host_name"):
        raise RuntimeError("C++ and Go chunks use different hosts")
    expected_names = expected_run_names()
    actual_names = set(rows_by_name)
    if actual_names != expected_names:
        missing = sorted(expected_names - actual_names)
        unexpected = sorted(actual_names - expected_names)
        raise RuntimeError(
            f"benchmark set mismatch: missing={missing}, unexpected={unexpected}"
        )

    family_names = sorted({name.split("/K:", 1)[0] for name in rows_by_name})
    family_indices = {name: index for index, name in enumerate(family_names)}
    raw_rows = []
    aggregates = []
    for run_name in sorted(rows_by_name):
        rows = rows_by_name[run_name]
        if len(rows) != args.repetitions:
            raise RuntimeError(
                f"{run_name} has {len(rows)} repetitions, want {args.repetitions}"
            )
        family_index = family_indices[run_name.split("/K:", 1)[0]]
        for repetition_index, row in enumerate(rows):
            row["name"] = run_name
            row["run_name"] = run_name
            row["run_type"] = "iteration"
            row["repetitions"] = args.repetitions
            row["repetition_index"] = repetition_index
            row["family_index"] = family_index
            row.pop("aggregate_name", None)
            row.pop("aggregate_unit", None)
            raw_rows.append(row)
        aggregates.extend(aggregate_rows(run_name, rows))

    context = dict(contexts["cpp"])
    for key, value in contexts["go"].items():
        context[f"klauspost_{key}"] = value
    context.update(
        {
            "cross_process_interleave": "randomized_chunks",
            "cross_process_seed": args.seed,
            "chunk_repetitions": args.chunk_repetitions,
        }
    )
    return {"context": context, "benchmarks": raw_rows + aggregates}


def main():
    args = parse_arguments()
    if args.cpu < 0:
        raise SystemExit("--cpu must be non-negative")
    if shutil.which("taskset") is None:
        raise SystemExit("taskset is required")
    if args.smoke:
        args.repetitions = 1
        args.chunk_repetitions = 1

    chunks = []
    for kind in ("cpp", "go"):
        for repetitions in repetition_chunks(
            args.repetitions, args.chunk_repetitions
        ):
            chunks.append((kind, repetitions))
    random.Random(args.seed).shuffle(chunks)

    if args.dry_run:
        for index, (kind, repetitions) in enumerate(chunks):
            print(f"{index + 1:02d}: {kind} x {repetitions}")
        return

    for binary in (args.cpp_binary, args.go_binary):
        if not binary.is_file():
            raise SystemExit(f"benchmark binary does not exist: {binary}")
    args.output.parent.mkdir(parents=True, exist_ok=True)

    completed = []
    with tempfile.TemporaryDirectory(
        prefix="rs-long-", dir=args.output.parent
    ) as temporary:
        temporary_path = Path(temporary)
        for index, (kind, repetitions) in enumerate(chunks):
            output = temporary_path / f"{index:03d}-{kind}.json"
            command = (
                cpp_command(args, repetitions, output)
                if kind == "cpp"
                else go_command(args, repetitions, output, index)
            )
            print(f"[{index + 1}/{len(chunks)}] {kind} x {repetitions}", flush=True)
            subprocess.run(command, check=True, stdout=subprocess.DEVNULL)
            completed.append((kind, output, repetitions))

        document = merge_chunks(completed, args)
        with args.output.open("w", encoding="utf-8") as destination:
            json.dump(document, destination, indent=2, sort_keys=True)
            destination.write("\n")


if __name__ == "__main__":
    main()
