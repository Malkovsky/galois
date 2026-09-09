#!/usr/bin/env python3
"""Frozen schema-1 coordinator for differential tests, never installed."""

import argparse
import array
import concurrent.futures
from contextlib import ExitStack
import ctypes
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import secrets
import sys
import struct
import time

SCHEMA = 1
FLOYD = "splitmix64 domain seeds; mt19937_64; rejection modulo; Floyd complement v1"
FISHER_YATES = "splitmix64 domain seeds; mt19937_64; rejection modulo; persistent Fisher-Yates complement v1; replay saved flips"
# Little endian, LSB-first bit = byte * 8 + bit_in_byte, row-major block.
# Per record: batch, trial, k, count, complement, 22 metrics, uint32 positions,
# SHA256(header + positions). File header binds format/version and run identity.
FLIP_MAGIC = b"RSFLIP01"
FLIP_HEADER = struct.Struct("<QQIIB7x22Q")
U64 = (1 << 64) - 1
INFORMATION_BYTES = 56896
METRICS = (
    "initial full block corrupted bits", "initial full block corrupted bytes",
    "initial information corrupted bits", "initial information corrupted bytes",
    "residual full block bits", "residual full block bytes",
    "residual information bits", "residual information bytes",
    "message failures", "full block failures", "zero syndrome outcomes",
    "zero syndrome wrong full blocks", "directional passes",
    "accepted bit changes", "accepted byte changes",
    "strong accepted bit changes", "strong accepted byte changes",
    "weak accepted bit changes", "weak accepted byte changes",
    "strong lines visited", "weak lines visited", "pass limit outcomes",
)


def timestamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def integer(text):
    if not text.isascii() or not text.isdecimal():
        raise argparse.ArgumentTypeError("expected an unsigned decimal integer")
    value = int(text)
    if value > U64:
        raise argparse.ArgumentTypeError("integer exceeds uint64")
    return value


def native():
    here = Path(__file__).resolve().parent
    candidates = [here / "product_monte_carlo_native.so",
                  here.parent / "lib" / "product_monte_carlo_native.so",
                  here.parent / "lib64" / "product_monte_carlo_native.so"]
    path = next((p for p in candidates if p.is_file()), None)
    if path is None:
        raise ValueError("native trial library not found beside CLI or in ../lib[64]")
    lib = ctypes.CDLL(str(path))
    lib.product_batch_k.argtypes = [ctypes.c_uint64] * 4
    lib.product_batch_k.restype = ctypes.c_uint64
    lib.product_trial.argtypes = ([ctypes.c_uint64] * 5 + [ctypes.c_int] * 2 +
                                 [ctypes.POINTER(ctypes.c_uint64)])
    lib.product_trial.restype = ctypes.c_int
    lib.product_trial_flips.argtypes = lib.product_trial.argtypes + [
        ctypes.c_int, ctypes.POINTER(ctypes.c_uint32)]
    lib.product_trial_flips.restype = ctypes.c_int
    lib.product_trial_reference.argtypes = lib.product_trial_flips.argtypes + [
        ctypes.c_int, ctypes.POINTER(ctypes.c_uint8)]
    lib.product_trial_reference.restype = ctypes.c_int
    lib.product_trials.argtypes = lib.product_trial.argtypes + [
        ctypes.c_uint64, ctypes.c_int, ctypes.POINTER(ctypes.c_uint32),
        ctypes.POINTER(ctypes.c_uint64)]
    lib.product_trials.restype = ctypes.c_int
    lib.product_interrupt_install.argtypes = []
    lib.product_interrupt_install.restype = ctypes.c_int
    lib.product_interrupted.argtypes = []
    lib.product_interrupted.restype = ctypes.c_int
    return lib


def strict_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_number(value):
    raise ValueError(f"noninteger JSON number: {value}")


def loads(text):
    return json.loads(text, object_pairs_hook=strict_object,
                      parse_float=reject_number, parse_constant=reject_number)


def canonical(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=True, separators=(",", ":"))


def identity(metadata):
    return hashlib.sha256(canonical(metadata).encode("ascii")).hexdigest()


def sync_directory(directory):
    fd = os.open(directory, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def atomic_json(path, value):
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="ascii") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, ensure_ascii=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(tmp, path)
    sync_directory(path.parent)


def empty_stats():
    return {name: {"sum": 0, "squared sum": 0} for name in METRICS}


def add_stats(destination, source):
    for name in METRICS:
        for field in ("sum", "squared sum"):
            destination[name][field] += source[name][field]


def stats_for(results):
    stats = empty_stats()
    for values in results:
        for name, value in zip(METRICS, values):
            stats[name]["sum"] += value
            stats[name]["squared sum"] += value * value
    return stats


class Aggregate:
    def __init__(self, run_identity):
        self.run_identity = run_identity
        self.by_k = {}
        self.count = 0
        self.stats = empty_stats()

    def add(self, record):
        k, count = record["flipped bit count"], record["trial count"]
        entry = self.by_k.setdefault(k, {"flipped bit count": k, "trial count": 0,
                                        "statistics": empty_stats()})
        entry["trial count"] += count
        add_stats(entry["statistics"], record["statistics"])
        self.count += count
        add_stats(self.stats, record["statistics"])

    def summary(self):
        return {"schema revision": SCHEMA, "run identity": self.run_identity,
                "overall": {"trial count": self.count, "statistics": self.stats},
                "by flipped bit count": [self.by_k[k] for k in sorted(self.by_k)]}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def natural(value, maximum=None):
    return type(value) is int and value >= 0 and (maximum is None or value <= maximum)


def flip_bytes(batch, trial, k, metrics, positions):
    header = FLIP_HEADER.pack(batch, trial, k, min(k, 524288 - k), k > 262144, *metrics)
    payload = bytes(positions)
    if sys.byteorder != "little":
        converted = array.array("I")
        converted.frombytes(payload)
        converted.byteswap()
        payload = converted.tobytes()
    body = header + payload
    return body + hashlib.sha256(body).digest()


def read_flips(stream, record, replay, settings, lib, codeword):
    require(stream.tell() == record["flip start"], "noncontiguous flip offsets")
    values = []
    for trial in range(record["first trial index"], record["past last trial index"]):
        header = stream.read(FLIP_HEADER.size)
        require(len(header) == FLIP_HEADER.size, "truncated flip header")
        batch, index, k, count, complement, *metrics = FLIP_HEADER.unpack(header)
        require((batch, index, k) == (record["batch id"], trial, record["flipped bit count"]),
                "flip trial identity mismatch")
        require(count == min(k, 524288 - k) and complement == int(k > 262144),
                "invalid flip count/complement")
        payload = stream.read(count * 4)
        checksum = stream.read(32)
        require(len(payload) == count * 4 and hashlib.sha256(header + payload).digest() == checksum,
                "truncated or corrupt flip record")
        positions = array.array("I")
        positions.frombytes(payload)
        if sys.byteorder != "little":
            positions.byteswap()
        require(all(p < 524288 for p in positions) and len(set(positions)) == count,
                "invalid or duplicate flip position")
        if replay:
            output = (ctypes.c_uint64 * len(METRICS))()
            buffer = (ctypes.c_uint32 * count)(*positions)
            status = lib.product_trial_reference(settings["root seed"], batch, trial, k,
                settings["maximum directional passes"], settings["anchors"], settings["binary image"],
                output, 2, buffer, codeword == "random", None)
            require(status == 0 and list(output) == metrics,
                    f"replay mismatch batch={batch} trial={trial} status={status}")
        values.append(metrics)
    require(stream.tell() == record["flip end"], "flip end offset mismatch")
    require(stats_for(values) == record["statistics"], "flip metrics disagree with journal")


def validate_settings(settings):
    expected = {"root seed", "batch size", "batches", "threads", "minimum flipped bits",
                "maximum flipped bits", "maximum directional passes", "anchors",
                "binary image", "checkpoint trials", "report seconds", "fsync seconds"}
    require(type(settings) is dict and set(settings) == expected, "incompatible settings")
    for key in expected - {"anchors", "binary image"}:
        require(natural(settings[key], U64), f"invalid {key}")
    require(type(settings["anchors"]) is bool and type(settings["binary image"]) is bool,
            "gates must be boolean")
    require(0 <= settings["minimum flipped bits"] <= settings["maximum flipped bits"] <= 524288,
            "require 0 <= minimum flipped bits <= maximum flipped bits <= 524288")
    for key, maximum in (("threads", 1024), ("checkpoint trials", 4096),
                         ("report seconds", 86400), ("fsync seconds", 86400)):
        require(1 <= settings[key] <= maximum, f"{key} must be in [1,{maximum}]")
    require(settings["batch size"] > 0, "batch size must be positive")
    require(2 <= settings["maximum directional passes"] <= 1000000,
            "maximum directional passes must be in [2,1000000]")


def validate_record(record, metadata, run_identity, expected_id, previous, lib):
    fields = {"schema revision", "run identity", "increment id", "batch id",
              "first trial index", "past last trial index", "trial count",
              "flipped bit count", "statistics"}
    if metadata["random algorithm"] == FISHER_YATES:
        fields |= {"flip start", "flip end"}
    require(type(record) is dict and set(record) == fields, "incompatible journal record")
    require(record["schema revision"] == SCHEMA and type(record["schema revision"]) is int and
            record["run identity"] == run_identity, "incompatible journal identity/schema")
    for key in fields - {"run identity", "statistics"}:
        require(natural(record[key]), f"invalid record {key}")
    require(record["increment id"] == expected_id, "duplicate or out-of-order increment id")
    s = metadata["settings"]
    batch, first, end = (record[key] for key in
                         ("batch id", "first trial index", "past last trial index"))
    require(batch <= U64 and (s["batches"] == 0 or batch < s["batches"]), "invalid batch id")
    require(0 <= first < end <= s["batch size"] and end - first == record["trial count"] and
            end - first <= s["checkpoint trials"], "invalid trial range/count")
    if previous is not None:
        old_batch, old_end = previous
        require(batch >= old_batch and (batch != old_batch or first >= old_end),
                "overlapping or out-of-order trial ranges")
    k = lib.product_batch_k(s["root seed"], batch, s["minimum flipped bits"], s["maximum flipped bits"])
    require(record["flipped bit count"] == k, "batch k disagrees with seed/settings")
    stats, count = record["statistics"], record["trial count"]
    require(type(stats) is dict and set(stats) == set(METRICS), "incompatible statistics")
    bounds = [524288, 65536, 455168, 56896, 524288, 65536, 455168, 56896,
              1, 1, 1, 1, s["maximum directional passes"]] + [
                  s["maximum directional passes"] * 524288] * 8 + [1]
    for name, bound in zip(METRICS, bounds):
        item = stats[name]
        require(type(item) is dict and set(item) == {"sum", "squared sum"}, "invalid moment fields")
        total, square = item["sum"], item["squared sum"]
        require(natural(total, count * bound) and natural(square, count * bound * bound),
                f"invalid moment: {name}")
        require(total * total <= count * square and total <= square <= bound * total,
                f"inconsistent moments: {name}")
    require(stats[METRICS[0]] == {"sum": count * k, "squared sum": count * k * k},
            "initial channel is not exact k")
    for total, strong, weak in ((13, 15, 17), (14, 16, 18)):
        require(stats[METRICS[total]]["sum"] == stats[METRICS[strong]]["sum"] +
                stats[METRICS[weak]]["sum"], "directional accepted totals disagree")
    return batch, end


def recover(directory, lib, replay=False):
    with (directory / "metadata.json").open(encoding="ascii") as stream:
        metadata = loads(stream.read(131073))
    require(type(metadata) is dict and set(metadata) - {"codeword"} == {
        "schema revision", "created at", "settings", "code", "random algorithm"}, "incompatible metadata")
    # Do not insert defaults into metadata: legacy run identities hash it as-is.
    codeword = metadata.get("codeword", "random")
    require(codeword in ("zero", "random"), "incompatible codeword convention")
    require(type(metadata["schema revision"]) is int and metadata["schema revision"] == SCHEMA and
            metadata["code"] == "RS256,224 x RS256,254 Cantor systematic row major" and
            metadata["random algorithm"] in (FLOYD, FISHER_YATES),
            "incompatible metadata schema/code/random algorithm")
    validate_settings(metadata["settings"])
    digest = identity(metadata)
    aggregate = Aggregate(digest)
    previous = None
    ignored = False
    recorded = metadata["random algorithm"] == FISHER_YATES
    require(not replay or recorded, "this run has no saved flips")
    with ExitStack() as stack:
        stream = stack.enter_context((directory / "journal.jsonl").open("rb"))
        flips = stack.enter_context((directory / "flips.bin").open("rb")) if recorded else None
        if flips:
            require(flips.read(40) == FLIP_MAGIC + bytes.fromhex(digest), "invalid flip file identity/version")
        index = 0
        while True:
            line = stream.readline(131073)
            if not line:
                break
            require(len(line) <= 131072, "journal line exceeds schema size limit")
            if not line.endswith(b"\n"):
                ignored = True
                break
            try:
                record = loads(line)
                previous = validate_record(record, metadata, digest, index, previous, lib)
                if flips:
                    read_flips(flips, record, replay, metadata["settings"], lib, codeword)
                aggregate.add(record)
            except (ValueError, TypeError, KeyError) as error:
                raise ValueError(f"journal line {index + 1}: {error}") from error
            index += 1
        if flips and flips.read(1):
            print("unreferenced flip tail ignored (not committed trials)", file=sys.stderr)
    atomic_json(directory / "summary.json", aggregate.summary())
    print(f"{timestamp()} regenerated {aggregate.count} trials; ignored incomplete tail={ignored}", file=sys.stderr)
    if replay:
        print(f"verified replay: {aggregate.count} trials, all 22 metrics match", file=sys.stderr)


def run(directory, settings, lib, sampler="floyd"):
    metadata = {"schema revision": SCHEMA, "created at": timestamp(), "settings": settings,
                "codeword": "zero",
                "code": "RS256,224 x RS256,254 Cantor systematic row major",
                "random algorithm": FISHER_YATES if sampler == "fisher-yates" else FLOYD}
    atomic_json(directory / "metadata.json", metadata)
    aggregate = Aggregate(identity(metadata))
    atomic_json(directory / "summary.json", aggregate.summary())
    if lib.product_interrupt_install() != 0:
        raise OSError("could not install native signal handlers")
    with ExitStack() as stack:
        journal = stack.enter_context((directory / "journal.jsonl").open("x", encoding="ascii"))
        log = stack.enter_context((directory / "progress.log").open("x", encoding="ascii"))
        flips = stack.enter_context((directory / "flips.bin").open("xb")) if sampler == "fisher-yates" else None
        if flips:
            flips.write(FLIP_MAGIC + bytes.fromhex(aggregate.run_identity))
            flips.flush()
            os.fsync(flips.fileno())
        tty = sys.stderr.isatty()
        bar_visible = False
        increment = batch = 0
        staged = None
        last_checkpoint = time.monotonic()

        def checkpoint():
            nonlocal staged, increment, last_checkpoint
            if staged is not None:
                if flips:
                    flips.flush()
                    os.fsync(flips.fileno())
                staged["increment id"] = increment
                journal.write(canonical(staged) + "\n")
                journal.flush()
                aggregate.add(staged)
                increment += 1
                staged = None
            last_checkpoint = time.monotonic()

        def progress(text, console=True):
            nonlocal bar_visible
            if console and bar_visible:
                print("\r\033[K", end="", file=sys.stderr)
                bar_visible = False
            line = f"{timestamp()} {text}"
            if console:
                print(line, file=sys.stderr, flush=True)
            print(line, file=log, flush=True)

        def durable():
            checkpoint()
            journal.flush()
            os.fsync(journal.fileno())
            log.flush()
            os.fsync(log.fileno())
            atomic_json(directory / "summary.json", aggregate.summary())

        progress(f"root seed={settings['root seed']} settings persisted before trials")
        durable()
        sync_directory(directory)
        sync_directory(directory.parent)
        last_sync = last_report = last_bar = time.monotonic()
        start = last_sync
        completed_total = 0
        error = None

        def throughput(now):
            # Whole simulation wall time, including coordinator work, not decoder time.
            elapsed = max(0.0, now - start)
            rate = completed_total / elapsed if elapsed > 0 else 0.0
            return (f"wall blocks/s={rate:.3f} information MiB/s={rate * INFORMATION_BYTES / (1 << 20):.3f} "
                    f"elapsed seconds={elapsed:.3f}")

        def bar(now, force=False):
            nonlocal bar_visible, last_bar
            if tty and (force or now - last_bar >= 0.2):
                fraction = batch_completed / settings["batch size"]
                width = int(30 * fraction)
                print(f"\r[{('#' * width).ljust(30)}] batch={batch} k={k} "
                      f"trials={batch_completed}/{settings['batch size']} {fraction:.1%} "
                      f"overall trials={completed_total} {throughput(now)}\033[K",
                      end="", file=sys.stderr, flush=True)
                bar_visible = True
                last_bar = now

        def trial(index, k):
            output = (ctypes.c_uint64 * len(METRICS))()
            arguments = (settings["root seed"], batch, index, k,
                         settings["maximum directional passes"],
                         settings["anchors"], settings["binary image"], output)
            if flips:
                positions = (ctypes.c_uint32 * min(k, 524288 - k))()
                status = lib.product_trial_flips(*arguments, 1, positions)
            else:
                status = lib.product_trial(*arguments)
            if status:
                raise RuntimeError(f"native trial batch={batch} index={index} failed: {status}")
            return (list(output), positions) if flips else list(output)

        # Amortize GIL/executor handoffs without starving workers at small caps.
        chunk = 4 if settings["checkpoint trials"] >= 4 * settings["threads"] else 1

        def trials(first, k):
            if chunk == 1:
                if lib.product_interrupted():
                    return [], None
                try:
                    return [trial(first, k)], None
                except Exception as exc:
                    return [], exc
            count = min(chunk, settings["batch size"] - first)
            output = (ctypes.c_uint64 * (len(METRICS) * count))()
            stride = min(k, 524288 - k)
            positions = (ctypes.c_uint32 * (stride * count))() if flips else None
            completed = ctypes.c_uint64()
            status = lib.product_trials(settings["root seed"], batch, first, k,
                settings["maximum directional passes"], settings["anchors"],
                settings["binary image"], output, count, bool(flips), positions,
                ctypes.byref(completed))
            results = []
            for i in range(completed.value):
                metrics = list(output[i * len(METRICS):(i + 1) * len(METRICS)])
                if flips:
                    saved = (ctypes.c_uint32 * stride).from_buffer(positions, i * stride * 4)
                    results.append((metrics, saved))
                else:
                    results.append(metrics)
            failure = None
            if status:
                failure = RuntimeError(
                    f"native trial batch={batch} index={first + completed.value} failed: {status}")
            return results, failure

        # Refill workers without wave barriers. Bound running plus out-of-order
        # results so a slow early trial cannot accumulate an unbounded flip log.
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=settings["threads"]) as pool:
                while not lib.product_interrupted() and not error and (
                        settings["batches"] == 0 or batch < settings["batches"]):
                    k = lib.product_batch_k(settings["root seed"], batch,
                                            settings["minimum flipped bits"], settings["maximum flipped bits"])
                    first = 0
                    cursor = 0
                    futures, completed = {}, {}
                    window = min(settings["checkpoint trials"], 2 * chunk * settings["threads"])
                    batch_completed = 0
                    progress(f"batch={batch} k={k} trials=0/{settings['batch size']} overall trials={completed_total}")
                    while futures or (first < settings["batch size"] and not lib.product_interrupted() and not error):
                        while (not error and not lib.product_interrupted() and
                               first < settings["batch size"] and
                               first - cursor < window and len(futures) < settings["threads"]):
                            end = min(first + chunk, settings["batch size"])
                            if end - cursor > window:
                                break
                            try:
                                futures[pool.submit(trials, first, k)] = (first, end)
                                first = end
                            except Exception as exc:
                                error = exc
                                break
                        if futures:
                            done, _ = concurrent.futures.wait(set(futures), timeout=0.1,
                                return_when=concurrent.futures.FIRST_COMPLETED)
                            for future in done:
                                index, end = futures.pop(future)
                                results = []
                                try:
                                    results, failure = future.result()
                                    if failure is not None:
                                        error = failure
                                except Exception as exc:
                                    error = exc
                                completed_total += len(results)
                                batch_completed += len(results)
                                for i in range(index, end):
                                    completed[i] = results[i - index] if i - index < len(results) else None
                            now = time.monotonic()
                            if now - last_checkpoint >= 1:
                                checkpoint()
                            if now - last_sync >= settings["fsync seconds"]:
                                durable()
                                last_sync = now
                            if now - last_report >= settings["report seconds"]:
                                progress(f"batch={batch} k={k} trials={batch_completed}/{settings['batch size']} "
                                         f"overall trials={completed_total} "
                                         f"persisted overall trials={aggregate.count} "
                                         f"message failures={aggregate.stats['message failures']['sum']} "
                                         f"full block failures={aggregate.stats['full block failures']['sum']} "
                                         f"{throughput(now)}", console=not tty)
                                last_report = now
                            bar(now)
                        ready = {}
                        while cursor in completed:
                            result = completed.pop(cursor)
                            if result is not None:
                                ready[cursor] = result
                            cursor += 1
                        indices = list(ready)
                        offset = 0
                        while offset < len(indices):
                            remaining = settings["checkpoint trials"] - (staged["trial count"] if staged else 0)
                            stop = offset + 1
                            while (stop < len(indices) and stop - offset < remaining and
                                   indices[stop] == indices[stop - 1] + 1):
                                stop += 1
                            flip_start = flips.tell() if flips else None
                            if flips:
                                for i in indices[offset:stop]:
                                    metrics, positions = ready[i]
                                    flips.write(flip_bytes(batch, i, k, metrics, positions))
                            record = {"schema revision": SCHEMA, "run identity": aggregate.run_identity,
                                      "increment id": increment, "batch id": batch,
                                      "first trial index": indices[offset], "past last trial index": indices[stop - 1] + 1,
                                      "trial count": stop - offset, "flipped bit count": k,
                                      "statistics": stats_for(ready[i][0] if flips else ready[i]
                                                              for i in indices[offset:stop])}
                            if flips:
                                record.update({"flip start": flip_start, "flip end": flips.tell()})
                            if staged is not None and staged["past last trial index"] != record["first trial index"]:
                                checkpoint()
                            if staged is None:
                                staged = record
                            else:
                                staged["past last trial index"] = record["past last trial index"]
                                staged["trial count"] += record["trial count"]
                                if flips:
                                    staged["flip end"] = record["flip end"]
                                add_stats(staged["statistics"], record["statistics"])
                            offset = stop
                            if staged["trial count"] >= settings["checkpoint trials"]:
                                checkpoint()
                        now = time.monotonic()
                        if (staged is not None and staged["trial count"] >= settings["checkpoint trials"]) or now - last_checkpoint >= 1:
                            checkpoint()
                        if now - last_sync >= settings["fsync seconds"]:
                            durable()
                            last_sync = now
                    checkpoint()
                    now = time.monotonic()
                    bar(now, force=True)
                    progress(f"batch={batch} k={k} completed trials={batch_completed}/{settings['batch size']} "
                             f"overall trials={completed_total} message failures={aggregate.stats['message failures']['sum']} "
                             f"full block failures={aggregate.stats['full block failures']['sum']} {throughput(now)}")
                    if batch == U64:
                        raise OverflowError("batch identity exhausted; start a new run")
                    batch += 1
        finally:
            checkpoint()
            progress(f"finalizing trials={completed_total} interrupted={bool(lib.product_interrupted())} "
                     f"error={error} {throughput(time.monotonic())}")
            durable()
        if error:
            raise error


def main():
    # Python 3.11's decimal conversion guard is not a statistical counter limit.
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, epilog=(
        "One uniform inclusive k per batch; each block has exactly k distinct flipped bits. "
        "Defaults: RS256,224 x RS256,254, indefinite batches until Ctrl+C. "
        "New trials send the all-zero codeword, bypassing message generation and encoding. "
        "Linear syndromes and BDD delta gates make residual errors codeword-translation invariant. "
        "Replay honors metadata; legacy records without codeword use the private random reference path. "
        "Journal increments aggregate up to checkpoint-trials blocks, flushed at least "
        "once per second subject to ordered in-flight completion. "
        "Crash loss: bounded in-flight window and staged increment plus writes since last fsync (default 5 seconds, "
        "subject to I/O scheduling). Ctrl+C drains in-flight blocks and fsyncs. "
        "Report mode ignores only a final non-newline-terminated journal fragment. "
        "Resume is not supported. Counts and squared sums are exact JSON integers. "
        "Fisher-Yates always records flips.bin (RSFLIP01): little-endian uint32 bit positions, "
        "bit=8*row-major byte index+LSB-first bit index. Dense records store N-k unflips "
        "with a complement flag. Records include batch/trial identity, 22 metrics, SHA256; "
        "journal flip start/end offsets reference data fsynced before journal publication. "
        "At k=2600 recording costs 10640 bytes/trial, plus a 40-byte file header. "
        "--replay verifies all committed trials; unreferenced crash tails are not trials."))
    parser.add_argument("--output", type=Path, help="new run directory; must not already exist")
    parser.add_argument("--report", type=Path, help="regenerate summary in an existing run, without trials")
    parser.add_argument("--replay", type=Path, help="verify every committed saved trial and its 22 metrics")
    parser.add_argument("--sampler", choices=("floyd", "fisher-yates"), default="floyd",
                        help="Fisher-Yates saves flips.bin; seeds alone cannot replay its worker history")
    parser.add_argument("--seed", type=integer, help="uint64 root seed; otherwise generated once and persisted")
    for flag, default, help_text in (
        ("batch-size", 1000, "blocks per random-k batch"),
        ("batches", 0, "batch limit, 0 means indefinite"),
        ("threads", 1, "parallel independent product blocks (1..1024)"),
        ("minimum-flipped-bits", 2500, "inclusive lower bound for batch k"),
        ("maximum-flipped-bits", 2700, "inclusive upper bound for batch k"),
        ("max-directional-passes", 16, "directional pass cap (2..1000000)"),
        ("checkpoint-trials", 64, "maximum trials per journal increment; also caps concurrency (1..4096)"),
        ("report-seconds", 2, "log/non-TTY progress snapshot interval (1..86400 seconds); TTY bar updates every 0.2 seconds"),
        ("fsync-seconds", 5, "durability and summary interval (1..86400 seconds)")):
        parser.add_argument("--" + flag, type=integer, default=default, help=help_text)
    parser.add_argument("--anchors", action=argparse.BooleanOptionalAction, default=True,
                        help="enable strong-success anchor protection")
    parser.add_argument("--binary-image", action=argparse.BooleanOptionalAction, default=True,
                        help="enable weak repair delta popcount <= 2 gate")
    args = parser.parse_args()
    require(sum(map(bool, (args.output, args.report, args.replay))) == 1,
            "specify exactly one of --output, --report or --replay")
    lib = native()
    if args.report or args.replay:
        flag = "--replay" if args.replay else "--report"
        directory = args.replay or args.report
        require((len(sys.argv) == 3 and sys.argv[1] == flag) or
                (len(sys.argv) == 2 and sys.argv[1].startswith(flag + "=")),
                f"{flag} accepts only a directory")
        with (directory / "run.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            recover(directory, lib, replay=bool(args.replay))
        return
    settings = {"root seed": args.seed if args.seed is not None else secrets.randbits(64),
                "batch size": args.batch_size, "batches": args.batches, "threads": args.threads,
                "minimum flipped bits": args.minimum_flipped_bits,
                "maximum flipped bits": args.maximum_flipped_bits,
                "maximum directional passes": args.max_directional_passes,
                "anchors": args.anchors, "binary image": args.binary_image,
                "checkpoint trials": args.checkpoint_trials,
                "report seconds": args.report_seconds, "fsync seconds": args.fsync_seconds}
    validate_settings(settings)
    args.output.mkdir()  # No exist_ok: never overwrite a previous experiment.
    with (args.output / "run.lock").open("x") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        run(args.output, settings, lib, args.sampler)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, OverflowError) as error:
        print(f"{timestamp()} error: {error}", file=sys.stderr)
        sys.exit(1)
