"""Legacy coordinator and C ABI regression tests (test-only reference)."""

import hashlib
import importlib.machinery
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import pty
import re
import select
import shutil
import signal
import subprocess
import struct
import sys
import tempfile
import threading
import time
import unittest
from unittest import mock

CLI = Path(sys.argv.pop(1)).resolve()
loader = importlib.machinery.SourceFileLoader("experiment", str(CLI))
spec = importlib.util.spec_from_loader(loader.name, loader)
experiment = importlib.util.module_from_spec(spec)
loader.exec_module(experiment)


class ExperimentTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="product-monte-carlo-")
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def invoke(self, *args, success=True):
        result = subprocess.run([sys.executable, str(CLI), *map(str, args)],
                                capture_output=True, text=True, timeout=90)
        self.assertEqual(result.returncode == 0, success, result.stderr)
        return result

    def run_case(self, name, *args):
        path = self.root / name
        self.invoke("--output", path, "--seed", "18446744073709551615",
                    "--batches", "2", "--batch-size", "3", "--checkpoint-trials", "3", *args)
        return path

    def read(self, path, name="summary.json"):
        return json.loads((path / name).read_text())

    def scientific(self, path):
        result = self.read(path)
        del result["run identity"]
        return result

    def test_reproducible_threads_and_recovery(self):
        args = ("--minimum-flipped-bits", "2590", "--maximum-flipped-bits", "2610")
        a = self.run_case("a", "--threads", "1", *args)
        b = self.run_case("b", "--threads", "3", *args)
        c = self.run_case("c", "--threads", "3", *args)
        self.assertEqual(self.scientific(a), self.scientific(b))
        self.assertEqual(self.scientific(a), self.scientific(c))
        before = self.read(b)
        self.invoke("--report", b)
        self.assertEqual(before, self.read(b))
        with (b / "journal.jsonl").open("ab") as stream:
            stream.write(b'{"incomplete crash tail')
        self.invoke("--report", b)
        self.assertEqual(before, self.read(b))
        self.assertEqual(before["overall"]["trial count"], 6)
        chunked = [self.run_case(f"chunked-{threads}", "--threads", threads,
                   "--checkpoint-trials", 64, "--batch-size", 17, *args)
                   for threads in (1, 4, 16)]
        self.assertEqual(self.scientific(chunked[0]), self.scientific(chunked[1]))
        self.assertEqual(self.scientific(chunked[0]), self.scientific(chunked[2]))

    def test_exact_channel_oracles_and_flags(self):
        for k in (0, 1, 524287, 524288):
            path = self.run_case(str(k), "--minimum-flipped-bits", k,
                                 "--maximum-flipped-bits", k, "--batches", "1",
                                 "--batch-size", "1", "--no-anchors", "--no-binary-image")
            self.invoke("--report", path)
            summary = self.read(path)["overall"]
            stats = summary["statistics"]
            self.assertEqual(stats["initial full block corrupted bits"], {"sum": k, "squared sum": k*k})
            if k <= 1:
                self.assertEqual(stats["accepted bit changes"]["sum"], k)
                self.assertEqual(stats["accepted byte changes"]["sum"], k)
                self.assertEqual(stats["residual full block bits"]["sum"], 0)
            else:
                self.assertEqual(stats["initial full block corrupted bytes"]["sum"], 65536)
                # The all-ones difference is itself a product codeword. With
                # one bit absent BDD repairs to that wrong valid codeword.
                self.assertEqual(stats["residual full block bits"]["sum"], 524288)
                self.assertEqual(stats["residual information bits"]["sum"], 455168)
                self.assertEqual(stats["accepted bit changes"]["sum"], 524288 - k)
                self.assertEqual(stats["zero syndrome wrong full blocks"]["sum"], 1)
            settings = self.read(path, "metadata.json")["settings"]
            self.assertFalse(settings["anchors"])
            self.assertFalse(settings["binary image"])
        for anchors in ("--anchors", "--no-anchors"):
            for binary in ("--binary-image", "--no-binary-image"):
                path = self.run_case(anchors + binary, anchors, binary,
                    "--minimum-flipped-bits", "0", "--maximum-flipped-bits", "0")
                self.assertEqual(self.read(path)["overall"]["statistics"]["directional passes"],
                                 {"sum": 12, "squared sum": 24})

    def test_native_known_answers_and_worker_reuse(self):
        lib = experiment.native()
        weights = (0, 1, 2590, 2600, 2610, 262143, 262144, 262145, 524287, 524288)

        def trials(_):
            # Reuse both the ctypes output and native thread-local scratch across
            # sparse/dense noise and different decoder gates.
            output = (experiment.ctypes.c_uint64 * len(experiment.METRICS))()
            values = []
            for index, k in enumerate(weights):
                for i in range(len(output)):
                    output[i] = experiment.U64
                status = lib.product_trial(42, 7, index, k, 16,
                                           index % 2, (index // 2) % 2, output)
                self.assertEqual(status, 0)
                self.assertEqual(output[0], k)
                values.extend(output)
            return hashlib.sha256(struct.pack("<" + "Q" * len(values), *values)).hexdigest()

        expected = "4a4257dff8b039967b33aba2a92dc387a8ccc5b05b39a7d91680910c6175c61d"
        self.assertEqual(trials(0), expected)
        with experiment.concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
            self.assertEqual(list(pool.map(trials, range(6))), [expected] * 6)

    def test_native_chunk_prefix_and_replay(self):
        lib = experiment.native()
        c = experiment.ctypes
        lib.product_interrupt_install()
        output = (c.c_uint64 * 88)()
        single = (c.c_uint64 * 22)()
        count = c.c_uint64()
        for k in (0, 2600, 524288):
            stride = min(k, 524288 - k)
            for sampler in (0, 1):
                positions = (c.c_uint32 * (4 * stride))()
                self.assertEqual(lib.product_trials(42, 0, 7, k, 16, 1, 1,
                    output, 4, sampler, positions, c.byref(count)), 0)
                self.assertEqual(count.value, 4)
                for i in range(4):
                    saved = (c.c_uint32 * stride).from_buffer(positions, i * stride * 4)
                    self.assertEqual(lib.product_trial_flips(42, 0, 7 + i, k, 16, 1, 1,
                        single, 2, saved), 0)
                    self.assertEqual(list(single), list(output[i * 22:(i + 1) * 22]))
        before = list(output)
        for first, size, passes in ((0, 5, 16), (experiment.U64, 2, 16), (0, 4, 1)):
            self.assertNotEqual(lib.product_trials(42, 0, first, 2600, passes, 1, 1,
                output, size, 0, None, c.byref(count)), 0)
            self.assertEqual(count.value, 0)
            self.assertEqual(list(output), before)

    def test_chunk_failure_preserves_prefix_and_later_successes(self):
        source = self.run_case("chunk-settings", "--minimum-flipped-bits", "0", "--maximum-flipped-bits", "0")
        settings = self.read(source, "metadata.json")["settings"]
        settings.update({"threads": 3, "checkpoint trials": 12, "batch size": 12, "batches": 1})
        lib = experiment.native()
        barrier = threading.Barrier(3)

        class FailChunk:
            def __getattr__(self, name):
                return getattr(lib, name)

            def product_trials(self, *args):
                barrier.wait(timeout=10)
                if args[2] == 4:
                    args = list(args)
                    args[8] = 1
                    self_status = lib.product_trials(*args)
                    assert self_status == 0
                    return 9
                return lib.product_trials(*args)

        for sampler in ("floyd", "fisher-yates"):
            path = self.root / sampler
            path.mkdir()
            with self.assertRaisesRegex(RuntimeError, "index=5 failed: 9"):
                experiment.run(path, settings, FailChunk(), sampler)
            records = [json.loads(line) for line in (path / "journal.jsonl").read_text().splitlines()]
            indices = [i for r in records for i in range(r["first trial index"], r["past last trial index"])]
            self.assertEqual(indices, list(range(5)) + list(range(8, 12)))
            self.invoke("--replay" if sampler == "fisher-yates" else "--report", path)

    def test_slow_first_trial_bounds_refill_and_orders_journal(self):
        source = self.run_case("window-settings", "--minimum-flipped-bits", "0", "--maximum-flipped-bits", "0")
        settings = self.read(source, "metadata.json")["settings"]
        settings.update({"threads": 3, "checkpoint trials": 6, "batch size": 13, "batches": 1})
        submitted = []
        released_first = False

        def submit(function, index, k):
            if not released_first:
                self.assertLess(index, 6)
            submitted.append(index)
            future = experiment.concurrent.futures.Future()
            future.index = index
            future.set_result(function(index, k))
            return future

        def wait(pending, **kwargs):
            nonlocal released_first
            latest = max(pending, key=lambda f: f.index)
            if not released_first and latest.index < 2:
                self.assertEqual(submitted, list(range(6)))
                latest = min(pending, key=lambda f: f.index)
                released_first = True
            return {latest}, pending - {latest}

        pool = mock.MagicMock()
        pool.__enter__.return_value.submit.side_effect = submit
        path = self.root / "window"
        path.mkdir()
        with mock.patch.object(experiment.concurrent.futures, "ThreadPoolExecutor", return_value=pool), \
                mock.patch.object(experiment.concurrent.futures, "wait", side_effect=wait):
            experiment.run(path, settings, experiment.native())
        self.assertTrue(released_first)
        self.assertEqual(submitted, list(range(13)))
        records = [json.loads(line) for line in (path / "journal.jsonl").read_text().splitlines()]
        self.assertTrue(all(r["trial count"] <= 6 for r in records))
        self.assertEqual([i for r in records for i in range(r["first trial index"], r["past last trial index"])], list(range(13)))
        self.invoke("--report", path)

    def test_invalid_inputs_and_no_overwrite(self):
        for flags in (("--minimum-flipped-bits", "5", "--maximum-flipped-bits", "4"),
                      ("--maximum-flipped-bits", "524289"), ("--seed", "18446744073709551616"),
                      ("--seed", "-1"), ("--threads", "0"), ("--threads", "1025"),
                      ("--max-directional-passes", "1"), ("--batch-size", "0"),
                      ("--checkpoint-trials", "0"), ("--fsync-seconds", "0")):
            self.invoke("--output", self.root / "invalid", *flags, success=False)
            self.assertFalse((self.root / "invalid").exists())
        path = self.run_case("existing", "--minimum-flipped-bits", "0", "--maximum-flipped-bits", "0")
        before = (path / "journal.jsonl").read_bytes()
        self.invoke("--output", path, success=False)
        self.assertEqual(before, (path / "journal.jsonl").read_bytes())

    def test_codeword_translation_equivariance(self):
        lib = experiment.native()
        c = experiment.ctypes
        output = (c.c_uint64 * 22)()
        residual = (c.c_uint8 * 65536)()
        # Fixed bit patterns: strong radius 16 and beyond; weak radius 1 and
        # beyond in columns that strong BDD cannot repair; parity boundaries.
        patterns = [[], [0]]
        patterns += [[8 * (r * 256 + col) for r in range(n) for col in range(width)]
                     for n in (16, 17, 33) for width in (1, 2, 3)]
        patterns += [[8 * (r * 256 + col) + bit for r, col, bit in
                      ((223, 253, 7), (223, 254, 0), (224, 253, 1), (255, 255, 7))]]
        failures = wrong = capped = False
        for anchors in (0, 1):
            for binary in (0, 1):
                for passes in (2, 3, 16):
                    cases = [(len(p), p) for p in patterns]
                    # Fixed Floyd realizations, including dense complement boundary.
                    for k in (2600, 4000, 262143, 262144, 262145, 524287, 524288):
                        positions = (c.c_uint32 * min(k, 524288-k))()
                        self.assertEqual(lib.product_trial_flips(
                            42, 7, 9, k, passes, anchors, binary, output, 0, positions), 0)
                        cases.append((k, positions))
                    for k, positions in cases:
                        positions = (c.c_uint32 * len(positions))(*positions)
                        expected = None
                        # Alternate specializations on the same native thread; do
                        # not let a stale random original contaminate zero trials.
                        for random in (0, 1, 0):
                            self.assertEqual(lib.product_trial_reference(
                                42, 7, 9, k, passes, anchors, binary, output,
                                2, positions, random, residual), 0)
                            actual = (list(output), bytes(residual))
                            if expected is None:
                                expected = actual
                            self.assertEqual(actual, expected,
                                             (anchors, binary, passes, k, random))
                        failures |= bool(output[9] and not output[10])
                        wrong |= bool(output[11])
                        capped |= bool(output[21])
                        if k == 524288:
                            self.assertEqual(bytes(residual), b"\xff" * 65536)
        self.assertTrue(failures)
        self.assertTrue(wrong)  # Dense all-ones valid-codeword miscorrection witness.
        self.assertTrue(capped)

    def test_reference_floyd_positions_and_transactionality(self):
        lib = experiment.native()
        c = experiment.ctypes
        outputs = [(c.c_uint64 * 22)() for _ in range(2)]
        positions = [(c.c_uint32 * 2600)() for _ in range(2)]
        for random in (0, 1):
            self.assertEqual(lib.product_trial_reference(
                experiment.U64, 3, 8, 2600, 16, 1, 1, outputs[random],
                0, positions[random], random, None), 0)
        self.assertEqual(list(outputs[0]), list(outputs[1]))
        self.assertEqual(list(positions[0]), list(positions[1]))
        output = (c.c_uint64 * 22)(*([experiment.U64] * 22))
        residual = (c.c_uint8 * 65536)(*([123] * 65536))
        for random in (0, 1, 2):
            for k, passes, sampler, p in ((524289, 16, 0, None), (0, 1, 0, None),
                    (0, 16, 3, None), (2, 16, 2, (c.c_uint32 * 2)(4, 4)),
                    (1, 16, 2, (c.c_uint32 * 1)(524288))):
                self.assertNotEqual(lib.product_trial_reference(
                    42, 0, 0, k, passes, 1, 1, output, sampler, p, random, residual), 0)
                self.assertEqual(list(output), [experiment.U64] * 22)
                self.assertEqual(bytes(residual), bytes([123]) * 65536)

    def test_legacy_random_replay_without_metadata_mutation(self):
        path = self.run_case("legacy", "--sampler", "fisher-yates",
                             "--minimum-flipped-bits", "2600", "--maximum-flipped-bits", "2600")
        metadata = self.read(path, "metadata.json")
        self.assertEqual(metadata.pop("codeword"), "zero")
        # Translate the identity only, preserving saved positions and counters.
        # Equivariance means these records are valid random-reference fixtures.
        digest = experiment.identity(metadata)
        experiment.atomic_json(path / "metadata.json", metadata)
        journal = [json.loads(line) for line in (path / "journal.jsonl").read_text().splitlines()]
        for record in journal:
            record["run identity"] = digest
        (path / "journal.jsonl").write_text("".join(experiment.canonical(r) + "\n" for r in journal))
        flips = (path / "flips.bin").read_bytes()
        (path / "flips.bin").write_bytes(experiment.FLIP_MAGIC + bytes.fromhex(digest) + flips[40:])
        before = (path / "metadata.json").read_bytes()
        lib = experiment.native()
        with mock.patch.object(lib, "product_trial_reference", wraps=lib.product_trial_reference) as trial:
            experiment.recover(path, lib, replay=True)
            self.assertEqual(trial.call_count, 6)
            self.assertTrue(all(call.args[-2] for call in trial.call_args_list))
        self.assertEqual(before, (path / "metadata.json").read_bytes())
        self.assertEqual(self.read(path)["run identity"], digest)
        self.invoke("--report", path)
        self.assertEqual(before, (path / "metadata.json").read_bytes())

    def test_saved_flips_replay_and_corruption(self):
        for k in (0, 1, 2600, 262144, 524287, 524288):
            path = self.run_case(f"fy-{k}", "--sampler", "fisher-yates", "--threads", "3",
                                 "--minimum-flipped-bits", k, "--maximum-flipped-bits", k)
            before = self.read(path)
            self.invoke("--replay", path)
            self.assertEqual(before, self.read(path))
            expected_size = 40 + 6 * (experiment.FLIP_HEADER.size + 4 * min(k, 524288-k) + 32)
            self.assertEqual((path / "flips.bin").stat().st_size, expected_size)
        data = (path / "flips.bin").read_bytes()
        for bad in (data[:-1], data[:45] + bytes([data[45] ^ 1]) + data[46:], b"bad"):
            (path / "flips.bin").write_bytes(bad)
            self.invoke("--report", path, success=False)
            self.assertEqual(before, self.read(path))
        (path / "flips.bin").write_bytes(data + b"uncommitted partial record")
        self.invoke("--replay", path)
        (path / "flips.bin").unlink()
        self.invoke("--report", path, success=False)

    def test_native_saved_positions_and_persistent_permutation(self):
        lib = experiment.native()
        output = (experiment.ctypes.c_uint64 * len(experiment.METRICS))()
        replay = type(output)()
        first = None
        for k in (2600, 2600, 0, 1, 262143, 262144, 262145, 524287, 524288):
            count = min(k, 524288-k)
            positions = (experiment.ctypes.c_uint32 * count)()
            self.assertEqual(lib.product_trial_flips(42, 7, 0, k, 16, 1, 1, output, 1, positions), 0)
            self.assertEqual(len(set(positions)), count)
            self.assertTrue(all(p < 524288 for p in positions))
            self.assertEqual(output[0], k)
            if k == 2600:
                if first is None:
                    first = list(positions)
                else:
                    self.assertNotEqual(first, list(positions))
            self.assertEqual(lib.product_trial_flips(42, 7, 0, k, 16, 1, 1, replay, 2, positions), 0)
            self.assertEqual(list(output), list(replay))
        positions = (experiment.ctypes.c_uint32 * 2)(3, 3)
        self.assertEqual(lib.product_trial_flips(42, 0, 0, 2, 16, 1, 1, replay, 2, positions), 6)
        positions[1] = 524288
        self.assertEqual(lib.product_trial_flips(42, 0, 0, 2, 16, 1, 1, replay, 2, positions), 6)

    def test_flip_sync_precedes_journal_publication(self):
        source = self.run_case("settings", "--minimum-flipped-bits", "1", "--maximum-flipped-bits", "1")
        settings = self.read(source, "metadata.json")["settings"]
        settings.update({"threads": 3, "checkpoint trials": 3})
        path = self.root / "sync-order"
        path.mkdir()
        real_sync = os.fsync
        synced_end = [40]

        def sync(fd):
            target = Path(os.readlink(f"/proc/self/fd/{fd}"))
            if target == path / "flips.bin":
                # Every already-published reference must have been covered by
                # the PREVIOUS flip sync, not this one.
                for line in (path / "journal.jsonl").read_text().splitlines():
                    self.assertLessEqual(json.loads(line)["flip end"], synced_end[0])
                synced_end[0] = target.stat().st_size
            if target == path / "journal.jsonl":
                for line in target.read_text().splitlines():
                    self.assertLessEqual(json.loads(line)["flip end"], synced_end[0])
            real_sync(fd)

        with mock.patch.object(experiment.os, "fsync", side_effect=sync):
            experiment.run(path, settings, experiment.native(), "fisher-yates")
        self.invoke("--replay", path)

    def test_corrupt_recovery_rejected_without_summary_replacement(self):
        source = self.run_case("source", "--minimum-flipped-bits", "1", "--maximum-flipped-bits", "1",
                               "--checkpoint-trials", "1")
        lines = (source / "journal.jsonl").read_text().splitlines(keepends=True)
        record = json.loads(lines[1])
        record["run identity"] = "wrong run"
        overlap = json.loads(lines[1])
        overlap["first trial index"] = 0
        overlap["past last trial index"] = 1
        moment = json.loads(lines[1])
        moment["statistics"]["initial full block corrupted bits"]["sum"] = 0
        for index, replacement in enumerate((lines[0], "garbage\n", json.dumps(record) + "\n",
                                             json.dumps(overlap) + "\n", json.dumps(moment) + "\n",
                                             '{"schema revision":1,"schema revision":1}\n')):
            path = self.root / f"bad-{index}"
            shutil.copytree(source, path)
            before = (path / "summary.json").read_bytes()
            (path / "journal.jsonl").write_text(lines[0] + replacement + "".join(lines[2:]))
            self.invoke("--report", path, success=False)
            self.assertEqual(before, (path / "summary.json").read_bytes())

    def test_big_integer_moments(self):
        aggregate = experiment.Aggregate("test")
        values = [10**30] * len(experiment.METRICS)
        stats = experiment.stats_for([values, values])
        aggregate.add({"flipped bit count": 7, "trial count": 2, "statistics": stats})
        output = json.loads(experiment.canonical(aggregate.summary()))
        self.assertEqual(output["overall"]["statistics"]["accepted bit changes"]["squared sum"], 2 * 10**60)

    def test_progress_success_counts_and_wall_throughput(self):
        for tty, scenario, step, expected in (
                (True, "success", 1.25, 12), (False, "success", 1.25, 12),
                (True, "success", 0.05, 12),
                (True, "short", 0.05, 1), (True, "short", 0.0, 1),
                (True, "failure", 1.25, 5), (False, "failure", 1.25, 5),
                (True, "interrupt", 1.25, 3), (False, "interrupt", 1.25, 3),
                (True, "submit failure", 1.25, 1)):
            with self.subTest(tty=tty, scenario=scenario, step=step):
                path = self.root / f"{tty}-{scenario}-{step}"
                path.mkdir()
                settings = {"root seed": 42, "batch size": 1 if scenario == "short" else 6,
                            "batches": 2 if scenario == "success" else 1, "threads": 3,
                            "minimum flipped bits": 0, "maximum flipped bits": 0,
                            "maximum directional passes": 16, "anchors": True, "binary image": True,
                            "checkpoint trials": 6, "report seconds": 1, "fsync seconds": 5}
                clock = [0.0]
                interrupted = [False]
                lib = mock.Mock()
                lib.product_interrupt_install.return_value = 0
                lib.product_interrupted.side_effect = lambda: interrupted[0]
                lib.product_batch_k.return_value = 0

                def trial(seed, batch, index, k, passes, anchors, binary, output):
                    if scenario == "failure" and index == 4:
                        return 1
                    return 0

                lib.product_trial.side_effect = trial

                def submit(function, index, k):
                    if scenario == "submit failure" and index == 1:
                        raise RuntimeError("injected submit failure")
                    future = experiment.concurrent.futures.Future()
                    try:
                        future.set_result(function(index, k))
                    except Exception as error:
                        future.set_exception(error)
                    future.index = index
                    return future

                def wait(pending, **kwargs):
                    # Deterministic completions with wall time independent of native execution.
                    clock[0] += step
                    if scenario == "interrupt":
                        interrupted[0] = True
                    done = {min(pending, key=lambda future: future.index)}
                    return done, pending - done

                console = io.StringIO()
                console.isatty = lambda: tty
                pool = mock.MagicMock()
                pool.__enter__.return_value.submit.side_effect = submit
                with mock.patch.object(experiment.sys, "stderr", console), \
                        mock.patch.object(experiment.time, "monotonic", lambda: clock[0]), \
                        mock.patch.object(experiment.concurrent.futures, "ThreadPoolExecutor", return_value=pool), \
                        mock.patch.object(experiment.concurrent.futures, "wait", side_effect=wait):
                    if "failure" in scenario:
                        with self.assertRaises(RuntimeError):
                            experiment.run(path, settings, lib)
                    else:
                        experiment.run(path, settings, lib)

                output = console.getvalue()
                log = (path / "progress.log").read_text()
                self.assertNotIn("\r", log)
                self.assertNotIn("\033", log)
                self.assertEqual(self.read(path)["overall"]["trial count"], expected)
                final = log.splitlines()[-1]
                self.assertIn(f"finalizing trials={expected} interrupted={scenario == 'interrupt'}", final)
                self.assertIn(final + "\n", output)
                batch_count = 6 if scenario == "success" else expected
                self.assertIn(f"completed trials={batch_count}/{settings['batch size']}", log)
                rate = expected / clock[0] if clock[0] else 0.0
                rates = re.search(r"wall blocks/s=([\d.]+) information MiB/s=([\d.]+)", final)
                self.assertIsNotNone(rates)
                for actual, wanted in zip(map(float, rates.groups()), (rate, rate * 56896 / 1048576)):
                    self.assertTrue(math.isfinite(actual))
                    self.assertAlmostEqual(actual, wanted, delta=0.00051)
                snapshots = [line for line in log.splitlines() if "persisted overall trials=" in line]
                if step >= 1:
                    self.assertTrue(snapshots)
                    # In-flight successes appear before the wave is journaled, exactly once.
                    self.assertIn("trials=1/6 overall trials=1 persisted overall trials=0", snapshots[0])
                    self.assertIn("wall blocks/s=0.800", snapshots[0])
                    for line in snapshots:
                        self.assertEqual(line in output, not tty)
                    if scenario == "success":
                        self.assertEqual([int(re.search(r"overall trials=(\d+) persisted", line)[1])
                                          for line in snapshots], list(range(1, 13)))
                if tty:
                    self.assertIn("\r[", output)
                    self.assertIn(f"trials={batch_count}/{settings['batch size']} "
                                  f"{batch_count / settings['batch size']:.1%}", output)
                    self.assertTrue(output.endswith("\n"))
                    if scenario == "short":
                        self.assertEqual(output.count("\r["), 1)
                    if scenario == "success" and step == 0.05:
                        # Two forced batch-end bars plus at most one timed update per 0.2s.
                        self.assertLessEqual(output.count("\r["), 2 + int(clock[0] / 0.2))
                else:
                    self.assertEqual(output, log)

    def test_worker_exception_drains_noncontiguous_successes(self):
        source = self.run_case("source", "--minimum-flipped-bits", "0", "--maximum-flipped-bits", "0")
        settings = self.read(source, "metadata.json")["settings"]
        settings.update({"threads": 3, "batch size": 3, "batches": 1})
        path = self.root / "failure"
        path.mkdir()
        lib = experiment.native()

        class FailOne:
            def __getattr__(self, name):
                return getattr(lib, name)

            def product_trial(self, seed, batch, trial, *args):
                if trial == 1:
                    raise RuntimeError("injected worker failure")
                return lib.product_trial(seed, batch, trial, *args)

        with self.assertRaisesRegex(RuntimeError, "injected worker failure"):
            experiment.run(path, settings, FailOne())
        summary = self.read(path)
        self.assertEqual(summary["overall"]["trial count"], 2)
        records = [json.loads(line) for line in (path / "journal.jsonl").read_text().splitlines()]
        self.assertEqual([(r["first trial index"], r["past last trial index"]) for r in records], [(0, 1), (2, 3)])
        self.invoke("--report", path)
        self.assertEqual(summary, self.read(path))
        recorded = self.root / "failure-recorded"
        recorded.mkdir()

        class FailSaved(FailOne):
            def product_trial_flips(self, seed, batch, trial, *args):
                if trial == 1:
                    raise RuntimeError("injected worker failure")
                return lib.product_trial_flips(seed, batch, trial, *args)

        with self.assertRaisesRegex(RuntimeError, "injected worker failure"):
            experiment.run(recorded, settings, FailSaved(), "fisher-yates")
        self.invoke("--replay", recorded)
        self.assertEqual(self.read(recorded)["overall"]["trial count"], 2)

    def test_weighted_aggregate_and_generated_seed(self):
        path = self.root / "generated"
        self.invoke("--output", path, "--batches", "2", "--batch-size", "2",
                    "--checkpoint-trials", "1",
                    "--minimum-flipped-bits", "1", "--maximum-flipped-bits", "1")
        seed = self.read(path, "metadata.json")["settings"]["root seed"]
        self.assertTrue(0 <= seed <= experiment.U64)
        self.assertIn(f"root seed={seed}", (path / "progress.log").read_text())
        # Retain a complete first batch plus one trial of the next batch.
        lines = (path / "journal.jsonl").read_text().splitlines(keepends=True)
        (path / "journal.jsonl").write_text("".join(lines[:3]))
        self.invoke("--report", path)
        entry = self.read(path)["by flipped bit count"][0]
        self.assertEqual(entry["trial count"], 3)
        self.assertEqual(entry["statistics"]["accepted bit changes"], {"sum": 3, "squared sum": 3})

    def test_interrupt_partial_batch(self):
        path = self.root / "signal"
        with (self.root / "stderr").open("w") as stderr:
            process = subprocess.Popen([sys.executable, str(CLI), "--output", str(path),
                "--sampler", "fisher-yates",
                "--checkpoint-trials", "4096",
                "--seed", "42", "--batch-size", "1000000", "--threads", "3",
                "--minimum-flipped-bits", "0", "--maximum-flipped-bits", "0"],
                stderr=stderr)
            try:
                deadline = time.monotonic() + 30
                journal = path / "journal.jsonl"
                while not journal.exists() or journal.stat().st_size == 0:
                    self.assertIsNone(process.poll())
                    self.assertLess(time.monotonic(), deadline)
                    time.sleep(0.01)
                # A live run cannot be reaggregated concurrently.
                self.invoke("--report", path, success=False)
                process.send_signal(signal.SIGINT)
                self.assertEqual(process.wait(timeout=30), 0)
            finally:
                if process.poll() is None:
                    process.kill()
                    process.wait()
        summary = self.read(path)
        self.assertGreater(summary["overall"]["trial count"], 0)
        self.assertLess(summary["overall"]["trial count"], 1000000)
        self.invoke("--report", path)
        self.assertEqual(summary, self.read(path))
        self.invoke("--replay", path)
        seen = set()
        for line in (path / "journal.jsonl").read_text().splitlines():
            record = json.loads(line)
            for index in range(record["first trial index"], record["past last trial index"]):
                self.assertNotIn(index, seen)
                seen.add(index)
        self.assertEqual(len(seen), summary["overall"]["trial count"])
        self.assertNotIn(b"\r", (path / "progress.log").read_bytes())

    def test_tty_progress_is_in_place_but_log_is_plain(self):
        path = self.root / "tty"
        master, slave = pty.openpty()
        process = subprocess.Popen([sys.executable, str(CLI), "--output", str(path),
            "--seed", "42", "--batch-size", "1000000", "--threads", "2",
            "--minimum-flipped-bits", "0", "--maximum-flipped-bits", "0"], stderr=slave)
        os.close(slave)
        output = b""
        try:
            deadline = time.monotonic() + 30
            while b"\r[" not in output:
                self.assertIsNone(process.poll())
                self.assertLess(time.monotonic(), deadline)
                ready, _, _ = select.select([master], [], [], 0.1)
                if ready:
                    output += os.read(master, 16384)
            process.send_signal(signal.SIGINT)
            self.assertEqual(process.wait(timeout=30), 0)
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()
            os.close(master)
        self.assertNotIn(b"\r", (path / "progress.log").read_bytes())
        self.assertIn("finalizing", (path / "progress.log").read_text())


if __name__ == "__main__":
    unittest.main()
