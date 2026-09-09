"""Black-box tests of the native executable; Python is not a runtime dependency."""
import copy
import fcntl
import hashlib
import json
import os
from pathlib import Path
import pty
import re
import select
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest

CLI = Path(sys.argv.pop(1)).resolve()
REFERENCE = CLI.parent / "product_monte_carlo_reference.py"
FAULT = CLI.parent / "product_monte_carlo_fault_cli"


def canonical(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=True, separators=(",", ":"))


def projected(row):
    s = row["statistics"]
    count = row["trial count"]
    return {"completed blocks": count, "total iterations": s["directional passes"]["sum"],
        "information bits": {"total bits": count * 455168,
            "raw corrupted bits": s["initial information corrupted bits"]["sum"],
            "post decoding corrupted bits": s["residual information bits"]["sum"]},
        "full-codeword bits": {"total bits": count * 524288,
            "raw corrupted bits": s["initial full block corrupted bits"]["sum"],
            "post decoding corrupted bits": s["residual full block bits"]["sum"]}}


class NativeTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="native-mc-")
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def invoke(self, *args, success=True, reference=False, fault=False, env=None):
        command = [sys.executable, str(REFERENCE)] if reference else [str(FAULT if fault else CLI)]
        environment = dict(os.environ) if env is None else dict(env)
        if reference and environment.get("MC_REFERENCE_LD_PRELOAD"):
            environment["LD_PRELOAD"] = environment["MC_REFERENCE_LD_PRELOAD"]
        p = subprocess.run([*command, *map(str, args)], capture_output=True,
                           text=True, timeout=90, env=environment)
        self.assertEqual(p.returncode == 0, success, p.stderr)
        return p

    def run_case(self, name, *args, **kwargs):
        path = self.root / name
        self.invoke("--output", path, "--seed", 42, "--batches", 2,
                    "--batch-size", 17, *args, **kwargs)
        return path

    def read(self, path, name="summary.json"):
        return json.loads((path / name).read_text())

    def records(self, path):
        return [json.loads(s) for s in (path / "journal.jsonl").read_text().splitlines()]

    def count(self, path):
        return self.read(path)["overall"]["statistics"]["completed blocks"]

    def test_floyd_matches_legacy_projection_and_thread_determinism(self):
        old = self.run_case("old", reference=True)
        expected = self.read(old)
        for threads in (1, 4, 16):
            path = self.run_case(str(threads), "--threads", threads)
            result = self.read(path)
            self.assertEqual(result["schema revision"], 2)
            self.assertEqual(result["overall"]["statistics"], projected(expected["overall"]))
            for a, b in zip(result["by flipped bit count"], expected["by flipped bit count"]):
                self.assertEqual(a["flipped bit count"], b["flipped bit count"])
                self.assertEqual(a["statistics"], projected(b))
            self.invoke("--report", path)
            self.assertEqual(result, self.read(path))
            self.assertNotIn("squared sum", (path / "summary.json").read_text())
        repeated = self.run_case("repeat", "--threads", 16)
        self.assertEqual(self.read(repeated)["overall"], result["overall"])

    def test_extremes_gates_and_pass_caps(self):
        for k in (0, 1, 524287, 524288):
            for sampler in ("floyd", "fisher-yates"):
                path = self.run_case(f"{k}-{sampler}", "--minimum-flipped-bits", k,
                    "--maximum-flipped-bits", k, "--sampler", sampler,
                    "--batches", 1, "--batch-size", 1, "--no-anchors", "--no-binary-image")
                s = self.read(path)["overall"]["statistics"]
                self.assertEqual(s["full-codeword bits"]["raw corrupted bits"], k)
                self.assertEqual(s["full-codeword bits"]["post decoding corrupted bits"],
                                 0 if k <= 1 else 524288)
                self.invoke("--replay" if sampler == "fisher-yates" else "--report", path)
        for anchors in ("--anchors", "--no-anchors"):
            for binary in ("--binary-image", "--no-binary-image"):
                args = (anchors, binary, "--max-directional-passes", 2,
                        "--minimum-flipped-bits", 2600, "--maximum-flipped-bits", 2600)
                a = self.run_case(anchors+binary, *args)
                b = self.run_case("old"+anchors+binary, *args, reference=True)
                self.assertEqual(self.read(a)["overall"]["statistics"], projected(self.read(b)["overall"]))

    def test_legacy_report_and_saved_replay_hash_defaults(self):
        for sampler in ("floyd", "fisher-yates"):
            path = self.run_case(sampler, "--sampler", sampler, reference=True)
            expected = self.read(path)
            self.invoke("--report", path)
            self.assertEqual(expected, self.read(path))
            if sampler == "fisher-yates":
                self.invoke("--replay", path)
                metadata = self.read(path, "metadata.json")
                del metadata["codeword"]
                # Legacy random-codeword fixture, identity hashed without a default.
                digest = hashlib.sha256(canonical(metadata).encode("ascii")).hexdigest()
                (path / "metadata.json").write_text(canonical(metadata))
                records = self.records(path)
                for r in records:
                    r["run identity"] = digest
                (path / "journal.jsonl").write_text("".join(canonical(r)+"\n" for r in records))
                data = (path / "flips.bin").read_bytes()
                (path / "flips.bin").write_bytes(data[:8]+bytes.fromhex(digest)+data[40:])
                before = (path / "metadata.json").read_bytes()
                self.invoke("--replay", path)
                self.assertEqual(before, (path / "metadata.json").read_bytes())
                self.assertEqual(self.read(path)["run identity"], digest)

    def test_saved_replay_and_corruption(self):
        for k in (0, 2600, 262144, 524288):
            path = self.run_case(f"k-{k}", "--sampler", "fisher-yates", "--threads", 4,
                "--batch-size", 3, "--minimum-flipped-bits", k, "--maximum-flipped-bits", k)
            before = self.read(path)
            p = self.invoke("--replay", path)
            self.assertIn("all 22 metrics match", p.stderr)
            self.assertEqual(before, self.read(path))
            self.assertEqual((path / "flips.bin").stat().st_size, 40+6*(240+4*min(k,524288-k)))
        data = (path / "flips.bin").read_bytes()
        for bad in (data[:-1], data[:45]+bytes([data[45]^1])+data[46:], b"bad"):
            (path / "flips.bin").write_bytes(bad)
            self.invoke("--report", path, success=False)
            self.assertEqual(before, self.read(path))
        (path / "flips.bin").write_bytes(data+b"uncommitted tail")
        self.invoke("--replay", path)

    def test_uint64_seed_and_canonical_unicode_identity(self):
        for reference in (False, True):
            path = self.run_case(str(reference), "--seed", 2**64-1,
                "--batches", 1, "--batch-size", 1, "--sampler", "fisher-yates",
                "--minimum-flipped-bits", 0, "--maximum-flipped-bits", 0,
                reference=reference)
            metadata = self.read(path, "metadata.json")
            self.assertEqual(metadata["settings"]["root seed"], 2**64-1)
            digest = hashlib.sha256(canonical(metadata).encode("ascii")).hexdigest()
            self.assertEqual(self.read(path)["run identity"], digest)
            self.assertEqual((path / "flips.bin").read_bytes()[8:40], bytes.fromhex(digest))
            metadata["created at"] = "unicode \u00e9 \U0001f600 / \x0f\n"
            digest = hashlib.sha256(canonical(metadata).encode("ascii")).hexdigest()
            # Exercise raw UTF-8 parsing as well as canonical ASCII serialization.
            (path / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False), encoding="utf-8")
            records = self.records(path)
            for record in records:
                record["run identity"] = digest
            (path / "journal.jsonl").write_text("".join(canonical(r)+"\n" for r in records))
            data = (path / "flips.bin").read_bytes()
            (path / "flips.bin").write_bytes(data[:8]+bytes.fromhex(digest)+data[40:])
            self.invoke("--replay", path)
            self.assertEqual(self.read(path)["run identity"], digest)

    def test_oversized_counters_and_legacy_moments_rejected(self):
        for reference in (False, True):
            path = self.run_case(str(reference), "--batches", 1, "--batch-size", 1,
                                 reference=reference)
            before = (path / "summary.json").read_bytes()
            original = self.records(path)[0]
            for value in (2**64, 10**100, 10**400, 1.0, -1):
                record = copy.deepcopy(original)
                if reference:
                    record["statistics"]["accepted bit changes"]["squared sum"] = value
                else:
                    record["statistics"]["total iterations"] = value
                (path / "journal.jsonl").write_text(canonical(record)+"\n")
                p = self.invoke("--report", path, success=False)
                self.assertRegex(p.stderr, "uint64|unsigned integer")
                self.assertEqual(before, (path / "summary.json").read_bytes())

    def test_counter_overflow_preserves_committed_prefix(self):
        for sampler in ("floyd", "fisher-yates"):
            path = self.root / sampler
            p = self.invoke("--output", path, "--seed", 42, "--batches", 1,
                "--batch-size", 3, "--checkpoint-trials", 1, "--sampler", sampler,
                "--minimum-flipped-bits", 0, "--maximum-flipped-bits", 0,
                fault=True, env=dict(os.environ, MC_TEST_COUNTER_OVERFLOW="1"), success=False)
            self.assertIn("uint64 counter addition overflow", p.stderr)
            self.assertEqual(len(self.records(path)), 1)
            self.assertEqual(self.count(path), 0)  # Last durable summary is unchanged.
            self.invoke("--replay" if sampler == "fisher-yates" else "--report", path)
            self.assertEqual(self.count(path), 1)

    def test_recovery_strictness_and_partial_tail(self):
        path = self.run_case("source", "--checkpoint-trials", 1)
        before = (path / "summary.json").read_bytes()
        lines = (path / "journal.jsonl").read_text().splitlines(keepends=True)
        bad_identity = json.loads(lines[1]); bad_identity["run identity"] = "wrong"
        overlap = json.loads(lines[1]); overlap["first trial index"] = 0
        bad_stats = json.loads(lines[1]); bad_stats["statistics"]["full-codeword bits"]["raw corrupted bits"] = 0
        for replacement in (lines[0], "garbage\n", "{}\n", canonical(bad_identity)+"\n",
                canonical(overlap)+"\n", canonical(bad_stats)+"\n",
                '{"schema revision":2,"schema revision":2}\n'):
            (path / "journal.jsonl").write_text(lines[0]+replacement+"".join(lines[2:]))
            self.invoke("--report", path, success=False)
            self.assertEqual(before, (path / "summary.json").read_bytes())
        for tail in ('{"incomplete', 'garbage', '{"complete but no newline":1}'):
            (path / "journal.jsonl").write_text("".join(lines)+tail)
            p = self.invoke("--report", path)
            self.assertIn("incomplete tail=True", p.stderr)
            self.assertEqual(json.loads(before), self.read(path))
        (path / "journal.jsonl").write_text("".join(lines)+"garbage\n")
        self.invoke("--report", path, success=False)

    def test_validation_seed_install_runtime_and_no_overwrite(self):
        for args in (("--seed", -1), ("--seed", 2**64), ("--seed", "1.0"),
                     ("--threads", 0), ("--threads", 1025), ("--batch-size", 0),
                     ("--minimum-flipped-bits", 5, "--maximum-flipped-bits", 4),
                     ("--max-directional-passes", 1), ("--checkpoint-trials", 4097),
                     ("--report-seconds", 0), ("--fsync-seconds", 0), ("--sampler", "unknown")):
            self.invoke("--output", self.root / "invalid", *args, success=False)
            self.assertFalse((self.root / "invalid").exists())
        path = self.root / "generated"
        env = dict(os.environ, PATH="/nonexistent")
        self.invoke("--output", path, "--batches", 1, "--batch-size", 1,
                    "--minimum-flipped-bits", 0, "--maximum-flipped-bits", 0, env=env)
        self.assertEqual(CLI.read_bytes()[:4], b"\x7fELF")
        seed = self.read(path, "metadata.json")["settings"]["root seed"]
        self.assertTrue(0 <= seed < 2**64)
        self.assertIn(f"root seed={seed}", (path / "progress.log").read_text())
        before = (path / "journal.jsonl").read_bytes()
        self.invoke("--output", path, success=False)
        self.assertEqual(before, (path / "journal.jsonl").read_bytes())
        self.invoke("--report", path, "--threads", 1, success=False)
        with (path / "run.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            self.invoke("--report", path, success=False)

    def test_worker_failure_drains_ordered_successes(self):
        for sampler in ("floyd", "fisher-yates"):
            path = self.root / sampler
            env = dict(os.environ, MC_TEST_FAIL="5", MC_TEST_SLOW_FIRST="200")
            p = self.invoke("--output", path, "--batches", 1, "--batch-size", 30,
                "--seed", 42, "--threads", 4, "--checkpoint-trials", 12,
                "--minimum-flipped-bits", 0, "--maximum-flipped-bits", 0,
                "--sampler", sampler, fault=True, env=env, success=False)
            self.assertIn("index=5 failed: injected worker failure", p.stderr)
            indices = [i for r in self.records(path) for i in range(r["first trial index"],r["past last trial index"])]
            self.assertEqual(indices, sorted(set(indices)))
            self.assertEqual(indices[:5], list(range(5)))
            self.assertNotIn(5, indices)
            self.assertTrue(any(i > 5 for i in indices))
            self.assertEqual(self.count(path), len(indices))
            self.invoke("--replay" if sampler == "fisher-yates" else "--report", path)

    def test_native_flip_sync_before_journal_publication(self):
        path = self.run_case("sync-order", "--sampler", "fisher-yates", "--threads", 4,
            "--checkpoint-trials", 3, fault=True, env=dict(os.environ, MC_TEST_SYNC_ORDER="1"))
        self.assertEqual(self.count(path), 34)
        self.invoke("--replay", path)

    def test_signal_drains_bounded_window_and_plain_progress(self):
        for sig in (signal.SIGINT, signal.SIGTERM):
            path = self.root / str(sig)
            err = self.root / f"stderr-{sig}"
            # Slow first block holds the ordered window; other workers must not
            # progress past 12 starts, and SIGTERM must still finish block zero.
            env = dict(os.environ, MC_TEST_SLOW_FIRST="1200")
            with err.open("w") as stream:
                p = subprocess.Popen([str(FAULT), "--output", str(path), "--seed", "42",
                    "--threads", "4", "--batch-size", "1000000", "--checkpoint-trials", "12",
                    "--minimum-flipped-bits", "0", "--maximum-flipped-bits", "0",
                    "--sampler", "fisher-yates", "--report-seconds", "1"], stderr=stream, env=env)
                try:
                    deadline = time.monotonic()+20
                    while not (path / "progress.log").exists() or "batch=0" not in (path / "progress.log").read_text():
                        self.assertIsNone(p.poll()); self.assertLess(time.monotonic(), deadline); time.sleep(.01)
                    time.sleep(.2)
                    self.invoke("--report", path, success=False)
                    p.send_signal(sig)
                    self.assertEqual(p.wait(timeout=30), 0)
                finally:
                    if p.poll() is None: p.kill(); p.wait()
            self.assertEqual(self.count(path), 12)
            indices = [i for r in self.records(path) for i in range(r["first trial index"],r["past last trial index"])]
            self.assertEqual(indices, list(range(12)))
            self.invoke("--replay", path)
            log = (path / "progress.log").read_text()
            self.assertEqual(log, err.read_text())
            self.assertNotIn("\x1b", log)
            self.assertIn("interrupted=True", log)
            final = log.splitlines()[-1]
            match = re.search(r"wall blocks/s=([\d.]+) information MiB/s=([\d.]+) elapsed seconds=([\d.]+)", final)
            self.assertIsNotNone(match)
            rate, mib, elapsed = map(float, match.groups())
            self.assertAlmostEqual(rate, 12/elapsed, delta=.02)
            self.assertAlmostEqual(mib, rate*56896/1048576, delta=.001)

    def test_tty_throttle_and_plain_log(self):
        path = self.root / "tty"
        master, slave = pty.openpty()
        p = subprocess.Popen([str(CLI), "--output", str(path), "--seed", "42",
            "--threads", "4", "--batch-size", "1000000"], stderr=slave)
        os.close(slave)
        output = b""
        start = time.monotonic()
        try:
            while output.count(b"\r[") < 3:
                self.assertIsNone(p.poll()); self.assertLess(time.monotonic()-start, 20)
                if select.select([master], [], [], .1)[0]: output += os.read(master, 65536)
            p.send_signal(signal.SIGINT)
            self.assertEqual(p.wait(timeout=30), 0)
        finally:
            if p.poll() is None: p.kill(); p.wait()
            os.close(master)
        self.assertLessEqual(output.count(b"\r["), int((time.monotonic()-start)/.2)+1)
        log = (path / "progress.log").read_bytes()
        self.assertNotIn(b"\r", log); self.assertNotIn(b"\x1b", log)
        self.assertIn(b"finalizing", log)


if __name__ == "__main__":
    unittest.main()
