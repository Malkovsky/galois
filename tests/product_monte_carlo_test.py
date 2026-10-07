"""Black-box tests of the native executable; Python is not a runtime dependency."""
import copy
import ctypes
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
import struct
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
        "stall patterns corrected": 0, "strong miscorrections detected": 0,
        "strong miscorrections corrected": 0,
        "information bits": {"total bits": count * 455168,
            "raw corrupted bits": s["initial information corrupted bits"]["sum"],
            "post decoding corrupted bits": s["residual information bits"]["sum"]},
        "full-codeword bits": {"total bits": count * 524288,
            "raw corrupted bits": s["initial full block corrupted bits"]["sum"],
            "post decoding corrupted bits": s["residual full block bits"]["sum"]}}


class NativeTest(unittest.TestCase):
    def test_postprocessing_supplemental_abi_known_stall_and_miscorrection(self):
        preload = os.environ.get("MC_REFERENCE_LD_PRELOAD")
        if preload and os.environ.get("LD_PRELOAD") != preload:
            child = subprocess.run([sys.executable, __file__, str(CLI),
                "NativeTest.test_postprocessing_supplemental_abi_known_stall_and_miscorrection"],
                env={**os.environ, "LD_PRELOAD": preload}, capture_output=True, text=True, timeout=90)
            self.assertEqual(child.returncode, 0, child.stderr)
            return
        lib = ctypes.CDLL(str(CLI.parent / "product_monte_carlo_native.so"))
        trial = lib.product_trial_postprocessing
        u64, u32 = ctypes.c_uint64, ctypes.c_uint32
        trial.argtypes = [u64] * 5 + [ctypes.c_int] * 2 + [ctypes.POINTER(u64),
            ctypes.c_int, ctypes.POINTER(u32), ctypes.c_int, ctypes.POINTER(ctypes.c_uint8)] + [u64] * 4 + [ctypes.c_int, ctypes.POINTER(u64)]
        for hidden in (False, True):
            indices = [8 * (row * 256 + col) + bit
                       for row in range(256 if hidden else 17)
                       for col in ([0] if hidden else [254, 255]) for bit in range(8)]
            positions = (u32 * len(indices))(*indices)
            metrics, supplemental = (u64 * 22)(), (u64 * 3)()
            residual = (ctypes.c_uint8 * 65536)()
            self.assertEqual(trial(42, 0, 0, len(indices), 2, 1, 1, metrics, 2,
                positions, 0, residual, 256, 224, 256, 254, 1, supplemental), 0)
            self.assertEqual(list(supplemental), [0, 1, 1] if hidden else [1, 0, 0])
            self.assertFalse(any(residual))
            self.assertEqual(metrics[13], len(indices))
            self.assertEqual(metrics[13], metrics[15] + metrics[17])
            self.assertEqual(metrics[14], metrics[16] + metrics[18])
            self.assertEqual(metrics[12], 2)

    def test_postprocessing_saved_replay_metadata_and_seed_reproducibility(self):
        extra = ("--postprocessing", "--n1", 8, "--k1", 4, "--n2", 5, "--k2", 3,
                 "--minimum-flipped-bits", 20, "--maximum-flipped-bits", 20)
        first = self.run_case("pp-floyd", *extra)
        second = self.run_case("pp-repeat", *extra, "--threads", 4)
        self.assertEqual(self.read(first)["overall"], self.read(second)["overall"])
        self.assertEqual(self.read(first)["by flipped bit count"], self.read(second)["by flipped bit count"])
        self.assertTrue(self.read(first, "metadata.json")["settings"]["postprocessing"])
        self.invoke("--report", first)
        saved = self.run_case("pp-saved", *extra, "--sampler", "fisher-yates")
        before = (saved / "summary.json").read_bytes()
        self.assertEqual((saved / "flips.bin").read_bytes()[:8], b"RSFLIP02")
        self.invoke("--report", saved)
        replay = self.invoke("--replay", saved)
        self.assertIn("all 3 postprocessing counters match", replay.stderr)
        self.assertEqual((saved / "summary.json").read_bytes(), before)
        stats = self.read(saved)["overall"]["statistics"]
        for name in ("stall patterns corrected", "strong miscorrections detected", "strong miscorrections corrected"):
            self.assertEqual(stats[name], sum(row["statistics"][name] for row in self.read(saved)["by flipped bit count"]))
        self.assertGreater(stats["stall patterns corrected"] + stats["strong miscorrections detected"], 0)
        self.invoke("--output", self.root / "badpp", "--postprocessing=false", success=False)

        # Alter only a supplemental counter, recompute checksum and matching
        # summary totals: report is structurally valid, replay must still fail.
        binary = bytearray((saved / "flips.bin").read_bytes())
        offset = 40
        count = struct.unpack_from("<I", binary, offset + 20)[0]
        previous = struct.unpack_from("<Q", binary, offset + 208)[0]
        replacement = 1 - previous
        struct.pack_into("<Q", binary, offset + 208, replacement)
        end = offset + 232 + 4 * count
        binary[end:end + 32] = hashlib.sha256(binary[offset:end]).digest()
        (saved / "flips.bin").write_bytes(binary)
        summary = self.read(saved)
        for row in [summary["overall"], *summary["by flipped bit count"]]:
            row["statistics"]["stall patterns corrected"] += replacement - previous
        (saved / "summary.json").write_text(canonical(summary))
        self.invoke("--report", saved)
        self.assertIn("replay mismatch", self.invoke("--replay", saved, success=False).stderr)

    def test_old_snapshot_missing_supplemental_counters_is_unchanged(self):
        path = self.run_case("old-minimal", "--sampler", "fisher-yates")
        summary = self.read(path)
        for row in [summary["overall"], *summary["by flipped bit count"]]:
            for name in ("stall patterns corrected", "strong miscorrections detected", "strong miscorrections corrected"):
                row["statistics"].pop(name)
        (path / "summary.json").write_text(canonical(summary))
        before = (path / "summary.json").read_bytes()
        self.invoke("--report", path)
        self.invoke("--replay", path)
        self.assertEqual((path / "summary.json").read_bytes(), before)

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
        if not kwargs.get("reference"):
            self.assertFalse((path / "journal.jsonl").exists())
            self.assertEqual(self.read(path, "metadata.json")["storage"], "atomic summary v1")
        return path

    def read(self, path, name="summary.json"):
        return json.loads((path / name).read_text())

    def records(self, path):
        return [json.loads(s) for s in (path / "journal.jsonl").read_text().splitlines()]

    def saved_indices(self, path):
        data = (path / "flips.bin").read_bytes()
        end = self.read(path)["checkpoint"]["flip end"]
        result, offset = [], 40
        while offset < end:
            batch, trial, k, count = struct.unpack_from("<QQII", data, offset)
            result.append(trial)
            offset += 240 + count*4
        self.assertEqual(offset, end)
        return result

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

    def test_shortened_dimensions_threads_report_replay_and_totals(self):
        common = ("--n1", 256, "--k1", 224, "--n2", 175, "--k2", 173,
                  "--minimum-flipped-bits", 1800, "--maximum-flipped-bits", 1810,
                  "--batches", 1, "--batch-size", 6)
        results = []
        for threads in (1, 3):
            path = self.run_case(f"short-{threads}", *common, "--threads", threads)
            result = self.read(path)
            self.assertEqual(result["code parameters"], {"n1": 256, "k1": 224, "n2": 175, "k2": 173})
            results.append(result["overall"])
            stats = result["overall"]["statistics"]
            self.assertEqual(stats["full-codeword bits"]["total bits"], 6*8*256*175)
            self.assertEqual(stats["information bits"]["total bits"], 6*8*224*173)
            metadata = self.read(path, "metadata.json")
            self.assertEqual([metadata["settings"][key] for key in ("n1", "k1", "n2", "k2")],
                             [256, 224, 175, 173])
            before = (path / "metadata.json").read_bytes()
            self.invoke("--report", path)
            self.assertEqual(result, self.read(path))
            self.assertEqual(before, (path / "metadata.json").read_bytes())
        self.assertEqual(results[0], results[1])
        for dims in ((256, 224, 175, 173), (4, 2, 5, 3)):
            n = 8*dims[0]*dims[2]
            for k in (0, 1, n//2, n-1, n):
                path = self.run_case(f"fy-{n}-{k}", "--sampler", "fisher-yates", "--threads", 3,
                    "--batches", 1, "--batch-size", 3,
                    *[arg for key, value in zip(("--n1", "--k1", "--n2", "--k2"), dims) for arg in (key, value)],
                    "--minimum-flipped-bits", k, "--maximum-flipped-bits", k)
                before = self.read(path)
                self.invoke("--replay", path)
                self.assertEqual(before, self.read(path))
                self.assertEqual((path / "flips.bin").stat().st_size, 40+3*(240+4*min(k,n-k)))
                self.assertEqual(before["overall"]["statistics"]["full-codeword bits"]["raw corrupted bits"], 3*k)

    def test_r4_dimensions_threads_report_and_replay(self):
        results = []
        for threads in (1, 3):
            path = self.run_case(f"r4-{threads}", "--n2", 256, "--k2", 252,
                                 "--threads", threads, "--batches", 1, "--batch-size", 6,
                                 "--sampler", "fisher-yates")
            result = self.read(path)
            self.assertEqual(result["code parameters"],
                             {"n1": 256, "k1": 224, "n2": 256, "k2": 252})
            stats = result["overall"]["statistics"]
            self.assertEqual(stats["information bits"]["total bits"], 6*8*224*252)
            self.assertEqual(stats["full-codeword bits"]["total bits"], 6*8*256*256)
            before = (path / "metadata.json").read_bytes()
            self.invoke("--report", path)
            self.invoke("--replay", path)
            self.assertEqual(before, (path / "metadata.json").read_bytes())
            self.assertEqual(result, self.read(path))
            results.append(result["overall"])
        # Fisher-Yates intentionally retains worker-local permutations, so
        # schedules may differ; saved-position replay above is the invariant.
        results = []
        for threads in (1, 3):
            path = self.run_case(f"r4-floyd-{threads}", "--n2", 256, "--k2", 252,
                                 "--threads", threads, "--batches", 1, "--batch-size", 6)
            results.append(self.read(path)["overall"])
        self.assertEqual(results[0], results[1])

    def test_dimension_and_dynamic_range_validation(self):
        for args in (("--n1", 255), ("--k1", 223), ("--n2", 175),
                     ("--n2", 175, "--k2", 171), ("--n2", 8, "--k2", 4),
                     ("--n2", 257, "--k2", 255), ("--n2", 3, "--k2", 1),
                     ("--n2", 175, "--k2", 173, "--maximum-flipped-bits", 358401),
                     ("--n1", 4, "--k1", 2, "--n2", 5, "--k2", 3)):
            self.invoke("--output", self.root / "invalid", *args, success=False)
            self.assertFalse((self.root / "invalid").exists())
        path = self.run_case("bad-range", "--n2", 175, "--k2", 173,
                             "--batches", 1, "--batch-size", 1)
        metadata = self.read(path, "metadata.json")
        metadata["settings"]["maximum flipped bits"] = 358401
        (path / "metadata.json").write_text(canonical(metadata))
        before = (path / "summary.json").read_bytes()
        self.invoke("--report", path, success=False)
        self.assertEqual(before, (path / "summary.json").read_bytes())

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

    def test_legacy_schema2_journal_regeneration_and_replay(self):
        for sampler in ("floyd", "fisher-yates"):
            path = self.run_case(f"legacy2-{sampler}", "--sampler", sampler,
                                 "--batches", 1, "--batch-size", 2, reference=True)
            metadata = self.read(path, "metadata.json")
            metadata["schema revision"] = 2
            digest = hashlib.sha256(canonical(metadata).encode("ascii")).hexdigest()
            (path / "metadata.json").write_text(canonical(metadata))
            records = self.records(path)
            for record in records:
                record["statistics"] = projected(record)
                record["schema revision"] = 2
                record["run identity"] = digest
            (path / "journal.jsonl").write_text("".join(canonical(r)+"\n" for r in records))
            if sampler == "fisher-yates":
                data = (path / "flips.bin").read_bytes()
                (path / "flips.bin").write_bytes(data[:8]+bytes.fromhex(digest)+data[40:])
            expected = projected(self.read(path)["overall"])
            (path / "summary.json").unlink()
            before = (path / "metadata.json").read_bytes()
            self.invoke("--replay" if sampler == "fisher-yates" else "--report", path)
            self.assertEqual(self.read(path)["overall"]["statistics"], expected)
            self.assertEqual((path / "metadata.json").read_bytes(), before)

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
            if reference:
                records = self.records(path)
                for record in records:
                    record["run identity"] = digest
                (path / "journal.jsonl").write_text("".join(canonical(r)+"\n" for r in records))
            else:
                summary = self.read(path)
                summary["run identity"] = digest
                (path / "summary.json").write_text(canonical(summary))
            data = (path / "flips.bin").read_bytes()
            (path / "flips.bin").write_bytes(data[:8]+bytes.fromhex(digest)+data[40:])
            self.invoke("--replay", path)
            self.assertEqual(self.read(path)["run identity"], digest)

    def test_oversized_counters_and_legacy_moments_rejected(self):
        for reference in (False, True):
            path = self.run_case(str(reference), "--batches", 1, "--batch-size", 1,
                                 reference=reference)
            before = (path / "summary.json").read_bytes()
            original = self.records(path)[0] if reference else self.read(path)
            for value in (2**64, 10**100, 10**400, 1.0, -1):
                record = copy.deepcopy(original)
                if reference:
                    record["statistics"]["accepted bit changes"]["squared sum"] = value
                else:
                    record["overall"]["statistics"]["total iterations"] = value
                target = path / ("journal.jsonl" if reference else "summary.json")
                target.write_text(canonical(record)+"\n")
                corrupt = target.read_bytes()
                p = self.invoke("--report", path, success=False)
                self.assertRegex(p.stderr, "uint64|unsigned integer")
                self.assertEqual(before if reference else corrupt, (path / "summary.json").read_bytes())

    def test_counter_overflow_preserves_committed_prefix(self):
        for sampler in ("floyd", "fisher-yates"):
            path = self.root / sampler
            p = self.invoke("--output", path, "--seed", 42, "--batches", 1,
                "--batch-size", 3, "--checkpoint-trials", 1, "--sampler", sampler,
                "--minimum-flipped-bits", 0, "--maximum-flipped-bits", 0,
                fault=True, env=dict(os.environ, MC_TEST_COUNTER_OVERFLOW="1"), success=False)
            self.assertIn("uint64 counter addition overflow", p.stderr)
            self.assertFalse((path / "journal.jsonl").exists())
            self.assertEqual(self.count(path), 0)  # Last durable summary is unchanged.
            self.invoke("--replay" if sampler == "fisher-yates" else "--report", path)
            self.assertEqual(self.count(path), 0)  # Uncommitted work is not recovered.

    def test_atomic_snapshot_failure_preserves_authoritative_summary(self):
        for sampler in ("floyd", "fisher-yates"):
            path = self.root / sampler
            p = self.invoke("--output", path, "--seed", 42, "--batches", 1,
                "--batch-size", 3, "--sampler", sampler, "--n2", 175, "--k2", 173,
                "--minimum-flipped-bits", 1, "--maximum-flipped-bits", 1,
                fault=True, env=dict(os.environ, MC_TEST_SNAPSHOT_FAIL="1"), success=False)
            self.assertIn("before atomic summary rename", p.stderr)
            self.assertFalse((path / "journal.jsonl").exists())
            self.assertEqual(self.count(path), 0)
            self.assertEqual(self.read(path, "summary.json.tmp")["overall"]["statistics"]["completed blocks"], 3)
            before = (path / "summary.json").read_bytes()
            self.invoke("--replay" if sampler == "fisher-yates" else "--report", path)
            self.assertEqual((path / "summary.json").read_bytes(), before)

    def test_snapshot_report_rejects_corruption_without_rewriting(self):
        for sampler in ("floyd", "fisher-yates"):
            path = self.run_case(sampler, "--batches", 1, "--batch-size", 2,
                                 "--sampler", sampler, "--n2", 175, "--k2", 173)
            original = self.read(path)
            mutations = [lambda s: s.update({"run identity": "wrong"}),
                lambda s: s["overall"]["statistics"].update({"completed blocks": 0}),
                lambda s: s["by flipped bit count"].append(copy.deepcopy(s["by flipped bit count"][0])),
                lambda s: s["code parameters"].update(k2=172)]
            if sampler == "fisher-yates":
                mutations.extend([lambda s: s["checkpoint"].update({"flip end": 40}),
                                  lambda s: s["checkpoint"].update({"flip end": 41}),
                                  lambda s: s["checkpoint"].update({"flip end": 2**64-1})])
            for mutation in mutations:
                broken = copy.deepcopy(original)
                mutation(broken)
                (path / "summary.json").write_text(canonical(broken))
                before = (path / "summary.json").read_bytes()
                self.invoke("--report", path, success=False)
                self.assertEqual((path / "summary.json").read_bytes(), before)
            (path / "summary.json").unlink()
            self.invoke("--report", path, success=False)
            self.assertFalse((path / "summary.json").exists())

    def test_periodic_snapshot_survives_abrupt_stop(self):
        path = self.root / "crash"
        p = subprocess.Popen([str(CLI), "--output", str(path), "--seed", "42",
            "--threads", "1", "--batch-size", "1000000", "--n2", "175", "--k2", "173",
            "--minimum-flipped-bits", "1800", "--maximum-flipped-bits", "1800",
            "--sampler", "fisher-yates", "--fsync-seconds", "1"], stderr=subprocess.DEVNULL)
        try:
            deadline = time.monotonic()+15
            while not (path / "summary.json").exists() or self.count(path) == 0:
                self.assertIsNone(p.poll())
                self.assertLess(time.monotonic(), deadline)
                time.sleep(.02)
            p.kill()
            p.wait(timeout=10)
        finally:
            if p.poll() is None:
                p.kill(); p.wait()
        self.assertFalse((path / "journal.jsonl").exists())
        before = (path / "summary.json").read_bytes()
        self.assertGreater(self.count(path), 0)
        self.invoke("--replay", path)
        self.assertEqual((path / "summary.json").read_bytes(), before)

    def test_recovery_strictness_and_partial_tail(self):
        path = self.run_case("source", "--checkpoint-trials", 1, reference=True)
        before = (path / "summary.json").read_bytes()
        lines = (path / "journal.jsonl").read_text().splitlines(keepends=True)
        bad_identity = json.loads(lines[1]); bad_identity["run identity"] = "wrong"
        overlap = json.loads(lines[1]); overlap["first trial index"] = 0
        bad_stats = json.loads(lines[1]); bad_stats["statistics"]["initial full block corrupted bits"]["sum"] = 0
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
        generated = self.root / "generated-again"
        self.invoke("--output", generated, "--batches", 1, "--batch-size", 1,
                    "--minimum-flipped-bits", 0, "--maximum-flipped-bits", 0, env=env)
        other_seed = self.read(generated, "metadata.json")["settings"]["root seed"]
        self.assertTrue(0 <= other_seed < 2**64)
        self.assertNotEqual(seed, other_seed)
        self.invoke("--report", path)
        self.assertEqual(seed, self.read(path, "metadata.json")["settings"]["root seed"])
        before = (path / "summary.json").read_bytes()
        self.invoke("--output", path, success=False)
        self.assertEqual(before, (path / "summary.json").read_bytes())
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
            self.assertFalse((path / "journal.jsonl").exists())
            self.assertGreater(self.count(path), 5)
            if sampler == "fisher-yates":
                indices = self.saved_indices(path)
                self.assertEqual(indices, sorted(set(indices)))
                self.assertEqual(indices[:5], list(range(5)))
                self.assertNotIn(5, indices)
                self.assertTrue(any(i > 5 for i in indices))
                self.assertEqual(self.count(path), len(indices))
            self.invoke("--replay" if sampler == "fisher-yates" else "--report", path)

    def test_native_flip_sync_before_snapshot_publication(self):
        path = self.run_case("sync-order", "--sampler", "fisher-yates", "--threads", 4,
            "--checkpoint-trials", 3, fault=True, env=dict(os.environ, MC_TEST_SYNC_ORDER="1"))
        self.assertEqual(self.count(path), 34)
        self.invoke("--replay", path)

    def test_signal_drains_bounded_window_and_plain_progress(self):
        for sig in (signal.SIGINT, signal.SIGTERM):
            dimensions = ["--n2", "175", "--k2", "173"] if sig == signal.SIGTERM else []
            path = self.root / str(sig)
            err = self.root / f"stderr-{sig}"
            # Slow first block holds the ordered window; other workers must not
            # progress past 12 starts, and SIGTERM must still finish block zero.
            env = dict(os.environ, MC_TEST_SLOW_FIRST="1200")
            with err.open("w") as stream:
                p = subprocess.Popen([str(FAULT), "--output", str(path), "--seed", "42",
                    "--threads", "4", "--batch-size", "1000000", "--checkpoint-trials", "12",
                    "--minimum-flipped-bits", "0", "--maximum-flipped-bits", "0",
                    "--sampler", "fisher-yates", "--report-seconds", "1", *dimensions], stderr=stream, env=env)
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
            self.assertFalse((path / "journal.jsonl").exists())
            indices = self.saved_indices(path)
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
            self.assertAlmostEqual(mib, rate*224*(173 if dimensions else 254)/1048576, delta=.001)

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
