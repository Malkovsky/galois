import copy
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "plot_product_monte_carlo.py"
SPEC = importlib.util.spec_from_file_location("plot_product", SCRIPT)
plot = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(plot)


class NumericalTest(unittest.TestCase):
    def test_binomial_identity_and_weighted_channel(self):
        n, p = 12, 0.17
        rows = {k: {"trials": 1, "full": k} for k in range(n + 1)}
        values, covered, missing = plot.evaluate(rows, p, n, {"full": n})
        self.assertAlmostEqual(math.exp(values["full"]), p, places=14)
        self.assertEqual(covered, 0)
        self.assertEqual(missing, -math.inf)
        for low in range(n + 1):
            for high in range(low, n + 1):
                expected = sum(math.comb(n, k) * p**k * (1 - p)**(n - k)
                               for k in range(low, high + 1))
                self.assertAlmostEqual(plot.log_interval(n, low, high, p), math.log(expected), places=12)

    def test_tiny_exact_endpoint_and_large_counters(self):
        n, p = plot.N, 0.0045
        count = 10**400
        rows = {n: {"trials": count, "full": count * n}}
        values, covered, missing = plot.evaluate(rows, p, n, {"full": n})
        expected = n * math.log(p)
        self.assertAlmostEqual(values["full"], expected, places=8)
        self.assertEqual(covered, expected)
        self.assertEqual(missing, 0)
        self.assertTrue(math.isfinite(values["full"] / math.log(10)))
        self.assertTrue(math.isnan(plot.plot_value(values["full"])))

    def test_nearly_complete_coverage_keeps_tiny_missing_tail(self):
        n, p = 100, 0.0045
        rows = {k: {"trials": 1, "full": 0} for k in range(n)}
        values, covered, missing = plot.evaluate(rows, p, n, {"full": n})
        self.assertAlmostEqual(missing, n * math.log(p), places=10)
        self.assertLess(covered, 0)
        self.assertEqual(values["full"], -math.inf)
        self.assertTrue(math.isnan(plot.plot_value(values["full"])))

    def test_empty_and_sparse_coverage_not_renormalized(self):
        values, covered, missing = plot.evaluate({}, 0.1, 10, {"full": 10})
        self.assertEqual((values["full"], covered, missing), (-math.inf, -math.inf, 0))
        rows = {2: {"trials": 2, "full": 10}}
        values, covered, missing = plot.evaluate(rows, 0.1, 10, {"full": 10})
        self.assertAlmostEqual(values["full"], covered + math.log(0.5))
        self.assertAlmostEqual(math.exp(covered) + math.exp(missing), 1)

    def test_disjoint_missing_intervals(self):
        n = 30
        for p in (0.0045, 0.5, 0.9955):
            for sampled in ({0, 2, 10, 29}, set(range(1, n + 1))):
                rows = {k: {"trials": 1, "full": 1} for k in sampled}
                _, covered, missing = plot.evaluate(rows, p, n, {"full": n})
                expected = sum(math.comb(n, k) * p**k * (1 - p)**(n - k)
                               for k in range(n + 1) if k not in sampled)
                self.assertAlmostEqual(missing, math.log(expected), places=11)
                self.assertAlmostEqual(math.exp(covered) + math.exp(missing), 1)


class ReportTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="plot-product-")
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def fixture(self, name="run", seed=1, trials=1, residual=1, passes=16):
        path = self.root / name
        path.mkdir()
        metadata = {"schema revision": 1, "created at": name, "code": plot.CODE,
                    "random algorithm": plot.RANDOM,
                    "settings": {"root seed": seed, "batch size": 1000, "batches": 0,
                                 "threads": 1, "minimum flipped bits": 0,
                                 "maximum flipped bits": plot.N,
                                 "maximum directional passes": passes, "anchors": True,
                                 "binary image": True, "checkpoint trials": 64,
                                 "report seconds": 2, "fsync seconds": 5}}
        stats = {metric: {"sum": residual * trials, "squared sum": residual**2 * trials}
                 for metric in plot.METRICS.values()}
        row = {"flipped bit count": 2600, "trial count": trials, "statistics": stats}
        summary = {"schema revision": 1, "run identity": plot.identity(metadata),
                   "overall": {"trial count": trials, "statistics": copy.deepcopy(stats)},
                   "by flipped bit count": [row]}
        self.save(path, metadata, summary)
        return path, metadata, summary

    def save(self, path, metadata, summary):
        (path / "metadata.json").write_text(json.dumps(metadata))
        (path / "summary.json").write_text(json.dumps(summary))

    def test_pool_counts_not_means_and_separate_configs(self):
        a, _, _ = self.fixture("a", trials=1, residual=10)
        b, _, _ = self.fixture("b", seed=2, trials=9, residual=0)
        c, _, _ = self.fixture("c", passes=32)
        groups = plot.pool_reports([a, b / "summary.json", c])
        self.assertEqual(len(groups), 2)
        row = groups[(16, True, True)][2600]
        self.assertEqual(row, {"trials": 10, "information": 10, "full": 10})

    def test_duplicate_identity_and_seed(self):
        a, _, _ = self.fixture("a")
        b, _, _ = self.fixture("b")
        with self.assertRaisesRegex(ValueError, "duplicate run identity"):
            plot.pool_reports([a, a / "summary.json"])
        with self.assertRaisesRegex(ValueError, "repeated seed"):
            plot.pool_reports([a, b])

    def test_schema2_pooling_exact_integers_and_reconciliation(self):
        old, _, _ = self.fixture("old", seed=1, trials=1, residual=10)
        path, metadata, summary = self.fixture("new", seed=2)
        n = 10**80
        metadata["schema revision"] = summary["schema revision"] = 2
        summary["run identity"] = plot.identity(metadata)
        stats = {"completed blocks": n, "total iterations": n*2,
                 "information bits": {"total bits": n*455168, "raw corrupted bits": n*2000,
                                      "post decoding corrupted bits": n*3},
                 "full-codeword bits": {"total bits": n*524288, "raw corrupted bits": n*2600,
                                        "post decoding corrupted bits": n*4}}
        summary["overall"] = {"statistics": copy.deepcopy(stats)}
        summary["by flipped bit count"] = [{"flipped bit count": 2600, "statistics": copy.deepcopy(stats)}]
        self.save(path, metadata, summary)
        rows = plot.pool_reports([old, path])[(16, True, True)]
        self.assertEqual(rows[2600], {"trials": n+1, "information": n*3+10, "full": n*4+10})
        for mutation in (
                lambda s: s["overall"]["statistics"].update({"completed blocks": n+1}),
                lambda s: s["by flipped bit count"][0]["statistics"]["information bits"].update({"total bits": 1}),
                lambda s: s["by flipped bit count"][0]["statistics"].update({"total iterations": True}),
                lambda s: s["by flipped bit count"][0]["statistics"]["full-codeword bits"].update({"raw corrupted bits": 0})):
            bad = copy.deepcopy(summary)
            mutation(bad)
            self.save(path, metadata, bad)
            with self.assertRaises(ValueError):
                plot.load_report(path)

    def test_fisher_yates_metadata_preserves_summary_format(self):
        path, metadata, summary = self.fixture()
        expected = plot.load_report(path)[1:]
        metadata["random algorithm"] = plot.RANDOM_FY
        summary["run identity"] = plot.identity(metadata)
        self.save(path, metadata, summary)
        self.assertEqual(plot.load_report(path)[1:], expected)

    def test_codeword_conventions_and_explicit_pooling(self):
        legacy, _, _ = self.fixture("legacy")
        zero, metadata, summary = self.fixture("zero", seed=2)
        metadata["codeword"] = "zero"
        summary["run identity"] = plot.identity(metadata)
        self.save(zero, metadata, summary)
        self.assertEqual(plot.load_report(legacy)[4], "random")
        self.assertEqual(plot.load_report(zero)[4], "zero")
        self.assertEqual(plot.pool_reports([zero])[(16, True, True)][2600]["trials"], 1)
        with self.assertRaisesRegex(ValueError, "mixed zero/random"):
            plot.pool_reports([legacy, zero])
        self.assertEqual(plot.pool_reports([legacy, zero], True)[(16, True, True)][2600]["trials"], 2)
        metadata["settings"]["root seed"] = 1
        summary["run identity"] = plot.identity(metadata)
        self.save(zero, metadata, summary)
        with self.assertRaisesRegex(ValueError, "repeated seed"):
            plot.pool_reports([legacy, zero], True)
        metadata["codeword"] = "unsupported"
        summary["run identity"] = plot.identity(metadata)
        self.save(zero, metadata, summary)
        with self.assertRaisesRegex(ValueError, "codeword convention"):
            plot.load_report(zero)

    def test_flags_separate_groups_and_range_does_not_bypass_seed_check(self):
        a, _, _ = self.fixture("a")
        paths = [a]
        for flag in ("anchors", "binary image"):
            path, metadata, summary = self.fixture(flag)
            metadata["settings"][flag] = False
            summary["run identity"] = plot.identity(metadata)
            self.save(path, metadata, summary)
            paths.append(path)
        self.assertEqual(len(plot.pool_reports(paths)), 3)
        path, metadata, summary = self.fixture("range")
        metadata["settings"]["minimum flipped bits"] = 2500
        summary["run identity"] = plot.identity(metadata)
        self.save(path, metadata, summary)
        with self.assertRaisesRegex(ValueError, "repeated seed"):
            plot.pool_reports([a, path])

    def test_malformed_reports(self):
        path, metadata, summary = self.fixture()
        mutations = [
            lambda s: s.update({"schema revision": True}),
            lambda s: s.update({"run identity": "wrong"}),
            lambda s: s["by flipped bit count"].append(copy.deepcopy(s["by flipped bit count"][0])),
            lambda s: s["by flipped bit count"][0].update({"flipped bit count": plot.N + 1}),
            lambda s: s["by flipped bit count"][0].update({"trial count": True}),
            lambda s: s["overall"].update({"trial count": 2}),
            lambda s: s["overall"]["statistics"][plot.METRICS["full"]].update({"sum": 0}),
            lambda s: s["by flipped bit count"][0]["statistics"][plot.METRICS["full"]].update(
                {"sum": plot.N + 1, "squared sum": (plot.N + 1)**2}),
        ]
        for mutate in mutations:
            with self.subTest(mutate=mutate):
                broken = copy.deepcopy(summary)
                mutate(broken)
                self.save(path, metadata, broken)
                with self.assertRaises(ValueError):
                    plot.load_report(path)
        for key in ("code", "random algorithm", "schema revision"):
            broken = dict(metadata, **{key: "unsupported"})
            self.save(path, broken, summary)
            with self.assertRaisesRegex(ValueError, "incompatible"):
                plot.load_report(path)

    def test_strict_json_numbers_and_duplicate_keys(self):
        path = self.root / "bad.json"
        for text in ('{"x":1,"x":2}', '{"x":1.0}', '{"x":NaN}'):
            path.write_text(text)
            with self.assertRaises(ValueError):
                plot.read_json(path)

    def test_large_json_counters(self):
        path, _, _ = self.fixture(trials=10**400)
        self.assertEqual(plot.load_report(path)[3][2600]["trials"], 10**400)

    def test_help_and_headless_outputs(self):
        result = subprocess.run([sys.executable, str(SCRIPT), "--help"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("No extrapolation", result.stdout)
        if importlib.util.find_spec("matplotlib") is None:
            self.skipTest("matplotlib not installed")
        path, _, _ = self.fixture(residual=0)
        metadata = plot.read_json(path / "metadata.json")
        summary = plot.read_json(path / "summary.json")
        metadata["codeword"] = "zero"
        summary["run identity"] = plot.identity(metadata)
        self.save(path, metadata, summary)
        output = self.root / "plot.svg"
        result = subprocess.run([sys.executable, str(SCRIPT), str(path), "--output", str(output),
                                 "--points", "3", "--metric", "both"],
                                capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("<svg", output.read_text())
        self.assertIn("-inf", output.with_suffix(".csv").read_text())
        self.assertIn("zero-observed strata", result.stderr)


if __name__ == "__main__":
    unittest.main()
