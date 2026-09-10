import copy
import csv
import importlib.util
import io
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "plot_product_monte_carlo.py"
SPEC = importlib.util.spec_from_file_location("plot_product", SCRIPT)
plot = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(plot)


class NumericalTest(unittest.TestCase):
    def test_dimension_specific_conditional_export_and_weighted_plot(self):
        dims = (4, 2, 5, 3)
        ds = plot.denominators(dims)
        n = ds["full"]
        config = (16, True, True) + dims
        rows = {k: {"trials": 2, "information": min(k, ds["information"])*2,
                    "full": 2*k} for k in range(n+1)}
        p = 0.17
        expected = sum(math.comb(n, k)*p**k*(1-p)**(n-k)*min(k, ds["information"])
                       / ds["information"] for k in rows)
        values, covered, missing = plot.evaluate(rows, p, n, ds)
        self.assertAlmostEqual(math.exp(values["full"]), p, places=13)
        self.assertAlmostEqual(math.exp(values["information"]), expected, places=13)
        self.assertEqual((covered, missing), (0, -math.inf))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "export.json"
            plot.export_conditional({config: rows}, ["information", "full"], path)
            exports = json.loads(path.read_text())
            self.assertIn("RS4,2 x RS5,3", exports[0]["configuration"])
            self.assertEqual(exports[0]["points"][1]["raw_ber"], 1/n)
            self.assertAlmostEqual(exports[0]["points"][1]["residual_ber"], 1/ds["information"])
            plt = mock.MagicMock()
            plt.subplots.return_value = (mock.MagicMock(), mock.MagicMock())
            plt.rcParams.__getitem__.return_value.by_key.return_value = {"color": ["blue"]}
            with mock.patch.object(plot, "evaluate", wraps=plot.evaluate) as evaluate, \
                    mock.patch("sys.stderr", io.StringIO()):
                plot.plot_results({config: rows}, ["full"], "both", [p],
                                  Path(directory)/"out.svg", Path(directory)/"out.csv", plt)
            evaluate.assert_called_once_with(rows, p, n, ds)

    def test_pure_conditional_export(self):
        groups = {(16, True, True): {
            1: {"trials": 10, "information": 0, "full": 0},
            3500: {"trials": 20, "information": 4, "full": 8}}}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "pure.json"
            plot.export_conditional(groups, ["information", "full"], path)
            with path.open() as stream:
                rows = json.load(stream, parse_constant=lambda value: self.fail(value))
        self.assertEqual(len(rows), 2)
        for row, metric in zip(rows, ("information", "full")):
            self.assertEqual(set(row), {"configuration", "metric", "points"})
            self.assertEqual(row["configuration"], plot.config_label((16, True, True)))
            self.assertEqual(row["metric"], metric)
            self.assertEqual(len(row["points"]), 2)
            for point in row["points"]:
                self.assertEqual(set(point), {"raw_ber", "residual_ber", "log10_residual_ber"})
        points = rows[0]["points"]
        self.assertEqual(points[0]["residual_ber"], 0)
        self.assertIsNone(points[0]["log10_residual_ber"])
        self.assertEqual(points[0]["raw_ber"], 1 / plot.N)
        self.assertAlmostEqual(points[1]["residual_ber"], 4 / (20 * 455168))
        self.assertAlmostEqual(rows[1]["points"][1]["residual_ber"], 8 / (20 * 524288))

    def test_thin_rejects_invalid_values_before_reading_inputs(self):
        for value in ("0", "-1", "nan", "inf", "-inf"):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "--thin"):
                plot.main(["nonexistent-input", "--thin=" + value])

    def test_thin_scales_markers_and_lines(self):
        for thin in (0.25, 0.5, 1.0):
            with self.subTest(thin=thin), tempfile.TemporaryDirectory() as directory:
                plt = mock.MagicMock()
                figure, axis = mock.MagicMock(), mock.MagicMock()
                plt.subplots.return_value = (figure, axis)
                plt.rcParams.__getitem__.return_value.by_key.return_value = {"color": ["blue"]}
                groups = {(16, True, True): {3500: {"trials": 10, "information": 20, "full": 30}}}
                with mock.patch.object(plot, "evaluate", return_value=(
                        {"information": -10, "full": -9}, -1, -1)), mock.patch("sys.stderr", io.StringIO()):
                    plot.plot_results(groups, ["information"], "both", [0.006],
                                      Path(directory) / "out.svg", Path(directory) / "out.csv", plt,
                                      thin=thin)
                self.assertEqual(axis.scatter.call_args.kwargs["s"], 18 * thin**2)
                self.assertEqual(axis.scatter.call_args.kwargs["linewidths"], thin)
                self.assertEqual(axis.plot.call_args.kwargs["linewidth"], 2 * thin)

    def test_conditional_large_counters_and_zero(self):
        count = 10**400
        rows = {3: {"trials": count, "full": 0},
                1: {"trials": count, "full": count * 3},
                2: {"trials": count, "full": 1}}
        values = list(plot.conditional_values(rows, "full"))
        for value, k, total, mean, log_mean in zip(
                values, (1, 2, 3), (count * 3, 1, 0), (3, 0, 0), (math.log10(3), -400, -math.inf)):
            self.assertEqual((value["k"], value["residual_bits_sum"], value["completed_blocks"],
                              value["mean"]), (k, total, count, mean))
            self.assertAlmostEqual(value["log10_mean"], log_mean)
            self.assertEqual(value["raw_ber"], k / plot.N)
            self.assertAlmostEqual(value["log10_conditional_ber"], log_mean - math.log10(plot.N))
        self.assertEqual(values[1]["conditional_ber"], 0)
        self.assertTrue(math.isfinite(values[1]["log10_conditional_ber"]))
        self.assertEqual(values[2]["conditional_ber"], 0)

    def test_conditional_metric_denominators(self):
        rows = {3500: {"trials": 10, "information": 20, "full": 30}}
        for metric, total in (("information", 20), ("full", 30)):
            value, = plot.conditional_values(rows, metric)
            self.assertEqual(value["raw_ber"], 3500 / 524288)
            expected = total / (10 * plot.DENOMINATORS[metric])
            self.assertAlmostEqual(value["conditional_ber"], expected, places=18)
            self.assertAlmostEqual(value["log10_conditional_ber"], math.log10(expected))

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

    def test_recursive_discovery_sorted_and_stops_at_runs(self):
        z, _, _ = self.fixture("z", seed=2)
        (self.root / "nested").mkdir()
        a, _, _ = self.fixture("nested/a")
        self.fixture("nested/a/internal", seed=3)
        (self.root / "outputs").mkdir()
        (self.root / "outputs" / "plot.svg").write_text("not a report")
        (self.root / "outputs" / "plot.csv").write_text("not a report")
        (self.root / "unrelated.json").write_text("not JSON")
        self.assertEqual(plot.discover_reports([self.root]),
                         [a / "summary.json", z / "summary.json"])
        self.assertEqual(plot.discover_reports([a]), [a / "summary.json"])
        self.assertEqual(plot.discover_reports([z / "summary.json"]), [z / "summary.json"])

    def test_discovery_does_not_follow_subdirectory_symlinks(self):
        path, _, _ = self.fixture()
        try:
            (self.root / "loop").symlink_to(self.root, target_is_directory=True)
            (self.root / "alias").symlink_to(path, target_is_directory=True)
        except OSError as error:
            self.skipTest(f"directory symlinks unavailable: {error}")
        self.assertEqual(plot.discover_reports([self.root]), [path / "summary.json"])
        with self.assertRaisesRegex(ValueError, "duplicate report source"):
            plot.discover_reports([path, self.root / "alias"])

    def test_empty_missing_and_partial_containers(self):
        with self.assertRaisesRegex(ValueError, "no runs found"):
            plot.discover_reports([self.root])
        with self.assertRaisesRegex(ValueError, "input does not exist"):
            plot.discover_reports([self.root / "missing"])
        self.fixture()
        partial = self.root / "partial"
        partial.mkdir()
        for name in ("metadata.json", "summary.json"):
            with self.subTest(name=name):
                file = partial / name
                file.write_text("{}")
                with self.assertRaisesRegex(ValueError, "partial: incomplete run"):
                    plot.discover_reports([self.root])
                file.unlink()

    def test_discovered_malformed_and_unsupported_reports_fail(self):
        path, metadata, summary = self.fixture()
        (path / "summary.json").write_text("not JSON")
        with self.assertRaisesRegex(ValueError, "summary.json"):
            plot.pool_reports([self.root])
        for key in ("code", "random algorithm"):
            broken = dict(metadata, **{key: "unsupported"})
            self.save(path, broken, summary)
            with self.assertRaisesRegex(ValueError, "incompatible schema/code/random algorithm"):
                plot.pool_reports([self.root])

    def test_container_overlap_and_distinct_duplicate_identities(self):
        a, metadata, summary = self.fixture("a")
        for inputs in ([self.root, a], [a, self.root], [self.root, self.root],
                       [self.root, a / "summary.json"]):
            with self.subTest(inputs=inputs):
                with self.assertRaisesRegex(ValueError, "overlapping inputs.*supply each run only once"):
                    plot.pool_reports(inputs)
        b, _, _ = self.fixture("b", seed=2)
        self.save(b, metadata, summary)
        with self.assertRaisesRegex(ValueError, "duplicate run identity"):
            plot.pool_reports([self.root])

    def test_container_pooling_keeps_configuration_and_seed_checks(self):
        self.fixture("a", trials=1, residual=10)
        b, metadata, summary = self.fixture("b", seed=2, trials=9, residual=0)
        metadata["settings"].update({"threads": 2, "minimum flipped bits": 2500})
        summary["run identity"] = plot.identity(metadata)
        self.save(b, metadata, summary)
        self.fixture("cap", passes=32)
        for flag in ("anchors", "binary image"):
            path, metadata, summary = self.fixture(flag)
            metadata["settings"][flag] = False
            summary["run identity"] = plot.identity(metadata)
            self.save(path, metadata, summary)
        groups = plot.pool_reports([self.root])
        self.assertEqual(len(groups), 4)
        self.assertEqual(groups[(16, True, True)][2600],
                         {"trials": 10, "information": 10, "full": 10})
        self.fixture("duplicate-seed")
        with self.assertRaisesRegex(ValueError, "repeated seed"):
            plot.pool_reports([self.root])

    def test_discovered_input_output_protection(self):
        path, _, _ = self.fixture()
        before = {p: p.read_bytes() for p in path.iterdir()}
        for protected in before:
            with self.subTest(protected=protected):
                with self.assertRaisesRegex(ValueError, "output paths must be distinct"):
                    plot.main([str(self.root), "--output", str(self.root / "plot.svg"),
                               "--csv", str(protected)])
        self.assertEqual(before, {p: p.read_bytes() for p in before})

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "matplotlib not installed")
    def test_container_conditional_smoke_preserves_sources(self):
        a, _, _ = self.fixture("a", residual=10)
        b, _, _ = self.fixture("b", seed=2, trials=9, residual=0)
        sources = [p for run in (a, b) for p in run.iterdir()]
        before = {p: p.read_bytes() for p in sources}
        output = self.root / "merged.svg"
        for mode in ("conditional", "ber", "both"):
            with self.subTest(mode=mode), mock.patch("sys.stderr", new_callable=io.StringIO) as stderr:
                plot.main([str(self.root), "--mode", mode, "--metric", "both",
                           "--points", "2", "--output", str(output)])
                self.assertIn("Validated 2 report(s)", stderr.getvalue())
                self.assertIn("<svg", output.read_text())
                with output.with_suffix(".csv").open() as stream:
                    rows = list(csv.DictReader(stream))
                self.assertEqual({row["metric"] for row in rows}, {"information", "full"})
                if mode == "conditional":
                    self.assertEqual(len(rows), 2)
                    self.assertTrue(all(row["completed_blocks"] == "10" and
                                        float(row["mean"]) == 1 for row in rows))
        self.assertEqual(before, {p: p.read_bytes() for p in sources})

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

    def test_r4_and_default_dimensions_never_pool(self):
        default, _, _ = self.fixture("default")
        r4, metadata, summary = self.fixture("r4")
        metadata["settings"].update(n1=256, k1=224, n2=256, k2=252)
        metadata["code"] = "RS256,224 x RS256,252 Cantor systematic row major"
        summary["run identity"] = plot.identity(metadata)
        self.save(r4, metadata, summary)
        groups = plot.pool_reports([default, r4])
        self.assertEqual(len(groups), 2)
        self.assertEqual({plot.config_dimensions(config) for config in groups},
                         {(256, 224, 256, 254), (256, 224, 256, 252)})

    def test_shortened_and_default_dimensions_never_pool(self):
        default, _, _ = self.fixture("default")
        short, metadata, summary = self.fixture("short")
        metadata["settings"].update(n1=256, k1=224, n2=175, k2=173)
        metadata["settings"]["maximum flipped bits"] = 8*256*175
        metadata["code"] = "RS256,224 x RS175,173 Cantor systematic row major"
        summary["run identity"] = plot.identity(metadata)
        self.save(short, metadata, summary)
        groups = plot.pool_reports([default, short])
        self.assertEqual(len(groups), 2)
        self.assertEqual({plot.config_dimensions(config) for config in groups},
                         {(256, 224, 256, 254), (256, 224, 175, 173)})
        before = (short / "metadata.json").read_bytes()
        self.assertEqual(plot.load_report(short)[0], summary["run identity"])
        self.assertEqual(before, (short / "metadata.json").read_bytes())
        metadata["schema revision"] = summary["schema revision"] = 2
        stats = {"completed blocks": 1, "total iterations": 2,
                 "information bits": {"total bits": 8*224*173, "raw corrupted bits": 2000,
                                      "post decoding corrupted bits": 1},
                 "full-codeword bits": {"total bits": 8*256*175, "raw corrupted bits": 2600,
                                        "post decoding corrupted bits": 1}}
        summary["overall"] = {"statistics": copy.deepcopy(stats)}
        summary["by flipped bit count"] = [{"flipped bit count": 2600, "statistics": stats}]
        summary["run identity"] = plot.identity(metadata)
        self.save(short, metadata, summary)
        self.assertEqual(plot.pool_reports([default, short]), groups)
        summary["code parameters"] = {"n1": 256, "k1": 224, "n2": 175, "k2": 173}
        self.save(short, metadata, summary)
        self.assertEqual(plot.pool_reports([default, short]), groups)
        summary["code parameters"]["k2"] = 172
        self.save(short, metadata, summary)
        with self.assertRaises(ValueError):
            plot.load_report(short)
        summary["code parameters"]["k2"] = 173
        for mutation in (lambda m: m["settings"].pop("k2"),
                         lambda m: m["settings"].update(n2=174),
                         lambda m: m["settings"].update({"maximum flipped bits": 358401}),
                         lambda m: m.update(code=plot.CODE)):
            broken = copy.deepcopy(metadata)
            mutation(broken)
            self.save(short, broken, summary)
            with self.assertRaises(ValueError):
                plot.load_report(short)

    def test_explicit_default_dimensions_preserve_legacy_group_and_hash(self):
        path, metadata, summary = self.fixture()
        implicit = plot.load_report(path)
        metadata["settings"].update(zip(("n1", "k1", "n2", "k2"), plot.DEFAULT_DIMENSIONS))
        summary["run identity"] = plot.identity(metadata)
        self.save(path, metadata, summary)
        explicit = plot.load_report(path)
        self.assertEqual(implicit[1:], explicit[1:])
        self.assertEqual(explicit[0], plot.identity(metadata))

    def test_fisher_yates_metadata_preserves_summary_format(self):
        path, metadata, summary = self.fixture()
        expected = plot.load_report(path)[1:]
        metadata["random algorithm"] = plot.RANDOM_FY
        summary["run identity"] = plot.identity(metadata)
        self.save(path, metadata, summary)
        self.assertEqual(plot.load_report(path)[1:], expected)

    def test_snapshot_storage_and_checkpoint_validation(self):
        path, metadata, summary = self.fixture()
        metadata["schema revision"] = summary["schema revision"] = 2
        metadata["storage"] = plot.SNAPSHOTS
        stats = {"completed blocks": 0, "total iterations": 0,
                 "information bits": {"total bits": 0, "raw corrupted bits": 0, "post decoding corrupted bits": 0},
                 "full-codeword bits": {"total bits": 0, "raw corrupted bits": 0, "post decoding corrupted bits": 0}}
        summary["overall"] = {"statistics": stats}
        summary["by flipped bit count"] = []
        for sampler in (plot.RANDOM, plot.RANDOM_FY):
            metadata["random algorithm"] = sampler
            summary["run identity"] = plot.identity(metadata)
            if sampler == plot.RANDOM_FY:
                summary["checkpoint"] = {"flip end": 40}
            self.save(path, metadata, summary)
            self.assertEqual(plot.load_report(path)[3], {})
        for checkpoint in ({"flip end": 39}, {"flip end": True}, {"flip end": 2**64}, {}):
            summary["checkpoint"] = checkpoint
            self.save(path, metadata, summary)
            with self.assertRaises(ValueError):
                plot.load_report(path)

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

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "matplotlib not installed")
    def test_conditional_pool_union_csv_and_no_binomial(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        a, _, _ = self.fixture("a", trials=1, residual=10)
        b, _, _ = self.fixture("b", seed=2, trials=9, residual=0)
        c, metadata, summary = self.fixture("c", seed=3, trials=7, residual=0)
        summary["by flipped bit count"][0]["flipped bit count"] = 2500
        self.save(c, metadata, summary)
        d, _, _ = self.fixture("d", passes=32, residual=2)
        output = self.root / "conditional.svg"
        with mock.patch.object(plot, "evaluate", side_effect=AssertionError("BER evaluation")), \
                mock.patch.object(plot, "log_binomial", side_effect=AssertionError("binomial")), \
                mock.patch.object(plt, "close", wraps=plt.close) as close:
            plot.main([str(a), str(b / "summary.json"), str(c), str(d),
                       "--mode", "conditional", "--metric", "both", "--points", "0",
                       "--output", str(output)])
        figure = close.call_args.args[0]
        axis = figure.axes[0]
        self.assertEqual((axis.get_xscale(), axis.get_yscale()), ("linear", "log"))
        self.assertEqual(axis.get_xlim(), (0.008, 0.0045))
        self.assertEqual(axis.get_ylim(), (1e-30, 1e-1))
        self.assertEqual(len(axis.collections), 4)
        for index, mean in ((0, 1), (2, 2)):
            x, y = axis.collections[index].get_offsets().tolist()[0]
            self.assertEqual(x, 2600 / plot.N)
            self.assertAlmostEqual(y, mean / plot.DENOMINATORS["information"], places=18)
        with output.with_suffix(".csv").open() as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 6)
        self.assertEqual([int(row["k"]) for row in rows[:2]], [2500, 2600])
        self.assertEqual(rows[0]["residual_bits_sum"], "0")
        self.assertEqual(rows[0]["completed_blocks"], "7")
        self.assertEqual(float(rows[0]["mean"]), 0)
        self.assertEqual(float(rows[0]["log10_mean"]), -math.inf)
        self.assertEqual(float(rows[0]["conditional_ber"]), 0)
        self.assertEqual(float(rows[0]["log10_conditional_ber"]), -math.inf)
        self.assertEqual(rows[1]["residual_bits_sum"], "10")
        self.assertEqual(rows[1]["completed_blocks"], "10")
        self.assertEqual(float(rows[1]["mean"]), 1)
        self.assertAlmostEqual(float(rows[1]["log10_mean"]), 0)
        self.assertEqual(rows[1]["record_type"], "conditional")
        self.assertEqual(rows[1]["p"], "")
        self.assertEqual({row["metric"] for row in rows}, {"information", "full"})
        self.assertIn("1/2 zero-observed strata omitted", output.read_text())

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "matplotlib not installed")
    def test_conditional_all_zero_and_empty_no_floor(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        path, metadata, summary = self.fixture(residual=0)
        for empty in (False, True):
            with self.subTest(empty=empty):
                if empty:
                    summary["by flipped bit count"] = []
                    summary["overall"]["trial count"] = 0
                    self.save(path, metadata, summary)
                output = self.root / "zero.svg"
                with mock.patch.object(plt, "close", wraps=plt.close) as close:
                    plot.main([str(path), "--mode", "conditional", "--output", str(output)])
                axis = close.call_args.args[0].axes[0]
                self.assertEqual(len(axis.collections[0].get_offsets()), 0)
                self.assertEqual(len(axis.lines), 0)
                self.assertEqual(axis.get_xlim(), (0.008, 0.0045))
                self.assertEqual(axis.get_ylim(), (1e-30, 1e-1))
                text = output.read_text()
                self.assertIn("No artificial floor is used", text)
                self.assertIn("No completed blocks" if empty else "All observed residual sums are zero", text)
                with output.with_suffix(".csv").open() as stream:
                    rows = list(csv.DictReader(stream))
                self.assertEqual(len(rows), 0 if empty else 1)

    def test_points_validation_for_weighted_modes(self):
        for mode in ("ber", "both"):
            for points in ("0", "1"):
                with self.subTest(mode=mode, points=points):
                    with self.assertRaisesRegex(ValueError, "--points must be at least 2"):
                        plot.main([str(self.root), "--mode", mode, "--points", points])

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "matplotlib not installed")
    def test_both_overlay_pooled_weights_and_unified_csv(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import to_rgba

        a, _, _ = self.fixture("a", residual=10)
        b, _, _ = self.fixture("b", seed=2, trials=9, residual=0)
        self.fixture("cap", passes=32, residual=2)
        output, csv_path = self.root / "both.svg", self.root / "explicit.csv"
        with mock.patch.object(plt, "close", wraps=plt.close) as close:
            plot.main([str(self.root), "--mode", "both", "--metric", "both", "--points", "3",
                       "--output", str(output), "--csv", str(csv_path)])
        axis = close.call_args.args[0].axes[0]
        self.assertEqual((axis.get_xscale(), axis.get_yscale()), ("linear", "log"))
        self.assertEqual(axis.get_xlim(), (0.008, 0.0045))
        self.assertEqual(axis.get_ylim(), (1e-30, 1e-1))
        self.assertEqual((len(axis.collections), len(axis.lines)), (4, 4))
        for index, (scatter, line) in enumerate(zip(axis.collections, axis.lines)):
            self.assertEqual(tuple(scatter.get_facecolors()[0]), to_rgba(line.get_color()))
            self.assertIn("Fixed-weight conditional BER", scatter.get_label())
            self.assertIn("Sampled-stratum BSC contribution", line.get_label())
            self.assertIn("information" if index % 2 == 0 else "full", line.get_label())
        self.assertEqual(axis.lines[0].get_color(), axis.lines[1].get_color())
        self.assertNotEqual(axis.lines[0].get_color(), axis.lines[2].get_color())
        self.assertNotEqual(axis.lines[0].get_linestyle(), axis.lines[1].get_linestyle())
        self.assertNotEqual(axis.collections[0].get_paths()[0].vertices.tolist(),
                            axis.collections[1].get_paths()[0].vertices.tolist())
        pooled = plot.pool_reports([a, b])[(16, True, True)]
        with csv_path.open() as stream:
            records = list(csv.DictReader(stream))
        self.assertEqual(len(records), 16)
        conditional = [r for r in records if r["record_type"] == "conditional"]
        weighted = [r for r in records if r["record_type"] == "ber"]
        self.assertEqual((len(conditional), len(weighted)), (4, 12))
        for row in conditional[:2]:
            self.assertEqual((row["k"], row["residual_bits_sum"], row["completed_blocks"]),
                             ("2600", "10", "10"))
            self.assertEqual(float(row["mean"]), 1)
            self.assertAlmostEqual(float(row["conditional_ber"]), 1 / plot.DENOMINATORS[row["metric"]])
            self.assertEqual(row["log10_contribution"], "")
        for row in weighted[:6]:
            p, metric = float(row["p"]), row["metric"]
            expected = plot.evaluate(pooled, p)[0][metric] / math.log(10)
            self.assertEqual(float(row["log10_contribution"]), expected)
            self.assertAlmostEqual(expected, (plot.log_binomial(plot.N, 2600, p) -
                                   math.log(plot.DENOMINATORS[metric])) / math.log(10))
            for field in ("k", "residual_bits_sum", "completed_blocks", "mean", "raw_ber",
                          "conditional_ber", "log10_conditional_ber"):
                self.assertEqual(row[field], "")
        self.assertFalse(output.with_suffix(".csv").exists())

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "matplotlib not installed")
    def test_conditional_underflow_csv_not_floored(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        output = self.root / "underflow.svg"
        groups = {(16, True, True): {3500: {"trials": 10**400, "information": 1}}}
        with mock.patch.object(plt, "close", wraps=plt.close) as close:
            plot.plot_results(groups, ["information"], "conditional", [], output,
                              output.with_suffix(".csv"), plt)
        axis = close.call_args.args[0].axes[0]
        self.assertEqual(len(axis.collections[0].get_offsets()), 0)
        with output.with_suffix(".csv").open() as stream:
            row, = csv.DictReader(stream)
        self.assertEqual(row["completed_blocks"], str(10**400))
        self.assertEqual(row["residual_bits_sum"], "1")
        self.assertEqual(float(row["conditional_ber"]), 0)
        self.assertAlmostEqual(float(row["log10_conditional_ber"]), -400 - math.log10(455168))
        self.assertIn("No representable positive BERs", output.read_text())

    def test_help_and_headless_outputs(self):
        result = subprocess.run([sys.executable, str(SCRIPT), "--help"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("No extrapolation", result.stdout)
        self.assertIn("recursively", result.stdout)
        self.assertIn("merged.svg", result.stdout)
        self.assertIn("Partial or malformed runs", " ".join(result.stdout.split()))
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
