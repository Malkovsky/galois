"""Portable CLI smoke checks, also usable with the tests-disabled preset."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

CLI = Path(sys.argv.pop(1)).resolve()
ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("plot_product", ROOT / "scripts/plot_product_monte_carlo.py")
plot = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(plot)


class MonteCarloSmokeTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="product-mc-smoke-")
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def invoke(self, *args):
        result = subprocess.run([str(CLI), *map(str, args)], capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result

    def run_case(self, name, k, *options, default_dimensions=False):
        path = self.root / name
        dims = () if default_dimensions else ("--n1", 64, "--k1", 56, "--n2", 175, "--k2", 173)
        self.invoke("--output", path, *dims, "--minimum-flipped-bits", k,
                    "--maximum-flipped-bits", k, "--batch-size", 24, "--batches", 1,
                    "--threads", 2, "--seed", 42, *options)
        metadata = json.loads((path / "metadata.json").read_text())
        summary = json.loads((path / "summary.json").read_text())
        _, _, config, rows, _ = plot.load_report(path)
        expected_dims = plot.DEFAULT_DIMENSIONS if default_dimensions else (64, 56, 175, 173)
        self.assertEqual(plot.config_dimensions(config), expected_dims)
        self.assertEqual(rows[k]["trials"], 24)
        stats = summary["overall"]["statistics"]
        self.assertEqual(stats["completed blocks"], 24)
        self.assertEqual(stats["full-codeword bits"]["raw corrupted bits"], 24*k)
        self.assertEqual(stats["full-codeword bits"]["total bits"], 24*8*expected_dims[0]*expected_dims[2])
        self.assertEqual(stats["information bits"]["total bits"], 24*8*expected_dims[1]*expected_dims[3])
        return path, metadata, summary

    def test_zero_and_single_bit_errors_recover(self):
        for k in (0, 1):
            with self.subTest(k=k):
                _, _, summary = self.run_case(str(k), k)
                stats = summary["overall"]["statistics"]
                for metric in ("information bits", "full-codeword bits"):
                    self.assertEqual(stats[metric]["post decoding corrupted bits"], 0)

    def test_default_dimensions_build_and_run(self):
        self.run_case("default", 0, default_dimensions=True)

    def test_seeded_results_do_not_depend_on_worker_count(self):
        _, _, first = self.run_case("serial", 320, "--threads", 1)
        _, _, second = self.run_case("parallel", 320, "--threads", 3)
        self.assertEqual(first["overall"], second["overall"])
        self.assertEqual(first["by flipped bit count"], second["by flipped bit count"])

    def test_saved_flips_postprocessing_report_and_replay(self):
        path, metadata, _ = self.run_case("saved", 448, "--postprocessing", "--sampler", "fisher-yates")
        self.assertTrue(metadata["settings"]["postprocessing"])
        self.assertEqual((path / "flips.bin").read_bytes()[:8], b"RSFLIP02")
        before = (path / "summary.json").read_bytes()
        self.invoke("--report", path)
        result = self.invoke("--replay", path)
        self.assertIn("all 3 postprocessing counters match", result.stderr)
        self.assertEqual((path / "summary.json").read_bytes(), before)

    def test_ideal_decoder_report_and_replay_with_both_repair_modes(self):
        for postprocessing in (False, True):
            with self.subTest(postprocessing=postprocessing):
                options = ("--postprocessing",) if postprocessing else ()
                path, metadata, summary = self.run_case(
                    f"ideal-{postprocessing}", 448, "--ideal-decoder",
                    "--sampler", "fisher-yates", *options)
                self.assertTrue(metadata["settings"]["ideal decoder"])
                _, _, config, _, _ = plot.load_report(path)
                self.assertIn("decoder=ideal", plot.config_label(config))
                self.assertEqual(plot.config_dimensions(config), (64, 56, 175, 173))
                before = (path / "summary.json").read_bytes()
                self.invoke("--report", path)
                self.invoke("--replay", path)
                self.assertEqual((path / "summary.json").read_bytes(), before)
                stats = summary["overall"]["statistics"]
                self.assertLessEqual(stats["full-codeword bits"]["post decoding corrupted bits"], 24*448)
                for name in ("strong miscorrections detected", "strong miscorrections corrected"):
                    self.assertEqual(stats[name], 0)
                _, _, serial = self.run_case(f"ideal-serial-{postprocessing}", 448,
                                             "--ideal-decoder", "--threads", 1, *options)
                _, _, parallel = self.run_case(f"ideal-parallel-{postprocessing}", 448,
                                               "--ideal-decoder", "--threads", 3, *options)
                self.assertEqual(serial["overall"], parallel["overall"])

    def test_ideal_decoder_zero_errors_and_full_error_complement(self):
        for k in (0, 8*64*175):
            with self.subTest(k=k):
                path, _, summary = self.run_case(f"ideal-extreme-{k}", k,
                                                "--ideal-decoder", "--sampler", "fisher-yates")
                stats = summary["overall"]["statistics"]
                self.assertEqual(stats["full-codeword bits"]["post decoding corrupted bits"], 24*k)
                self.invoke("--replay", path)


if __name__ == "__main__":
    unittest.main()
