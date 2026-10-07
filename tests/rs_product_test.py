"""Stdlib integration tests for the native fixture/search tool."""
import json
import pathlib
import subprocess
import sys
import tempfile
import unittest

EXE = sys.argv.pop(1)


def invoke(*args, text=None, status=0):
    result = subprocess.run([EXE, *args], input=text, text=True,
                            capture_output=True, timeout=30)
    assert result.returncode == status, (result.returncode, result.stderr, result.stdout[:500])
    if status == 2:
        assert not result.stdout
        return result
    return json.loads(result.stdout)


def run(fixture, status=0):
    return invoke("run", "-", text=json.dumps(fixture), status=status)


def column_fixture(n, k, column):
    return dict(strong_n=n, strong_k=k, weak_n=4, weak_k=2,
                column_masks=[dict(column=0, xor_hex=column.hex())])


def miscorrection_witness(word, t):
    positions = sorted((i for i, b in enumerate(word) if b),
                       key=lambda i: (bin(word[i]).count("1"), i))
    received = bytearray(word)
    for i in positions[-t:]:
        received[i] = 0
    return bytes(received)


class ProductTest(unittest.TestCase):
    def sample(self, scenario, *args, details=True):
        result = subprocess.run(
            [EXE, "sample", scenario, "--strong-n", "8", "--strong-k", "4",
             "--weak-n", "4", "--weak-k", "2",
             *(["--details", "1"] if details else []), *args],
            capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('\n  "code_params":', result.stdout)
        return json.loads(result.stdout)

    def test_sampling_scenarios_and_replay(self):
        for scenario in ("none", "stall", "miscorrection", "mixed",
                         "miscorrection-stall", "two-miscorrections"):
            with self.subTest(scenario=scenario):
                records = self.sample(scenario, "--random-bits", "7",
                                      "--samples", "2", "--seed", "13")
                self.assertEqual(records, self.sample(
                    scenario, "--random-bits", "7", "--samples", "2", "--seed", "13"))
                samples = records["samples"]
                decoded = self.decode_samples(records)
                self.assertEqual(len(samples), 2)
                for sample, result in zip(samples, decoded["results"]):
                    fixture = sample["fixture"]
                    bits = sum(bin(b).count("1") for mask in fixture["row_masks"]
                               for b in bytes.fromhex(mask["xor_hex"]))
                    self.assertEqual(bits, 7)
                    self.assertEqual(sample["received_weight"]["bits"],
                                     sample["planted_weight"]["bits"] + 7 -
                                     2 * sample["cancelled_planted_bits"])
                    columns = [m["column"] for m in fixture.get("column_masks", [])]
                    self.assertEqual(len(columns), len(set(columns)))
                    self.assertEqual(len(columns), dict(none=0, stall=2,
                                     miscorrection=1, mixed=2,
                                     **{"miscorrection-stall": 3,
                                        "two-miscorrections": 2})[scenario])
                    self.assertNotIn("expected", fixture)
                    fixture["options"]["use_postprocessing"] = True
                    replay = run(fixture)
                    self.assertEqual(result["outcome"], replay["outcome"])
                    self.assertEqual(set(result), {"sample_index", "outcome", "patterns"})
                expected = {outcome: sum(s["outcome"] == outcome for s in decoded["results"])
                            for outcome in ("recovered", "undetected-failure", "detected-failure")}
                self.assertEqual(decoded["outcomes"], expected)
                self.assertEqual(decoded["schema_version"], 5)

    def decode_samples(self, records):
        result = subprocess.run([EXE, "test", "-"],
                                input=json.dumps(records, indent=2),
                                text=True, capture_output=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('\n  "outcomes":', result.stdout)
        return json.loads(result.stdout)

    def test_compact_sampling_matches_detailed_results(self):
        for scenario in ("none", "stall", "miscorrection", "mixed",
                         "miscorrection-stall", "two-miscorrections"):
            args = ("--random-bits", "7", "--samples", "2", "--seed", "13")
            compact = self.sample(scenario, *args, details=False)
            detailed = self.sample(scenario, *args)
            header = compact
            self.assertEqual(set(header), {"schema_version", "code_params",
                                           "error_generation", "error_format", "sample_count", "samples"})
            self.assertEqual(header["schema_version"], 4)
            self.assertEqual(header["code_params"],
                             dict(strong=dict(n=8, k=4), weak=dict(n=4, k=2)))
            self.assertEqual(header["error_generation"]["scenario"], scenario)
            self.assertEqual(header["error_generation"]["random_bits"], 7)
            self.assertEqual(header["error_generation"]["seed"], 13)
            self.assertEqual(len(compact["samples"]), 2)
            for small, full in zip(compact["samples"], detailed["samples"]):
                self.assertEqual(small, dict(type="sample", sample_index=full["sample_index"],
                                            error_hex=full["error_hex"]))
                expected = bytearray(32)
                for mask in full["fixture"].get("column_masks", []):
                    for r, b in enumerate(bytes.fromhex(mask["xor_hex"])):
                        expected[r * 4 + mask["column"]] ^= b
                for mask in full["fixture"].get("row_masks", []):
                    for c, b in enumerate(bytes.fromhex(mask["xor_hex"])):
                        expected[mask["row"] * 4 + c] ^= b
                self.assertEqual(bytes.fromhex(small["error_hex"]), expected)
                self.assertNotIn("baseline", full)
                self.assertNotIn("postprocessed", full)
            self.assertEqual(compact["sample_count"], 2)
            self.assertEqual(self.decode_samples(compact), self.decode_samples(detailed))
            self.assertLess(len(json.dumps(compact)), len(json.dumps(detailed)))

    def test_sampling_zero_noise_and_full_noise(self):
        expected = {"none": "recovered", "stall": "recovered",
                    "miscorrection": "recovered", "mixed": "detected-failure",
                    "miscorrection-stall": "undetected-failure"}
        for scenario, outcome in expected.items():
            records = self.sample(scenario, "--random-bits", "0")
            result = self.decode_samples(records)["results"][0]
            self.assertEqual(result["outcome"], outcome)
            markers = {"none": [], "stall": ["stall-repaired"],
                       "miscorrection": ["strong-miscorrection-detected", "strong-miscorrection-repaired"],
                       "mixed": ["strong-miscorrection-detected"],
                       "miscorrection-stall": ["stall-repaired"]}
            self.assertEqual(result["patterns"], markers[scenario])
        sample = self.sample("stall", "--random-bits", "256")["samples"][0]
        self.assertEqual(sample["cancelled_planted_bits"], sample["planted_weight"]["bits"])
        self.assertEqual(sample["received_weight"]["bits"],
                         256 - sample["planted_weight"]["bits"])
        for mask in sample["fixture"]["row_masks"]:
            self.assertEqual(mask["xor_hex"], "ffffffff")

    def test_sampling_invalid_arguments(self):
        for args in (("unknown", "--random-bits", "0"),
                     ("none",), ("none", "--random-bits", "-1"),
                     ("none", "--random-bits", "524289"),
                     ("none", "--random-bits", "0", "--samples", "0"),
                     ("none", "--random-bits", "0", "--details", "2"),
                     ("none", "--random-bits", "0", "--columns", "0"),
                     ("stall", "--random-bits", "0", "--codewords", "unused")):
            invoke("sample", *args, status=2)

    def test_noise_excludes_miscorrection_columns(self):
        for scenario, protected in (("miscorrection", 1), ("mixed", 1),
                                    ("miscorrection-stall", 1),
                                    ("two-miscorrections", 2)):
            eligible = 8 * 8 * (4 - protected)
            for flips in (7, eligible):
                doc = self.sample(scenario, "--random-bits", str(flips),
                                  "--samples", "3", "--seed", "42")
                self.assertEqual(doc["error_generation"]["eligible_noise_bits"], eligible)
                self.assertEqual(doc["error_generation"]["excluded_miscorrection_columns"], protected)
                for s in doc["samples"]:
                    masks = s["fixture"]["column_masks"]
                    columns = [masks[0]["column"]]
                    if protected == 2:
                        columns.append(masks[1]["column"])
                    noise = bytearray(32)
                    for row in s["fixture"].get("row_masks", []):
                        noise[4*row["row"]:4*(row["row"]+1)] = bytes.fromhex(row["xor_hex"])
                    self.assertEqual(sum(b.bit_count() for b in noise), flips)
                    received = bytes.fromhex(s["error_hex"])
                    for c in columns:
                        self.assertEqual(noise[c::4], bytes(8))
                        planted = next(m for m in masks if m["column"] == c)
                        self.assertEqual(received[c::4], bytes.fromhex(planted["xor_hex"]))
                    if flips == eligible:
                        for c in set(range(4)) - set(columns):
                            self.assertEqual(noise[c::4], bytes([255]) * 8)
            invoke("sample", scenario, "--strong-n", "8", "--strong-k", "4",
                   "--weak-n", "4", "--weak-k", "2", "--random-bits",
                   str(eligible + 1), status=2)

    def test_two_miscorrections_construction(self):
        records = self.sample("two-miscorrections", "--random-bits", "0",
                              "--samples", "4", "--seed", "42")
        for sample in records["samples"]:
            self.assertEqual(len(sample["construction"]), 2)
            for planted in sample["construction"]:
                self.assertEqual(planted["scenario"], "miscorrection")
                checks = planted["checks"]
                self.assertTrue(checks["witness_exact_miscorrection"])
                self.assertTrue(checks["verified"])
                self.assertEqual(checks["strong_outcomes"][0]["status"], "miscorrected")
            # Two wrong columns cannot be a nonzero weak codeword (d=3).
            fixture = sample["fixture"]
            baseline = run(fixture)
            self.assertEqual(baseline["outcome"], "detected-failure")
            self.assertEqual(baseline["final"]["components"]["invalid_columns"], [])
            self.assertTrue(baseline["final"]["components"]["invalid_rows"])

    def test_sample_stream_validation(self):
        records = self.sample("none", "--random-bits", "0", details=False)
        invalid = [json.dumps(records)[:-1], json.dumps(records) + '{}',
                   json.dumps(dict(records, sample_count=2)),
                   json.dumps(dict(records, samples=[]))]
        for field, value in (("error_hex", "00"), ("error_hex", "gg" * 32),
                             ("sample_index", 1), ("type", "result")):
            changed = json.loads(json.dumps(records))
            changed["samples"][0][field] = value
            invalid.append(json.dumps(changed))
        for data in invalid:
            result = subprocess.run([EXE, "test", "-"],
                                    input=data,
                                    text=True, capture_output=True, timeout=30)
            self.assertEqual(result.returncode, 2)
            self.assertEqual(result.stdout, "")

    def check_pattern(self, report, name):
        self.assertEqual(report["scenario"], name)
        self.assertTrue(report["expected_matched"])
        self.assertTrue(report["assertion_passed"])
        construction = report["construction"]
        self.assertTrue(construction["verified"])
        self.assertLessEqual(construction["attempts_used"], construction["attempt_limit"])
        self.assertTrue(construction["baseline_matches_initial_strong_output"])
        baseline, post = report["baseline"], report["postprocessed"]
        self.assertEqual(baseline["outcome"], "detected-failure")
        self.assertEqual(baseline["correction"]["weak_changed_symbols"], 0)
        self.assertEqual(baseline["postprocessing"], "none")
        self.assertEqual(post["postprocessing"], "enabled")
        self.assertEqual(post["outcome"], "detected-failure" if name == "mixed" else "recovered")
        columns = construction["columns"]
        expected_bad = columns if name == "stall" else columns[1:]
        self.assertEqual(baseline["final"]["components"]["invalid_columns"], sorted(expected_bad))
        t = construction["strong_radius"]
        for i, outcome in enumerate(construction["strong_outcomes"]):
            self.assertTrue(outcome["verified"])
            self.assertGreater(outcome["received_weight"]["symbols"], t)
            wrong = name != "stall" and i == 0
            self.assertEqual(outcome["status"], "miscorrected" if wrong else "uncorrectable")
            self.assertEqual(outcome["corrections"], t if wrong else 0)
            self.assertEqual(outcome["unchanged"], not wrong)
        masks = baseline["fixture"]["column_masks"]
        if name == "stall":
            self.assertEqual(masks[0]["xor_hex"], masks[1]["xor_hex"])
        else:
            word = bytes.fromhex(construction["wrong_codeword_hex"])
            witness = bytes.fromhex(masks[0]["xor_hex"])
            self.assertEqual(witness, miscorrection_witness(word, t))
            self.assertEqual(bytes.fromhex(baseline["final"]["block_hex"])[columns[0]::baseline["fixture"]["weak_n"]], word)
            self.assertTrue(construction["witness_exact_miscorrection"])
            self.assertEqual(construction["witness_corrections"], t)
            if construction["source"]["kind"] == "seeded-single-information-bit":
                self.assertEqual(construction["wrong_codeword_weight"]["symbols"], 2 * t + 1)
        expected_counters = (1, 0, 0) if name == "stall" else (0, 1, int(name == "miscorrection"))
        for key, value in zip(("stall_patterns_corrected", "strong_miscorrections_detected",
                               "strong_miscorrections_corrected"), expected_counters):
            self.assertEqual(baseline["correction"][key], 0)
            self.assertEqual(post["correction"][key], value)
        if name == "mixed":
            self.assertTrue(construction["known_column_erasure_inconsistent"])
            self.assertLess(max(construction["weak_proposal_support"].values(), default=0),
                            construction["consensus_threshold"])
            self.assertEqual(post["final"], baseline["final"])
            self.assertGreater(post["final"]["full_residual"]["symbols"], 0)
            second = bytes.fromhex(masks[1]["xor_hex"])
            self.assertEqual([i for i, b in enumerate(second) if b], [i for i, b in enumerate(witness) if b])
            self.assertEqual(sum(a != b for a, b in zip(second, witness)), 1)
        else:
            self.assertEqual(post["final"]["full_residual"], dict(symbols=0, bits=0))
        for branch in (baseline, post):
            self.assertEqual(run(branch["fixture"]), branch)

    def test_patterns_small_repeat_and_parity(self):
        for name in ("stall", "miscorrection", "mixed"):
            for columns in ("0" if name == "miscorrection" else "0,1",
                            "5" if name == "miscorrection" else "5,4"):
                with self.subTest(name=name, columns=columns):
                    args = ("pattern", name, "--strong-n", "8", "--strong-k", "4",
                            "--weak-n", "6", "--weak-k", "4", "--seed", "42",
                            "--columns", columns)
                    report = invoke(*args)
                    self.assertEqual(report, invoke(*args))
                    self.check_pattern(report, name)

    def test_patterns_default_smoke(self):
        for name in ("stall", "miscorrection", "mixed"):
            with self.subTest(name=name):
                report = invoke("pattern", name)
                self.assertEqual([report["baseline"]["fixture"][key] for key in
                                  ("strong_n", "strong_k", "weak_n", "weak_k")],
                                 [256, 224, 175, 173])
                self.check_pattern(report, name)

    def test_patterns_other_dimensions(self):
        for n, k in ((16, 12), (32, 24), (64, 32), (256, 128)):
            for name in ("stall", "miscorrection", "mixed"):
                with self.subTest(n=n, k=k, name=name):
                    report = invoke("pattern", name, "--strong-n", str(n), "--strong-k", str(k),
                                    "--weak-n", "175", "--weak-k", "173",
                                    "--columns", "174" if name == "miscorrection" else "174,0")
                    self.check_pattern(report, name)

    def test_patterns_compact_generator_source(self):
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "words.json"
            for dimensions in ((), ("--strong-n", "8", "--strong-k", "4")):
                source = invoke("generate", *dimensions, "--iterations", "12", "--retain", "2")
                path.write_text(json.dumps(source))
                before = path.read_bytes()
                for name in ("miscorrection", "mixed"):
                    report = invoke("pattern", name, *dimensions, "--codewords", str(path), "--index", "1")
                    self.check_pattern(report, name)
                    self.assertEqual(report["construction"]["wrong_codeword_hex"], source["codewords"][1]["codeword"])
                    self.assertEqual(report["construction"]["source"], dict(kind="compact-codewords", index=1))
                self.assertEqual(before, path.read_bytes())

    def test_patterns_invalid_options_and_sources(self):
        base = ("pattern", "mixed", "--strong-n", "8", "--strong-k", "4",
                "--weak-n", "4", "--weak-k", "2")
        for args in (("--attempts", "0"), ("--attempts", "10001"), ("--seed", "-1"),
                     ("--seed", "18446744073709551616"), ("--columns", "0,0"),
                     ("--columns", "0,4"), ("--columns", "-1,0"), ("--columns", "0"),
                     ("--columns", "0,1,2"), ("--columns", "0,"), ("--columns", ""),
                     ("--columns", "0, 1"), ("--index", "0"), ("--unknown", "1"),
                     ("--seed", "1", "--seed", "2"), ("--seed",)):
            invoke(*base, *args, status=2)
        invoke("pattern", status=2)
        invoke("pattern", "other", status=2)
        invoke("pattern", "miscorrection", "--columns", "0,1", status=2)
        invoke("pattern", "stall", "--codewords", "unused", status=2)
        error = invoke("pattern", "stall", "--weak-n", "256", "--weak-k", "252", status=2)
        self.assertIn("unsupported for weak R=4", error.stderr)
        error = invoke("pattern", "stall", "--strong-n", "8", "--strong-k", "6", status=2)
        self.assertIn("strong R>=4", error.stderr)
        invoke("pattern", "stall", "--strong-n", "9", status=2)
        source = invoke("generate", "--strong-n", "8", "--strong-k", "4", "--iterations", "1", "--retain", "1")
        invalid = [[], {}, dict(source, unknown=0), dict(source, codewords=[]),
                   dict(source, codewords={}), dict(source, **{"code params": dict(n=16, k=12)}),
                   dict(source, codewords=[dict(input="00000000", codeword="00" * 8)]),
                   dict(source, codewords=[dict(input="ff" * 4, codeword=source["codewords"][0]["codeword"])]),
                   dict(source, codewords=[dict(input="00000000", codeword="00" * 7 + "01")]),
                   dict(source, codewords=[dict(input="zz" * 4, codeword="00" * 8)]),
                   dict(source, codewords=[dict(input="00", codeword="00" * 8)]),
                   dict(source, codewords=[dict(source["codewords"][0], extra=0)])]
        for value in invalid:
            invoke(*base, "--codewords", "-", text=json.dumps(value), status=2)
        for text in ('{"codewords":[],"codewords":[]}', json.dumps(source) + "\0", "{"):
            invoke(*base, "--codewords", "-", text=text, status=2)
        invoke(*base, "--codewords", "-", "--index", "1", text=json.dumps(source), status=2)
        invoke(*base, "--codewords", "/nonexistent/rs-product-test-codewords.json", status=2)

    def test_replayed_counter_mismatch(self):
        report = invoke("pattern", "stall", "--strong-n", "8", "--strong-k", "4",
                        "--weak-n", "4", "--weak-k", "2")
        fixture = report["postprocessed"]["fixture"]
        fixture["expected_correction"]["stall_patterns_corrected"] = 0
        self.assertFalse(run(fixture, status=1)["assertion_passed"])
        for value in (True, -1, 257, "1"):
            fixture["expected_correction"]["stall_patterns_corrected"] = value
            run(fixture, status=2)
        fixture["expected_correction"] = dict(unknown=1)
        run(fixture, status=2)

    def test_pattern_attempt_exhaustion(self):
        args = ("pattern", "mixed", "--strong-n", "8", "--strong-k", "4",
                "--weak-n", "4", "--weak-k", "2", "--seed", "3280")
        error = invoke(*args, "--attempts", "1", status=2)
        self.assertIn("construction failed: exhausted 1 attempts", error.stderr)
        self.assertIn("postprocessing was not run", error.stderr)
        report = invoke(*args, "--attempts", "2")
        self.assertEqual(report["construction"]["attempts_used"], 2)
        self.check_pattern(report, "mixed")

    def test_search_and_replay(self):
        args = ("generate", "--strong-n", "8", "--strong-k", "4",
                "--iterations", "80", "--restart-interval", "13",
                "--retain", "32", "--seed", "42")
        result = invoke(*args)
        self.assertEqual(result, invoke(*args))
        self.assertEqual(set(result), {"code params", "codewords"})
        self.assertEqual(result["code params"], dict(n=8, k=4))
        self.assertEqual(len(result["codewords"]), 32)
        self.assertTrue(any(sum(b != 0 for b in bytes.fromhex(c["input"])) > 1
                            for c in result["codewords"]))
        words = [bytes.fromhex(c["codeword"]) for c in result["codewords"]]
        self.assertEqual(len(set(words)), len(words))
        self.assertEqual(words, sorted(words, key=lambda w: (sum(bin(b).count("1") for b in w), w)))
        for candidate in result["codewords"]:
            self.assertEqual(set(candidate), {"input", "codeword"})
            self.assertRegex(candidate["input"], r"^[0-9a-f]{8}$")
            self.assertRegex(candidate["codeword"], r"^[0-9a-f]{16}$")
            word = bytes.fromhex(candidate["codeword"])
            self.assertEqual(bytes.fromhex(candidate["input"]), word[:4])
            self.assertTrue(any(word))
            witness = miscorrection_witness(word, 2)
            self.assertEqual(sum(a != b for a, b in zip(word, witness)), 2)
            self.assertGreater(sum(b != 0 for b in witness), 2)
            post = run(column_fixture(8, 4, word))
            self.assertEqual(post["initial"]["components"]["invalid_columns"], [])
            replay = run(column_fixture(8, 4, witness))
            self.assertEqual(replay["correction"]["strong_changed_symbols"], 2)
            self.assertEqual(replay["final"]["block_hex"], post["initial"]["block_hex"])
            self.assertEqual(replay["outcome"], "detected-failure")

    def test_default_exact_witness(self):
        result = invoke("generate", "--iterations", "12", "--retain", "1")
        self.assertEqual(set(result), {"code params", "codewords"})
        self.assertEqual(result["code params"], dict(n=256, k=224))
        self.assertEqual(len(result["codewords"]), 1)
        candidate = result["codewords"][0]
        self.assertEqual(set(candidate), {"input", "codeword"})
        word = bytes.fromhex(candidate["codeword"])
        received = miscorrection_witness(word, 16)
        self.assertEqual(len(word), 256)
        self.assertEqual(len(bytes.fromhex(candidate["input"])), 224)
        self.assertEqual(bytes.fromhex(candidate["input"]), word[:224])
        self.assertEqual(sum(a != b for a, b in zip(word, received)), 16)
        positions = sorted((i for i, b in enumerate(word) if b),
                           key=lambda i: (bin(word[i]).count("1"), i))
        self.assertEqual([i for i, b in enumerate(received) if b], sorted(positions[:-16]))
        self.assertGreater(sum(b != 0 for b in received), 16)
        post = run(column_fixture(256, 224, word))
        self.assertEqual(post["initial"]["components"]["invalid_columns"], [])
        replay = run(column_fixture(256, 224, received))
        self.assertEqual(replay["correction"]["strong_changed_symbols"], 16)
        self.assertEqual(bytes.fromhex(replay["final"]["block_hex"])[::4], word)

    def test_recovered_parity_and_explicit_transmission(self):
        base = dict(strong_n=8, strong_k=4, weak_n=4, weak_k=2)
        fixture = dict(base, errors=[dict(row=7, column=3, xor_hex="ff")], expected="recovered")
        report = run(fixture)
        self.assertEqual(report["initial"]["full_residual"], dict(symbols=1, bits=8))
        self.assertEqual(report["initial"]["information_residual"]["bits"], 0)
        self.assertEqual(report["initial"]["residual_errors"][0]["column"], 3)
        self.assertTrue(report["final"]["components"]["valid"])
        # Constant nonzero data evaluates to a constant full component word.
        explicit = dict(base, transmitted_hex="01" * 32, expected="recovered")
        self.assertEqual(run(explicit)["outcome"], "recovered")
        wrong = dict(base, row_masks=[dict(row=r, xor_hex="01010101") for r in range(8)],
                     expected="undetected-failure")
        self.assertEqual(run(wrong)["outcome"], "undetected-failure")
        wrong["expected"] = "recovered"
        self.assertFalse(run(wrong, status=1)["assertion_passed"])
        cancel = dict(base, errors=[dict(row=0, column=0, xor_hex="03")] * 2)
        self.assertEqual(run(cancel)["initial"]["full_residual"]["bits"], 0)

    def test_invalid_input_and_no_overwrite(self):
        base = dict(strong_n=8, strong_k=4, weak_n=4, weak_k=2)
        invalid = [dict(strong_n=257), dict(strong_k=8), dict(strong_n=True),
                   dict(options={"use_anchors": 1}), dict(options={"max_directional_passes": 1}),
                   dict(expected="future-recovery"), dict(transmitted_hex="00"),
                   dict(transmitted_hex="01" + "00" * 31), dict(unknown=0)]
        for update in invalid:
            run(dict(base, **update), status=2)
        for row, col, mask in [(8, 0, "01"), (0, 4, "01"), (-1, 0, "01"),
                                (0, 0, "zz"), (0, 0, "1"), (0, 0, " 1")]:
            run(dict(base, errors=[dict(row=row, column=col, xor_hex=mask)]), status=2)
        run(dict(base, column_masks=[dict(column=4, xor_hex="00" * 8)]), status=2)
        run(dict(base, row_masks=[dict(row=0, xor_hex="00")]), status=2)
        for text in ["{", '{"strong_n":8,"strong_n":8}', "[]", '{} garbage']:
            invoke("run", "-", text=text, status=2)
        for text in [json.dumps(base) + "\0 garbage", json.dumps(base) + "\0",
                     "\0" + json.dumps(base)]:
            result = invoke("run", "-", text=text, status=2)
            self.assertIn("raw NUL byte", result.stderr)
        for args in [("--seed", "-1"), ("--iterations", "0"), ("--retain", "33"),
                     ("--max-perturb-bits", "1"), ("--unknown", "2")]:
            invoke("generate", *args, status=2)
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "bad.json"
            path.write_text('{"transmitted_hex":"ff"}')
            before = path.read_bytes()
            invoke("run", str(path), status=2)
            self.assertEqual(before, path.read_bytes())

    def test_help(self):
        result = subprocess.run([EXE, "--help"], text=True, capture_output=True, check=True)
        self.assertIn("NOT a global minimum", result.stdout)
        self.assertIn('"code params":{"n":N,"k":K}', result.stdout)
        self.assertIn('"codewords":[{"input":"<K-byte hex>","codeword":"<N-byte hex>"}', result.stdout)
        self.assertIn("No metadata, weights, witnesses, or fixtures are embedded", result.stdout)
        self.assertNotIn("post_miscorrection_fixture", result.stdout)

    def test_examples_and_options(self):
        for name in ["product_genuine_witness.json", "product_post_miscorrection.json"]:
            path = pathlib.Path(__file__).parent / "fixtures" / name
            report = invoke("run", str(path))
            self.assertEqual(report["outcome"], "detected-failure")
            fixture = json.loads(path.read_text())
            fixture["options"] = dict(use_anchors=False, use_binary_image=False,
                                      max_directional_passes=4)
            fixture["expected"] = "recovered"
            self.assertEqual(run(fixture)["outcome"], "recovered")

    def test_postprocessing_fixture_option_and_stats(self):
        fixture = column_fixture(256, 224, bytes([255]) * 256)
        baseline = run(fixture)
        self.assertEqual(baseline["outcome"], "detected-failure")
        fixture["options"] = dict(use_postprocessing=True, max_directional_passes=2)
        report = run(fixture)
        self.assertEqual(report["outcome"], "recovered")
        self.assertEqual(report["correction"]["strong_miscorrections_detected"], 1)
        self.assertEqual(report["correction"]["strong_miscorrections_corrected"], 1)
        self.assertEqual(report["correction"]["stall_patterns_corrected"], 0)
        self.assertEqual(report["correction"]["weak_changed_symbols"], 256)
        fixture["options"]["use_postprocessing"] = 1
        run(fixture, status=2)


if __name__ == "__main__":
    unittest.main()
