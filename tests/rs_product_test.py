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
                     expected="valid-wrong")
        self.assertEqual(run(wrong)["outcome"], "valid-wrong")
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


if __name__ == "__main__":
    unittest.main()
