"""Public example discovery, immutable references, and retained run outcomes."""
from __future__ import annotations

import hashlib
import importlib.util
import contextlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
RUNNER = EXAMPLES / "run_examples.py"


class FeatureExamplesTests(unittest.TestCase):
    def test_catalog_links_and_reference_payloads(self):
        catalog = json.loads((EXAMPLES / "catalog.json").read_text())
        self.assertEqual(catalog["schema_version"], 1)
        ids = [item["id"] for item in catalog["features"]]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertEqual(set(ids), {"fukui-chern", "kubo-curvature", "hall-valley",
                                    "z2", "circular-dichroism", "wavefunction"})
        for item in catalog["features"]:
            with self.subTest(feature=item["id"]):
                for name in ("readme", "runner", "input"):
                    self.assertTrue((EXAMPLES / item[name]).is_file(), item[name])
                reference = EXAMPLES / item["reference"]
                for name in item["required_outputs"]:
                    self.assertTrue((reference / name).is_file(), str(reference / name))
                result = json.loads((reference / "result.json").read_text())
                self.assertEqual(result["status"], "PASS")
                self.assertEqual(result["schema_version"], 1)
                self.assertEqual(result["feature_id"], item["id"])
                self.assertTrue(result["workflow_mode"])
                # Validate the stored payload itself; this is not a claim that
                # the current runner has just recomputed it.
                for name, digest in result["output_sha256"].items():
                    self.assertEqual(hashlib.sha256((reference / name).read_bytes()).hexdigest(), digest)
        modes = {item["id"]: item["default_mode"] for item in catalog["features"]}
        self.assertEqual(modes["z2"], "stored_reference_validation")

    def run_cli(self, *args):
        return subprocess.run([sys.executable, str(RUNNER), *map(str, args)],
                              cwd=ROOT, capture_output=True, text=True)

    def test_fresh_fukui_batch_has_checked_results_and_logs(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "run"
            run = self.run_cli("fukui-chern", "--output-dir", out)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            manifest = json.loads((out / "run_manifest.json").read_text())
            self.assertEqual(manifest["status"], "PASS")
            self.assertEqual(len(manifest["features"]), 1)
            feature = manifest["features"][0]
            self.assertEqual(feature["exit_code"], 0)
            self.assertEqual(feature["missing_outputs"], [])
            self.assertTrue((out / feature["stdout"]).is_file())
            self.assertTrue((out / feature["stderr"]).is_file())
            result = json.loads((out / "fukui-chern/result.json").read_text())
            self.assertEqual(result["status"], "PASS")
            self.assertEqual([case["expected_chern"] for case in result["cases"]], [1, -1, 0])
            for case in result["cases"]:
                self.assertLess(case["absolute_error"], 1e-10)
                self.assertAlmostEqual(case["lower_band_sigma_mu0_e2h"], -case["chern"], places=10)
            before = (out / "run_manifest.json").read_bytes()
            repeated = self.run_cli("fukui-chern", "--output-dir", out)
            self.assertNotEqual(repeated.returncode, 0)
            self.assertEqual((out / "run_manifest.json").read_bytes(), before)

    def test_failed_feature_retains_manifest_and_stderr(self):
        with tempfile.TemporaryDirectory() as tmp:
            # A regular non-executable file passes the existence preflight and
            # forces an actual subprocess failure inside the selected workflow.
            binary = Path(tmp) / "not-executable"
            binary.write_text("not a program\n")
            out = Path(tmp) / "failed"
            run = self.run_cli("wavefunction", "--binary", binary, "--output-dir", out)
            self.assertNotEqual(run.returncode, 0)
            manifest = json.loads((out / "run_manifest.json").read_text())
            self.assertEqual(manifest["status"], "FAIL")
            feature = manifest["features"][0]
            self.assertEqual(feature["status"], "FAIL")
            self.assertNotEqual(feature["exit_code"], 0)
            self.assertTrue((out / feature["stderr"]).read_text().strip())

    def test_unknown_feature_is_rejected_before_creating_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "unused"
            run = self.run_cli("typo-feature", "--output-dir", out)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn("unknown feature", run.stderr)
            self.assertFalse(out.exists())

    def test_zero_exit_rejected_result_is_not_reported_as_pass(self):
        spec = importlib.util.spec_from_file_location("public_batch_runner", RUNNER)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as tmp:
            examples = Path(tmp)
            catalog = {"features": [{"id": "rejected", "requires_fortran": False,
                "default_mode": "test_fixture", "runner": "fake.py",
                "required_outputs": ["result.json", "summary.csv", "figure.png"]}]}
            (examples / "catalog.json").write_text(json.dumps(catalog))
            # A workflow may exit normally while explicitly rejecting its
            # scientific result. Its declaration must override process exit 0.
            (examples / "fake.py").write_text(
                "import json,sys\nfrom pathlib import Path\n"
                "out=Path(sys.argv[2]);out.mkdir()\n"
                "(out/'summary.csv').write_text('value\\n')\n"
                "(out/'figure.png').write_bytes(b'placeholder')\n"
                "(out/'result.json').write_text(json.dumps({'schema_version':1,"
                "'feature_id':'rejected','status':'FAIL'}))\n")
            out = examples / "results"
            with patch.object(module, "EXAMPLES", examples), \
                    patch.object(sys, "argv", [str(RUNNER), "rejected", "--output-dir", str(out)]), \
                    contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(module.main(), 1)
            manifest = json.loads((out / "run_manifest.json").read_text())
            self.assertEqual(manifest["status"], "FAIL")
            self.assertEqual(manifest["features"][0]["exit_code"], 0)
            self.assertIn("status=PASS", manifest["features"][0]["error"])


if __name__ == "__main__":
    unittest.main()
