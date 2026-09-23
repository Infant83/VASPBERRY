"""Public example discovery, immutable references, and retained run outcomes."""
from __future__ import annotations

import hashlib
import importlib.util
import contextlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
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
        self.assertEqual(set(ids), {"fukui-berry-curvature", "fukui-chern", "kubo-curvature", "kubo-hall", "hall-valley",
                                    "z2", "circular-dichroism", "wavefunction"})
        for item in catalog["features"]:
            with self.subTest(feature=item["id"]):
                for name in ("readme", "runner"):
                    self.assertTrue((EXAMPLES / item[name]).is_file(), item[name])
                if item.get("input_preparation"):
                    self.assertTrue((EXAMPLES / item["input_preparation"]).is_file())
                    self.assertFalse(item["input_bundled"])
                else:
                    self.assertTrue((EXAMPLES / item["input"]).is_file(), item["input"])
                reference = EXAMPLES / item["reference"]
                for name in item["required_outputs"]:
                    self.assertTrue((reference / name).is_file(), str(reference / name))
                result = json.loads((reference / "result.json").read_text())
                self.assertEqual(result["status"], "PASS")
                self.assertEqual(result["schema_version"], 1)
                self.assertEqual(result["feature_id"], item["id"])
                self.assertEqual(result["workflow_mode"], "vasp_wavecar_calculation")
                self.assertEqual(item["default_mode"], "vasp_wavecar_calculation")
                self.assertEqual(Path(item["input"]).name, "WAVECAR")
                self.assertFalse((EXAMPLES / item["runner"]).with_name("input.json").exists())
                # Validate the stored payload itself; this is not a claim that
                # the current runner has just recomputed it.
                for name, digest in result["output_sha256"].items():
                    self.assertEqual(hashlib.sha256((reference / name).read_bytes()).hexdigest(), digest)
                if "raw_archive_sha256" in result:
                    archive = reference / "raw-output.tar.gz"
                    self.assertEqual(hashlib.sha256(archive.read_bytes()).hexdigest(),
                                     result["raw_archive_sha256"])
                    with tarfile.open(archive) as raw:
                        for name, digest in result["raw_output_sha256"].items():
                            self.assertEqual(hashlib.sha256(raw.extractfile(name).read()).hexdigest(), digest)

    def run_cli(self, *args):
        return subprocess.run([sys.executable, str(RUNNER), *map(str, args)],
                              cwd=ROOT, capture_output=True, text=True)

    @unittest.skipUnless(shutil.which("gfortran"), "gfortran required for actual VASP example")
    def test_fresh_mos2_kubo_batch_has_checked_results_and_logs(self):
        with tempfile.TemporaryDirectory() as tmp:
            binary = Path(tmp) / "vaspberry"
            built = subprocess.run([
                "gfortran", "-cpp", "-O0", "-ffixed-line-length-none",
                "-fallow-argument-mismatch", str(ROOT / "vaspberry.f"),
                "-o", str(binary), "-llapack", "-lblas",
            ], capture_output=True, text=True)
            self.assertEqual(built.returncode, 0, built.stderr)
            out = Path(tmp) / "run"
            run = self.run_cli("kubo-curvature", "--binary", binary, "--output-dir", out)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            manifest = json.loads((out / "run_manifest.json").read_text())
            self.assertEqual(manifest["status"], "PASS")
            self.assertEqual(len(manifest["features"]), 1)
            feature = manifest["features"][0]
            self.assertEqual(feature["exit_code"], 0)
            self.assertEqual(feature["missing_outputs"], [])
            self.assertTrue((out / feature["stdout"]).is_file())
            self.assertTrue((out / feature["stderr"]).is_file())
            result = json.loads((out / "kubo-curvature/result.json").read_text())
            self.assertEqual(result["status"], "PASS")
            self.assertEqual(result["workflow_mode"], "vasp_wavecar_calculation")
            checks = result["numerical_checks"]
            self.assertEqual(checks["native_rows"], 96)
            self.assertEqual(checks["valid_points_per_band"], {"17": 44, "18": 44})
            self.assertTrue(checks["reference_validity_mask_matches"])
            self.assertLess(checks["reference_comparison_max_abs_A2"], 1e-5)
            before = (out / "run_manifest.json").read_bytes()
            repeated = self.run_cli("kubo-curvature", "--binary", binary, "--output-dir", out)
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

    def test_lfs_pointer_never_falls_back_to_stored_or_model_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            pointer = Path(tmp) / "WAVECAR"
            pointer.write_text("version https://git-lfs.github.com/spec/v1\noid sha256:dummy\n")
            out = Path(tmp) / "unused"
            run = self.run_cli("hall-valley", "--bi-wavecar", pointer, "--output-dir", out)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn("LFS pointer", run.stderr)
            self.assertFalse(out.exists())

    def test_missing_full_mos2_mesh_points_to_vasp_preparation(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "unused"
            # Use an existing file for the executable preflight; it must not run.
            run = self.run_cli("fukui-berry-curvature", "--binary", RUNNER,
                               "--mos2-mesh-wavecar", Path(tmp) / "WAVECAR",
                               "--output-dir", out)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn("fukui-berry-curvature/inputs/README.md", run.stderr)
            self.assertFalse(out.exists())

    def test_zero_exit_rejected_result_is_not_reported_as_pass(self):
        spec = importlib.util.spec_from_file_location("public_batch_runner", RUNNER)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as tmp:
            examples = Path(tmp)
            catalog = {"features": [{"id": "rejected", "requires_fortran": False,
                "dataset": "mos2", "default_mode": "test_fixture", "runner": "fake.py",
                "required_outputs": ["result.json", "summary.csv", "figure.png"]}]}
            (examples / "catalog.json").write_text(json.dumps(catalog))
            wavecar = examples / "WAVECAR"
            wavecar.write_bytes(b"batch orchestration fixture")
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
                    patch.object(sys, "argv", [str(RUNNER), "rejected", "--output-dir", str(out),
                                               "--mos2-wavecar", str(wavecar)]), \
                    contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(module.main(), 1)
            manifest = json.loads((out / "run_manifest.json").read_text())
            self.assertEqual(manifest["status"], "FAIL")
            self.assertEqual(manifest["features"][0]["exit_code"], 0)
            self.assertIn("status=PASS", manifest["features"][0]["error"])


if __name__ == "__main__":
    unittest.main()
