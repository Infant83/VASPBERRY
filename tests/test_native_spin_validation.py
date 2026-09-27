"""Focused checks for the prebuilt-binary spin CI runner and its loop oracle."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

from test_fortran_spin_kubo_math import model

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "native_spin_validation", ROOT / ".github/scripts/validate_native_spin.py")
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)


class NativeSpinValidationTests(unittest.TestCase):
    def test_prebuilt_output_must_match_source_release_version(self):
        version = (ROOT / "VERSION").read_text().strip()
        runner.require_current_producer({"producer": f"VASPBERRY {version}"})
        for metadata in ({}, {"producer": "VASPBERRY 0.0.0"}):
            with self.assertRaisesRegex(AssertionError, "does not match VERSION"):
                runner.require_current_producer(metadata)

    def test_overlap_oracle_matches_known_conserved_spin_model(self):
        k = np.array([.4, .7, .2])
        c, dc = model(k, mixing=0.)
        x, y, _ = k
        d = np.array([np.sin(x), np.sin(y), -1+np.cos(x)+np.cos(y)])
        dx = np.array([np.cos(x), 0, -np.sin(x)])
        dy = np.array([0, np.cos(y), -np.sin(y)])
        expected = np.zeros((3, 3))
        expected[:2, 2] = np.dot(d, np.cross(dx, dy))/(2*np.linalg.norm(d)**3)*np.array([1, -1])
        np.testing.assert_allclose(runner.loop_curvature(c, dc, 2e-4), expected, atol=2e-8)

    def test_loop_oracle_keeps_spin_mixing_and_nonorthogonal_frame(self):
        c, dc = model((.4, .7, .2))
        normal = runner.loop_curvature(c, dc, 2e-4)
        matrix = np.array([[1.4, .2+.3j], [-.1j, .65]])
        transformed = runner.loop_curvature(c @ matrix, dc @ matrix, 2e-4)
        np.testing.assert_allclose(transformed, normal, atol=3e-8, rtol=1e-6)
        self.assertGreater(np.max(abs(normal[:, :2])), 1e-3)
        np.testing.assert_allclose(normal[:2].sum(axis=0), normal[2], atol=3e-8)

    def test_launch_failure_retains_receipt_and_streams(self):
        with tempfile.TemporaryDirectory() as directory:
            case = Path(directory)
            command = [sys.executable, "-c", "import sys; print('started'); sys.exit(7)"]
            with self.assertRaisesRegex(AssertionError, "exit 7"):
                runner.invoke(command, case)
            receipt = json.loads((case / "command.json").read_text())
            self.assertEqual(receipt["command"], command)
            self.assertEqual(receipt["exit_code"], 7)
            self.assertEqual(receipt["status"], "FAIL")
            self.assertIn("started", (case / "stdout.log").read_text())
            self.assertTrue((case / "stderr.log").is_file())

    def test_timeout_is_failure_with_saved_receipt(self):
        with tempfile.TemporaryDirectory() as directory:
            case = Path(directory)
            with self.assertRaises(runner.subprocess.TimeoutExpired):
                runner.invoke([sys.executable, "-c", "import time; time.sleep(5)"], case, timeout=.02)
            receipt = json.loads((case / "command.json").read_text())
            self.assertEqual(receipt["status"], "FAIL")
            self.assertIn("TimeoutExpired", receipt["error"])

    def test_expected_rejection_requires_nonzero_exit_message_and_no_output(self):
        with tempfile.TemporaryDirectory() as directory:
            case = Path(directory)
            command = [sys.executable, "-c", "import sys; sys.stderr.write('bad gap'); sys.exit(1)"]
            runner.invoke(command, case, expected_error="bad gap")
            self.assertEqual(json.loads((case / "command.json").read_text())["status"], "EXPECTED_REJECTION")
            (case / "SPIN_KUBO.csv").write_text("stale result")
            with self.assertRaisesRegex(AssertionError, "produced result files"):
                runner.invoke(command, case, expected_error="bad gap")

    def test_existing_output_directory_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            sentinel = Path(directory) / "summary.json"
            sentinel.write_text("keep")
            with self.assertRaises(FileExistsError):
                runner.main(["--serial", sys.executable, "--mpi", sys.executable,
                             "--output-dir", directory])
            self.assertEqual(sentinel.read_text(), "keep")

    def test_validation_failure_still_records_partial_outputs_and_hashes(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "new"

            def fail(args, folder):
                (folder / "partial.txt").write_text("diagnostic")
                raise AssertionError("reference mismatch")

            with mock.patch.object(runner, "validate", side_effect=fail):
                with self.assertRaisesRegex(AssertionError, "reference mismatch"):
                    runner.main(["--serial", sys.executable, "--mpi", sys.executable,
                                 "--output-dir", str(output)])
            self.assertEqual(json.loads((output / "summary.json").read_text())["status"], "FAIL")
            hashes = json.loads((output / "checksums.json").read_text())
            self.assertEqual(hashes["partial.txt"]["sha256"], runner.sha256(output / "partial.txt"))
            self.assertIn("summary.json", hashes)

    def test_uses_only_supplied_binaries_and_preserves_launcher_arguments(self):
        with tempfile.TemporaryDirectory() as directory:
            args = argparse.Namespace(serial=sys.executable, mpi=sys.executable,
                                      mpiexec=sys.executable, mpi_arg=["--oversubscribe"], timeout=120)
            calls = []
            with mock.patch.object(runner, "invoke", side_effect=lambda command, *a: calls.append(command)):
                # No outputs are synthesized: validation must fail after launching
                # the two supplied executables, rather than compile a substitute.
                with self.assertRaises(FileNotFoundError):
                    runner.validate(args, Path(directory))
            binary = str(Path(sys.executable).resolve())
            self.assertEqual(calls[0][0], binary)
            self.assertEqual(calls[1][:5], [binary, "--oversubscribe", "-n", "2", binary])
            self.assertEqual(len(calls), 2)


if __name__ == "__main__":
    unittest.main()
