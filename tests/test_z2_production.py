"""Public production-path tests for the guarded Wilson-loop implementation.

The fixtures below only construct occupied eigenspaces of a time-reversal
doubled Qi-Wu-Zhang Hamiltonian.  Wilson links, spectra, largest-gap parity,
and every validity guard are evaluated by :mod:`wavecar_z2` itself.
"""

from __future__ import annotations

import copy
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import wavecar_z2 as z2  # noqa: E402


SIGMA_X = np.array(((0.0, 1.0), (1.0, 0.0)), dtype=complex)
SIGMA_Y = np.array(((0.0, -1.0j), (1.0j, 0.0)), dtype=complex)
SIGMA_Z = np.array(((1.0, 0.0), (0.0, -1.0)), dtype=complex)
ZERO_2 = np.zeros((2, 2), dtype=complex)
IDENTITY_2 = np.eye(2, dtype=complex)


def doubled_qwz_occupied_frames(
    mass: float, mesh: int = 32
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return frames, energies, and the unitary part of spinful Theta.

    In the basis ``(up orbital 1/2, down orbital 1/2)``, the Hamiltonian is
    ``diag(h(k), h*(-k))``.  Thus ``Theta = (i s_y) K`` squares to -1 and the
    two occupied bands form one Kramers pair.  ``mass=-3`` is trivial, while
    ``mass=-1`` contains two opposite-Chern QWZ blocks and has Z2 parity one.
    """

    frames = np.empty((mesh, mesh, 4, 2), dtype=complex)
    energies = np.empty((mesh, mesh, 4), dtype=float)

    def qwz(kx: float, ky: float) -> np.ndarray:
        return (
            np.sin(kx) * SIGMA_X
            + np.sin(ky) * SIGMA_Y
            + (mass + np.cos(kx) + np.cos(ky)) * SIGMA_Z
        )

    for ix in range(mesh):
        kx = 2.0 * np.pi * ix / mesh
        for iy in range(mesh):
            ky = 2.0 * np.pi * iy / mesh
            hamiltonian = np.block(
                [[qwz(kx, ky), ZERO_2], [ZERO_2, qwz(-kx, -ky).conj()]]
            )
            eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)
            energies[ix, iy] = eigenvalues
            frames[ix, iy] = eigenvectors[:, :2]

    theta = np.block([[ZERO_2, -IDENTITY_2], [IDENTITY_2, ZERO_2]])
    return frames, energies, theta


class Z2ProductionPathTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = {}
        for mass in (-3.0, -1.0):
            frames, energies, theta = doubled_qwz_occupied_frames(mass)
            cls.results[mass] = z2.analyze_analytic_frames(
                frames, energies, theta, axes=(0, 1)
            )

    def test_trivial_and_topological_doubled_qwz_on_both_loop_axes(self):
        for mass, expected in ((-3.0, 0), (-1.0, 1)):
            with self.subTest(mass=mass):
                result = self.results[mass]
                self.assertTrue(result["valid"], result["failed_guards"])
                self.assertEqual(result["z2"], expected)
                self.assertEqual(result["candidate_z2"], expected)
                self.assertEqual(result["failed_guards"], [])
                self.assertEqual(
                    {report["candidate_z2"] for report in result["axes"].values()},
                    {expected},
                )

    def test_production_time_reversal_and_kramers_guards_pass(self):
        for mass, result in self.results.items():
            with self.subTest(mass=mass):
                checks = {item["name"]: item for item in result["checks"]}
                self.assertTrue(checks["band_energy_time_reversal"]["pass"])
                self.assertTrue(checks["occupied_projector_time_reversal"]["pass"])
                for axis in ("x", "y"):
                    self.assertTrue(checks[f"{axis}_loop_wcc_time_reversal_set"]["pass"])
                    self.assertTrue(checks[f"{axis}_loop_endpoint_kramers_pairing"]["pass"])

                self.assertLess(result["metrics"]["max_projector_tr_residual"], 1.0e-12)
                self.assertLess(result["metrics"]["max_flux_odd_residual_rad"], 1.0e-12)

    def test_broken_projector_time_reversal_is_not_reportable(self):
        frames, energies, theta = doubled_qwz_occupied_frames(-3.0)
        # Mix one occupied vector with a genuinely unoccupied basis direction
        # at only one non-TRIM k point.  This changes its projector, rather than
        # merely choosing a different gauge inside the occupied subspace.
        altered = frames.copy()
        local_frame = altered[1, 2]
        complement = np.eye(4, dtype=complex) - local_frame @ local_frame.conj().T
        column = int(np.argmax(np.linalg.norm(complement, axis=0)))
        unoccupied_direction = complement[:, column]
        unoccupied_direction /= np.linalg.norm(unoccupied_direction)
        vector = altered[1, 2, :, 0] + 0.2 * unoccupied_direction
        altered[1, 2, :, 0] = vector / np.linalg.norm(vector)

        result = z2.analyze_analytic_frames(
            altered, energies, theta, axes=(0, 1)
        )

        self.assertFalse(result["valid"])
        self.assertIsNone(result["z2"])
        self.assertIn("occupied_projector_time_reversal", result["failed_guards"])

    def _run_cli_with_result(
        self, input_path: Path, output_dir: Path, result: dict[str, object]
    ) -> int:
        argv = [
            str(input_path),
            "--nx",
            "32",
            "--ny",
            "32",
            "--occupied-bands",
            "2",
            "--axis",
            "both",
            "--output-dir",
            str(output_dir),
        ]
        with (
            mock.patch.object(z2, "Wavecar", return_value=object()),
            mock.patch.object(
                z2, "analyze_wavecar", side_effect=lambda *args: copy.deepcopy(result)
            ),
            redirect_stdout(io.StringIO()),
            redirect_stderr(io.StringIO()),
        ):
            return z2.main(argv)

    def test_valid_cli_result_writes_csv_and_json(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            wavecar = root / "WAVECAR"
            wavecar.write_bytes(b"synthetic test sentinel")
            output = root / "result"

            return_code = self._run_cli_with_result(
                wavecar, output, self.results[-1.0]
            )

            self.assertEqual(return_code, 0)
            self.assertEqual(wavecar.read_bytes(), b"synthetic test sentinel")
            csv_path = output / z2.DEFAULT_CSV
            json_path = output / z2.DEFAULT_JSON
            self.assertTrue(csv_path.is_file())
            self.assertTrue(json_path.is_file())
            self.assertFalse((output / z2.DEFAULT_PLOT).exists())
            self.assertEqual(
                csv_path.read_text(encoding="utf-8").splitlines()[0],
                "loop_axis,pump_axis,pump_fractional,pump_signed,wcc_index,wcc_cycles",
            )
            payload = json.loads(json_path.read_text(encoding="utf-8"))
            self.assertTrue(payload["valid"])
            self.assertEqual(payload["z2"], 1)
            self.assertEqual(payload["candidate_z2"], 1)
            self.assertEqual(payload["input"], "WAVECAR")

    def test_invalid_cli_result_keeps_diagnostics_but_not_a_reportable_z2(self):
        invalid = copy.deepcopy(self.results[-1.0])
        invalid["valid"] = False
        invalid["z2"] = None
        invalid["failed_guards"] = ["x_loop_wcc_position_converged"]
        for check in invalid["checks"]:
            if check["name"] == "x_loop_wcc_position_converged":
                check["pass"] = False

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            wavecar = root / "WAVECAR"
            wavecar.write_bytes(b"synthetic test sentinel")
            output = root / "result"

            return_code = self._run_cli_with_result(wavecar, output, invalid)

            self.assertEqual(return_code, 2)
            self.assertTrue((output / z2.DEFAULT_CSV).is_file())
            json_path = output / z2.DEFAULT_JSON
            self.assertTrue(json_path.is_file())
            payload = json.loads(json_path.read_text(encoding="utf-8"))
            self.assertFalse(payload["valid"])
            self.assertIsNone(payload["z2"])
            self.assertEqual(payload["candidate_z2"], 1)
            self.assertEqual(
                payload["failed_guards"], ["x_loop_wcc_position_converged"]
            )

    def test_cli_output_collision_never_overwrites_the_input(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            wavecar = root / "WAVECAR"
            original = b"input bytes must survive an output collision"
            wavecar.write_bytes(original)
            argv = [
                str(wavecar),
                "--nx",
                "32",
                "--ny",
                "32",
                "--occupied-bands",
                "2",
                "--output-dir",
                str(root),
                "--csv-name",
                wavecar.name,
            ]
            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                return_code = z2.main(argv)

            self.assertEqual(return_code, 1)
            self.assertEqual(wavecar.read_bytes(), original)
            self.assertFalse((root / z2.DEFAULT_JSON).exists())


if __name__ == "__main__":
    unittest.main()
