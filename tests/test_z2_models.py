import importlib.util
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "wavecar_z2", ROOT / "tools" / "wavecar_z2.py"
)
Z2 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = Z2
SPEC.loader.exec_module(Z2)


def doubled_qwz(mass: float, size: int = 32):
    """Time-reversed QWZ blocks with one occupied state per block."""

    sigma_x = np.asarray(((0.0, 1.0), (1.0, 0.0)), dtype=complex)
    sigma_y = np.asarray(((0.0, -1.0j), (1.0j, 0.0)), dtype=complex)
    sigma_z = np.asarray(((1.0, 0.0), (0.0, -1.0)), dtype=complex)
    frames = np.empty((size, size, 4, 2), dtype=complex)
    energies = np.empty((size, size, 4), dtype=float)

    def block(kx, ky):
        return (
            math.sin(kx) * sigma_x
            + math.sin(ky) * sigma_y
            + (mass + math.cos(kx) + math.cos(ky)) * sigma_z
        )

    for ix in range(size):
        kx = 2.0 * math.pi * ix / size
        for iy in range(size):
            ky = 2.0 * math.pi * iy / size
            energy_up, vector_up = np.linalg.eigh(block(kx, ky))
            energy_down, vector_down = np.linalg.eigh(block(-kx, -ky).conj())
            frames[ix, iy, :, 0] = np.r_[vector_up[:, 0], np.zeros(2)]
            frames[ix, iy, :, 1] = np.r_[np.zeros(2), vector_down[:, 0]]
            energies[ix, iy] = np.sort(
                np.r_[energy_up[0], energy_down[0], energy_up[1], energy_down[1]]
            )
    theta = np.block(
        [[np.zeros((2, 2)), -np.eye(2)], [np.eye(2), np.zeros((2, 2))]]
    ).astype(complex)
    return frames, energies, theta


def random_unitary(rank: int, rng: np.random.Generator) -> np.ndarray:
    trial = rng.normal(size=(rank, rank)) + 1.0j * rng.normal(size=(rank, rank))
    unitary, diagonal = np.linalg.qr(trial)
    phases = np.diag(diagonal).copy()
    phases /= np.where(np.abs(phases) > 0.0, np.abs(phases), 1.0)
    return unitary @ np.diag(phases.conj())


class Z2ModelTests(unittest.TestCase):
    def test_version_is_release_version(self):
        self.assertEqual(Z2.__version__, "1.1.1")

    def test_polar_unitary_is_unitary_and_preserves_singular_values(self):
        matrix = np.asarray(((1.0 + 0.2j, 0.3), (-0.4j, 0.8 - 0.1j)))
        polar, singular = Z2.polar_unitary(matrix)
        np.testing.assert_allclose(polar.conj().T @ polar, np.eye(2), atol=1e-14)
        np.testing.assert_allclose(singular, np.linalg.svd(matrix, compute_uv=False))

    def test_skinny_projector_distance_matches_dense_reference(self):
        rng = np.random.default_rng(51)
        left, _ = np.linalg.qr(
            rng.normal(size=(9, 3)) + 1.0j * rng.normal(size=(9, 3))
        )
        right, _ = np.linalg.qr(
            rng.normal(size=(9, 3)) + 1.0j * rng.normal(size=(9, 3))
        )
        dense = np.linalg.norm(
            left @ left.conj().T - right @ right.conj().T, 2
        )
        self.assertAlmostEqual(Z2.projector_distance(left, right), dense, places=12)

    def test_skinny_projector_distance_handles_large_plane_wave_basis(self):
        # A dense 10000x10000 complex projector would require about 1.6 GB.
        # The production implementation must retain only the 10000x4 frames.
        rng = np.random.default_rng(52)
        left, _ = np.linalg.qr(
            rng.normal(size=(10_000, 4)) + 1.0j * rng.normal(size=(10_000, 4))
        )
        right = left.copy()
        right[:, 0] += 1.0e-4 * (
            rng.normal(size=10_000) + 1.0j * rng.normal(size=10_000)
        )
        distance = Z2.projector_distance(left, right)
        self.assertTrue(np.isfinite(distance))
        self.assertGreater(distance, 0.0)
        self.assertLess(distance, 0.1)

    def test_unordered_circular_set_matcher_handles_periodic_boundary(self):
        left = (0.99, 0.24, 0.51)
        right = (0.01, 0.50, 0.25)
        self.assertAlmostEqual(Z2.circular_set_distance(left, right), 0.02)

    def test_z2pack_matcher_cuts_at_union_largest_gap(self):
        left = np.asarray((0.02, 0.31, 0.67, 0.96))
        right = np.asarray((0.99, 0.29, 0.65, 0.04))
        cut = Z2.largest_circular_gap(np.r_[left, right])[0]
        displacement = np.abs(
            np.sort(np.mod(left - cut, 1.0))
            - np.sort(np.mod(right - cut, 1.0))
        )
        expected = float(np.max(np.minimum(displacement, 1.0 - displacement)))
        self.assertAlmostEqual(Z2.z2pack_max_move(left, right), expected)

        # A documented regression: circular bottleneck matching and Z2Pack's
        # union-gap-cut matching are not interchangeable for arbitrary sets.
        regression_left = (0.02831967, 0.12428328)
        regression_right = (0.53649929, 0.58338584)
        self.assertNotAlmostEqual(
            Z2.circular_set_distance(regression_left, regression_right),
            Z2.z2pack_max_move(regression_left, regression_right),
        )

    def test_z2pack_gap_ratio_checks_both_neighbor_directions(self):
        left = np.asarray((0.10, 0.40, 0.70))
        right = np.asarray((0.12, 0.39, 0.73))
        lc, lg, _ = Z2.largest_circular_gap(left)
        rc, rg, _ = Z2.largest_circular_gap(right)
        expected = min(
            np.min(np.abs(right - lc)) / lg,
            np.min(np.abs(left - rc)) / rg,
        )
        self.assertAlmostEqual(Z2.z2pack_gap_ratio(left, right), expected)

    def test_z2pack_fixed_grid_threshold_equalities_fail_strictly(self):
        self.assertFalse(Z2._check("pos", 0.01, 0.01, "<")["pass"])
        self.assertFalse(Z2._check("move", 0.30, 0.30, "<")["pass"])
        self.assertFalse(Z2._check("gap", 0.30, 0.30, ">")["pass"])

        result = Z2.analyze_analytic_frames(*doubled_qwz(-3.0, 32))
        checks = {item["name"]: item for item in result["checks"]}
        self.assertEqual(checks["x_loop_wcc_position_converged"]["relation"], "<")
        self.assertEqual(checks["x_loop_wcc_move_neighbor_converged"]["relation"], "<")
        self.assertEqual(checks["x_loop_wcc_gap_neighbor_converged"]["relation"], ">")

    def test_largest_gap_parity_uses_official_linear_interval(self):
        # Largest-gap centers are 0.1 then 0.9.  The right row has three WCCs
        # strictly inside (0.1, 0.9), hence odd parity.  Replacing this by the
        # shortest periodic arc would incorrectly return zero.
        wcc = np.asarray(((0.3, 0.5, 0.7, 0.9), (0.1, 0.3, 0.5, 0.7)))
        parity, centers, _ = Z2.z2pack_largest_gap_parity(wcc)
        np.testing.assert_allclose(centers, (0.1, 0.9), atol=1e-14)
        self.assertEqual(parity, 1)

    def test_doubled_qwz_phase_diagram_is_0_1_1_0_on_both_axes(self):
        observed = []
        for mass in (-3.0, -1.0, 1.0, 3.0):
            result = Z2.analyze_analytic_frames(*doubled_qwz(mass, 32))
            self.assertTrue(result["valid"], (mass, result["failed_guards"]))
            self.assertEqual(
                result["axes"]["x"]["candidate_z2"],
                result["axes"]["y"]["candidate_z2"],
            )
            observed.append(result["z2"])
        self.assertEqual(observed, [0, 1, 1, 0])

    def test_random_occupied_gauge_does_not_change_topological_phase(self):
        frames, energies, theta = doubled_qwz(-1.0, 32)
        rng = np.random.default_rng(20260831)
        transformed = frames.copy()
        for ix in range(frames.shape[0]):
            for iy in range(frames.shape[1]):
                transformed[ix, iy] = frames[ix, iy] @ random_unitary(2, rng)
        result = Z2.analyze_analytic_frames(transformed, energies, theta)
        self.assertTrue(result["valid"], result["failed_guards"])
        self.assertEqual(result["z2"], 1)
        self.assertLess(result["metrics"]["max_projector_tr_residual"], 1e-12)

    def test_underresolved_4x4_candidate_is_not_reported_as_z2(self):
        result = Z2.analyze_analytic_frames(*doubled_qwz(-1.0, 4))
        self.assertEqual(result["candidate_z2"], 1)
        self.assertFalse(result["valid"])
        self.assertIsNone(result["z2"])
        self.assertTrue(
            any("converged" in guard for guard in result["failed_guards"]),
            result["failed_guards"],
        )

    def test_single_axis_run_remains_diagnostic_only(self):
        result = Z2.analyze_analytic_frames(
            *doubled_qwz(-3.0, 32), axes=(0,)
        )
        self.assertFalse(result["valid"])
        self.assertEqual(result["candidate_z2"], 0)
        self.assertIsNone(result["z2"])
        self.assertIn("both_loop_axes_evaluated", result["failed_guards"])

    def test_non_time_reversal_invariant_kz_plane_is_rejected(self):
        wavecar = SimpleNamespace(
            header=SimpleNamespace(ispin=1, nbands=4), spinor_components=2
        )
        make_grid = lambda kz: Z2.Grid(
            nx=4,
            ny=4,
            index=np.arange(16).reshape(4, 4),
            offset=np.asarray((0.0, 0.0, kz)),
        )
        with self.assertRaisesRegex(ValueError, "kz=0"):
            Z2._validate_wavecar_inputs(wavecar, make_grid(0.1), occupied=2)
        # kz=1 differs from Gamma by a reciprocal lattice vector.
        Z2._validate_wavecar_inputs(wavecar, make_grid(1.0), occupied=2)

    def test_breaking_projector_time_reversal_fails_hard_gate(self):
        frames, energies, theta = doubled_qwz(-3.0, 32)
        broken = frames.copy()
        vector = broken[1, 2, :, 0]
        vector = vector + 0.2 * np.asarray((0.0, 0.0, 1.0, 0.0))
        broken[1, 2, :, 0] = vector / np.linalg.norm(vector)
        result = Z2.analyze_analytic_frames(broken, energies, theta)
        self.assertFalse(result["valid"])
        self.assertIn("occupied_projector_time_reversal", result["failed_guards"])

    def test_output_collision_never_deletes_input(self):
        with tempfile.TemporaryDirectory() as directory:
            wavecar = Path(directory) / "WAVECAR"
            payload = b"not-a-wavecar-but-must-survive"
            wavecar.write_bytes(payload)
            status = Z2.main(
                (
                    str(wavecar),
                    "--nx",
                    "4",
                    "--ny",
                    "4",
                    "--occupied-bands",
                    "2",
                    "--output-dir",
                    directory,
                    "--json-name",
                    "WAVECAR",
                )
            )
            self.assertEqual(status, 1)
            self.assertEqual(wavecar.read_bytes(), payload)

    def test_hard_error_removes_stale_planned_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            (output / Z2.DEFAULT_CSV).write_text("stale", encoding="utf-8")
            (output / Z2.DEFAULT_JSON).write_text("stale", encoding="utf-8")
            status = Z2.main(
                (
                    str(output / "missing-WAVECAR"),
                    "--nx",
                    "4",
                    "--ny",
                    "4",
                    "--occupied-bands",
                    "2",
                    "--output-dir",
                    str(output),
                )
            )
            self.assertEqual(status, 1)
            self.assertFalse((output / Z2.DEFAULT_CSV).exists())
            self.assertFalse((output / Z2.DEFAULT_JSON).exists())

    def test_public_json_writer_is_strict_and_candidate_is_separate(self):
        payload = {
            "valid": False,
            "z2": None,
            "candidate_z2": 0,
            "metric": np.inf,
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "diagnostics.json"
            Z2._atomic_write_json(path, payload)
            raw = path.read_text(encoding="utf-8")
            self.assertNotIn("Infinity", raw)
            parsed = json.loads(raw, parse_constant=lambda value: self.fail(value))
            self.assertIsNone(parsed["z2"])
            self.assertEqual(parsed["candidate_z2"], 0)
            self.assertIsNone(parsed["metric"])


if __name__ == "__main__":
    unittest.main()
