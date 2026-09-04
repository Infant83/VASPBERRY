"""Regression tests for folded reciprocal-space time reversal at TRIM."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import wavecar_z2 as z2  # noqa: E402
from wavecar_fukui import Grid  # noqa: E402


TRIM_NAMES = ("Gamma", "M2", "M1", "M3")
KPOINTS = np.asarray(
    (
        (0.0, 0.0, 0.0),
        (0.0, 0.5, 0.0),
        (0.5, 0.0, 0.0),
        (0.5, 0.5, 0.0),
    )
)
GRID = Grid(
    nx=2,
    ny=2,
    index=np.asarray(((0, 1), (2, 3)), dtype=np.int64),
    offset=np.zeros(3),
)
EXPECTED_SHIFT = {
    "Gamma": (0, 0, 0),
    "M1": (-1, 0, 0),
    "M2": (0, -1, 0),
    "M3": (-1, -1, 0),
}


def _theta_spinor(values: np.ndarray) -> np.ndarray:
    return np.asarray((-values[1].conjugate(), values[0].conjugate()))


class SyntheticTrimWavecar:
    spinor_components = 2
    header = SimpleNamespace(ispin=1, nbands=3)

    def __init__(
        self,
        *,
        legacy_indices: frozenset[int] = frozenset(),
        reverse_g_order: bool = False,
        missing_target_index: int | None = None,
    ) -> None:
        self.kpoints = KPOINTS.copy()
        self.energies = np.tile(np.asarray((-1.0, -1.0, 1.0)), (4, 1))
        self._g: list[np.ndarray] = []
        self._coefficients: list[np.ndarray] = []
        spinor = np.asarray((1.0 + 2.0j, -0.3 + 0.7j), dtype=np.complex128)
        spinor /= np.linalg.norm(spinor)
        for index, kpoint in enumerate(self.kpoints):
            shift = np.rint(-2.0 * kpoint).astype(np.int64)
            origin = np.zeros(3, dtype=np.int64)
            if np.any(shift):
                gvectors = np.stack((origin, shift))
            else:
                gvectors = origin[None, :]
            coefficients = np.zeros(
                (2, 2, gvectors.shape[0]), dtype=np.complex128
            )
            origin_index = int(
                np.flatnonzero(np.all(gvectors == origin, axis=1))[0]
            )
            exact_target_index = int(
                np.flatnonzero(np.all(gvectors == shift, axis=1))[0]
            )
            coefficients[0, :, origin_index] = spinor
            partner_index = (
                origin_index if index in legacy_indices else exact_target_index
            )
            coefficients[1, :, partner_index] = _theta_spinor(spinor)
            if index == missing_target_index:
                keep = np.arange(gvectors.shape[0]) != exact_target_index
                gvectors = gvectors[keep]
                coefficients = coefficients[:, :, keep]
            if reverse_g_order and gvectors.shape[0] > 1:
                permutation = np.arange(gvectors.shape[0] - 1, -1, -1)
                gvectors = gvectors[permutation]
                coefficients = coefficients[:, :, permutation]
            self._g.append(gvectors)
            self._coefficients.append(coefficients)

    def g_vectors(self, ik: int) -> np.ndarray:
        return self._g[ik]

    def coefficients(self, ik: int, bands: list[int]) -> np.ndarray:
        return self._coefficients[ik][np.asarray(bands) - 1]


class TrimReciprocalMappingTests(unittest.TestCase):
    def test_complex64_raw_norm_is_invariant_to_g_vector_order(self) -> None:
        rng = np.random.default_rng(20260904)
        coefficients = (
            rng.standard_normal((2, 2, 10_000))
            + 1j * rng.standard_normal((2, 2, 10_000))
        ).astype(np.complex64)
        coefficients *= np.geomspace(1.0e-5, 1.0, 10_000).astype(np.float32)
        permutation = rng.permutation(coefficients.shape[-1])

        reference = z2._raw_band_norm(coefficients)
        permuted = z2._raw_band_norm(coefficients[:, :, permutation])

        self.assertEqual(reference.dtype, np.dtype(np.float64))
        np.testing.assert_allclose(permuted, reference, rtol=1.0e-14, atol=0.0)
        self.assertGreaterEqual(float(np.min(permuted / reference)), 0.999999)

    def test_exact_mapping_at_gamma_and_all_three_m_points(self) -> None:
        wavecar = SyntheticTrimWavecar()
        for index, name in enumerate(TRIM_NAMES):
            with self.subTest(trim=name):
                kpoint = wavecar.kpoints[index]
                target_g = np.rint(-2.0 * kpoint).astype(np.int64)
                self.assertEqual(tuple(target_g), EXPECTED_SHIFT[name])
                np.testing.assert_allclose(
                    kpoint + target_g + kpoint,
                    np.zeros(3),
                    atol=1.0e-14,
                    rtol=0.0,
                )
        diagnostics = z2.wavecar_time_reversal_diagnostics(
            wavecar, GRID, [1, 2]
        )
        self.assertLess(diagnostics["max_projector_tr_residual"], 1.0e-12)
        self.assertEqual(diagnostics["min_tr_g_basis_coverage"], 1.0)
        self.assertEqual(diagnostics["min_tr_raw_norm_coverage"], 1.0)

    def test_shift_free_legacy_mapping_fails_at_each_nonzero_trim(self) -> None:
        gamma = z2.wavecar_time_reversal_diagnostics(
            SyntheticTrimWavecar(legacy_indices=frozenset((0,))),
            GRID,
            [1, 2],
        )
        self.assertLess(gamma["max_projector_tr_residual"], 1.0e-12)
        for index, name in ((1, "M2"), (2, "M1"), (3, "M3")):
            with self.subTest(trim=name):
                diagnostics = z2.wavecar_time_reversal_diagnostics(
                    SyntheticTrimWavecar(
                        legacy_indices=frozenset((index,))
                    ),
                    GRID,
                    [1, 2],
                )
                self.assertGreater(
                    diagnostics["max_projector_tr_residual"], 0.5
                )

    def test_mapping_is_invariant_to_g_vector_order(self) -> None:
        reference = z2.wavecar_time_reversal_diagnostics(
            SyntheticTrimWavecar(), GRID, [1, 2]
        )
        permuted = z2.wavecar_time_reversal_diagnostics(
            SyntheticTrimWavecar(reverse_g_order=True), GRID, [1, 2]
        )
        self.assertLess(permuted["max_projector_tr_residual"], 1.0e-12)
        for metric in (
            "max_projector_tr_residual",
            "min_tr_g_basis_coverage",
            "min_tr_raw_norm_coverage",
        ):
            self.assertAlmostEqual(
                permuted[metric], reference[metric], places=14
            )

    def test_missing_folded_target_g_fails_hard(self) -> None:
        with self.assertRaisesRegex(
            ValueError,
            "time-reversal projector comparison has no common plane waves",
        ):
            z2.wavecar_time_reversal_diagnostics(
                SyntheticTrimWavecar(missing_target_index=2),
                GRID,
                [1, 2],
            )


if __name__ == "__main__":
    unittest.main()
