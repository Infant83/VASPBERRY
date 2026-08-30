"""Regression tests for all-band cumulative occupied-subspace transport."""

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "wavecar_fukui_full_test", ROOT / "tools" / "wavecar_fukui.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def stable_links(shape=(2, 2), minimum=1.0):
    ones = np.ones(shape)
    return MODULE.LinkSet(
        phase=np.zeros(shape),
        logabs=np.zeros(shape),
        min_singular=ones * minimum,
        max_singular=ones,
        coverage_left=ones,
        coverage_right=ones,
    )


class FullCumulativeFormulaTests(unittest.TestCase):
    def test_vertex_counts_select_cumulative_subspaces(self):
        phases = np.asarray([[[0.0]], [[1.0]], [[2.0]]])
        energies = np.asarray(
            [
                [[[0.0, 1.0, 2.0, 3.0]]],
                [[[0.5, 1.5, 2.5, 3.5]]],
            ]
        )
        effective, counts = MODULE.full_cumulative_t0_effective_phi(
            phases, energies, 1.25
        )
        self.assertAlmostEqual(float(effective[0, 0]), (2.0 + 1.0) / 4.0)
        self.assertEqual(counts[0, 0].tolist(), [2, 1, 0, 0])

    def test_exactly_degenerate_group_enters_together(self):
        phases = np.asarray([[[0.0]], [[100.0]], [[2.0]]])
        energies = np.ones((2, 1, 1, 4))
        effective, counts = MODULE.full_cumulative_t0_effective_phi(
            phases, energies, 1.0
        )
        self.assertEqual(counts[0, 0].tolist(), [2, 2, 2, 2])
        self.assertEqual(float(effective[0, 0]), 2.0)

    def test_formula_rejects_inverted_energy_order_and_nonzero_empty_phase(self):
        phases = np.zeros((3, 1, 1))
        energies = np.asarray(
            [
                [[[1.0, 1.0, 1.0, 1.0]]],
                [[[0.0, 2.0, 2.0, 2.0]]],
            ]
        )
        with self.assertRaisesRegex(ValueError, "nondecreasing"):
            MODULE.full_cumulative_t0_effective_phi(phases, energies, 0.5)
        phases[0, 0, 0] = 1.0
        energies[1, 0, 0, 0] = 2.0
        with self.assertRaisesRegex(ValueError, "empty occupied subspace"):
            MODULE.full_cumulative_t0_effective_phi(phases, energies, 0.5)


class BatchedCumulativeMapTests(unittest.TestCase):
    def test_batched_leading_subspaces_match_individual_map_path(self):
        grid = MODULE.Grid(
            2, 2, np.asarray([[0, 1], [2, 3]]), np.zeros(3)
        )
        fake = SimpleNamespace(
            header=SimpleNamespace(nbands=2),
            kpoints=np.asarray(
                [
                    (0.0, 0.0, 0.0),
                    (0.0, 0.5, 0.0),
                    (0.5, 0.0, 0.0),
                    (0.5, 0.5, 0.0),
                ]
            ),
        )
        matrix = np.asarray(
            [[np.exp(0.1j), 0.03j], [0.02, 0.95 * np.exp(-0.2j)]],
            dtype=np.complex128,
        )

        def fixed_link(*args, **_kwargs):
            band_count = len(args[3])
            return matrix[:band_count, :band_count].copy(), 0.98, 0.97

        with mock.patch.object(MODULE, "link_matrix", side_effect=fixed_link):
            batched = MODULE.compute_cumulative_band_maps(fake, grid, 2)
            self.assertEqual(MODULE.link_matrix.call_count, 2 * grid.nx * grid.ny)
            for band in (1, 2):
                direct = MODULE.compute_band_map(fake, grid, list(range(1, band + 1)))
                for got, expected in zip(batched[band], direct):
                    if isinstance(got, np.ndarray):
                        np.testing.assert_allclose(got, expected)
                    else:
                        np.testing.assert_allclose(got.phase, expected.phase)
                        np.testing.assert_allclose(
                            got.min_singular, expected.min_singular
                        )
                        np.testing.assert_allclose(
                            got.coverage_left, expected.coverage_left
                        )


class FullCumulativeWriterTests(unittest.TestCase):
    def setUp(self):
        self.grid = MODULE.Grid(
            2, 2, np.asarray([[0, 1], [2, 3]]), np.zeros(3)
        )
        self.wavecar = SimpleNamespace(
            header=SimpleNamespace(nbands=3, reciprocal=np.eye(3)),
            energies=np.tile(np.asarray([[-1.0, 0.0, 10.0]]), (4, 1)),
            kpoints=np.asarray(
                [
                    (0.0, 0.0, 0.0),
                    (0.0, 0.5, 0.0),
                    (0.5, 0.0, 0.0),
                    (0.5, 0.5, 0.0),
                ]
            ),
        )
        stable = stable_links()
        self.maps = {
            1: (np.zeros((2, 2)), stable, stable),
            2: (np.zeros((2, 2)), stable, stable),
        }

    def call_writer(self, output, *, allow_invalid=False, max_phi=0.8 * np.pi):
        with mock.patch.object(
            MODULE, "compute_cumulative_band_maps", return_value=self.maps
        ):
            MODULE.write_full_t0_transport(
                output,
                self.wavecar,
                self.grid,
                2,
                -2.0,
                0.5,
                3,
                np.asarray((0.0, 0.0, 0.0)),
                np.asarray((0.25, 0.25, 0.0)),
                None,
                1.0e-3,
                1.0e-5,
                0.9,
                max_phi,
                allow_invalid,
            )

    def test_empty_through_fully_occupied_window_is_written(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            self.call_writer(output)
            with (output / "transport_full_t0.csv").open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 3)
            self.assertEqual(rows[0]["occupied_band_count_max"], "0")
            self.assertEqual(rows[0]["worst_active_link_sv"], "")
            self.assertNotIn("inf", (output / "transport_full_t0.csv").read_text())
            self.assertEqual(rows[-1]["occupied_band_count_min"], "2")
            self.assertIn(
                "sigma_valley_contrast_relative_to_full_max_bundle_e2_over_h",
                rows[-1],
            )
            self.assertEqual(
                float(
                    rows[-1][
                        "sigma_valley_contrast_relative_to_full_max_bundle_e2_over_h"
                    ]
                ),
                0.0,
            )
            self.assertTrue(all(row["quality"] == "PASS" for row in rows))
            diagnostics = json.loads(
                (output / "transport_full_t0_diagnostics.json").read_text()
            )
            self.assertTrue(diagnostics["validated"])
            self.assertEqual(diagnostics["sentinel_band"], 3)
            self.assertTrue(diagnostics["positive_indirect_gap_above_max_bundle"])
            self.assertTrue(diagnostics["max_bundle_plateau_in_requested_range"])
            self.assertEqual(diagnostics["max_bundle_plateau_sample_count"], 1)
            self.assertEqual(
                diagnostics["max_bundle_plateau_validated_sample_count"], 1
            )
            self.assertTrue(diagnostics["max_bundle_plateau_quality_pass"])
            self.assertEqual(diagnostics["max_bundle_plateau_max_chern_error"], 0.0)
            self.assertTrue(diagnostics["max_bundle_reference_validated"])
            self.assertTrue(diagnostics["max_bundle_quality"]["quality_pass"])
            self.assertEqual(
                diagnostics["max_bundle_chern_nearest_integer_residual"], 0.0
            )
            figure = output / "full.png"
            MODULE.plot_full_transport_csv(
                output / "transport_full_t0.csv",
                figure,
                output / "transport_full_t0_diagnostics.json",
            )
            self.assertTrue(figure.exists())
            self.assertGreater(figure.stat().st_size, 1000)

    def test_full_window_rejects_occupied_sentinel(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "reaches sentinel band 3"):
                with mock.patch.object(
                    MODULE, "compute_cumulative_band_maps", return_value=self.maps
                ):
                    MODULE.write_full_t0_transport(
                        Path(directory), self.wavecar, self.grid, 2,
                        -2.0, 10.0, 3,
                        np.zeros(3), np.asarray((0.25, 0.25, 0.0)), None,
                        1.0e-3, 1.0e-5, 0.9, 0.8 * np.pi, False,
                    )

    def test_invalid_reason_is_recorded_and_csv_refused(self):
        stable = stable_links()
        self.maps = {
            1: (np.zeros((2, 2)), stable, stable),
            2: (np.full((2, 2), 0.9 * np.pi), stable, stable),
        }
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            with self.assertRaisesRegex(RuntimeError, "validation failed"):
                self.call_writer(output)
            self.assertFalse((output / "transport_full_t0.csv").exists())
            diagnostics = json.loads(
                (output / "transport_full_t0_diagnostics.json").read_text()
            )
            self.assertGreater(diagnostics["invalid_mu_count"], 0)
            self.assertGreater(
                diagnostics["invalid_mu"][-1]["phase_failure_cell_count"], 0
            )

    def test_event_audit_detects_narrow_invalid_interval_between_mu_samples(self):
        self.wavecar.energies = np.tile(
            np.asarray([[0.0, 1.0e-7, 10.0]]), (4, 1)
        )
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            with self.assertRaisesRegex(RuntimeError, "occupation-event subinterval"):
                self.call_writer(output)
            diagnostics = json.loads(
                (output / "transport_full_t0_diagnostics.json").read_text()
            )
            self.assertTrue(diagnostics["requested_mu_points_validated"])
            self.assertFalse(diagnostics["continuous_mu_range_validated"])
            self.assertGreater(diagnostics["continuous_invalid_event_count"], 0)
            event = diagnostics["continuous_invalid_events"][0]
            self.assertIn(event["vertex"], range(4))
            self.assertTrue(event["lower_bound_inclusive"])
            self.assertFalse(event["upper_bound_inclusive"])
            self.assertEqual(
                event["occupation_interval_convention"], "[E_n,E_(n+1))"
            )
            intervals = diagnostics["merged_invalid_mu_intervals_ev"]
            self.assertTrue(any(lo <= 0.0 and hi >= 1.0e-7 for lo, hi in intervals))

    def test_exact_degenerate_pair_skips_unsafe_partial_subspace(self):
        self.wavecar.energies = np.tile(
            np.asarray([[0.0, 0.0, 10.0]]), (4, 1)
        )
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            self.call_writer(output)
            diagnostics = json.loads(
                (output / "transport_full_t0_diagnostics.json").read_text()
            )
            self.assertTrue(diagnostics["validated"])
            self.assertEqual(diagnostics["continuous_invalid_event_count"], 0)

    def test_opposite_branch_safe_phases_are_not_rejected_by_span(self):
        self.wavecar.energies = np.asarray(
            [
                [0.0, 0.1, 10.0],
                [0.0, 1.0, 10.0],
                [0.0, 1.0, 10.0],
                [0.0, 1.0, 10.0],
            ]
        )
        stable = stable_links()
        self.maps = {
            1: (np.full((2, 2), 0.7 * np.pi), stable, stable),
            2: (np.full((2, 2), -0.7 * np.pi), stable, stable),
        }
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            self.call_writer(output)
            diagnostics = json.loads(
                (output / "transport_full_t0_diagnostics.json").read_text()
            )
            self.assertTrue(diagnostics["validated"])
            self.assertEqual(diagnostics["continuous_invalid_event_count"], 0)
            header = (output / "transport_full_t0.csv").read_text().splitlines()[0]
            self.assertNotIn("phase_span", header)

    def test_active_singular_rank_is_blank_and_explicitly_invalid(self):
        singular = stable_links(minimum=0.0)
        singular.phase[:] = np.nan
        singular.logabs[:] = -np.inf
        self.maps = {
            1: (np.full((2, 2), np.nan), singular, singular),
            2: self.maps[2],
        }
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            self.call_writer(output, allow_invalid=True)
            diagnostics = json.loads(
                (output / "transport_full_t0_diagnostics.json").read_text()
            )
            self.assertFalse(diagnostics["validated"])
            self.assertTrue(
                any(
                    event["nonfinite_or_singular_failure"]
                    and event["occupied_band_count"] == 1
                    for event in diagnostics["continuous_invalid_events"]
                )
            )
            self.assertEqual(
                diagnostics["cumulative_map_numeric_issues"][0][
                    "occupied_band_count"
                ],
                1,
            )
            csv_text = (output / "transport_full_t0.csv").read_text()
            self.assertNotIn("nan", csv_text.lower())
            rows = list(csv.DictReader(csv_text.splitlines()))
            active_row = next(
                row for row in rows if row["occupied_band_count_min"] == "1"
            )
            self.assertEqual(active_row["quality"], "INVALID")
            self.assertGreater(
                int(active_row["nonfinite_or_singular_active_cell_count"]), 0
            )
            self.assertEqual(active_row["sigma_total_e2_over_h"], "")
            figure = output / "singular-diagnostic.png"
            MODULE.plot_full_transport_csv(
                output / "transport_full_t0.csv",
                figure,
                output / "transport_full_t0_diagnostics.json",
            )
            self.assertTrue(figure.exists())

    def test_inactive_singular_rank_is_tolerated(self):
        singular = stable_links(minimum=0.0)
        singular.phase[:] = np.nan
        singular.logabs[:] = -np.inf
        self.maps = {
            1: (np.full((2, 2), np.nan), singular, singular),
            2: self.maps[2],
        }
        self.wavecar.energies = np.tile(
            np.asarray([[0.0, 0.0, 10.0]]), (4, 1)
        )
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            self.call_writer(output)
            diagnostics = json.loads(
                (output / "transport_full_t0_diagnostics.json").read_text()
            )
            self.assertTrue(diagnostics["validated"])
            self.assertEqual(diagnostics["continuous_invalid_event_count"], 0)
            self.assertEqual(
                diagnostics["cumulative_map_numeric_issues"][0][
                    "occupied_band_count"
                ],
                1,
            )
            self.assertNotIn(
                "nan", (output / "transport_full_t0.csv").read_text().lower()
            )

    def test_undefined_max_bundle_baseline_is_always_refused(self):
        singular = stable_links(minimum=0.0)
        singular.phase[:] = np.nan
        singular.logabs[:] = -np.inf
        self.maps = {
            1: self.maps[1],
            2: (np.full((2, 2), np.nan), singular, singular),
        }
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            with self.assertRaisesRegex(RuntimeError, "MAX_BAND reference phase"):
                self.call_writer(output, allow_invalid=True)
            diagnostics = json.loads(
                (output / "transport_full_t0_diagnostics.json").read_text()
            )
            self.assertFalse(diagnostics["validated"])
            self.assertFalse(diagnostics["max_bundle_quality"]["phase_finite"])
            self.assertFalse((output / "transport_full_t0.csv").exists())


class FullCumulativeCliTests(unittest.TestCase):
    def test_full_and_narrow_transport_options_are_mutually_exclusive(self):
        completed = __import__("subprocess").run(
            [
                sys.executable,
                str(ROOT / "tools" / "wavecar_fukui.py"),
                "missing-WAVECAR",
                "--nx", "2", "--ny", "2",
                "--energy-band", "2",
                "--output-dir", "unused",
                "--transport-t0",
                "--transport-full-t0", "2",
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("mutually exclusive", completed.stderr)
        self.assertNotIn("Traceback", completed.stderr)

    def test_full_only_mode_does_not_validate_irrelevant_energy_band(self):
        reciprocal = np.eye(3)
        fake = SimpleNamespace(
            header=SimpleNamespace(
                logical_recl=256,
                stride_bytes=256,
                ispin=1,
                rtag=45200,
                nkpoints=4,
                nbands=3,
                encut_ev=400.0,
                lattice=np.eye(3),
                reciprocal=reciprocal,
            ),
            spin=1,
            spinor_components=1,
            nbmax=(1, 1, 1),
            kpoints=np.asarray(
                [
                    (0.0, 0.0, 0.0),
                    (0.0, 0.5, 0.0),
                    (0.5, 0.0, 0.0),
                    (0.5, 0.5, 0.0),
                ]
            ),
            energies=np.tile(np.asarray([[-1.0, 0.0, 10.0]]), (4, 1)),
        )
        grid = MODULE.Grid(
            2, 2, np.asarray([[0, 1], [2, 3]]), np.zeros(3)
        )
        with tempfile.TemporaryDirectory() as directory, mock.patch.object(
            sys,
            "argv",
            [
                "wavecar_fukui.py", "missing-WAVECAR",
                "--nx", "2", "--ny", "2",
                "--energy-band", "999",
                "--output-dir", directory,
                "--transport-full-t0", "2",
                "--mu-min", "-2", "--mu-max", "0.5", "--mu-num", "3",
                "--valley-k", "0,0,0",
                "--valley-kp", "0.25,0.25,0",
            ],
        ), mock.patch.object(MODULE, "Wavecar", return_value=fake), mock.patch.object(
            MODULE, "infer_uniform_grid", return_value=grid
        ), mock.patch.object(MODULE, "write_full_t0_transport") as writer, mock.patch(
            "builtins.print"
        ):
            MODULE.main()
            writer.assert_called_once()
            diagnostics = json.loads(
                (Path(directory) / "diagnostics.json").read_text()
            )
            self.assertIsNone(diagnostics["energy_diagnostics"])


if __name__ == "__main__":
    unittest.main()
