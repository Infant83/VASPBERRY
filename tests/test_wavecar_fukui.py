import argparse
import importlib.util
import json
import math
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "wavecar_fukui", ROOT / "tools" / "wavecar_fukui.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class WavecarFukuiUnitTests(unittest.TestCase):
    def test_vasp_signed_integer_order(self):
        self.assertEqual(
            MODULE.signed_vasp_sequence(3), [0, 1, 2, 3, -3, -2, -1]
        )

    def test_legacy_four_byte_recl_is_detected(self):
        logical_recl = 256
        stride = 4 * logical_recl
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "WAVECAR"
            data = bytearray(4 * stride)
            struct.pack_into("<3d", data, 0, logical_recl, 1.0, 45200.0)
            lattice = (3.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 15.0)
            struct.pack_into("<12d", data, stride, 1.0, 1.0, 400.0, *lattice)
            struct.pack_into(
                "<7d", data, 2 * stride, 2.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0
            )
            path.write_bytes(data)
            parsed = MODULE.Wavecar(path, spinor_components=2)
            self.assertEqual(parsed.header.logical_recl, logical_recl)
            self.assertEqual(parsed.header.stride_bytes, stride)
            self.assertEqual(parsed.header.nkpoints, 1)
            self.assertEqual(parsed.header.nbands, 1)

    def test_ispin_two_requires_both_spin_blocks_and_reads_second_channel(self):
        stride = 256
        lattice = (3.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 15.0)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "WAVECAR"
            truncated = bytearray(4 * stride)
            struct.pack_into("<3d", truncated, 0, stride, 2.0, 45200.0)
            struct.pack_into("<12d", truncated, stride, 1.0, 1.0, 400.0, *lattice)
            struct.pack_into(
                "<7d", truncated, 2 * stride,
                1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0,
            )
            path.write_bytes(truncated)
            with self.assertRaisesRegex(ValueError, "physical RECL stride"):
                MODULE.Wavecar(path, spinor_components=1, spin=1)

            complete = bytearray(6 * stride)
            complete[: len(truncated)] = truncated
            struct.pack_into(
                "<7d", complete, 4 * stride,
                1.0, 0.0, 0.0, 0.0, 2.5, 0.0, 0.0,
            )
            path.write_bytes(complete)
            parsed = MODULE.Wavecar(path, spinor_components=1, spin=2)
            self.assertEqual(parsed.record_number(0, spin=2), 5)
            self.assertAlmostEqual(float(parsed.energies[0, 0]), 2.5)

    def test_invalid_ispin_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "WAVECAR"
            data = bytearray(512)
            struct.pack_into("<3d", data, 0, 256.0, 3.0, 45200.0)
            path.write_bytes(data)
            with self.assertRaisesRegex(ValueError, "invalid ISPIN"):
                MODULE.Wavecar(path, spinor_components=1)

    def test_wavecar_rejects_nonpositive_nplane_and_nonfinite_headers(self):
        stride = 256
        lattice = (3.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 15.0)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "WAVECAR"
            data = bytearray(4 * stride)
            struct.pack_into("<3d", data, 0, stride, 1.0, 45200.0)
            struct.pack_into("<12d", data, stride, 1.0, 1.0, 400.0, *lattice)
            struct.pack_into(
                "<7d", data, 2 * stride,
                0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0,
            )
            path.write_bytes(data)
            with self.assertRaisesRegex(ValueError, "NPLANE must be positive"):
                MODULE.Wavecar(path, spinor_components=1)

            struct.pack_into(
                "<7d", data, 2 * stride,
                1.0, 0.0, 0.0, 0.0, np.nan, 0.0, 1.0,
            )
            path.write_bytes(data)
            with self.assertRaisesRegex(ValueError, "band energies must all be finite"):
                MODULE.Wavecar(path, spinor_components=1)

    def test_nonfinite_coefficient_record_is_rejected(self):
        stride = 256
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "WAVECAR"
            data = bytearray(3 * stride)
            struct.pack_into("<ff", data, 2 * stride, np.nan, 0.0)
            path.write_bytes(data)
            wavecar = object.__new__(MODULE.Wavecar)
            wavecar.path = path
            wavecar.spinor_components = 1
            wavecar.spin = 1
            wavecar.header = SimpleNamespace(stride_bytes=stride)
            wavecar.nplane = np.asarray([1])
            wavecar.g_vectors = lambda _ik: np.zeros((1, 3), dtype=int)
            wavecar.record_number = lambda _ik, _band, _spin: 3
            with self.assertRaisesRegex(ValueError, "nonfinite coefficient record"):
                wavecar.coefficients(0, [1])

    def test_json_writer_replaces_nonfinite_values_and_emits_strict_json(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "strict.json"
            MODULE.write_strict_json(
                path, {"finite": 1.0, "infinite": np.inf, "nan": np.nan}
            )
            raw = path.read_text(encoding="utf-8")
            self.assertNotIn("Infinity", raw)
            self.assertNotIn("NaN", raw)
            parsed = json.loads(raw, parse_constant=lambda token: self.fail(token))
            self.assertIsNone(parsed["infinite"])
            self.assertIsNone(parsed["nan"])

    def test_uniform_signed_mesh_and_neighbors(self):
        kpoints = np.asarray(
            [
                (0.0, 0.0, 0.0),
                (0.0, -0.5, 0.0),
                (-0.5, 0.0, 0.0),
                (-0.5, -0.5, 0.0),
            ]
        )
        grid = MODULE.infer_uniform_grid(kpoints, 2, 2)
        self.assertEqual(grid.index.tolist(), [[0, 1], [2, 3]])
        self.assertEqual(MODULE.vertex_indices(grid, 1, 1), (3, 1, 0, 2))

    def test_uniform_shifted_monkhorst_mesh_preserves_offset(self):
        kpoints = np.asarray(
            [
                (0.125, 0.125, 0.0),
                (0.125, -0.375, 0.0),
                (-0.375, 0.125, 0.0),
                (-0.375, -0.375, 0.0),
            ]
        )
        grid = MODULE.infer_uniform_grid(kpoints, 2, 2)
        np.testing.assert_allclose(grid.offset, (0.125, 0.125, 0.0))
        self.assertEqual(grid.index.tolist(), [[0, 1], [2, 3]])

    def test_uniform_grid_requires_two_points_per_plaquette_axis(self):
        with self.assertRaisesRegex(ValueError, "both be at least 2"):
            MODULE.infer_uniform_grid(np.zeros((2, 3)), 1, 2)
        kpoints = np.asarray(
            [(0.0, 0.0, 0.0), (0.0, 0.5, 0.0),
             (0.5, 0.0, 0.0), (np.nan, 0.5, 0.0)]
        )
        with self.assertRaisesRegex(ValueError, "must all be finite"):
            MODULE.infer_uniform_grid(kpoints, 2, 2)

    def test_boundary_g_matching(self):
        left = np.asarray([(0, 0, 0), (1, 0, 0), (-1, 0, 0)])
        right = np.asarray([(0, 0, 0), (1, 0, 0), (2, 0, 0)])
        il, ir = MODULE.common_g_indices(left, right, np.asarray((1, 0, 0)))
        pairs = [(tuple(left[i]), tuple(right[j] - (1, 0, 0))) for i, j in zip(il, ir)]
        self.assertEqual(pairs, [((0, 0, 0), (0, 0, 0)), ((1, 0, 0), (1, 0, 0)), ((-1, 0, 0), (-1, 0, 0))])

    def test_cumulative_t0_vertex_selection_formula(self):
        phi_v = np.asarray([[0.4]])
        phi_v1 = np.asarray([[1.2]])
        phi_v2 = np.asarray([[2.0]])
        energy1 = np.asarray([[[0.0, 1.0, 2.0, 3.0]]])
        energy2 = np.asarray([[[0.5, 1.5, 2.5, 3.5]]])
        effective, n0, n1, n2 = MODULE.cumulative_t0_effective_phi(
            phi_v, phi_v1, phi_v2, energy1, energy2, 1.25
        )
        # Vertices select V+2, V+1, V, V respectively.
        self.assertAlmostEqual(float(effective[0, 0]), (2.0 + 1.2 + 0.4 + 0.4) / 4)
        self.assertEqual((int(n0[0, 0]), int(n1[0, 0]), int(n2[0, 0])), (2, 1, 1))

    def test_cumulative_state_count_invariant_rejects_inverted_energies(self):
        phi = np.zeros((1, 1))
        energy1 = np.asarray([[[1.0, 1.0, 2.0, 3.0]]])
        energy2 = np.asarray([[[0.5, 1.5, 2.5, 3.5]]])
        with self.assertRaisesRegex(RuntimeError, "state-count invariant"):
            MODULE.cumulative_t0_effective_phi(
                phi, phi, phi, energy1, energy2, 0.75
            )
        energy1[0, 0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "must be finite"):
            MODULE.cumulative_t0_effective_phi(
                phi, phi, phi, energy1, energy2, 0.75
            )

    def test_default_maps_follow_requested_energy_band(self):
        maps = dict(MODULE.default_map_specifications(5))
        self.assertEqual(maps["V"], [1, 2, 3, 4])
        self.assertEqual(maps["VC"], [1, 2, 3, 4, 5])
        self.assertEqual(maps["VCC"], [1, 2, 3, 4, 5, 6])
        self.assertEqual(maps["C"], [5])
        self.assertEqual(maps["Cpair"], [5, 6])
        with self.assertRaisesRegex(ValueError, "at least 2"):
            MODULE.default_map_specifications(1)

    def test_explicit_maps_cannot_bypass_target_or_cumulative_bounds(self):
        with self.assertRaisesRegex(ValueError, "at least 2"):
            MODULE.validate_map_specifications([("X", [1])], 1, 4, False)
        # Raw export may target NBANDS-1 because it only reads target+1.
        MODULE.validate_map_specifications([("X", [1])], 3, 4, False)
        cumulative = [
            ("V", [1, 2]),
            ("VC", [1, 2, 3]),
            ("VCC", [1, 2, 3, 4]),
        ]
        MODULE.validate_map_specifications(cumulative, 3, 4, True)
        with self.assertRaisesRegex(ValueError, r"transport-t0 requires target\+2"):
            MODULE.validate_map_specifications(
                cumulative, 3, 4, True, require_transport_headroom=True
            )
        with self.assertRaisesRegex(ValueError, r"leave target\+1"):
            MODULE.validate_map_specifications([("X", [1])], 4, 4, False)
        with self.assertRaisesRegex(ValueError, "requires cumulative maps"):
            MODULE.validate_map_specifications([("X", [2])], 2, 4, True)
        with self.assertRaisesRegex(ValueError, "exceeds WAVECAR"):
            MODULE.validate_map_specifications([("X", [1, 2, 3, 4, 5])], 2, 4, False)
        for unsafe in ("../escape=1", "foo/bar=1"):
            with self.subTest(unsafe=unsafe), self.assertRaises(argparse.ArgumentTypeError):
                MODULE.parse_map(unsafe)

    def test_expected_cli_error_has_no_traceback(self):
        completed = subprocess.run(
            [
                sys.executable,
                str(ROOT / "tools" / "wavecar_fukui.py"),
                "missing-WAVECAR",
                "--nx", "1", "--ny", "1",
                "--energy-band", "1",
                "--output-dir", "unused",
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("error:", completed.stderr)
        self.assertNotIn("Traceback", completed.stderr)

    def test_valley_radius_analytic_overlap_is_rejected(self):
        fake = SimpleNamespace(
            header=SimpleNamespace(reciprocal=np.eye(3)),
            kpoints=np.zeros((1, 3)),
        )
        grid = MODULE.Grid(1, 1, np.asarray([[0]]), np.zeros(3))
        with self.assertRaisesRegex(ValueError, "overlap analytically"):
            MODULE.make_valley_partition(
                fake,
                grid,
                np.asarray((0.0, 0.0, 0.0)),
                np.asarray((0.4, 0.0, 0.0)),
                0.2,
            )

    def test_nonpositive_valley_radius_is_rejected(self):
        fake = SimpleNamespace(
            header=SimpleNamespace(reciprocal=np.eye(3)),
            kpoints=np.zeros((1, 3)),
        )
        grid = MODULE.Grid(1, 1, np.asarray([[0]]), np.zeros(3))
        with self.assertRaisesRegex(ValueError, "must be positive"):
            MODULE.make_valley_partition(
                fake,
                grid,
                np.asarray((0.0, 0.0, 0.0)),
                np.asarray((0.4, 0.0, 0.0)),
                0.0,
            )

    def test_valley_inputs_must_be_finite_and_periodically_distinct(self):
        fake = SimpleNamespace(
            header=SimpleNamespace(reciprocal=np.eye(3)),
            kpoints=np.zeros((1, 3)),
        )
        grid = MODULE.Grid(1, 1, np.asarray([[0]]), np.zeros(3))
        with self.assertRaisesRegex(ValueError, "finite three-component"):
            MODULE.make_valley_partition(
                fake, grid, np.asarray((np.nan, 0.0, 0.0)),
                np.asarray((0.4, 0.0, 0.0)), None,
            )
        with self.assertRaisesRegex(ValueError, "periodically distinct"):
            MODULE.make_valley_partition(
                fake, grid, np.zeros(3), np.asarray((1.0, 0.0, 0.0)), None,
            )
        with self.assertRaisesRegex(ValueError, "periodically distinct"):
            MODULE.make_valley_partition(
                fake, grid, np.zeros(3), np.asarray((2.0, 0.0, 0.0)), None,
            )
        with self.assertRaisesRegex(ValueError, "must be finite"):
            MODULE.make_valley_partition(
                fake, grid, np.zeros(3), np.asarray((0.4, 0.0, 0.0)), np.inf,
            )
        with self.assertRaises(argparse.ArgumentTypeError):
            MODULE.parse_fractional_vector("nan,0,0")

    def test_voronoi_exact_tie_is_outside(self):
        fake = SimpleNamespace(
            header=SimpleNamespace(reciprocal=np.eye(3)),
            kpoints=np.asarray(
                [(0.0, 0.0, 0.0), (0.0, 0.5, 0.0),
                 (0.5, 0.0, 0.0), (0.5, 0.5, 0.0)]
            ),
        )
        grid = MODULE.Grid(2, 2, np.asarray([[0, 1], [2, 3]]), np.zeros(3))
        mask_k, mask_kp, outside, _, _ = MODULE.make_valley_partition(
            fake,
            grid,
            np.asarray((0.0, 0.25, 0.0)),
            np.asarray((0.5, 0.25, 0.0)),
            None,
        )
        self.assertTrue(outside[0, 0])
        self.assertFalse(mask_k[0, 0])
        self.assertFalse(mask_kp[0, 0])
        self.assertTrue(np.all(mask_k | mask_kp | outside))
        self.assertTrue(
            np.all(
                mask_k.astype(int) + mask_kp.astype(int) + outside.astype(int)
                == 1
            )
        )

    def test_skew_reciprocal_basis_finds_distant_integer_image(self):
        reciprocal = np.asarray(
            [(1.0, 0.0, 0.0), (10.0, 0.1, 0.0), (0.0, 0.0, 1.0)]
        )
        distance = MODULE.periodic_reciprocal_distance(
            np.asarray((0.0, 0.2, 0.0)), np.zeros(3), reciprocal
        )
        # The closest image is (-2, 0), outside the legacy +/-1 stencil.
        self.assertAlmostEqual(distance, 0.02, places=12)
        np.testing.assert_allclose(
            MODULE.closest_periodic_image(
                np.zeros(3), np.asarray((0.0, 0.2, 0.0)), reciprocal
            ),
            (-2.0, 0.2, 0.0),
            atol=1.0e-12,
        )

        fake = SimpleNamespace(
            header=SimpleNamespace(reciprocal=reciprocal),
            kpoints=np.asarray(
                [(-0.25, -0.05, 0.0), (-0.25, 0.45, 0.0),
                 (0.25, -0.05, 0.0), (0.25, 0.45, 0.0)]
            ),
        )
        grid = MODULE.Grid(2, 2, np.asarray([[0, 1], [2, 3]]), np.zeros(3))
        mask_k, mask_kp, _, _, _ = MODULE.make_valley_partition(
            fake,
            grid,
            np.zeros(3),
            np.asarray((0.1, 0.2, 0.0)),
            None,
        )
        self.assertTrue(mask_k[0, 0])
        self.assertFalse(mask_kp[0, 0])

        high_shear = np.asarray(
            [(1.0, 0.0, 0.0), (1.0e8, 1.0, 0.0), (0.0, 0.0, 1.0)]
        )
        self.assertAlmostEqual(
            MODULE.periodic_reciprocal_distance(
                np.asarray((0.0, 0.25, 0.0)), np.zeros(3), high_shear
            ),
            0.25,
            places=12,
        )
        tiny = np.diag((1.0e-20, 2.0e-20, 1.0))
        self.assertTrue(
            np.isclose(
                MODULE.periodic_reciprocal_distance(
                    np.asarray((0.2, 0.3, 0.0)), np.zeros(3), tiny
                ),
                math.sqrt(40.0) * 1.0e-21,
                rtol=1.0e-12,
                atol=0.0,
            )
        )

    def test_radius_partition_requires_cells_for_both_valleys(self):
        fake = SimpleNamespace(
            header=SimpleNamespace(reciprocal=np.eye(3)),
            kpoints=np.asarray(
                [(0.0, 0.0, 0.0), (0.0, 0.5, 0.0),
                 (0.5, 0.0, 0.0), (0.5, 0.5, 0.0)]
            ),
        )
        grid = MODULE.Grid(2, 2, np.asarray([[0, 1], [2, 3]]), np.zeros(3))
        with self.assertRaisesRegex(ValueError, "selects no plaquette centers"):
            MODULE.make_valley_partition(
                fake, grid, np.asarray((0.0, 0.0, 0.0)),
                np.asarray((0.4, 0.0, 0.0)), 1.0e-6,
            )

    def test_transport_energy_window_rejects_unrepresented_band(self):
        fake = SimpleNamespace(
            header=SimpleNamespace(nbands=4),
            energies=np.asarray([[-1.0, 0.0, 0.5, 0.6]]),
        )
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "reaches band 4"):
                MODULE.write_t0_transport(
                    Path(directory), fake, None, {}, 2, 0.0, 0.7, 2,
                    np.zeros(3), np.ones(3) * 0.5, None,
                    1e-3, 1e-5, 0.8 * np.pi, False,
                )

    def test_failed_transport_removes_stale_csv_diagnostics_and_plot(self):
        fake = SimpleNamespace(
            header=SimpleNamespace(nbands=4),
            energies=np.asarray([[-1.0, 0.0, 0.5, 0.6]]),
        )
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            for name in MODULE.TRANSPORT_OUTPUT_NAMES:
                (output / name).write_text("stale", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "reaches band 4"):
                MODULE.write_t0_transport(
                    output, fake, None, {}, 2, 0.0, 0.7, 2,
                    np.zeros(3), np.ones(3) * 0.5, None,
                    1e-3, 1e-5, 0.8 * np.pi, False,
                )
            for name in MODULE.TRANSPORT_OUTPUT_NAMES:
                self.assertFalse((output / name).exists(), name)

    def test_wavecar_input_cannot_collide_with_any_planned_output(self):
        maps = [("V", [1]), ("VC", [1, 2])]
        names = MODULE.planned_output_names(maps, plot=True, transport_t0=True)
        self.assertTrue(
            {"diagnostics.json", "fukui_V.csv", "fukui_VC.csv"}.issubset(names)
        )
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            for name in names:
                collision = output / name
                collision.write_bytes(b"input must survive")
                with self.subTest(name=name), self.assertRaisesRegex(
                    ValueError, "resolves to planned output"
                ):
                    MODULE.validate_input_output_collision(collision, output, names)
                self.assertEqual(collision.read_bytes(), b"input must survive")

            collision = output / "transport_t0.csv"
            fake = SimpleNamespace(path=collision)
            with self.assertRaisesRegex(ValueError, "refusing to overwrite"):
                MODULE.write_t0_transport(
                    output, fake, None, {}, 2, 0.0, 0.1, 2,
                    np.zeros(3), np.ones(3), None,
                    1e-3, 1e-5, 0.8 * np.pi, False,
                )
            self.assertEqual(collision.read_bytes(), b"input must survive")

    def test_transport_refuses_low_quality_active_cell(self):
        fake = SimpleNamespace(
            header=SimpleNamespace(nbands=4, reciprocal=np.eye(3)),
            energies=np.asarray([[-1.0, 0.0, 0.5, 10.0]]),
            kpoints=np.zeros((1, 3)),
        )
        grid = MODULE.Grid(1, 1, np.asarray([[0]]), np.zeros(3))

        def links(minimum):
            one = np.ones((1, 1))
            return MODULE.LinkSet(
                phase=np.zeros((1, 1)),
                logabs=np.zeros((1, 1)),
                min_singular=one * minimum,
                max_singular=one,
                coverage_left=one,
                coverage_right=one,
            )

        stable, unstable = links(1.0), links(1.0e-6)
        maps = {
            "V": ([1], np.asarray([[0.0]]), stable, stable),
            "V1": ([1, 2], np.asarray([[0.2]]), unstable, unstable),
            "V2": ([1, 2, 3], np.asarray([[0.4]]), stable, stable),
        }
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            with self.assertRaisesRegex(RuntimeError, "transport validation failed"):
                MODULE.write_t0_transport(
                    output, fake, grid, maps, 2, 0.1, 0.2, 2,
                    np.zeros(3), np.ones(3) * 0.5, None,
                    1e-3, 1e-5, 0.8 * np.pi, False,
                )
            self.assertTrue((output / "transport_t0_diagnostics.json").exists())
            self.assertFalse((output / "transport_t0.csv").exists())

    def test_v1_flux_checks_gap_at_all_four_vertices(self):
        fake = SimpleNamespace(
            header=SimpleNamespace(nbands=4, reciprocal=np.eye(3)),
            energies=np.asarray(
                [
                    [-1.0, 0.4, 0.8, 10.0],
                    [-1.0, 1.0, 1.000001, 10.0],
                    [-1.0, 1.0, 2.0, 10.0],
                    [-1.0, 1.0, 2.0, 10.0],
                ]
            ),
            kpoints=np.asarray(
                [(0.0, 0.0, 0.0), (0.0, 0.5, 0.0), (0.5, 0.0, 0.0), (0.5, 0.5, 0.0)]
            ),
        )
        grid = MODULE.Grid(2, 2, np.asarray([[0, 1], [2, 3]]), np.zeros(3))
        one = np.ones((2, 2))
        stable = MODULE.LinkSet(
            phase=np.zeros((2, 2)),
            logabs=np.zeros((2, 2)),
            min_singular=one,
            max_singular=one,
            coverage_left=one,
            coverage_right=one,
        )
        maps = {
            "V": ([1], np.zeros((2, 2)), stable, stable),
            "V1": ([1, 2], np.zeros((2, 2)), stable, stable),
            "V2": ([1, 2, 3], np.zeros((2, 2)), stable, stable),
        }
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            with self.assertRaisesRegex(RuntimeError, "active plaquettes"):
                MODULE.write_t0_transport(
                    output, fake, grid, maps, 2, 0.5, 0.51, 2,
                    np.zeros(3), np.asarray((0.25, 0.25, 0.0)), None,
                    1e-3, 1e-5, 0.8 * np.pi, False,
                )
            diagnostics = json.loads(
                (output / "transport_t0_diagnostics.json").read_text()
            )
            self.assertAlmostEqual(
                diagnostics["invalid_mu"][0]["worst_neighbor_gap_ev"], 1e-6
            )
            self.assertFalse((output / "transport_t0.csv").exists())

    def test_valence_baseline_link_gap_and_phase_are_hard_gates(self):
        grid = MODULE.Grid(1, 1, np.asarray([[0]]), np.zeros(3))

        def links(minimum):
            one = np.ones((1, 1))
            return MODULE.LinkSet(
                phase=np.zeros((1, 1)),
                logabs=np.zeros((1, 1)),
                min_singular=one * minimum,
                max_singular=one,
                coverage_left=one,
                coverage_right=one,
            )

        stable = links(1.0)
        cases = (
            ("link", np.asarray([[-1.0, 1.0, 2.0, 10.0]]), links(1e-6), 0.0),
            ("gap", np.asarray([[0.0, 1e-6, 2.0, 10.0]]), stable, 0.0),
            ("phase", np.asarray([[-1.0, 1.0, 2.0, 10.0]]), stable, 0.9 * np.pi),
        )
        for name, energies, baseline_links, baseline_phi in cases:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as directory:
                fake = SimpleNamespace(
                    header=SimpleNamespace(nbands=4, reciprocal=np.eye(3)),
                    energies=energies,
                    kpoints=np.zeros((1, 3)),
                )
                maps = {
                    "V": ([1], np.asarray([[baseline_phi]]), baseline_links, baseline_links),
                    "V1": ([1, 2], np.asarray([[0.2]]), stable, stable),
                    "V2": ([1, 2, 3], np.asarray([[0.4]]), stable, stable),
                }
                output = Path(directory)
                with self.assertRaisesRegex(RuntimeError, "valence baseline failed"):
                    MODULE.write_t0_transport(
                        output, fake, grid, maps, 2, 0.1, 0.2, 2,
                        np.zeros(3), np.ones(3) * 0.5, None,
                        1e-3, 1e-5, 0.8 * np.pi, True,
                    )
                diagnostics = json.loads(
                    (output / "transport_t0_diagnostics.json").read_text()
                )
                self.assertFalse(diagnostics["baseline_quality"]["pass"])
                self.assertFalse((output / "transport_t0.csv").exists())

    def test_transport_lower_window_must_stay_above_valence_maximum(self):
        fake = SimpleNamespace(
            header=SimpleNamespace(nbands=4),
            energies=np.asarray([[0.2, 1.0, 2.0, 10.0]]),
        )
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "fixed valence baseline would contain holes"):
                MODULE.write_t0_transport(
                    Path(directory), fake, None, {}, 2, 0.1, 0.3, 2,
                    np.zeros(3), np.ones(3) * 0.5, None,
                    1e-3, 1e-5, 0.8 * np.pi, False,
                )

    def test_periodic_bilinear_sampling_respects_shifted_center_grid(self):
        grid = MODULE.Grid(2, 2, np.asarray([[0, 1], [2, 3]]), np.asarray((0.125, 0.125, 0.0)))
        values = np.asarray([[1.0, 2.0], [3.0, 4.0]])
        sampled = MODULE.periodic_bilinear_sample(
            values, np.asarray([[0.875, 0.375, 0.0]]), grid, centered=True
        )
        self.assertAlmostEqual(float(sampled[0]), 3.0)

    def test_shifted_plot_domain_uses_periodic_marker_image(self):
        offset = np.asarray((0.125, 0.25, 0.0))
        self.assertEqual(
            MODULE.fractional_plot_extent(offset), (0.125, 1.125, 0.25, 1.25)
        )
        marker = MODULE.periodic_image_in_plot_domain(
            np.asarray((0.0, 0.1, 0.0)), offset
        )
        np.testing.assert_allclose(marker, (1.0, 1.1, 0.0))
        self.assertTrue(offset[0] <= marker[0] < offset[0] + 1.0)
        self.assertTrue(offset[1] <= marker[1] < offset[1] + 1.0)
        self.assertEqual(
            MODULE.symmetric_finite_color_limits(
                np.asarray([[-2.0, 5.0], [np.nan, 1.0]])
            ),
            (-5.0, 5.0),
        )
        self.assertEqual(
            MODULE.symmetric_finite_color_limits(np.zeros((2, 2))), (-1.0, 1.0)
        )

    def test_invalid_transport_plot_is_written(self):
        header = (
            "mu_eV,quality,delta_sigma_e2_over_h,delta_sigma_K_e2_over_h,"
            "delta_sigma_Kp_e2_over_h,delta_sigma_valley_contrast_e2_over_h\n"
        )
        body = "0.0,PASS,0,0,0,0\n0.1,INVALID,0.1,0.2,-0.1,0.3\n"
        with tempfile.TemporaryDirectory() as directory:
            csv_path = Path(directory) / "transport.csv"
            png_path = Path(directory) / "transport.png"
            csv_path.write_text(header + body, encoding="utf-8")
            MODULE.plot_transport_csv(csv_path, png_path)
            self.assertTrue(png_path.exists())
            self.assertGreater(png_path.stat().st_size, 1000)

    def test_checked_in_line_wavecar_byte_recl_and_g_count(self):
        path = ROOT / "1H-MoS2" / "KPATH" / "2.band" / "WAVECAR"
        if not path.exists():
            self.skipTest("checked-in WAVECAR fixture is absent")
        parsed = MODULE.Wavecar(path, spinor_components=2)
        self.assertEqual(parsed.header.logical_recl, 38160)
        self.assertEqual(parsed.header.stride_bytes, 38160)
        self.assertEqual(parsed.header.nkpoints, 48)
        self.assertEqual(parsed.header.nbands, 32)
        self.assertEqual(parsed.nbmax, (6, 6, 25))
        self.assertEqual(2 * len(parsed.g_vectors(0)), int(parsed.nplane[0]))
        with self.assertRaisesRegex(ValueError, "off the requested uniform grid"):
            MODULE.infer_uniform_grid(parsed.kpoints, 8, 6)


if __name__ == "__main__":
    unittest.main()
