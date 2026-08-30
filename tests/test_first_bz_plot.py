import importlib.util
import inspect
import math
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "wavecar_fukui_first_bz", ROOT / "tools" / "wavecar_fukui.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class FirstBrillouinZonePlotTests(unittest.TestCase):
    @staticmethod
    def hexagonal_reciprocal() -> np.ndarray:
        return np.asarray(
            (
                (1.0, 0.0, 0.0),
                (-0.5, math.sqrt(3.0) / 2.0, 0.0),
                (0.0, 0.0, 1.0),
            )
        )

    def test_hexagonal_wigner_seitz_cell_has_six_vertices_and_correct_area(self):
        reciprocal = self.hexagonal_reciprocal()
        polygon = MODULE.wigner_seitz_polygon_2d(reciprocal)
        self.assertEqual(polygon.shape, (6, 2))
        self.assertAlmostEqual(
            MODULE.polygon_signed_area(polygon),
            float(np.linalg.norm(np.cross(reciprocal[0], reciprocal[1]))),
            places=12,
        )
        radii = np.linalg.norm(polygon, axis=1)
        np.testing.assert_allclose(radii, np.full(6, radii[0]), rtol=1e-12, atol=1e-12)

    def test_wigner_seitz_cell_is_basis_independent_for_large_integer_shear(self):
        reciprocal = np.asarray(
            ((3.0, 0.0, 0.0), (8.0, 1.0, 0.0), (0.0, 0.0, 2.0))
        )
        polygon = MODULE.wigner_seitz_polygon_2d(reciprocal)
        self.assertIn(polygon.shape[0], (4, 6))
        self.assertAlmostEqual(MODULE.polygon_signed_area(polygon), 3.0, places=11)

    def test_tilted_reciprocal_plane_uses_an_intrinsic_cartesian_frame(self):
        reciprocal = self.hexagonal_reciprocal()
        angle = 0.37
        rotation = np.asarray(
            (
                (math.cos(angle), 0.0, math.sin(angle)),
                (0.0, 1.0, 0.0),
                (-math.sin(angle), 0.0, math.cos(angle)),
            )
        )
        tilted = reciprocal @ rotation.T
        polygon = MODULE.wigner_seitz_polygon_2d(tilted)
        self.assertEqual(polygon.shape, (6, 2))
        self.assertAlmostEqual(
            MODULE.polygon_signed_area(polygon),
            float(np.linalg.norm(np.cross(tilted[0], tilted[1]))),
            places=12,
        )

    def test_k_and_kprime_markers_fold_onto_first_bz_boundary(self):
        reciprocal = self.hexagonal_reciprocal()
        polygon = MODULE.wigner_seitz_polygon_2d(reciprocal)
        folded_k = MODULE.fold_fractional_to_first_bz(
            np.asarray((2.0 / 3.0, 1.0 / 3.0, 0.0)), reciprocal
        )
        folded_kp = MODULE.fold_fractional_to_first_bz(
            np.asarray((1.0 / 3.0, 2.0 / 3.0, 0.0)), reciprocal
        )
        self.assertTrue(bool(MODULE.point_in_convex_polygon_2d(folded_k, polygon)))
        self.assertTrue(bool(MODULE.point_in_convex_polygon_2d(folded_kp, polygon)))
        self.assertGreater(float(np.linalg.norm(folded_k - folded_kp)), 0.1)

        k_images = MODULE.periodic_images_in_first_bz(
            np.asarray((2.0 / 3.0, 1.0 / 3.0, 0.0)), reciprocal, polygon
        )
        kp_images = MODULE.periodic_images_in_first_bz(
            np.asarray((1.0 / 3.0, 2.0 / 3.0, 0.0)), reciprocal, polygon
        )
        self.assertEqual(k_images.shape, (3, 2))
        self.assertEqual(kp_images.shape, (3, 2))
        labeled = [(point, "K") for point in k_images] + [
            (point, "Kp") for point in kp_images
        ]
        labeled.sort(key=lambda item: math.atan2(item[0][1], item[0][0]))
        labels = [label for _point, label in labeled]
        self.assertTrue(
            all(labels[index] != labels[(index + 1) % 6] for index in range(6))
        )
        for point, _label in labeled:
            self.assertLess(
                float(np.min(np.linalg.norm(polygon - point, axis=1))), 1.0e-11
            )

    def test_boundary_marker_enumeration_is_reciprocal_scale_invariant(self):
        for scale in (1.0e-8, 1.0e8):
            with self.subTest(scale=scale):
                reciprocal = self.hexagonal_reciprocal()
                reciprocal[:2] *= scale
                polygon = MODULE.wigner_seitz_polygon_2d(reciprocal)
                expected_area = float(
                    np.linalg.norm(np.cross(reciprocal[0], reciprocal[1]))
                )
                self.assertTrue(
                    math.isclose(
                        MODULE.polygon_signed_area(polygon),
                        expected_area,
                        rel_tol=2.0e-12,
                    )
                )
                k_images = MODULE.periodic_images_in_first_bz(
                    np.asarray((2.0 / 3.0, 1.0 / 3.0, 0.0)),
                    reciprocal,
                    polygon,
                )
                kp_images = MODULE.periodic_images_in_first_bz(
                    np.asarray((1.0 / 3.0, 2.0 / 3.0, 0.0)),
                    reciprocal,
                    polygon,
                )
                self.assertEqual(k_images.shape, (3, 2))
                self.assertEqual(kp_images.shape, (3, 2))

    def test_first_bz_raster_preserves_only_source_plaquette_values(self):
        reciprocal = self.hexagonal_reciprocal()
        polygon = MODULE.wigner_seitz_polygon_2d(reciprocal)
        grid = MODULE.Grid(
            nx=2,
            ny=2,
            index=np.asarray(((0, 1), (2, 3))),
            offset=np.zeros(3),
        )
        source = np.asarray(((1.0, 2.0), (3.0, 4.0)))
        raster, extent = MODULE.first_bz_display_grid(
            source, grid, reciprocal, polygon, maximum_pixels=80
        )
        finite = raster[np.isfinite(raster)]
        self.assertGreater(finite.size, 100)
        self.assertTrue(set(np.unique(finite)).issubset({1.0, 2.0, 3.0, 4.0}))
        self.assertLess(extent[0], 0.0)
        self.assertGreater(extent[1], 0.0)
        self.assertLess(extent[2], 0.0)
        self.assertGreater(extent[3], 0.0)

        x_limits, y_limits = MODULE.padded_cartesian_plot_limits(extent)
        self.assertLess(x_limits[0], extent[0])
        self.assertGreater(x_limits[1], extent[1])
        self.assertLess(y_limits[0], extent[2])
        self.assertGreater(y_limits[1], extent[3])
        self.assertAlmostEqual(
            x_limits[1] - x_limits[0], 1.1 * (extent[1] - extent[0])
        )
        self.assertAlmostEqual(
            y_limits[1] - y_limits[0], 1.1 * (extent[3] - extent[2])
        )

    def test_first_bz_plot_is_opt_in_and_writes_a_map(self):
        signature = inspect.signature(MODULE.plot_k_resolved_maps)
        self.assertEqual(signature.parameters["plot_domain"].default, "fractional")
        self.assertIn("Shortest periodic-image", MODULE.SHORTEST_K_KP_LINE_TITLE)
        self.assertIn("Gamma", MODULE.SHORTEST_K_KP_LINE_TITLE)
        reciprocal = self.hexagonal_reciprocal()
        grid = MODULE.Grid(
            nx=2,
            ny=2,
            index=np.asarray(((0, 1), (2, 3))),
            offset=np.zeros(3),
        )
        shape = (2, 2)
        links = MODULE.LinkSet(
            phase=np.zeros(shape),
            logabs=np.zeros(shape),
            min_singular=np.ones(shape),
            max_singular=np.ones(shape),
            coverage_left=np.ones(shape),
            coverage_right=np.ones(shape),
        )
        maps = {
            "V": ([1], np.zeros(shape), links, links),
            "VC": ([1, 2], np.zeros(shape), links, links),
            "VCC": ([1, 2, 3], np.zeros(shape), links, links),
        }
        fake = SimpleNamespace(
            header=SimpleNamespace(reciprocal=reciprocal),
            energies=np.asarray(
                (
                    (-1.0, 0.1, 1.0),
                    (-1.0, 0.2, 1.1),
                    (-1.0, 0.3, 1.2),
                    (-1.0, 0.4, 1.3),
                )
            ),
        )
        with tempfile.TemporaryDirectory() as directory:
            with mock.patch.object(
                MODULE,
                "first_bz_display_grid",
                side_effect=AssertionError("legacy default must not use first-BZ raster"),
            ) as first_bz_raster:
                default_output = MODULE.plot_k_resolved_maps(
                    Path(directory),
                    fake,
                    grid,
                    maps,
                    2,
                    np.asarray((2.0 / 3.0, 1.0 / 3.0, 0.0)),
                    np.asarray((1.0 / 3.0, 2.0 / 3.0, 0.0)),
                )
            first_bz_raster.assert_not_called()
            self.assertEqual(default_output.name, "wavecar_fukui_kresolved.png")
            self.assertTrue(default_output.exists())
            output = MODULE.plot_k_resolved_maps(
                Path(directory),
                fake,
                grid,
                maps,
                2,
                np.asarray((2.0 / 3.0, 1.0 / 3.0, 0.0)),
                np.asarray((1.0 / 3.0, 2.0 / 3.0, 0.0)),
                plot_domain="first-bz",
            )
            self.assertTrue(output.exists())
            self.assertGreater(output.stat().st_size, 1000)
            line_output = MODULE.plot_k_to_kp_line(
                Path(directory),
                fake,
                grid,
                maps,
                2,
                np.asarray((2.0 / 3.0, 1.0 / 3.0, 0.0)),
                np.asarray((1.0 / 3.0, 2.0 / 3.0, 0.0)),
            )
            self.assertEqual(line_output.name, "wavecar_fukui_line_K_Kp.png")
            self.assertTrue(line_output.exists())
            self.assertGreater(line_output.stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
