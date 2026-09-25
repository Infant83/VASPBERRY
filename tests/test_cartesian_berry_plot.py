"""Physical first-BZ geometry and unsmoothed native plaquette rendering."""
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from plot_berry_curvature import (first_bz_polygon, polygon_area, plaquettes_in_first_bz,
                                 read_native_curvature, reciprocal_from_poscar, uniform_plaquettes)


def reciprocal(basis):
    result = np.eye(3)
    result[:2, :2] = basis
    return result


def grid(nx=6, ny=8):
    x, y = np.meshgrid((np.arange(nx) + .5) / nx, (np.arange(ny) + .5) / ny, indexing="ij")
    return np.column_stack((x.ravel(), y.ravel(), np.zeros(nx * ny)))


class CartesianFirstBZTests(unittest.TestCase):
    def test_square_brillouin_zone_has_physical_cartesian_edges(self):
        vertices = first_bz_polygon(reciprocal([[2, 0], [0, 4]]))
        self.assertEqual(len(vertices), 4)
        self.assertTrue(np.allclose(np.max(abs(vertices), axis=0), [1, 2]))
        self.assertAlmostEqual(polygon_area(vertices), 8)

    def test_triangular_reciprocal_lattice_has_regular_hexagonal_bz(self):
        basis = np.array([[1., 0], [-.5, np.sqrt(3) / 2]])
        vertices = first_bz_polygon(reciprocal(basis))
        self.assertEqual(len(vertices), 6)
        self.assertTrue(np.allclose(np.linalg.norm(vertices, axis=1), 1 / np.sqrt(3)))
        self.assertAlmostEqual(polygon_area(vertices), np.sqrt(3) / 2)

    def test_unimodular_skew_basis_does_not_change_physical_first_bz(self):
        basis = np.array([[1., 0], [-.5, np.sqrt(3) / 2]])
        original = first_bz_polygon(reciprocal(basis))
        skewed = first_bz_polygon(reciprocal(np.array([[1, 17], [0, 1]]) @ basis))
        distances = np.linalg.norm(original[:, None, :] - skewed[None, :, :], axis=2)
        self.assertEqual(len(original), len(skewed))
        self.assertLess(np.max(np.min(distances, axis=1)), 1e-10)

    def test_rigid_rotation_rotates_physical_bz(self):
        angle = .39
        rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        basis = np.array([[2., 0], [0, 3.]])
        expected = first_bz_polygon(reciprocal(basis)) @ rotation
        observed = first_bz_polygon(reciprocal(basis @ rotation))
        distances = np.linalg.norm(expected[:, None, :] - observed[None, :, :], axis=2)
        self.assertLess(np.max(np.min(distances, axis=1)), 1e-10)

    def test_clipping_preserves_area_and_area_weighted_native_values(self):
        q = grid()
        values = np.sin(2 * np.pi * q[:, 0]) + .2 * np.cos(2 * np.pi * q[:, 1])
        basis = reciprocal([[1., 0], [-.5, np.sqrt(3) / 2]])
        polygons, colors, boundary, mesh = plaquettes_in_first_bz(q, values, basis)
        areas = np.array([polygon_area(polygon) for polygon in polygons])
        area = polygon_area(boundary)
        self.assertEqual(mesh, (6, 8))
        self.assertAlmostEqual(areas.sum(), area, places=10)
        self.assertAlmostEqual(np.dot(areas, colors), area * values.mean(), places=10)
        self.assertTrue(set(colors).issubset(set(values)))  # No interpolated values.

    def test_reversed_plane_basis_still_tiles_once(self):
        q = grid(4, 4)
        polygons, colors, boundary, _ = plaquettes_in_first_bz(q, np.ones(len(q)),
                                                              reciprocal([[0, 1], [1, 0]]))
        self.assertAlmostEqual(sum(polygon_area(polygon) for polygon in polygons), 1)
        self.assertTrue(np.all(colors == 1))

    def test_inconsistent_periodic_display_values_are_rejected(self):
        q = grid(4, 4)
        tiled = np.vstack((q, q + [1, 0, 0]))
        values = np.ones(len(tiled))
        values[-1] = 2
        with self.assertRaisesRegex(ValueError, "inconsistent values"):
            uniform_plaquettes(tiled, values)

    def test_missing_plaquette_is_rejected(self):
        q = grid(4, 4)[:-1]
        with self.assertRaisesRegex(ValueError, "complete uniform"):
            uniform_plaquettes(q, np.ones(len(q)))

    def test_tilted_reciprocal_plane_cannot_be_mislabelled_cartesian_xy(self):
        basis = reciprocal([[1, 0], [0, 1]])
        basis[0, 2] = .2
        with self.assertRaisesRegex(ValueError, "Cartesian xy plane"):
            first_bz_polygon(basis)

    def test_actual_mos2_native_data_match_lattice_and_cover_hexagon(self):
        lattice = reciprocal_from_poscar(ROOT / "examples/1H-MoS2/POSCAR")
        source = ROOT / "examples/1H-MoS2/BERRYCURV.dat"
        q, values = read_native_curvature(source, lattice)
        polygons, colors, boundary, mesh = plaquettes_in_first_bz(q, values, lattice)
        self.assertEqual(mesh, (12, 12))
        self.assertEqual(len(boundary), 6)
        self.assertAlmostEqual(sum(polygon_area(polygon) for polygon in polygons),
                               abs(np.linalg.det(lattice[:2, :2])), places=8)
        with self.assertRaisesRegex(ValueError, "disagree"):
            read_native_curvature(source, 1.01 * lattice)

    def test_native_oriented_curvature_is_converted_to_cartesian_z(self):
        basis = reciprocal([[0, 1], [1, 0]])
        q = grid(4, 4)
        raw = np.column_stack((q @ basis, np.ones(len(q)), q))
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "BERRYCURV.dat"
            np.savetxt(path, raw)
            _, curvature = read_native_curvature(path, basis)
            self.assertTrue(np.all(curvature == -1))
            path.write_text("# KUBO point data\n" + path.read_text())
            with self.assertRaisesRegex(ValueError, "not Kubo point data"):
                read_native_curvature(path, basis)


if __name__ == "__main__":
    unittest.main()
