"""Regression checks for corrupted coordinates in native Bi Fukui output.

Run with ``python3 -m unittest discover -s tests -p test_real_fukui_example.py``.
These checks use archived public output/geometry and do not need the LFS payload.
"""
from pathlib import Path
import importlib.util
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples/features/fukui-chern"
sys.path.insert(0, str(EXAMPLE))
SPEC = importlib.util.spec_from_file_location("real_bi_fukui_example", EXAMPLE / "run.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
validate_field = MODULE.validate_field


class NativePlaquetteValidation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        here = EXAMPLE
        cls.raw = np.loadtxt(here / "reference/BERRYCURV.dat")
        cls.reference = np.loadtxt(here / "reference/plaquettes.csv", delimiter=",", skiprows=1)
        lines = (ROOT / "examples/Bi_Z2/archive-2016-run/EIGENVAL").read_text().splitlines()
        cls.kpoints = np.array([[float(x) for x in lines[7 + ik * 20].split()[:3]] for ik in range(144)])
        poscar = (ROOT / "examples/Bi_Z2/inputs/POSCAR").read_text().splitlines()
        lattice = float(poscar[1]) * np.array([[float(x) for x in line.split()] for line in poscar[2:5]])
        cls.reciprocal = 2 * np.pi * np.linalg.inv(lattice).T

    def validate(self, raw=None, kpoints=None, reference=None):
        return validate_field(self.raw if raw is None else raw,
                              self.kpoints if kpoints is None else kpoints,
                              self.reciprocal,
                              self.reference if reference is None else reference)

    def test_actual_native_reference_covers_input_mesh(self):
        field, error = self.validate()
        self.assertEqual(field.shape, (144, 7))
        self.assertEqual(error, 0)

    def test_repeated_144_first_bz_coordinates_are_rejected(self):
        malformed = self.raw.copy()
        selected = np.flatnonzero((malformed[:, 4] >= 0) & (malformed[:, 4] < 1)
                                  & (malformed[:, 5] >= 0) & (malformed[:, 5] < 1))
        malformed[selected] = malformed[selected[0]]
        with self.assertRaisesRegex(ValueError, "duplicate plaquette"):
            self.validate(raw=malformed)

    def test_off_mesh_fractional_coordinates_are_rejected(self):
        malformed = self.raw.copy()
        malformed[0, 4] += .01
        malformed[0, :3] = malformed[0, 4:] @ self.reciprocal
        with self.assertRaisesRegex(ValueError, "coordinates do not match"):
            self.validate(raw=malformed)

    def test_cartesian_coordinates_must_match_actual_lattice(self):
        malformed = self.raw.copy()
        malformed[0, 0] += .01
        with self.assertRaisesRegex(ValueError, "Cartesian and fractional"):
            self.validate(raw=malformed)

    def test_map_value_outside_native_precision_is_rejected(self):
        malformed = self.raw.copy()
        selected = np.flatnonzero((malformed[:, 4] >= 0) & (malformed[:, 4] < 1)
                                  & (malformed[:, 5] >= 0) & (malformed[:, 5] < 1))
        malformed[selected[0], 3] = .0001
        with self.assertRaisesRegex(ValueError, "curvature map disagrees"):
            self.validate(raw=malformed)

    def test_duplicate_input_mesh_point_is_rejected(self):
        malformed = self.kpoints.copy()
        malformed[1] = malformed[0]
        with self.assertRaisesRegex(ValueError, "duplicate or missing mesh"):
            self.validate(kpoints=malformed)

    def test_duplicate_reference_coordinate_is_rejected(self):
        malformed = self.reference.copy()
        malformed[1] = malformed[0]
        with self.assertRaisesRegex(ValueError, "reference plaquettes"):
            self.validate(reference=malformed)


if __name__ == "__main__":
    unittest.main()
