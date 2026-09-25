"""Synthetic schema/physics tests; optional local generated-output crosscheck.

No VASP source, POTCAR, or private calculation outputs are embedded here.
"""
import sys
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from procar_projection import (ProcarFormatError, aggregate_groups, parse_procar,
                               validate_state_alignment)


def synthetic_procar(*, nk=2, nb=2, shuffle=False, raw=None):
    """Two-ion, two-orbital joint spin example generated independently."""
    if raw is None:
        raw = np.zeros((2, 2, 4))
        raw[0, 0] = [.2, -.02, .01, .2]
        raw[1, 1] = [.6, .03, -.04, -.6]
    raw = np.asarray(raw)
    rows = ["PROCAR lm decomposed",
            f"# of k-points: {nk} # of bands: {nb} # of ions: 2", ""]
    k_ids = list(range(1, nk+1))
    b_ids = list(range(1, nb+1))
    ions = [1, 2]
    if shuffle:
        k_ids.reverse(); b_ids.reverse(); ions.reverse()
    for ik in k_ids:
        rows.append(f" k-point {ik} : {(ik-1)/nk:.8f} 0.00000000 0.00000000 weight = {1/nk:.8f}")
        for ib in b_ids:
            rows.extend([f"band {ib} # energy {ik+.1*ib:.8f} # occ. {float(ib == 1):.8f}",
                         "ion s p tot"])
            for component in range(4):
                block = raw[:, :, component]
                for ion in ions:
                    values = [*block[ion-1], block[ion-1].sum()]
                    rows.append(str(ion)+" "+" ".join(f"{v:.3f}" for v in values))
                values = [*block.sum(axis=0), block.sum()]
                rows.append("tot "+" ".join(f"{v:.3f}" for v in values))
            rows.append("")
    return "\n".join(rows)+"\n"


class ProcarParserTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name)/"PROCAR"

    def read(self, text=None, **kwargs):
        self.path.write_text(synthetic_procar() if text is None else text)
        defaults = dict(noncollinear=True, spin_to_cartesian=np.eye(3))
        defaults.update(kwargs)
        return parse_procar(self.path, **defaults)

    def test_ids_shapes_raw_components_and_empty_band_projection(self):
        p = self.read()
        self.assertEqual((p.nkpoints, p.nbands, p.nions), (2, 2, 2))
        np.testing.assert_array_equal(p.kpoint_ids, [1, 2])
        np.testing.assert_array_equal(p.band_ids, [1, 2])
        np.testing.assert_array_equal(p.ion_ids, [1, 2])
        self.assertEqual(p.orbital_names, ("s", "p"))
        self.assertEqual(p.ion_orbital.shape, (2, 2, 2, 2, 4))
        self.assertEqual(p.ion_totals.shape, (2, 2, 2, 4))
        self.assertEqual(p.orbital_totals.shape, (2, 2, 2, 4))
        np.testing.assert_allclose(p.state_totals[0, 0], [.8, .01, -.03, -.4])
        self.assertEqual(p.occupations[0, 1], 0)
        self.assertAlmostEqual(p.state_totals[0, 1, 0], .8)
        np.testing.assert_array_equal(p.weights, p.kpoint_weights)

    def test_shuffled_ids_are_canonicalized(self):
        expected = self.read()
        shuffled = self.read(synthetic_procar(shuffle=True))
        for name in ("kpoints", "energies_eV", "occupations", "ion_orbital", "ion_totals"):
            np.testing.assert_array_equal(getattr(shuffled, name), getattr(expected, name))

    def test_adjacent_fixed_width_negative_coordinates_and_exponents(self):
        text = synthetic_procar().replace(
            "0.50000000 0.00000000 0.00000000", "0.14285714-0.00000000-0.25000000")
        p = self.read(text)
        np.testing.assert_array_equal(p.kpoints[1], [.14285714, -0., -.25])
        self.assertTrue(np.signbit(p.kpoints[1, 1]))
        text = text.replace("0.14285714-0.00000000-0.25000000", "1.4285714D-1-0.0D+0-2.5D-1")
        np.testing.assert_array_equal(self.read(text).kpoints[1], p.kpoints[1])
        with self.assertRaises(ProcarFormatError):
            self.read(text.replace("1.4285714D-1-0.0D+0-2.5D-1", "0.142857140.00000000 0.00000000"))

    def test_joint_spin_is_not_group_charge_times_whole_spin(self):
        p = self.read()
        groups = aggregate_groups(p, {"left": [1], "right": [2], "all": [1, 2]}, axis_cartesian=[0, 0, 1])
        left, right, total = (groups[key] for key in ("left", "right", "all"))
        np.testing.assert_allclose(left.p_plus, .2)
        np.testing.assert_allclose(right.p_plus, 0)
        np.testing.assert_allclose(total.q, .8)  # never force PAW charge to one
        np.testing.assert_allclose(left.p_plus+right.p_plus, total.p_plus)
        np.testing.assert_allclose(left.p_minus+right.p_minus, total.p_minus)
        naive = left.q*total.p_plus/total.q
        np.testing.assert_allclose(naive, .05)
        self.assertGreater(float(np.max(abs(left.p_plus-naive))), .1)
        # m is Pauli weight, with no extra 1/2 until physical S/hbar is needed.
        np.testing.assert_allclose(left.m_axis, .2)

    def test_small_negative_joint_weight_is_preserved(self):
        raw = np.zeros((2, 2, 4)); raw[0, 0] = [.100, 0, 0, .101]
        p = self.read(synthetic_procar(raw=raw))
        g = aggregate_groups(p, {"one": [1]}, axis_cartesian=[0, 0, 1])["one"]
        np.testing.assert_allclose(g.p_minus, -.0005)
        np.testing.assert_allclose(g.p_plus+g.p_minus, g.q)

    def test_printed_totals_are_retained_within_rounding_bounds(self):
        text = synthetic_procar().replace("1 0.200 0.000 0.200", "1 0.200 0.000 0.201")
        p = self.read(text)
        self.assertAlmostEqual(p.ion_orbital[0, 0, 0, :, 0].sum(), .2)
        self.assertAlmostEqual(p.ion_totals[0, 0, 0, 0], .201)
        self.assertAlmostEqual(p.state_totals[0, 0, 0], .8)
        g = aggregate_groups(p, {"one": [1]}, axis_cartesian=[0, 0, 1])["one"]
        self.assertAlmostEqual(g.q[0, 0], .201)
        self.assertAlmostEqual(p.rounding_diagnostics["max_orbital_sum_vs_ion_total"], .001)

    def test_explicit_spin_frame_and_unit_axis(self):
        rotation = np.array([[0., 0, 1], [0, 1, 0], [-1, 0, 0]])
        p = self.read(spin_to_cartesian=rotation)
        g = aggregate_groups(p, {"one": [1]}, axis_cartesian=[1, 0, 0])["one"]
        np.testing.assert_allclose(g.m_axis, .2)
        np.testing.assert_allclose(g.m_cartesian[0, 0], [.2, .01, .02])
        np.testing.assert_allclose(g.p_plus, .2)
        for matrix in (None, np.eye(2), np.diag([1., 1., -1.]), 2*np.eye(3), np.full((3, 3), np.nan)):
            with self.subTest(matrix=matrix), self.assertRaises(ValueError):
                self.read(spin_to_cartesian=matrix)
        with self.assertRaises(TypeError):
            parse_procar(self.path, noncollinear=True)
        with self.assertRaises(TypeError):
            aggregate_groups(p, {"one": [1]})
        for axis in ([0, 0, 0], [0, 0, 2], [np.nan, 0, 1], [0, 1]):
            with self.subTest(axis=axis), self.assertRaises(ValueError):
                aggregate_groups(p, {"one": [1]}, axis_cartesian=axis)

    def test_explicit_noncollinear_mode_required(self):
        for mode in (False, None, "True", 1):
            with self.subTest(mode=mode), self.assertRaises(ValueError):
                self.read(noncollinear=mode)
        with self.assertRaises(TypeError):
            parse_procar(self.path, spin_to_cartesian=np.eye(3))

    def test_groups_require_nonempty_unique_one_based_integers(self):
        p = self.read()
        for groups in ({}, {"": [1]}, {"a": []}, {"a": [0]}, {"a": [3]},
                       {"a": [1, 1]}, {"a": [1.]}, {"a": [True]}, {"a": 1}):
            with self.subTest(groups=groups), self.assertRaises(ValueError):
                aggregate_groups(p, groups, axis_cartesian=[0, 0, 1])
        g = aggregate_groups(p, {"one": [np.int64(1)]}, axis_cartesian=[0, 0, 1])["one"]
        self.assertEqual(g.ion_ids, (1,))

    def test_truncated_missing_and_collinear_blocks_fail(self):
        text = synthetic_procar(nk=1, nb=1)
        nonblank = [line for line in text.splitlines() if line.strip()]
        for changed in ("\n".join(nonblank[:-1]), "\n".join(nonblank[:8]),
                        "\n".join(nonblank[:9]+nonblank[10:]),
                        text.replace("ion s p tot", "spin component 1\nion s p tot")):
            with self.subTest(changed=changed[:60]), self.assertRaises(ProcarFormatError):
                self.read(changed)

    def test_duplicate_and_out_of_range_ids_fail(self):
        text = synthetic_procar()
        for changed in (text.replace("k-point 2", "k-point 1"), text.replace("k-point 2", "k-point 3"),
                        text.replace("band 2", "band 1"), text.replace("band 2", "band 3"),
                        text.replace("2 0.000 0.600 0.600", "1 0.000 0.600 0.600", 1),
                        text.replace("2 0.000 0.600 0.600", "3 0.000 0.600 0.600", 1)):
            with self.subTest(changed=changed[:60]), self.assertRaises(ProcarFormatError):
                self.read(changed)

    def test_extra_blocks_changed_orbitals_and_corrupt_totals_fail(self):
        text = synthetic_procar()
        for changed in (text+"tot 0 0 0\n", text.replace("ion s p tot", "ion s s tot", 1),
                        text.replace("ion s p tot", "ion s d tot", 1),
                        text.replace("1 0.200 0.000 0.200", "1 0.200 0.000 0.300", 1),
                        text.replace("1 -0.020 0.000 -0.020", "ion s p tot\n1 -0.020 0.000 -0.020", 1)):
            with self.subTest(changed=changed[:60]), self.assertRaises(ProcarFormatError):
                self.read(changed)

    def test_nonfinite_header_and_projection_values_fail(self):
        text = synthetic_procar()
        for old, new in (("energy 1.10000000", "energy nan"), ("occ. 1.00000000", "occ. inf"),
                         ("weight = 0.50000000", "weight = nan"),
                         (": 0.50000000", ": nan"), ("1 0.200", "1 nan"), ("1 0.200", "1 *****")):
            with self.subTest(old=old), self.assertRaises(ProcarFormatError):
                self.read(text.replace(old, new, 1))

    def test_alignment_complete_ids_energy_coordinates_and_weights(self):
        p = self.read()
        result = validate_state_alignment(p, p.kpoints+[1, -2, 0], p.energies_eV+4e-9,
                                          kpoint_ids=[1, 2], band_ids=[1, 2], weights=p.weights)
        self.assertEqual(result["max_k_fractional_delta"], 0)
        self.assertLess(result["max_energy_delta_eV"], 5.1e-9)
        self.assertEqual(result["max_weight_delta"], 0)
        validate_state_alignment(p, p.kpoints+4e-9, p.energies_eV)
        for kw in ({"band_ids": [2, 1]}, {"kpoint_ids": [2, 1]}, {"weights": [1, 0]},
                   {"weights": [1]}, {"weights": [np.nan, 0]}, {"k_atol": -1}, {"energy_atol_eV": np.inf}):
            with self.subTest(kw=kw), self.assertRaises(ValueError):
                validate_state_alignment(p, p.kpoints, p.energies_eV, **kw)
        for points, energies in ((p.kpoints+.001, p.energies_eV), (p.kpoints, p.energies_eV+1e-8),
                                 (p.kpoints, p.energies_eV[:, :1]), (p.kpoints[:1], p.energies_eV),
                                 (p.kpoints*np.nan, p.energies_eV), (p.kpoints, p.energies_eV*np.inf)):
            with self.subTest(shape=np.shape(energies)), self.assertRaises(ValueError):
                validate_state_alignment(p, points, energies)



if __name__ == "__main__":
    unittest.main()
