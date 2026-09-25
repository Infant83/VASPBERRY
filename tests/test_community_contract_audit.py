"""Independent point-curvature contract and oriented transport edge cases."""
import copy
from dataclasses import replace
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from berry_data import (CurvatureData, base_metadata, hall_spectrum, read_curvature,
                       regions_from_spec, validate, write_curvature)
from exported_matrix_kubo import MatrixContractError, DegenerateBandError, berry_curvature


def fixture(reciprocal=None):
    reciprocal = (2*np.pi*np.eye(3) if reciprocal is None
                  else np.asarray(reciprocal, dtype=np.float64))
    q = np.array([[0., 0., 0.], [0., .5, 0.], [.5, 0., 0.], [.5, .5, 0.]])
    omega = np.zeros((4, 2, 3))
    omega[:, 0] = [.2, -.3, .7]
    omega[:, 1] = -omega[:, 0]
    meta = base_metadata(
        source_nbands=2, spin_multiplicity=1, method='analytic_contract_fixture',
        energy_reference='fixture energy zero',
        sampling={'kind': 'uniform_full_2d', 'mesh': [2, 2], 'plane_axes': [0, 1]},
        source_operator={'kind': 'analytic_contract_fixture', 'accuracy_status': 'validated'},
        provenance={'fixture': 'test_community_contract_audit.py'})
    return CurvatureData(meta, np.arange(1, 5, dtype=np.int64), np.array([1, 2], np.int64),
        np.array([1, 2], np.int64), q, np.full(4, .25), np.tile([-1., 2.], (4, 1)),
        omega, np.ones((4, 2), bool), np.full((4, 2), 3.),
        2*np.pi*np.linalg.inv(reciprocal).T, reciprocal)


def total_at_zero(data, **kwargs):
    rows, _ = hall_spectrum(data, [0.], [0.], mu_reference=0., **kwargs)
    return next(r['sigma_e2_over_h'] for r in rows if r['region'] == 'total' and r['band_id'] == 0)


class IndependentCommunityContractTests(unittest.TestCase):
    def test_known_spin_metadata_cannot_double_a_channel_or_spinor(self):
        data = fixture()
        invalid = [
            {'source_nspin': 2, 'spin_channel_1based': 1, 'spin_multiplicity': 2},
            {'source_nspin': 1, 'spin_channel_1based': 2},
            {'source_nspin': True, 'spin_channel_1based': 1},
            {'source_nspin': 2, 'spin_channel_1based': True},
            {'source_nspin': 2}, {'spin_channel_1based': 1},
            {'spinor_components': 2, 'spin_multiplicity': 2},
            {'spinor_components': True},
        ]
        for fields in invalid:
            with self.subTest(fields=fields):
                with self.assertRaises(MatrixContractError):
                    validate(replace(data, metadata=dict(data.metadata, **fields)))
        for fields in (
            {'source_nspin': 2, 'spin_channel_1based': 2},
            {'source_nspin': 1, 'spin_channel_1based': 1, 'spinor_components': 2},
            {'source_nspin': 1, 'spin_channel_1based': 1, 'spinor_components': 1, 'spin_multiplicity': 2},
            {'spin_multiplicity': 2},  # Unknown source representation stays explicitly user-declared.
        ):
            validate(replace(data, metadata=dict(data.metadata, **fields)))

    def test_delta_hall_keeps_small_signal_above_large_filled_baseline(self):
        data = fixture()
        e = np.tile([-1000., .1, 1000.], (4, 1))
        omega = np.zeros((4, 3, 3)); omega[:, :, 2] = [1e20, .2, 0.]
        data = replace(data, metadata=dict(data.metadata, source_nbands=3),
            band_ids=np.array([1, 2, 3], np.int64), intermediate_band_ids=np.array([1, 2, 3], np.int64),
            energies_eV=e, omega_A2=omega, valid_nondegenerate=np.ones((4, 3), bool),
            min_gap_eV=np.full((4, 3), 999.9))
        rows, _ = hall_spectrum(data, [.2], [0., 300.], mu_reference=0., band_resolved=True, mu_chunk=1)
        for t in (0., 300.):
            df = 1. if t == 0 else 1/(1+np.exp(-.1/(8.617333262145e-5*t)))-1/(1+np.exp(.1/(8.617333262145e-5*t)))
            expected = -2*np.pi*.2*df
            for band in (0, 2):
                row = next(r for r in rows if r['temperature_K'] == t and r['region'] == 'total' and r['band_id'] == band)
                self.assertAlmostEqual(row['delta_sigma_e2_over_h'], expected, places=12)
                self.assertAlmostEqual(row['delta_electrons_per_cell'], df, places=12)

    def test_legacy_import_normalizes_once_and_restores_unsorted_gap_order(self):
        from vaspberry_kubo import import_legacy
        data = fixture()
        energies = np.tile([3., 0., 1.], (4, 1))
        wave = SimpleNamespace(energies=energies, kpoints=data.kpoints_fractional,
            header=SimpleNamespace(ispin=1, lattice=data.lattice_A, reciprocal=data.reciprocal_inv_A),
            coefficients=lambda *args: np.zeros((1, 1)))
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); wave_path = root/'WAVECAR'; wave_path.write_bytes(b'fixture')
            for convention, expected in (('legacy-double', 4.), ('physical', 8.)):
                csv_path = root/(convention+'.csv')
                column = 'omega_legacy_A2' if convention == 'legacy-double' else 'omega_z_A2'
                lines = ['k_index,band,kx_frac,ky_frac,kz_frac,energy_eV,min_gap_eV,'+column]
                for k, q in enumerate(data.kpoints_fractional, 1):
                    for b in (1, 2, 3):
                        lines.append(','.join(map(str, [k, b, *q, energies[k-1, b-1], 10., 8.])))
                csv_path.write_text('\n'.join(lines)+'\n')
                args = SimpleNamespace(csv=csv_path, wavecar=wave_path, output_dir=root/convention,
                    spin=1, spinor_components=1, spin_multiplicity=1, normalization=convention,
                    degeneracy_threshold_eV=1e-8, energy_reference='fixture zero', mesh=[2, 2], plane_axes=[0, 1])
                with patch('vaspberry_kubo.Wavecar', return_value=wave):
                    import_legacy(args)
                imported = read_curvature(args.output_dir)
                np.testing.assert_array_equal(imported.omega_A2[:, :, 2], expected)
                np.testing.assert_array_equal(imported.min_gap_eV, np.tile([2., 1., 1.], (4, 1)))
                self.assertEqual(imported.metadata['gap_validation']['max_abs_difference_reported_vs_energies_eV'], 9.)
                self.assertEqual(imported.metadata.get('source_nspin'), 1)
                self.assertEqual(imported.metadata.get('spin_channel_1based'), 1)
                self.assertEqual(imported.metadata.get('spinor_components'), 1)
                if convention == 'legacy-double':
                    csv_path.write_text('# STANDARD_MINUS_TWO_IM\n'+csv_path.read_text())
                    args.output_dir = root/'double-division-refused'
                    with patch('vaspberry_kubo.Wavecar', return_value=wave):
                        with self.assertRaises(MatrixContractError):
                            import_legacy(args)
                    self.assertFalse(args.output_dir.exists())

    def test_legacy_csv_gap_cannot_hide_exact_wavecar_degeneracy(self):
        from vaspberry_kubo import import_legacy
        data = fixture()
        wave = SimpleNamespace(energies=np.zeros((4, 2)), kpoints=data.kpoints_fractional,
            header=SimpleNamespace(ispin=1, lattice=data.lattice_A, reciprocal=data.reciprocal_inv_A),
            coefficients=lambda *args: np.zeros((1, 1)))
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); csv_path = root/'legacy.csv'; wave_path = root/'WAVECAR'
            wave_path.write_bytes(b'fixture provenance only; Wavecar object is mocked')
            lines = ['k_index,band,kx_frac,ky_frac,kz_frac,energy_eV,min_gap_eV,omega_legacy_A2']
            for k, q in enumerate(data.kpoints_fractional, 1):
                for b in (1, 2):
                    lines.append(','.join(map(str, [k, b, *q, 0., 1., .5])))
            csv_path.write_text('\n'.join(lines)+'\n')
            args = SimpleNamespace(csv=csv_path, wavecar=wave_path, output_dir=root/'out',
                spin=1, spinor_components=1, spin_multiplicity=1, normalization='legacy-double',
                degeneracy_threshold_eV=1e-8, energy_reference='fixture zero', mesh=[2, 2], plane_axes=[0, 1])
            with patch('vaspberry_kubo.Wavecar', return_value=wave):
                try:
                    import_legacy(args)
                except MatrixContractError:
                    return  # Refusing inconsistent source gap metadata is also valid.
            imported = read_curvature(args.output_dir)
            self.assertFalse(imported.valid_nondegenerate.any())
            np.testing.assert_array_equal(imported.min_gap_eV, 0.)
            self.assertEqual(imported.metadata['gap_validation']['max_abs_difference_reported_vs_energies_eV'], 1.)

    def test_m_subset_cannot_hide_exported_degenerate_partner(self):
        from test_exported_matrix_kubo import model_data
        data = model_data([0., 0., 2.], np.zeros((3, 3, 3), complex))
        kwargs = dict(interest_band_ids=[1], intermediate_band_ids=[1, 3],
                      degeneracy_threshold_eV=1e-8, allow_experimental=True)
        with self.assertRaises(DegenerateBandError):
            berry_curvature(data, **kwargs)
        result = berry_curvature(data, degeneracy_policy='mask', **kwargs)
        self.assertFalse(result.valid_nondegenerate.any())
        self.assertTrue(np.isnan(result.omega_A2).all())
        self.assertEqual(result.min_gap_eV[0, 0], 0.)
        self.assertEqual(result.diagnostics['degenerate_pairs'][0]['m_band'], 2)
        self.assertTrue(result.diagnostics['all_source_band_energies_available'])

    def test_masked_state_may_lack_matrix_but_valid_state_may_not(self):
        from test_exported_matrix_kubo import model_data
        data = model_data([[0., 0., 2.], [0., 1., 2.]], np.zeros((2, 3, 3, 3), complex))
        matrix = data.D_eVA.copy(); covered = data.coverage.copy()
        matrix[0, :, 0, 2] = matrix[0, :, 2, 0] = complex(np.nan, np.nan)
        covered[0, :, 0, 2] = covered[0, :, 2, 0] = False
        kwargs = dict(interest_band_ids=[1], intermediate_band_ids=[1, 3],
                      degeneracy_threshold_eV=1e-8, allow_experimental=True, degeneracy_policy='mask')
        result = berry_curvature(replace(data, D_eVA=matrix, coverage=covered), **kwargs)
        np.testing.assert_array_equal(result.valid_nondegenerate[:, 0], [False, True])
        self.assertTrue(np.isnan(result.omega_A2[0]).all())
        self.assertTrue(np.isfinite(result.omega_A2[1]).all())
        matrix[1, :, 0, 2] = matrix[1, :, 2, 0] = complex(np.nan, np.nan)
        covered[1, :, 0, 2] = covered[1, :, 2, 0] = False
        with self.assertRaises(MatrixContractError):
            berry_curvature(replace(data, D_eVA=matrix, coverage=covered), **kwargs)

    def test_gap_scope_uses_known_energies_and_reports_unknown_source_bands(self):
        from test_exported_matrix_kubo import model_data
        data = model_data([0., .1, 2.], np.zeros((3, 3, 3), complex))
        kwargs = dict(interest_band_ids=[1], intermediate_band_ids=[1, 3],
                      degeneracy_threshold_eV=1e-8, allow_experimental=True)
        result = berry_curvature(data, **kwargs)
        self.assertEqual(result.min_gap_eV[0, 0], .1)
        self.assertTrue(result.diagnostics['all_source_band_energies_available'])
        keep = [0, 2]
        cropped = replace(data, band_ids=data.band_ids[keep], energies_eV=data.energies_eV[:, keep],
                          D_eVA=data.D_eVA[:, :, keep, :][:, :, :, keep],
                          coverage=data.coverage[:, :, keep, :][:, :, :, keep])
        result = berry_curvature(cropped, **kwargs)
        self.assertEqual(result.min_gap_eV[0, 0], 2.)
        self.assertFalse(result.diagnostics['all_source_band_energies_available'])

    def test_oriented_skew_plane_uses_all_axial_components_and_not_integer_gate(self):
        data = fixture(2*np.pi*np.array([[1., .3, .4], [0., 1., .2], [0., 0., 1.]]))
        expected = -np.dot(data.omega_A2[0, 0], np.cross(*data.reciprocal_inv_A[:2]))/(2*np.pi)
        self.assertAlmostEqual(total_at_zero(data), expected, places=12)
        self.assertGreater(abs(expected-round(expected)), .01)
        reverse_meta = copy.deepcopy(data.metadata)
        reverse_meta['sampling']['plane_axes'] = [1, 0]
        self.assertAlmostEqual(total_at_zero(replace(data, metadata=reverse_meta)), -expected, places=12)

    def test_shuffled_global_band_and_k_storage_preserve_result(self):
        data = fixture()
        k, b = np.array([2, 0, 3, 1]), np.array([1, 0])
        shuffled = replace(data, k_ids=data.k_ids[k], band_ids=data.band_ids[b],
            kpoints_fractional=data.kpoints_fractional[k], weights=data.weights[k],
            energies_eV=data.energies_eV[k][:, b], omega_A2=data.omega_A2[k][:, b],
            valid_nondegenerate=data.valid_nondegenerate[k][:, b], min_gap_eV=data.min_gap_eV[k][:, b])
        self.assertAlmostEqual(total_at_zero(shuffled), total_at_zero(data), places=12)

    def test_available_xy_only_is_rejected_for_tilted_plane(self):
        data = fixture(2*np.pi*np.array([[1., 0., .4], [0., 1., 0.], [0., 0., 1.]]))
        meta = copy.deepcopy(data.metadata); meta['available_components'] = [False, False, True]
        omega = data.omega_A2.copy(); omega[:, :, :2] = np.nan
        with self.assertRaises(MatrixContractError):
            total_at_zero(replace(data, metadata=meta, omega_A2=omega))

    def test_point_source_and_duplicate_endpoint_cannot_become_full_bz_by_weight_sum(self):
        data = fixture()
        meta = copy.deepcopy(data.metadata); meta['sampling'] = {'kind': 'points'}
        validate(replace(data, metadata=meta))
        with self.assertRaises(MatrixContractError):
            total_at_zero(replace(data, metadata=meta))
        q = data.kpoints_fractional.copy(); q[3] = q[0]+[1., 0., 0.]
        with self.assertRaises(ValueError):
            total_at_zero(replace(data, kpoints_fractional=q))

    def test_partial_n_window_and_occupied_top_need_explicit_partial_scope(self):
        data = fixture()
        reduced = replace(data, band_ids=data.band_ids[:1], energies_eV=data.energies_eV[:, :1],
            omega_A2=data.omega_A2[:, :1], valid_nondegenerate=data.valid_nondegenerate[:, :1],
            min_gap_eV=data.min_gap_eV[:, :1])
        for d, mu in ((reduced, 0.), (data, 2.)):
            with self.assertRaises(MatrixContractError):
                hall_spectrum(d, [mu], [0.], mu_reference=mu)
            _, meta = hall_spectrum(d, [mu], [0.], mu_reference=mu, allow_partial_bands=True)
            self.assertEqual(meta['scope'], 'partial_band_contribution')

    def test_equivalent_third_axis_center_representatives_select_same_region(self):
        data = fixture()
        def masks(z):
            return regions_from_spec(data, {'regions': [{'name': 'valley',
                'center_fractional': [0., 0., z], 'radius_inv_A': .2}]})
        reference = masks(0.)
        shifted = masks(1.)
        for key in reference:
            np.testing.assert_array_equal(shifted[key], reference[key])

    def test_geometric_circle_overlap_is_rejected_between_sampled_points(self):
        data = fixture()
        # Each disk contains only its own coarse-mesh point, but their interiors
        # intersect between x=0 and x=1/2. A disjoint-disk declaration is false.
        spec = {'regions': [
            {'name': 'A', 'center_fractional': [0., 0., 0.], 'radius_inv_A': .26*2*np.pi},
            {'name': 'B', 'center_fractional': [.5, 0., 0.], 'radius_inv_A': .26*2*np.pi}]}
        with self.assertRaises(MatrixContractError):
            regions_from_spec(data, spec)

    def test_mixed_id_and_circle_sampled_overlap_remains_rejected(self):
        data = fixture()
        circle = {'name': 'circle', 'center_fractional': [0., 0., 0.], 'radius_inv_A': .2}
        ids = {'name': 'ids', 'k_ids': [1]}
        for entries in ([circle, ids], [ids, circle]):
            with self.assertRaises(MatrixContractError):
                regions_from_spec(data, {'regions': entries})

    def test_exact_zero_gap_cannot_claim_valid_individual_band(self):
        data = fixture()
        gaps = data.min_gap_eV.copy(); gaps[0, 0] = 0.
        with self.assertRaises(MatrixContractError):
            validate(replace(data, min_gap_eV=gaps))

    def test_schema_version_requires_integer_not_bool_or_float(self):
        data = fixture()
        for version in (True, 1.0):
            with self.subTest(version=version):
                meta = copy.deepcopy(data.metadata); meta['version'] = version
                with self.assertRaises(MatrixContractError):
                    validate(replace(data, metadata=meta))

    def test_nonobject_json_is_a_clean_contract_error(self):
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)/'curvature'
            write_curvature(out, fixture())
            (out/'curvature.json').write_text(json.dumps([]))
            with self.assertRaises(MatrixContractError):
                read_curvature(out)

    def test_malformed_sampling_and_array_types_are_clean_contract_errors(self):
        data = fixture()
        meta = copy.deepcopy(data.metadata); meta['sampling'] = []
        cases = [replace(data, metadata=meta),
                 replace(data, valid_nondegenerate=data.valid_nondegenerate.tolist()),
                 replace(data, min_gap_eV=data.min_gap_eV.tolist()),
                 replace(data, omega_A2=data.omega_A2.tolist())]
        for malformed in cases:
            with self.subTest(malformed=type(malformed.metadata.get('sampling'))):
                with self.assertRaises(MatrixContractError):
                    validate(malformed)


if __name__ == '__main__':
    unittest.main()
