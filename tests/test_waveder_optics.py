"""Circular PAW transition invariants and standard-run adapter contract."""
import contextlib
import csv
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import waveder_optics as optics
import waveder_hall as hall
import vaspberry_kubo as cli
import test_waveder_hall as fixtures


class CircularPhysicsTests(unittest.TestCase):
    def setUp(self):
        self.e = np.tile([-1., -1., 2., 2., 5.], (2, 1))
        self.q = np.array([[1/3, 1/3, 0.], [-1/3, -1/3, 0.]])
        self.rng = np.random.default_rng(1846)
        self.c = self.rng.normal(size=(2, 3, 5, 2))+1j*self.rng.normal(size=(2, 3, 5, 2))
        self.options = dict(occupied=2, initial=[1, 2], final=[3, 4],
                            photon_eV=np.linspace(2., 4., 201), sigma_eV=.05)

    def tables(self, c=None, **options):
        return optics.transition_tables(self.e, self.c if c is None else c, self.q, np.array([.2, .8]),
                                        **dict(self.options, **options))

    def test_both_complete_multiplets_are_unitary_invariant(self):
        expected, expected_spectrum, _ = self.tables()
        changed = self.c.copy()
        for k in range(2):
            ui = np.linalg.qr(self.rng.normal(size=(2, 2))+1j*self.rng.normal(size=(2, 2)))[0]
            uf = np.linalg.qr(self.rng.normal(size=(2, 2))+1j*self.rng.normal(size=(2, 2)))[0]
            changed[k, :, 2:4, :] = uf.conj().T@changed[k, :, 2:4, :]@ui
        actual, actual_spectrum, _ = self.tables(changed)
        for name in expected:
            np.testing.assert_allclose(actual[name], expected[name], atol=3e-13)
        for name in expected_spectrum:
            np.testing.assert_allclose(actual_spectrum[name], expected_spectrum[name], atol=5e-13)

    def test_time_reversal_exchanges_helicities_at_opposite_k(self):
        self.c[1] = self.c[0].conj()
        transitions, spectra, _ = self.tables()
        self.assertAlmostEqual(transitions['I_plus_A2'][0], transitions['I_minus_A2'][1])
        self.assertAlmostEqual(transitions['eta'][0], -transitions['eta'][1])
        nk = len(self.options['photon_eV'])
        np.testing.assert_allclose(spectra['I_plus_A2_per_eV'][:nk], spectra['I_minus_A2_per_eV'][nk:])

    def test_circular_sum_is_transverse_linear_strength_and_beam_reversal_swaps(self):
        e1, e2, beam = optics.polarization_frame([1, 2, 3], [1, .5, 0])
        block = self.c[0, :, 2:4, :]
        p, m = optics.circular_strength(block, e1, e2)
        linear = sum(np.linalg.norm(np.einsum('a,aij->ij', direction, block))**2 for direction in (e1, e2))
        self.assertAlmostEqual(p+m, linear, places=12)
        r1, r2, _ = optics.polarization_frame(-beam, e1)
        reverse = optics.circular_strength(block, r1, r2)
        np.testing.assert_allclose(reverse, [m, p], atol=2e-14)

    def test_pure_channels_pin_nonconjugated_bra_ket_helicity_sign(self):
        block = np.array([1, 1j, 0.]).reshape(3, 1, 1)
        np.testing.assert_allclose(optics.circular_strength(block, [1, 0, 0], [0, 1, 0]), [0, 2], atol=1e-15)
        np.testing.assert_allclose(optics.circular_strength(block.conj(), [1, 0, 0], [0, 1, 0]), [2, 0], atol=1e-15)

    def test_centroid_group_strength_integrates_to_gaussian_area(self):
        self.e[:] = [-1.0002, -.9998, 1.9997, 2.0003, 5.]
        photons = np.linspace(0., 6., 6001)
        transitions, spectra, meta = self.tables(photon_eV=photons)
        np.testing.assert_allclose(transitions['transition_eV'], 3., atol=1e-15)
        np.testing.assert_allclose(transitions['transition_min_eV'], 2.9995)
        self.assertAlmostEqual(meta['degeneracy']['maximum_group_energy_spread_eV'], .0006)
        for channel in ('plus', 'minus'):
            values = spectra[f'I_{channel}_A2_per_eV'].reshape(2, -1)
            areas = np.sum((values[:, 1:]+values[:, :-1])*.5*np.diff(photons), axis=1)
            np.testing.assert_allclose(areas, transitions[f'I_{channel}_A2'], rtol=2e-14)
        self.assertFalse(meta['weighting']['k_weight_applied'])

    def test_zero_and_weak_channels_have_explicit_eta_mask(self):
        self.c[:] = 0
        with np.errstate(all='raise'):
            transitions, spectra, _ = self.tables()
        for result in (transitions, spectra):
            self.assertFalse(result['eta_valid'].any())
            self.assertTrue(np.isnan(result['eta']).all())
        eta, valid = optics.contrast(np.array([1., 1e-12, 0.]), np.zeros(3), 1e-10)
        np.testing.assert_array_equal(valid, [True, False, False])
        self.assertEqual(eta[0], 1.)

    def test_initial_and_final_cut_groups_and_unobserved_top_boundary_reject(self):
        for option in ({'initial': [2, 2]}, {'final': [3, 3]}):
            with self.subTest(option=option), self.assertRaisesRegex(ValueError, 'cuts a degenerate group'):
                self.tables(**option)
        with self.assertRaisesRegex(ValueError, 'source band above'):
            self.tables(final=[3, 5])
        self.e[:] = [-1., -1., -.999, 2., 5.]
        with self.assertRaisesRegex(ValueError, 'occupied/empty boundary'):
            self.tables()

    def test_transitive_clusters_use_complete_stored_axis(self):
        e = np.array([0., .0015, .003, 1.])
        with self.assertRaisesRegex(ValueError, 'cuts a degenerate group'):
            optics.complete_groups(e, [1, 2], .002)
        groups = optics.complete_groups(e, [1, 3], .002)
        self.assertEqual(len(groups), 1)
        np.testing.assert_array_equal(groups[0], [0, 1, 2])

    def test_bad_geometry_threshold_and_photon_grids_reject(self):
        for options in ({'beam': [0, 0, 0]}, {'axis': [0, 0, 1]},
                        {'degeneracy_threshold_eV': .001}, {'sigma_eV': 0.},
                        {'photon_eV': [1., 1.]}, {'photon_eV': [-1.]},
                        {'relative_intensity_floor': 1.}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                self.tables(**options)


class StandardOpticalRunTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.WaveDerHallTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        f = self.fixture
        f.energies = np.tile([-1., -1., 2., 4.], (1, 4, 1))
        f.occupations = np.tile([.99999, .99998, .00003, 0.], (1, 4, 1))
        f.c = np.zeros((1, 4, 3, 4, 2), dtype=np.complex64)
        f.c[0, :, 0, 2] = [1., 2.]; f.c[0, :, 1, 2] = [1j, 2j]
        f.c[0, :, :, :2, :] = 9999+777j  # unused occupied block cannot affect optics
        f.write_run()
        self.options = dict(occupied=2, spin=1, spinor_components=2, energy_reference='fixture unchanged VASP zero',
            initial=[1, 2], final=[3, 3], photon_eV=np.linspace(2., 4., 101), sigma_eV=.05)

    def run_api(self, **options):
        return optics.waveder_optics(self.fixture.run, **dict(self.options, **options))

    def arguments(self, directory):
        return ['--run-dir', str(self.fixture.run), '--occupied', '2', '--spinor-components', '2',
            '--energy-reference', 'fixture zero', '--initial', '1', '2', '--final', '3', '3',
            '--photon-min', '2', '--photon-max', '4', '--photon-num', '101',
            '--output-dir', str(self.fixture.root/directory)]

    def test_real_binary_parser_lower_block_orientation_and_hall_helicity_identity(self):
        transitions, _, meta = self.run_api()
        np.testing.assert_allclose(transitions['I_plus_A2'], 0., atol=1e-15)
        np.testing.assert_allclose(transitions['I_minus_A2'], 10., atol=1e-14)
        np.testing.assert_allclose(transitions['eta'], -1.)
        # For the same stored occupied-empty blocks, I+ - I- = occupied Ω_xy.
        curvature = hall.occupied_curvature(self.fixture.c[0], 2)
        np.testing.assert_allclose(transitions['I_plus_A2']-transitions['I_minus_A2'], curvature[:, 2], atol=1e-14)
        self.assertEqual(meta['source_nbands'], 4)
        self.assertEqual(meta['source_ndbands'], 2)
        self.assertIn('cannot be independently authenticated', meta['same_run_association']['limitation'])

    def test_path_zero_weights_allowed_without_relaxing_hall_uniform_guard(self):
        p = self.fixture.run/'OUTCAR'
        p.write_text(p.read_text().replace('0.250\n', '0.000\n'))
        transitions, spectra, _ = self.run_api()
        self.assertTrue(np.all(transitions['source_k_weight'] == 0))
        np.testing.assert_allclose(transitions['I_minus_A2'], 10.)
        self.assertGreater(spectra['I_minus_A2_per_eV'].max(), 0.)
        with self.assertRaisesRegex(ValueError, 'not uniform'):
            self.fixture.spectrum()

    def test_same_state_mismatch_and_unsupported_optical_branch_reject(self):
        p = self.fixture.run/'OUTCAR'; original = p.read_text()
        for before, after in [('2.0000', '2.1000'), ('LNABLA = F', 'LNABLA = T'),
                              ('vasp.5.4.4', 'vasp.6.5.0')]:
            p.write_text(original.replace(before, after))
            with self.subTest(after=after), self.assertRaises(ValueError):
                self.run_api()
        p.write_text(original)
        with self.assertRaisesRegex(ValueError, 'NELECT'):
            self.run_api(occupied=1)

    def test_source_changes_during_analysis_are_rejected(self):
        original = optics.transition_tables
        def changing_source(*args, **kwargs):
            result = original(*args, **kwargs)
            p = self.fixture.run/'OUTCAR'
            p.write_text(p.read_text()+'source changed during optical analysis\n')
            return result
        with patch.object(optics, 'transition_tables', side_effect=changing_source), \
             self.assertRaisesRegex(ValueError, 'source files changed during optical calculation'):
            self.run_api()

    def test_cli_formats_share_columns_values_and_metadata_and_prevent_overwrite(self):
        with contextlib.redirect_stdout(io.StringIO()) as stdout:
            cli.main(['waveder-optics', *self.arguments('optical')])
        self.assertEqual(json.loads(stdout.getvalue())['schema'], 'vaspberry.circular-transition-strength')
        out = self.fixture.root/'optical'
        meta = json.loads((out/'optical.json').read_text())
        self.assertEqual(len(meta['output_sha256']), 6)
        for name in ('transitions', 'spectra'):
            with (out/f'{name}.csv').open() as f:
                rows = list(csv.DictReader(f))
            with (out/f'{name}.dat').open() as f:
                lines = f.readlines(); lines[0] = lines[0].removeprefix('# ')
                dat_rows = list(csv.DictReader(lines, delimiter='\t'))
            self.assertEqual(rows, dat_rows)
            with np.load(out/f'{name}.npz', allow_pickle=False) as data:
                self.assertEqual(list(rows[0]), data.files)
                for key in data.files:
                    expected = [r[key] == 'True' for r in rows] if data[key].dtype == bool else [float(r[key]) for r in rows]
                    np.testing.assert_allclose(expected, data[key], equal_nan=True)
        before = (out/'optical.json').read_bytes()
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            cli.main(['waveder-optics', *self.arguments('optical')])
        self.assertEqual((out/'optical.json').read_bytes(), before)

    def test_standalone_and_composable_output_match_and_choose_formats(self):
        with contextlib.redirect_stdout(io.StringIO()):
            optics.main([*self.arguments('standalone'), '--formats', 'npz'])
            cli.main(['waveder-optics', *self.arguments('composable'), '--formats', 'npz'])
        with np.load(self.fixture.root/'standalone/transitions.npz') as a, np.load(self.fixture.root/'composable/transitions.npz') as b:
            for key in a.files:
                np.testing.assert_equal(a[key], b[key])
        self.assertFalse((self.fixture.root/'composable/transitions.csv').exists())


if __name__ == '__main__':
    unittest.main()
