"""Independent analytic character attribution, same-run checks and CLI outputs."""
import contextlib
from dataclasses import replace
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import procar_character as pc
from kubo_pairs import import_native_pairs, write_pairs, pair_hall_spectrum
from procar_projection import parse_procar

spec = importlib.util.spec_from_file_location('procar_fixture', ROOT/'examples/features/procar-character/make_fixture.py')
fixture = importlib.util.module_from_spec(spec); spec.loader.exec_module(fixture)


class CharacterTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.inputs = self.root/'fixture'; fixture.create(self.inputs)
        self.args = pc.parser().parse_args(['project', '--procar', str(self.inputs/'PROCAR'),
            '--wavecar', str(self.inputs/'WAVECAR'), '--outcar', str(self.inputs/'OUTCAR'),
            '--groups', str(self.inputs/'groups.json'), '--axis', '1', '0', '0',
            '--output-dir', str(self.root/'character')])
        pc.project_command(self.args)
        self.character = pc.read_characters(self.args.output_dir)
        self.pairs = import_native_pairs(self.inputs/'PAIRS.csv', self.inputs/'WAVECAR', spin=1,
            spinor_components=2, spin_multiplicity=1, sampling={'kind': 'uniform_full_2d', 'mesh': [2, 2], 'plane_axes': [0, 1]},
            energy_reference='synthetic analytic fixture zero')

    def scan(self, bands=(1,), mus=(0.,), ts=(0.,), **kwargs):
        return pc.character_hall(self.character, self.pairs, list(bands), mus, ts, mu_reference=0., **kwargs)

    def test_project_rotation_orbitals_raw_joint_weights_and_roundtrip(self):
        d = self.character
        np.testing.assert_allclose(d.characters[:, :, 0], np.broadcast_to([.2, .2, .2, 0], (4, 2, 4)))
        np.testing.assert_allclose(d.characters[:, :, 1], np.broadcast_to([.6, -.6, 0, .6], (4, 2, 4)))
        np.testing.assert_array_equal(d.characters[:, :, 0], d.characters[:, :, 3])
        np.testing.assert_allclose(d.all_ions_cartesian[..., 0], .8)
        np.testing.assert_allclose(d.all_ions_cartesian[..., 1], -.4)
        self.assertTrue((self.args.output_dir/'projection_diagnostics.csv').is_file())
        self.assertEqual(d.metadata['all_ions_charge_range'], [.8, .8])
        self.assertLess(d.metadata['alignment']['max_energy_delta_eV'], 1e-12)

    def test_analytic_hall_sign_no_extra_half_and_residual(self):
        rows, meta, omega, gaps = self.scan()
        got = {(r['group'], r['component']): r['attribution_e2_over_h'] for r in rows}
        np.testing.assert_allclose(omega[:, 0, 2], 1.)
        np.testing.assert_allclose(gaps, 2.)
        self.assertAlmostEqual(got['lower', 'charge'], -.4*np.pi)
        self.assertAlmostEqual(got['lower', 'plus'], -.4*np.pi)
        self.assertAlmostEqual(got['lower', 'minus'], 0.)
        self.assertAlmostEqual(got['upper', 'minus'], -1.2*np.pi)
        self.assertAlmostEqual(got['$unweighted', 'charge'], -2*np.pi)
        self.assertAlmostEqual(got['$all_projected', 'charge'], -1.6*np.pi)
        self.assertAlmostEqual(got['$unprojected_residual', 'charge'], -.4*np.pi)
        self.assertIn('selected_band', meta['scope'])
        self.assertEqual(meta['virtual_band_window'], [1, 2])

    def test_pair_reversal_total_agrees_with_existing_kernel(self):
        rows, _, omega, _ = self.scan(bands=(1, 2), mus=(-.5, .2), ts=(0., 300.))
        np.testing.assert_allclose(omega[:, 0, 2], -omega[:, 1, 2])
        expected, _ = pair_hall_spectrum(self.pairs, [-.5, .2], [0., 300.], mu_reference=0.)
        baseline = [r for r in rows if r['group'] == '$unweighted']
        for r in baseline:
            ref = next(e for e in expected if e['mu_eV'] == r['mu_eV'] and e['temperature_K'] == r['temperature_K'])
            self.assertAlmostEqual(r['attribution_e2_over_h'], ref['sigma_e2_over_h'], places=13)
            self.assertAlmostEqual(r['delta_attribution_e2_over_h'], ref['delta_sigma_e2_over_h'], places=13)

    def test_finite_temperature_direct_fermi_and_region_closure(self):
        regions = {'regions': [{'name': 'first', 'k_ids': [1]}]}
        rows, _, _, _ = self.scan(mus=(-1., -.8), ts=(300.,), region_spec=regions)
        selected = [r for r in rows if r['group'] == 'lower' and r['component'] == 'charge' and r['mu_eV'] == -1.]
        values = {r['region']: r['attribution_e2_over_h'] for r in selected}
        self.assertAlmostEqual(values['total'], -.2*np.pi)  # f(E=mu,T)=1/2
        self.assertAlmostEqual(values['first'], values['total']/4)
        self.assertAlmostEqual(values['total'], values['first']+values['rest'])
        ref_f = 1/(1+np.exp(-1/(8.617333262145e-5*300)))
        total = next(r for r in selected if r['region'] == 'total')
        self.assertAlmostEqual(total['delta_attribution_e2_over_h'], -.4*np.pi*(.5-ref_f))

    def test_unresolved_selected_states_fail_even_equal_occupation(self):
        e = np.tile([-1., -1.], (4, 1))
        character = replace(self.character, energies_eV=e)
        pairs = replace(self.pairs, energies_eV=e)
        with self.assertRaisesRegex(ValueError, 'unresolved individual-state'):
            pc.character_hall(character, pairs, [1], [0.], [0.], mu_reference=0.)
        for bands, threshold in [([0], 1e-5), ([1, 1], 1e-5), ([3], 1e-5), ([1], 0.), ([1], np.nan)]:
            with self.subTest(bands=bands, threshold=threshold), self.assertRaises(ValueError):
                pc.selected_curvature(self.pairs, bands, threshold)

    def test_mismatched_hash_states_and_nonuniform_weights_fail(self):
        cases = [replace(self.character, metadata=dict(self.character.metadata, wavecar_sha256='f'*64)),
                 replace(self.character, energies_eV=self.character.energies_eV+.01),
                 replace(self.character, weights=np.array([.1, .2, .3, .4]))]
        for changed in cases:
            with self.subTest(changed=changed.metadata['wavecar_sha256']), self.assertRaises(ValueError):
                pc.character_hall(changed, self.pairs, [1], [0.], [0.], mu_reference=0.)

    def test_project_rejects_bad_actual_outcar_and_state_mismatch(self):
        original = (self.inputs/'OUTCAR').read_text()
        variants = [original.replace('LNONCOLLINEAR = T', 'LNONCOLLINEAR = F'),
                    original.replace('LORBIT = 11', 'LORBIT = 12'),
                    original.replace('NSW = 0', 'NSW = 1'),
                    original.replace('1 -1.0000 1.0000', '1 -1.1000 1.0000'),
                    original.replace('transformation matrix from SAXIS', 'missing transformation from SAXIS'),
                    original.replace('0.0000000 m_x 0.0000000 m_y 1.0000000 m_z', '0.0000000 m_x 0.0000000 m_y 2.0000000 m_z')]
        self.args.output_dir = self.root/'bad'
        for text in variants:
            (self.inputs/'OUTCAR').write_text(text)
            with self.assertRaises(ValueError): pc.project_command(self.args)
            self.assertFalse(self.args.output_dir.exists())
        (self.inputs/'OUTCAR').write_text(original)
        procar = self.inputs/'PROCAR'; procar.write_text(procar.read_text().replace('energy -1.00000000', 'energy -1.10000000'))
        with self.assertRaisesRegex(ValueError, 'states mismatch'): pc.project_command(self.args)

    def test_invalid_group_axis_and_orbital_selectors_fail(self):
        d = parse_procar(self.inputs/'PROCAR', noncollinear=True, spin_to_cartesian=np.eye(3))
        cases = [{}, {'groups': []}, {'groups': [{'name': 'a', 'ions': [0]}]},
                 {'groups': [{'name': 'a', 'ions': [1, 1]}]},
                 {'groups': [{'name': '$unweighted', 'ions': [1]}]},
                 {'groups': [{'name': 'a', 'ions': [1], 'orbitals': ['dxy']}]},
                 {'groups': [{'name': 'a', 'ions': [1], 'orbitals': ['s', 's']}]}]
        for groups in cases:
            with self.subTest(groups=groups), self.assertRaises(ValueError): pc.group_characters(d, groups, [0, 0, 1])
        with self.assertRaisesRegex(ValueError, 'unit vector'):
            pc.group_characters(d, {'groups': [{'name': 'a', 'ions': [1]}]}, [0, 0, 2])

    def test_saved_cache_rejects_reserved_names_and_nonunit_axis(self):
        for metadata in (dict(self.character.metadata, axis_cartesian=[1.000001, 0, 0]),
                         dict(self.character.metadata, group_names=['$unweighted', 'upper', 'all', 'lower_s'])):
            with self.subTest(metadata=metadata['axis_cartesian']), self.assertRaises(ValueError):
                pc.validate_characters(replace(self.character, metadata=metadata))

    def test_outcar_weight_rounding_does_not_relax_hall_fullmesh_check(self):
        procar = self.inputs/'PROCAR'
        procar.write_text(procar.read_text().replace('weight = 0.25000000', 'weight = 0.25030000'))
        self.args.output_dir = self.root/'rounded'
        pc.project_command(self.args)
        data = pc.read_characters(self.args.output_dir)
        self.assertAlmostEqual(data.metadata['alignment']['max_outcar_weight_delta'], .0003)
        with self.assertRaisesRegex(ValueError, 'same full uniform mesh'):
            pc.character_hall(data, self.pairs, [1], [0.], [0.], mu_reference=0.)

    def test_cli_saved_outputs_plots_checksums_and_no_clobber(self):
        write_pairs(self.root/'pairs', self.pairs)
        argv = ['hall', '--character-dir', str(self.root/'character'), '--pairs-dir', str(self.root/'pairs'),
                '--bands', '1', '--mu-min', '-1.2', '--mu-max', '.2', '--mu-num', '3',
                '--mu-reference', '0', '--temperatures', '0', '300', '--output-dir', str(self.root/'hall')]
        with contextlib.redirect_stdout(io.StringIO()): self.assertEqual(pc.main(argv), 0)
        with np.load(self.root/'hall/character_hall.npz', allow_pickle=False) as saved:
            self.assertEqual(len(saved['mu_eV']), 3*2*2*(4*4+3))  # total and automatic rest
            self.assertEqual(saved['group'].dtype.kind, 'U')
        with contextlib.redirect_stdout(io.StringIO()):
            pc.main(['plot', '--character-dir', str(self.root/'character'), '--hall-dir', str(self.root/'hall'),
                     '--group', 'lower', '--band', '1', '--temperature', '300', '--delta',
                     '--output-dir', str(self.root/'figures')])
        self.assertEqual(len(list((self.root/'figures').glob('*.png'))), 2)
        for path in (self.root/'figures').iterdir(): self.assertGreater(path.stat().st_size, 100)
        with self.assertRaisesRegex(ValueError, 'exists'): pc.project_command(self.args)
        path = self.root/'character/character.npz'; path.write_bytes(path.read_bytes()+b'corrupt')
        with self.assertRaisesRegex(ValueError, 'checksum'): pc.read_characters(self.root/'character')


if __name__ == '__main__':
    unittest.main()
