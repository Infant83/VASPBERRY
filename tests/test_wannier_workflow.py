"""Public ASCII import through full-connection Hall workflow and failure guards."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
sys.path.insert(0, str(ROOT/'tests'))
import wannier_workflow as ww
import wannier_bands as wb
from berry_data import CONDUCTANCE_QUANTUM_S
from exported_matrix_kubo import sha256
from plot_hall import read_table
from test_wannier_hall_physics import qwz, analytic_lower_curvature


def write_ascii(directory, model):
    directory.mkdir()
    hh, aa, poscar = (directory/name for name in ('model_HH_R.dat', 'model_AA_R.dat', 'POSCAR'))
    n = model.hamiltonian_eV.shape[-1]
    with hh.open('w') as h, aa.open('w') as a:
        h.write(f'Analytic test Hamiltonian\n{n}\n{len(model.irvec)}\n')
        a.write('Analytic test full position connection\n')
        for r, hr, ar in zip(model.irvec, model.hamiltonian_eV, model.position_A):
            for i in range(n):
                for j in range(n):
                    prefix = ''.join(f'{v:5d}' for v in (*r, i+1, j+1))
                    h.write(prefix+f'{hr[i,j].real:12.5E}{hr[i,j].imag:12.5E}\n')
                    a.write(prefix+''.join(f'{v:12.5E}' for c in range(3)
                            for v in (ar[c,i,j].real, ar[c,i,j].imag))+'\n')
    poscar.write_text('Analytic test lattice\n1.0\n'+'\n'.join(' '.join(map(str, row))
                         for row in model.lattice_A)+'\nX\n1\nDirect\n0 0 0\n')
    return hh, aa, poscar


class WannierWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)
        self.env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
                        MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', MPLBACKEND='Agg')

    def cache(self, name='operators', model=None):
        hh, aa, poscar = write_ascii(self.path/(name+'-input'), qwz() if model is None else model)
        out = self.path/name
        ww.import_command(SimpleNamespace(hh=hh, aa=aa, poscar=poscar, spinor_components=1,
            spin_multiplicity=1, energy_reference='analytic unchanged energy zero', output_dir=out))
        return out

    def args(self, operators, output='result', **overrides):
        values = dict(operators=operators, occupied=1, mesh=[8, 8], plane_axes=[0, 1],
            refine=1, refine_radius=None, refine_center=[], mu_min=-.2, mu_max=.2,
            mu_num=3, mu_reference=.1, gap_threshold=1e-8, workers=1, batch_size=4,
            time_limit=60., formats=['csv', 'dat', 'npz'], output_dir=self.path/output)
        values.update(overrides)
        return SimpleNamespace(**values)

    def test_actual_cli_import_then_hall_all_formats_and_metadata(self):
        hh, aa, poscar = write_ascii(self.path/'input', qwz())
        cache, output = self.path/'operators', self.path/'hall-result'
        command = [sys.executable, str(ROOT/'tools/vaspberry_kubo.py')]
        imported = subprocess.run(command+['wannier-import', '--hh', str(hh), '--aa', str(aa),
            '--poscar', str(poscar), '--spinor-components', '1', '--spin-multiplicity', '1',
            '--energy-reference', 'analytic zero', '--output-dir', str(cache)],
            env=self.env, capture_output=True, text=True, timeout=30)
        self.assertEqual(imported.returncode, 0, imported.stderr)
        self.assertEqual(json.loads(imported.stdout)['schema'], 'vaspberry.wannier-operators')
        response = subprocess.run(command+['wannier-hall', '--operators', str(cache), '--occupied', '1',
            '--mesh', '8', '8', '--mu-min', '-.2', '--mu-max', '.2', '--mu-num', '3',
            '--mu-reference', '.1', '--formats', 'csv', 'dat', 'npz', '--output-dir', str(output)],
            env=self.env, capture_output=True, text=True, timeout=30)
        self.assertEqual(response.returncode, 0, response.stderr)
        self.assertEqual(json.loads(response.stdout)['schema'], 'vaspberry.hall-spectrum')
        tables = [read_table(output/'hall'/('conductivity.'+fmt)) for fmt in ('csv', 'dat', 'npz')]
        for table in tables[1:]:
            for key in tables[0]: np.testing.assert_array_equal(table[key], tables[0][key])
        with np.load(output/'curvature.npz') as sample:
            expected = -2*np.pi*np.mean([analytic_lower_curvature(*(2*np.pi*q[:2]))
                                        for q in sample['kpoints_fractional']])
        np.testing.assert_allclose(tables[0]['sigma_e2_over_h'], expected, atol=2e-13)
        np.testing.assert_allclose(tables[0]['sigma_S'], expected*CONDUCTANCE_QUANTUM_S, atol=1e-17)
        meta = json.loads((output/'hall/conductivity.json').read_text())
        self.assertEqual(meta['method'], 'vaspberry_wannier_full_connection_T0')
        self.assertEqual(meta['term_order'], ['J0', 'J1', 'J2'])
        self.assertFalse(meta['integer_rounding_applied'])
        self.assertEqual(meta['source_metadata']['lattice_source']['sha256'], sha256(poscar))
        self.assertEqual(json.loads((output/'run.json').read_text())['status'], 'PASS')
        self.assertFalse(output.with_name(output.name+'.partial').exists())

    def test_npz_only_requested_mu_range_and_explicit_model_count(self):
        args = self.args(self.cache(), formats=['npz'], mu_min=-.1, mu_max=.15,
                         mu_num=4, mu_reference=.5)
        meta = ww.hall_command(args)
        files = {p.name for p in (args.output_dir/'hall').iterdir()}
        self.assertEqual(files, {'conductivity.npz', 'conductivity.json'})
        table = read_table(args.output_dir/'hall/conductivity.npz')
        np.testing.assert_array_equal(table['mu_eV'], np.linspace(-.1, .15, 4))
        np.testing.assert_array_equal(table['temperature_K'], np.zeros(4))
        np.testing.assert_array_equal(table['electrons_per_cell'], np.ones(4))
        np.testing.assert_array_equal(table['delta_sigma_e2_over_h'], np.zeros(4))
        self.assertIn('model', meta['electron_count_scope'])

    def test_refined_periodic_full_torus_weights_on_skew_cell(self):
        lattice = np.array([[2., 0., 0.], [1.7, .6, 0.], [0., 0., 4.]])
        center = np.array([.97, .02]); radius = .9; mesh = [6, 4]; refine = 3
        q, w, parent, reciprocal, refined = ww.quadrature(lattice, mesh, [0, 1], refine, radius, [center])
        # Independent finite periodic-image search for selected parent centers.
        offsets = np.array([(i, j) for i in range(-5, 6) for j in range(-5, 6)])
        selected = []
        for iy in range(mesh[1]):
            for ix in range(mesh[0]):
                delta = np.array([ix/mesh[0], iy/mesh[1]])-center
                distance = np.linalg.norm((delta[None]+offsets)@reciprocal[:2], axis=1).min()
                selected.append(distance <= radius)
        counts = np.bincount(parent)
        np.testing.assert_array_equal(counts, np.where(selected, refine**2, 1))
        np.testing.assert_allclose(np.bincount(parent, weights=w), 1/24, atol=2e-16)
        self.assertEqual(refined, sum(selected))
        self.assertAlmostEqual(w.sum(), 1., places=14)
        self.assertEqual(len(np.unique(q, axis=0)), len(q))
        self.assertTrue(np.all(q[:, 2] == 0))
        self.assertTrue(np.all((q >= 0) & (q < 1)))

    def test_nonorthogonal_inplane_geometry_and_orientation(self):
        lattice = np.array([[2., 0., 0.], [.7, 3., 0.], [0., 0., 5.]])
        cache = self.cache(model=qwz(lattice=lattice))
        forward = ww.hall_command(self.args(cache, 'forward', mesh=[24, 24]))
        reverse = ww.hall_command(self.args(cache, 'reverse', mesh=[24, 24], plane_axes=[1, 0]))
        a = read_table(self.path/'forward/hall/conductivity.npz')['sigma_e2_over_h']
        b = read_table(self.path/'reverse/hall/conductivity.npz')['sigma_e2_over_h']
        np.testing.assert_allclose(a, -b, atol=2e-14)
        np.testing.assert_allclose(a, -1., atol=2e-6)
        np.testing.assert_allclose(forward['normal'], -np.asarray(reverse['normal']), atol=1e-15)
        self.assertAlmostEqual(forward['bz_area_inv_A2'], (2*np.pi)**2/6, places=12)

    def test_serial_and_threaded_full_connection_are_identical(self):
        cache = self.cache(model=qwz(embedded=True))
        ww.hall_command(self.args(cache, 'serial', workers=1))
        ww.hall_command(self.args(cache, 'threaded', workers=2))
        with np.load(self.path/'serial/curvature.npz') as a, np.load(self.path/'threaded/curvature.npz') as b:
            for key in a.files: np.testing.assert_array_equal(a[key], b[key])
            self.assertGreater(abs(a['omega_terms_A2'][:, :2]).max(), .05)
        np.testing.assert_array_equal(read_table(self.path/'serial/hall/conductivity.npz')['sigma_e2_over_h'],
                                     read_table(self.path/'threaded/hall/conductivity.npz')['sigma_e2_over_h'])

    def test_invalid_options_reject_before_outputs_and_no_clobber(self):
        cache = self.cache()
        invalid = [{'mu_min': float('nan')}, {'mu_num': 0}, {'workers': 0}, {'batch_size': 0},
                   {'formats': ['npz', 'npz']}, {'refine': 2}, {'refine_center': [[0., 0.]]},
                   {'plane_axes': [0, 0]}, {'time_limit': 0}, {'memory_limit_mib': 1e-9}]
        for i, overrides in enumerate(invalid):
            args = self.args(cache, 'invalid'+str(i), **overrides)
            with self.subTest(overrides=overrides), self.assertRaises(ValueError): ww.hall_command(args)
            self.assertFalse(args.output_dir.exists())
            self.assertFalse(args.output_dir.with_name(args.output_dir.name+'.partial').exists())
        args = self.args(cache)
        args.output_dir.mkdir(); sentinel = args.output_dir/'keep'; sentinel.write_text('user data')
        with self.assertRaisesRegex(ValueError, 'exists'): ww.hall_command(args)
        self.assertEqual(sentinel.read_text(), 'user data')

    def test_outside_gap_and_gap_closure_preserve_failed_partial(self):
        for name, model, changes in [('outside', qwz(), {'mu_max': 5.}),
                                     ('closed', qwz(mass=-2.), {})]:
            args = self.args(self.cache(name+'-cache', model), name, **changes)
            with self.subTest(name=name), self.assertRaises(ValueError): ww.hall_command(args)
            partial = args.output_dir.with_name(args.output_dir.name+'.partial')
            record = json.loads((partial/'run.json').read_text())
            self.assertEqual(record['status'], 'FAILED')
            self.assertFalse(args.output_dir.exists()); self.assertFalse((partial/'hall').exists())
            self.assertTrue(record['error'])

    def test_timeout_preserves_failed_partial_and_blocks_reuse(self):
        args = self.args(self.cache(), time_limit=1e-30, batch_size=1)
        with self.assertRaisesRegex(ValueError, 'time'): ww.hall_command(args)
        partial = args.output_dir.with_name(args.output_dir.name+'.partial')
        record = json.loads((partial/'run.json').read_text())
        self.assertEqual(record['status'], 'FAILED')
        self.assertGreater(record['completed_points'], 0)
        self.assertLess(record['completed_points'], record['total_points'])
        self.assertFalse(args.output_dir.exists()); self.assertFalse((partial/'hall').exists())
        args.time_limit = 60.
        with self.assertRaisesRegex(ValueError, 'unfinished'): ww.hall_command(args)

    def test_source_mutation_during_compute_cannot_claim_pass(self):
        args = self.args(self.cache())
        original = ww.occupied_curvature; changed = False
        def mutate(*positional, **keywords):
            nonlocal changed
            if not changed:
                with (args.operators/'operators.npz').open('ab') as handle: handle.write(b'changed during run')
                changed = True
            return original(*positional, **keywords)
        with patch.object(ww, 'occupied_curvature', side_effect=mutate), self.assertRaisesRegex(ValueError, 'changed'):
            ww.hall_command(args)
        self.assertFalse(args.output_dir.exists())
        partial = args.output_dir.with_name(args.output_dir.name+'.partial')
        self.assertEqual(json.loads((partial/'run.json').read_text())['status'], 'FAILED')

    def test_t0_cli_rejects_finite_temperature_option(self):
        command = [sys.executable, str(ROOT/'tools/vaspberry_kubo.py'), 'wannier-hall',
            '--operators', 'unused', '--occupied', '1', '--mesh', '4', '4', '--mu-min', '0',
            '--mu-max', '0', '--mu-num', '1', '--mu-reference', '0', '--temperatures', '300',
            '--output-dir', str(self.path/'unsupported')]
        result = subprocess.run(command, env=self.env, capture_output=True, text=True, timeout=30)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('unrecognized arguments', result.stderr)
        self.assertFalse((self.path/'unsupported').exists())

    def test_band_path_cartesian_distance_and_unique_internal_junctions(self):
        lattice = np.array([[2., 0., 0.], [.7, 3., 0.], [0., 0., 5.]])
        reciprocal = 2*np.pi*np.linalg.inv(lattice).T
        vertices = np.array([[0., 0., 0.], [.5, 0., 0.], [.5, .5, 0.], [0., 0., 0.]])
        q, distance, ticks = wb.path_grid(vertices, 5, reciprocal)
        self.assertEqual(len(q), 13)
        np.testing.assert_array_equal(q[[0, 4, 8, 12]], vertices)
        self.assertTrue(np.all(np.diff(distance) > 0))
        np.testing.assert_allclose(distance[[0, 4, 8, 12]], ticks, atol=1e-15)
        expected = np.linalg.norm(np.diff(vertices, axis=0)@reciprocal, axis=1)
        np.testing.assert_allclose(np.diff(ticks), expected, atol=1e-15)
        for vertex in vertices[1:-1]:
            self.assertEqual(np.count_nonzero(np.all(q == vertex, axis=1)), 1)
        with self.assertRaisesRegex(ValueError, 'differ'):
            wb.path_grid(np.array([[0., 0., 0.], [0., 0., 0.]]), 5, reciprocal)

    def test_actual_band_cli_analytic_energies_and_csv_npz_parity(self):
        cache = self.cache(); output = self.path/'bands'
        command = [sys.executable, str(ROOT/'tools/vaspberry_kubo.py'), 'wannier-bands',
            '--operators', str(cache), '--vertices', '0', '0', '0', '.5', '0', '0', '.5', '.5', '0',
            '--labels', 'G', 'X', 'M', '--points-per-segment', '7', '--formats', 'csv', 'npz',
            '--output-dir', str(output)]
        result = subprocess.run(command, env=self.env, capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout)['schema'], 'vaspberry.wannier-bands')
        with (output/'bands.csv').open() as handle: rows = list(csv.DictReader(handle))
        with np.load(output/'bands.npz') as data:
            q = data['kpoints_fractional']; energies = data['energies_eV']
            kx, ky = (2*np.pi*q[:, :2]).T
            norm = np.sqrt(np.sin(kx)**2+np.sin(ky)**2+(-1.2+np.cos(kx)+np.cos(ky))**2)
            np.testing.assert_allclose(energies, np.column_stack((-norm, norm)), atol=1e-14)
            np.testing.assert_array_equal(np.array([float(r['energy_eV']) for r in rows]).reshape(-1, 2), energies)
            np.testing.assert_array_equal(np.array([float(r['distance_inv_A']) for r in rows])[::2], data['distance_inv_A'])
        meta = json.loads((output/'bands.json').read_text())
        self.assertEqual(meta['point_count'], 13)
        self.assertEqual(meta['labels'], ['G', 'X', 'M'])


if __name__ == '__main__':
    unittest.main()
