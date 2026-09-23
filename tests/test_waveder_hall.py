"""Standard WAVEDER occupied-bundle physics and same-run validation contract."""
import contextlib
import io
import json
from pathlib import Path
import struct
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import waveder_hall as wh
import vaspberry_kubo as cli


def write_waveder(path, c):
    ns, nk, _, nb, nd = c.shape
    def record(blob):
        return struct.pack('<i', len(blob))+blob+struct.pack('<i', len(blob))
    arrays = [struct.pack('<4i', nb, nd, nk, ns), struct.pack('<d', 0.),
              np.zeros((3, 3), dtype='<f8').tobytes(),
              np.asarray(c.transpose(3, 4, 1, 0, 2), dtype='<c8').tobytes(order='F')]
    path.write_bytes(b''.join(record(v) for v in arrays))


class WaveDerHallTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name); self.run = self.root/'run'; self.run.mkdir()
        self.q = np.array([[0., 0., 0.], [0., .5, 0.], [.5, 0., 0.], [.5, .5, 0.]])
        self.energies = np.tile([-1., -1., 2.], (1, 4, 1))
        self.occupations = np.tile([.99999, .99998, .00003], (1, 4, 1))
        rng = np.random.default_rng(551)
        self.c = (rng.normal(size=(1, 4, 3, 3, 3))+1j*rng.normal(size=(1, 4, 3, 3, 3))).astype(np.complex64)
        self.spinors = 2; self.nelect = 2
        self.fake_sources = {}
        self.options = dict(occupied=2, spin=1, spinor_components=2, spin_multiplicity=1,
            sampling={'kind': 'uniform_full_2d', 'mesh': [2, 2], 'plane_axes': [0, 1]},
            energy_reference='unchanged fixture VASP zero', mu_reference=.25)
        self.write_run()
        self.patcher = patch.object(wh, 'Wavecar', side_effect=self.fake_wavecar)
        self.patcher.start(); self.addCleanup(self.patcher.stop)

    def fake_wavecar(self, path, spinor_components, spin):
        e, o, q = self.fake_sources.get(Path(path).parent.resolve(), (self.energies, self.occupations, self.q))
        return SimpleNamespace(energies=e[spin-1].copy(), occupations=o[spin-1].copy(),
            kpoints=q.copy(), header=SimpleNamespace(ispin=len(e), lattice=np.eye(3),
            reciprocal=2*np.pi*np.eye(3)), coefficients=lambda k, b: np.ones((1, spinor_components, 1)))

    def write_run(self):
        ns, nk, nb = self.energies.shape
        (self.run/'WAVECAR').write_bytes(b'fixture WAVECAR, reader mocked')
        write_waveder(self.run/'WAVEDER', self.c)
        (self.run/'INCAR').write_text('SYSTEM = fixture\nLOPTICS = .TRUE.; LPEAD=.FALSE.\nLNABLA = .FALSE.\n'
                                    'ICHARG = 11; LCHARG = .FALSE.\n')
        out = [' vasp.5.4.4.18Apr17 complex', ' NSPIN fixture',
               '| So try LREAL= Auto in the INCAR file. |',
               '| reciprocal projection scheme (i.e. LREAL=.FALSE.) |',
               f' ISPIN = {ns}', f' k-points NKPTS = {nk} k-points in BZ NKDIM = {nk} number of bands NBANDS = {nb}',
               f' LNONCOLLINEAR = {"T" if self.spinors == 2 else "F"}', ' LSORBIT = F',
               ' LNABLA = F', ' LREAL = F', ' LHFCALC = F', ' METAGGA = F', ' LEPSILON = F', ' LVEL = F',
               ' ISYM = -1', ' NSW = 0', ' ICHARG = 11', ' LCHARG = F',
               f' NELECT = {self.nelect}', ' DEG_THRESHOLD=0.2000000E-02',
               ' direct lattice vectors                 reciprocal lattice vectors',
               '1.000000000 0.000000000 0.000000000 1 0 0',
               '0.000000000 1.000000000 0.000000000 0 1 0',
               '0.000000000 0.000000000 1.000000000 0 0 1',
               ' k-points in reciprocal lattice and weights:']
        out += [' '.join(f'{v:.8f}' for v in q)+f' {1/nk:.3f}' for q in self.q]
        factor = 2 if ns == 1 and self.spinors == 1 else 1
        for s in range(ns):
            if ns == 2:
                out.append(f' spin component {s+1}')
            for k, q in enumerate(self.q):
                out += [f' k-point {k+1} : '+''.join(f'{v:10.4f}' for v in q),
                        '  band No.  band energies     occupation ']
                out += [f'{b+1:7d} {e:10.4f} {o*factor:10.5f}'
                        for b, (e, o) in enumerate(zip(self.energies[s, k], self.occupations[s, k]))]
        out += ['aborting loop because EDIFF is reached',
                'frequency dependent IMAGINARY DIELECTRIC FUNCTION (independent particle) density-density',
                'OPTICS:  cpu time 1.0', 'General timing and accounting informations for this job:']
        (self.run/'OUTCAR').write_text('\n'.join(out)+'\n')

    def make_chunks(self, selections, prefix='chunk'):
        original = self.run, self.energies, self.occupations, self.q, self.c
        runs = []
        try:
            for i, selected in enumerate(selections):
                self.run = self.root/f'{prefix}-{i}'; self.run.mkdir()
                self.energies = original[1][:, selected].copy()
                self.occupations = original[2][:, selected].copy()
                self.q = original[3][selected].copy(); self.c = original[4][:, selected].copy()
                self.write_run()
                for name in ('CHGCAR', 'POTCAR', 'POSCAR'):
                    (self.run/name).write_bytes(('same retained fixture input '+name).encode())
                self.fake_sources[self.run.resolve()] = self.energies.copy(), self.occupations.copy(), self.q.copy()
                runs.append(self.run)
        finally:
            self.run, self.energies, self.occupations, self.q, self.c = original
        return runs

    def spectrum(self, mus=(1., -.5, 0.), **options):
        return wh.waveder_hall_spectrum(self.run, mus, **dict(self.options, **options))

    def cli_arguments(self, output='hall'):
        return ['--run-dir', str(self.run), '--occupied', '2', '--spinor-components', '2',
                '--spin-multiplicity', '1', '--mesh', '2', '2', '--energy-reference', 'fixture zero',
                '--mu-min', '-.2', '--mu-max', '.2', '--mu-num', '3', '--mu-reference', '.25',
                '--output-dir', str(self.root/output)]

    def test_lower_occupied_empty_trace_sign_factor_and_source_precision(self):
        rows, meta = self.spectrum()
        expected = []
        for k in range(4):
            value = 0.
            for v in range(2):
                x, y = complex(self.c[0, k, 0, 2, v]), complex(self.c[0, k, 1, 2, v])
                value += -2*(x.conjugate()*y).imag
            expected.append(value)
        sigma = -2*np.pi*np.mean(expected)
        self.assertEqual([r['mu_eV'] for r in rows if r['region'] == 'total'], [1., -.5, 0.])
        for row in rows:
            self.assertAlmostEqual(row['sigma_e2_over_h'], sigma, places=12)
            self.assertEqual(row['delta_sigma_e2_over_h'], 0.)
            self.assertEqual(row['electrons_per_cell'], 2.)
        self.assertEqual(meta['global_gap_eV'], 3.)
        self.assertFalse(meta['integer_rounding_applied'])
        self.assertIn('complex64', meta['source_metadata']['source_precision'])
        self.assertIn('cannot be independently authenticated', meta['source_metadata']['same_run_association']['limitation'])

    def test_covariant_trace_and_selection_rule_zeros_do_not_need_internal_channels(self):
        before = wh.occupied_curvature(self.c[0], 2)
        rng = np.random.default_rng(785)
        rotation = np.linalg.qr(rng.normal(size=(2, 2))+1j*rng.normal(size=(2, 2)))[0]
        changed = self.c[0].astype(complex)
        changed[:, :, 2:, :2] = changed[:, :, 2:, :2]@rotation
        # Occupied-occupied and unused upper blocks have no role in this trace.
        changed[:, :, :2, :] = 9999+123j
        np.testing.assert_allclose(wh.occupied_curvature(changed, 2), before, atol=1e-13)
        self.c[:] = 0; self.write_run()
        rows, _ = self.spectrum()
        self.assertTrue(all(row['sigma_e2_over_h'] == 0 for row in rows))

    def test_two_level_sign_and_all_three_components(self):
        c = np.zeros((1, 3, 2, 1), complex)
        c[0, :, 1, 0] = [.5j, -.5, .25+.3j]
        # Dx=sigma_x, Dy=sigma_y, energies=(-1,+1): occupied Ωxy=-1/2.
        omega = wh.occupied_curvature(c, 1)
        np.testing.assert_allclose(omega, [[.3, -.25, -.5]], atol=1e-15)

    def test_global_gap_cluster_boundary_filling_and_top_band_guards(self):
        for mus in [(-1.,), (2.,), (3.,)]:
            with self.subTest(mus=mus), self.assertRaisesRegex(ValueError, 'global insulating gap'):
                self.spectrum(mus)
        with self.assertRaisesRegex(ValueError, 'count disagrees'):
            self.spectrum(occupied=1)
        self.nelect = 3; self.write_run()
        with self.assertRaisesRegex(ValueError, 'NELECT'):
            self.spectrum()
        self.nelect = 2; self.energies[:] = [-.5, 0., .002]; self.write_run()
        with self.assertRaisesRegex(ValueError, '0.002 eV'):
            self.spectrum([.001], mu_reference=.001)

    def test_collinear_spin_filling_and_scalar_spin_factor(self):
        self.spinors = 1; self.nelect = 4; self.write_run()
        one, _ = self.spectrum(spinor_components=1, spin_multiplicity=1)
        two, _ = self.spectrum(spinor_components=1, spin_multiplicity=2)
        self.assertAlmostEqual(two[0]['sigma_e2_over_h'], 2*one[0]['sigma_e2_over_h'])
        self.energies = np.concatenate((self.energies, np.tile([-2., 1., 3.], (1, 4, 1))))
        self.occupations = np.concatenate((self.occupations, np.tile([1., 0., 0.], (1, 4, 1))))
        self.c = np.concatenate((self.c, self.c)); self.nelect = 3; self.write_run()
        rows, meta = self.spectrum([0.], spinor_components=1, spin=2, occupied=1)
        self.assertEqual(meta['source_metadata']['occupied_counts_by_spin'], [2, 1])
        self.assertEqual(rows[0]['electrons_per_cell'], 1.)
        with self.assertRaisesRegex(ValueError, 'multiplicity one'):
            self.spectrum([0.], spinor_components=1, spin_multiplicity=2)

    def test_settings_unsupported_or_inconsistent_are_rejected(self):
        original = (self.run/'OUTCAR').read_text()
        for old, new in [('vasp.5.4.4', 'vasp.6.5.0'), ('NSW = 0', 'NSW = 1'),
                         ('LNABLA = F', 'LNABLA = T'), ('LHFCALC = F', 'LHFCALC = T'),
                         ('ISYM = -1', 'ISYM = 0'), ('DEG_THRESHOLD=0.2000000E-02', 'DEG_THRESHOLD=1E-10'),
                         ('aborting loop because EDIFF is reached', 'not converged')]:
            with self.subTest(new=new):
                (self.run/'OUTCAR').write_text(original.replace(old, new))
                with self.assertRaises(ValueError):
                    self.spectrum()
        (self.run/'OUTCAR').write_text(original)
        (self.run/'INCAR').write_text('LOPTICS=T; LPEAD=T; LNABLA=F\n')
        with self.assertRaisesRegex(ValueError, 'LPEAD'):
            self.spectrum()
        (self.run/'INCAR').write_text('LOPTICS=T; LPEAD=F; LNABLA=F; NELECT=3\n')
        with self.assertRaisesRegex(ValueError, 'INCAR/OUTCAR'):
            self.spectrum()
        (self.run/'INCAR').write_text('LOPTICS=T; LPEAD=F; LNABLA=F; NBANDS=2\n')
        self.spectrum()  # Effective NBANDS may be increased by MPI distribution.
        (self.run/'INCAR').write_text('LOPTICS=T; LPEAD=F; LNABLA=F; NBANDS=4\n')
        with self.assertRaisesRegex(ValueError, 'NBANDS exceeds'):
            self.spectrum()

    def test_final_energy_k_lattice_occupation_and_dimension_mismatches_reject(self):
        original = (self.run/'OUTCAR').read_text()
        for old, new in [('-1.0000', '-1.0100'), ('0.99999', '0.89999'),
                         ('0.50000000', '0.51000000'), ('1.000000000', '1.100000000')]:
            with self.subTest(new=new):
                (self.run/'OUTCAR').write_text(original.replace(old, new))
                with self.assertRaisesRegex(ValueError, 'mismatch'):
                    self.spectrum()
        (self.run/'OUTCAR').write_text(original)
        write_waveder(self.run/'WAVEDER', self.c[:, :3])
        with self.assertRaisesRegex(ValueError, 'dimensions disagree'):
            self.spectrum()
        write_waveder(self.run/'WAVEDER', self.c[:, :, :, :, :1])
        with self.assertRaisesRegex(ValueError, 'coverage'):
            self.spectrum()

    def test_mesh_regions_orientation_and_nonfinite_checks(self):
        rows, _ = self.spectrum([0.], region_spec={'regions': [{'name': 'left', 'k_ids': [1, 2]},
                                                       {'name': 'right', 'k_ids': [3, 4]}]},
                                differences=[('contrast', 'left', 'right')])
        r = {row['region']: row['sigma_e2_over_h'] for row in rows}
        self.assertAlmostEqual(r['total'], r['left']+r['right'])
        self.assertAlmostEqual(r['contrast'], r['left']-r['right'])
        reverse, _ = self.spectrum([0.], sampling={'kind': 'uniform_full_2d', 'mesh': [2, 2], 'plane_axes': [1, 0]})
        self.assertAlmostEqual(reverse[0]['sigma_e2_over_h'], -r['total'])
        with self.assertRaises(ValueError):
            self.spectrum(sampling={'kind': 'uniform_full_2d', 'mesh': [3, 3], 'plane_axes': [0, 1]})
        with self.assertRaisesRegex(ValueError, 'finite'):
            self.spectrum([np.nan])

    def test_cli_all_formats_no_clobber_and_no_hand_authored_manifest(self):
        args = self.cli_arguments()
        with contextlib.redirect_stdout(io.StringIO()):
            wh.main(args)
        out = self.root/'hall'
        self.assertEqual(set(p.name for p in out.iterdir()),
                         {'conductivity.csv', 'conductivity.dat', 'conductivity.npz', 'conductivity.json'})
        with np.load(out/'conductivity.npz', allow_pickle=False) as data:
            np.testing.assert_array_equal(data['mu_eV'][data['region'] == 'total'], [-.2, 0., .2])
        meta = json.loads((out/'conductivity.json').read_text())
        self.assertEqual(set(meta['source_metadata']['source_sha256']), {'INCAR', 'OUTCAR', 'WAVEDER', 'WAVECAR'})
        before = (out/'conductivity.json').read_bytes()
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            wh.main(args)
        self.assertEqual((out/'conductivity.json').read_bytes(), before)

    def test_composable_cli_reaches_same_api_and_npz_results(self):
        with contextlib.redirect_stdout(io.StringIO()):
            wh.main([*self.cli_arguments('standalone'), '--formats', 'npz'])
        with patch.object(wh, 'waveder_hall_spectrum', wraps=wh.waveder_hall_spectrum) as api, \
             contextlib.redirect_stdout(io.StringIO()) as stdout:
            cli.main(['waveder-hall', *self.cli_arguments('composable'), '--formats', 'npz'])
        api.assert_called_once()
        self.assertEqual(json.loads(stdout.getvalue())['schema'], 'vaspberry.hall-spectrum')
        with np.load(self.root/'standalone/conductivity.npz', allow_pickle=False) as one, \
             np.load(self.root/'composable/conductivity.npz', allow_pickle=False) as two:
            self.assertEqual(one.files, two.files)
            for key in one.files:
                np.testing.assert_array_equal(one[key], two[key])
        self.assertFalse((self.root/'composable/conductivity.csv').exists())

    def test_composable_cli_rejects_unsupported_temperature_and_source_branch(self):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as exc:
            cli.main(['waveder-hall', *self.cli_arguments('temperature'), '--temperatures', '300'])
        self.assertEqual(exc.exception.code, 2)
        self.assertFalse((self.root/'temperature').exists())
        (self.run/'INCAR').write_text('LOPTICS=T; LPEAD=T; LNABLA=F\n')
        with contextlib.redirect_stderr(io.StringIO()) as stderr, self.assertRaises(SystemExit) as exc:
            cli.main(['waveder-hall', *self.cli_arguments('pead')])
        self.assertEqual(exc.exception.code, 2)
        self.assertIn('LPEAD', stderr.getvalue())
        self.assertFalse((self.root/'pead').exists())

    def test_genuine_fixed_charge_chunks_equal_one_full_run_and_cli_is_reusable(self):
        expected, _ = self.spectrum([-.2, 0., .2])
        runs = self.make_chunks([[0, 2], [1, 3]])
        rows, meta = wh.waveder_hall_spectrum(runs, [-.2, 0., .2], **self.options)
        np.testing.assert_allclose([r['sigma_e2_over_h'] for r in rows],
                                   [r['sigma_e2_over_h'] for r in expected], atol=1e-13)
        self.assertEqual(meta['source_metadata']['source_run_count'], 2)
        self.assertEqual(meta['source_metadata']['source_nkpoints'], 4)
        self.assertEqual(len(meta['source_metadata']['chunk_input_records']), 2)
        self.assertEqual(meta['source_metadata']['source_runs'][0]['producer_settings']['NKPTS'], 2)
        args = self.cli_arguments('chunk-hall')
        args[1:2] = [str(p) for p in runs]
        with contextlib.redirect_stdout(io.StringIO()):
            cli.main(['waveder-hall', *args, '--formats', 'npz'])
        saved = json.loads((self.root/'chunk-hall/conductivity.json').read_text())
        self.assertEqual(saved['provenance']['command']['run_dir'], [str(p) for p in runs])
        self.assertFalse((self.root/'chunk-hall/WAVEDER').exists())
        self.assertFalse((self.root/'chunk-hall/OUTCAR').exists())

    def test_chunks_reject_mixed_hamiltonians_or_nonretained_charge(self):
        for name in ('CHGCAR', 'POTCAR', 'POSCAR'):
            runs = self.make_chunks([[0, 2], [1, 3]], prefix=name)
            (runs[1]/name).write_bytes(b'different Hamiltonian input')
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, 'different CHGCAR/POTCAR/POSCAR'):
                wh.waveder_hall_spectrum(runs, [0.], **self.options)
        runs = self.make_chunks([[0, 2], [1, 3]], prefix='mode')
        p = runs[1]/'OUTCAR'; p.write_text(p.read_text().replace(' ICHARG = 11', ' ICHARG = 2'))
        with self.assertRaisesRegex(ValueError, 'effective ICHARG'):
            wh.waveder_hall_spectrum(runs, [0.], **self.options)
        runs = self.make_chunks([[0, 2], [1, 3]], prefix='parameters')
        p = runs[1]/'INCAR'; p.write_text(p.read_text()+'LDAUU = 5.0\n')
        with self.assertRaisesRegex(ValueError, 'identical INCAR'):
            wh.waveder_hall_spectrum(runs, [0.], **self.options)

    def test_chunks_reject_duplicate_incomplete_union_spin_or_occupation(self):
        for name, selected in [('duplicate', [[0, 1], [0, 1]]), ('incomplete', [[0, 1], [2]])]:
            runs = self.make_chunks(selected, prefix=name)
            with self.subTest(name=name), self.assertRaises(ValueError):
                wh.waveder_hall_spectrum(runs, [0.], **self.options)
        runs = self.make_chunks([[0, 2], [1, 3]], prefix='spin')
        p = runs[1]/'OUTCAR'; p.write_text(p.read_text().replace(' ISPIN = 1', ' ISPIN = 2'))
        with self.assertRaises(ValueError):
            wh.waveder_hall_spectrum(runs, [0.], **self.options)
        runs = self.make_chunks([[0, 2], [1, 3]], prefix='filling')
        with self.assertRaisesRegex(ValueError, 'fixed occupied count'):
            wh.waveder_hall_spectrum(runs, [-2.], **dict(self.options, mu_reference=-2.))
        with self.assertRaisesRegex(ValueError, 'distinct nonempty'):
            wh.waveder_hall_spectrum([runs[0], runs[0]], [0.], **self.options)

    def test_chunk_files_must_remain_unchanged_until_integration_finishes(self):
        runs = self.make_chunks([[0, 2], [1, 3]], prefix='mutating')
        original = wh.read_validated_run
        def changing_run(run, *args, **kwargs):
            data = original(run, *args, **kwargs)
            if Path(run) == runs[1]:
                path = runs[0]/'OUTCAR'
                path.write_text(path.read_text()+'source changed during integration\n')
            return data
        with patch.object(wh, 'read_validated_run', side_effect=changing_run), \
             self.assertRaisesRegex(ValueError, 'chunk source files changed'):
            wh.waveder_hall_spectrum(runs, [0.], **self.options)
