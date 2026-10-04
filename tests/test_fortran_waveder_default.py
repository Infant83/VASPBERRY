"""Compiled native standard-WAVEDER route and explicit approximation barrier.

The four-record optical fixture is generated independently of production
readers; its occupied-empty contraction is known without plane-wave velocities.
"""
import csv
import io
import itertools
import math
import os
from pathlib import Path
import runpy
import shutil
import struct
import subprocess
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def write_sources(path, *, components=2, channels=1, nd=2, smear=False,
                  q=((0., 0., 0.), (.5, 0., 0.)), empty_shift=0., gap=None):
    path.mkdir(parents=True, exist_ok=True)
    q = np.asarray(q, dtype=float)
    nk, nb, recl, cutoff = len(q), 4, 2048, 8.
    lattice = np.eye(3) * (2 * math.pi)
    energies = np.tile([-2., -2., 1.+empty_shift, 3.+empty_shift], (channels, nk, 1))
    if gap is not None:
        energies[..., 2] = -2.+gap
    occ = np.tile([.99, .98, .02, .01] if smear else [1., 1., 0., 0.], (channels, nk, 1))
    data = bytearray(recl * (2 + channels*nk*(nb+1)))
    struct.pack_into('<3d', data, 0, recl, channels, 45200.)
    struct.pack_into('<12d', data, recl, nk, nb, cutoff, *lattice.ravel())
    for s in range(channels):
        for k, point in enumerate(q):
            ng = sum(sum((g[i]+point[i])**2 for i in range(3)) / .262465831 < cutoff
                     for g in itertools.product(range(-3, 4), repeat=3))
            rec = 2 + (s*nk+k)*(nb+1)
            values = [components*ng, *point]
            for e, o in zip(energies[s, k], occ[s, k]):
                values.extend((e, 0., o))
            struct.pack_into('<'+str(len(values))+'d', data, rec*recl, *values)
            for band in range(nb):
                struct.pack_into('<2f', data, (rec+1+band)*recl + band*8, 1., 0.)
    (path/'WAVECAR').write_bytes(data)
    # Direct storage layout m,n,k,spin,direction, all columns remain rectangular.
    c = np.zeros((nb, nd, nk, channels, 3), dtype=np.complex64)
    for s in range(channels):
        for k in range(nk):
            for n in range(min(nd, 2)):
                for m in range(2, nb):
                    c[m, n, k, s] = [1+.25j*(n+1), .2*(m+1)+1j*(k+s+1), -.3j]
            # Occupied internal elements are intentionally huge and not Hermitian:
            # they must not enter the occupied trace or require square completion.
            if nd >= 2:
                c[0, 1, k, s] = [1e4, 2e4j, 0.]
    def record(payload):
        return struct.pack('<i', len(payload))+payload+struct.pack('<i', len(payload))
    der = record(struct.pack('<4i', nb, nd, nk, channels))
    der += record(struct.pack('<d', 0.)) + record(np.zeros((3, 3), '<f8').tobytes())
    der += record(c.tobytes(order='F'))
    (path/'WAVEDER').write_bytes(der)
    (path/'INCAR').write_text('LOPTICS=T; LPEAD=F; LNABLA=F\nLREAL=F; ISYM=-1; NSW=0\n')
    mult = 2 if channels == 1 and components == 1 else 1
    out = [' vasp.5.4.4.18Apr17 complex',
           f' ISPIN = {channels}',
           f' k-points NKPTS = {nk} k-points in BZ NKDIM = {nk} number of bands NBANDS = {nb}',
           f' LNONCOLLINEAR = {"T" if components == 2 else "F"}', ' LSORBIT = F',
           ' LNABLA = F', ' LREAL = F', ' LHFCALC = F', ' METAGGA = F',
           ' LEPSILON = F', ' LVEL = F', ' ISYM = -1', ' NSW = 0',
           f' NELECT = {channels*2*mult}', ' DEG_THRESHOLD = 0.002',
           ' direct lattice vectors                 reciprocal lattice vectors']
    out.extend(' '.join(f'{v:.12f}' for v in row)+' 0 0 0' for row in lattice)
    out.append(' k-points in reciprocal lattice and weights:')
    out.extend(' '.join(f'{v:.8f}' for v in point)+f' {1/nk:.8f}' for point in q)
    out.extend((' k-point  1 :   0.0000 0.0000 0.0000  plane waves: 17',
                ' k-point  2 :   0.0000-0.3333 0.0000  plane waves: 17'))
    for s in range(channels):
        if channels == 2:
            out.append(f' spin component {s+1}')
        for k, point in enumerate(q):
            out.extend((f' k-point {k+1} : '+' '.join(f'{v:.4f}' for v in point),
                        ' band No.  band energies     occupation'))
            out.extend(f'{i+1} {e:.4f} {o*mult:.5f}'
                       for i, (e, o) in enumerate(zip(energies[s, k], occ[s, k])))
    out.extend(('aborting loop because EDIFF is reached',
                'frequency dependent IMAGINARY DIELECTRIC FUNCTION',
                'OPTICS:  cpu time 1.0',
                'General timing and accounting informations for this job:'))
    (path/'OUTCAR').write_text('\n'.join(out)+'\n')
    expected = -2*np.imag((c[2:, :2, ..., 0].astype(complex).conj() *
                           c[2:, :2, ..., 1]).sum(axis=(0, 1))).T
    return expected


def table(path):
    text = path.read_text()
    metadata = dict(line[2:].split('=', 1) for line in text.splitlines()
                    if line.startswith('# ') and '=' in line)
    rows = list(csv.DictReader(io.StringIO('\n'.join(
        line for line in text.splitlines() if not line.startswith('#')))))
    return metadata, rows


@unittest.skipUnless(shutil.which('gfortran'), 'gfortran required')
class NativeWavederTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix='vaspberry-waveder-')
        cls.work = Path(cls.temp.name)
        cls.binary = cls.work/'vaspberry'
        cls.mpi = cls.work/'vaspberry-mpi'
        flags = ['-cpp', '-O0', '-g', '-fcheck=all', '-ffixed-line-length-none',
                 '-fallow-argument-mismatch']
        cls.invoke(['gfortran', *flags, str(ROOT/'vaspberry.f'), '-llapack', '-lblas',
                    '-o', str(cls.binary)], cls.work)
        cls.has_mpi = bool(shutil.which('mpifort') and shutil.which('mpiexec'))
        if cls.has_mpi:
            cls.invoke(['mpifort', *flags, '-DMPI_USE', str(ROOT/'vaspberry.f'),
                        '-llapack', '-lblas', '-o', str(cls.mpi)], cls.work)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    @staticmethod
    def invoke(command, cwd, success=True):
        env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT='1', OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1',
                   OMPI_MCA_rmaps_base_oversubscribe='1')
        result = subprocess.run(command, cwd=cwd, capture_output=True, text=True, env=env, timeout=90)
        if success and result.returncode:
            raise AssertionError(f'{command}\n{result.stdout}\n{result.stderr}')
        return result

    def case(self, name, **kwargs):
        path = self.work/(self._testMethodName+'-'+name)
        expected = write_sources(path, **kwargs)
        return path, expected

    def run_case(self, path, flags=(), *, success=True, mpi=False):
        command = [str(self.mpi if mpi else self.binary), '--task', 'kubo', *flags]
        if mpi:
            command = ['mpiexec', '-n', '2', *command]
        result = self.invoke(command, path, success=success)
        if not success:
            self.assertNotEqual(result.returncode, 0, result.stdout+result.stderr)
            self.assertFalse((path/'KUBO_WAVEDER.csv').exists())
        return result

    def test_rectangular_trace_matches_direct_contraction_and_not_momentum(self):
        path, expected = self.case('normal')
        self.run_case(path)
        metadata, rows = table(path/'KUBO_WAVEDER.csv')
        self.assertEqual(metadata['schema'], 'VASPBERRY_WAVEDER_KUBO_OCCUPIED_V1')
        self.assertEqual(metadata['source_ndbands'], '2')
        self.assertEqual(metadata['source_nbands'], '4')
        self.assertEqual(metadata['result_status'], 'PASS')
        self.assertEqual(metadata['integration'], 'NONE; source k-point trace only')
        np.testing.assert_allclose([float(r['omega_z_A2']) for r in rows], expected.ravel(), rtol=1e-14)
        self.assertFalse((path/'KUBO.csv').exists())
        self.assertFalse(list(path.glob('BERRYCURV*')))

    def test_no_extra_denominator_and_smeared_source_still_defines_t0_gap(self):
        values = []
        for shift in (0., 4.):
            path, expected = self.case(str(shift), empty_shift=shift, smear=True)
            self.run_case(path, ['--bands', '1:2'])
            _, rows = table(path/'KUBO_WAVEDER.csv')
            values.append([float(r['omega_z_A2']) for r in rows])
        self.assertEqual(*values)

    def test_scalar_spin_factor_and_integral_sign(self):
        path, expected = self.case('scalar', components=1)
        self.run_case(path, ['--mesh', '2,1'])
        metadata, _ = table(path/'KUBO_WAVEDER.csv')
        target_chern = 2*expected.sum()/2/(2*math.pi)
        self.assertAlmostEqual(float(metadata['total_chern']), target_chern, delta=1e-13)
        self.assertAlmostEqual(float(metadata['sheet_hall_e2_over_h']), -target_chern, delta=1e-13)
        self.assertEqual(metadata['physical_spin_multiplicity'], '2')

    def test_legacy_mesh_sizes_do_not_request_standard_integration(self):
        path, _ = self.case('legacy-sizes')
        self.run_case(path, ['--bands', '1:2', '-kx', '2', '-ky', '1'])
        metadata, rows = table(path/'KUBO_WAVEDER.csv')
        self.assertEqual(len(rows), 2)
        self.assertEqual(metadata['integration'], 'NONE; source k-point selected geometry only')
        self.assertNotIn('total_chern', metadata)
        integral, _ = self.case('integral-needs-modern-mesh')
        result = self.invoke([str(self.binary), '--task', 'kubo-integral',
                              '--bands', '1:2', '-kx', '2', '-ky', '1'], integral, success=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('requires explicit --mesh', result.stderr)
        self.assertFalse((integral/'KUBO_WAVEDER.csv').exists())

    def test_standard_ignored_controls_fail_before_reading_source_or_writing_output(self):
        path = self.work/self._testMethodName
        path.mkdir()
        for flag, value in (('--kpoint', '2'), ('--theta', '90'), ('-sigma', '0.1')):
            result = self.invoke([str(self.binary), '--task', 'kubo',
                                  '--bands', '1:2', flag, value], path, success=False)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('not used by standard WAVEDER Kubo', result.stderr)
            self.assertNotIn('Error opening', result.stdout + result.stderr)
            self.assertFalse(list(path.iterdir()))

    def test_two_spin_channels_have_complete_distinct_rows(self):
        path, expected = self.case('collinear', components=1, channels=2)
        self.run_case(path, ['--mesh', '2,1'])
        metadata, rows = table(path/'KUBO_WAVEDER.csv')
        self.assertEqual(metadata['expected_rows'], '4')
        self.assertEqual(metadata['physical_spin_multiplicity'], '1')
        np.testing.assert_allclose([float(r['omega_z_A2']) for r in rows], expected.ravel(), rtol=1e-14)

    def test_missing_or_corrupt_waveder_never_silently_falls_back(self):
        for kind in ('missing', 'marker', 'dimension', 'truncated', 'nan', 'trailing', 'coverage'):
            path, _ = self.case(kind, nd=1 if kind == 'coverage' else 2)
            blob = bytearray((path/'WAVEDER').read_bytes())
            if kind == 'missing':
                (path/'WAVEDER').unlink()
            else:
                if kind == 'marker': struct.pack_into('<i', blob, 0, -16)
                if kind == 'dimension': struct.pack_into('<i', blob, 12, 999)
                if kind == 'truncated': del blob[-8:]
                if kind == 'trailing': blob.extend(b'junk')
                if kind == 'nan': struct.pack_into('<f', blob, 124, float('nan'))
                (path/'WAVEDER').write_bytes(blob)
            result = self.run_case(path, success=False)
            self.assertIn('WAVEDER', result.stderr)
            self.assertNotIn('WARNING:', result.stderr)
            self.assertFalse((path/'KUBO.csv').exists())

    def test_scope_metadata_and_mesh_guards(self):
        changes = {
            'branch': lambda p: (p/'INCAR').write_text('LOPTICS=T; LPEAD=T; LNABLA=F\n'),
            'effective': lambda p: self.replace_outcar(p, ' LNABLA = F', ' LNABLA = T'),
            'energy': lambda p: self.replace_outcar(p, '1 -2.0000 1.00000', '1 -2.1000 1.00000'),
            'version': lambda p: self.replace_outcar(p, 'vasp.5.4.4.', 'vasp.6.4.3.'),
            'incomplete': lambda p: self.replace_outcar(p, 'General timing and accounting informations for this job:', ''),
        }
        for name, change in changes.items():
            path, _ = self.case(name);change(path)
            self.run_case(path, success=False)
        for name, kwargs, flags in (
            ('cluster', {'gap': .001}, []),
            ('mesh', {'q': ((0., 0., 0.), (.2, 0., 0.))}, ['--mesh', '2,1']),
            ('partial-bundle', {}, ['--bands', '2']),
            ('per-band', {}, ['--bands', '1:2', '--per-band', '1']),
            ('ignored-prefix', {}, ['--output', 'old-prefix']),
        ):
            path, _ = self.case(name, **kwargs)
            self.run_case(path, flags, success=False)

    @staticmethod
    def replace_outcar(path, before, after):
        out = path/'OUTCAR';out.write_text(out.read_text().replace(before, after))

    def test_explicit_wavecar_is_warned_even_when_waveder_exists(self):
        path, _ = self.case('opt-in')
        result = self.run_case(path, ['--kubo-source', 'wavecar', '--bands', '1:2'])
        self.assertIn('WARNING:', result.stderr)
        self.assertIn('bare-momentum', result.stderr)
        self.assertTrue((path/'KUBO.csv').exists())
        self.assertFalse((path/'KUBO_WAVEDER.csv').exists())
        meta, _ = table(path/'KUBO.csv')
        self.assertIn('BARE_MOMENTUM', meta['schema'])

    def test_operator_comparison_example_keeps_its_canonical_pair_source(self):
        path, _ = self.case('operator-comparison')
        example = runpy.run_path(str(ROOT/'examples/features/kubo-hall/operator-comparison/run.py'))
        command = example['native_export_command'](self.binary, path/'WAVECAR')
        result = self.invoke(command, path)
        self.assertIn('WARNING:', result.stderr)
        self.assertTrue((path/'PAIRS.csv').is_file())
        meta, _ = table(path/'PAIRS.csv')
        self.assertEqual(meta['schema'], 'VASPBERRY_BARE_MOMENTUM_KUBO_PAIRS_V1')
        self.assertEqual(meta['result_status'], 'PASS')
        self.assertFalse((path/'KUBO_WAVEDER.csv').exists())

    def test_input_directory_explicit_paths_legacy_and_integral_dispatch(self):
        path, _ = self.case('input-dir')
        run = path/'source';run.mkdir()
        for name in ('WAVECAR', 'WAVEDER', 'INCAR', 'OUTCAR'):
            (path/name).rename(run/name)
        # A WAVECAR override alone never redirects the missing cwd auxiliaries.
        failed = self.run_case(path, ['--wavecar', 'source/WAVECAR'], success=False)
        self.assertIn('# Input WAVECAR: '+str((run/'WAVECAR').resolve()), failed.stdout)
        self.assertIn('# Input WAVEDER: '+str((path/'WAVEDER').resolve()), failed.stdout)
        self.run_case(path, ['--input-dir', 'source', '--curvature-csv', 'answer.csv'])
        self.assertTrue((path/'answer.csv').exists())
        self.assertFalse((run/'answer.csv').exists())
        self.invoke([str(self.binary), '-kubo', '2', '--input-dir', 'source',
                     '-kubo_csv', 'legacy.csv'], path)
        self.assertTrue((path/'legacy.csv').exists())
        result = self.invoke([str(self.binary), '--task', 'kubo-integral', '--input-dir', 'source'],
                             path, success=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('explicit --mesh', result.stderr)
        self.invoke([str(self.binary), '--task', 'kubo-integral', '--input-dir', 'source',
                     '--mesh', '2,1', '--curvature-csv', 'integral.csv'], path)
        self.assertIn('total_chern', table(path/'integral.csv')[0])

    def test_wavecar_override_keeps_cwd_auxiliaries_and_reports_before_open(self):
        path, expected = self.case('cwd-aux')
        other = path/'other';other.mkdir()
        (path/'WAVECAR').rename(other/'different.wave')
        for name in ('WAVEDER', 'INCAR', 'OUTCAR'):
            (other/name).write_text('invalid sibling must never be selected')
        result = self.run_case(path, ['--wavecar', 'other/different.wave'])
        self.assertIn('# Input WAVECAR: '+str((other/'different.wave').resolve()), result.stdout)
        for name in ('WAVEDER', 'INCAR', 'OUTCAR'):
            self.assertIn('# Input '+name+': '+str((path/name).resolve()), result.stdout)
        _, rows = table(path/'KUBO_WAVEDER.csv')
        np.testing.assert_allclose([float(r['omega_z_A2']) for r in rows], expected.ravel(), rtol=1e-14)
        (path/'KUBO_WAVEDER.csv').unlink()
        failed = self.run_case(path, ['--wavecar', 'missing.wave'], success=False)
        self.assertIn('# Input WAVECAR: '+str((path/'missing.wave').resolve()), failed.stdout)
        for name in ('WAVEDER', 'INCAR', 'OUTCAR'):
            self.assertIn('# Input '+name+': '+str((path/name).resolve()), failed.stdout)
        self.assertLess(failed.stdout.index('# Input OUTCAR:'), failed.stdout.index('# File reading...'))

    def test_each_file_override_is_cwd_relative_and_order_independent(self):
        path, expected = self.case('overrides')
        default = path/'input space';default.mkdir()
        overrides = path/'overrides';overrides.mkdir()
        names = ('WAVECAR', 'WAVEDER', 'INCAR', 'OUTCAR')
        for name in names:
            (path/name).rename(default/name)
            shutil.copy2(default/name, overrides/('custom.'+name))
        for index, name in enumerate(names):
            original = (default/name).read_bytes()
            (default/name).write_text('invalid default; explicit override required')
            for order in (0, 1):
                output = f'one-{index}-{order}.csv'
                flags = ['--'+name.lower(), 'overrides/custom.'+name]
                directory = ['--input-dir', 'input space']
                flags = [*flags, *directory] if order else [*directory, *flags]
                result = self.run_case(path, [*flags, '--curvature-csv', output])
                self.assertIn('# Input '+name+': '+str((overrides/('custom.'+name)).resolve()), result.stdout)
                _, rows = table(path/output)
                np.testing.assert_allclose([float(r['omega_z_A2']) for r in rows], expected.ravel(), rtol=1e-14)
            (default/name).write_bytes(original)
        # Every file may be overridden independently while the directory remains validated.
        flags = []
        for name in names:
            flags.extend(('--'+name.lower(), str(overrides/('custom.'+name))))
        self.run_case(path, ['--input-dir', 'input space', *flags])
        meta, _ = table(path/'KUBO_WAVEDER.csv')
        for name in names:
            self.assertEqual(meta[name.lower()], str(overrides/('custom.'+name)))

    def test_directory_and_path_errors_fail_even_with_all_file_overrides(self):
        path, _ = self.case('invalid-dir')
        files = []
        for name in ('WAVECAR', 'WAVEDER', 'INCAR', 'OUTCAR'):
            files.extend(('--'+name.lower(), str(path/name)))
        for directory in ('missing', 'WAVECAR'):
            failed = self.run_case(path, ['--input-dir', directory, *files], success=False)
            self.assertIn('must be an existing directory', failed.stderr)
        for option in ('--input-dir', '--wavecar', '--waveder', '--incar', '--outcar', '-f'):
            for value in ('', 'x'*257):
                self.run_case(path, [option, value], success=False)
        # A valid directory may fit while its default filename would overflow.
        longdir = path/('a'*(251-len(str(path))-1))
        longdir.mkdir()
        self.assertEqual(len(str(longdir)), 251)
        failed = self.run_case(path, ['--input-dir', str(longdir)], success=False)
        self.assertIn('exceeds 256 characters', failed.stderr)
        self.run_case(path, ['--input-dir', str(longdir), *files])

    def test_input_directory_serial_mpi_output_and_failure_agree(self):
        if not self.has_mpi:
            self.skipTest('MPI tools required')
        path, _ = self.case('mpi-dir')
        source = path/'source';source.mkdir()
        for name in ('WAVECAR', 'WAVEDER', 'INCAR', 'OUTCAR'):
            (path/name).rename(source/name)
        flags = ['--input-dir', 'source']
        serial = self.run_case(path, flags)
        before = (path/'KUBO_WAVEDER.csv').read_bytes()
        (path/'KUBO_WAVEDER.csv').unlink()
        parallel = self.run_case(path, flags, mpi=True)
        self.assertEqual(before, (path/'KUBO_WAVEDER.csv').read_bytes())
        self.assertEqual(parallel.stdout.count('# Input WAVECAR:'), 1)
        (path/'KUBO_WAVEDER.csv').unlink()
        failed = self.run_case(path, ['--input-dir', 'not-a-directory'], mpi=True, success=False)
        self.assertIn('must be an existing directory', failed.stderr)

    def test_duplicate_or_partial_latest_outcar_tables_are_rejected(self):
        for duplicate in (False, True):
            path, _ = self.case(str(duplicate))
            out = path/'OUTCAR'
            text = out.read_text()
            start = text.index(' k-point 1 :')
            end = text.index(' k-point 2 :', start)
            k1 = text[start:end]
            footer = text.index('aborting loop because EDIFF is reached')
            text = text[:footer] + k1*(2 if duplicate else 1) + text[footer:]
            out.write_text(text)
            self.run_case(path, success=False)

    def test_source_option_mismatches_are_rejected(self):
        path, _ = self.case('source-mix')
        for flags in (['--kubo-source', 'wavecar', '--outcar', 'DOES_NOT_EXIST'],
                      ['--kubo-source', 'wavecar', '--waveder', 'WAVEDER']):
            self.run_case(path, flags, success=False)

    def test_spin_and_pair_kubo_need_explicit_approximation_choice(self):
        path, _ = self.case('other-kubo')
        for task, flags in (('spin-kubo', ['--bands', '1:2']),
                            ('kubo-pairs', ['--pairs-csv', 'PAIRS.csv'])):
            result = self.invoke([str(self.binary), '--task', task, *flags], path, success=False)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('--kubo-source wavecar', result.stderr)
        self.assertFalse((path/'PAIRS.csv').exists())

    def test_existing_output_is_preserved_and_mpi_matches_serial(self):
        path, _ = self.case('protected')
        self.run_case(path)
        before = (path/'KUBO_WAVEDER.csv').read_bytes()
        result = self.invoke([str(self.binary), '--task', 'kubo'], path, success=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual((path/'KUBO_WAVEDER.csv').read_bytes(), before)
        if self.has_mpi:
            (path/'KUBO_WAVEDER.csv').unlink()
            self.run_case(path, mpi=True)
            self.assertEqual((path/'KUBO_WAVEDER.csv').read_bytes(), before)
            (path/'KUBO_WAVEDER.csv').unlink()
            (path/'WAVEDER').unlink()
            self.run_case(path, mpi=True, success=False)


if __name__ == '__main__':
    unittest.main()
