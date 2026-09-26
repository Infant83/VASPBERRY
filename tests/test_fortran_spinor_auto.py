"""Exercise automatic WAVECAR layout detection through compiled VASPBERRY."""
import csv
import io
import itertools
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def write_wavecar(path, components=2, channels=1, counts=None, rtag=45200,
                  kpoints=((0.0, 0.0, 0.0), (0.17, 0.09, 0.0))):
    """Independent full-complex fixture with two orthogonal plane-wave states."""
    recl, bands, cutoff = 1024, 2, 8.0
    lattice = [2 * math.pi if i == j else 0.0 for i in range(3) for j in range(3)]
    data = bytearray(recl * (2 + channels * len(kpoints) * (bands + 1)))
    struct.pack_into('3d', data, 0, recl, channels, rtag)
    struct.pack_into('12d', data, recl, len(kpoints), bands, cutoff, *lattice)
    full_counts = []
    for spin in range(channels):
        for ik, k in enumerate(kpoints):
            ng = sum(sum((g[i] + k[i]) ** 2 for i in range(3)) / 0.262465831 < cutoff
                     for g in itertools.product(range(-4, 5), repeat=3))
            ncoeff = components * ng if counts is None else counts[spin][ik](ng)
            full_counts.append(ng)
            rec = 2 + (spin * len(kpoints) + ik) * (bands + 1)
            struct.pack_into('10d', data, rec * recl, ncoeff, *k, -1., 0., 1., 1., 0., 0.)
            for band in range(bands):
                struct.pack_into('2f', data, (rec + 1 + band) * recl + band * 8, 1., 0.)
    path.write_bytes(data)
    return full_counts


def table(path):
    text = path.read_text()
    metadata = dict(line[2:].split('=', 1) for line in text.splitlines()
                    if line.startswith('# ') and '=' in line)
    rows = list(csv.DictReader(io.StringIO('\n'.join(
        line for line in text.splitlines() if not line.startswith('#')))))
    return metadata, rows


@unittest.skipUnless(shutil.which('gfortran'), 'gfortran required')
class NativeSpinorAutoTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix='vaspberry-spinor-')
        cls.work = Path(cls.temp.name)
        cls.serial = cls.work / 'vaspberry'
        cls.parallel = cls.work / 'vaspberry-mpi'
        flags = ['-cpp', '-O1', '-g', '-fcheck=all', '-ffixed-line-length-none',
                 '-fallow-argument-mismatch']
        cls.invoke(['gfortran', *flags, str(ROOT / 'vaspberry.f'), '-llapack', '-lblas',
                    '-o', str(cls.serial)], cls.work)
        cls.has_mpi = bool(shutil.which('mpifort') and shutil.which('mpiexec'))
        if cls.has_mpi:
            cls.invoke(['mpifort', *flags, '-DMPI_USE', str(ROOT / 'vaspberry.f'),
                        '-llapack', '-lblas', '-o', str(cls.parallel)], cls.work)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    @staticmethod
    def invoke(command, cwd, success=True):
        env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT='1', OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1')
        result = subprocess.run(command, cwd=cwd, capture_output=True, text=True, env=env, timeout=90)
        if success and result.returncode:
            raise AssertionError(f'{command}\n{result.stdout}\n{result.stderr}')
        return result

    def run_case(self, name, components=2, channels=1, flags=(), mpi=False,
                 success=True, **fixture):
        case = self.work / (self._testMethodName + '-' + name)
        case.mkdir()
        write_wavecar(case / 'WAVECAR', components=components, channels=channels, **fixture)
        command = [str(self.parallel if mpi else self.serial), '--task', 'kubo-pairs',
                   '--wavecar', 'WAVECAR', '--pairs-csv', 'PAIRS.csv', *flags]
        if mpi:
            command = ['mpiexec', '-n', '2', *command]
        result = self.invoke(command, case, success=success)
        if not success:
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((case / 'PAIRS.csv').exists())
        return case, result

    def test_auto_scalar_spinor_and_collinear_match_explicit_layout(self):
        for components, channels in ((1, 1), (2, 1), (1, 2)):
            with self.subTest(components=components, channels=channels):
                name = f'{components}-{channels}'
                auto, result = self.run_case(name + '-auto', components, channels)
                explicit, _ = self.run_case(name + '-explicit', components, channels,
                                            flags=['-s' if components == 1 else '--spinor', str(components)])
                self.assertEqual((auto / 'PAIRS.csv').read_bytes(), (explicit / 'PAIRS.csv').read_bytes())
                metadata, rows = table(auto / 'PAIRS.csv')
                self.assertEqual(metadata['spinor_components'], str(components))
                self.assertEqual(metadata['source_nspin'], str(channels))
                self.assertEqual(len(rows), 2 * channels)
                self.assertNotIn('LSORBIT', result.stdout)

    def test_explicit_auto_and_legacy_s_alias(self):
        baseline, _ = self.run_case('default')
        for name, flags in [('auto', ['--spinor', 'auto']), ('legacy-auto', ['-s', 'auto']),
                            ('legacy-explicit', ['-s', '2'])]:
            case, _ = self.run_case(name, flags=flags)
            self.assertEqual((baseline / 'PAIRS.csv').read_bytes(), (case / 'PAIRS.csv').read_bytes())

    def test_explicit_mismatch_fails_including_collinear_override(self):
        for components, channels in ((1, 1), (2, 1), (1, 2)):
            for option in ('--spinor', '-s'):
                _, failed = self.run_case(f'{components}-{channels}-{option}', components, channels,
                                          flags=[option, str(3 - components)], success=False)
                self.assertIn('disagrees with WAVECAR', failed.stderr)

    def test_all_k_and_spin_records_are_validated_before_export(self):
        cases = [
            ('noninteger-count', 1, [[lambda n: n + .25, lambda n: n]], 'noninteger WAVECAR coefficient count'),
            ('nonfinite-count', 1, [[lambda n: float('nan'), lambda n: n]], 'nonfinite WAVECAR coefficient layout'),
            ('zero-count', 1, [[lambda n: 0, lambda n: n]], 'invalid WAVECAR coefficient count'),
            ('invalid-count', 1, [[lambda n: n, lambda n: n + 1]], 'unsupported WAVECAR layout'),
            ('changed-layout', 1, [[lambda n: n, lambda n: 2 * n]], 'inconsistent WAVECAR spinor layout'),
            ('later-spin', 2, [[lambda n: n, lambda n: n], [lambda n: n, lambda n: n + 1]],
             'unsupported WAVECAR layout'),
            ('collinear-spinor', 2, [[lambda n: 2*n, lambda n: 2*n], [lambda n: 2*n, lambda n: 2*n]],
             'ISPIN=2 needs scalar coefficients'),
        ]
        for name, channels, counts, message in cases:
            with self.subTest(name=name):
                _, failed = self.run_case(name, channels=channels, counts=counts, success=False)
                self.assertIn(message, failed.stderr)

    def test_unsupported_precision_and_gamma_tags_are_not_guessed(self):
        for rtag in (45210, 53300, 53310, 12345):
            _, failed = self.run_case(str(rtag), rtag=rtag, success=False)
            self.assertIn('unsupported WAVECAR RTAG', failed.stderr)

    def test_mpi_matches_serial_and_mismatch_aborts_all_ranks(self):
        if not self.has_mpi:
            self.skipTest('MPI compiler and launcher required')
        for components, channels in ((1, 1), (2, 1), (1, 2)):
            name = f'{components}-{channels}'
            serial, _ = self.run_case(name + '-serial', components, channels)
            parallel, _ = self.run_case(name + '-mpi', components, channels, mpi=True)
            self.assertEqual((serial / 'PAIRS.csv').read_bytes(), (parallel / 'PAIRS.csv').read_bytes())
        _, failed = self.run_case('mismatch-mpi', flags=['-s', '1'], mpi=True, success=False)
        self.assertIn('disagrees with WAVECAR', failed.stderr)
        _, failed = self.run_case('later-k-mpi', mpi=True, success=False,
                                  counts=[[lambda n: 2*n, lambda n: n + 1]])
        self.assertIn('unsupported WAVECAR layout', failed.stderr)


if __name__ == '__main__':
    unittest.main()
