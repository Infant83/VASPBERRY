"""Effective full-operator import: orientation, residual sums and rejection guards."""
from dataclasses import replace
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import wannier_operators as wo
from exported_matrix_kubo import sha256


def field(value):
    # Explicit decimals prevent the Fortran F12.6 implied-decimal rule.
    return f'{value:12.5E}'


def row(r, i, j, values):
    return ''.join(f'{v:5d}' for v in (*r, i, j))+''.join(field(v) for v in values)+'\n'


def fixture(path, order=(0, 1, 2)):
    r = np.array([[0, 0, 0], [1, -2, 0], [-1, 2, 0]], dtype=np.int64)[list(order)]
    h0 = np.array([[1., .2+.3j], [.2-.3j, 2.]])
    h1 = np.array([[.11+.12j, .13+.14j], [.15+.16j, .17+.18j]])
    h = np.array([h0, h1, h1.conj().T], dtype=np.complex128)[list(order)]
    a0 = np.array([[.4, .5+.6j], [.5-.6j, .7]])
    a1 = np.array([[.21+.22j, .23+.24j], [.25+.26j, .27+.28j]])
    a = np.array([[a0*(c+1) for c in range(3)], [a1*(c+1) for c in range(3)],
                  [a1.conj().T*(c+1) for c in range(3)]], dtype=np.complex128)[list(order)]
    hh, aa = path/'test_HH_R.dat', path/'test_AA_R.dat'
    hs, az = 'effective H\n2\n3\n', 'effective A\n'
    for k, v in enumerate(r):
        for i in range(2):
            for j in range(2):
                hs += row(v, i+1, j+1, [h[k, i, j].real, h[k, i, j].imag])
                az += row(v, i+1, j+1, np.column_stack((a[k, :, i, j].real, a[k, :, i, j].imag)).ravel())
    hh.write_text(hs); aa.write_text(az)
    return hh, aa, r, h, a


class EffectiveOperatorsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)
        self.hh, self.aa, self.r, self.h, self.a = fixture(self.path)

    def load(self, **kwargs):
        options = dict(spinor_components=2, spin_multiplicity=1, energy_reference='fixture energy zero')
        options.update(kwargs)
        return wo.load_effective_operators(self.hh, self.aa, np.diag([2., 3., 5.]), **options)

    def test_complex_orientation_and_cartesian_order(self):
        data = self.load()
        np.testing.assert_array_equal(data.irvec, self.r)
        np.testing.assert_allclose(data.hamiltonian_eV, self.h, atol=1e-16)
        np.testing.assert_allclose(data.position_A, self.a, atol=1e-15)
        # Independently Fourier-sum an off-diagonal element: row/column reversal
        # or conjugating the file changes this deliberately complex quantity.
        q = np.array([.17, -.13, .09]); phase = np.exp(2j*np.pi*(self.r@q))
        actual = np.einsum('r,rij->ij', phase, data.hamiltonian_eV)
        expected = self.h[0]+phase[1]*self.h[1]+phase[2]*self.h[2]
        np.testing.assert_allclose(actual, expected, atol=1e-15)
        self.assertEqual(data.metadata['position_components'], ['x', 'y', 'z'])
        self.assertEqual(data.metadata['source_files'], {'HH_R': self.hh.name, 'AA_R': self.aa.name})

    def test_nonlexical_common_r_order_is_preserved(self):
        self.hh, self.aa, r, h, a = fixture(self.path, (1, 0, 2))
        data = self.load()
        np.testing.assert_array_equal(data.irvec, r)
        np.testing.assert_array_equal(data.hamiltonian_eV, h)
        np.testing.assert_allclose(data.position_A, a, atol=1e-15)

    def test_repeated_residual_entries_add_across_parse_chunks(self):
        lines = self.hh.read_text().splitlines(keepends=True)
        lines.insert(4, row([0, 0, 0], 1, 1, [2e-11, 0.]))
        self.hh.write_text(''.join(lines))
        from unittest.mock import patch
        with patch.object(wo, 'PARSE_ROWS', 1):
            data = self.load()
        self.assertEqual(data.hamiltonian_eV[0, 0, 0], 1.+2e-11)
        self.assertEqual(data.metadata['source_rows']['HH_R'], 13)

    def test_fortran_fixed_fields_and_implied_decimal(self):
        # This fixed-width form has no whitespace between adjacent real fields.
        self.assertEqual(wo._real('123e2'), .0123)
        self.assertAlmostEqual(wo._real('123'), .000123)
        self.assertEqual(wo._real('1.0D+2'), 100.)
        self.assertEqual(wo._real('1.0-2'), .01)
        line = row([0, 0, 0], 1, 1, [1., 0.])
        self.assertEqual(len(line.rstrip('\n')), 49)
        self.assertTrue(np.isfinite(self.load().hamiltonian_eV).all())

    def test_sparse_zeros_are_supported_without_discarding_nonzero(self):
        self.hh.write_text('sparse H\n2\n1\n'+row([0, 0, 0], 1, 1, [1., 0.]))
        self.aa.write_text('sparse A\n'+row([0, 0, 0], 1, 1, [1e-18, 0., 0., 0., 0., 0.]))
        data = self.load()
        self.assertEqual(np.count_nonzero(data.hamiltonian_eV), 1)
        self.assertEqual(data.position_A[0, 0, 0, 0], 1e-18)

    def test_missing_or_reordered_blocks_rejected(self):
        original = self.aa.read_text().splitlines(keepends=True)
        for bad in (original[:1]+original[1:5]+original[9:]+original[5:9], original[:-4]):
            with self.subTest(rows=len(bad)):
                self.aa.write_text(''.join(bad))
                with self.assertRaisesRegex(ValueError, 'R block'):
                    self.load()
        self.aa.write_text(''.join(original))
        lines = self.hh.read_text().splitlines(keepends=True)
        self.hh.write_text(''.join(lines+lines[3:4]))
        with self.assertRaisesRegex(ValueError, 'noncontiguous'):
            self.load()

    def test_opposite_r_and_origin_required(self):
        self.hh.write_text('H\n2\n1\n'+row([1, 0, 0], 1, 1, [1., 0.]))
        self.aa.write_text('A\n'+row([1, 0, 0], 1, 1, [0.]*6))
        with self.assertRaisesRegex(ValueError, 'origin'):
            self.load()
        self.hh.write_text('H\n2\n2\n'+row([0, 0, 0], 1, 1, [1., 0.])+row([1, 0, 0], 1, 1, [1., 0.]))
        self.aa.write_text('A\n'+row([0, 0, 0], 1, 1, [0.]*6)+row([1, 0, 0], 1, 1, [0.]*6))
        with self.assertRaisesRegex(ValueError, 'opposite'):
            self.load()

    def test_bad_indices_dimensions_and_rows_rejected(self):
        original = self.hh.read_text(); lines = original.splitlines(keepends=True)
        badlines = [lines[3][:-3]+'\n', lines[3].rstrip('\n')+' x\n', '\n',
                    row([0, 0, 0], 0, 1, [1., 0.]), row([0, 0, 0], 3, 1, [1., 0.]),
                    '  0.5'+lines[3][5:], row([0, 0, 0], 1, 1, [float('nan'), 0.]),
                    row([0, 0, 0], 1, 1, [float('inf'), 0.])]
        for bad in badlines:
            with self.subTest(bad=bad):
                self.hh.write_text(''.join(lines[:3]+[bad]+lines[4:]))
                with self.assertRaises(ValueError):
                    self.load()
        for dimension in ('2.0', '0', '-1', '99999999'):
            with self.subTest(dimension=dimension):
                self.hh.write_text(original.replace('\n2\n', '\n'+dimension+'\n', 1))
                with self.assertRaises(ValueError):
                    self.load()

    def test_both_operators_require_hermiticity_without_repair(self):
        for path, components in ((self.hh, 1), (self.aa, 3)):
            original = path.read_text(); lines = original.splitlines(keepends=True)
            index = 3 if components == 1 else 1
            fields = [1., .1]+[0.]*(2*components-2)
            lines[index] = row([0, 0, 0], 1, 1, fields)
            path.write_text(''.join(lines))
            with self.assertRaisesRegex(ValueError, 'Hermiticity'):
                self.load()
            path.write_text(original)
        lines = self.hh.read_text().splitlines(keepends=True)
        lines[3] = row([0, 0, 0], 1, 1, [1., 1e-12])
        self.hh.write_text(''.join(lines))
        self.assertEqual(self.load().hamiltonian_eV[0, 0, 0].imag, 1e-12)

    def test_spin_lattice_and_energy_declarations(self):
        for options in ({'spin_multiplicity': 2}, {'spinor_components': True},
                        {'spin_multiplicity': 1.0}, {'energy_reference': ''}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                self.load(**options)
        self.assertEqual(self.load(spinor_components=1, spin_multiplicity=2).metadata['spin_multiplicity'], 2)
        for lattice in (np.zeros((3, 3)), np.diag([1., 1., -1.]), np.full((3, 3), np.nan)):
            with self.assertRaises(ValueError):
                wo.load_effective_operators(self.hh, self.aa, lattice, spinor_components=1,
                                            spin_multiplicity=1, energy_reference='zero')

    def test_cache_roundtrip_checksums_and_no_overwrite(self):
        data = self.load(); out = self.path/'cache'
        metadata = wo.write_operators(out, data)
        restored = wo.read_operators(out)
        self.assertEqual(metadata, restored.metadata)
        for name in wo.ARRAYS:
            np.testing.assert_array_equal(getattr(restored, name), getattr(data, name))
        with self.assertRaisesRegex(ValueError, 'exists'):
            wo.write_operators(out, data)
        with (out/'operators.npz').open('ab') as handle:
            handle.write(b'corrupted')
        with self.assertRaisesRegex(ValueError, 'checksum'):
            wo.read_operators(out)

    def test_cache_bad_metadata_and_payload_cannot_bypass_validation(self):
        data = self.load(); out = self.path/'cache'
        wo.write_operators(out, data)
        meta = json.loads((out/'operators.json').read_text()); meta['fourier_phase'] = 'opposite'
        (out/'operators.json').write_text(json.dumps(meta))
        with self.assertRaisesRegex(ValueError, 'conventions'):
            wo.read_operators(out)
        meta['fourier_phase'] = wo.PHASE
        np.savez_compressed(out/'operators.npz', **{key: getattr(data, key).astype(float)
             if key == 'irvec' else getattr(data, key) for key in wo.ARRAYS})
        meta['data_npz_sha256'] = sha256(out/'operators.npz')
        (out/'operators.json').write_text(json.dumps(meta))
        with self.assertRaisesRegex(ValueError, 'integer R'):
            wo.read_operators(out)
        invalid = replace(data, hamiltonian_eV=data.hamiltonian_eV*1j)
        with self.assertRaisesRegex(ValueError, 'Hermiticity'):
            wo.write_operators(self.path/'invalid', invalid)
        self.assertFalse((self.path/'invalid').exists())


if __name__ == '__main__':
    unittest.main()
