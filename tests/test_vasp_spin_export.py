"""Portable tests of the opt-in PAW spin/velocity interchange boundary."""
import json
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

import vasp_spin_export as export


def stream():
    nk, nb = 2, 2
    lattice = np.array([[3., 0, 0], [1., 4., 0], [0, 0, 12.]])
    reciprocal = 2*np.pi*np.linalg.inv(lattice).T
    pseudo = np.zeros((nb, nb, nk, 1, 4), dtype=np.complex128)
    pseudo[:, :, :, 0, 0] = np.eye(nb)[:, :, None]
    pseudo[:, :, :, 0, 1] = np.array([[0, 1j], [-1j, 0]])[:, :, None]
    parts = [export.MAGIC, struct.pack('<7i', 1, nk, nb, 1, 3, 8, 16)]
    for a in [np.array([.3]), lattice.T, reciprocal.T,
              np.array([[0, 0, 0], [.25, .375, 0]]).T, np.array([.5, .5]),
              np.array([[-1., -2.], [1., 2.]]).reshape(nb, nk, 1), np.ones((nb, nk, 1))]:
        parts.append(np.asarray(a, dtype='<f8').tobytes(order='F'))
    for a in [pseudo, np.zeros_like(pseudo)] + [np.zeros((nb, nb, nk, 1, 3), dtype=np.complex128)]*3:
        parts.append(np.asarray(a, dtype='<c16').tobytes(order='F'))
    return b''.join(parts)+export.FOOTER


class SpinVelocityExportTests(unittest.TestCase):
    def test_fortran_order_bra_ket_complex_and_nonorthogonal_lattice(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'stream'; path.write_bytes(stream())
            data = export.read_spin_velocity(path)
        np.testing.assert_array_equal(data['energies_eV'], [[-1, 1], [-2, 2]])
        self.assertEqual(data['pseudo'][1, 1, 0, 1], 1j)
        self.assertEqual(data['pseudo'][1, 1, 1, 0], -1j)
        np.testing.assert_array_equal(data['lattice_A'], [[3, 0, 0], [1, 4, 0], [0, 0, 12]])

    def test_bad_header_footer_size_and_nonfinite_reject(self):
        good = stream()
        bad = [b'x'+good[1:], good[:-1], good[:-32]+b'x'*32, good+b'x',
               good[:32]+struct.pack('<7i', 1, 2**30, 2**30, 1, 3, 8, 16)+good[60:],
               good[:60]+struct.pack('<d', float('nan'))+good[68:]]
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'stream'
            for value in bad:
                with self.subTest(size=len(value)):
                    path.write_bytes(value)
                    with self.assertRaises(ValueError): export.read_spin_velocity(path)

    def test_wrong_reciprocal_and_negative_weight_reject(self):
        good = stream()
        # Header68 + direct lattice72; reciprocal begins at140.
        corrupt = bytearray(good); corrupt[140:148] = struct.pack('<d', 200.)
        # k coordinates follow both lattices; two weights start at260.
        negative = bytearray(good); negative[260:268] = struct.pack('<d', -.5)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'stream'
            for value in [corrupt, negative]:
                path.write_bytes(value)
                with self.assertRaises(ValueError): export.read_spin_velocity(path)

    def test_failed_or_incomplete_run_never_reaches_matrix_audit(self):
        with tempfile.TemporaryDirectory() as folder:
            case = Path(folder)
            for n in ['SPIN_VELOCITY.bin', 'BERRY_CONNECTION.bin', 'WAVECAR', 'OUTCAR', 'INCAR', 'POSCAR', 'KPOINTS']:
                (case/n).write_text('test')
            for status, footer in [('TIMEOUT', True), ('FINISHED', False), ('FAILED_OUTPUT', False)]:
                (case/'run.json').write_text(json.dumps(dict(status=status, returncode=0, converged=True, normal_footer=footer)))
                with self.assertRaisesRegex(ValueError, 'successful producer'):
                    export.audit_producer(case)

    def test_validation_before_output_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder)/'output'
            with patch.object(export, 'audit_producer', side_effect=ValueError('invalid producer')):
                with self.assertRaisesRegex(ValueError, 'invalid producer'): export.convert_run('unused', out)
            self.assertFalse(out.exists())
            out.mkdir(); (out/'keep').write_text('unchanged')
            with patch.object(export, 'audit_producer') as audit:
                with self.assertRaisesRegex(ValueError, 'exists'): export.convert_run('unused', out)
                audit.assert_not_called()
            self.assertEqual((out/'keep').read_text(), 'unchanged')

    def test_supported_settings_and_optional_unsupported_branches(self):
        incar = '\n'.join(f'{key} = {value}' for key, value in {
            'LOPTICS': 'T', 'LBERRY_EXPORT': 'T', 'LSPIN_EXPORT': 'T',
            'LPEAD': 'F', 'LNABLA': 'F', 'LREAL': 'F', 'LSORBIT': 'T',
            'SAXIS': '0 0 1'}.items())+'\n'
        outcar = ('vasp.5.4.4\nEDIFF is reached\nGeneral timing and accounting\n'
                  'Euler angles ALPHA= 0.0 BETA= 0.0\n')
        effective = {'ISYM': -1, 'ISPIN': 1, 'NSW': 0,
                     'LSORBIT': True, 'LNONCOLLINEAR': True}
        with tempfile.TemporaryDirectory() as folder:
            case = Path(folder)
            with patch.object(export, 'outcar_value', side_effect=lambda text, key, cast: effective.get(key, False)):
                (case/'INCAR').write_text(incar)
                (case/'OUTCAR').write_text(outcar)
                export._settings(case)  # VASP 5.4.4 omits these false defaults.
                for extra in ('LDAU = T', 'LSPIRAL = T', 'QSPIRAL = .1 0 0',
                              'QSPIRAL = invalid', 'QSPIRAL = 0 0'):
                    with self.subTest(incar=extra):
                        (case/'INCAR').write_text(incar+extra+'\n')
                        with self.assertRaises(ValueError): export._settings(case)
                (case/'INCAR').write_text(incar)
                for extra in ('LDAU = T', 'LSPIRAL = .TRUE.'):
                    with self.subTest(outcar=extra):
                        (case/'OUTCAR').write_text(outcar+extra+'\n')
                        with self.assertRaises(ValueError): export._settings(case)


if __name__ == '__main__':
    unittest.main()
