"""Binary fixtures independently exercise full-basis WAVECAR layout detection."""
from itertools import product
from pathlib import Path
import struct
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from wavecar_fukui import KINETIC_C, Wavecar


class AutomaticSpinorTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.path = Path(temporary.name) / 'WAVECAR'

    def write(self, components=1, channels=1, *, override=None, rtag=45200, legacy=False):
        stride, nb = 4096, 2
        q = np.array([[0., 0., 0.], [.25, 0., 0.]])
        data = bytearray(stride * (2 + channels * len(q) * (nb + 1)))
        struct.pack_into('<3d', data, 0, stride // 4 if legacy else stride, channels, rtag)
        struct.pack_into('<12d', data, stride, len(q), nb, 160., *np.eye(3).ravel())
        for channel in range(channels):
            for k, point in enumerate(q):
                # Small independently bounded brute-force basis, not reader helpers.
                count = sum(np.linalg.norm((np.array(g) + point) * (2 * np.pi))**2 / KINETIC_C < 160.
                            for g in product(range(-3, 4), repeat=3)) * components
                count = (override or {}).get((channel, k), count)
                offset = stride * (2 + (channel * len(q) + k) * (nb + 1))
                struct.pack_into('<10d', data, offset, count, *point, -1., 0., 1., 2., 0., 0.)
                for band in range(nb):
                    np.frombuffer(data, dtype='<c8', count=count, offset=offset + stride * (band + 1))[:] = 1 + 2j
        self.path.write_bytes(data)

    def test_auto_scalar_spinor_and_collinear_with_checked_override(self):
        for components, channels, multiplicity in [(1, 1, 2), (2, 1, 1), (1, 2, 1)]:
            with self.subTest(components=components, channels=channels):
                self.write(components, channels)
                wave = Wavecar(self.path, spin=channels)
                explicit = Wavecar(self.path, spinor_components=components, spin=channels)
                self.assertEqual((wave.spinor_components, wave.spin_multiplicity), (components, multiplicity))
                for k in range(2):
                    np.testing.assert_array_equal(wave.coefficients(k, [1, 2]), explicit.coefficients(k, [1, 2]))
                with self.assertRaisesRegex(ValueError, 'declared spinor components'):
                    Wavecar(self.path, spinor_components=3-components)

    def test_later_k_and_unselected_channel_are_checked(self):
        self.write(1, override={(0, 1): 3})
        with self.assertRaisesRegex(ValueError, 'spin=1, k=2'):
            Wavecar(self.path)
        self.write(1, 2, override={(1, 1): 3})
        with self.assertRaisesRegex(ValueError, 'spin=2, k=2'):
            Wavecar(self.path, spin=1)

    def test_mixed_layout_and_ispin_two_spinors_are_rejected(self):
        self.write(1)
        count = int(Wavecar(self.path).nplane[1])
        self.write(1, override={(0, 1): 2 * count})
        with self.assertRaisesRegex(ValueError, 'inconsistent spinor layout'):
            Wavecar(self.path)
        self.write(2, 2)
        with self.assertRaisesRegex(ValueError, 'ISPIN=2'):
            Wavecar(self.path)

    def test_spin_channels_must_share_coordinates_even_when_basis_counts_match(self):
        self.write(1, 2)
        data = bytearray(self.path.read_bytes())
        struct.pack_into('<d', data, 4096 * 8 + 8, .01)
        self.path.write_bytes(data)
        for selected in (1, 2):
            with self.subTest(selected=selected), self.assertRaisesRegex(ValueError, 'disagree between spin channels'):
                Wavecar(self.path, spin=selected)

    def test_scalar_multiplicity_override_is_preserved(self):
        self.write(1)
        wave = Wavecar(self.path)
        self.assertEqual(wave.resolve_spin_multiplicity(), 2)
        self.assertEqual(wave.resolve_spin_multiplicity(1), 1)
        self.write(2)
        with self.assertRaisesRegex(ValueError, 'multiplicity one'):
            Wavecar(self.path).resolve_spin_multiplicity(2)

    def test_legacy_stride_and_unsupported_precision(self):
        self.write(2, legacy=True)
        self.assertEqual(Wavecar(self.path).spinor_components, 2)
        self.write(rtag=45210)
        with self.assertRaisesRegex(ValueError, 'RTAG=45210 is unsupported'):
            Wavecar(self.path)

    def test_fractional_and_nonfinite_header_fields_are_rejected(self):
        for field, value in [(0, 4096.25), (1, 1.1), (2, 45200.1), (2, float('nan'))]:
            self.write()
            data = bytearray(self.path.read_bytes())
            struct.pack_into('<d', data, 8 * field, value)
            self.path.write_bytes(data)
            with self.subTest(field=field, value=value), self.assertRaisesRegex(ValueError, 'finite integers'):
                Wavecar(self.path)


if __name__ == '__main__':
    unittest.main()
