"""Native Fortran layout, physical conversion, corruption, and no-repair tests."""
from dataclasses import replace
import json
from pathlib import Path
import struct
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from vasp_optical_export import (FOOTER, MAGIC, OpticalExportError, compare_waveder_lower,
                                 hermiticity_stats, main, read_connection, read_waveder, _provenance,
                                 write_bundles)
from exported_matrix_kubo import MatrixContractError, berry_curvature, read_matrix_bundle


def synthetic_arrays(ns=2, nk=2, nb=3):
    # A non-orthogonal lattice exposes accidental A-vs-A.T conventions.
    lattice = np.array([[3., .2, 0.], [.4, 4., 0.], [.1, .3, 12.]])
    reciprocal = 2*np.pi*np.linalg.inv(lattice).T
    kpoints = np.array([[i/(nk+1), .2*i, .1*i] for i in range(nk)])
    weights = np.arange(1, nk+1, dtype=float); weights /= weights.sum()
    energies = np.empty((ns, nk, nb))
    occupations = np.zeros_like(energies)
    derivative = np.empty((ns, nk, 3, nb))
    c = np.empty((ns, nk, 3, nb, nb), dtype=np.complex128)
    for s in range(ns):
        for k in range(nk):
            energies[s, k] = np.arange(nb)*1.3 + .1*s + .02*k
            occupations[s, k, 0] = 1
            for a in range(3):
                for m in range(nb):
                    derivative[s, k, a, m] = 1+s+2*k+3*a+.01*m
                    for n in range(nb):
                        c[s, k, a, m, n] = 10000*s+1000*k+100*a+10*m+n + 1j*(m-n)
    return dict(lattice_A=lattice, reciprocal_inv_A=reciprocal, kpoints_fractional=kpoints,
                weights=weights, energies_eV=energies, occupations=occupations,
                energy_der_eVA=derivative, C_A=c)


def write_connection(path, arrays, threshold=1e-10, fermi=.35):
    ns, nk, nb = arrays['energies_eV'].shape
    chunks = [MAGIC, struct.pack('<7i', 1, nk, nb, ns, 3, 8, 16), struct.pack('<2d', threshold, fermi)]
    # Explicit source-to-Fortran orientation; distinct sentinels above make every
    # axis observable independently of the parser implementation.
    for key, axes, dtype in [('lattice_A', (1, 0), '<f8'), ('reciprocal_inv_A', (1, 0), '<f8'),
                             ('kpoints_fractional', (1, 0), '<f8'), ('weights', (0,), '<f8'),
                             ('energies_eV', (2, 1, 0), '<f8'), ('occupations', (2, 1, 0), '<f8'),
                             ('energy_der_eVA', (3, 1, 0, 2), '<f8'), ('C_A', (3, 4, 1, 0, 2), '<c16')]:
        chunks.append(np.asarray(arrays[key].transpose(axes), dtype=dtype).tobytes(order='F'))
    path.write_bytes(b''.join(chunks)+FOOTER)


def write_waveder(path, c):
    ns, nk, _, nb, nd = c.shape
    def record(data):
        return struct.pack('<i', len(data))+data+struct.pack('<i', len(data))
    records = [struct.pack('<4i', nb, nd, nk, ns), struct.pack('<d', -.2),
               np.arange(9, dtype='<f8').tobytes(),
               np.asarray(c.transpose(3, 4, 1, 0, 2), dtype='<c8').tobytes(order='F')]
    path.write_bytes(b''.join(record(data) for data in records))


def two_band_arrays():
    arrays = synthetic_arrays(ns=1, nk=1, nb=2)
    arrays['energies_eV'][:] = [-1., 1.]
    arrays['C_A'][:] = 0
    # Independent exact two-level oracle: Dx=sigma_x, Dy=sigma_y.
    arrays['C_A'][0, 0, 0] = [[0, -.5j], [.5j, 0]]
    arrays['C_A'][0, 0, 1] = [[0, -.5], [-.5, 0]]
    arrays['energy_der_eVA'][:] = 0
    return arrays


class OpticalExportTests(unittest.TestCase):
    def test_all_axes_and_nonorthogonal_lattice_are_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'connection.bin'
            arrays = synthetic_arrays()
            write_connection(path, arrays)
            data = read_connection(path)
            for key, value in arrays.items():
                np.testing.assert_array_equal(getattr(data, key), value, err_msg=key)
            self.assertEqual(data.C_A[1, 1, 2, 1, 2], 11212-1j)
            self.assertEqual(data.degeneracy_threshold_eV, 1e-10)
            self.assertEqual((data.nspin, data.nkpoints, data.nbands), (2, 2, 3))
            self.assertEqual(len(data.source_sha256), 64)

    def test_velocity_sign_diagonal_and_kubo_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'connection.bin'; arrays = two_band_arrays()
            write_connection(path, arrays)
            data = read_connection(path)
            np.testing.assert_allclose(data.D_eVA[0, 0, 0], [[0, 1], [1, 0]], atol=0)
            np.testing.assert_allclose(data.D_eVA[0, 0, 1], [[0, -1j], [1j, 0]], atol=0)
            provenance = {'run_id': 'analytic-serialization-test', 'exporter_revision': 'fixture',
                          'source_hashes': {'connection': data.source_sha256}}
            out = Path(tmp)/'converted'
            report = write_bundles(data, out, provenance)
            self.assertTrue(report['D_eVA_hermiticity']['passed_default_kubo_tolerance'])
            parsed = read_matrix_bundle(out/'matrix-spin1.npz', out/'matrix-spin1.json', allow_experimental=True)
            result = berry_curvature(parsed, interest_band_ids=[1, 2], intermediate_band_ids=[1, 2],
                                    degeneracy_threshold_eV=1e-10, allow_experimental=True)
            np.testing.assert_allclose(result.omega_A2[0, :, 2], [-.5, .5], atol=1e-15)
            with self.assertRaisesRegex(MatrixContractError, 'experimental'):
                read_matrix_bundle(out/'matrix-spin1.npz', out/'matrix-spin1.json')
            with self.assertRaisesRegex(OpticalExportError, 'new'):
                write_bundles(data, out, provenance)
            changed = replace(data, energy_der_eVA=np.ones_like(data.energy_der_eVA)*7)
            np.testing.assert_array_equal(np.diagonal(changed.D_eVA, axis1=-2, axis2=-1), 7)
            np.testing.assert_array_equal(np.diagonal(changed.raw_D_eVA, axis1=-2, axis2=-1), 0)

    def test_exact_degenerate_pairs_missing_but_diagonal_present(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'connection.bin'; arrays = synthetic_arrays(ns=1, nk=1)
            arrays['energies_eV'][0, 0] = [0., 1e-11, 1.]
            write_connection(path, arrays)
            data = read_connection(path)
            self.assertFalse(data.coverage[0, 0, :, 0, 1].any())
            self.assertFalse(data.coverage[0, 0, :, 1, 0].any())
            self.assertTrue(np.isnan(data.D_eVA[0, 0, :, 0, 1].real).all())
            self.assertTrue(np.isnan(data.D_eVA[0, 0, :, 0, 1].imag).all())
            self.assertTrue(data.coverage[0, 0, :, 0, 0].all())
            self.assertTrue(np.isfinite(data.raw_D_eVA).all())
            np.testing.assert_array_equal(data.C_A, arrays['C_A'])

    def test_transitive_native_clusters_are_unavailable_per_spin_and_k(self):
        threshold = 1e-10
        arrays = synthetic_arrays(ns=2, nk=2, nb=3)
        arrays['energies_eV'][:] = np.array([
            [[0., .75*threshold, 1.5*threshold], [0., 1.01*threshold, 3*threshold]],
            [[0., threshold, 2*threshold], [0., .5*threshold, 4*threshold]],
        ])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'connection.bin'; write_connection(path, arrays)
            data = read_connection(path)
            expected = np.array([[np.eye(3, dtype=bool), np.ones((3, 3), dtype=bool)],
                                 [np.eye(3, dtype=bool), [[True, False, True],
                                                        [False, True, True], [True, True, True]]]])
            np.testing.assert_array_equal(data.coverage, np.broadcast_to(expected[:, :, None], data.C_A.shape))
            self.assertTrue(np.isnan(data.D_eVA[~data.coverage]).all())
            self.assertTrue(np.isfinite(data.raw_D_eVA).all())
            np.testing.assert_array_equal(data.C_A, arrays['C_A'])
            e = arrays['energies_eV']
            np.testing.assert_array_equal(data.raw_D_eVA,
                1j*(e[:, :, None, :]-e[:, :, :, None])[:, :, None]*arrays['C_A'])
            np.testing.assert_array_equal(np.diagonal(data.D_eVA, axis1=-2, axis2=-1), arrays['energy_der_eVA'])

    def test_native_denominator_literal_precision_and_strict_boundary(self):
        threshold = 1e-10
        promoted = float(np.float32(1e-10))
        gaps = np.array([threshold, (threshold+promoted)/2, promoted, np.nextafter(promoted, np.inf)])
        arrays = synthetic_arrays(ns=1, nk=4, nb=2)
        arrays['energies_eV'][0, :, 0] = 0.
        arrays['energies_eV'][0, :, 1] = gaps
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'connection.bin'; write_connection(path, arrays)
            data = read_connection(path)
            np.testing.assert_array_equal(data.coverage[0, :, 0, 0, 1], [False, False, True, True])
            np.testing.assert_array_equal(data.coverage, data.coverage.swapaxes(-1, -2))
            self.assertTrue(np.diagonal(data.coverage, axis1=-2, axis2=-1).all())

    def test_sparse_band_selection_retains_full_native_cluster_mask(self):
        arrays = synthetic_arrays(ns=1, nk=1, nb=3)
        arrays['energies_eV'][0, 0] = [0., .75e-10, 1.5e-10]
        arrays['C_A'][:] = 0  # Producer zeroed the whole transitive cluster.
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'connection.bin'; write_connection(path, arrays)
            data = read_connection(path); out = Path(tmp)/'converted'
            write_bundles(data, out, {'run_id': 'cluster-test', 'exporter_revision': 'fixture',
                                      'source_hashes': {'connection': data.source_sha256}})
            full = read_matrix_bundle(out/'matrix-spin1.npz', out/'matrix-spin1.json', allow_experimental=True)
            selected = np.array([0, 2])
            sparse = replace(full, band_ids=full.band_ids[selected], energies_eV=full.energies_eV[:, selected],
                D_eVA=full.D_eVA[:, :, selected][:, :, :, selected],
                coverage=full.coverage[:, :, selected][:, :, :, selected])
            self.assertFalse(sparse.coverage[0, :, 0, 1].any())
            with self.assertRaisesRegex(MatrixContractError, 'missing required element'):
                berry_curvature(sparse, interest_band_ids=[1], intermediate_band_ids=[3],
                                degeneracy_threshold_eV=1e-10, allow_experimental=True)
            with np.load(out/'raw-connection.npz') as raw:
                np.testing.assert_array_equal(raw['C_A'], arrays['C_A'])
                np.testing.assert_array_equal(raw['raw_D_eVA'], data.raw_D_eVA)

    def test_nonhermitian_evidence_is_not_mirrored_or_averaged(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'connection.bin'; arrays = two_band_arrays()
            arrays['C_A'][0, 0, 0, 0, 1] += .07j
            write_connection(path, arrays)
            data = read_connection(path)
            stats = hermiticity_stats(data.D_eVA, data.coverage)
            self.assertAlmostEqual(stats['max_absolute_residual'], .14)
            self.assertFalse(stats['repair_applied'])
            self.assertFalse(stats['passed_default_kubo_tolerance'])
            provenance = {'run_id': 'asymmetry-test', 'exporter_revision': 'fixture',
                          'source_hashes': {'connection': data.source_sha256}}
            out = Path(tmp)/'converted'; write_bundles(data, out, provenance)
            with np.load(out/'raw-connection.npz') as raw:
                np.testing.assert_array_equal(raw['C_A'], arrays['C_A'])
                np.testing.assert_array_equal(raw['raw_D_eVA'], data.raw_D_eVA)
            with self.assertRaisesRegex(MatrixContractError, 'not Hermitian'):
                read_matrix_bundle(out/'matrix-spin1.npz', out/'matrix-spin1.json', allow_experimental=True)

    def test_header_byte_length_footer_and_endian_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'connection.bin'; write_connection(path, two_band_arrays())
            valid = path.read_bytes()
            cases = [valid[:20], valid[:-1], valid+b'\0', b'X'+valid[1:], valid[:-1]+b'X']
            for offset, value in [(32, 2), (36, 0), (40, -1), (44, 3), (48, 4), (52, 4), (56, 8)]:
                bad = bytearray(valid); struct.pack_into('<i', bad, offset, value); cases.append(bad)
            endian = bytearray(valid); struct.pack_into('>7i', endian, 32, 1, 1, 2, 1, 3, 8, 16); cases.append(endian)
            for bad in cases:
                with self.subTest(prefix=bytes(bad[:60])):
                    path.write_bytes(bad)
                    with self.assertRaises(OpticalExportError):
                        read_connection(path)

    def test_nonfinite_values_and_bad_normalization_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'connection.bin'
            for key in synthetic_arrays():
                arrays = synthetic_arrays(); arrays[key].flat[0] = np.nan
                write_connection(path, arrays)
                with self.subTest(key=key), self.assertRaises(OpticalExportError):
                    read_connection(path)
            arrays = synthetic_arrays(); arrays['weights'] *= 2; write_connection(path, arrays)
            with self.assertRaisesRegex(OpticalExportError, 'weights'):
                read_connection(path)
            write_connection(path, synthetic_arrays(), threshold=np.inf)
            with self.assertRaisesRegex(OpticalExportError, 'nonfinite'):
                read_connection(path)
            write_connection(path, synthetic_arrays(), threshold=-1)
            with self.assertRaisesRegex(OpticalExportError, 'negative'):
                read_connection(path)
            for threshold in (0., .001, .0021):
                write_connection(path, synthetic_arrays(), threshold=threshold)
                with self.assertRaisesRegex(OpticalExportError, 'threshold must equal'):
                    read_connection(path)

    def test_spin_mode_threshold_masks_complete_clusters_and_preserves_raw_values(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'connection.bin'
            arrays = synthetic_arrays(ns=1, nk=1, nb=4)
            arrays['energies_eV'][0, 0] = [0., .0015, .003, 1.]
            write_connection(path, arrays, threshold=.002)
            data = read_connection(path)
            self.assertEqual(data.degeneracy_threshold_eV, .002)
            # The outer pair spans 3meV but belongs to the same transitive group.
            self.assertFalse(data.coverage[0, 0, :, 0, 2].any())
            self.assertTrue(data.coverage[0, 0, :, 0, 3].all())
            self.assertTrue(np.isnan(data.D_eVA[0, 0, :, 0, 2]).all())
            np.testing.assert_array_equal(data.C_A, arrays['C_A'])
            np.testing.assert_array_equal(np.diagonal(data.D_eVA, axis1=-2, axis2=-1), arrays['energy_der_eVA'])
            out = Path(tmp)/'new'
            report = write_bundles(data, out, {'run_id': 'threshold fixture'})
            self.assertEqual(report['unavailable_offdiagonal_gap_threshold_eV'], .002)
            saved = json.loads((out/'matrix-spin1.json').read_text())
            self.assertEqual(saved['missing_offdiagonal_gap_threshold_eV'], .002)
            # Legacy LBERRY-only producer remains readable with its original mask.
            write_connection(path, arrays, threshold=1e-10)
            legacy = read_connection(path)
            self.assertTrue(legacy.coverage[0, 0, :, 0, 2].all())

    def test_native_waveder_layout_lower_only_and_complex64_rounding(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'connection.bin'; wpath = Path(tmp)/'WAVEDER'
            arrays = synthetic_arrays(); write_connection(path, arrays)
            data = read_connection(path)
            native = arrays['C_A'][..., :2].copy()
            native[..., 0, 1] = 99+100j  # Simulate altered native upper triangle.
            write_waveder(wpath, native)
            wave = read_waveder(wpath)
            self.assertEqual(wave.C_A.dtype, np.dtype('complex64'))
            np.testing.assert_array_equal(wave.C_A, native.astype(np.complex64))
            stats = compare_waveder_lower(data, wave)
            self.assertEqual(stats['max_absolute_difference_A'], 0)
            self.assertEqual(stats['compared_elements'], 36)
            zero_reference = replace(wave, C_A=np.zeros_like(wave.C_A))
            zero_stats = compare_waveder_lower(data, zero_reference)
            self.assertGreater(zero_stats['max_absolute_difference_A'], 0)
            self.assertIsNone(zero_stats['relative_frobenius_difference'])
            bad = replace(wave, C_A=wave.C_A[:, :1])
            with self.assertRaisesRegex(OpticalExportError, 'dimensions'):
                compare_waveder_lower(data, bad)

    def test_waveder_corrupt_record_markers_truncation_and_extra_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'WAVEDER'; write_waveder(path, two_band_arrays()['C_A'])
            valid = path.read_bytes()
            cases = [valid[:3], valid[:-1], valid+b'junk']
            for offset, value in [(0, -16), (0, 32), (20, 17)]:
                bad = bytearray(valid); struct.pack_into('<i', bad, offset, value); cases.append(bad)
            for bad in cases:
                path.write_bytes(bad)
                with self.assertRaises(OpticalExportError):
                    read_waveder(path)

    def test_cli_run_provenance_and_multi_spin_outputs(self):
        from contextlib import redirect_stdout
        import io
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)/'run'; run.mkdir()
            write_connection(run/'BERRY_CONNECTION.bin', synthetic_arrays())
            data = read_connection(run/'BERRY_CONNECTION.bin')
            record = {'binary_sha256': 'a'*64, 'input_sha256': {'INCAR': 'b'*64}, 'nbands': 3, 'nkpoints': 2,
                      'completed_export': True, 'returncode': 0, 'normal_timing_footer': True,
                      'outputs': {'BERRY_CONNECTION.bin': {'sha256': data.source_sha256}}}
            (run/'run.json').write_text(json.dumps(record))
            out = Path(tmp)/'converted'
            with redirect_stdout(io.StringIO()):
                main(['--run-dir', str(run), '--output-dir', str(out)])
            self.assertTrue((out/'matrix-spin2.npz').is_file())
            meta = json.loads((out/'matrix-spin2.json').read_text())
            self.assertEqual(meta['spin_channel_1based'], 2)
            self.assertEqual(meta['operator']['accuracy_status'], 'experimental')
            self.assertEqual(meta['provenance']['source_hashes']['binary_sha256'], 'a'*64)
            self.assertFalse(meta['repair_applied'])
            unbound = dict(record, outputs={})
            (run/'run.json').write_text(json.dumps(unbound))
            with self.assertRaisesRegex(OpticalExportError, 'producer-recorded SHA256'):
                _provenance(run, data)
            for key, bad_value in (('completed_export', False), ('returncode', 1),
                                   ('normal_timing_footer', False), ('timed_out', True)):
                bad = dict(record); bad[key] = bad_value
                (run/'run.json').write_text(json.dumps(bad))
                with self.subTest(key=key), self.assertRaises(OpticalExportError):
                    _provenance(run, data)
            record['outputs']['BERRY_CONNECTION.bin']['sha256'] = '0'*64
            (run/'run.json').write_text(json.dumps(record))
            with self.assertRaisesRegex(OpticalExportError, 'output SHA256 mismatch'):
                _provenance(run, data)


if __name__ == '__main__':
    unittest.main()
