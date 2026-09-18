"""Analytic sign, gauge, coverage and legacy-canonical-matrix checks."""
import copy
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from exported_matrix_kubo import (ARRAYS, SCHEMA, UNITS, BerryResult, DegenerateBandError,
                                 MatrixContractError, MatrixDataset, berry_curvature,
                                 read_matrix_bundle, validate_matrix_data)


def model_data(energies, matrices, *, status='experimental'):
    energies = np.asarray(energies, dtype=np.float64)
    if energies.ndim == 1:
        energies = energies[None]
    nk, nb = energies.shape
    matrices = np.asarray(matrices, dtype=np.complex128)
    if matrices.ndim == 3:
        matrices = matrices[None]
    metadata = {
        'schema': SCHEMA, 'version': 1, 'complete': True, 'source_nkpoints': nk, 'source_nbands': nb,
        'units': UNITS.copy(), 'matrix_axes': ['k', 'cartesian', 'bra_band', 'ket_band'],
        'matrix_element_convention': '<n|D_a|m>', 'cartesian_components': ['x', 'y', 'z'],
        'reciprocal_convention': '2pi', 'weights_convention': 'sum_one', 'diagonal_status': 'present',
        'operator': {'kind': 'analytic_model_hbar_velocity', 'definition': 'Explicit derivative of model H(k).',
                     'hermitian': True, 'accuracy_status': status, 'included_terms': ['entire analytic model'],
                     'missing_terms': [], 'validation_evidence': ['independent two-band Dirac oracle']},
        'provenance': {'run_id': 'synthetic-test', 'exporter_revision': 'synthetic-python-v1',
                       'source_hashes': {'model': hashlib.sha256(b'analytic-model').hexdigest()}},
    }
    return MatrixDataset(metadata, np.arange(1, nk+1, dtype=np.int64), np.arange(1, nb+1, dtype=np.int64),
                         np.zeros((nk, 3)), np.full(nk, 1/nk), energies, matrices,
                         np.ones(matrices.shape, dtype=bool), np.eye(3)*10, np.eye(3)*2*np.pi/10)


def dirac_model(kx=.3, ky=-.2, mass=.7):
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.diag([1., -1.])
    h = kx*sx+ky*sy+mass*sz
    energy, u = np.linalg.eigh(h)
    d = np.stack((u.conj().T@sx@u, u.conj().T@sy@u, np.zeros((2, 2), complex)))
    oracle = mass/(2*(kx*kx+ky*ky+mass*mass)**1.5)
    return model_data(energy, d), np.array([oracle, -oracle])


def curvature(data, n=None, m=None, **kwargs):
    defaults = dict(interest_band_ids=data.band_ids if n is None else n,
                    intermediate_band_ids=data.band_ids if m is None else m,
                    degeneracy_threshold_eV=1e-10, allow_experimental=True)
    defaults.update(kwargs)
    return berry_curvature(data, **defaults)


class MatrixKuboTests(unittest.TestCase):
    def test_dirac_sign_and_units_for_both_mass_signs(self):
        for mass in (.7, -.7):
            data, oracle = dirac_model(mass=mass)
            out = curvature(data)
            self.assertIsInstance(out, BerryResult)
            np.testing.assert_allclose(out.omega_A2[0, :, 2], oracle, atol=1e-14, rtol=0)
            np.testing.assert_array_equal(out.omega_A2[0, :, :2], 0)
            np.testing.assert_allclose(out.omega_A2.sum(axis=1), 0, atol=1e-14)
            # If derivatives scale by L and energies remain fixed, Omega scales by L^2.
            scaled = replace(data, D_eVA=data.D_eVA*3)
            np.testing.assert_allclose(curvature(scaled).omega_A2, out.omega_A2*9, atol=1e-14)

    def test_independent_band_phase_gauge_and_complex_parts(self):
        data, _ = dirac_model()
        phase = np.exp(1j*np.array([.123, -1.932]))
        gauge = data.D_eVA*phase.conj()[None, None, :, None]*phase[None, None, None, :]
        np.testing.assert_allclose(curvature(replace(data, D_eVA=gauge)).omega_A2,
                                   curvature(data).omega_A2, atol=1e-14)
        # Real-only data would destroy the Dirac curvature; it is not an equivalent format.
        with self.assertRaises(MatrixContractError):
            curvature(replace(data, D_eVA=data.D_eVA.real))
        transposed = replace(data, D_eVA=data.D_eVA.swapaxes(-1, -2).copy())
        np.testing.assert_allclose(curvature(transposed).omega_A2, -curvature(data).omega_A2, atol=1e-14)

    def test_cartesian_cyclic_component_order(self):
        data, _ = dirac_model()
        original = curvature(data).omega_A2
        # x->y, y->z changes the nonzero curvature component from z to x.
        rotated = replace(data, D_eVA=data.D_eVA[:, [2, 0, 1]])
        np.testing.assert_allclose(curvature(rotated).omega_A2[:, :, 0], original[:, :, 2], atol=1e-14)
        swapped = replace(data, D_eVA=data.D_eVA[:, [1, 0, 2]])
        np.testing.assert_allclose(curvature(swapped).omega_A2[:, :, 2], -original[:, :, 2], atol=1e-14)

    def test_hermiticity_fail_does_not_silently_symmetrize(self):
        data, _ = dirac_model()
        bad = data.D_eVA.copy(); bad[0, 0, 0, 1] += .01j
        with self.assertRaisesRegex(MatrixContractError, 'not Hermitian'):
            curvature(replace(data, D_eVA=bad))
        bad = data.D_eVA.copy(); bad[0, 0, 0, 0] += .01j
        with self.assertRaisesRegex(MatrixContractError, 'not Hermitian'):
            curvature(replace(data, D_eVA=bad))
        # E_m*dS is not automatically Hermitian even for Hermitian dH and dS.
        energies = np.array([-1., 2.]); overlap_derivative = np.array([[0, 1j], [-1j, 0]])
        unsym = data.D_eVA.copy(); unsym[0, 0] -= overlap_derivative*energies[None, :]
        with self.assertRaisesRegex(MatrixContractError, 'not Hermitian'):
            curvature(replace(data, D_eVA=unsym))

    def test_interest_and_intermediate_windows_are_independent(self):
        matrix = np.zeros((3, 3, 3), complex)
        matrix[0, 0, 1] = matrix[0, 1, 0] = 1
        matrix[1, 0, 1] = -1j; matrix[1, 1, 0] = 1j
        matrix[0, 0, 2] = matrix[0, 2, 0] = 2
        matrix[1, 0, 2] = -2j; matrix[1, 2, 0] = 2j
        data = model_data([-1., 1., 3.], matrix)
        full = curvature(data, n=[1], m=[1, 2, 3])
        partial = curvature(data, n=[1], m=[1, 2])
        np.testing.assert_allclose(full.omega_A2[0, 0, 2], -1.)
        np.testing.assert_allclose(partial.omega_A2[0, 0, 2], -.5)
        self.assertTrue(full.diagnostics['intermediate_contains_all_source_bands'])
        self.assertFalse(partial.diagnostics['intermediate_contains_all_source_bands'])
        for n, m in (([1], [1]), ([1], [2, 2]), ([4], [1, 2]), ([1.], [1, 2]), ([True], [1, 2])):
            with self.subTest(n=n, m=m), self.assertRaises(MatrixContractError):
                curvature(data, n=n, m=m)

    def test_degeneracy_error_or_nan_with_explicit_threshold(self):
        data, _ = dirac_model()
        data = replace(data, energies_eV=np.array([[0., .001529]]))
        self.assertTrue(np.isfinite(curvature(data).omega_A2).all())
        with self.assertRaises(DegenerateBandError):
            curvature(data, degeneracy_threshold_eV=.002)
        masked = curvature(data, degeneracy_threshold_eV=.002, degeneracy_policy='mask')
        self.assertFalse(masked.valid_nondegenerate.any())
        self.assertTrue(np.isnan(masked.omega_A2).all())
        self.assertEqual(masked.diagnostics['invalid_individual_band_count'], 2)
        self.assertEqual(masked.diagnostics['degenerate_pairs'][0]['gap_eV'], .001529)
        exact = replace(data, energies_eV=np.zeros((1, 2)))
        with self.assertRaises(DegenerateBandError):
            curvature(exact, degeneracy_threshold_eV=0)

    def test_missing_diagonal_allowed_but_missing_interband_never_zeroed(self):
        data, _ = dirac_model()
        matrix, coverage = data.D_eVA.copy(), data.coverage.copy()
        for n in range(2):
            matrix[:, :, n, n] = complex(np.nan, np.nan); coverage[:, :, n, n] = False
        meta = copy.deepcopy(data.metadata); meta['diagonal_status'] = 'unavailable'
        nodiag = replace(data, metadata=meta, D_eVA=matrix, coverage=coverage)
        np.testing.assert_allclose(curvature(nodiag).omega_A2, curvature(data).omega_A2, atol=1e-14)
        matrix[:, :, 0, 1] = matrix[:, :, 1, 0] = complex(np.nan, np.nan)
        coverage[:, :, 0, 1] = coverage[:, :, 1, 0] = False
        with self.assertRaisesRegex(MatrixContractError, 'missing required element'):
            curvature(replace(nodiag, D_eVA=matrix, coverage=coverage))
        matrix[~coverage] = 0
        with self.assertRaisesRegex(MatrixContractError, 'explicit complex NaN'):
            curvature(replace(nodiag, D_eVA=matrix, coverage=coverage))

    def test_experimental_operator_requires_explicit_diagnostic_override(self):
        data, _ = dirac_model()
        with self.assertRaisesRegex(MatrixContractError, 'experimental operator rejected'):
            curvature(data, allow_experimental=False)
        meta = copy.deepcopy(data.metadata)
        meta['operator']['kind'] = 'derived_from_vasp_paw_optical_connection'
        result = curvature(replace(data, metadata=meta))
        self.assertEqual(result.diagnostics['operator_accuracy_status'], 'experimental')
        self.assertEqual(result.diagnostics['operator_kind'], 'derived_from_vasp_paw_optical_connection')
        validated = copy.deepcopy(meta); validated['operator']['accuracy_status'] = 'validated'
        curvature(replace(data, metadata=validated), allow_experimental=False)
        validated['operator']['validation_evidence'] = []
        with self.assertRaises(MatrixContractError):
            curvature(replace(data, metadata=validated), allow_experimental=False)

    def test_canonical_plane_wave_kernel_matches_half_legacy_circular_formula(self):
        rng = np.random.default_rng(9321)
        nk, nb, ng = 2, 4, 7
        energies = np.array([[-2., -.3, .7, 2.], [-1.8, -.1, .8, 2.2]])
        # Spinor coefficients intentionally need not be orthonormal, like PAW pseudo states.
        coeff = rng.normal(size=(nk, 2*ng, nb))+1j*rng.normal(size=(nk, 2*ng, nb))
        coeff /= np.sqrt(np.sum(abs(coeff)**2, axis=1))[:, None, :]
        kg = rng.normal(size=(nk, ng, 3))
        unit = 58.063843708134961  # Actual legacy Fortran runtime U; eV^2 Angstrom^4.
        matrices = np.empty((nk, 3, nb, nb), complex)
        legacy = np.zeros((nk, nb))
        for k in range(nk):
            for a in range(3):
                g = np.tile(kg[k, :, a], 2)
                matrices[k, a] = np.sqrt(unit)*(coeff[k].conj().T @ (g[:, None]*coeff[k]))
            # Independent scalar loops reproduce legacy spinor plane-wave accumulation.
            for n in range(nb):
                for m in range(nb):
                    if n == m:
                        continue
                    x, y = 0j, 0j
                    for g in range(ng):
                        overlap = (coeff[k, g, n].conj()*coeff[k, g, m]
                                   +coeff[k, g+ng, n].conj()*coeff[k, g+ng, m])
                        x += kg[k, g, 0]*overlap; y += kg[k, g, 1]*overlap
                    legacy[k, n] -= (abs(x+1j*y)**2-abs(x-1j*y)**2)*unit/(energies[k, n]-energies[k, m])**2
        result = curvature(model_data(energies, matrices))
        np.testing.assert_allclose(result.omega_A2[:, :, 2], legacy/2, atol=2e-14, rtol=1e-13)
        np.testing.assert_allclose(result.omega_A2.sum(axis=1), 0, atol=1e-13)

    def test_bundle_roundtrip_checksum_and_complete_array_contract(self):
        data, _ = dirac_model()
        with tempfile.TemporaryDirectory() as tmp:
            path, meta_path = Path(tmp)/'matrix.npz', Path(tmp)/'matrix.json'
            np.savez_compressed(path, **{key: getattr(data, key) for key in ARRAYS})
            meta = copy.deepcopy(data.metadata)
            meta['matrix_npz_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
            meta_path.write_text(json.dumps(meta))
            read = read_matrix_bundle(path, meta_path, allow_experimental=True)
            np.testing.assert_array_equal(read.D_eVA, data.D_eVA)
            np.testing.assert_array_equal(curvature(read).omega_A2, curvature(data).omega_A2)
            with self.assertRaises(MatrixContractError):
                read_matrix_bundle(path, meta_path)
            meta['matrix_npz_sha256'] = '0'*64; meta_path.write_text(json.dumps(meta))
            with self.assertRaisesRegex(MatrixContractError, 'SHA256'):
                read_matrix_bundle(path, meta_path, allow_experimental=True)
            # Rehash a deliberately incomplete archive: checksum alone is insufficient.
            np.savez_compressed(path, **{key: getattr(data, key) for key in ARRAYS if key != 'coverage'})
            meta['matrix_npz_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest(); meta_path.write_text(json.dumps(meta))
            with self.assertRaisesRegex(MatrixContractError, 'missing required arrays'):
                read_matrix_bundle(path, meta_path, allow_experimental=True)

    def test_bad_units_precision_ids_weights_and_nonfinite_values_fail(self):
        data, _ = dirac_model()
        for field, value in (('complete', False), ('version', 2), ('units', {'matrix': 'm/s'}),
                             ('matrix_element_convention', '<m|D_a|n>'), ('reciprocal_convention', 'no2pi')):
            meta = copy.deepcopy(data.metadata); meta[field] = value
            with self.subTest(field=field), self.assertRaises(MatrixContractError):
                curvature(replace(data, metadata=meta))
        for replacement in ({'D_eVA': data.D_eVA.astype(np.complex64)}, {'band_ids': np.array([1, 1], dtype=np.int64)},
                            {'energies_eV': np.array([[np.nan, 1.]])}, {'weights': np.array([.5])},
                            {'reciprocal_inv_A': np.eye(3)}, {'coverage': np.ones((1, 3, 2, 2), dtype=float)}):
            with self.subTest(replacement=list(replacement)), self.assertRaises(MatrixContractError):
                curvature(replace(data, **replacement))


if __name__ == '__main__':
    unittest.main()
