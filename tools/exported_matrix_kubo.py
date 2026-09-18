#!/usr/bin/env python3
"""Strict, exporter-independent Hermitian interband matrix Kubo calculator.

Interchange v1 uses JSON metadata and an NPZ with arrays:
  k_ids:int64[K], band_ids:int64[B], kpoints_fractional:float64[K,3],
  weights:float64[K] (sum one), energies_eV:float64[K,B],
  D_eVA:complex128[K,3,B,B], coverage:bool[K,3,B,B],
  lattice_A:float64[3,3], reciprocal_inv_A:float64[3,3] (includes 2*pi).
Here D[k,a,n,m]=<n|D_a|m> in eV*Angstrom, with Cartesian a=x,y,z.
Covered elements are finite; unavailable elements are explicit complex NaN.

Only a declared Hermitian operator is accepted. In particular, an unsymmetrized
dH-E_m*dS array is not automatically a Hermitian physical hbar-velocity vertex.
The metadata accuracy label is preserved, not established by this calculator.
Experimental operators require an explicit diagnostic override.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

SCHEMA = 'vaspberry.interband-matrix'
SCHEMA_VERSION = 1
ARRAYS = ('k_ids', 'band_ids', 'kpoints_fractional', 'weights', 'energies_eV',
          'D_eVA', 'coverage', 'lattice_A', 'reciprocal_inv_A')
UNITS = {'energy': 'eV', 'matrix': 'eV*Angstrom', 'lattice': 'Angstrom', 'reciprocal': '1/Angstrom'}


class MatrixContractError(ValueError):
    """Unsupported, incomplete or internally inconsistent matrix data."""


class DegenerateBandError(ValueError):
    """Individual-band Abelian curvature is undefined at a declared cluster."""


@dataclass(frozen=True)
class MatrixDataset:
    metadata: dict
    k_ids: np.ndarray
    band_ids: np.ndarray
    kpoints_fractional: np.ndarray
    weights: np.ndarray
    energies_eV: np.ndarray
    D_eVA: np.ndarray
    coverage: np.ndarray
    lattice_A: np.ndarray
    reciprocal_inv_A: np.ndarray


@dataclass(frozen=True)
class BerryResult:
    k_ids: np.ndarray
    band_ids: np.ndarray
    intermediate_band_ids: np.ndarray
    omega_A2: np.ndarray               # (k, interest band, Cartesian x/y/z)
    min_gap_eV: np.ndarray             # minimum |En-El| over all exported l != n
    valid_nondegenerate: np.ndarray    # (k, interest band)
    diagnostics: dict


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def _fail(condition, message):
    if condition:
        raise MatrixContractError(message)


def validate_matrix_data(data: MatrixDataset, *, allow_experimental=False,
                         hermitian_atol_eVA=1e-10, hermitian_rtol=1e-8):
    """Validate storage, coverage and Hermiticity without repairing any values."""
    meta = data.metadata
    _fail(not isinstance(meta, dict) or meta.get('schema') != SCHEMA
          or type(meta.get('version')) is not int or meta.get('version') != SCHEMA_VERSION,
          f'expected {SCHEMA} version {SCHEMA_VERSION}')
    _fail(meta.get('complete') is not True, 'matrix export is not marked complete')
    _fail(meta.get('units') != UNITS, f'unit contract must be exactly {UNITS}')
    _fail(meta.get('matrix_axes') != ['k', 'cartesian', 'bra_band', 'ket_band'], 'ambiguous matrix axes')
    _fail(meta.get('matrix_element_convention') != '<n|D_a|m>', 'ambiguous bra/ket convention')
    _fail(meta.get('cartesian_components') != ['x', 'y', 'z'], 'Cartesian component order must be x,y,z')
    _fail(meta.get('reciprocal_convention') != '2pi', 'reciprocal lattice must include 2*pi')
    _fail(meta.get('weights_convention') != 'sum_one', 'normalized k-weight convention required')
    operator = meta.get('operator', {})
    _fail(not isinstance(operator, dict) or operator.get('hermitian') is not True,
          'explicitly Hermitian operator definition required; no silent symmetrization')
    for key in ('kind', 'definition'):
        _fail(not isinstance(operator.get(key), str) or not operator[key].strip(), f'operator.{key} required')
    for key in ('included_terms', 'missing_terms', 'validation_evidence'):
        value = operator.get(key)
        _fail(not isinstance(value, list) or any(not isinstance(x, str) or not x.strip() for x in value),
              f'operator.{key} must be an explicit list of text descriptions')
    status = operator.get('accuracy_status')
    _fail(status not in ('experimental', 'validated'), 'operator accuracy_status must be experimental or validated')
    _fail(status == 'experimental' and allow_experimental is not True,
          'experimental operator rejected; allow_experimental=True permits labeled diagnostic calculation only')
    _fail(status == 'validated' and not operator['validation_evidence'],
          'validated accuracy label requires explicit validation evidence/provenance')
    provenance = meta.get('provenance', {})
    _fail(not isinstance(provenance, dict), 'provenance mapping required')
    for key in ('run_id', 'exporter_revision'):
        _fail(not isinstance(provenance.get(key), str) or not provenance[key].strip(), f'provenance.{key} required')
    hashes = provenance.get('source_hashes')
    _fail(not isinstance(hashes, dict) or not hashes, 'nonempty provenance.source_hashes required')
    for key, value in hashes.items():
        _fail(not isinstance(key, str) or not isinstance(value, str) or len(value) != 64
              or any(c not in '0123456789abcdef' for c in value.lower()), 'source hashes must be named SHA256 hex strings')
    for label, ids, source_count in (('k', data.k_ids, meta.get('source_nkpoints')),
                                     ('band', data.band_ids, meta.get('source_nbands'))):
        _fail(not isinstance(source_count, int) or isinstance(source_count, bool) or source_count <= 0,
              f'positive integer source_{"nkpoints" if label == "k" else "nbands"} required')
        _fail(not isinstance(ids, np.ndarray) or ids.ndim != 1 or ids.dtype.kind != 'i' or ids.dtype.itemsize != 8
              or not len(ids), f'{label}_ids must be nonempty signed int64 array')
        _fail(np.any(ids < 1) or np.any(ids > source_count) or len(np.unique(ids)) != len(ids),
              f'duplicate or out-of-range global 1-based {label} IDs')
    nk, nb = len(data.k_ids), len(data.band_ids)
    shapes = {'kpoints_fractional': (nk, 3), 'weights': (nk,), 'energies_eV': (nk, nb),
              'lattice_A': (3, 3), 'reciprocal_inv_A': (3, 3)}
    for name, shape in shapes.items():
        a = getattr(data, name)
        _fail(not isinstance(a, np.ndarray) or a.shape != shape or a.dtype.kind != 'f' or a.dtype.itemsize != 8,
              f'{name} must be float64 with shape {shape}')
        _fail(not np.isfinite(a).all(), f'{name} contains nonfinite/missing values')
    _fail(np.any(data.weights < 0) or not np.isclose(data.weights.sum(), 1, atol=1e-12, rtol=0),
          'k weights must be nonnegative and sum to one; no automatic normalization')
    _fail(np.linalg.det(data.lattice_A) <= 0 or not np.allclose(data.lattice_A @ data.reciprocal_inv_A.T,
          2*np.pi*np.eye(3), atol=1e-8, rtol=1e-10), 'inconsistent or non-right-handed lattice/2pi reciprocal pair')
    shape = (nk, 3, nb, nb)
    _fail(not isinstance(data.D_eVA, np.ndarray) or data.D_eVA.shape != shape or data.D_eVA.dtype.kind != 'c'
          or data.D_eVA.dtype.itemsize != 16, f'D_eVA must be complex128 with complete band dimensions {shape}')
    _fail(not isinstance(data.coverage, np.ndarray) or data.coverage.shape != shape or data.coverage.dtype != np.bool_,
          'explicit boolean coverage with the full matrix shape is required')
    _fail(not np.isfinite(data.D_eVA[data.coverage]).all(), 'covered matrix elements must be finite')
    unavailable = data.D_eVA[~data.coverage]
    _fail(not (np.isnan(unavailable.real) & np.isnan(unavailable.imag)).all(),
          'unavailable elements must be explicit complex NaN, not zero or infinity')
    _fail(not np.array_equal(data.coverage, data.coverage.swapaxes(-1, -2)),
          'coverage must include both bra/ket directions for each Hermitian matrix element')
    diagonal = np.diagonal(data.coverage, axis1=-2, axis2=-1)
    diagonal_status = meta.get('diagonal_status')
    _fail(diagonal_status not in ('present', 'unavailable', 'coverage_mask'), 'explicit diagonal_status required')
    _fail(diagonal_status == 'present' and not diagonal.all(), 'diagonal_status=present but diagonal data missing')
    _fail(diagonal_status == 'unavailable' and diagonal.any(), 'diagonal_status=unavailable but covered diagonal data present')
    _fail(not np.isfinite([hermitian_atol_eVA, hermitian_rtol]).all()
          or min(hermitian_atol_eVA, hermitian_rtol) < 0, 'finite nonnegative Hermiticity tolerances required')
    adjoint = data.D_eVA.swapaxes(-1, -2).conj()
    error = np.abs(data.D_eVA[data.coverage]-adjoint[data.coverage])
    scale = np.maximum(np.abs(data.D_eVA[data.coverage]), np.abs(adjoint[data.coverage]))
    _fail(np.any(error > hermitian_atol_eVA+hermitian_rtol*scale),
          'matrix is not Hermitian within declared numerical tolerances; check phases, energy convention and missing terms')
    return {'max_hermitian_residual_eVA': float(error.max(initial=0)),
            'hermitian_atol_eVA': float(hermitian_atol_eVA), 'hermitian_rtol': float(hermitian_rtol),
            'covered_elements': int(data.coverage.sum()), 'missing_elements': int((~data.coverage).sum()),
            'full_source_band_coverage': bool(nb == meta['source_nbands']),
            'full_source_k_coverage': bool(nk == meta['source_nkpoints']),
            'operator_kind': operator['kind'], 'operator_accuracy_status': status}


def read_matrix_bundle(npz_path, metadata_path, *, allow_experimental=False,
                       hermitian_atol_eVA=1e-10, hermitian_rtol=1e-8):
    """Read only a complete NPZ whose exact bytes are bound by the JSON hash."""
    npz_path, metadata_path = Path(npz_path), Path(metadata_path)
    metadata = json.loads(metadata_path.read_text())
    _fail(not isinstance(metadata, dict), 'matrix metadata must be a JSON object')
    _fail(metadata.get('matrix_npz_sha256') != sha256(npz_path), 'matrix NPZ SHA256 mismatch or missing checksum')
    with np.load(npz_path, allow_pickle=False) as archive:
        _fail(any(key not in archive.files for key in ARRAYS), 'matrix NPZ is missing required arrays')
        data = MatrixDataset(metadata, **{key: archive[key] for key in ARRAYS})
    validate_matrix_data(data, allow_experimental=allow_experimental,
                         hermitian_atol_eVA=hermitian_atol_eVA, hermitian_rtol=hermitian_rtol)
    return data


def _select_ids(requested, available, label):
    a = np.asarray(requested)
    _fail(a.ndim != 1 or not a.size or a.dtype.kind not in 'iu' or len(np.unique(a)) != len(a),
          f'{label} must be a nonempty, duplicate-free list of integer global band IDs')
    locations = {int(value): index for index, value in enumerate(available)}
    _fail(any(int(value) not in locations for value in a), f'{label} contains bands absent from the export')
    return a.astype(np.int64), np.array([locations[int(value)] for value in a], dtype=int)


def berry_curvature(data: MatrixDataset, *, interest_band_ids, intermediate_band_ids,
                    degeneracy_threshold_eV, degeneracy_policy='error', allow_experimental=False,
                    hermitian_atol_eVA=1e-10, hermitian_rtol=1e-8) -> BerryResult:
    """Compute Omega=-2 Im sum_m D_a(n,m)D_b(m,n)/(En-Em)^2.

    Output xyz uses (a,b)=(y,z),(z,x),(x,y), consistent with A=i<u|grad u>.
    n and m windows are mandatory and independent. A small-gap pair never
    disappears silently: error is default, mask returns NaN for the affected
    individual-band curvature and preserves the explicit pair list. There is
    no eta regulator, lifetime broadening, or degenerate-subspace algorithm.
    """
    diagnostics = validate_matrix_data(data, allow_experimental=allow_experimental,
                                       hermitian_atol_eVA=hermitian_atol_eVA, hermitian_rtol=hermitian_rtol)
    _fail(not np.isfinite(degeneracy_threshold_eV) or degeneracy_threshold_eV < 0,
          'explicit finite nonnegative degeneracy threshold in eV required')
    _fail(degeneracy_policy not in ('error', 'mask'), 'degeneracy_policy must be error or mask')
    n_ids, n_indices = _select_ids(interest_band_ids, data.band_ids, 'interest n window')
    m_ids, m_indices = _select_ids(intermediate_band_ids, data.band_ids, 'intermediate m window')
    nk = len(data.k_ids)
    omega = np.zeros((nk, len(n_ids), 3), dtype=float)
    gaps = np.empty((nk, len(n_ids)), dtype=float)
    valid = np.ones((nk, len(n_ids)), dtype=bool)
    pairs = []
    for j, (n_id, n) in enumerate(zip(n_ids, n_indices)):
        others = m_indices[m_ids != n_id]
        _fail(not len(others), f'interest band {n_id}: m window has no other band')
        delta = data.energies_eV[:, n, None]-data.energies_eV[:, others]
        # Truncating the intermediate sum must never hide an exported band
        # that makes the individual state degenerate. Unknown omitted source
        # energies remain a separate limitation, recorded below.
        known_others = np.flatnonzero(data.band_ids != n_id)
        known_delta = data.energies_eV[:, n, None]-data.energies_eV[:, known_others]
        gaps[:, j] = np.min(abs(known_delta), axis=1)
        known_near = abs(known_delta) <= degeneracy_threshold_eV
        near = abs(delta) <= degeneracy_threshold_eV
        for k, other in np.argwhere(known_near):
            pairs.append({'k_id': int(data.k_ids[k]), 'n_band': int(n_id),
                          'm_band': int(data.band_ids[known_others[other]]), 'gap_eV': float(abs(known_delta[k, other]))})
        if known_near.any() and degeneracy_policy == 'error':
            first = pairs[-int(known_near.sum())]
            raise DegenerateBandError(f'band curvature undefined within threshold {degeneracy_threshold_eV:g} eV: '
                                      f'k={first["k_id"]}, n={first["n_band"]}, m={first["m_band"]}, '
                                      f'gap={first["gap_eV"]:.12g} eV; use an appropriate complete-subspace method or explicit mask')
        valid[:, j] = ~known_near.any(axis=1)
        needed = np.take(data.coverage[:, :, n, :], others, axis=-1)
        # Explicitly masked states need no reconstructed missing matrix data;
        # every vertex used for a valid state is still mandatory.
        missing = ~needed & valid[:, j, None, None]
        if missing.any():
            k, component, other = np.argwhere(missing)[0]
            raise MatrixContractError(f'missing required element k={data.k_ids[k]} axis={"xyz"[component]} '
                                      f'n={n_id} m={data.band_ids[others[other]]}; regenerate complete interband data')
        # Both directions are stored and checked; do not silently reconstruct
        # the reverse vertex, renormalize states, or drop imaginary components.
        forward = np.take(data.D_eVA[:, :, n, :], others, axis=-1)
        reverse = np.take(data.D_eVA[:, :, :, n], others, axis=-1)
        safe_denominator = np.where(near, np.nan, delta**2)
        for component, (a, b) in enumerate(((1, 2), (2, 0), (0, 1))):
            terms = -2*np.imag(forward[:, a, :]*reverse[:, b, :])/safe_denominator
            omega[:, j, component] = terms.sum(axis=1)
        omega[~valid[:, j], j] = np.nan
    finite_values = omega[valid]
    _fail(not np.isfinite(finite_values).all(), 'overflow/nonfinite curvature outside declared degenerate bands')
    diagnostics.update({'formula': '-2 Im sum_m D_a(n,m) D_b(m,n)/(En-Em)^2',
                        'berry_connection_convention': 'A=i<u|grad_k u>',
                        'output_components': ['Omega_yz', 'Omega_zx', 'Omega_xy'],
                        'curvature_units': 'Angstrom^2', 'degeneracy_threshold_eV': float(degeneracy_threshold_eV),
                        'degeneracy_policy': degeneracy_policy, 'degenerate_pairs': pairs,
                        'invalid_individual_band_count': int((~valid).sum()),
                        'interest_band_ids': n_ids.tolist(), 'intermediate_band_ids': m_ids.tolist(),
                        'intermediate_contains_all_source_bands': bool(len(m_ids) == data.metadata['source_nbands']),
                        'isolation_check_scope': 'all exported band energies, independently of requested m window',
                        'all_source_band_energies_available': bool(len(data.band_ids) == data.metadata['source_nbands']),
                        'physical_accuracy_claim': 'None added: numerical checks preserve the input operator accuracy label.'})
    return BerryResult(data.k_ids.copy(), n_ids, m_ids, omega, gaps, valid, diagnostics)


def _band_range(text):
    result = []
    for block in text.split(','):
        ends = block.split(':')
        if len(ends) == 1:
            result.append(int(ends[0]))
        elif len(ends) == 2 and int(ends[0]) <= int(ends[1]):
            result.extend(range(int(ends[0]), int(ends[1])+1))
        else:
            raise argparse.ArgumentTypeError('use 1:40 or 1,3,5:8 for global band IDs')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrices', required=True, type=Path)
    parser.add_argument('--metadata', required=True, type=Path)
    parser.add_argument('--n-bands', required=True, type=_band_range)
    parser.add_argument('--m-bands', required=True, type=_band_range)
    parser.add_argument('--degeneracy-threshold-eV', required=True, type=float)
    parser.add_argument('--degeneracy-policy', choices=['error', 'mask'], default='error')
    parser.add_argument('--allow-experimental', action='store_true')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error('output directory exists; preserve prior results with a fresh path')
    data = read_matrix_bundle(args.matrices, args.metadata, allow_experimental=args.allow_experimental)
    result = berry_curvature(data, interest_band_ids=args.n_bands, intermediate_band_ids=args.m_bands,
                             degeneracy_threshold_eV=args.degeneracy_threshold_eV,
                             degeneracy_policy=args.degeneracy_policy, allow_experimental=args.allow_experimental)
    args.output_dir.mkdir(parents=True)
    _, indices = _select_ids(result.band_ids, data.band_ids, 'output')
    np.savez_compressed(args.output_dir/'curvature.npz', k_ids=result.k_ids, band_ids=result.band_ids,
                        intermediate_band_ids=result.intermediate_band_ids,
                        kpoints_fractional=data.kpoints_fractional, weights=data.weights,
                        energies_eV=data.energies_eV[:, indices], omega_A2=result.omega_A2,
                        min_gap_eV=result.min_gap_eV, valid_nondegenerate=result.valid_nondegenerate,
                        lattice_A=data.lattice_A, reciprocal_inv_A=data.reciprocal_inv_A)
    with (args.output_dir/'curvature.csv').open('w', newline='') as handle:
        fields = ['k_index', 'band', 'kx_frac', 'ky_frac', 'kz_frac', 'energy_eV', 'k_weight',
                  'omega_x_A2', 'omega_y_A2', 'omega_z_A2', 'min_gap_eV', 'valid_nondegenerate']
        writer = csv.writer(handle); writer.writerow(fields)
        for k, j in np.ndindex(result.min_gap_eV.shape):
            writer.writerow([result.k_ids[k], result.band_ids[j], *data.kpoints_fractional[k],
                             data.energies_eV[k, indices[j]], data.weights[k],
                             *[float(x) if np.isfinite(x) else '' for x in result.omega_A2[k, j]],
                             result.min_gap_eV[k, j], bool(result.valid_nondegenerate[k, j])])
    record = dict(result.diagnostics)
    record.update({'source_matrix_npz': {'path': str(args.matrices.resolve()), 'sha256': sha256(args.matrices)},
                   'source_metadata': {'path': str(args.metadata.resolve()), 'sha256': sha256(args.metadata)},
                   'source_operator': data.metadata['operator'], 'source_provenance': data.metadata['provenance']})
    (args.output_dir/'diagnostics.json').write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'output': str(args.output_dir), 'nkpoints': len(result.k_ids), 'n_interest_bands': len(result.band_ids),
                      'operator_accuracy_status': record['operator_accuracy_status'],
                      'invalid_individual_band_count': record['invalid_individual_band_count']}, indent=2))


if __name__ == '__main__':
    main()
