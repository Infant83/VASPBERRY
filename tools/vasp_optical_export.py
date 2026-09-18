#!/usr/bin/env python3
"""Read the experimental VASP 5.4.4 optical-connection stream without repairing it.

Native format v1 is little endian, with explicit 32-byte ASCII sentinels and
Fortran-order arrays. Public arrays retain a leading spin dimension. C_A has
axes (spin, k, Cartesian, bra, ket), energies have (spin, k, band), and direct
and reciprocal lattice vectors are rows. D=i*(E_ket-E_bra)*C; its diagonal is
supplied independently by ENERGY_DER. Off-diagonal elements zeroed by the
producer's degenerate clusters or denominator guard are explicit complex NaNs
in the interchange bundle. Raw C and D
remain unchanged in raw-connection.npz, including any Hermiticity mismatch.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import struct

import numpy as np

try:
    from .exported_matrix_kubo import SCHEMA, SCHEMA_VERSION, UNITS, sha256
except ImportError:
    from exported_matrix_kubo import SCHEMA, SCHEMA_VERSION, UNITS, sha256


MAGIC = b'VBERRY_PAW_CONNECTION_V1'.ljust(32, b' ')
FOOTER = b'VBERRY_CONNECTION_COMPLETE_V1'.ljust(32, b' ')
ZERO_GAP_THRESHOLD_EV = 1e-10
# Native v1's current VASP 5.4.4 producer uses default-real 1E-10 (not
# 1E-10_q) in EDDIAG_LR's strict denominator comparison. The audited build
# promotes that float32 literal to q; its cluster threshold is 1E-10_q.
PRODUCER_DENOMINATOR_ZERO_THRESHOLD_EV = float(np.float32(1e-10))
REVISION = 'vasp-optical-stream-adapter-v2'
COVERAGE_RULE = ('Full native band-order adjacent gaps <= producer threshold form transitive '
                 'clusters independently per spin/k; off-diagonal entries require different '
                 'clusters, gap > producer threshold and gap >= promoted denominator cutoff; '
                 'diagonal ENERGY_DER entries remain available. Slice this full mask with the bands.')


class OpticalExportError(ValueError):
    """Malformed, unsupported or inconsistent native export."""


def _require(condition, message):
    if not condition:
        raise OpticalExportError(message)


def _finite(name, value):
    _require(np.isfinite(value).all(), f'{name} contains nonfinite values')


@dataclass(frozen=True)
class ConnectionData:
    path: Path
    source_sha256: str
    degeneracy_threshold_eV: float
    fermi_eV: float
    lattice_A: np.ndarray
    reciprocal_inv_A: np.ndarray
    kpoints_fractional: np.ndarray
    weights: np.ndarray
    energies_eV: np.ndarray             # spin, k, band
    occupations: np.ndarray             # spin, k, band; native FERTOT, no rescaling
    energy_der_eVA: np.ndarray           # spin, k, Cartesian, band
    C_A: np.ndarray                     # spin, k, Cartesian, bra, ket

    @property
    def nspin(self):
        return self.energies_eV.shape[0]

    @property
    def nkpoints(self):
        return self.energies_eV.shape[1]

    @property
    def nbands(self):
        return self.energies_eV.shape[2]

    @property
    def raw_D_eVA(self):
        """Algebraic i*(E_ket-E_bra)*C, retaining computed entries unchanged."""
        e = self.energies_eV
        return 1j*(e[:, :, None, :]-e[:, :, :, None])[:, :, None]*self.C_A

    @property
    def coverage(self):
        """Reproduce native-v1 zeroing on the complete source band axis.

        FIND_DEG_CLUSTERS joins adjacent band indices transitively without
        sorting. Recomputing clusters after a band selection could wrongly
        restore an outer pair whose connecting middle band was omitted.
        """
        e = self.energies_eV
        labels = np.zeros(e.shape, dtype=np.int64)
        labels[:, :, 1:] = np.cumsum(abs(np.diff(e, axis=-1)) > self.degeneracy_threshold_eV, axis=-1)
        gap = abs(e[:, :, None, :]-e[:, :, :, None])
        mask = ((labels[:, :, None, :] != labels[:, :, :, None])
                & (gap > self.degeneracy_threshold_eV)
                & (gap >= PRODUCER_DENOMINATOR_ZERO_THRESHOLD_EV))
        diagonal = np.arange(self.nbands)
        mask[:, :, diagonal, diagonal] = True
        return np.broadcast_to(mask[:, :, None], self.C_A.shape).copy()

    @property
    def D_eVA(self):
        d = self.raw_D_eVA
        diagonal = np.arange(self.nbands)
        d[:, :, :, diagonal, diagonal] = self.energy_der_eVA
        d[~self.coverage] = complex(np.nan, np.nan)
        return d


@dataclass(frozen=True)
class WaveDerData:
    path: Path
    source_sha256: str
    dielectric_node: float
    plasmon: np.ndarray
    C_A: np.ndarray                     # spin, k, Cartesian, bra, ket; complex64


def read_connection(path) -> ConnectionData:
    """Strict native-v1 reader; does not symmetrize, mirror, or normalize data."""
    path = Path(path)
    blob = path.read_bytes()
    _require(len(blob) >= 108, 'truncated optical-connection header/footer')
    _require(blob[:32] == MAGIC, 'invalid optical-connection magic or padding')
    version, nk, nb, ns, nc, real_bytes, complex_bytes = struct.unpack_from('<7i', blob, 32)
    _require(version == 1, f'unsupported connection version {version}; little-endian v1 required')
    _require(nk > 0 and nb > 0 and ns in (1, 2), 'invalid nkpoints/nbands/spin header')
    _require((nc, real_bytes, complex_bytes) == (3, 8, 16), 'unsupported components or scalar precision')
    # Fixed header 76, two 3x3 lattices, k coordinates/weights, energy/occupation,
    # diagonal derivatives, full complex C, and 32-byte completion footer.
    expected = 108 + 144 + 32*nk + 40*nb*nk*ns + 48*nb*nb*nk*ns
    _require(len(blob) == expected, f'connection byte length {len(blob)} != expected {expected}')
    _require(blob[-32:] == FOOTER, 'missing or damaged connection completion footer')
    threshold, fermi = struct.unpack_from('<2d', blob, 60)
    _finite('threshold/Fermi energy', np.array([threshold, fermi]))
    _require(threshold >= 0, 'negative producer degeneracy threshold')
    _require(threshold == ZERO_GAP_THRESHOLD_EV,
             'v1 producer degeneracy threshold must equal 1e-10 eV; other cluster-zeroing policies are unsupported')
    offset = 76

    def take(dtype, shape):
        nonlocal offset
        count = int(np.prod(shape))
        result = np.frombuffer(blob, dtype=dtype, count=count, offset=offset).reshape(shape, order='F')
        offset += count*np.dtype(dtype).itemsize
        return result

    lattice = take('<f8', (3, 3)).T.copy()
    reciprocal = take('<f8', (3, 3)).T.copy()
    kpoints = take('<f8', (3, nk)).T.copy()
    weights = take('<f8', (nk,)).copy()
    energies = take('<f8', (nb, nk, ns)).transpose(2, 1, 0).copy()
    occupations = take('<f8', (nb, nk, ns)).transpose(2, 1, 0).copy()
    derivatives = take('<f8', (nb, nk, ns, 3)).transpose(2, 1, 3, 0).copy()
    c = take('<c16', (nb, nb, nk, ns, 3)).transpose(3, 2, 4, 0, 1).copy()
    _require(offset == len(blob)-32, 'internal connection array-length mismatch')
    for name, value in (('lattice', lattice), ('reciprocal', reciprocal), ('kpoints', kpoints),
                        ('weights', weights), ('energies', energies), ('occupations', occupations),
                        ('energy derivatives', derivatives), ('connection', c)):
        _finite(name, value)
    _require(np.linalg.det(lattice) > 0 and np.allclose(lattice @ reciprocal.T,
             2*np.pi*np.eye(3), atol=1e-8, rtol=1e-10), 'inconsistent/right-handed lattice and reciprocal pair')
    _require((weights >= 0).all() and np.isclose(weights.sum(), 1, atol=1e-12, rtol=0),
             'k weights must be nonnegative and sum to one; no normalization performed')
    data = ConnectionData(path.resolve(), hashlib.sha256(blob).hexdigest(), threshold, fermi,
                          lattice, reciprocal, kpoints, weights, energies, occupations, derivatives, c)
    _finite('derived connection velocity', data.raw_D_eVA)
    return data


def read_waveder(path) -> WaveDerData:
    """Read native complex64 WAVEDER: four records with 4-byte LE markers.

    Split/continuation records and gamma-real variants are explicitly unsupported.
    Both stored triangles are retained; the comparison helper uses only the lower.
    """
    path = Path(path)
    blob = path.read_bytes()
    offset = 0

    def record(expected=None):
        nonlocal offset
        _require(offset+8 <= len(blob), 'truncated WAVEDER record marker')
        length = struct.unpack_from('<i', blob, offset)[0]
        _require(length >= 0, 'split/negative WAVEDER record marker unsupported')
        _require(expected is None or length == expected,
                 f'WAVEDER record length {length} != expected {expected}')
        _require(offset+8+length <= len(blob), 'truncated WAVEDER record payload')
        trailer = struct.unpack_from('<i', blob, offset+4+length)[0]
        _require(trailer == length, 'WAVEDER leading/trailing record markers disagree')
        payload = memoryview(blob)[offset+4:offset+4+length]
        offset += length+8
        return payload

    nb, nd, nk, ns = struct.unpack('<4i', record(16))
    _require(nb > 0 and 0 < nd <= nb and nk > 0 and ns in (1, 2), 'invalid WAVEDER dimensions')
    node = struct.unpack('<d', record(8))[0]
    plasmon = np.frombuffer(record(72), '<f8').reshape((3, 3), order='F').copy()
    payload = record(8*nb*nd*nk*ns*3)
    c = np.frombuffer(payload, '<c8').reshape((nb, nd, nk, ns, 3), order='F').transpose(3, 2, 4, 0, 1).copy()
    _require(offset == len(blob), 'trailing bytes/extra records after WAVEDER matrix')
    _finite('WAVEDER scalar metadata', np.array([node]))
    _finite('WAVEDER plasmon', plasmon)
    _finite('WAVEDER connection', c)
    return WaveDerData(path.resolve(), hashlib.sha256(blob).hexdigest(), node, plasmon, c)


def hermiticity_stats(matrix, coverage=None):
    """Report residuals without changing either matrix triangle."""
    matrix = np.asarray(matrix)
    _require(matrix.ndim >= 2 and matrix.shape[-1] == matrix.shape[-2], 'square matrix axes required')
    mask = np.isfinite(matrix) if coverage is None else np.asarray(coverage, dtype=bool)
    _require(mask.shape == matrix.shape, 'Hermiticity mask shape mismatch')
    mask = mask & mask.swapaxes(-1, -2)
    lhs = matrix[mask]
    rhs = matrix.swapaxes(-1, -2).conj()[mask]
    _finite('covered Hermiticity entries', lhs)
    error = abs(lhs-rhs)
    denominator = np.linalg.norm(lhs)
    diagonal = np.diagonal(matrix, axis1=-2, axis2=-1)
    finite_diagonal = diagonal[np.isfinite(diagonal)]
    return {'compared_elements': int(mask.sum()), 'max_absolute_residual': float(error.max(initial=0)),
            'rms_absolute_residual': float(np.sqrt(np.mean(error**2))) if error.size else 0.0,
            'relative_frobenius_residual': float(np.linalg.norm(lhs-rhs)/denominator) if denominator else 0.0,
            'maximum_diagonal_imaginary': float(abs(finite_diagonal.imag).max(initial=0)),
            'passed_default_kubo_tolerance': bool(np.all(error <= 1e-10+1e-8*np.maximum(abs(lhs), abs(rhs)))),
            'repair_applied': False}


def compare_waveder_lower(connection: ConnectionData, waveder: WaveDerData):
    source = waveder.C_A
    ns, nk, nc, nb, nd = source.shape
    _require((ns, nk, nc, nb) == (connection.nspin, connection.nkpoints, 3, connection.nbands),
             'WAVEDER and connection dimensions disagree')
    mask = np.arange(nb)[:, None] > np.arange(nd)[None, :]
    raw = connection.C_A[..., :nd]
    error = abs(raw[..., mask]-source[..., mask])
    denominator = np.linalg.norm(source[..., mask])
    return {'comparison': 'strict lower triangle (bra_band > ket_band); no mirroring',
            'connection_units': 'Angstrom', 'waveder_precision': 'complex64',
            'native_shape': list(source.shape), 'compared_elements': int(error.size),
            'max_absolute_difference_A': float(error.max(initial=0)),
            'rms_absolute_difference_A': float(np.sqrt(np.mean(error**2))) if error.size else 0.0,
            'relative_frobenius_difference': (float(np.linalg.norm(error)/denominator) if denominator
                                                else (0.0 if not error.any() else None)),
            'waveder_sha256': waveder.source_sha256}


def _provenance(run_dir, data):
    path = run_dir/'run.json'
    _require(path.is_file(), 'run.json is required for a CLI conversion with provenance')
    record = json.loads(path.read_text())
    _require(isinstance(record, dict), 'run.json must contain an object')
    _require(record.get('completed_export') is True,
             'run.json must confirm completed_export=true before CLI conversion')
    _require(type(record.get('returncode')) is int and record['returncode'] == 0
             and record.get('normal_timing_footer') is True and record.get('timed_out') is not True,
             'run.json must confirm returncode=0 and normal timing footer without timeout')
    outputs = record.get('outputs', {})
    _require(isinstance(outputs, dict), 'run.json outputs must be a mapping')
    if 'BERRY_CONNECTION.bin' in outputs:
        recorded = outputs['BERRY_CONNECTION.bin']
        _require(isinstance(recorded, dict) and recorded.get('sha256') == data.source_sha256,
                 'run.json BERRY_CONNECTION.bin output SHA256 mismatch')
        if 'bytes' in recorded:
            _require(recorded['bytes'] == data.path.stat().st_size, 'run.json connection output byte count mismatch')
    output_hashes = record.get('output_sha256', {})
    _require(isinstance(output_hashes, dict), 'run.json output_sha256 must be a mapping')
    _require('BERRY_CONNECTION.bin' in outputs or 'BERRY_CONNECTION.bin' in output_hashes,
             'run.json must bind BERRY_CONNECTION.bin to its producer-recorded SHA256')
    if 'BERRY_CONNECTION.bin' in output_hashes:
        _require(output_hashes['BERRY_CONNECTION.bin'] == data.source_sha256,
                 'run.json BERRY_CONNECTION.bin output SHA256 mismatch')
    hashes = {'BERRY_CONNECTION.bin': data.source_sha256, 'run.json': sha256(path),
              'adapter_source': sha256(Path(__file__))}
    for name in ('binary_sha256', 'patch_sha256', 'modified_source_sha256'):
        if name in record:
            hashes[name] = record[name]
    for group in ('input_sha256', 'source_hashes', 'source_input_sha256'):
        value = record.get(group, {})
        _require(isinstance(value, dict), f'run.json {group} must be a mapping')
        for name, digest in value.items():
            hashes[f'{group}/{name}'] = digest
    for name, digest in hashes.items():
        _require(isinstance(digest, str) and len(digest) == 64 and all(c in '0123456789abcdef' for c in digest.lower()),
                 f'invalid SHA256 in provenance: {name}')
    for field, actual in (('nkpoints', data.nkpoints), ('nbands', data.nbands)):
        if field in record:
            _require(record[field] == actual, f'run.json {field} disagrees with exported header')
    return {'run_id': str(record.get('run_id') or run_dir.name), 'exporter_revision': REVISION,
            'source_hashes': hashes, 'run_directory': str(run_dir.resolve()), 'run_record': record}


def write_bundles(data: ConnectionData, output_dir, provenance, *, waveder=None):
    """Write a fresh diagnostic export even when Hermiticity residuals fail.

    The strict downstream Kubo reader still rejects failed Hermiticity. Metadata
    describes the intended Hermitian operator and preserves experimental status;
    it does not certify that every measured entry satisfies that definition.
    """
    output = Path(output_dir)
    _require(not output.exists(), 'output directory must be new; existing results are never overwritten')
    d = data.D_eVA
    coverage = data.coverage
    report = {'schema': 'vaspberry.optical-stream-conversion', 'version': 1,
              'source_connection_sha256': data.source_sha256, 'source_shape': list(data.C_A.shape),
              'producer_degeneracy_threshold_eV': data.degeneracy_threshold_eV,
              'unavailable_offdiagonal_gap_threshold_eV': ZERO_GAP_THRESHOLD_EV,
              'producer_denominator_zero_threshold_eV': PRODUCER_DENOMINATOR_ZERO_THRESHOLD_EV,
              'coverage_rule': COVERAGE_RULE,
              'C_A_hermiticity': hermiticity_stats(data.C_A),
              'D_eVA_hermiticity': hermiticity_stats(d, coverage),
              'missing_offdiagonal_elements': int((~coverage).sum()),
              'physical_accuracy_status': 'experimental', 'repair_applied': False, 'bundles': []}
    if waveder is not None:
        report['native_waveder_lower_triangle'] = compare_waveder_lower(data, waveder)
    output.mkdir(parents=True)
    raw_path = output/'raw-connection.npz'
    np.savez_compressed(raw_path, C_A=data.C_A, raw_D_eVA=data.raw_D_eVA,
                        energy_der_eVA=data.energy_der_eVA, energies_eV=data.energies_eV,
                        occupations=data.occupations, coverage=coverage,
                        kpoints_fractional=data.kpoints_fractional, weights=data.weights,
                        lattice_A=data.lattice_A, reciprocal_inv_A=data.reciprocal_inv_A)
    report['raw_npz_sha256'] = sha256(raw_path)
    for spin in range(data.nspin):
        npz_path = output/f'matrix-spin{spin+1}.npz'
        json_path = npz_path.with_suffix('.json')
        np.savez_compressed(npz_path,
                            k_ids=np.arange(1, data.nkpoints+1, dtype=np.int64),
                            band_ids=np.arange(1, data.nbands+1, dtype=np.int64),
                            kpoints_fractional=data.kpoints_fractional, weights=data.weights,
                            energies_eV=data.energies_eV[spin], D_eVA=d[spin], coverage=coverage[spin],
                            lattice_A=data.lattice_A, reciprocal_inv_A=data.reciprocal_inv_A,
                            occupations=data.occupations[spin])
        metadata = {'schema': SCHEMA, 'version': SCHEMA_VERSION, 'complete': True,
                    'source_nkpoints': data.nkpoints, 'source_nbands': data.nbands,
                    'source_nspin': data.nspin, 'spin_channel_1based': spin+1,
                    'units': UNITS.copy(), 'matrix_axes': ['k', 'cartesian', 'bra_band', 'ket_band'],
                    'matrix_element_convention': '<n|D_a|m>', 'cartesian_components': ['x', 'y', 'z'],
                    'reciprocal_convention': '2pi', 'weights_convention': 'sum_one',
                    'diagonal_status': 'present', 'fermi_eV': data.fermi_eV,
                    'producer_degeneracy_threshold_eV': data.degeneracy_threshold_eV,
                    'missing_offdiagonal_gap_threshold_eV': ZERO_GAP_THRESHOLD_EV,
                    'producer_denominator_zero_threshold_eV': PRODUCER_DENOMINATOR_ZERO_THRESHOLD_EV,
                    'coverage_rule': COVERAGE_RULE,
                    'operator': {'kind': 'derived_from_vasp_paw_optical_connection',
                                 'definition': 'D_mn=i*(E_n-E_m)*C_mn; diagonal=ENERGY_DER, from unmirrored VASP response.',
                                 'hermitian': True, 'accuracy_status': 'experimental',
                                 'included_terms': ['VASP analytic optical-response PAW connection and independent ENERGY_DER'],
                                 'missing_terms': ['undefined producer-zeroed degenerate-cluster or small-gap off-diagonal reconstruction'],
                                 'validation_evidence': ['binary completeness/finite values checked; measured Hermiticity preserved',
                                                         'full physical operator and material convergence not established']},
                    'hermiticity_diagnostics': hermiticity_stats(d[spin], coverage[spin]),
                    'provenance': provenance, 'matrix_npz_sha256': sha256(npz_path),
                    'native_connection_convention': 'C_mn=<u_m|-i partial_k|u_n>',
                    'occupation_note': 'Native FERTOT retained separately; no spin multiplicity or occupation weighting applied.',
                    'repair_applied': False}
        json_path.write_text(json.dumps(metadata, indent=2, allow_nan=False)+'\n')
        report['bundles'].append({'spin_channel_1based': spin+1, 'npz': str(npz_path), 'metadata': str(json_path)})
    (output/'conversion-report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        run = args.run_dir.resolve()
        data = read_connection(run/'BERRY_CONNECTION.bin')
        provenance = _provenance(run, data)
        waveder = read_waveder(run/'WAVEDER') if (run/'WAVEDER').is_file() else None
        report = write_bundles(data, args.output_dir, provenance, waveder=waveder)
    except (OpticalExportError, OSError, json.JSONDecodeError) as error:
        parser.exit(2, f'error: {error}\n')
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
