"""Spin operators in the unchanged WAVECAR eigenvector gauge.

``spin_pauli[k,a,n,m]`` is <u_n|sigma_a|u_m>, dimensionless; physical
spin is (hbar/2) times this matrix. WAVECAR coefficients alone yield the
pseudo operator and pseudo overlap, without normalization. A separately
audited, matching PAW correction is required for a physical PAW operator.
The WAVECAR does not encode SAXIS: callers must declare the spin basis.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import re
import zipfile

import numpy as np

from berry_data import require
from exported_matrix_kubo import sha256
from wavecar_fukui import Wavecar

SCHEMA = 'vaspberry.spin-operators'
AUGMENTATION_SCHEMA = 'vaspberry.spin-augmentation'
UNITS = {'spin': 'dimensionless Pauli', 'overlap': 'dimensionless'}
GAUGE = 'same_eigenvectors_as_source_wavecar'
MAX_ARRAY_BYTES = 2 * 1024**3
HERMITICITY_ATOL = 1e-10
PHYSICAL_METRIC_ATOL = 5e-5
ARRAYS = ('kpoints_fractional', 'energies_eV', 'lattice_A', 'band_indices',
          'spin_pauli', 'overlap', 'pseudo_spin_pauli', 'pseudo_overlap')
AUGMENTATION_ARRAYS = ('kpoints_fractional', 'energies_eV', 'lattice_A',
                       'band_indices', 'delta_spin_pauli', 'delta_overlap')


@dataclass(frozen=True)
class SpinOperators:
    kpoints_fractional: np.ndarray
    energies_eV: np.ndarray
    lattice_A: np.ndarray
    band_indices: np.ndarray
    spin_pauli: np.ndarray
    overlap: np.ndarray
    pseudo_spin_pauli: np.ndarray
    pseudo_overlap: np.ndarray
    metadata: dict


def _finite(array, shape, dtype, name):
    require(isinstance(array, np.ndarray) and array.dtype == dtype
            and array.shape == shape and np.isfinite(array).all(),
            f'{name} must be finite {np.dtype(dtype)} with shape {shape}')


def _hermitian(array, name):
    error = float(np.max(abs(array-array.conj().swapaxes(-1, -2))))
    require(error <= HERMITICITY_ATOL, name+' is not Hermitian; no repair is applied')
    return error


def _hash(value):
    return isinstance(value, str) and re.fullmatch('[0-9a-f]{64}', value) is not None


def _basis(value):
    require(isinstance(value, str) and bool(value.strip()), 'explicit nonempty spin_basis required; SAXIS cannot be inferred from WAVECAR')


def pauli_from_coefficients(coefficients):
    """Return raw (Pauli[3,N,N], overlap[N,N]) from [N,2,G] coefficients.

    No band normalization, orthogonalization, or PAW correction is applied.
    The first spinor block is up and the second down in the declared basis.
    """
    c = np.asarray(coefficients, dtype=np.complex128)
    require(c.ndim == 3 and c.shape[0] > 0 and c.shape[1] == 2
            and c.shape[2] > 0 and np.isfinite(c).all(),
            'finite [band,2,G] spinor coefficients required')
    up, down = c[:, 0], c[:, 1]
    uu, ud = up.conj() @ up.T, up.conj() @ down.T
    du, dd = down.conj() @ up.T, down.conj() @ down.T
    return np.stack((ud+du, -1j*ud+1j*du, uu-dd)), uu+dd


def projector_spin_correction(projectors, metric):
    """Contract p* Q sigma p for a spin-independent projector metric.

    ``projectors[k,n,s,i]`` is <p_i|u_ns>; ``metric[k,i,j]`` is the
    AE-minus-pseudo overlap of partial waves. A shared [i,j] metric is
    accepted. Projector ordering and gauge must be established externally.
    A spin-dependent Q is deliberately unsupported. Q itself may be
    indefinite: it is a correction, not the complete physical metric.
    """
    p = np.asarray(projectors, dtype=np.complex128)
    q = np.asarray(metric, dtype=np.complex128)
    require(p.ndim == 4 and min(p.shape) > 0 and p.shape[2] == 2
            and np.isfinite(p).all(), 'finite projectors [k,band,2,projector] required')
    nk, nb, _, ni = p.shape
    if q.shape == (ni, ni):
        q = np.broadcast_to(q, (nk, ni, ni))
    require(q.shape == (nk, ni, ni) and np.isfinite(q).all(),
            'spin-independent projector metric [projector,projector] or [k,projector,projector] required')
    _hermitian(q, 'projector correction metric')
    blocks = np.einsum('knsi,kij,kmtj->kstnm', p.conj(), q, p, optimize=True)
    uu, ud, du, dd = blocks[:, 0, 0], blocks[:, 0, 1], blocks[:, 1, 0], blocks[:, 1, 1]
    return np.stack((ud+du, -1j*ud+1j*du, uu-dd), axis=1), uu+dd


def validate_spin(data):
    """Check source declarations, shapes, Hermiticity and metric semantics."""
    m = data.metadata
    require(isinstance(m, dict) and m.get('schema') == SCHEMA and m.get('version') == 1
            and type(m.get('version')) is int and m.get('complete') is True,
            'complete version1 spin cache required')
    require(m.get('units') == UNITS and m.get('spin_components') == ['x', 'y', 'z']
            and m.get('spinor_components') == 2 and m.get('spin_multiplicity') == 1
            and m.get('spinor_storage') == 'all_up_then_all_down' and m.get('gauge') == GAUGE,
            'spin units, components, spinor storage or gauge disagree')
    _basis(m.get('spin_basis'))
    require(_hash(m.get('source_wavecar_sha256')), 'source WAVECAR checksum required')
    require(m.get('normalization_applied') is False, 'spin normalization must not be applied')
    nk, nb = m.get('nkpoints'), m.get('nbands')
    require(type(nk) is int and nk > 0 and type(nb) is int and nb > 0
            and 8*nk*nb*nb*16 <= MAX_ARRAY_BYTES, 'positive bounded spin dimensions required')
    _finite(data.kpoints_fractional, (nk, 3), np.float64, 'kpoints_fractional')
    _finite(data.energies_eV, (nk, nb), np.float64, 'energies_eV')
    _finite(data.lattice_A, (3, 3), np.float64, 'lattice_A')
    require(abs(np.linalg.det(data.lattice_A)) > 1e-12, 'nonsingular lattice required')
    _finite(data.band_indices, (nb,), np.int64, 'band_indices')
    require(np.all(data.band_indices > 0) and len(set(data.band_indices.tolist())) == nb
            and type(m.get('source_nbands')) is int and np.all(data.band_indices <= m['source_nbands']),
            'unique 1-based band indices within source band count required')
    residual = {}
    for key, shape in [('spin_pauli', (nk, 3, nb, nb)), ('overlap', (nk, nb, nb)),
                       ('pseudo_spin_pauli', (nk, 3, nb, nb)), ('pseudo_overlap', (nk, nb, nb))]:
        array = getattr(data, key)
        _finite(array, shape, np.complex128, key)
        residual[key] = _hermitian(array, key)
    require(np.all(np.diagonal(data.pseudo_overlap, axis1=-2, axis2=-1).real > 0),
            'pseudo overlap must have positive band norms')
    metric_error = float(np.max(abs(data.overlap-np.eye(nb))))
    pseudo_error = float(np.max(abs(data.pseudo_overlap-np.eye(nb))))
    scope = m.get('operator_scope')
    if scope == 'raw_pseudo':
        require(m.get('physical_paw_validated') is False
                and np.array_equal(data.spin_pauli, data.pseudo_spin_pauli)
                and np.array_equal(data.overlap, data.pseudo_overlap),
                'raw pseudo arrays must remain unchanged and cannot claim PAW completeness')
    elif scope == 'paw_augmented':
        require(m.get('physical_paw_validated') is True and isinstance(m.get('augmentation'), dict),
                'audited augmentation provenance required')
        _validate_augmentation_metadata(m['augmentation'], data)
        require(metric_error <= PHYSICAL_METRIC_ATOL,
                'PAW overlap is not identity within physical metric tolerance; no normalization is applied')
        # For any represented state, |<sigma_a>| <= <1>. This is a
        # matrix inequality, not only a check of diagonal expectation values.
        # It catches unit/sign/contraction failures in otherwise Hermitian
        # supplied corrections without forcing eigenvalues to +/-1.
        bound_min = min(float(np.linalg.eigvalsh(data.overlap[:, None]+sign*data.spin_pauli).min())
                        for sign in (-1, 1))
        require(bound_min >= -PHYSICAL_METRIC_ATOL,
                'PAW Pauli operator violates the overlap-relative spin bound')
    else:
        raise ValueError('operator_scope must be raw_pseudo or paw_augmented')
    return {'hermiticity_max_abs': residual, 'overlap_identity_max_abs': metric_error,
            'pseudo_overlap_identity_max_abs': pseudo_error,
            'physical_metric_tolerance': PHYSICAL_METRIC_ATOL}


def extract_wavecar_spin(wavecar_path, *, spin_basis, bands=None):
    """Extract all complex Pauli matrices from a two-component WAVECAR."""
    _basis(spin_basis)
    path = Path(wavecar_path)
    digest = sha256(path)
    wave = Wavecar(path, spinor_components=2, spin=1)
    source_n = wave.header.nbands
    if bands is None:
        bands = list(range(1, source_n+1))
    require(isinstance(bands, (list, tuple, np.ndarray)) and len(bands) > 0
            and all(isinstance(b, (int, np.integer)) and not isinstance(b, (bool, np.bool_)) for b in bands),
            'nonempty integer band indices required')
    band_indices = np.asarray(bands, dtype=np.int64)
    require(np.all((band_indices >= 1) & (band_indices <= source_n))
            and len(set(band_indices.tolist())) == len(band_indices), 'unique band indices within WAVECAR required')
    nk, nb = wave.header.nkpoints, len(band_indices)
    require(8*nk*nb*nb*16 <= MAX_ARRAY_BYTES, 'spin output exceeds bounded storage limit')
    spin = np.empty((nk, 3, nb, nb), dtype=np.complex128)
    overlap = np.empty((nk, nb, nb), dtype=np.complex128)
    for ik in range(nk):
        spin[ik], overlap[ik] = pauli_from_coefficients(wave.coefficients(ik, band_indices.tolist()))
    meta = dict(schema=SCHEMA, version=1, complete=True, nkpoints=nk, nbands=nb,
        source_nbands=source_n, source_wavecar_sha256=digest, source_wavecar_file=path.name,
        source_format='VASP_WAVECAR', units=UNITS.copy(), spin_components=['x', 'y', 'z'],
        spinor_components=2, spin_multiplicity=1, spinor_storage='all_up_then_all_down',
        spin_basis=spin_basis, spin_basis_association='Caller declaration; SAXIS is not encoded in WAVECAR.',
        gauge=GAUGE, operator_scope='raw_pseudo', physical_paw_validated=False,
        normalization_applied=False, energy_reference='Unchanged WAVECAR eigenvalue zero in eV',
        matrix_convention='[n,m] = bra(n) operator ket(m)', physical_spin='S=(hbar/2)*spin_pauli',
        limitation='Pseudo coefficients only; PAW augmentation is absent. Raw pseudo overlap is preserved.')
    data = SpinOperators(wave.kpoints.copy(), wave.energies[:, band_indices-1].copy(),
        wave.header.lattice.copy(), band_indices, spin, overlap, spin, overlap, meta)
    meta['validation'] = validate_spin(data)
    require(digest == sha256(path), 'WAVECAR changed during extraction')
    return data


def _validate_augmentation_metadata(m, source):
    require(isinstance(m, dict) and m.get('schema') == AUGMENTATION_SCHEMA
            and type(m.get('version')) is int and m['version'] == 1 and m.get('complete') is True
            and m.get('producer_status') == 'PASS', 'complete successful version1 augmentation producer required')
    require(m.get('source_wavecar_sha256') == source.metadata['source_wavecar_sha256']
            and m.get('spin_basis') == source.metadata['spin_basis']
            and m.get('gauge') == GAUGE and m.get('units') == UNITS,
            'augmentation source WAVECAR, spin basis, gauge or units mismatch')
    require(isinstance(m.get('producer'), str) and bool(m['producer'].strip()), 'augmentation producer required')
    hashes = m.get('source_sha256')
    require(isinstance(hashes, dict) and len(hashes) > 0
            and all(isinstance(k, str) and k and _hash(v) for k, v in hashes.items()),
            'augmentation audited source checksums required')


def add_paw_augmentation(source, delta_spin_pauli, delta_overlap, *, metadata,
                         kpoints_fractional, energies_eV, lattice_A, band_indices):
    """Attach a same-gauge PAW correction; never normalize the resulting states."""
    validate_spin(source)
    require(source.metadata['operator_scope'] == 'raw_pseudo', 'PAW augmentation can only be applied once')
    _validate_augmentation_metadata(metadata, source)
    for name, actual, expected, atol in (
        ('kpoints', kpoints_fractional, source.kpoints_fractional, 1e-10),
        ('energies', energies_eV, source.energies_eV, 1e-7),
        ('lattice', lattice_A, source.lattice_A, 1e-10)):
        a = np.asarray(actual)
        require(a.shape == expected.shape and np.isfinite(a).all()
                and np.allclose(a, expected, rtol=0, atol=atol), 'augmentation '+name+' mismatch')
    require(np.array_equal(band_indices, source.band_indices), 'augmentation band indices mismatch')
    ds, do = np.asarray(delta_spin_pauli, dtype=np.complex128), np.asarray(delta_overlap, dtype=np.complex128)
    _finite(ds, source.spin_pauli.shape, np.complex128, 'delta_spin_pauli')
    _finite(do, source.overlap.shape, np.complex128, 'delta_overlap')
    _hermitian(ds, 'delta_spin_pauli'); _hermitian(do, 'delta_overlap')
    m = dict(source.metadata, operator_scope='paw_augmented', physical_paw_validated=True,
             augmentation=json.loads(json.dumps(metadata, allow_nan=False)),
             limitation='PAW spin supplied by the declared producer; finite band coverage and source DFT convergence remain separate.')
    data = SpinOperators(source.kpoints_fractional.copy(), source.energies_eV.copy(),
        source.lattice_A.copy(), source.band_indices.copy(), source.spin_pauli+ds,
        source.overlap+do, source.pseudo_spin_pauli.copy(), source.pseudo_overlap.copy(), m)
    m['validation'] = validate_spin(data)
    return data


def _load_arrays(path, names):
    with zipfile.ZipFile(path) as archive:
        entries = archive.infolist()
        require(len(entries) == len(names) and {e.filename for e in entries} == {k+'.npy' for k in names},
                'spin archive entries disagree')
        require(sum(e.file_size for e in entries) <= MAX_ARRAY_BYTES+1024*1024,
                'spin archive exceeds bounded uncompressed storage limit')
    with np.load(path, allow_pickle=False) as arrays:
        require(set(arrays.files) == set(names), 'spin archive arrays disagree')
        return {key: arrays[key] for key in names}


def import_paw_augmentation(source, npz_path, metadata_path):
    """Import the explicit delta-matrix NPZ/JSON producer interchange format."""
    npz_path, metadata_path = Path(npz_path), Path(metadata_path)
    digests = (sha256(npz_path), sha256(metadata_path))
    m = json.loads(metadata_path.read_text())
    require(m.get('data_npz_sha256') == digests[0], 'augmentation NPZ checksum mismatch')
    arrays = _load_arrays(npz_path, AUGMENTATION_ARRAYS)
    m['imported_files'] = {'npz': npz_path.name, 'metadata': metadata_path.name}
    m['imported_sha256'] = {'npz': digests[0], 'metadata': digests[1]}
    data = add_paw_augmentation(source, metadata=m, **arrays)
    require(digests == (sha256(npz_path), sha256(metadata_path)), 'augmentation source changed while reading')
    return data


def write_spin(directory, data):
    """Write a checked portable cache into a new directory."""
    validate_spin(data)
    out = Path(directory)
    require(not out.exists(), 'spin output directory exists')
    out.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(out/'spin.npz', **{name:getattr(data, name) for name in ARRAYS})
    m = dict(data.metadata, data_npz_sha256=sha256(out/'spin.npz'))
    (out/'spin.json').write_text(json.dumps(m, indent=2, allow_nan=False)+'\n')
    return m


def read_spin(directory):
    """Read a cache after integrity, conventions and matrix validation."""
    path = Path(directory)
    digests = (sha256(path/'spin.npz'), sha256(path/'spin.json'))
    m = json.loads((path/'spin.json').read_text())
    require(m.get('data_npz_sha256') == digests[0], 'spin NPZ checksum mismatch')
    data = SpinOperators(**_load_arrays(path/'spin.npz', ARRAYS), metadata=m)
    validate_spin(data)
    require(digests == (sha256(path/'spin.npz'), sha256(path/'spin.json')), 'spin cache changed while reading')
    return data
