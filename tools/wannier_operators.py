"""Hamiltonian and position operators in the Wannier90 effective-model format.

This is the paired ``*_HH_R.dat`` / ``*_AA_R.dat`` interface, not the usual
``*_hr.dat`` Hamiltonian alone. Real-space and pair-translation degeneracies
must already be absorbed. Sparse omitted entries are zero; repeated entries
are additive. No nonzero entries are discarded or repaired by this importer.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import islice
import json
from pathlib import Path
import re
import zipfile

import numpy as np

from berry_data import require
from exported_matrix_kubo import sha256

SCHEMA = 'vaspberry.wannier-operators'
ARRAYS = ('irvec', 'hamiltonian_eV', 'position_A', 'lattice_A')
UNITS = {'hamiltonian': 'eV', 'position': 'Angstrom', 'lattice': 'Angstrom'}
PHASE = 'exp(+2*pi*i*q_fractional.dot(R_integer))'
HERMITICITY_ATOL = 1e-8
HERMITICITY_RTOL = 1e-10
MAX_OPERATOR_BYTES = 2*1024**3
PARSE_ROWS = 32768


@dataclass(frozen=True)
class WannierOperators:
    irvec: np.ndarray
    hamiltonian_eV: np.ndarray
    position_A: np.ndarray
    lattice_A: np.ndarray
    metadata: dict


def _integer(field):
    value = field.strip()
    require(re.fullmatch(r'[+-]?\d+', value) is not None, 'integer operator field required')
    return int(value)


def _real(field):
    """Fortran F12.6 input, including implied decimals and D exponents."""
    value = field.strip().replace('D', 'E').replace('d', 'e')
    require(bool(value), 'empty operator real field')
    try:
        result = float(value)
    except ValueError:
        # Fortran also permits an exponent without its E letter.
        match = re.fullmatch(r'([+-]?(?:\d+(?:\.\d*)?|\.\d+))([+-]\d+)', value)
        require(match is not None, 'invalid operator real field')
        value = match[1]+'e'+match[2]
        result = float(value)
    if '.' not in re.split('[Ee]', value)[0]:
        result *= 1e-6
    return result


def _read_blocks(handle, n, nr, components, expected=None):
    """Stream bounded text buffers into one dense array; sum residual rows."""
    result = np.zeros((nr, components, n, n), dtype=np.complex128)
    width = 25+24*components
    vectors = []; seen = set(); current = None; block = -1; rows = 0
    while True:
        lines = list(islice(handle, PARSE_ROWS))
        if not lines:
            break
        ids = np.empty((len(lines), 3), dtype=np.int64)
        values = np.empty((len(lines), components), dtype=np.complex128)
        for ri, raw in enumerate(lines):
            line = raw.rstrip('\r\n')
            require(len(line) >= width and not line[width:].strip(),
                    'operator row has incorrect fixed-width field count')
            integer = [_integer(line[i:i+5]) for i in range(0, 25, 5)]
            r = tuple(integer[:3]); row, column = integer[3:]
            require(1 <= row <= n and 1 <= column <= n, 'operator orbital index outside declared dimension')
            if r != current:
                require(r not in seen, 'noncontiguous repeated R block')
                block += 1
                require(block < nr, 'too many operator R blocks')
                require(expected is None or r == tuple(expected[block]), 'HH/AA R block order or coverage disagrees')
                seen.add(r); vectors.append(r); current = r
            ids[ri] = block, row-1, column-1
            real = [_real(line[i:i+12]) for i in range(25, width, 12)]
            values[ri] = [complex(real[i], real[i+1]) for i in range(0, len(real), 2)]
        require(np.isfinite(values).all(), 'nonfinite operator values')
        for component in range(components):
            np.add.at(result[:, component], tuple(ids.T), values[:, component])
        rows += len(lines)
    require(block+1 == nr, 'missing operator R blocks')
    require(np.isfinite(result).all(), 'nonfinite accumulated operator values')
    return np.asarray(vectors, dtype=np.int64), result, rows


def validate_operators(data):
    m = data.metadata
    require(isinstance(m, dict) and m.get('schema') == SCHEMA and type(m.get('version')) is int
            and m['version'] == 1 and m.get('complete') is True, 'complete version1 Wannier operator cache required')
    require(m.get('source_format') == 'wannier90_effective_HH_R_AA_R'
            and m.get('units') == UNITS and m.get('fourier_phase') == PHASE
            and m.get('real_space_degeneracy') == 'already_absorbed'
            and m.get('position_components') == ['x', 'y', 'z'], 'operator units or Fourier conventions disagree')
    require(isinstance(m.get('energy_reference'), str) and bool(m['energy_reference'].strip()), 'energy reference required')
    require(isinstance(m.get('source_operator'), dict)
            and m['source_operator'].get('kind') == 'Hamiltonian_and_position_connection'
            and m['source_operator'].get('scope') == 'represented finite Wannier subspace',
            'full Hamiltonian/position finite-subspace provenance required')
    require(isinstance(m.get('source_sha256'), dict) and set(m['source_sha256']) == {'HH_R', 'AA_R'}
            and all(isinstance(v, str) and re.fullmatch('[0-9a-f]{64}', v)
                    for v in m['source_sha256'].values()), 'HH/AA source checksums required')
    spinor, multiplicity = m.get('spinor_components'), m.get('spin_multiplicity')
    require(type(spinor) is int and spinor in (1, 2) and type(multiplicity) is int
            and multiplicity in (1, 2) and (spinor == 1 or multiplicity == 1), 'invalid explicit spinor/multiplicity declaration')
    require(isinstance(data.lattice_A, np.ndarray) and data.lattice_A.dtype == np.float64
            and data.lattice_A.shape == (3, 3) and np.isfinite(data.lattice_A).all()
            and np.linalg.det(data.lattice_A) > 1e-12, 'finite right-handed nonsingular lattice required')
    r = data.irvec; h = data.hamiltonian_eV; a = data.position_A
    require(isinstance(r, np.ndarray) and r.dtype == np.int64 and r.ndim == 2
            and r.shape[1] == 3 and len(r) > 0, 'integer R vectors required')
    require(np.all((r >= -9999) & (r <= 99999)), 'R vector outside source I5 field range')
    require(isinstance(h, np.ndarray) and h.dtype == np.complex128 and h.ndim == 3
            and h.shape[0] == len(r) and h.shape[1] == h.shape[2] and h.shape[1] > 0,
            'complex128 Hamiltonian shape must be [R,orbital,orbital]')
    n = h.shape[1]
    require(isinstance(a, np.ndarray) and a.dtype == np.complex128 and a.shape == (len(r), 3, n, n),
            'complex128 position shape must be [R,Cartesian,orbital,orbital]')
    require(h.nbytes+a.nbytes <= MAX_OPERATOR_BYTES and np.isfinite(h).all() and np.isfinite(a).all(),
            'finite operator arrays within 2 GiB storage limit required')
    require(type(m.get('num_wann')) is int and m['num_wann'] == n
            and type(m.get('nrpts')) is int and m['nrpts'] == len(r), 'operator dimensions disagree with metadata')
    lookup = {tuple(v): i for i, v in enumerate(r)}
    require(len(lookup) == len(r) and (0, 0, 0) in lookup, 'unique R vectors including origin required')
    require(all(tuple(-v) in lookup for v in r), 'missing opposite R vector')
    residuals = {'hamiltonian_eV': 0., 'position_A': 0.}
    for i, v in enumerate(r):
        j = lookup[tuple(-v)]
        for key, array in (('hamiltonian_eV', h), ('position_A', a)):
            left, right = array[i], array[j].conj().swapaxes(-1, -2)
            error = abs(left-right)
            residuals[key] = max(residuals[key], float(error.max()))
            require(np.all(error <= HERMITICITY_ATOL+HERMITICITY_RTOL*np.maximum(abs(left), abs(right))),
                    key+' violates R/-R Hermiticity; input is not symmetrized')
    return residuals


def load_effective_operators(hh_path, aa_path, lattice_A, *, spinor_components,
                             spin_multiplicity, energy_reference):
    """Load full operators; lattice rows are direct vectors in Angstrom.

    Both files must contain the same contiguous R blocks. In this effective
    convention each block has unit degeneracy, and repeated matrix entries ADD.
    The caller declares spin and energy zero; these are not present in the files.
    """
    hh_path, aa_path = Path(hh_path), Path(aa_path)
    require(hh_path.resolve() != aa_path.resolve(), 'distinct HH and AA input files required')
    hashes = {'HH_R': sha256(hh_path), 'AA_R': sha256(aa_path)}
    with hh_path.open(encoding='ascii') as handle:
        require(bool(handle.readline().strip()), 'HH comment header required')
        n, nr = _integer(handle.readline()), _integer(handle.readline())
        require(n > 0 and nr > 0 and 4*nr*n*n*16 <= MAX_OPERATOR_BYTES,
                'positive dimensions within 2 GiB operator storage limit required')
        r, h, hrows = _read_blocks(handle, n, nr, 1)
    with aa_path.open(encoding='ascii') as handle:
        require(bool(handle.readline().strip()), 'AA comment header required')
        _, a, arows = _read_blocks(handle, n, nr, 3, expected=r)
    metadata = dict(schema=SCHEMA, version=1, complete=True,
        source_format='wannier90_effective_HH_R_AA_R', num_wann=n, nrpts=nr,
        units=UNITS.copy(), fourier_phase=PHASE, real_space_degeneracy='already_absorbed',
        position_components=['x', 'y', 'z'], spinor_components=spinor_components,
        spin_multiplicity=spin_multiplicity, energy_reference=energy_reference,
        source_operator={'kind': 'Hamiltonian_and_position_connection',
            'scope': 'represented finite Wannier subspace',
            'physical_accuracy': 'not established by format validation; retain producer provenance'},
        source_sha256=hashes, source_files={'HH_R': hh_path.name, 'AA_R': aa_path.name},
        source_rows={'HH_R': hrows, 'AA_R': arows},
        sparse_semantics='omitted entries are zero; duplicate entries summed without a magnitude cutoff',
        lattice_association='Caller-supplied direct lattice rows; operator files do not encode the lattice.',
        hermiticity_tolerance={'absolute_eV_or_Angstrom': HERMITICITY_ATOL, 'relative': HERMITICITY_RTOL})
    data = WannierOperators(r, h[:, 0], a, np.asarray(lattice_A, dtype=np.float64).copy(), metadata)
    metadata['hermiticity_residual'] = validate_operators(data)
    require(hashes == {'HH_R': sha256(hh_path), 'AA_R': sha256(aa_path)}, 'operator source changed while reading')
    return data


def write_operators(directory, data):
    """Write a checked normalized NPZ/JSON cache to a new directory."""
    validate_operators(data)
    out = Path(directory)
    require(not out.exists(), 'operator output directory exists')
    out.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(out/'operators.npz', **{key: getattr(data, key) for key in ARRAYS})
    meta = dict(data.metadata, data_npz_sha256=sha256(out/'operators.npz'))
    (out/'operators.json').write_text(json.dumps(meta, indent=2, allow_nan=False)+'\n')
    return meta


def read_operators(directory):
    """Validate cache integrity, conventions, shapes and real-space Hermiticity."""
    directory = Path(directory)
    meta = json.loads((directory/'operators.json').read_text())
    require(meta.get('data_npz_sha256') == sha256(directory/'operators.npz'), 'operator NPZ checksum mismatch')
    with zipfile.ZipFile(directory/'operators.npz') as archive:
        entries = archive.infolist()
        require(len(entries) == len(ARRAYS) and {e.filename for e in entries} == {k+'.npy' for k in ARRAYS},
                'operator cache archive entries disagree')
        require(sum(e.file_size for e in entries) <= MAX_OPERATOR_BYTES+MAX_OPERATOR_BYTES//2+4096,
                'operator cache exceeds bounded uncompressed storage limit')
    with np.load(directory/'operators.npz', allow_pickle=False) as arrays:
        require(len(arrays.files) == len(ARRAYS) and set(arrays.files) == set(ARRAYS), 'operator cache arrays disagree')
        data = WannierOperators(**{key: arrays[key] for key in ARRAYS}, metadata=meta)
    validate_operators(data)
    return data
