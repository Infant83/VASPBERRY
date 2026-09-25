"""Occupation-weighted Kubo pairs and insulating bundles from actual VASP data.

Store undivided antisymmetric matrix products. Equal-occupation internal
pairs cancel before evaluating a denominator, including exact degeneracies.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
import json
from itertools import islice
from pathlib import Path

import numpy as np

from berry_data import CONDUCTANCE_QUANTUM_S, regions_from_spec, require
from exported_matrix_kubo import sha256
from vaspberry_transport import fermi_dirac
from wavecar_fukui import Wavecar, infer_uniform_grid

SCHEMA = 'vaspberry.kubo-pairs'
NATIVE_SCHEMA = 'VASPBERRY_BARE_MOMENTUM_KUBO_PAIRS_V1'
OPERATOR = 'WAVECAR_BARE_MOMENTUM_NO_PAW_NONLOCAL_VELOCITY'
CONNECTION = 'A_i=i<u|d/dk_i u>'
ARRAYS = ('k_ids', 'band_ids', 'pair_n', 'pair_m', 'kpoints_fractional',
          'weights', 'energies_eV', 'numerator_eV2_A2', 'lattice_A', 'reciprocal_inv_A')


@dataclass(frozen=True)
class PairData:
    metadata: dict
    k_ids: np.ndarray
    band_ids: np.ndarray
    pair_n: np.ndarray
    pair_m: np.ndarray
    kpoints_fractional: np.ndarray
    weights: np.ndarray
    energies_eV: np.ndarray
    numerator_eV2_A2: np.ndarray
    lattice_A: np.ndarray
    reciprocal_inv_A: np.ndarray


def metadata_lines(path):
    """Read comments through the end: native completion is a final footer."""
    result = {}
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith('#') and '=' in line:
                key, value = line[1:].strip().split('=', 1)
                require(key not in result or result[key] == value, 'conflicting metadata: '+key)
                result[key] = value
    return result


def validate_geometry(data):
    nk = len(data.k_ids)
    require(isinstance(data.k_ids, np.ndarray) and data.k_ids.dtype.kind in 'iu'
            and nk > 0 and np.array_equal(data.k_ids, np.arange(1, nk+1)), 'consecutive integer k IDs required')
    for key, shape in (('kpoints_fractional', (nk, 3)), ('weights', (nk,)),
                       ('lattice_A', (3, 3)), ('reciprocal_inv_A', (3, 3))):
        a = getattr(data, key)
        require(isinstance(a, np.ndarray) and a.dtype == np.float64 and a.shape == shape
                and np.isfinite(a).all(), key+' requires finite float64 values and correct shape')
    require(np.linalg.det(data.lattice_A) > 0 and np.allclose(
        data.lattice_A @ data.reciprocal_inv_A.T, 2*np.pi*np.eye(3), atol=1e-8, rtol=1e-10),
        'right-handed real/reciprocal lattice mismatch')
    require(np.all(data.weights >= 0) and np.isclose(data.weights.sum(), 1., atol=1e-12, rtol=0),
            'weights must sum to one')
    spec = data.metadata.get('sampling', {})
    require(spec.get('kind') in ('points', 'uniform_full_2d'), 'sampling kind required')
    if spec['kind'] == 'uniform_full_2d':
        axes, mesh = spec.get('plane_axes'), spec.get('mesh')
        require(isinstance(axes, list) and len(axes) == 2 and len(set(axes)) == 2
                and all(type(a) is int and a in (0, 1, 2) for a in axes), 'two distinct plane axes required')
        require(isinstance(mesh, list) and len(mesh) == 2
                and all(type(n) is int and n >= 2 for n in mesh), 'two mesh dimensions >=2 required')
        other = next(a for a in range(3) if a not in axes)
        infer_uniform_grid(data.kpoints_fractional[:, axes+[other]], *mesh)
        require(np.allclose(data.weights, 1/nk, atol=1e-12, rtol=0), 'uniform full-mesh weights required')
    mult = data.metadata.get('spin_multiplicity')
    require(type(mult) is int and mult in (1, 2), 'explicit spin multiplicity required')
    if data.metadata.get('source_nspin') == 2 or data.metadata.get('spinor_components') == 2:
        require(mult == 1, 'explicit spin channels and SOC spinors have multiplicity one')
    return data


def validate_pairs(data):
    m = data.metadata
    require(m.get('schema') == SCHEMA and m.get('version') == 1 and m.get('complete') is True,
            'complete vaspberry.kubo-pairs version1 required')
    require(m.get('normalization') == 'numerator=-2Im(D_a_nm*D_b_mn)'
            and m.get('components') == ['yz', 'zx', 'xy']
            and m.get('units') == {'energy': 'eV', 'numerator': 'eV^2 Angstrom^2',
                                   'lattice': 'Angstrom', 'reciprocal': '1/Angstrom'},
            'pair conventions or units disagree')
    validate_geometry(data)
    nk, nb = len(data.k_ids), len(data.band_ids)
    require(all(isinstance(a, np.ndarray) and a.dtype.kind in 'iu' for a in
                (data.band_ids, data.pair_n, data.pair_m)), 'integer band and pair indices required')
    require(nb >= 2 and m.get('source_nbands') == nb
            and np.array_equal(data.band_ids, np.arange(1, nb+1)), 'all source bands required')
    n, p = np.triu_indices(nb, 1)
    require(np.array_equal(data.pair_n, n) and np.array_equal(data.pair_m, p),
            'complete unique lexicographic n<m pairs required')
    require(data.energies_eV.dtype == np.float64 and data.energies_eV.shape == (nk, nb)
            and np.isfinite(data.energies_eV).all()
            and np.all(np.diff(data.energies_eV, axis=1) >= -1e-8), 'finite energy-ordered bands required')
    require(data.numerator_eV2_A2.dtype == np.float64
            and data.numerator_eV2_A2.shape == (nk, len(n), 3)
            and np.isfinite(data.numerator_eV2_A2).all(), 'complete finite undivided pair numerators required')
    require(isinstance(m.get('source_operator'), dict) and bool(m.get('energy_reference')),
            'operator provenance and energy-reference description required')
    return data


def native_wavecar(csv_path, wavecar_path, spin, spinor_components, spin_multiplicity):
    w = Wavecar(wavecar_path, spin=spin, spinor_components=spinor_components)
    w.coefficients(0, [1])  # Resolve actual scalar/spinor layout.
    require(spin_multiplicity in (1, 2), 'spin multiplicity must be1 or2')
    require(spin_multiplicity == 1 or (w.header.ispin == 1 and spinor_components == 1),
            'SOC spinors or explicit collinear channels have multiplicity one')
    meta = metadata_lines(csv_path)
    require(meta.get('normalization') == 'STANDARD_MINUS_TWO_IM'
            and meta.get('operator') == OPERATOR and meta.get('berry_connection') == CONNECTION
            and meta.get('result_status') == 'PASS', 'current complete native physical output required')
    return w, meta


def import_native_pairs(csv_path, wavecar_path, *, spin, spinor_components, spin_multiplicity,
                        sampling, energy_reference):
    w, source = native_wavecar(csv_path, wavecar_path, spin, spinor_components, spin_multiplicity)
    nk, nb = w.energies.shape
    require(source.get('schema') == NATIVE_SCHEMA and source.get('pair_order') == 'n_lt_m'
            and source.get('result_kind') == 'UNORDERED_INTERBAND_NUMERATORS'
            and source.get('components') == 'yz,zx,xy'
            and source.get('numerator_units') == 'eV^2*Angstrom^2'
            and source.get('gap_definition') == 'ABS_EN_MINUS_EM'
            and source.get('index_base') == '1'
            and source.get('occupation_weighting') == 'NONE'
            and source.get('denominator_weighting') == 'NONE', 'undivided native unordered-pair CSV required')
    for key, expected in (('source_nbands', nb), ('source_nkpoints', nk),
                          ('source_nspin', w.header.ispin), ('spinor_components', spinor_components),
                          ('pairs_per_k', nb*(nb-1)//2), ('expected_rows', nk*nb*(nb-1)//2*w.header.ispin)):
        require(int(source[key]) == expected, key+' disagrees with WAVECAR')
    require(source.get('reciprocal_convention') == '2pi', 'native reciprocal convention required')
    for prefix, expected, atol in (('lattice_A', w.header.lattice, 1e-8),
                                   ('reciprocal_inv_A', w.header.reciprocal, 1e-6)):
        actual = np.array([[float(x) for x in source[f'{prefix}_{a+1}'].split(',')] for a in range(3)])
        require(actual.shape == (3, 3) and np.allclose(actual, expected, atol=atol, rtol=1e-7),
                prefix+' disagrees with WAVECAR')
    n, m = np.triu_indices(nb, 1)
    index = np.full((nb, nb), -1, dtype=int); index[n, m] = np.arange(len(n))
    values = np.full((nk, len(n), 3), np.nan)
    seen = np.zeros((nk, len(n)), dtype=bool)
    columns = ['spin', 'k_index', 'n_band', 'm_band', 'kx_frac', 'ky_frac', 'kz_frac',
               'energy_n_eV', 'energy_m_eV', 'numerator_yz_eV2_A2', 'numerator_zx_eV2_A2',
               'numerator_xy_eV2_A2', 'gap_eV']
    total_rows = 0
    with Path(csv_path).open() as handle:
        reader = csv.reader(line for line in handle if not line.startswith('#'))
        require(next(reader, None) == columns, 'native pair column names/order disagree')
        while True:
            # Bounded parsing buffer; validate coordinates/energies in arrays,
            # avoiding one NumPy call for every one of millions of pair rows.
            chunk = list(islice(reader, 32768))
            if not chunk:
                break
            block = np.asarray(chunk, dtype=float)
            require(block.ndim == 2 and block.shape[1] == len(columns) and np.isfinite(block).all(),
                    'complete finite native pair rows required')
            total_rows += len(block)
            ids = block[:, :4]
            require(np.all(ids == np.floor(ids)) and np.all(ids[:, 0] >= 1)
                    and np.all(ids[:, 0] <= w.header.ispin), 'integer indices and valid spin channels required')
            block = block[ids[:, 0] == spin]
            if not len(block):
                continue
            ids = block[:, 1:4]
            require(np.all((ids[:, 0] >= 1) & (ids[:, 0] <= nk)
                           & (ids[:, 1] >= 1) & (ids[:, 1] < ids[:, 2]) & (ids[:, 2] <= nb)),
                    'out-of-range k/pair')
            k, a, b = (ids.astype(np.int64)-1).T
            p = index[a, b]
            flat = k*len(n)+p
            require(len(np.unique(flat)) == len(flat) and not seen[k, p].any(), 'duplicate k/pair row')
            require(np.allclose(block[:, 4:7], w.kpoints[k], atol=1e-8, rtol=0),
                    'pair coordinates disagree with WAVECAR')
            e = np.column_stack((w.energies[k, a], w.energies[k, b]))
            require(np.allclose(block[:, 7:9], e, atol=1e-7, rtol=0), 'pair energies disagree with WAVECAR')
            require(np.allclose(block[:, 12], abs(block[:, 7]-block[:, 8]), atol=1e-7, rtol=0),
                    'pair gap disagrees')
            values[k, p] = block[:, 9:12]
            seen[k, p] = True
    require(total_rows == int(source['expected_rows']), 'native total row count disagrees')
    require(seen.all(), 'missing k/pair rows; do not integrate a partial export')
    meta = dict(schema=SCHEMA, version=1, complete=True, source_nbands=nb,
        normalization='numerator=-2Im(D_a_nm*D_b_mn)', components=['yz', 'zx', 'xy'],
        units={'energy': 'eV', 'numerator': 'eV^2 Angstrom^2', 'lattice': 'Angstrom', 'reciprocal': '1/Angstrom'},
        sampling=sampling, spin_multiplicity=spin_multiplicity, spin_channel_1based=spin,
        source_nspin=int(w.header.ispin), spinor_components=spinor_components,
        energy_reference=energy_reference,
        source_operator={'kind': 'canonical_momentum', 'accuracy_status': 'approximation',
                         'missing_terms': ['PAW augmentation and nonlocal/SOC velocity corrections']},
        native_metadata=source,
        source_sha256={str(csv_path): sha256(csv_path), str(wavecar_path): sha256(wavecar_path)})
    return validate_pairs(PairData(meta, np.arange(1, nk+1), np.arange(1, nb+1), n, m,
                         w.kpoints.copy(), np.full(nk, 1/nk), w.energies.copy(), values,
                         w.header.lattice.copy(), w.header.reciprocal.copy()))


def write_pairs(directory, data):
    validate_pairs(data)
    out = Path(directory); require(not out.exists(), 'pair output directory exists')
    out.mkdir(parents=True)
    np.savez_compressed(out/'pairs.npz', **{key: getattr(data, key) for key in ARRAYS})
    meta = dict(data.metadata, data_npz_sha256=sha256(out/'pairs.npz'))
    (out/'pairs.json').write_text(json.dumps(meta, indent=2, allow_nan=False)+'\n')
    return meta


def read_pairs(directory):
    d = Path(directory); meta = json.loads((d/'pairs.json').read_text())
    require(meta.get('data_npz_sha256') == sha256(d/'pairs.npz'), 'pair NPZ checksum mismatch')
    with np.load(d/'pairs.npz', allow_pickle=False) as arrays:
        require(all(k in arrays for k in ARRAYS), 'missing pair arrays')
        data = PairData(meta, **{k: arrays[k] for k in ARRAYS})
    return validate_pairs(data)


def integration_setup(data, mus, temperatures, mu_reference, region_spec, differences):
    validate_geometry(data)
    require(data.metadata['sampling']['kind'] == 'uniform_full_2d', 'Hall integration needs a full uniform 2D mesh')
    mus, ts = np.asarray(mus, dtype=float), np.asarray(temperatures, dtype=float)
    require(mus.ndim == ts.ndim == 1 and len(mus) and len(ts) and np.isfinite(mus).all()
            and np.isfinite(ts).all() and np.all(ts >= 0) and np.isfinite(mu_reference), 'finite mu/T required')
    require(len(np.unique(mus)) == len(mus) and len(np.unique(ts)) == len(ts), 'duplicate mu/T values')
    axes = data.metadata['sampling']['plane_axes']
    vector = np.cross(*data.reciprocal_inv_A[axes]); area = float(np.linalg.norm(vector)); normal = vector/area
    regions = regions_from_spec(data, region_spec)
    for name, left, right in differences:
        require(name and name not in regions and left in regions and right in regions and left != right,
                'difference requires unique name and two known distinct regions')
        regions[name] = regions[left].astype(float)-regions[right].astype(float)
    return mus, ts, area, normal, regions


def coalesce_energies(energies, threshold):
    """Explicit numerical-degeneracy approximation; keep each group span bounded."""
    effective = energies.copy(); groups = []
    first = 0
    while first < len(energies):
        last = first+1
        while last < len(energies) and energies[last]-energies[first] <= threshold:
            last += 1
        if last-first > 1:
            effective[first:last] = np.mean(energies[first:last])
            groups.append((first, last, float(np.min(energies[first:last])), float(np.max(energies[first:last]))))
        first = last
    return effective, groups


def pair_hall_spectrum(data, mus, temperatures, *, mu_reference, region_spec=None, differences=(),
                       degeneracy_threshold_eV=1e-7, degeneracy_policy='error',
                       allow_partial_bands=False, mu_chunk=32, pair_band_max=None):
    validate_pairs(data)
    require(np.isfinite(degeneracy_threshold_eV) and degeneracy_threshold_eV >= 0,
            'finite nonnegative degeneracy threshold required')
    require(degeneracy_policy in ('error', 'coalesce') and type(mu_chunk) is int and mu_chunk > 0,
            'valid degeneracy policy and positive mu chunk required')
    nb = len(data.band_ids)
    cap = nb if pair_band_max is None else pair_band_max
    require(type(cap) is int and 2 <= cap <= nb, 'pair band maximum must be an integer in 2:source_nbands')
    mus, ts, area, normal, regions = integration_setup(data, mus, temperatures, mu_reference, region_spec, differences)
    all_mus = np.unique(np.r_[mus, mu_reference])
    effective = data.energies_eV.copy(); grouped = []
    if degeneracy_policy == 'coalesce':
        for k, e in enumerate(effective):
            effective[k], groups = coalesce_energies(e, degeneracy_threshold_eV)
            grouped += [(k, *g) for g in groups]
        if np.any(ts == 0):
            for k, first, last, lo, hi in grouped:
                require(not np.any((all_mus >= lo) & (all_mus < hi)),
                        f'T=0 mu cuts unresolved energy group at k={k+1}, bands={first+1}:{last}')
    raw_cutoff_gap = float(np.min(data.energies_eV[:, cap]-data.energies_eV[:, cap-1])) if cap < nb else None
    effective_cutoff_gap = float(np.min(effective[:, cap]-effective[:, cap-1])) if cap < nb else None
    cutoff_gap = min(raw_cutoff_gap, effective_cutoff_gap) if cap < nb else None
    require(cutoff_gap is None or cutoff_gap > degeneracy_threshold_eV,
            'pair band cutoff splits an unresolved group; include the complete group')
    cutoff_occ = max(float(fermi_dirac(effective[:, cap-1], float(all_mus.max()), float(t)).max()) for t in ts)
    partial = cutoff_occ > 1e-8
    require(not partial or allow_partial_bands,
            ('highest source band occupied' if cap == nb else 'highest selected band occupied')
            + '; extend band window or explicitly allow partial bands')
    masks = np.array(list(regions.values()), dtype=float)
    # Per k and per mu chunk: no allocation spanning every mu*k*pair.
    values = np.zeros((len(ts), len(regions), len(all_mus), 4))
    selected = data.pair_m < cap
    n, m = data.pair_n[selected], data.pair_m[selected]
    max_occ_change = 0.
    mult = data.metadata['spin_multiplicity']
    for k, e in enumerate(effective):
        gap = abs(e[n]-e[m]); safe = gap > degeneracy_threshold_eV
        numerator = data.numerator_eV2_A2[k, selected] @ normal
        coeff = np.zeros(len(n)); coeff[safe] = numerator[safe]/gap[safe]**2
        factor = -area/(2*np.pi)*mult*data.weights[k]
        for ti, t in enumerate(ts):
            ref = fermi_dirac(e, mu_reference, float(t))
            ref_difference = ref[n]-ref[m]
            require(not np.any(ref_difference[~safe] != 0), f'unequal occupation in unresolved reference pair at k={k+1}')
            for start in range(0, len(all_mus), mu_chunk):
                stop = min(start+mu_chunk, len(all_mus))
                occ = fermi_dirac(e[None, :]-all_mus[start:stop, None], 0., float(t))
                difference = occ[:, n]-occ[:, m]
                require(not np.any(difference[:, ~safe] != 0),
                        f'unequal occupations in unresolved pair at k={k+1}; refine states or explicitly coalesce numerical degeneracies')
                delta_occ = occ-ref[None, :]
                delta_difference = delta_occ[:, n]-delta_occ[:, m]
                sigma = factor*(difference @ coeff)
                delta = factor*(delta_difference @ coeff)
                count = mult*data.weights[k]*occ[:, :cap].sum(axis=1)
                delta_count = mult*data.weights[k]*delta_occ[:, :cap].sum(axis=1)
                block = np.column_stack((sigma, count, delta, delta_count))
                values[ti, :, start:stop] += masks[:, k, None, None]*block[None, :, :]
                if grouped:
                    original = fermi_dirac(data.energies_eV[k][None, :]-all_mus[start:stop, None], 0., float(t))
                    max_occ_change = max(max_occ_change, float(np.max(abs(original-occ))))
    require(np.isfinite(values).all(), 'nonfinite pair Hall result')
    rows = result_rows(values, mus, ts, all_mus, regions, mu_reference)
    meta = hall_metadata(data, area, normal, region_spec, differences, regions, mu_reference)
    meta.update(method='occupation_weighted_interband_pairs',
        scope='partial_band_contribution' if partial else ('truncated_pair_space' if cap < nb else 'all_source_bands'),
        source_nbands=nb, pair_band_window=[1, cap], minimum_pair_cutoff_gap_eV=cutoff_gap,
        raw_pair_cutoff_gap_eV=raw_cutoff_gap, effective_pair_cutoff_gap_eV=effective_cutoff_gap,
        occupation_window_max_band=cap,
        band_window_interpretation='Both pair endpoints lie in the selected window; source eigenstates and provenance are retained. Converge the virtual-state sum separately from eigenstate accuracy.',
        highest_band_max_occupation=cutoff_occ,
        formula='sigma/(e^2/h)=-A_BZ/(2*pi)*g_s*sum_k w_k sum_n<m (f_n-f_m)*N_nm/(E_n-E_m)^2',
        delta_formula='replace pair occupations by [f_n(mu)-f_n(ref)]-[f_m(mu)-f_m(ref)] before summation',
        degeneracy={'policy': degeneracy_policy, 'threshold_eV': degeneracy_threshold_eV,
                    'coalesced_groups': len(grouped), 'max_energy_shift_eV': float(np.max(abs(effective-data.energies_eV))),
                    'max_occupation_shift': max_occ_change,
                    'interpretation': 'Coalesce is an explicit numerical-degeneracy approximation; raw energies remain unchanged.'},
        algorithm={'k_streamed': True, 'mu_chunk': mu_chunk, 'internal_equal_occupation_pairs': 'cancelled before division'})
    return rows, meta


def result_rows(values, mus, temperatures, all_mus, regions, mu_reference):
    rows = []; chosen = np.searchsorted(all_mus, mus)
    for ti, t in enumerate(temperatures):
        for ri, name in enumerate(regions):
            for mu, v in zip(mus, values[ti, ri, chosen]):
                rows.append(dict(mu_eV=float(mu), mu_minus_reference_eV=float(mu-mu_reference),
                    temperature_K=float(t), region=name, band_id=0, sigma_e2_over_h=float(v[0]),
                    sigma_S=float(v[0]*CONDUCTANCE_QUANTUM_S), delta_sigma_e2_over_h=float(v[2]),
                    electrons_per_cell=float(v[1]), delta_electrons_per_cell=float(v[3])))
    return rows


def hall_metadata(data, area, normal, region_spec, differences, regions, mu_reference):
    return dict(schema='vaspberry.hall-spectrum', version=1, complete=True,
        units={'sigma_e2_over_h': 'e^2/h', 'sigma_S': 'S', 'energy': 'eV', 'temperature': 'K'},
        energy_reference=data.metadata['energy_reference'], mu_reference_eV=float(mu_reference),
        source_metadata=data.metadata, spin_multiplicity=data.metadata['spin_multiplicity'],
        bz_area_inv_A2=area, normal=normal.tolist(), regions=region_spec or {'regions': []},
        differences=[list(v) for v in differences], region_points={n: int(np.count_nonzero(v)) for n, v in regions.items()},
        band_id_zero='total represented subspace; no gauge-dependent individual-band decomposition',
        zero_T_equal_energy_occupation=1., convergence_status='not_established_by_format_validation',
        limitations=['Intrinsic clean-limit charge sheet response; no scattering or lifetime broadening.',
                     'Region differences depend on the supplied partition; not a conserved valley-current operator.',
                     'Fixed-band chemical-potential scan; no self-consistent doping or thermal magnetism.',
                     'Operator approximations are preserved in source_metadata.'])


def bundle_hall_spectrum(csv_path, wavecar_path, mus, *, occupied, spin, spinor_components,
                         spin_multiplicity, sampling, energy_reference, mu_reference,
                         region_spec=None, differences=()):
    """Only the zero-temperature common insulating gap is valid for a fixed bundle."""
    from types import SimpleNamespace
    w, source = native_wavecar(csv_path, wavecar_path, spin, spinor_components, spin_multiplicity)
    nk, nb = w.energies.shape
    require(0 < occupied < nb and source.get('schema') == 'VASPBERRY_BARE_MOMENTUM_KUBO_BUNDLE_V1',
            'native occupied bundle with stored empty states required')
    require(tuple(int(source[k]) for k in ('band_min', 'band_max', 'band_rank', 'source_nbands'))
            == (1, occupied, occupied, nb), 'bundle must contain exactly the occupied leading bands')
    require(source.get('intermediate_bands') == 'EXTERNAL_TO_SELECTED_BUNDLE_WITHIN_SOURCE_NBANDS',
            'bundle external-state convention required')
    omega = np.full(nk, np.nan); seen = set()
    with Path(csv_path).open() as f:
        for row in csv.DictReader(line for line in f if not line.startswith('#')):
            if int(row['spin']) != spin:
                continue
            k = int(row['k_index'])-1
            require(0 <= k < nk and k not in seen, 'duplicate/out-of-range bundle point')
            q = np.array([float(row['k'+a+'_frac']) for a in 'xyz'])
            require(np.allclose(q, w.kpoints[k], atol=1e-8, rtol=0), 'bundle coordinates disagree with WAVECAR')
            gap = float(np.min(abs(w.energies[k, :occupied, None]-w.energies[k, None, occupied:])))
            require(gap > 1e-5 and np.isclose(float(row['min_external_gap_eV']), gap, atol=1e-7, rtol=0),
                    'bundle external gap disagrees with WAVECAR or is unresolved')
            omega[k] = float(row['omega_z_A2']); seen.add(k)
    require(len(seen) == nk and np.isfinite(omega).all(), 'complete finite bundle grid required')
    vbm, cbm = float(w.energies[:, :occupied].max()), float(w.energies[:, occupied:].min())
    requested = np.r_[mus, mu_reference]
    require(cbm > vbm and np.all((requested > vbm) & (requested < cbm)),
            'fixed occupied bundle supports mu strictly inside its positive global insulating gap only')
    metadata = dict(sampling=sampling, spin_multiplicity=spin_multiplicity, energy_reference=energy_reference,
                    source_nspin=int(w.header.ispin), spin_channel_1based=spin, spinor_components=spinor_components,
                    native_metadata=source, source_sha256={str(csv_path): sha256(csv_path), str(wavecar_path): sha256(wavecar_path)})
    data = SimpleNamespace(metadata=metadata, k_ids=np.arange(1, nk+1), kpoints_fractional=w.kpoints.copy(),
                           weights=np.full(nk, 1/nk), lattice_A=w.header.lattice.copy(), reciprocal_inv_A=w.header.reciprocal.copy())
    mus, ts, area, normal, regions = integration_setup(data, mus, [0.], mu_reference, region_spec, differences)
    require(np.max(abs(normal[:2])) < 1e-12, 'native bundle has only xy curvature; selected plane must be parallel to xy')
    all_mus = np.unique(np.r_[mus, mu_reference])
    vals = np.zeros((1, len(regions), len(all_mus), 4))
    for ri, mask in enumerate(regions.values()):
        vals[0, ri, :, 0] = -area/(2*np.pi)*spin_multiplicity*np.sum(data.weights*omega*normal[2]*mask)
        vals[0, ri, :, 1] = occupied*spin_multiplicity*np.sum(data.weights*mask)
    rows = result_rows(vals, mus, ts, all_mus, regions, mu_reference)
    meta = hall_metadata(data, area, normal, region_spec, differences, regions, mu_reference)
    meta.update(method='occupied_bundle_insulating_T0', scope='fixed_occupied_bundle', occupied_bands=[1, occupied],
                global_gap_eV=cbm-vbm, valence_max_eV=vbm, conduction_min_eV=cbm,
                formula='sigma/(e^2/h)=-A_BZ/(2*pi)*g_s*sum_k w_k Omega_bundle_normal',
                delta_formula='zero inside the same global insulating gap at T=0')
    return rows, meta
