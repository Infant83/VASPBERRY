"""Versioned, physical point-curvature data and bounded-memory 2D Hall quadrature.

This format is distinct from a Fukui plaquette flux. No reader divides an
already normalized curvature by two, invents missing states, or repairs data.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np

from exported_matrix_kubo import MatrixContractError, sha256
from vaspberry_transport import fermi_dirac
from wavecar_fukui import infer_uniform_grid, periodic_reciprocal_distance

SCHEMA = 'vaspberry.band-curvature'
VERSION = 1
NORMALIZATION = 'physical_Omega=-2Im;A=i<u|grad_k u>'
UNITS = {'energy': 'eV', 'curvature': 'Angstrom^2', 'lattice': 'Angstrom',
         'reciprocal': '1/Angstrom'}
ARRAYS = ('k_ids', 'band_ids', 'intermediate_band_ids', 'kpoints_fractional',
          'weights', 'energies_eV', 'omega_A2', 'valid_nondegenerate',
          'min_gap_eV', 'lattice_A', 'reciprocal_inv_A')
CONDUCTANCE_QUANTUM_S = 1.602176634e-19**2 / 6.62607015e-34


@dataclass(frozen=True)
class CurvatureData:
    metadata: dict
    k_ids: np.ndarray
    band_ids: np.ndarray
    intermediate_band_ids: np.ndarray
    kpoints_fractional: np.ndarray
    weights: np.ndarray
    energies_eV: np.ndarray
    omega_A2: np.ndarray
    valid_nondegenerate: np.ndarray
    min_gap_eV: np.ndarray
    lattice_A: np.ndarray
    reciprocal_inv_A: np.ndarray


def require(condition, message):
    if not condition:
        raise MatrixContractError(message)


def validate(data):
    m = data.metadata
    require(isinstance(m, dict) and m.get('schema') == SCHEMA
            and type(m.get('version')) is int and m['version'] == VERSION,
            'expected vaspberry.band-curvature version 1; plaquette flux is a different format')
    require(m.get('complete') is True and m.get('normalization') == NORMALIZATION,
            'complete physical curvature with explicit normalization required')
    require(m.get('units') == UNITS and m.get('weights_convention') == 'sum_one'
            and m.get('reciprocal_convention') == '2pi', 'incompatible units or weights')
    require(m.get('axes') == ['k', 'band', 'cartesian'] and m.get('components') == ['yz', 'zx', 'xy'],
            'explicit curvature axes and Cartesian component order required')
    require(type(m.get('source_nbands')) is int and m['source_nbands'] > 0, 'source_nbands required')
    require(type(m.get('spin_multiplicity')) is int and m['spin_multiplicity'] in (1, 2),
            'explicit spin_multiplicity 1 (spinor or spin channel) or 2 (scalar spin degenerate) required')
    if 'source_nspin' in m or 'spin_channel_1based' in m:
        require(type(m.get('source_nspin')) is int and m['source_nspin'] in (1, 2)
                and type(m.get('spin_channel_1based')) is int
                and 1 <= m['spin_channel_1based'] <= m['source_nspin'],
                'optional source_nspin and spin_channel_1based require a consistent integer pair')
        require(m['source_nspin'] != 2 or m['spin_multiplicity'] == 1,
                'an explicit collinear spin channel has multiplicity one')
    if 'spinor_components' in m:
        require(type(m['spinor_components']) is int and m['spinor_components'] in (1, 2),
                'optional spinor_components must be integer 1 or 2')
        require(m['spinor_components'] != 2 or m['spin_multiplicity'] == 1,
                'two-component spinors have multiplicity one')
    for key in ('method', 'energy_reference'):
        require(isinstance(m.get(key), str) and bool(m[key].strip()), key+' required')
    require(isinstance(m.get('source_operator'), dict) and bool(m['source_operator'].get('kind')),
            'source_operator.kind required')
    require(isinstance(m.get('provenance'), dict) and bool(m['provenance']), 'provenance required')
    for key in ('k_ids', 'band_ids', 'intermediate_band_ids'):
        a = getattr(data, key)
        require(isinstance(a, np.ndarray) and a.ndim == 1 and a.dtype == np.int64 and len(a) > 0,
                key+' must be nonempty int64')
        require(np.all(a > 0) and len(np.unique(a)) == len(a), key+' must be positive and unique')
        if key != 'k_ids':
            require(a.max() <= m['source_nbands'], key+' exceeds source_nbands')
    nk, nb = len(data.k_ids), len(data.band_ids)
    for key, shape in {'kpoints_fractional': (nk, 3), 'weights': (nk,), 'energies_eV': (nk, nb),
                       'lattice_A': (3, 3), 'reciprocal_inv_A': (3, 3)}.items():
        a = getattr(data, key)
        require(isinstance(a, np.ndarray) and a.dtype == np.float64 and a.shape == shape
                and np.isfinite(a).all(), key+' must be finite float64 with shape '+str(shape))
    require(np.all(data.weights >= 0) and np.isclose(data.weights.sum(), 1, atol=1e-12, rtol=0),
            'weights must sum to one; no automatic normalization')
    require(np.linalg.det(data.lattice_A) > 0 and np.allclose(
            data.lattice_A @ data.reciprocal_inv_A.T, 2*np.pi*np.eye(3), atol=1e-8, rtol=1e-10),
            'inconsistent right-handed lattice and 2pi reciprocal lattice')
    require(isinstance(data.valid_nondegenerate, np.ndarray)
            and data.valid_nondegenerate.dtype == np.bool_ and data.valid_nondegenerate.shape == (nk, nb),
            'valid_nondegenerate must be boolean (k,band)')
    require(isinstance(data.min_gap_eV, np.ndarray)
            and data.min_gap_eV.dtype == np.float64 and data.min_gap_eV.shape == (nk, nb)
            and np.isfinite(data.min_gap_eV).all() and np.all(data.min_gap_eV >= 0),
            'finite nonnegative minimum gaps required')
    require(np.all(data.min_gap_eV[data.valid_nondegenerate] > 0),
            'valid individual-band curvature requires a positive known minimum gap')
    available = m.get('available_components')
    require(isinstance(available, list) and len(available) == 3
            and all(type(x) is bool for x in available) and any(available), 'available_components required')
    require(isinstance(data.omega_A2, np.ndarray)
            and data.omega_A2.dtype == np.float64 and data.omega_A2.shape == (nk, nb, 3),
            'omega_A2 must be float64 (k,band,cartesian)')
    covered = data.valid_nondegenerate[:, :, None] & np.array(available)[None, None, :]
    require(np.isfinite(data.omega_A2[covered]).all() and np.isnan(data.omega_A2[~covered]).all(),
            'curvature must be finite where valid/available and NaN elsewhere')
    sampling = m.get('sampling', {})
    require(isinstance(sampling, dict), 'sampling must be a metadata object')
    require(sampling.get('kind') in ('points', 'uniform_full_2d'), 'sampling kind required')
    if sampling['kind'] == 'uniform_full_2d':
        mesh, axes = sampling.get('mesh'), sampling.get('plane_axes')
        require(isinstance(mesh, list) and len(mesh) == 2 and all(type(n) is int and n >= 2 for n in mesh),
                'uniform_full_2d requires two mesh sizes >=2')
        require(isinstance(axes, list) and len(axes) == 2 and len(set(axes)) == 2
                and all(type(a) is int and a in (0, 1, 2) for a in axes), 'two distinct plane_axes required')
        other = next(i for i in range(3) if i not in axes)
        q = data.kpoints_fractional[:, axes+[other]].copy()
        # Existing tested grid validator checks a full tensor grid, arbitrary offset and duplicates.
        infer_uniform_grid(q, *mesh)
        require(np.allclose(data.weights, 1/nk, atol=1e-12, rtol=0),
                'uniform full mesh requires uniform weights (not irreducible mesh weights)')
    return data


def write_curvature(directory, data):
    validate(data)
    out = Path(directory)
    require(not out.exists(), 'output directory exists; choose a new directory')
    out.mkdir(parents=True)
    np.savez_compressed(out/'curvature.npz', **{key: getattr(data, key) for key in ARRAYS})
    meta = dict(data.metadata, data_npz_sha256=sha256(out/'curvature.npz'))
    (out/'curvature.json').write_text(json.dumps(meta, indent=2, allow_nan=False)+'\n')
    with (out/'curvature.csv').open('w', newline='') as handle:
        w = csv.writer(handle)
        w.writerow(['k_id', 'band_id', 'kx_frac', 'ky_frac', 'kz_frac', 'weight', 'energy_eV',
                    'omega_x_A2', 'omega_y_A2', 'omega_z_A2', 'min_gap_eV', 'valid_nondegenerate'])
        for k, n in np.ndindex(data.energies_eV.shape):
            w.writerow([data.k_ids[k], data.band_ids[n], *data.kpoints_fractional[k], data.weights[k],
                        data.energies_eV[k, n], *[v if np.isfinite(v) else '' for v in data.omega_A2[k, n]],
                        data.min_gap_eV[k, n], bool(data.valid_nondegenerate[k, n])])
    return meta


def read_curvature(directory):
    d = Path(directory)
    meta = json.loads((d/'curvature.json').read_text())
    require(isinstance(meta, dict), 'curvature metadata must be a JSON object')
    require(meta.get('data_npz_sha256') == sha256(d/'curvature.npz'), 'curvature NPZ SHA256 mismatch')
    with np.load(d/'curvature.npz', allow_pickle=False) as a:
        require(all(key in a for key in ARRAYS), 'missing curvature arrays')
        data = CurvatureData(meta, **{key: a[key] for key in ARRAYS})
    return validate(data)


def base_metadata(*, source_nbands, spin_multiplicity, method, energy_reference,
                  sampling, source_operator, provenance, available_components=(True, True, True)):
    return dict(schema=SCHEMA, version=VERSION, complete=True, normalization=NORMALIZATION,
                units=UNITS.copy(), weights_convention='sum_one', reciprocal_convention='2pi',
                axes=['k', 'band', 'cartesian'], components=['yz', 'zx', 'xy'],
                available_components=list(available_components), source_nbands=int(source_nbands),
                spin_multiplicity=spin_multiplicity, method=method, energy_reference=energy_reference,
                sampling=sampling, source_operator=source_operator, provenance=provenance,
                convergence_status='not_established_by_format_validation')


def regions_from_spec(data, spec=None):
    """User-named disjoint periodic circles OR explicit global k-ID sets; rest is automatic."""
    regions = {'total': np.ones(len(data.k_ids), dtype=bool)}
    used = np.zeros(len(data.k_ids), dtype=bool)
    circles = []
    require(spec is None or (isinstance(spec, dict) and set(spec) == {'regions'}
            and isinstance(spec['regions'], list)), 'region JSON must contain a regions list')
    for item in (spec or {'regions': []})['regions']:
        require(isinstance(item, dict), 'region entry must be an object')
        name = item.get('name')
        require(isinstance(name, str) and name.strip() and name not in regions and name != 'rest',
                'region names must be unique; total/rest reserved')
        if 'k_ids' in item:
            require(set(item) == {'name', 'k_ids'}, 'k-ID region accepts name and k_ids only')
            ids = item['k_ids']
            require(isinstance(ids, list) and len(ids)>0 and all(type(v) is int for v in ids)
                    and len(set(ids)) == len(ids) and set(ids) <= set(data.k_ids), 'invalid region k_ids')
            mask = np.isin(data.k_ids, ids)
        else:
            require(set(item) == {'name', 'center_fractional', 'radius_inv_A'},
                    'circle region requires name, center_fractional and radius_inv_A')
            center = np.asarray(item['center_fractional'], dtype=float).copy()
            radius = item['radius_inv_A']
            require(center.shape == (3,) and np.isfinite(center).all()
                    and type(radius) in (int, float) and np.isfinite(radius) and radius > 0,
                    'finite center and positive radius required')
            # In a 2D integral, periodic images are restricted to the chosen plane.
            require(data.metadata['sampling']['kind'] == 'uniform_full_2d',
                    'circle regions require a declared uniform full 2D sampling plane')
            axes = data.metadata['sampling']['plane_axes']
            other = next(a for a in range(3) if a not in axes)
            require(abs((center[other]-data.kpoints_fractional[0, other]+.5)%1-.5)<1e-8,
                    'region center must lie in the sampled reciprocal plane')
            # The plane test is periodic, so use the sampled representative for
            # its fixed coordinate before the in-plane closest-image search.
            center[other] = data.kpoints_fractional[0, other]
            reciprocal = data.reciprocal_inv_A[axes+[other]]
            for previous_name, previous_center, previous_radius in circles:
                separation = periodic_reciprocal_distance(
                    center[axes+[other]], previous_center[axes+[other]], reciprocal)
                require(radius+previous_radius < separation,
                        'circle regions '+previous_name+' and '+name+' overlap or touch geometrically')
            circles.append((name, center, radius))
            mask = np.array([periodic_reciprocal_distance(q[axes+[other]], center[axes+[other]], reciprocal)
                             <= radius for q in data.kpoints_fractional])
        require(mask.any(), 'region '+name+' has no sampled points')
        require(not np.any(mask & used), 'region masks overlap; provide a disjoint partition')
        regions[name] = mask
        used |= mask
    regions['rest'] = ~used
    return regions


def hall_spectrum(data, mus, temperatures, *, mu_reference, region_spec=None,
                  differences=(), allow_partial_bands=False, band_resolved=False, mu_chunk=64):
    """Occupation-resolved sheet response; zero-T uses sorted cumulative sums.

    At E=mu the existing fermi_dirac convention is full occupation. Finite-T
    memory is O(mu_chunk*K*B); there is no allocation spanning all mu and T.
    Counts refer to represented bands and include the explicit spin multiplicity.
    """
    validate(data)
    require(data.metadata['sampling']['kind'] == 'uniform_full_2d',
            'Hall integration requires a validated uniform full 2D mesh, not points/path/irreducible samples')
    require(data.valid_nondegenerate.all(), 'invalid/degenerate band curvature cannot be integrated; use a subspace method')
    require(type(mu_chunk) is int and mu_chunk > 0, 'positive mu_chunk required')
    mus, temperatures = np.asarray(mus, dtype=float), np.asarray(temperatures, dtype=float)
    require(mus.ndim == temperatures.ndim == 1 and mus.size > 0 and temperatures.size > 0
            and np.isfinite(mus).all() and np.isfinite(temperatures).all()
            and np.all(temperatures >= 0) and np.isfinite(mu_reference), 'finite mu and nonnegative temperatures required')
    require(len(np.unique(mus)) == len(mus) and len(np.unique(temperatures)) == len(temperatures),
            'duplicate chemical potentials or temperatures')
    axes = data.metadata['sampling']['plane_axes']
    area_vector = np.cross(*data.reciprocal_inv_A[axes])
    area = float(np.linalg.norm(area_vector)); normal = area_vector/area
    needed = np.abs(normal) > 1e-12
    require(np.all(np.array(data.metadata['available_components'])[needed]),
            'curvature components for requested plane are unavailable')
    omega = np.einsum('knc,c->kn', data.omega_A2[:, :, needed], normal[needed])
    order_b = np.argsort(data.band_ids)
    require(np.all(np.diff(data.energies_eV[:, order_b], axis=1) >= -1e-8),
            'global band IDs must follow energy order for occupation-window checks')
    all_mus = np.unique(np.r_[mus, mu_reference])
    highest = data.energies_eV[:, order_b[-1]]
    cutoff_occ = max(float(fermi_dirac(highest, float(all_mus.max()), float(t)).max()) for t in temperatures)
    missing_bands = not np.array_equal(np.sort(data.band_ids), np.arange(1, data.metadata['source_nbands']+1))
    partial = missing_bands or cutoff_occ > 1e-8
    require(not partial or allow_partial_bands,
            'incomplete occupied-band window or occupied highest exported band; extend n bands or explicitly --allow-partial-bands')
    regions = regions_from_spec(data, region_spec)
    for name, left, right in differences:
        require(name not in regions and name and left in regions and right in regions and left != right,
                'difference must have unique name and two known distinct regions')
        regions[name] = regions[left].astype(float)-regions[right].astype(float)
    mult = data.metadata['spin_multiplicity']
    state_weights = np.broadcast_to(mult*data.weights[:, None], data.energies_eV.shape)
    factors = -area/(2*np.pi)*state_weights*omega
    groups = [(0, np.ones(len(data.band_ids), bool))]
    if band_resolved:
        groups += [(int(b), data.band_ids == b) for b in data.band_ids]
    rows = []
    # Sort each represented state at most twice (total and its own band), then
    # reuse the ordering across regions. Band-resolved work does not scale B^2.
    orderings = {}
    for band, bands in groups:
        energies = data.energies_eV[:, bands].ravel()
        order = np.argsort(energies, kind='stable')
        sorted_e = energies[order]
        orderings[band] = (order, np.searchsorted(sorted_e, all_mus, side='right'),
                          int(np.searchsorted(sorted_e, mu_reference, side='right')))
    chosen = np.searchsorted(all_mus, mus)
    channels = [(name, band, (mask, bands))
                for name, mask in regions.items() for band, bands in groups]
    for temperature in temperatures:
        values = np.empty((len(channels), len(all_mus), 4))
        if temperature == 0:
            for c, (name, band, (region_mask, band_mask)) in enumerate(channels):
                coeff = np.stack(((factors[:, band_mask]*region_mask[:, None]).ravel(),
                                  (state_weights[:, band_mask]*region_mask[:, None]).ravel()), axis=1)
                order, right, at_ref = orderings[band]
                coeff = coeff[order]
                cumulative = np.vstack((np.zeros((1, 2)), np.cumsum(coeff, axis=0)))
                values[c, :, :2] = cumulative[right]
                # Delta starts at the reference occupation boundary, so a huge
                # filled baseline cannot erase a small newly occupied signal.
                above = np.vstack((np.zeros((1, 2)), np.cumsum(coeff[at_ref:], axis=0)))
                below = np.vstack((np.zeros((1, 2)), np.cumsum(coeff[:at_ref][::-1], axis=0)))
                hi = right >= at_ref
                values[c, hi, 2:] = above[right[hi]-at_ref]
                values[c, ~hi, 2:] = -below[at_ref-right[~hi]]
        else:
            # Compute the Fermi function once per chunk, shared by every
            # region and band. No mu*K*B*regions*bands temporary is formed.
            ref_occ = fermi_dirac(data.energies_eV, mu_reference, float(temperature))
            for start in range(0, len(all_mus), mu_chunk):
                end = min(start+mu_chunk, len(all_mus))
                occ = fermi_dirac(data.energies_eV[None, :, :]-all_mus[start:end, None, None],
                                 0., float(temperature))
                for offset in (0, 2):
                    for r, (name, mask) in enumerate(regions.items()):
                        sigma = np.einsum('qkn,kn->qn', occ, factors*mask[:, None])
                        count = np.einsum('qkn,kn->qn', occ, state_weights*mask[:, None])
                        c = r*len(groups)
                        values[c, start:end, offset] = sigma.sum(axis=1)
                        values[c, start:end, offset+1] = count.sum(axis=1)
                        if band_resolved:
                            values[c+1:c+len(groups), start:end, offset] = sigma.T
                            values[c+1:c+len(groups), start:end, offset+1] = count.T
                    if offset == 0:
                        occ -= ref_occ[None, :, :]
        require(np.isfinite(values).all(), 'nonfinite Hall result; inspect input magnitudes and numerical overflow')
        for (name, band, mask), vals in zip(channels, values):
            for mu, v in zip(mus, vals[chosen]):
                rows.append(dict(mu_eV=float(mu), mu_minus_reference_eV=float(mu-mu_reference),
                                 temperature_K=float(temperature), region=name, band_id=band,
                                 sigma_e2_over_h=float(v[0]), sigma_S=float(v[0]*CONDUCTANCE_QUANTUM_S),
                                 delta_sigma_e2_over_h=float(v[2]),
                                 electrons_per_cell=float(v[1]), delta_electrons_per_cell=float(v[3])))
    metadata = dict(schema='vaspberry.hall-spectrum', version=1, complete=True,
                    formula='sigma/(e^2/h)=-A_BZ/(2*pi)*sum_kn w_k f_kn Omega_normal',
                    delta_formula='same integrand with f(mu,T)-f(mu_reference,T); no subtraction of filled-band baselines',
                    units={'sigma_e2_over_h': 'e^2/h', 'sigma_S': 'S', 'energy': 'eV', 'temperature': 'K'},
                    band_id_zero='sum of all represented bands; not a physical band ID',
                    scope='partial_band_contribution' if partial else 'all_represented_bands',
                    spin_multiplicity=mult, mu_reference_eV=float(mu_reference),
                    energy_reference=data.metadata['energy_reference'],
                    source_curvature_metadata=data.metadata, bz_area_inv_A2=area, normal=normal.tolist(),
                    regions=region_spec or {'regions': []}, differences=[list(v) for v in differences],
                    region_points={name: int(np.count_nonzero(mask)) for name, mask in regions.items()},
                    highest_band_max_occupation=cutoff_occ,
                    zero_T_equal_energy_occupation=1.0, algorithm={'zero_T': 'sorted cumulative sum', 'mu_chunk': mu_chunk},
                    convergence_status='not_established_by_format_validation',
                    limitations=['Intrinsic clean-limit charge response; no scattering or lifetime broadening.',
                                 'Regional differences depend on the supplied partition; not a conserved valley-current operator.',
                                 'State spin/layer character does not define a spin/layer current response.',
                                 'mu scan uses fixed electronic states; self-consistent doping changes are not included.'])
    return rows, metadata


def write_hall(directory, rows, metadata):
    out = Path(directory)
    require(not out.exists(), 'output directory exists; choose a new directory')
    out.mkdir(parents=True)
    with (out/'conductivity.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    # Long-form arrays make both CSV and NPZ independently usable without reshaping guesses.
    np.savez_compressed(out/'conductivity.npz', **{key: np.array([r[key] for r in rows]) for key in rows[0]})
    meta = dict(metadata, output_sha256={name: sha256(out/name) for name in ('conductivity.csv', 'conductivity.npz')})
    (out/'conductivity.json').write_text(json.dumps(meta, indent=2, allow_nan=False)+'\n')
    return meta
