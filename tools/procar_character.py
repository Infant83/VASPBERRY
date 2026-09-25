#!/usr/bin/env python3
"""Reusable atom/orbital/spin characters and projected charge-Hall attribution.

This is not a spin-current, layer-current or orbital-current Kubo operator.
Use matching noncollinear LORBIT=11 PROCAR/WAVECAR/OUTCAR from one static run.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
from pathlib import Path
import re

import numpy as np

from berry_data import require
from exported_matrix_kubo import sha256
from kubo_pairs import integration_setup, read_pairs, validate_pairs
from procar_projection import parse_procar, validate_state_alignment, _rotation
from vaspberry_transport import fermi_dirac, __version__
from waveder_hall import outcar_arrays, outcar_value, boolean
from wavecar_fukui import Wavecar

COMPONENTS = ('charge', 'pauli_axis', 'plus', 'minus')
CHAR_ARRAYS = ('kpoints_fractional', 'energies_eV', 'occupations', 'weights',
               'lattice_A', 'reciprocal_inv_A', 'raw_cartesian', 'characters',
               'all_ions_cartesian', 'printed_state_cartesian')
DIAGNOSTICS = ('$unweighted', '$all_projected', '$unprojected_residual')
ASSOCIATION = ('PROCAR, WAVECAR and OUTCAR are asserted to come from one unmodified static run. '
               'Coordinates, energies, occupations, lattice and dimensions are checked; these '
               'checks cannot independently authenticate the run or a degenerate-state gauge.')


@dataclass(frozen=True)
class CharacterData:
    metadata: dict
    kpoints_fractional: np.ndarray
    energies_eV: np.ndarray
    occupations: np.ndarray
    weights: np.ndarray
    lattice_A: np.ndarray
    reciprocal_inv_A: np.ndarray
    raw_cartesian: np.ndarray
    characters: np.ndarray
    all_ions_cartesian: np.ndarray
    printed_state_cartesian: np.ndarray


def json_write(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')


def write_rows(path, rows):
    require(bool(rows), 'empty output table')
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def outcar_rotation(text):
    """Read the actual printed forward rotation, without assuming SAXIS=z."""
    blocks = list(re.finditer(r'transformation matrix from SAXIS to cartesian coordinates\s*\n[-\s]+\n', text))
    require(len(blocks) == 1, 'OUTCAR needs one unambiguous SAXIS-to-Cartesian matrix')
    rows = text[blocks[0].end():].splitlines()[:3]
    matrix = []
    for line in rows:
        tokens = line.split()
        require(len(tokens) == 6 and tokens[1::2] == ['m_x', 'm_y', 'm_z'], 'malformed OUTCAR spin rotation')
        matrix.append([float(v) for v in tokens[::2]])
    return _rotation(np.asarray(matrix))


def group_characters(data, spec, axis):
    axis = np.asarray(axis, dtype=float)
    require(axis.shape == (3,) and np.isfinite(axis).all()
            and np.isclose(np.linalg.norm(axis), 1., atol=1e-10, rtol=0),
            '--axis requires a Cartesian unit vector')
    require(isinstance(spec, dict) and set(spec) == {'groups'} and isinstance(spec['groups'], list)
            and bool(spec['groups']), 'groups JSON requires a nonempty groups list')
    names, raw = [], []
    for item in spec['groups']:
        require(isinstance(item, dict) and set(item) in ({'name', 'ions'}, {'name', 'ions', 'orbitals'}),
                'group accepts name, ions and optional orbitals')
        name, ions = item['name'], item['ions']
        require(isinstance(name, str) and name.strip() and name not in names and name not in DIAGNOSTICS,
                'unique nonempty group names required; $ diagnostic names are reserved')
        require(isinstance(ions, list) and bool(ions) and all(type(i) is int and 1 <= i <= data.nions for i in ions)
                and len(set(ions)) == len(ions), 'group ions require unique 1-based integer IDs')
        if 'orbitals' in item:
            labels = item['orbitals']
            require(isinstance(labels, list) and bool(labels) and all(isinstance(s, str) for s in labels)
                    and len(set(labels)) == len(labels) and set(labels) <= set(data.orbital_names),
                    'group orbitals must be unique names printed in PROCAR: '+', '.join(data.orbital_names))
            block = np.take(data.ion_orbital, np.array(ions)-1, axis=2)
            values = np.take(block, [data.orbital_names.index(s) for s in labels], axis=3).sum(axis=(2, 3))
        else:
            values = data.ion_totals[:, :, np.array(ions)-1].sum(axis=2)
        values = values.copy()
        values[..., 1:] = data.pauli_cartesian(values)
        names.append(name); raw.append(values)
    raw = np.stack(raw, axis=2)
    q = raw[..., 0]; m = np.einsum('kngi,i->kng', raw[..., 1:], axis)
    characters = np.stack((q, m, (q+m)/2, (q-m)/2), axis=-1)
    return names, raw, characters


def validate_characters(data):
    m = data.metadata
    require(m.get('schema') == 'vaspberry.procar-character' and m.get('version') == 1
            and m.get('complete') is True and m.get('components') == list(COMPONENTS),
            'complete version 1 PROCAR character cache required')
    names = m.get('group_names')
    require(isinstance(names, list) and bool(names) and all(isinstance(x, str) and x and x not in DIAGNOSTICS for x in names)
            and len(set(names)) == len(names), 'unique named character groups required')
    require(data.energies_eV.ndim == 2 and min(data.energies_eV.shape) > 0, 'nonempty (k,band) energies required')
    nk, nb = data.energies_eV.shape
    shapes = {'kpoints_fractional': (nk, 3), 'energies_eV': (nk, nb), 'occupations': (nk, nb),
              'weights': (nk,), 'lattice_A': (3, 3), 'reciprocal_inv_A': (3, 3),
              'raw_cartesian': (nk, nb, len(names), 4), 'characters': (nk, nb, len(names), 4),
              'all_ions_cartesian': (nk, nb, 4), 'printed_state_cartesian': (nk, nb, 4)}
    for name, shape in shapes.items():
        a = getattr(data, name)
        require(a.dtype == np.float64 and a.shape == shape and np.isfinite(a).all(), 'invalid character array '+name)
    require(np.all(data.weights >= 0), 'nonnegative PROCAR weights required')
    require(np.allclose(data.lattice_A @ data.reciprocal_inv_A.T, 2*np.pi*np.eye(3), atol=1e-8, rtol=1e-10),
            'character lattice mismatch')
    q, m_axis, plus, minus = np.moveaxis(data.characters, -1, 0)
    axis = np.asarray(m.get('axis_cartesian'), dtype=float)
    require(axis.shape == (3,) and np.isfinite(axis).all() and np.isclose(np.linalg.norm(axis), 1., atol=1e-10, rtol=0),
            'invalid saved spin axis')
    require(np.allclose(q, data.raw_cartesian[..., 0], atol=1e-14, rtol=1e-13)
            and np.allclose(m_axis, data.raw_cartesian[..., 1:] @ axis, atol=1e-14, rtol=1e-13)
            and np.allclose(plus+minus, q, atol=1e-14, rtol=1e-13)
            and np.allclose(plus-minus, m_axis, atol=1e-14, rtol=1e-13), 'saved joint-spin identities disagree')
    require(isinstance(m.get('wavecar_sha256'), str) and len(m['wavecar_sha256']) == 64, 'WAVECAR provenance required')
    return data


def write_characters(directory, data):
    validate_characters(data)
    out = Path(directory); require(not out.exists(), 'output directory exists; choose a fresh directory')
    out.mkdir(parents=True)
    np.savez_compressed(out/'character.npz', **{name: getattr(data, name) for name in CHAR_ARRAYS})
    meta = dict(data.metadata, data_npz_sha256=sha256(out/'character.npz'))
    rows = []
    for k, b, g in np.ndindex(data.characters.shape[:3]):
        q, m, plus, minus = data.characters[k, b, g]
        rows.append(dict(k_id=k+1, band_id=b+1, group=meta['group_names'][g],
            kx_frac=data.kpoints_fractional[k, 0], ky_frac=data.kpoints_fractional[k, 1],
            kz_frac=data.kpoints_fractional[k, 2], energy_eV=data.energies_eV[k, b],
            occupation=data.occupations[k, b], weight=data.weights[k], charge=q, pauli_axis=m,
            plus=plus, minus=minus, mx=data.raw_cartesian[k, b, g, 1],
            my=data.raw_cartesian[k, b, g, 2], mz=data.raw_cartesian[k, b, g, 3]))
    write_rows(out/'character.csv', rows)
    diagnostic_rows = []
    for k, b in np.ndindex(data.energies_eV.shape):
        q = data.all_ions_cartesian[k, b, 0]
        diagnostic_rows.append(dict(k_id=k+1, band_id=b+1, all_ions_charge=q, charge_residual=1-q,
            mx=data.all_ions_cartesian[k, b, 1], my=data.all_ions_cartesian[k, b, 2],
            mz=data.all_ions_cartesian[k, b, 3], printed_state_charge=data.printed_state_cartesian[k, b, 0]))
    write_rows(out/'projection_diagnostics.csv', diagnostic_rows)
    meta['data_csv_sha256'] = sha256(out/'character.csv')
    meta['diagnostics_csv_sha256'] = sha256(out/'projection_diagnostics.csv')
    json_write(out/'character.json', meta)
    return meta


def read_characters(directory):
    d = Path(directory); meta = json.loads((d/'character.json').read_text())
    require(meta.get('data_npz_sha256') == sha256(d/'character.npz'), 'character NPZ checksum mismatch')
    with np.load(d/'character.npz', allow_pickle=False) as saved:
        require(set(CHAR_ARRAYS) <= set(saved.files), 'missing character arrays')
        data = CharacterData(meta, **{name: saved[name] for name in CHAR_ARRAYS})
    return validate_characters(data)


def project_command(args):
    require(not args.output_dir.exists(), 'output directory exists; choose a fresh directory')
    text = args.outcar.read_text()
    require('General timing and accounting informations for this job:' in text, 'completed static OUTCAR required')
    require(outcar_value(text, 'LNONCOLLINEAR', boolean), 'noncollinear OUTCAR required')
    for tag, value in [('ISPIN', 1), ('NSW', 0), ('LORBIT', 11)]:
        require(outcar_value(text, tag, int) == value, 'OUTCAR requires '+tag+'='+str(value))
    rotation = outcar_rotation(text)
    data = parse_procar(args.procar, noncollinear=True, spin_to_cartesian=rotation)
    wave = Wavecar(args.wavecar, spin=1, spinor_components=2)
    wave.coefficients(0, [1])
    require(wave.header.ispin == 1, 'noncollinear ISPIN=1 WAVECAR required')
    checks = validate_state_alignment(data, wave.kpoints, wave.energies)
    require(np.allclose(data.occupations, wave.occupations, atol=5.1e-9, rtol=0), 'PROCAR/WAVECAR occupations mismatch')
    nk, nb = wave.energies.shape
    require(outcar_value(text, 'NKPTS', int) == nk and outcar_value(text, 'NBANDS', int) == nb,
            'OUTCAR/WAVECAR dimensions mismatch')
    oq, oe, oo, lattice, grid = outcar_arrays(text, nk, nb, 1)
    require(np.allclose(lattice, wave.header.lattice, atol=6e-9, rtol=0), 'OUTCAR/WAVECAR lattice mismatch')
    require(np.allclose(grid[:, :3], wave.kpoints, atol=6e-9, rtol=0)
            and np.allclose(oq[0], wave.kpoints, atol=5.1e-5, rtol=0), 'OUTCAR/WAVECAR k coordinates mismatch')
    require(np.allclose(oe[0], wave.energies, atol=5.1e-5, rtol=0)
            and np.allclose(oo[0], wave.occupations, atol=5.1e-5, rtol=0), 'OUTCAR/WAVECAR final states mismatch')
    # VASP 5.4.4 prints this OUTCAR weight column with only three decimals;
    # PROCAR keeps eight. Never normalize either printed array to hide errors.
    require(np.allclose(grid[:, 3], data.weights, atol=5.1e-4, rtol=0), 'OUTCAR/PROCAR weights mismatch')
    checks['max_outcar_weight_delta'] = float(np.max(abs(grid[:, 3]-data.weights)))
    spec = json.loads(args.groups.read_text())
    names, raw, chars = group_characters(data, spec, args.axis)
    all_ions = data.ion_totals.sum(axis=2)
    all_ions[..., 1:] = data.pauli_cartesian(all_ions)
    printed_state = data.state_totals.copy()
    printed_state[..., 1:] = data.pauli_cartesian(printed_state)
    inputs = {name: sha256(getattr(args, name)) for name in ('procar', 'wavecar', 'outcar', 'groups')}
    meta = dict(schema='vaspberry.procar-character', version=1, complete=True, vaspberry_version=__version__,
        group_names=names, groups=spec, components=list(COMPONENTS),
        axis_cartesian=list(args.axis), spin_to_cartesian=rotation.tolist(),
        raw_components=['charge', 'mx', 'my', 'mz'], orbital_names=list(data.orbital_names),
        wavecar_sha256=inputs['wavecar'], source_sha256=inputs, alignment=checks,
        rounding_diagnostics=data.rounding_diagnostics, association=ASSOCIATION,
        all_ions_charge_range=[float(all_ions[..., 0].min()), float(all_ions[..., 0].max())],
        charge_residual_definition='1 - sum_ions printed ion charge totals; not clipped; no residual spin inferred',
        normalization='raw PAW projector weights; no clipping or normalization to unity',
        spin_definition='m is projected Pauli weight; p_plus/minus=(q +/- axis.m)/2; S_axis/hbar=m/2',
        limitations=['Overlapping atom/orbital groups do not form a partition.',
                     'Interstitial weight is not supplied by PROCAR.',
                     'Individual degenerate-state characters depend on the eigenvector gauge.'],
        provenance=provenance(args))
    return write_characters(args.output_dir, CharacterData(meta, wave.kpoints.copy(), wave.energies.copy(),
        wave.occupations.copy(), data.weights.copy(), wave.header.lattice.copy(), wave.header.reciprocal.copy(),
        raw, chars, all_ions, printed_state))


def match_pairs(character, pairs):
    validate_characters(character); validate_pairs(pairs)
    require(pairs.metadata.get('spinor_components') == 2 and pairs.metadata.get('spin_multiplicity') == 1,
            'character attribution requires a matching noncollinear spinor pair cache')
    hashes = pairs.metadata.get('source_sha256', {})
    require(character.metadata['wavecar_sha256'] in hashes.values(), 'pair cache and characters have different WAVECAR hashes')
    for key in ('kpoints_fractional', 'energies_eV', 'lattice_A', 'reciprocal_inv_A'):
        a, b = getattr(character, key), getattr(pairs, key)
        require(a.shape == b.shape and np.allclose(a, b, atol=1e-10, rtol=0), 'pair/character '+key+' mismatch')
    require(np.allclose(character.weights, pairs.weights, atol=5.1e-9, rtol=0),
            'PROCAR weights are not the same full uniform mesh as the pair cache')


def selected_curvature(pairs, bands, threshold):
    """Physical -2Im numerator divided once; reverse pair has opposite sign."""
    require(np.isfinite(threshold) and threshold > 0, 'positive finite degeneracy threshold required')
    require(bool(bands) and all(type(b) is int and 1 <= b <= len(pairs.band_ids) for b in bands)
            and len(set(bands)) == len(bands), 'unique selected band IDs within source range required')
    result = np.zeros((len(pairs.k_ids), len(bands), 3))
    minimum = np.full((len(pairs.k_ids), len(bands)), np.inf)
    for i, band in enumerate(bands):
        n = band-1
        keep = (pairs.pair_n == n) | (pairs.pair_m == n)
        gap = abs(pairs.energies_eV[:, pairs.pair_n[keep]]-pairs.energies_eV[:, pairs.pair_m[keep]])
        minimum[:, i] = gap.min(axis=1)
        bad = np.argwhere(gap <= threshold)
        require(not len(bad), f'band {band} has an unresolved individual-state gap; '
                'character-weighted Hall cannot coalesce states; use isolated bands or ordinary pair-hall')
        sign = np.where(pairs.pair_n[keep] == n, 1., -1.)
        result[:, i] = np.einsum('kpj,p,kp->kj', pairs.numerator_eV2_A2[:, keep], sign, 1/gap**2)
    require(np.isfinite(result).all(), 'nonfinite selected-band curvature')
    return result, minimum


def character_hall(character, pairs, bands, mus, temperatures, *, mu_reference,
                   region_spec=None, degeneracy_threshold_eV=1e-5):
    match_pairs(character, pairs)
    mus, ts, area, normal, regions = integration_setup(pairs, mus, temperatures, mu_reference, region_spec, ())
    omega, gaps = selected_curvature(pairs, bands, degeneracy_threshold_eV)
    omega_normal = omega @ normal
    chosen = np.array(bands)-1
    e = pairs.energies_eV[:, chosen]
    chars = character.characters[:, chosen]
    factor = -area/(2*np.pi)
    rows = []
    for temperature in ts:
        ref = fermi_dirac(e, float(mu_reference), float(temperature))
        for mu in mus:
            occ = fermi_dirac(e, float(mu), float(temperature))
            for region, mask in regions.items():
                for gi, group in enumerate(character.metadata['group_names']):
                    coeff = pairs.weights[:, None, None]*mask[:, None, None]*omega_normal[:, :, None]*chars[:, :, gi]
                    sigma = factor*np.einsum('kb,kbc->c', occ, coeff)
                    delta = factor*np.einsum('kb,kbc->c', occ-ref, coeff)
                    for ci, component in enumerate(COMPONENTS):
                        rows.append(dict(mu_eV=float(mu), mu_minus_reference_eV=float(mu-mu_reference),
                            temperature_K=float(temperature), region=region, group=group, component=component,
                            attribution_e2_over_h=float(sigma[ci]), delta_attribution_e2_over_h=float(delta[ci])))
                q = character.all_ions_cartesian[:, chosen, 0]
                diagnostic_chars = (np.ones_like(q), q, 1-q)
                for name, value in zip(DIAGNOSTICS, diagnostic_chars):
                    coeff = pairs.weights[:, None]*mask[:, None]*omega_normal*value
                    rows.append(dict(mu_eV=float(mu), mu_minus_reference_eV=float(mu-mu_reference),
                        temperature_K=float(temperature), region=region, group=name, component='charge',
                        attribution_e2_over_h=float(factor*np.sum(occ*coeff)),
                        delta_attribution_e2_over_h=float(factor*np.sum((occ-ref)*coeff))))
    meta = dict(schema='vaspberry.character-hall', version=1, complete=True, vaspberry_version=__version__,
        scope='selected_band_character_weighted_charge_Hall_attribution', selected_bands=bands,
        virtual_band_window=[1, len(pairs.band_ids)], min_selected_gap_eV=float(gaps.min()),
        highest_source_band_max_occupation=max(float(fermi_dirac(pairs.energies_eV[:, -1],
            float(max(np.max(mus), mu_reference)), float(t)).max()) for t in ts),
        degeneracy_threshold_eV=degeneracy_threshold_eV, degeneracy_policy='reject',
        group_names=character.metadata['group_names'], components=list(COMPONENTS),
        diagnostic_groups=dict(zip(DIAGNOSTICS, ['unweighted selected-band charge Hall',
            'charge Hall weighted by sum of all ion charge projections',
            'difference from unweighted selected-band charge Hall; no residual spin inferred'])),
        source_operator=pairs.metadata['source_operator'], energy_reference=pairs.metadata['energy_reference'],
        mu_reference_eV=float(mu_reference), bz_area_inv_A2=area, normal=normal.tolist(),
        regions=region_spec or {'regions': []}, region_points={n: int(np.count_nonzero(v)) for n, v in regions.items()},
        formula='-A_BZ/(2*pi) sum_k,n_selected w_k f_n(mu,T) Omega_n.normal c_n,group,component',
        delta_formula='replace f_n(mu,T) with f_n(mu,T)-f_n(mu_reference,T) before summation',
        normalization='Omega_n=sum_m!=n N_nm/(E_n-E_m)^2; native STANDARD_MINUS_TWO_IM, no extra half',
        character_metadata=character.metadata, pair_metadata=pairs.metadata,
        limitations=['Charge-Hall attribution weighted by state character; not a spin/layer/orbital-current conductivity.',
                     'Selected bands only: omitted bands can contribute, even when fully occupied.',
                     'Raw PAW projections are incomplete and possibly overlapping; no unity normalization.',
                     'Converge mesh, number of virtual bands and operator approximation separately.',
                     'Rigid-band mu/T scan, not self-consistent doping or thermal magnetism.'])
    return rows, meta, omega, gaps


def hall_command(args):
    require(not args.output_dir.exists(), 'output directory exists; choose a fresh directory')
    require(args.mu_num >= 2 and np.isfinite([args.mu_min, args.mu_max]).all() and args.mu_min < args.mu_max,
            'increasing finite mu range and at least two samples required')
    character, pairs = read_characters(args.character_dir), read_pairs(args.pairs_dir)
    spec = json.loads(args.regions.read_text()) if args.regions else None
    rows, meta, omega, gaps = character_hall(character, pairs, args.bands,
        np.linspace(args.mu_min, args.mu_max, args.mu_num), args.temperatures,
        mu_reference=args.mu_reference, region_spec=spec, degeneracy_threshold_eV=args.degeneracy_threshold_eV)
    out = args.output_dir; out.mkdir(parents=True)
    write_rows(out/'character_hall.csv', rows)
    columns = {name: np.array([row[name] for row in rows]) for name in rows[0]}
    np.savez_compressed(out/'character_hall.npz', **columns)
    np.savez_compressed(out/'selected_curvature.npz', k_ids=pairs.k_ids, band_ids=np.array(args.bands),
                        omega_A2=omega, min_gap_eV=gaps, kpoints_fractional=pairs.kpoints_fractional,
                        energies_eV=pairs.energies_eV[:, np.array(args.bands)-1])
    meta.update(provenance=provenance(args), data_csv_sha256=sha256(out/'character_hall.csv'),
                data_npz_sha256=sha256(out/'character_hall.npz'), selected_curvature_npz_sha256=sha256(out/'selected_curvature.npz'))
    json_write(out/'character_hall.json', meta)
    return meta


def plot_command(args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    require(not args.output_dir.exists(), 'output directory exists; choose a fresh directory')
    data = read_characters(args.character_dir)
    require(args.group in data.metadata['group_names'], 'unknown character group')
    require(1 <= args.band <= data.energies_eV.shape[1], 'band outside character cache')
    gi = data.metadata['group_names'].index(args.group)
    meta = json.loads((args.hall_dir/'character_hall.json').read_text())
    require(meta.get('data_npz_sha256') == sha256(args.hall_dir/'character_hall.npz'), 'Hall NPZ checksum mismatch')
    require(meta['character_metadata'].get('data_npz_sha256') == data.metadata.get('data_npz_sha256'),
            'plots require Hall and characters from the same cache')
    with np.load(args.hall_dir/'character_hall.npz', allow_pickle=False) as src:
        hall = {key: src[key] for key in src.files}
    take = (hall['group'] == args.group) & (hall['temperature_K'] == args.temperature) & (hall['region'] == args.region)
    require(np.any(take), 'requested group/temperature/region absent from Hall table')
    out = args.output_dir; out.mkdir(parents=True)
    fig, axes = plt.subplots(2, 2, figsize=(9, 7), constrained_layout=True)
    for ci, ax in enumerate(axes.flat):
        q = data.kpoints_fractional
        a, b = meta['pair_metadata']['sampling']['plane_axes']
        artist = ax.scatter(q[:, a], q[:, b], c=data.characters[:, args.band-1, gi, ci], s=35)
        fig.colorbar(artist, ax=ax)
        ax.set(title=COMPONENTS[ci], xlabel=f'fractional k{a+1}', ylabel=f'fractional k{b+1}')
    fig.suptitle(f'{args.group}, band {args.band}: raw projected character')
    for suffix in ('png', 'pdf', 'svg'): fig.savefig(out/('character.'+suffix), dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    field = 'delta_attribution_e2_over_h' if args.delta else 'attribution_e2_over_h'
    for component in COMPONENTS:
        mask = take & (hall['component'] == component)
        order = np.argsort(hall['mu_minus_reference_eV'][mask])
        ax.plot(hall['mu_minus_reference_eV'][mask][order], hall[field][mask][order], label=component)
    ax.set(xlabel='mu - reference (eV)', ylabel=('Delta ' if args.delta else '')+'charge-Hall attribution (e²/h)',
           title=f'{args.group}; {args.region}; {args.temperature:g} K; selected bands {meta["selected_bands"]}')
    ax.axhline(0, color='0.7', lw=.7); ax.legend()
    for suffix in ('png', 'pdf', 'svg'): fig.savefig(out/('character_hall.'+suffix), dpi=180)
    plt.close(fig)
    result = dict(schema='vaspberry.character-plots', version=1, complete=True, provenance=provenance(args),
                  figure_sha256={p.name: sha256(p) for p in sorted(out.iterdir())})
    json_write(out/'plots.json', result)
    return result


def provenance(args):
    return dict(command={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        implementation_sha256={name: sha256(Path(__file__).with_name(name)) for name in
            ('procar_character.py', 'procar_projection.py', 'kubo_pairs.py', 'berry_data.py', 'vaspberry_transport.py')})


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version', action='version', version=__version__)
    sub = p.add_subparsers(dest='command', required=True)
    project = sub.add_parser('project', help='same-run PROCAR/WAVECAR/OUTCAR -> reusable character CSV/NPZ')
    for key in ('procar', 'wavecar', 'outcar', 'groups'):
        project.add_argument('--'+key, type=Path, required=True)
    project.add_argument('--axis', nargs=3, type=float, required=True, metavar=('X', 'Y', 'Z'), help='Cartesian unit spin-analysis axis')
    hall = sub.add_parser('hall', help='character cache + native pair cache -> selected-band charge-Hall attribution')
    hall.add_argument('--character-dir', type=Path, required=True)
    hall.add_argument('--pairs-dir', type=Path, required=True)
    hall.add_argument('--bands', type=int, nargs='+', required=True, help='1-based individual bands; each must be isolated at every k')
    hall.add_argument('--mu-min', type=float, required=True)
    hall.add_argument('--mu-max', type=float, required=True)
    hall.add_argument('--mu-num', type=int, default=241)
    hall.add_argument('--mu-reference', type=float, required=True)
    hall.add_argument('--temperatures', type=float, nargs='+', default=[0.])
    hall.add_argument('--regions', type=Path, help='same named region JSON as vaspberry_kubo.py pair-hall')
    hall.add_argument('--degeneracy-threshold-eV', type=float, default=1e-5)
    plot = sub.add_parser('plot', help='saved tables -> raw character map and attributed Hall PNG/PDF/SVG')
    plot.add_argument('--character-dir', type=Path, required=True)
    plot.add_argument('--hall-dir', type=Path, required=True)
    plot.add_argument('--group', required=True)
    plot.add_argument('--band', type=int, required=True)
    plot.add_argument('--temperature', type=float, required=True)
    plot.add_argument('--region', default='total')
    plot.add_argument('--delta', action='store_true')
    for cmd in (project, hall, plot): cmd.add_argument('--output-dir', type=Path, required=True)
    return p


def main(argv=None):
    p = parser(); args = p.parse_args(argv)
    try:
        meta = {'project': project_command, 'hall': hall_command, 'plot': plot_command}[args.command](args)
    except (ValueError, OSError, KeyError) as exc:
        p.error(str(exc))
    print(json.dumps({'schema': meta['schema'], 'version': meta['version'], 'output_dir': str(args.output_dir)}, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
