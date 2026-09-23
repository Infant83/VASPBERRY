#!/usr/bin/env python3
"""T=0 insulating sheet Hall response from a standard VASP 5.4.4 WAVEDER run.

Use WAVEDER, WAVECAR, INCAR and OUTCAR from the same completed static run.
Several fixed-charge runs may be supplied when their union is a full mesh;
their retained CHGCAR/POTCAR/POSCAR and physical settings must match.
WAVEDER has no energies or k coordinates: consistency checks cannot authenticate
its association with the other files. The supplied run directory asserts that
association, which is recorded together with automatically computed file hashes.

Only the non-PEAD, longitudinal PAW optical branch is supported. Occupied-to-empty
matrix elements form the occupied-bundle trace; internal degenerate occupied
states never require division. This is not a metallic or finite-temperature path.
"""
from __future__ import annotations

import argparse
from bisect import bisect_right
import json
from pathlib import Path
import re
from types import SimpleNamespace

import numpy as np

from berry_data import require, write_hall
from exported_matrix_kubo import sha256
from kubo_pairs import hall_metadata, integration_setup, result_rows
from vasp_optical_export import read_waveder
from wavecar_fukui import Wavecar

PRODUCER_THRESHOLD_EV = .002
NUMBER = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'
ASSOCIATION = ('The supplied run directory asserts a single unmodified standard VASP run. '
               'OUTCAR/WAVECAR consistency and WAVEDER dimensions are checked; WAVEDER has '
               'no energies, coordinates or producer hashes, so its same-run association '
               'cannot be independently authenticated from these files.')


def number(value):
    return float(value.replace('D', 'E').replace('d', 'e'))


def boolean(value):
    require(value.upper().strip('.') in ('T', 'F', 'TRUE', 'FALSE'), 'invalid VASP logical value')
    return value.upper().strip('.') in ('T', 'TRUE')


def incar_values(text):
    values = {}
    for line in text.splitlines():
        for part in re.split('[;]', re.split('[!#]', line, maxsplit=1)[0]):
            if not part.strip():
                continue
            match = re.fullmatch(r'\s*([A-Za-z_][A-Za-z_0-9]*)\s*=\s*(.*?)\s*', part)
            require(match is not None, 'unsupported/ambiguous INCAR assignment')
            name, value = match.group(1).upper(), match.group(2).strip()
            require(value and (name not in values or values[name] == value), 'conflicting INCAR tag '+name)
            values[name] = value
    return values


def outcar_value(text, name, convert):
    # Advisory boxes may quote alternatives such as "try LREAL= Auto".
    # Read the parameter report, not arbitrary occurrences in prose.
    prefix = r'^\s*k-points[^\n]*?\b' if name in ('NKPTS', 'NBANDS') else r'^\s*'
    tokens = re.findall(prefix+re.escape(name)+r'\s*=\s*([^\s;]+)', text, re.M)
    require(bool(tokens), 'OUTCAR lacks effective '+name)
    values = [convert(token) for token in tokens]
    require(all(value == values[0] for value in values), 'ambiguous OUTCAR '+name)
    return values[0]


def read_run_settings(incar_text, outcar_text):
    """Restrict operator semantics to the audited standard optical branch."""
    versions = re.findall(r'^\s*vasp\.(\d+\.\d+\.\d+)\b', outcar_text, re.M)
    require(versions and set(versions) == {'5.4.4'}, 'only standard VASP 5.4.4 WAVEDER is supported')
    require('General timing and accounting informations for this job:' in outcar_text
            and 'aborting loop because EDIFF is reached' in outcar_text,
            'completed, electronically converged static OUTCAR required')
    require('frequency dependent IMAGINARY DIELECTRIC FUNCTION' in outcar_text
            and re.search(r'\bOPTICS:\s+cpu time', outcar_text), 'OUTCAR lacks completed optical calculation')
    inp = incar_values(incar_text)
    for tag, expected in (('LOPTICS', True), ('LPEAD', False), ('LNABLA', False)):
        require(tag in inp and boolean(inp[tag]) == expected, 'explicit supported INCAR '+tag+' required')
    require(not boolean(inp.get('LBERRY_EXPORT', 'F')), 'use standard WAVEDER output without custom export mode')
    settings = {'version': '5.4.4', 'LOPTICS': True, 'LPEAD': False}
    for name in ('LNABLA', 'LNONCOLLINEAR', 'LSORBIT', 'LREAL', 'LHFCALC', 'METAGGA', 'LEPSILON', 'LVEL'):
        settings[name] = outcar_value(outcar_text, name, boolean)
    for name in ('ISPIN', 'ISYM', 'NSW', 'NKPTS', 'NBANDS'):
        settings[name] = outcar_value(outcar_text, name, int)
    for name in ('NELECT', 'DEG_THRESHOLD'):
        settings[name] = outcar_value(outcar_text, name, number)
    require(not any(settings[name] for name in ('LNABLA', 'LREAL', 'LHFCALC', 'METAGGA', 'LEPSILON', 'LVEL')),
            'unsupported PAW optical branch; require LNABLA/LREAL/LHFCALC/METAGGA/LEPSILON/LVEL false')
    require(settings['ISYM'] == -1 and settings['NSW'] == 0, 'static NSW=0 and full-mesh ISYM=-1 required')
    require(settings['DEG_THRESHOLD'] == PRODUCER_THRESHOLD_EV, 'standard producer DEG_THRESHOLD=0.002 eV required')
    require(settings['ISPIN'] in (1, 2) and settings['NKPTS'] > 0 and settings['NBANDS'] > 1
            and np.isfinite(settings['NELECT']) and settings['NELECT'] > 0, 'invalid OUTCAR dimensions/electron count')
    # Where OUTCAR reports a setting, an explicitly supplied INCAR value must agree.
    for name, effective in settings.items():
        if name not in inp or name in ('version', 'LOPTICS', 'LPEAD'):
            continue
        converted = boolean(inp[name]) if isinstance(effective, bool) else number(inp[name])
        if name == 'NBANDS':
            # VASP may increase a requested NBANDS for MPI band distribution.
            require(converted.is_integer() and 0 < converted <= effective,
                    'INCAR NBANDS exceeds effective OUTCAR band count')
            continue
        require(converted == effective, 'INCAR/OUTCAR setting mismatch: '+name)
    return settings


def outcar_arrays(text, nk, nb, ns):
    """Read final printed eigenstates, preserving OUTCAR rounding tolerances."""
    def lines_at(start, count):
        lines = []
        for _ in range(count):
            end = text.find('\n', start)
            require(end >= 0, 'truncated OUTCAR table')
            lines.append(text[start:end]); start = end+1
        return lines
    header = re.compile(r'^\s*k-point\s+(\d+)\s*:\s*('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+
                        r')\s*\n\s*band No\.\s+band energies\s+occupation\s*\n', re.M)
    matches = list(header.finditer(text))
    require(len(matches) >= nk*ns, 'OUTCAR lacks complete final eigenvalue tables')
    matches = matches[-nk*ns:]
    spins = [(m.start(), int(m.group(1))) for m in re.finditer(r'^\s*spin component\s+(\d+)\s*$', text, re.M)]
    starts = [v[0] for v in spins]
    q = np.empty((ns, nk, 3)); energy = np.empty((ns, nk, nb)); occ = np.empty_like(energy)
    for i, match in enumerate(matches):
        s, k = divmod(i, nk)
        require(int(match.group(1)) == k+1, 'OUTCAR final k-point order/coverage mismatch')
        if ns == 2:
            si = bisect_right(starts, match.start())-1
            require(si >= 0 and spins[si][1] == s+1, 'OUTCAR final spin-channel order mismatch')
        q[s, k] = [number(v) for v in match.groups()[1:]]
        lines = lines_at(match.end(), nb)
        require(len(lines) == nb, 'truncated OUTCAR final bands')
        for b, line in enumerate(lines):
            values = line.split()
            require(len(values) == 3 and values[0] == str(b+1), 'malformed OUTCAR final band row')
            energy[s, k, b], occ[s, k, b] = [number(v) for v in values[1:]]
    lattice_blocks = list(re.finditer(r'direct lattice vectors\s+reciprocal lattice vectors\s*\n', text))
    require(bool(lattice_blocks), 'OUTCAR lattice required')
    lattice_rows = lines_at(lattice_blocks[-1].end(), 3)
    lattice = np.array([[number(v) for v in row.split()[:3]] for row in lattice_rows])
    grid_blocks = list(re.finditer(r'k-points in reciprocal lattice and weights:[^\n]*\n', text))
    require(bool(grid_blocks), 'OUTCAR reciprocal k-point list required')
    grid_lines = lines_at(grid_blocks[-1].end(), nk)
    grid = np.array([[number(v) for v in row.split()] for row in grid_lines])
    require(lattice.shape == (3, 3) and grid.shape == (nk, 4), 'incomplete OUTCAR lattice/k-point table')
    require(all(np.isfinite(v).all() for v in (q, energy, occ, lattice, grid)), 'nonfinite OUTCAR data')
    return q, energy, occ, lattice, grid


def occupied_curvature(connection, occupied):
    """Trace from lower occupied-to-empty blocks; zeros remain valid zeros."""
    require(connection.ndim == 4 and connection.shape[1] == 3
            and 0 < occupied < connection.shape[2] and occupied <= connection.shape[3],
            'WAVEDER must cover every occupied ket and stored empty bra')
    c = np.asarray(connection[:, :, occupied:, :occupied], dtype=np.complex128)
    require(np.isfinite(c).all(), 'nonfinite occupied-empty WAVEDER block')
    return np.stack([-2*np.imag(np.sum(c[:, a].conj()*c[:, b], axis=(1, 2)))
                     for a, b in ((1, 2), (2, 0), (0, 1))], axis=1)


def read_validated_run(run_dir, mus, *, occupied, spin, spinor_components, spin_multiplicity,
                       energy_reference, mu_reference):
    """Read genuine per-run states and optical curvature without a mesh assumption."""
    run = Path(run_dir)
    paths = {name: run/name for name in ('WAVEDER', 'WAVECAR', 'INCAR', 'OUTCAR')}
    require(all(p.is_file() for p in paths.values()), 'run directory requires WAVEDER, WAVECAR, INCAR and OUTCAR')
    hashes = {name: sha256(p) for name, p in paths.items()}
    text = paths['OUTCAR'].read_text()
    settings = read_run_settings(paths['INCAR'].read_text(), text)
    expected_spinors = 2 if settings['LNONCOLLINEAR'] else 1
    require(spinor_components == expected_spinors and (not settings['LSORBIT'] or expected_spinors == 2),
            'spinor declaration disagrees with OUTCAR')
    require(type(spin) is int and 1 <= spin <= settings['ISPIN'], 'selected spin channel outside OUTCAR')
    require(type(occupied) is int and 0 < occupied < settings['NBANDS'], 'occupied leading bundle and stored empty bands required')
    require(isinstance(energy_reference, str) and bool(energy_reference.strip()), 'describe the unchanged VASP energy reference')
    ns, nk, nb = (settings[name] for name in ('ISPIN', 'NKPTS', 'NBANDS'))
    printed_q, printed_e, printed_occ, lattice, grid = outcar_arrays(text, nk, nb, ns)
    require(np.allclose(grid[:, 3], 1/nk, atol=.00051, rtol=0), 'OUTCAR full-grid weights are not uniform')
    waves = []
    physical_spin_factor = 2 if ns == 1 and spinor_components == 1 else 1
    for channel in range(1, ns+1):
        w = Wavecar(paths['WAVECAR'], spin=channel, spinor_components=spinor_components)
        w.coefficients(0, [1])
        require(w.header.ispin == ns and w.energies.shape == (nk, nb), 'OUTCAR/WAVECAR dimensions disagree')
        require(np.all(np.diff(w.energies, axis=1) >= -1e-8), 'energy-ordered bands required')
        require(np.allclose(w.header.lattice, lattice, atol=6e-9, rtol=0), 'OUTCAR/WAVECAR lattice mismatch')
        require(np.allclose(w.kpoints, grid[:, :3], atol=6e-9, rtol=0)
                and np.allclose(w.kpoints, printed_q[channel-1], atol=5.1e-5, rtol=0),
                'OUTCAR/WAVECAR k-point mismatch')
        require(np.allclose(w.energies, printed_e[channel-1], atol=5.1e-5, rtol=0),
                'OUTCAR/WAVECAR final energies mismatch; no energy shifts are inferred')
        require(np.allclose(w.occupations*physical_spin_factor, printed_occ[channel-1], atol=5.1e-6, rtol=0),
                'OUTCAR/WAVECAR occupations mismatch')
        waves.append(w)
    w = waves[spin-1]
    source = read_waveder(paths['WAVEDER'])
    require(source.C_A.shape[:4] == (ns, nk, 3, nb), 'WAVEDER/WAVECAR spin/k/band dimensions disagree')
    require(source.C_A.shape[4] >= occupied, 'WAVEDER occupied-ket coverage is incomplete')
    require(source.source_sha256 == hashes['WAVEDER'], 'WAVEDER changed while reading')
    requested = np.asarray([*mus, mu_reference], dtype=float)
    require(np.isfinite(requested).all(), 'finite chemical potentials required')
    fillings = []; vbms = []; cbms = []
    for channel, state in enumerate(waves):
        counts = np.sum(state.energies < mu_reference, axis=1)
        require(np.all(counts == counts[0]) and 0 < counts[0] < nb,
                'T=0 reference must give a fixed occupied count and stored empty states in every spin channel')
        count = int(counts[0]); fillings.append(count)
        if channel == spin-1:
            require(count == occupied, 'selected occupied count disagrees with T=0 reference filling')
        vbms.append(float(state.energies[:, :count].max()))
        cbms.append(float(state.energies[:, count:].min()))
    require(abs(sum(fillings)*physical_spin_factor-settings['NELECT']) < 5.1e-5,
            'T=0 insulating band filling disagrees with NELECT')
    vbm, cbm = max(vbms), min(cbms)
    require(cbm > vbm and np.all((requested > vbm) & (requested < cbm)),
            'T=0 WAVEDER scan requires mu strictly inside the common global insulating gap')
    cross_gap = w.energies[:, occupied]-w.energies[:, occupied-1]
    require(np.all(cross_gap > PRODUCER_THRESHOLD_EV), 'occupied-empty separation must exceed producer 0.002 eV cluster threshold')
    omega = occupied_curvature(source.C_A[spin-1], occupied)
    metadata = dict(sampling={'kind': 'points'}, spin_multiplicity=spin_multiplicity, energy_reference=energy_reference,
        source_nspin=ns, spinor_components=spinor_components, spin_channel_1based=spin,
        source_nbands=nb, source_ndbands=source.C_A.shape[4], source_precision='complex64 WAVEDER; complex128 accumulation',
        source_operator={'kind': 'vasp_longitudinal_paw_optical_occupied_empty',
                         'accuracy_status': 'supported_operator_scope; convergence_not_established',
                         'coverage': 'selected occupied kets to all stored empty bras; no internal occupied terms'},
        producer_settings=settings, producer_cluster_threshold_eV=PRODUCER_THRESHOLD_EV,
        same_run_association={'status': 'user_supplied_run_directory_consistency_checked', 'limitation': ASSOCIATION},
        source_sha256=hashes, validation_tolerances={'OUTCAR_energy_eV': 5.1e-5, 'OUTCAR_occupation': 5.1e-6},
        source_occupations='Consistency checked including smearing; T=0 response uses explicit insulating filling.',
        physical_spin_factor_for_NELECT=physical_spin_factor, occupied_counts_by_spin=fillings,
        run_directory=str(run.resolve()))
    data = SimpleNamespace(metadata=metadata, k_ids=np.arange(1, nk+1), kpoints_fractional=w.kpoints.copy(),
        weights=np.full(nk, 1/nk), lattice_A=w.header.lattice.copy(), reciprocal_inv_A=w.header.reciprocal.copy(),
        omega_A2=omega, energies_eV=w.energies.copy(), valence_max_eV=vbm, conduction_min_eV=cbm,
        minimum_cross_gap_eV=float(cross_gap.min()))
    require(all(sha256(p) == hashes[name] for name, p in paths.items()), 'source files changed during validation')
    return data


def shared_chunk_inputs(runs):
    """Bind chunk Hamiltonians through retained, identical fixed-charge inputs."""
    records = []; first_parameters = None
    for run in runs:
        validation_hashes = {name: sha256(run/name) for name in ('INCAR', 'OUTCAR')}
        inp = incar_values((run/'INCAR').read_text())
        text = (run/'OUTCAR').read_text()
        require(inp.get('ICHARG') is not None and number(inp['ICHARG']) == 11
                and 'LCHARG' in inp and not boolean(inp['LCHARG'])
                and outcar_value(text, 'ICHARG', int) == 11
                and not outcar_value(text, 'LCHARG', boolean),
                'multiple runs require effective ICHARG=11 and LCHARG=false in every run')
        parameters = {key: value for key, value in inp.items() if key != 'SYSTEM'}
        require(first_parameters is None or parameters == first_parameters,
                'multiple runs require identical INCAR parameters except SYSTEM')
        first_parameters = parameters
        paths = {name: run/name for name in ('CHGCAR', 'POTCAR', 'POSCAR')}
        require(all(p.is_file() and p.stat().st_size > 0 for p in paths.values()),
                'multiple runs require retained CHGCAR, POTCAR and POSCAR inputs')
        hashes = {name: sha256(p) for name, p in paths.items()}
        require(not records or hashes == records[0]['input_sha256'],
                'multiple runs have different CHGCAR/POTCAR/POSCAR inputs; do not mix Hamiltonians')
        records.append({'run_directory': str(run.resolve()), 'input_sha256': hashes,
                        'validation_record_sha256': validation_hashes})
    return records


def waveder_hall_spectrum(run_dir, mus, *, occupied, spin, spinor_components, spin_multiplicity,
                         sampling, energy_reference, mu_reference, region_spec=None, differences=()):
    """Integrate one full run or a validated union of genuine fixed-charge runs."""
    runs = [Path(run_dir)] if isinstance(run_dir, (str, Path)) else [Path(v) for v in run_dir]
    require(bool(runs) and len({p.resolve() for p in runs}) == len(runs), 'distinct nonempty run directories required')
    mus = np.asarray(mus, dtype=float)
    require(mus.ndim == 1 and len(mus) and np.isfinite(mus).all(), 'finite chemical-potential array required')
    chunk_inputs = shared_chunk_inputs(runs) if len(runs) > 1 else None
    chunks = [read_validated_run(run, mus, occupied=occupied, spin=spin, spinor_components=spinor_components,
        spin_multiplicity=spin_multiplicity, energy_reference=energy_reference, mu_reference=mu_reference) for run in runs]
    first = chunks[0]
    common_settings = {k: v for k, v in first.metadata['producer_settings'].items() if k != 'NKPTS'}
    for data in chunks[1:]:
        require({k: v for k, v in data.metadata['producer_settings'].items() if k != 'NKPTS'} == common_settings,
                'chunk effective settings, spin channels or band counts disagree')
        require(data.metadata['occupied_counts_by_spin'] == first.metadata['occupied_counts_by_spin'],
                'chunk T=0 spin fillings disagree')
        require(np.allclose(data.lattice_A, first.lattice_A, atol=1e-9, rtol=0)
                and np.allclose(data.reciprocal_inv_A, first.reciprocal_inv_A, atol=1e-9, rtol=0),
                'chunk lattice/reciprocal geometry mismatch')
    q = np.concatenate([c.kpoints_fractional for c in chunks]); nk = len(q)
    omega = np.concatenate([c.omega_A2 for c in chunks])
    metadata = dict(first.metadata, sampling=sampling, source_nkpoints=nk, source_run_count=len(runs))
    if len(runs) > 1:
        metadata.pop('source_sha256'); metadata.pop('run_directory')
        metadata.pop('source_ndbands')
        metadata['source_ndbands_by_run'] = [c.metadata['source_ndbands'] for c in chunks]
        metadata['producer_settings'] = dict(common_settings, NKPTS=nk)
        metadata['producer_settings_scope'] = 'common effective settings; NKPTS is the validated union size'
        metadata['source_runs'] = [c.metadata for c in chunks]
        metadata['chunk_input_records'] = chunk_inputs
        metadata['chunk_combination'] = ('Concatenated genuine per-run energies, coordinates and optical curvatures; '
            'no synthetic WAVEDER/OUTCAR; uniform weights assigned only after full union mesh validation.')
        metadata['same_run_association'] = {'status': 'multiple_fixed_charge_runs_consistency_checked',
            'limitation': ASSOCIATION+' The runs also share retained fixed-charge inputs and INCAR parameters.'}
    data = SimpleNamespace(metadata=metadata, k_ids=np.arange(1, nk+1), kpoints_fractional=q,
        weights=np.full(nk, 1/nk), lattice_A=first.lattice_A, reciprocal_inv_A=first.reciprocal_inv_A,
        energies_eV=np.concatenate([c.energies_eV for c in chunks]))
    vbm = max(c.valence_max_eV for c in chunks); cbm = min(c.conduction_min_eV for c in chunks)
    mus, ts, area, normal, regions = integration_setup(data, mus, [0.], mu_reference, region_spec, differences)
    all_mus = np.unique(np.r_[mus, mu_reference]); values = np.zeros((1, len(regions), len(all_mus), 4))
    for ri, mask in enumerate(regions.values()):
        values[0, ri, :, 0] = -area/(2*np.pi)*spin_multiplicity*np.sum(data.weights*(omega@normal)*mask)
        values[0, ri, :, 1] = occupied*spin_multiplicity*np.sum(data.weights*mask)
    rows = result_rows(values, mus, ts, all_mus, regions, mu_reference)
    meta = hall_metadata(data, area, normal, region_spec, differences, regions, mu_reference)
    meta.update(method='standard_waveder_occupied_bundle_T0', scope='fixed_occupied_bundle',
        occupied_bands=[1, occupied], valence_max_eV=vbm, conduction_min_eV=cbm,
        global_gap_eV=cbm-vbm, minimum_selected_cross_gap_eV=min(c.minimum_cross_gap_eV for c in chunks),
        formula='Omega_ab=-2 Im sum_v_occ,c_empty conj(C_a_cv)*C_b_cv; sigma/(e^2/h)=-A_BZ/(2*pi)*g_s*sum_k w_k Omega_normal',
        delta_formula='zero inside the same global insulating gap at T=0',
        integer_rounding_applied=False, zero_matrix_elements='valid selection-rule zeros; not a coverage mask',
        references=['https://vasp.at/wiki/WAVEDER', 'https://vasp.at/wiki/LPEAD', 'https://vasp.at/wiki/LNABLA'])
    if chunk_inputs is not None:
        for run, chunk, record in zip(runs, chunks, chunk_inputs):
            expected = {**chunk.metadata['source_sha256'], **record['input_sha256'],
                        **record['validation_record_sha256']}
            require(all(sha256(run/name) == value for name, value in expected.items()),
                    'chunk source files changed during validation')
    return rows, meta


def add_arguments(p):
    """Shared argument contract for the standalone and composable commands."""
    p.add_argument('--run-dir', type=Path, nargs='+', required=True,
                   help='one full optical run, or fixed-charge chunks whose union is the full mesh')
    p.add_argument('--occupied', type=int, required=True)
    p.add_argument('--spin', type=int, default=1)
    p.add_argument('--spinor-components', type=int, choices=[1, 2], required=True)
    p.add_argument('--spin-multiplicity', type=int, choices=[1, 2], required=True)
    p.add_argument('--mesh', type=int, nargs=2, required=True)
    p.add_argument('--plane-axes', type=int, nargs=2, default=[0, 1])
    p.add_argument('--energy-reference', required=True)
    p.add_argument('--mu-min', type=float, required=True); p.add_argument('--mu-max', type=float, required=True)
    p.add_argument('--mu-num', type=int, default=241); p.add_argument('--mu-reference', type=float, required=True)
    p.add_argument('--regions', type=Path); p.add_argument('--difference', action='append', default=[])
    p.add_argument('--formats', nargs='+', choices=['csv', 'dat', 'npz'], default=['csv', 'dat', 'npz'])
    p.add_argument('--output-dir', type=Path, required=True)
    return p


def add_command(subparsers):
    return add_arguments(subparsers.add_parser('waveder-hall',
        help='standard VASP 5.4.4 optical run -> insulating T=0 PAW sheet Hall', description=__doc__))


def command(args):
    require(not args.output_dir.exists(), 'output directory exists; choose a new directory')
    require(args.mu_num >= 1 and np.isfinite([args.mu_min, args.mu_max]).all() and args.mu_max >= args.mu_min
            and (args.mu_num > 1 or args.mu_min == args.mu_max), 'valid chemical-potential scan required')
    differences = [tuple(value.split(':')) for value in args.difference]
    require(all(len(v) == 3 for v in differences), 'difference syntax is NAME:LEFT:RIGHT')
    rows, meta = waveder_hall_spectrum(args.run_dir, np.unique(np.linspace(args.mu_min, args.mu_max, args.mu_num)),
        occupied=args.occupied, spin=args.spin, spinor_components=args.spinor_components,
        spin_multiplicity=args.spin_multiplicity, sampling={'kind': 'uniform_full_2d', 'mesh': args.mesh,
        'plane_axes': args.plane_axes}, energy_reference=args.energy_reference, mu_reference=args.mu_reference,
        region_spec=json.loads(args.regions.read_text()) if args.regions else None, differences=differences)
    meta['provenance'] = {'command': {k: [str(p) for p in v] if k == 'run_dir' else
                                    str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                          'adapter_sha256': sha256(Path(__file__))}
    return write_hall(args.output_dir, rows, meta, formats=args.formats)


def main(argv=None):
    p = add_arguments(argparse.ArgumentParser(description=__doc__))
    args = p.parse_args(argv)
    try:
        meta = command(args)
    except (ValueError, OSError, KeyError, TypeError) as exc:
        p.error(str(exc))
    print(json.dumps({'output': str(args.output_dir), 'schema': meta['schema'], 'version': meta['version']}))


if __name__ == '__main__':
    main()
