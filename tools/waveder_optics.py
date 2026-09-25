#!/usr/bin/env python3
"""Independent-particle circular transition strengths from standard VASP WAVEDER.

Reads one completed, unmodified VASP 5.4.4 longitudinal PAW optical run.
Outputs k-resolved strengths and Gaussian spectral densities, not absolute
absorption or photoluminescence. A path is allowed; no BZ integral is inferred.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from berry_data import require
from exported_matrix_kubo import sha256
from waveder_hall import ASSOCIATION, PRODUCER_THRESHOLD_EV, read_validated_optical_states


def polarization_frame(beam, axis):
    """Right-handed real frame (e1, e2, beam); e2 = beam cross e1."""
    beam, axis = np.asarray(beam, float), np.asarray(axis, float)
    require(beam.shape == axis.shape == (3,) and np.isfinite([beam, axis]).all(),
            'finite Cartesian beam and polarization-axis vectors required')
    require(np.linalg.norm(beam) > 0 and np.linalg.norm(axis) > 0, 'nonzero polarization vectors required')
    beam = beam/np.linalg.norm(beam); axis = axis/np.linalg.norm(axis)
    e1 = axis-np.dot(axis, beam)*beam
    require(np.linalg.norm(e1) > 1e-10, 'polarization axis must not be parallel to beam')
    e1 /= np.linalg.norm(e1)
    return e1, np.cross(beam, e1), beam


def complete_groups(energy, selected, threshold):
    """Transitive adjacent-energy clusters on the complete stored band axis."""
    energy = np.asarray(energy, float)
    require(energy.ndim == 1 and len(energy) > 1 and np.isfinite(energy).all()
            and np.all(np.diff(energy) >= -1e-8), 'finite energy-ordered bands required')
    lo, hi = selected
    require(type(lo) is int and type(hi) is int and 1 <= lo <= hi <= len(energy),
            'one-based inclusive band range outside source states')
    cuts = np.r_[0, np.flatnonzero(np.diff(energy) > threshold)+1, len(energy)]
    require(lo-1 in cuts and hi in cuts, 'band range cuts a degenerate group; include all partners')
    return [np.arange(a, b) for a, b in zip(cuts[:-1], cuts[1:]) if lo-1 <= a < b <= hi]


def circular_strength(connection, e1, e2):
    """Sum complete final-bra/initial-ket blocks; no complex conjugation of epsilon.

    C has Cartesian, final, initial axes. Its common phase relative to the
    electric-dipole matrix cancels in |epsilon dot C|^2.
    """
    c = np.asarray(connection, complex)
    require(c.ndim == 3 and c.shape[0] == 3 and np.isfinite(c).all(), 'finite Cartesian transition block required')
    x, y = np.einsum('a,aij->ij', e1, c), np.einsum('a,aij->ij', e2, c)
    plus = float(np.sum(np.abs((x+1j*y)/np.sqrt(2))**2))
    minus = float(np.sum(np.abs((x-1j*y)/np.sqrt(2))**2))
    return plus, minus


def contrast(plus, minus, relative_floor):
    total = plus+minus
    valid = total > relative_floor*total.max(initial=0.)
    eta = np.full(total.shape, np.nan)
    np.divide(plus-minus, total, out=eta, where=valid)
    return eta, valid


def transition_tables(energies, connection, kpoints, weights, *, occupied, initial, final,
                      photon_eV, sigma_eV, beam=(0, 0, 1), axis=(1, 0, 0),
                      degeneracy_threshold_eV=PRODUCER_THRESHOLD_EV, relative_intensity_floor=1e-10):
    """Group-summed point spectra with centroid energies for each near-degenerate group.

    A stored empty band beyond the selected final range is required, so the
    final range cannot silently end inside an unobserved degenerate multiplet.
    Adjacent splittings <= threshold are unresolved and represented by their
    group mean. Energy ranges and maximum group spread are retained.
    """
    e, c = np.asarray(energies, float), np.asarray(connection)
    q, w = np.asarray(kpoints, float), np.asarray(weights, float)
    photons = np.asarray(photon_eV, float)
    require(e.ndim == 2 and e.shape[1] > 2 and np.isfinite(e).all(), 'finite k/band energies required')
    nk, nb = e.shape
    require(c.ndim == 4 and c.shape[:3] == (nk, 3, nb), 'connection axes must be k, Cartesian, final bra, initial ket')
    require(q.shape == (nk, 3) and w.shape == (nk,) and np.isfinite(q).all()
            and np.isfinite(w).all() and np.all(w >= 0), 'finite point coordinates and nonnegative source weights required')
    require(type(occupied) is int and 0 < occupied < nb and c.shape[3] >= occupied,
            'occupied leading bundle and complete WAVEDER occupied-ket coverage required')
    require(len(initial) == len(final) == 2 and 1 <= initial[0] <= initial[1] <= occupied
            and occupied < final[0] <= final[1] < nb,
            'select occupied initial and empty final ranges, retaining a source band above the final range')
    require(photons.ndim == 1 and len(photons) and np.isfinite(photons).all()
            and np.all(photons >= 0) and np.all(np.diff(photons) > 0), 'increasing nonnegative photon-energy grid required')
    require(np.isfinite(sigma_eV) and sigma_eV > 0, 'positive finite Gaussian sigma required')
    require(np.isfinite(degeneracy_threshold_eV) and degeneracy_threshold_eV >= PRODUCER_THRESHOLD_EV,
            'degeneracy threshold must be at least the WAVEDER producer threshold of 0.002 eV')
    require(np.isfinite(relative_intensity_floor) and 0 <= relative_intensity_floor < 1,
            'relative intensity floor must lie in [0,1)')
    require(np.all(e[:, occupied]-e[:, occupied-1] > degeneracy_threshold_eV),
            'occupied/empty boundary must separate complete degenerate groups')
    e1, e2, beam = polarization_frame(beam, axis)
    records = []; plus = np.zeros((nk, len(photons))); minus = np.zeros_like(plus)
    max_spread = 0.
    for k in range(nk):
        before = len(records)
        ig = complete_groups(e[k], initial, degeneracy_threshold_eV)
        fg = complete_groups(e[k], final, degeneracy_threshold_eV)
        for i in ig:
            for f in fg:
                block = c[k][:, f][:, :, i]
                p, m = circular_strength(block, e1, e2)
                de = float(e[k, f].mean()-e[k, i].mean())
                gap_min, gap_max = float(e[k, f].min()-e[k, i].max()), float(e[k, f].max()-e[k, i].min())
                max_spread = max(max_spread, float(np.ptp(e[k, f])), float(np.ptp(e[k, i])))
                gaussian = np.exp(-.5*((photons-de)/sigma_eV)**2)/(np.sqrt(2*np.pi)*sigma_eV)
                plus[k] += p*gaussian; minus[k] += m*gaussian
                records.append((k+1, *q[k], w[k], int(i[0]+1), int(i[-1]+1), int(f[0]+1), int(f[-1]+1),
                                de, gap_min, gap_max, p, m))
        local = np.asarray(records[before:])
        eta, valid = contrast(local[:, -2], local[:, -1], relative_intensity_floor)
        records[before:] = [(*row, float(v), bool(ok)) for row, v, ok in zip(records[before:], eta, valid)]
    names = ['k_index', 'kx_frac', 'ky_frac', 'kz_frac', 'source_k_weight',
             'initial_first', 'initial_last', 'final_first', 'final_last',
             'transition_eV', 'transition_min_eV', 'transition_max_eV', 'I_plus_A2', 'I_minus_A2', 'eta', 'eta_valid']
    integer = {'k_index', 'initial_first', 'initial_last', 'final_first', 'final_last'}
    transitions = {key: np.array([row[j] for row in records], dtype=bool if key == 'eta_valid' else int if key in integer else float)
                   for j, key in enumerate(names)}
    eta = np.empty_like(plus); valid = np.empty_like(plus, dtype=bool)
    for k in range(nk):
        eta[k], valid[k] = contrast(plus[k], minus[k], relative_intensity_floor)
    spectra = {'k_index': np.repeat(np.arange(1, nk+1), len(photons))}
    spectra.update({f'k{a}_frac': np.repeat(q[:, j], len(photons)) for j, a in enumerate('xyz')})
    spectra.update(source_k_weight=np.repeat(w, len(photons)), photon_eV=np.tile(photons, nk),
        I_plus_A2_per_eV=plus.ravel(), I_minus_A2_per_eV=minus.ravel(), eta=eta.ravel(), eta_valid=valid.ravel())
    meta = dict(schema='vaspberry.circular-transition-strength', version=1, complete=True,
        scope='independent_particle_k_resolved_occupied_to_empty',
        units={'transition_strength': 'Angstrom^2', 'spectral_density': 'Angstrom^2/eV', 'energy': 'eV'},
        polarization={'beam_cartesian': beam.tolist(), 'e1_cartesian': e1.tolist(), 'e2_cartesian': e2.tolist(),
            'epsilon_plus': '(e1+i*e2)/sqrt(2)', 'epsilon_minus': '(e1-i*e2)/sqrt(2)',
            'field_convention': 'Re[epsilon * exp(i*(q.r-omega*t))]; q parallel to beam',
            'matrix_contraction': 'abs(sum_a epsilon_a C_a[final_bra,initial_ket]))^2; epsilon is not conjugated',
            'handedness_labels': 'plus/minus defined algebraically; no observer-dependent left/right label'},
        initial_bands_1based=list(initial), final_bands_1based=list(final), occupied_count=occupied,
        degeneracy={'threshold_eV': degeneracy_threshold_eV, 'source_top_band_retained_as_boundary_guard': True,
            'grouping': 'transitive adjacent gaps on all stored states; both selected boundaries must close groups',
            'transition_energy': 'difference of unweighted group mean energies; sum strengths before eta',
            'maximum_group_energy_spread_eV': max_spread},
        broadening={'kind': 'normalized_Gaussian', 'sigma_eV': sigma_eV,
            'normalization': 'integral over the whole energy axis is one; finite output grid is not renormalized'},
        eta_mask={'relative_intensity_floor': relative_intensity_floor,
            'rule': 'I_plus+I_minus > floor * maximum total intensity at the same k; zero total always invalid',
            'scales': 'transitions and spectra use their own per-k maxima', 'invalid_eta': 'NaN with eta_valid=false'},
        weighting={'k_weight_applied': False, 'source_weights': 'OUTCAR printed values, retained without normalization',
            'spin_multiplicity_applied': 1, 'occupation_factors': 'one occupied-to-empty channel, fixed integer filling'},
        limitations=['Not absolute absorption, dielectric response or photoluminescence polarization.',
            'No excitons, scattering, recombination, nonequilibrium populations or BZ integration.',
            'Source NBANDS, selected final window, k sampling and Gaussian width require convergence.',
            'Near-degenerate groups use centroid energies within the recorded threshold.',
            'One explicit spin channel; a scalar spin-degenerate total would require a factor of two.'])
    return transitions, spectra, meta


def waveder_optics(run_dir, *, occupied, spin, spinor_components, energy_reference, **options):
    state = read_validated_optical_states(run_dir, spin=spin, spinor_components=spinor_components,
                                         energy_reference=energy_reference)
    fillings = []
    for wave in state.waves:
        counts = np.sum(wave.occupations > .5, axis=1)
        require(np.all(counts == counts[0]) and 0 < counts[0] < wave.energies.shape[1],
                'fixed integer occupied/empty partition required; metallic occupations are unsupported')
        count = int(counts[0]); fillings.append(count)
        require(np.all(wave.occupations[:, :count] > .5) and np.all(wave.occupations[:, count:] < .5),
                'occupied states must form the leading band bundle')
    require(fillings[spin-1] == occupied and abs(sum(fillings)*state.physical_spin_factor-state.settings['NELECT']) < 5.1e-5,
            'occupied count and source spin fillings must agree with NELECT')
    transitions, spectra, meta = transition_tables(state.wave.energies, state.source.C_A[spin-1],
        state.wave.kpoints, state.grid[:, 3], occupied=occupied, **options)
    meta.update(energy_reference=energy_reference,
        source_operator={'kind': 'vasp_longitudinal_paw_optical_occupied_empty',
            'accuracy_status': 'supported_operator_scope; convergence_not_established',
            'matrix_axes': ['k', 'cartesian', 'final_bra', 'initial_ket'],
            'source_conversion': 'standard lower occupied-to-empty WAVEDER C, in Angstrom; no second energy denominator'},
        source_nbands=state.settings['NBANDS'], source_nkpoints=state.settings['NKPTS'],
        source_ndbands=state.source.C_A.shape[4], source_nspin=state.settings['ISPIN'],
        spin_channel_1based=spin, spinor_components=spinor_components,
        source_precision='complex64 WAVEDER; complex128 accumulation', producer_settings=state.settings,
        occupied_counts_by_spin=fillings, lattice_A=state.wave.header.lattice.tolist(),
        reciprocal_inv_A=state.wave.header.reciprocal.tolist(), reciprocal_convention='2pi',
        same_run_association={'status': 'user_supplied_run_directory_consistency_checked', 'limitation': ASSOCIATION},
        source_sha256=state.hashes, run_directory=str(state.run.resolve()))
    require(all(sha256(p) == state.hashes[name] for name, p in state.paths.items()), 'source files changed during optical calculation')
    return transitions, spectra, meta


def write_optics(directory, transitions, spectra, metadata, formats=('csv', 'dat', 'npz')):
    formats = tuple(formats)
    require(formats and len(set(formats)) == len(formats) and set(formats) <= {'csv', 'dat', 'npz'},
            'choose distinct CSV, DAT or NPZ output formats')
    out = Path(directory)
    require(not out.exists(), 'output directory exists; choose a new directory')
    out.mkdir(parents=True)
    hashes = {}
    for name, table in [('transitions', transitions), ('spectra', spectra)]:
        for fmt in formats:
            path = out/f'{name}.{fmt}'
            if fmt == 'npz':
                np.savez_compressed(path, **table)
            else:
                with path.open('w', newline='') as stream:
                    if fmt == 'dat':
                        stream.write('# ')
                    writer = csv.writer(stream, delimiter='\t' if fmt == 'dat' else ',')
                    writer.writerow(table); writer.writerows(zip(*table.values()))
            hashes[path.name] = sha256(path)
    meta = dict(metadata, output_formats=list(formats), output_sha256=hashes)
    (out/'optical.json').write_text(json.dumps(meta, indent=2, allow_nan=False)+'\n')
    return meta


def add_arguments(p):
    p.add_argument('--run-dir', type=Path, required=True)
    p.add_argument('--occupied', type=int, required=True)
    p.add_argument('--spin', type=int, default=1)
    p.add_argument('--spinor-components', type=int, choices=[1, 2], required=True)
    p.add_argument('--energy-reference', required=True)
    p.add_argument('--initial', type=int, nargs=2, metavar=('FIRST', 'LAST'), required=True)
    p.add_argument('--final', type=int, nargs=2, metavar=('FIRST', 'LAST'), required=True,
                   help='complete empty groups, with at least one stored band above LAST')
    p.add_argument('--beam-vector', type=float, nargs=3, default=[0, 0, 1])
    p.add_argument('--polarization-axis', type=float, nargs=3, default=[1, 0, 0])
    p.add_argument('--photon-min', type=float, required=True); p.add_argument('--photon-max', type=float, required=True)
    p.add_argument('--photon-num', type=int, default=401)
    p.add_argument('--sigma-eV', type=float, default=.05, help='Gaussian standard deviation in eV')
    p.add_argument('--degeneracy-threshold-eV', type=float, default=PRODUCER_THRESHOLD_EV)
    p.add_argument('--relative-intensity-floor', type=float, default=1e-10)
    p.add_argument('--formats', nargs='+', choices=['csv', 'dat', 'npz'], default=['csv', 'dat', 'npz'])
    p.add_argument('--output-dir', type=Path, required=True)
    return p


def add_command(subparsers):
    return add_arguments(subparsers.add_parser('waveder-optics', help='standard PAW optical run -> circular transition strengths', description=__doc__))


def command(args):
    require(not args.output_dir.exists(), 'output directory exists; choose a new directory')
    require(args.photon_num >= 1 and np.isfinite([args.photon_min, args.photon_max]).all()
            and 0 <= args.photon_min <= args.photon_max
            and (args.photon_num > 1 or args.photon_min == args.photon_max)
            and (args.photon_min < args.photon_max or args.photon_num == 1), 'valid photon range and sample count required')
    transitions, spectra, meta = waveder_optics(args.run_dir, occupied=args.occupied, spin=args.spin,
        spinor_components=args.spinor_components, energy_reference=args.energy_reference, initial=args.initial,
        final=args.final, photon_eV=np.linspace(args.photon_min, args.photon_max, args.photon_num),
        sigma_eV=args.sigma_eV, beam=args.beam_vector, axis=args.polarization_axis,
        degeneracy_threshold_eV=args.degeneracy_threshold_eV, relative_intensity_floor=args.relative_intensity_floor)
    meta['provenance'] = {'command': {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        'implementation_sha256': {name: sha256(Path(__file__).with_name(name))
            for name in ('waveder_optics.py', 'waveder_hall.py', 'vasp_optical_export.py')}}
    return write_optics(args.output_dir, transitions, spectra, meta, formats=args.formats)


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
