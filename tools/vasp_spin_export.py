"""Validate and import an opt-in VASP PAW spin/velocity producer stream.

This is an original binary reader and validation layer, not VASP source.
The supported stream preserves complete band matrices, including diagonal
and degenerate velocity blocks. Its producer recipe is documented separately.
All Pauli matrices are dimensionless; the velocity matrix is hbar*v in eV A.
No normalization, Hermitian averaging, or missing-block filling is performed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import struct

import numpy as np

from berry_data import require
from exported_matrix_kubo import sha256
from spin_operators import pauli_from_coefficients
from vasp_optical_export import read_connection
from waveder_hall import boolean, incar_values, outcar_value
from wavecar_fukui import Wavecar

MAGIC = b'VBERRY_PAW_SPIN_VELOCITY_V1'.ljust(32, b' ')
FOOTER = b'VBERRY_SPIN_VELOCITY_COMPLETE_V1'.ljust(32, b' ')
MAX_STREAM_BYTES = 2 * 1024**3
SCOPE = ('VASP 5.4.4 serial SOC spinors, spin-independent PAW metric, '
         'Cartesian +z SAXIS, semilocal functional without Hubbard U, '
         'no spin spiral, LREAL/LNABLA/LPEAD false; '
         'complete predivision commutator plus projector derivative and '
         'one-center dipole correction. No Hermitian repair.')
LIMITATIONS = [
    'Spin current uses the anticommutator of finite source-band projected spin and velocity; source-band closure must converge.',
    'Finite plane-wave and PAW partial-wave basis; semilocal functional and spin-independent PAW metric only.',
    'Exporter consistency checks do not establish k-mesh, DFT, or source-band convergence.',
]


def read_spin_velocity(path):
    """Read the bounded little-endian, Fortran-order producer version 1.

    Arrays use [k, Cartesian, bra, ket]. ``pseudo`` and ``delta`` instead
    have four components: overlap, Pauli x, y, z. ``xi_restored`` already
    includes the removed energy derivative times the full PAW metric.
    """
    path = Path(path)
    size = path.stat().st_size
    require(244 <= size <= MAX_STREAM_BYTES, 'spin-velocity stream size is unsupported')
    with path.open('rb') as f:
        require(f.read(32) == MAGIC, 'invalid spin-velocity magic')
        version, nk, nb, ns, nc, rb, cb = struct.unpack('<7i', f.read(28))
        require((version, ns, nc, rb, cb) == (1, 1, 3, 8, 16)
                and nk > 0 and nb > 0, 'unsupported spin-velocity header')
        expected = 244 + 32*nk + 16*nk*nb + 17*16*nb*nb*nk
        require(size == expected, 'spin-velocity dimensions or byte length disagree')

        def take(shape, dtype='<f8'):
            count = int(np.prod(shape))
            a = np.fromfile(f, dtype=dtype, count=count)
            require(a.size == count and np.isfinite(a).all(), 'truncated or nonfinite spin-velocity array')
            return a.reshape(shape, order='F')

        result = {'fermi_eV': float(take((1,))[0]),
                  'lattice_A': take((3, 3)).T.copy(),
                  'reciprocal_inv_A': take((3, 3)).T.copy(),
                  'kpoints_fractional': take((3, nk)).T.copy(),
                  'weights': take((nk,)),
                  'energies_eV': take((nb, nk, ns))[:, :, 0].T.copy(),
                  'occupations': take((nb, nk, ns))[:, :, 0].T.copy()}
        for name, components in [('pseudo', 4), ('delta', 4),
                                 ('xi_restored', 3), ('correction_A', 3),
                                 ('velocity_eVA', 3)]:
            result[name] = take((nb, nb, nk, ns, components), '<c16')[:, :, :, 0, :].transpose(2, 3, 0, 1)
        require(f.read(32) == FOOTER and not f.read(1), 'missing completion footer or trailing bytes')
    require(abs(np.linalg.det(result['lattice_A'])) > 1e-12, 'singular source lattice')
    require(np.all(result['weights'] > 0)
            and abs(float(result['weights'].sum())-1) <= 1e-10, 'positive normalized producer weights required')
    reciprocal = 2*np.pi*np.linalg.inv(result['lattice_A']).T
    require(np.allclose(result['reciprocal_inv_A'], reciprocal, atol=1e-10, rtol=1e-10),
            'producer reciprocal lattice must include 2*pi')
    return result


def _hermiticity(a):
    return float(np.max(abs(a-a.conj().swapaxes(-1, -2))))


def _settings(case):
    out = (case/'OUTCAR').read_text(errors='replace')
    inp = incar_values((case/'INCAR').read_text())
    require(out.lstrip().startswith('vasp.5.4.4'), 'audited VASP 5.4.4 producer required')
    require('General timing and accounting' in out and 'EDIFF is reached' in out,
            'completed, electronically converged OUTCAR required')
    for key, expected in [('LOPTICS', True), ('LBERRY_EXPORT', True), ('LSPIN_EXPORT', True),
                          ('LPEAD', False), ('LNABLA', False), ('LREAL', False), ('LSORBIT', True)]:
        require(key in inp and boolean(inp[key]) == expected, 'explicit supported INCAR '+key+' required')
    for key in ('LREAL', 'LNABLA', 'LHFCALC', 'METAGGA', 'LEPSILON', 'LVEL'):
        require(not outcar_value(out, key, boolean), 'unsupported effective OUTCAR '+key)
    # This audited producer revision defaults these optional branches to false;
    # its OUTCAR does not print the false defaults. Reject explicit activation
    # in either the associated input or an effective-output declaration.
    for key in ('LDAU', 'LSPIRAL'):
        require(not boolean(inp.get(key, '.FALSE.')), 'unsupported INCAR '+key)
        values = re.findall(r'(?<!\w)'+key+r'\s*=\s*(\.?[TF](?:RUE|ALSE)?\.?)\b', out, re.I)
        require(all(not boolean(value) for value in values), 'unsupported effective OUTCAR '+key)
    try:
        qspiral = np.array([float(x) for x in inp.get('QSPIRAL', '0 0 0').split()])
    except ValueError as exc:
        raise ValueError('numeric zero QSPIRAL required') from exc
    require(qspiral.shape == (3,) and np.array_equal(qspiral, [0, 0, 0]),
            'finite-Q spin spirals are unsupported')
    require(outcar_value(out, 'ISYM', int) == -1 and outcar_value(out, 'ISPIN', int) == 1
            and outcar_value(out, 'NSW', int) == 0, 'static full-mesh spinor producer required')
    require(outcar_value(out, 'LSORBIT', boolean) and outcar_value(out, 'LNONCOLLINEAR', boolean),
            'effective OUTCAR SOC spinors required')
    try:
        saxis = np.array([float(x) for x in inp.get('SAXIS', '').split()])
    except ValueError as exc:
        raise ValueError('numeric SAXIS required') from exc
    require(saxis.shape == (3,) and np.array_equal(saxis, [0, 0, 1]), 'explicit SAXIS=0 0 1 required')
    angles = re.findall(r'Euler angles ALPHA=\s*([\d.Ee+\-]+)\s+BETA=\s*([\d.Ee+\-]+)', out)
    require(angles and all(abs(float(a))+abs(float(b)) < 1e-10 for a,b in angles),
            'effective OUTCAR spin basis is not Cartesian +z')


def audit_producer(run_dir):
    """Check completed same-run files and all exported operator identities."""
    case = Path(run_dir)
    names = ['SPIN_VELOCITY.bin', 'BERRY_CONNECTION.bin', 'WAVECAR', 'OUTCAR', 'INCAR', 'POSCAR', 'KPOINTS', 'run.json']
    digests = {n: sha256(case/n) for n in names}
    run = json.loads((case/'run.json').read_text())
    require(run.get('status') == 'FINISHED' and run.get('returncode') == 0
            and run.get('converged') is True and run.get('normal_footer') is True,
            'successful producer execution record required; exit zero alone is insufficient')
    require(not run.get('timed_out', False) and run.get('stop_reason') is None,
            'producer timeout or stop reason contradicts successful completion')
    for name in ('SPIN_VELOCITY.bin', 'BERRY_CONNECTION.bin', 'WAVECAR', 'OUTCAR'):
        require(run.get('output_sha256', {}).get(name) == digests[name], 'producer output association failed: '+name)
    for name in ('INCAR', 'POSCAR', 'KPOINTS'):
        require(run.get('input_sha256', {}).get(name) == digests[name], 'producer input association failed: '+name)
    require(isinstance(run.get('binary_sha256'), str) and re.fullmatch('[0-9a-f]{64}', run['binary_sha256']) is not None,
            'producer binary checksum required')
    _settings(case)
    data = read_spin_velocity(case/'SPIN_VELOCITY.bin')
    optical = read_connection(case/'BERRY_CONNECTION.bin')
    wave = Wavecar(case/'WAVECAR', spinor_components=2, spin=1)
    nk, nb = data['energies_eV'].shape
    require((wave.header.nkpoints, wave.header.nbands, wave.header.ispin) == (nk, nb, 1),
            'WAVECAR dimensions disagree')
    for left, right, label in [(data['energies_eV'], wave.energies, 'energies'),
                               (data['kpoints_fractional'], wave.kpoints, 'k points'),
                               (data['lattice_A'], wave.header.lattice, 'lattice')]:
        require(np.allclose(left, right, atol=1e-10, rtol=0), 'WAVECAR '+label+' disagree')
    require(optical.nspin == 1 and optical.energies_eV.shape == (1, nk, nb), 'optical dimensions disagree')
    require(np.allclose(optical.energies_eV[0], data['energies_eV'], atol=1e-10, rtol=0)
            and np.allclose(optical.kpoints_fractional, data['kpoints_fractional'], atol=1e-10, rtol=0)
            and np.allclose(optical.lattice_A, data['lattice_A'], atol=1e-10, rtol=0),
            'optical stream association failed')
    require(np.allclose(optical.weights, data['weights'], atol=1e-12, rtol=0),
            'optical and spin-stream k weights disagree')
    out = (case/'OUTCAR').read_text(errors='replace')
    require((outcar_value(out, 'NKPTS', int), outcar_value(out, 'NBANDS', int)) == (nk, nb),
            'effective OUTCAR dimensions disagree')
    pseudo_error = 0.
    for ik in range(nk):
        s, o = pauli_from_coefficients(wave.coefficients(ik, list(range(1, nb+1))))
        pseudo_error = max(pseudo_error, float(np.max(abs(data['pseudo'][ik]-np.concatenate((o[None], s))))))
    full = data['pseudo']+data['delta']
    velocity = data['velocity_eVA']
    gap = data['energies_eV'][:, None, :]-data['energies_eV'][:, :, None]
    old = optical.D_eVA[0]
    mask = (abs(gap[:, None]) > 1e-7) & np.isfinite(old)
    require(np.any(mask), 'no separated optical pair for consistency comparison')
    errors = {'metric_maxabs_error': float(np.max(abs(full[:, 0]-np.eye(nb)))),
              'pseudo_vs_saved_wavecar_maxabs': pseudo_error,
              'spin_hermiticity_maxabs': _hermiticity(full[:, 1:]),
              'augmentation_hermiticity_maxabs': _hermiticity(data['delta']),
              'velocity_hermiticity_eVA': _hermiticity(velocity),
              'optical_consistency_eVA': float(np.max(abs(velocity[mask]-old[mask]))),
              'diagonal_consistency_eVA': float(np.max(abs(np.diagonal(velocity, axis1=-2, axis2=-1)-optical.energy_der_eVA[0]))),
              'commutator_identity_eVA': float(np.max(abs(velocity-1j*(data['xi_restored']+gap[:, None]*data['correction_A']))))}
    require(errors['metric_maxabs_error'] < 1e-8 and pseudo_error < 1e-6, 'PAW metric or source gauge validation failed')
    require(max(errors[k] for k in ('spin_hermiticity_maxabs', 'augmentation_hermiticity_maxabs', 'velocity_hermiticity_eVA')) < 1e-8,
            'non-Hermitian physical operator; no repair is applied')
    require(max(errors[k] for k in ('optical_consistency_eVA', 'diagonal_consistency_eVA', 'commutator_identity_eVA')) < 1e-7,
            'producer velocity consistency failed')
    bound = min(float(np.linalg.eigvalsh(full[:, 0, None]+sign*full[:, 1:]).min()) for sign in (-1, 1))
    require(bound >= -5e-5, 'PAW spin exceeds the physical metric bound')
    require(digests == {n: sha256(case/n) for n in names}, 'producer files changed during validation')
    report = dict(status='PASS', nkpoints=nk, source_nbands=nb, source_sha256=digests,
                  producer_binary_sha256=run['binary_sha256'], pauli_metric_min_eigenvalue=bound,
                  optical_check_scope='Internal consistency with the same optical derivation; not an independent physical benchmark.', **errors)
    return data, report


def convert_run(run_dir, output_dir):
    """Create normalized physical and augmentation bundles after all checks."""
    out = Path(output_dir)
    require(not out.exists(), 'output directory exists')
    data, report = audit_producer(run_dir)
    nk, nb = data['energies_eV'].shape
    band_indices = np.arange(1, nb+1, dtype=np.int64)
    out.mkdir(parents=True, exist_ok=False)
    (out/'producer-audit.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    full = data['pseudo']+data['delta']
    np.savez_compressed(out/'physical-matrices.npz', energies_eV=data['energies_eV'],
                        kpoints_fractional=data['kpoints_fractional'], weights=data['weights'],
                        lattice_A=data['lattice_A'], spin_pauli=full[:, 1:], overlap=full[:, 0],
                        velocity_eVA=data['velocity_eVA'], band_indices=band_indices)
    metadata = {'schema': 'vaspberry.spin-velocity', 'version': 1, 'complete': True, 'producer_status': 'PASS',
                'data_npz_sha256': sha256(out/'physical-matrices.npz'),
                'source_wavecar_sha256': report['source_sha256']['WAVECAR'],
                'source_nbands': nb, 'nbands': nb, 'nkpoints': nk,
                'units': {'spin': 'dimensionless Pauli', 'overlap': 'dimensionless', 'velocity': 'eV angstrom = hbar*v', 'energies': 'eV', 'lattice': 'angstrom'},
                'spin_basis': 'Cartesian axes; SAXIS=(0,0,1)', 'spin_basis_id': 'cartesian', 'velocity_basis': 'cartesian',
                'matrix_axes': ['k', 'cartesian', 'bra', 'ket'], 'gauge': 'same_eigenvectors_as_source_wavecar',
                'operator_accuracy': 'paw_full_velocity', 'spin_operator_accuracy': 'paw_augmented', 'physical_paw_validated': True,
                'normalization_applied': False, 'spinor_components': 2, 'spin_multiplicity': 1,
                'full_source_band_coverage': True, 'diagonal_available': True, 'degenerate_blocks_available': True,
                'producer': 'Audited opt-in VASP5.4.4 spin and velocity instrumentation',
                'operator_scope': SCOPE, 'spin_current_construction': 'finite_band_projected_product',
                'method_limitations': LIMITATIONS, 'source_sha256': report['source_sha256'],
                'source_files': {n: n for n in report['source_sha256']}, 'producer_binary_sha256': report['producer_binary_sha256'],
                'producer_audit_sha256': sha256(out/'producer-audit.json'), 'metric_maxabs_error': report['metric_maxabs_error'], 'metric_tolerance': 5e-5,
                'source_energy_reference': 'Raw eigenvalue zero of matching WAVECAR; no shift applied.'}
    (out/'physical-matrices.json').write_text(json.dumps(metadata, indent=2, allow_nan=False)+'\n')
    np.savez_compressed(out/'augmentation.npz', delta_spin_pauli=data['delta'][:, 1:], delta_overlap=data['delta'][:, 0],
                        kpoints_fractional=data['kpoints_fractional'], energies_eV=data['energies_eV'],
                        lattice_A=data['lattice_A'], band_indices=band_indices)
    augmentation = {'schema': 'vaspberry.spin-augmentation', 'version': 1, 'complete': True,
                    'source_wavecar_sha256': metadata['source_wavecar_sha256'], 'spin_basis': metadata['spin_basis'],
                    'gauge': metadata['gauge'], 'units': {'spin': 'dimensionless Pauli', 'overlap': 'dimensionless'},
                    'producer': metadata['producer'], 'producer_status': 'PASS', 'source_sha256': report['source_sha256'],
                    'source_files': metadata['source_files'], 'data_npz_sha256': sha256(out/'augmentation.npz')}
    (out/'augmentation.json').write_text(json.dumps(augmentation, indent=2, allow_nan=False)+'\n')
    return metadata


def add_arguments(parser):
    parser.add_argument('--run-dir', required=True, type=Path, help='Completed instrumented VASP run with associated WAVECAR and run.json')
    parser.add_argument('--output-dir', required=True, type=Path, help='New normalized output directory')


def command(args):
    return convert_run(args.run_dir, args.output_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    print(json.dumps(command(parser.parse_args()), indent=2))
