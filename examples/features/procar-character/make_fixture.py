#!/usr/bin/env python3
"""Create tiny analytic format fixtures, not a VASP material calculation.

The two-level curvature and projector characters are assigned independently
for testing the join/integration conventions. They are not derived from the
artificial plane-wave coefficient records stored in this toy WAVECAR.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import struct

import numpy as np


def create(directory):
    out = Path(directory)
    if out.exists():
        raise ValueError('choose a fresh fixture directory')
    out.mkdir(parents=True)
    q = np.array([[0., 0., 0.], [0., .5, 0.], [.5, 0., 0.], [.5, .5, 0.]])
    stride, nk, nb = 512, 4, 2
    blob = bytearray(stride*(2+nk*(nb+1)))
    struct.pack_into('<3d', blob, 0, stride, 1, 45200)
    struct.pack_into('<12d', blob, stride, nk, nb, 120., *np.eye(3).ravel())
    for ik, point in enumerate(q):
        count = sum(sum((2*math.pi*(point[j]+g[j]))**2 for j in range(3))/.262465831 < 120.
                    for g in [(x, y, z) for x in range(-2, 3) for y in range(-2, 3) for z in range(-2, 3)])
        npw = 2*count
        record = 2+ik*(nb+1)
        struct.pack_into('<10d', blob, record*stride, npw, *point, -1., 0., 1., 1., 0., 0.)
        for b in range(nb):
            coeff = np.exp(2j*np.pi*b*np.arange(npw)/npw)/math.sqrt(npw)
            payload = np.asarray(coeff, dtype='<c8').tobytes()
            blob[(record+1+b)*stride:(record+1+b)*stride+len(payload)] = payload
    (out/'WAVECAR').write_bytes(blob)
    # Spin frame z points along Cartesian +x. Lower atom is +x, upper -x.
    rotation = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])
    raw = np.zeros((2, 2, 4)); raw[0, 0] = [.2, 0, 0, .2]; raw[1, 1] = [.6, 0, 0, -.6]
    lines = ['PROCAR lm decomposed', '# of k-points: 4 # of bands: 2 # of ions: 2']
    for k, point in enumerate(q):
        lines.append(f'k-point {k+1} : '+ ' '.join(f'{v:.8f}' for v in point)+' weight = 0.25000000')
        for b in range(nb):
            lines += [f'band {b+1} # energy {(-1+2*b):.8f} # occ. {1-b:.8f}', 'ion s p tot']
            for component in range(4):
                block = raw[..., component]
                for ion in range(2):
                    lines.append(str(ion+1)+' '+' '.join(f'{v:.3f}' for v in [*block[ion], block[ion].sum()]))
                lines.append('tot '+' '.join(f'{v:.3f}' for v in [*block.sum(axis=0), block.sum()]))
    (out/'PROCAR').write_text('\n'.join(lines)+'\n')
    lines = ['vasp.5.4.4 SYNTHETIC FORMAT FIXTURE; not a production run',
             ' LNONCOLLINEAR = T', ' ISPIN = 1', ' NSW = 0', ' LORBIT = 11',
             ' k-points NKPTS = 4 NBANDS = 2',
             ' direct lattice vectors                 reciprocal lattice vectors']
    for row in np.eye(3): lines.append(' '.join(f'{v:.10f}' for v in [*row, *row]))
    lines += [' k-points in reciprocal lattice and weights:']
    for point in q: lines.append(' '.join(f'{v:.8f}' for v in [*point, .25]))
    lines += [' transformation matrix from SAXIS to cartesian coordinates', ' ---------------------------------------------------------']
    for row in rotation:
        lines.append(' '.join(f'{v:.7f} m_{a}' for v, a in zip(row, 'xyz')))
    for k, point in enumerate(q):
        lines += [f' k-point {k+1} : '+' '.join(f'{v:.4f}' for v in point),
                  ' band No.  band energies     occupation', ' 1 -1.0000 1.0000', ' 2 1.0000 0.0000']
    lines += [' General timing and accounting informations for this job:']
    (out/'OUTCAR').write_text('\n'.join(lines)+'\n')
    groups = {'groups': [{'name': 'lower', 'ions': [1]}, {'name': 'upper', 'ions': [2]},
                         {'name': 'all', 'ions': [1, 2]}, {'name': 'lower_s', 'ions': [1], 'orbitals': ['s']}]}
    (out/'groups.json').write_text(json.dumps(groups, indent=2)+'\n')
    # This deliberately uses the current native pair schema for importer tests;
    # the fixture label is carried into the user-specified energy reference.
    meta = dict(schema='VASPBERRY_BARE_MOMENTUM_KUBO_PAIRS_V1', result_kind='UNORDERED_INTERBAND_NUMERATORS',
        normalization='STANDARD_MINUS_TWO_IM', operator='WAVECAR_BARE_MOMENTUM_NO_PAW_NONLOCAL_VELOCITY',
        berry_connection='A_i=i<u|d/dk_i u>', occupation_weighting='NONE', denominator_weighting='NONE',
        pair_order='n_lt_m', gap_definition='ABS_EN_MINUS_EM', numerator_units='eV^2*Angstrom^2',
        components='yz,zx,xy', reciprocal_convention='2pi', index_base=1, source_nbands=nb,
        source_nkpoints=nk, source_nspin=1, spinor_components=2, pairs_per_k=1, expected_rows=nk,
        fixture='ANALYTIC_ASSIGNED_NUMERATOR_NOT_A_MATERIAL_RESULT')
    for i in range(3):
        meta[f'lattice_A_{i+1}'] = ','.join(map(str, np.eye(3)[i]))
        meta[f'reciprocal_inv_A_{i+1}'] = ','.join(map(str, (2*np.pi*np.eye(3))[i]))
    with (out/'PAIRS.csv').open('w', newline='') as f:
        for key, value in meta.items(): f.write(f'# {key}={value}\n')
        w = csv.writer(f)
        w.writerow(['spin', 'k_index', 'n_band', 'm_band', 'kx_frac', 'ky_frac', 'kz_frac',
                    'energy_n_eV', 'energy_m_eV', 'numerator_yz_eV2_A2', 'numerator_zx_eV2_A2', 'numerator_xy_eV2_A2', 'gap_eV'])
        for k, point in enumerate(q): w.writerow([1, k+1, 1, 2, *point, -1., 1., 0., 0., 4., 2.])
        f.write('# result_status=PASS\n')
    (out/'expected.json').write_text(json.dumps(dict(fixture=True, band1_omega_z_A2=1.,
        band2_omega_z_A2=-1., at_T0_mu0_band1={'lower_charge': -.4*np.pi, 'lower_plus': -.4*np.pi,
        'lower_minus': 0., 'upper_charge': -1.2*np.pi, 'upper_plus': 0., 'upper_minus': -1.2*np.pi,
        'all_charge': -1.6*np.pi}), indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    create(parser.parse_args().output_dir)
