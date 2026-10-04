#!/usr/bin/env python3
"""Create declared synthetic optical files for selected-WAVEDER software examples.

This is a format/contraction fixture, not VASP output or a material model.
No electronic-structure or convergence claim follows from its response.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import runpy
import struct

import numpy as np

ROOT = Path(__file__).resolve().parents[3]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_waveder(path, c):
    nb, nd, nk, ns, _ = c.shape
    def record(payload):
        return struct.pack('<i', len(payload)) + payload + struct.pack('<i', len(payload))
    path.write_bytes(record(struct.pack('<4i', nb, nd, nk, ns)) +
                     record(struct.pack('<d', 0.)) +
                     record(np.zeros((3, 3), '<f8').tobytes()) +
                     record(np.asarray(c, dtype='<c8').tobytes(order='F')))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True,
                        help='new fixture directory, with full and rectangular optical inputs')
    args = parser.parse_args()
    out = args.output_dir.resolve()
    if out.exists():
        parser.error('output directory exists; use a new directory')
    out.mkdir(parents=True)
    helper = ROOT / 'tests/test_fortran_waveder_default.py'
    write_sources = runpy.run_path(str(helper))['write_sources']
    q = ((0., 0., 0.), (0., .5, 0.), (.5, 0., 0.), (.5, .5, 0.))
    c = np.zeros((4, 4, 4, 1, 3), dtype=np.complex64)
    base = {(2, 0): (.2+.1j, .3-.2j, .1j),
            (2, 1): (.1+.2j, .25+.3j, -.1j),
            (3, 0): (.15-.1j, -.25+.2j, .05j),
            (3, 1): (.1-.1j, .2+.25j, -.05j),
            (3, 2): (.3+.25j, -.15+.2j, .2j)}
    for k in range(4):
        for (m, n), vector in base.items():
            c[m, n, k, 0] = np.asarray(vector) * (1 + k / 8)
            c[n, m, k, 0] = c[m, n, k, 0].conj()
    for name, nd in [('full', 4), ('rectangular', 2)]:
        path = out / name
        write_sources(path, q=q, nd=nd)
        write_waveder(path / 'WAVEDER', c[:, :nd])
        with (path / 'OUTCAR').open('a') as stream:
            stream.write('SYNTHETIC SOFTWARE FIXTURE; NOT AN EXECUTED VASP CALCULATION\n')
    config = Path(__file__).with_name('selected.ini')
    (out/'analysis.ini').write_text(config.read_text())
    # Independent direct target-band sums over every m in the full matrix.
    omega = -2 * np.imag((c[..., 0].astype(complex).conj() *
                          c[..., 1].astype(complex)).sum(axis=0))[:, :, 0].T
    expected = {'kind': 'synthetic_format_and_contraction_fixture',
                'not_material_data': True, 'mesh': [2, 2], 'source_nbands': 4,
                'energies_eV': [-2., -2., 1., 3.],
                'source_occupations': [1., 1., 0., 0.],
                'source_spinor_components': 2, 'source_ispin': 1,
                'full_matrix_shape_m_n_k_spin_direction': list(c.shape),
                'omega_by_k_band_A2': omega.tolist(),
                'selected_3_4_trace_A2': omega[:, 2:4].sum(axis=1).tolist(),
                'occupied_1_2_trace_A2': omega[:, :2].sum(axis=1).tolist(),
                'source_helper_sha256': sha(helper),
                'generator_sha256': sha(Path(__file__)),
                'analysis_ini_sha256': sha(out/'analysis.ini'),
                'inputs_sha256': {str(p.relative_to(out)): sha(p)
                                  for name in ('full', 'rectangular')
                                  for p in sorted((out/name).iterdir())}}
    (out/'expected.json').write_text(json.dumps(expected, indent=2)+'\n')
    print(f'Created synthetic software fixture in {out}; no VASP/material result is claimed.')


if __name__ == '__main__':
    main()
