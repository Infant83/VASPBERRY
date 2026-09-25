#!/usr/bin/env python3
"""Prepare fixed-geometry MoS2 stacking examples for ordinary VASP calculations."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil

import numpy as np

HERE = Path(__file__).resolve().parent


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def structure(path):
    """Read VASP5 structure headers, including the rounded CHGCAR header."""
    with Path(path).open() as stream:
        lines = [next(stream) for _ in range(8)]
        scale = float(lines[1].split()[0])
        lattice = np.array([[float(v) for v in line.split()[:3]] for line in lines[2:5]])
        if scale <= 0 or not np.isfinite(lattice).all() or abs(np.linalg.det(lattice)) < 1e-8:
            raise ValueError('This example requires a positive-scale finite VASP lattice')
        lattice *= scale
        species = lines[5].split(); counts = [int(v) for v in lines[6].split()]
        mode = lines[7].strip().lower()
        if mode.startswith('s'):
            mode = next(stream).strip().lower()
        if not mode or mode[0] not in 'dck':
            raise ValueError('Unrecognized structure coordinate mode')
        positions = np.array([[float(v) for v in next(stream).split()[:3]] for _ in range(sum(counts))])
        if mode[0] != 'd':
            positions = scale*positions @ np.linalg.inv(lattice)
    return species, counts, lattice, positions


def check_source(source, template, potential_hash):
    text = (source/'OUTCAR').read_text(errors='replace')
    if not ('General timing and accounting informations for this job:' in text
            and 'aborting loop because EDIFF is reached' in text):
        raise ValueError('A completed electronically converged SCF OUTCAR is required')
    if sha(source/'POTCAR') != potential_hash or sha(source/'POSCAR') != sha(template/'POSCAR'):
        raise ValueError('SCF structure/potential differs from this reference case')
    reference = structure(template/'POSCAR'); charge = structure(source/'CHGCAR')
    delta = charge[3]-reference[3]; delta -= np.rint(delta)
    if (reference[:2] != charge[:2] or not np.allclose(reference[2], charge[2], rtol=0, atol=1e-6)
            or np.max(abs(delta)) > 1e-6):
        raise ValueError('SCF charge density has a different structure')
    # Reference reproduction permits a changed source band count or k mesh,
    # while keeping the Hamiltonian and SCF settings of the material fixed.
    def settings(path):
        values = {}
        for line in path.read_text().splitlines():
            line = re.split('[!#]', line, maxsplit=1)[0]
            if '=' in line:
                key, value = line.split('=', 1); values[key.strip().upper()] = value.strip().upper()
        return {k:v for k,v in values.items() if k not in {'SYSTEM', 'NBANDS'}}
    if settings(source/'INCAR') != settings(template/'INCAR'):
        raise ValueError('SCF settings differ from the reference; use matching fresh-SCF inputs')


def main():
    meta = json.loads((HERE/'inputs/provenance.json').read_text())
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=tuple(meta['cases']), required=True)
    parser.add_argument('--stage', choices=['scf', 'scf-dipole', 'path', 'mesh'], required=True)
    parser.add_argument('--potcar', type=Path, required=True,
                        help='Licensed concatenated Mo/S PAW datasets listed in PSEUDOPOTENTIAL.md')
    parser.add_argument('--scf-dir', type=Path, help='Completed matching fresh SCF run; required for path/mesh')
    parser.add_argument('--mesh', type=int, nargs=2, default=[6, 6], metavar=('NX', 'NY'))
    parser.add_argument('--nbands', type=int, help='Default from each supplied stage template')
    parser.add_argument('--path-points-per-segment', type=int, default=13)
    parser.add_argument('--without-optics', action='store_true',
                        help='Prepare WAVECAR only; omit standard LOPTICS output in path/mesh stage')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(); case = meta['cases'][args.case]
    template = HERE/'inputs'/args.case/args.stage
    if not template.is_dir(): parser.error('This case has no reference template for the requested stage')
    template_text = (template/'INCAR').read_text()
    reference_nbands = int(re.search(r'(?m)^NBANDS\s*=\s*(\d+)', template_text).group(1))
    nbands = args.nbands if args.nbands is not None else reference_nbands
    if args.output_dir.exists(): parser.error('Choose a new output directory')
    if min(args.mesh) < 2 or args.path_points_per_segment < 3:
        parser.error('Mesh sizes >=2 and path points per segment >=3 required')
    if nbands <= case['occupied_spinor_bands_expected']+4:
        parser.error('Retain occupied states, empty optical bands and a boundary guard band')
    if not args.potcar.is_file() or sha(args.potcar) != meta['potcar_sha256']:
        parser.error('POTCAR does not match the documented reference Mo/S datasets')
    if (args.stage == 'scf') == (args.scf_dir is not None):
        parser.error('Use --scf-dir for every continuation stage; omit it for the initial scf stage')
    if args.scf_dir is not None:
        parent_stage = 'scf' if args.stage == 'scf-dipole' else case.get('response_scf_stage', 'scf')
        check_source(args.scf_dir, HERE/'inputs'/args.case/parent_stage, meta['potcar_sha256'])
    incar, count = re.subn(r'(?m)^NBANDS\s*=.*$', f'NBANDS = {nbands}', template_text)
    if count != 1: raise ValueError('Expected exactly one NBANDS in reference template')
    if args.without_optics:
        incar = re.sub(r'(?m)^(LOPTICS|LPEAD|LNABLA|CSHIFT|NEDOS)\s*=.*\n', '', incar)
    nx, ny = args.mesh
    if args.stage in {'scf', 'scf-dipole'}:
        kpoints = f'Gamma-centered SCF mesh\n0\nGamma\n{nx} {ny} 1\n0 0 0\n'
        sampling = {'kind':'scf_mesh', 'mesh':[nx, ny, 1]}
    else:
        if args.stage == 'mesh':
            points = np.array([[i/nx, j/ny, 0.] for j in range(ny) for i in range(nx)])
            sampling = {'kind':'uniform_full_2d', 'mesh':[nx, ny, 1]}
        else:
            nodes = np.array(meta['path']['vertices_fractional']); n = args.path_points_per_segment
            points = np.vstack([np.linspace(a, b, n)[:-1] for a,b in zip(nodes[:-1],nodes[1:])]+[nodes[-1:]])
            sampling = {'kind':'path', 'labels':meta['path']['labels'],
                        'node_indices_zero_based':[i*(n-1) for i in range(len(nodes))]}
        kpoints = f'MoS2 {args.stage}\n{len(points)}\nReciprocal\n'
        kpoints += ''.join(' '.join(f'{v:.16f}' for v in q)+' 1\n' for q in points)
    args.output_dir.mkdir(parents=True)
    (args.output_dir/'INCAR').write_text(incar); (args.output_dir/'KPOINTS').write_text(kpoints)
    shutil.copyfile(template/'POSCAR', args.output_dir/'POSCAR')
    shutil.copyfile(args.potcar, args.output_dir/'POTCAR')
    if args.scf_dir is not None: shutil.copyfile(args.scf_dir/'CHGCAR', args.output_dir/'CHGCAR')
    if args.stage == 'scf-dipole':
        shutil.copyfile(args.scf_dir/'WAVECAR', args.output_dir/'WAVECAR')
    names = ['INCAR','KPOINTS','POSCAR','POTCAR']+(['CHGCAR'] if args.scf_dir is not None else [])
    if args.stage == 'scf-dipole': names.append('WAVECAR')
    record = {'status':'PREPARED', 'case':args.case, 'stage':args.stage,
              'sampling':sampling, 'nbands':nbands, 'fixed_geometry':True,
              'operator_route':'ordinary VASP; no source instrumentation',
              'charge_source':'fresh matching SCF' if args.scf_dir else 'atomic initialization',
              'input_sha256':{name:sha(args.output_dir/name) for name in names},
              'preparation_script_sha256':sha(Path(__file__))}
    (args.output_dir/'preparation.json').write_text(json.dumps(record, indent=2)+'\n')
    print(f'Prepared {args.case} {args.stage}: {args.output_dir}. Run your licensed noncollinear VASP there.')


if __name__ == '__main__': main()
