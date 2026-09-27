#!/usr/bin/env python3
"""Prepare actual VASP points around both graphene valleys at fixed SCF density."""
from pathlib import Path
import argparse
import csv
import importlib.util
import json
import math
import re
import shutil
import numpy as np

HERE = Path(__file__).resolve().parent
PARENT = HERE.parents[1]
spec = importlib.util.spec_from_file_location('graphene_prepare', PARENT/'prepare_vasp.py')
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scf-dir', type=Path, required=True)
    parser.add_argument('--potcar', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--ediff', type=float, default=1e-9)
    parser.add_argument('--nelmin', type=int, default=20)
    parser.add_argument('--angles', type=int, default=8, help='Even counterclockwise angles per circle; at least 8')
    parser.add_argument('--radii', default='1e-9,3e-9,1e-8,3e-8,1e-7,3e-7,1e-6,3e-6,1e-5', help='Comma-separated Cartesian radii in inverse Angstrom')
    args = parser.parse_args()
    radii = [float(x) for x in args.radii.split(',')]
    if args.output_dir.exists(): parser.error('Choose a fresh output directory')
    if args.angles < 8 or args.angles % 4: parser.error('--angles must be a multiple of 4, at least 8')
    if not 0 < args.ediff <= 1e-8 or not 20 <= args.nelmin <= 200: parser.error('Require 0<EDIFF<=1e-8 and 20<=NELMIN<=200')
    if not radii or radii != sorted(set(radii)) or any(not math.isfinite(x) or x < 1e-10 or x > 1e-2 for x in radii): parser.error('Radii must increase strictly within [1e-10,1e-2] inverse Angstrom')
    potcar = json.loads((PARENT/'pseudopotential.json').read_text())
    if base.sha(args.potcar) != potcar['sha256']: parser.error('POTCAR differs from the parent example')
    poscar = PARENT/'inputs/POSCAR'
    source = base.validate_scf_source(args.scf_dir, potcar['sha256'], base.sha(poscar))
    lattice = np.array([[float(x) for x in row.split()] for row in poscar.read_text().splitlines()[2:5]])
    reciprocal = 2*np.pi*np.linalg.inv(lattice).T
    inverse = np.linalg.inv(reciprocal)
    points = []
    for valley, center in [('K', np.array([1/3,1/3,0.])), ('Kprime', np.array([-1/3,-1/3,0.]))]:
        points.append(dict(valley=valley, radius_Ainv=0., angle_index=-1, theta_rad=0., qx_Ainv=0., qy_Ainv=0., k=center))
        for radius in radii:
            for j in range(args.angles):
                theta = 2*np.pi*j/args.angles
                q = radius*np.array([math.cos(theta), math.sin(theta), 0.])
                points.append(dict(valley=valley, radius_Ainv=radius, angle_index=j, theta_rad=theta, qx_Ainv=q[0], qy_Ainv=q[1], k=center+q@inverse))
    args.output_dir.mkdir(parents=True)
    shutil.copy2(poscar, args.output_dir/'POSCAR')
    shutil.copy2(args.potcar, args.output_dir/'POTCAR')
    repair = base.restore_exact_charge_header(args.scf_dir/'CHGCAR', args.output_dir/'CHGCAR', poscar)
    (args.output_dir/'charge-header-repair.json').write_text(json.dumps(repair, indent=2)+'\n')
    incar = (PARENT/'inputs/scf/INCAR').read_text()
    for key, value in dict(ENCUT='520', EDIFF=f'{args.ediff:.12g}', NELMIN=str(args.nelmin), ICHARG='11', ISYM='-1', LCHARG='.FALSE.').items():
        incar, count = re.subn(r'(?m)^'+key+r'\s*=.*$', f'{key} = {value}', incar)
        if count != 1: raise ValueError('Expected exactly one '+key)
    (args.output_dir/'INCAR').write_text(incar)
    text = f'Graphene intrinsic SOC; Cartesian radial convergence points\n{len(points)}\nReciprocal\n'
    text += ''.join(' '.join(f'{v:.16f}' for v in row['k'])+' 1.0\n' for row in points)
    (args.output_dir/'KPOINTS').write_text(text)
    with (args.output_dir/'points.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['k_index','valley','radius_Ainv','angle_index','theta_rad','qx_Ainv','qy_Ainv','k1','k2','k3'])
        for i, row in enumerate(points,1): writer.writerow([i]+[row[k] for k in ['valley','radius_Ainv','angle_index','theta_rad','qx_Ainv','qy_Ainv']]+row['k'].tolist())
    record = dict(stage='local_valley_convergence', ediff_eV=args.ediff, nelmin=args.nelmin, angles=args.angles, radii_Ainv=radii,
                  points=len(points), encut_eV=520, source_bands=16, occupied_bands=[1,8], reciprocal_Ainv=reciprocal.tolist(),
                  source_scf=source, preparation_sha256=base.sha(Path(__file__)),
                  scope='Actual ordinary VASP fixed-density local points; not a full BZ integration mesh',
                  inputs={n:base.sha(args.output_dir/n) for n in ['POSCAR','INCAR','KPOINTS','CHGCAR','POTCAR','points.csv']})
    (args.output_dir/'input_manifest.json').write_text(json.dumps(record, indent=2)+'\n')
    print(f'Prepared {len(points)} actual VASP points in {args.output_dir}')


if __name__ == '__main__': main()
