#!/usr/bin/env python3
"""Prepare ordinary VASP inputs for the fixed planar graphene SOC example."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil

HERE=Path(__file__).resolve().parent

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024**2),b''):h.update(block)
    return h.hexdigest()

def restore_exact_charge_header(source, destination, poscar):
    """Undo six-decimal structure serialization; preserve the entire density payload."""
    original=source.read_bytes()
    lines=original.splitlines(keepends=True)
    exact=poscar.read_bytes()
    target=exact.splitlines()
    if len(lines)<10 or len(target)!=10:
        raise ValueError('Expected complete two-carbon structure headers')
    if lines[5].split()!=[b'C'] or lines[6].split()!=[b'2'] or not lines[7].strip().lower().startswith(b'd'):
        raise ValueError('Expected the documented two-carbon direct-coordinate CHGCAR header')
    for i in (1,2,3,4,8,9):
        old=[float(x) for x in lines[i].split()]
        new=[float(x) for x in target[i].split()]
        if not old or not all(math.isfinite(x) for x in old+new):
            raise ValueError('Nonfinite or missing geometry value')
        if len(old)!=len(new) or any(min(abs(x-y),abs(x-float(f'{y:.6f}')))>1e-12 for x,y in zip(old,new)):
            raise ValueError('CHGCAR geometry mismatch exceeds known writer rounding')
    body=b''.join(lines[10:])
    destination.write_bytes(exact+body)
    return dict(operation='Restore exact generating POSCAR header; all density and augmentation bytes unchanged',
                original_sha256=hashlib.sha256(original).hexdigest(),
                output_sha256=sha(destination),
                unchanged_density_payload_sha256=hashlib.sha256(body).hexdigest(),
                source_structure_sha256=sha(poscar))

def validate_scf_source(directory, expected_potential, expected_structure):
    """Bind the converged source and actual inputs to its preparation/run record."""
    names=['OUTCAR','CHGCAR','POTCAR','POSCAR','INCAR','KPOINTS']
    if any(not (directory/n).is_file() for n in names):
        raise ValueError('Incomplete SCF source')
    hashes={n:sha(directory/n) for n in names}
    if hashes['POTCAR']!=expected_potential or hashes['POSCAR']!=expected_structure:
        raise ValueError('SCF structure/potential differs from the documented graphene example')
    outcar=(directory/'OUTCAR').read_text(errors='replace')
    if 'EDIFF is reached' not in outcar or 'General timing and accounting' not in outcar:
        raise ValueError('SCF must have converged and ended normally')
    manifest=directory/'input_manifest.json'
    runfile=directory/'run.json'
    if manifest.is_file():
        record=json.loads(manifest.read_text())
        if record.get('stage')!='scf':raise ValueError('Source manifest is not an SCF preparation')
        recorded=record.get('inputs',{})
        record_file=manifest
    elif runfile.is_file():
        record=json.loads(runfile.read_text())
        if record.get('status')!='PASS':raise ValueError('Source run record is not PASS')
        recorded=record.get('input_sha256',{})
        record_file=runfile
        for n in ['OUTCAR','CHGCAR']:
            if record.get('outputs',{}).get(n,{}).get('sha256')!=hashes[n]:
                raise ValueError('Source output differs from run record: '+n)
    else:
        raise ValueError('Source requires input_manifest.json from the SCF preparation')
    for n in ['POSCAR','POTCAR','INCAR','KPOINTS']:
        if recorded.get(n)!=hashes[n]:raise ValueError('Source input differs from recorded preparation: '+n)
    incar=(directory/'INCAR').read_text()
    if not re.search(r'(?m)^ICHARG\s*=\s*2\s*(?:$|[!#])',incar):
        raise ValueError('Source must be the self-consistent ICHARG=2 stage')
    return dict(record_file=record_file.name,record_sha256=sha(record_file),files_sha256=hashes)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage',choices=['scf','mesh','gap','path'],required=True)
    p.add_argument('--potcar',type=Path,required=True,help='Matching locally licensed C POTCAR; see pseudopotential.json')
    p.add_argument('--scf-dir',type=Path,help='Completed matching SCF directory, required outside the scf stage')
    p.add_argument('--mesh',type=int,default=9,help='Gamma-centered size; multiple of three so both valleys are present')
    p.add_argument('--encut',type=float,default=520)
    p.add_argument('--ediff',type=float,help='Default 1e-11 for SCF, 1e-9 for fixed-charge runs')
    p.add_argument('--no-soc',action='store_true',help='Matched spinor control for gap diagnostics; intrinsic SOC otherwise')
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args()
    if a.ediff is None:a.ediff=1e-11 if a.stage=='scf' else 1e-9
    if a.output_dir.exists():p.error('Choose a new output directory')
    if a.mesh<3 or a.mesh%3:p.error('--mesh must be a positive multiple of three (at least 3)')
    if not math.isfinite(a.encut) or a.encut<400 or not 0<a.ediff<=1e-8:p.error('Require finite ENCUT >=400 eV and 0<EDIFF<=1e-8 eV')
    expected=json.loads((HERE/'pseudopotential.json').read_text())
    if not a.potcar.is_file() or sha(a.potcar)!=expected['sha256']:
        p.error('POTCAR must match the documented licensed carbon dataset')
    if a.stage!='scf':
        if a.scf_dir is None:p.error('--scf-dir is required for fixed-charge runs')
        try:source_record=validate_scf_source(a.scf_dir,expected['sha256'],sha(HERE/'inputs/POSCAR'))
        except (ValueError,KeyError,TypeError) as error:p.error(str(error))
    incar=(HERE/'inputs/scf/INCAR').read_text()
    changes={'ENCUT':f'{a.encut:.12g}','EDIFF':f'{a.ediff:.12g}'}
    if a.stage!='scf':changes.update(ICHARG='11',ISYM='-1',LCHARG='.FALSE.',NELMIN='20')
    if a.no_soc:changes['LSORBIT']='.FALSE.'
    for key,value in changes.items():
        incar,count=re.subn(r'(?m)^'+key+r'\s*=.*$',f'{key} = {value}',incar)
        if count!=1:raise ValueError('Expected one '+key+' in input template')
    if a.no_soc:incar+='LNONCOLLINEAR = .TRUE.\n'
    if a.stage=='scf':
        kpoints=f'Graphene SOC SCF, Gamma centered {a.mesh}x{a.mesh}\n0\nGamma\n{a.mesh} {a.mesh} 1\n0 0 0\n'
    else:
        if a.stage=='mesh':
            points=[(i/a.mesh,j/a.mesh,0.) for i in range(a.mesh) for j in range(a.mesh)]
        elif a.stage=='gap':points=[(1/3,1/3,0.),(-1/3,-1/3,0.),(0.,0.,0.)]
        else:
            vertices=[(0.,0.,0.),(1/3,1/3,0.),(.5,0.,0.),(0.,0.,0.)]
            points=[]
            for index,(first,last) in enumerate(zip(vertices[:-1],vertices[1:])):
                points += [tuple((1-t/16)*x+t/16*y for x,y in zip(first,last)) for t in range(16)]
            points.append(vertices[-1])
        kpoints=f'Graphene {a.stage}; full explicit reciprocal coordinates\n{len(points)}\nReciprocal\n'
        kpoints+=''.join(' '.join(f'{x:.16f}' for x in point)+' 1.0\n' for point in points)
    a.output_dir.mkdir(parents=True)
    for name in ['POSCAR']:
        shutil.copy2(HERE/'inputs'/name,a.output_dir/name)
    shutil.copy2(a.potcar,a.output_dir/'POTCAR')
    (a.output_dir/'INCAR').write_text(incar);(a.output_dir/'KPOINTS').write_text(kpoints)
    if a.stage!='scf':
        repair=restore_exact_charge_header(a.scf_dir/'CHGCAR',a.output_dir/'CHGCAR',a.scf_dir/'POSCAR')
        (a.output_dir/'charge-header-repair.json').write_text(json.dumps(repair,indent=2)+'\n')
    record=dict(stage=a.stage,mesh=a.mesh,encut_eV=a.encut,ediff_eV=a.ediff,intrinsic_soc=not a.no_soc,
                fixed_geometry=True,occupied_spinor_bands=8,source_bands=16,
                inputs={n:sha(a.output_dir/n) for n in ['POSCAR','INCAR','KPOINTS','POTCAR','CHGCAR'] if (a.output_dir/n).is_file()},
                preparation_sha256=sha(Path(__file__)))
    if a.stage!='scf':record['source_scf']=source_record
    (a.output_dir/'input_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    print('Prepared '+str(a.output_dir)+'. Run ordinary noncollinear VASP in this directory.')

if __name__=='__main__':main()
