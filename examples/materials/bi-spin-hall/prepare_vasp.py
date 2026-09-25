#!/usr/bin/env python3
"""Prepare the actual two-atom Bi SCF, ordinary WAVECAR or PAW spin calculation."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import shutil

import numpy as np

HERE = Path(__file__).resolve().parent


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024**2),b''):h.update(block)
    return h.hexdigest()


def restore_charge(target, metadata):
    parts = []
    for item in metadata['charge_archive_parts']:
        p = HERE/'inputs/scf'/item['file']
        if p.stat().st_size != item['bytes'] or sha(p) != item['sha256']:
            raise ValueError('Supplied charge archive part failed integrity check: '+p.name)
        parts.append(p.read_bytes())
    payload = b''.join(parts)
    if hashlib.sha256(payload).hexdigest() != metadata['compressed_charge_sha256']:
        raise ValueError('Combined charge archive failed integrity check')
    with gzip.GzipFile(fileobj=io.BytesIO(payload)) as source,target.open('wb') as out:
        shutil.copyfileobj(source,out)
    if sha(target) != metadata['charge_sha256']:
        raise ValueError('Restored SCF charge density failed integrity check')


def structure_header(path):
    """Read only the small VASP5 structure header, including CHGCAR headers."""
    with Path(path).open() as f:
        lines=[next(f) for _ in range(8)]
        scale=float(lines[1].split()[0])
        lattice=np.array([[float(x) for x in line.split()[:3]] for line in lines[2:5]])
        if not np.isfinite(lattice).all() or abs(np.linalg.det(lattice)) < 1e-12 or scale == 0:
            raise ValueError('Invalid structure lattice')
        factor=scale if scale > 0 else (-scale/abs(np.linalg.det(lattice)))**(1/3)
        lattice*=factor
        species=lines[5].split();counts=[int(x) for x in lines[6].split()]
        if len(species) != len(counts) or sum(counts) != 2:
            raise ValueError('This example requires the matching two-atom Bi structure')
        mode=lines[7].strip().lower()
        if mode.startswith('s'):mode=next(f).strip().lower()
        if not mode or mode[0] not in 'dck':raise ValueError('Invalid structure coordinate mode')
        positions=np.array([[float(x) for x in next(f).split()[:3]] for _ in range(sum(counts))])
        if mode[0] != 'd':positions=(factor*positions)@np.linalg.inv(lattice)
        if not np.isfinite(positions).all():raise ValueError('Nonfinite structure coordinates')
    return species,counts,lattice,positions


def check_charge_geometry(path):
    a=structure_header(HERE/'inputs/POSCAR');b=structure_header(path)
    delta=a[3]-b[3];delta-=np.rint(delta)
    # This supported VASP CHGCAR writer prints the structure to six decimals.
    if a[:2] != b[:2] or not np.allclose(a[2],b[2],rtol=0,atol=1e-6) or np.max(abs(delta))>1e-6:
        raise ValueError('CHGCAR structure differs from the supplied Bi geometry')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage',choices=['scf','wavecar','spin'],default='spin',
                   help='wavecar uses ordinary VASP for native Z2/Chern; spin adds optional physical-matrix exports')
    p.add_argument('--potcar',type=Path,required=True,help='Licensed PAW-PBE Bi 08Apr2002 potential; see PSEUDOPOTENTIAL.md')
    p.add_argument('--charge',type=Path,help='Own completed matching SCF CHGCAR; default restore the supplied fresh SCF density')
    p.add_argument('--mesh',nargs=2,type=int,default=[12,12],metavar=('NX','NY'))
    p.add_argument('--chunks',type=int,default=1,help='Number of contiguous fixed-charge k-point partitions; default 1')
    p.add_argument('--chunk-index',type=int,default=0,help='Zero-based partition index; default 0')
    p.add_argument('--nbands',type=int,help='Default 32 for SCF, 48 for WAVECAR/spin export; occupied count is 10')
    p.add_argument('--output-dir',type=Path,required=True)
    a = p.parse_args()
    nb = a.nbands if a.nbands is not None else 32 if a.stage == 'scf' else 48
    if min(a.mesh) < 2 or nb <= 10:p.error('mesh sizes >=2 and NBANDS>10 required')
    nx,ny = a.mesh
    if not 1 <= a.chunks <= nx*ny:p.error('--chunks must be between 1 and NX*NY')
    if not 0 <= a.chunk_index < a.chunks:p.error('--chunk-index must satisfy 0 <= index < chunks')
    if a.stage == 'scf' and a.chunks != 1:p.error('SCF density generation cannot be partitioned; use --chunks 1')
    if a.stage == 'scf' and a.charge is not None:p.error('--charge is used only for fixed-charge wavecar/spin stages')
    indices = np.array_split(np.arange(nx*ny),a.chunks)[a.chunk_index]
    meta = json.loads((HERE/'inputs/provenance.json').read_text())
    if not a.potcar.is_file() or sha(a.potcar) != meta['potcar_sha256']:
        p.error('POTCAR must match the dataset in PSEUDOPOTENTIAL.md')
    if a.output_dir.exists():p.error('output directory exists; choose a new directory')
    if a.charge is not None and not a.charge.is_file():p.error('CHGCAR does not exist')
    if a.charge is not None:
        try:check_charge_geometry(a.charge)
        except (ValueError,StopIteration,IndexError) as exc:p.error(str(exc))
    template = HERE/'inputs'/('scf/INCAR' if a.stage == 'scf' else 'INCAR.nscf')
    incar,count = re.subn(r'(?m)^NBANDS\s*=.*$',f'NBANDS = {nb}',template.read_text())
    if count != 1:raise ValueError('Reference INCAR must contain one NBANDS setting')
    if a.stage == 'wavecar':
        # The electronic Hamiltonian and fixed charge are unchanged. Ordinary
        # VASP only needs WAVECAR; physical operators and optics are optional.
        incar = re.sub(r'(?m)^(?:LOPTICS|LPEAD|LNABLA|LBERRY_EXPORT|LSPIN_EXPORT)\s*=.*\n?', '', incar)
        incar = re.sub(r'(?m)^SYSTEM\s*=.*$', 'SYSTEM = Bi full-mesh SOC wavefunctions', incar)
    a.output_dir.mkdir(parents=True)
    (a.output_dir/'INCAR').write_text(incar)
    shutil.copyfile(HERE/'inputs/POSCAR',a.output_dir/'POSCAR')
    shutil.copyfile(a.potcar,a.output_dir/'POTCAR')
    if a.stage == 'scf':
        kpoints=f'Bi Gamma-centered SCF mesh\n0\nGamma\n{nx} {ny} 1\n0 0 0\n'
    else:
        # The ordered union of all partitions is the original complete mesh.
        kind='spin' if a.stage == 'spin' else 'wavefunction'
        title=f'Bi full Gamma-centered SOC {kind} mesh' if a.chunks == 1 else f'Bi Gamma-centered SOC {kind} mesh: chunk {a.chunk_index+1}/{a.chunks}'
        kpoints=f'{title}\n{len(indices)}\nReciprocal\n'
        kpoints+=''.join(f'{(n//ny)/nx:.16f} {(n%ny)/ny:.16f} 0.0 1.0\n' for n in indices)
    (a.output_dir/'KPOINTS').write_text(kpoints)
    if a.stage != 'scf':
        if a.charge is None:restore_charge(a.output_dir/'CHGCAR',meta)
        else:shutil.copyfile(a.charge,a.output_dir/'CHGCAR')
    names=['INCAR','KPOINTS','POSCAR','POTCAR']+(['CHGCAR'] if a.stage != 'scf' else [])
    result=dict(status='PREPARED',stage=a.stage,mesh=[nx,ny,1],nbands=nb,occupied=10,
        fixed_geometry=True,spinor_components=2,
        chunk_count=a.chunks,chunk_index=a.chunk_index,npoints=len(indices),full_mesh_point_count=nx*ny,
        point_count_scope='explicit KPOINTS rows' if a.stage != 'scf' else 'requested full mesh before VASP symmetry reduction',
        mesh_partition=dict(ordering='qx outer, qy inner',index_base=0,
            index_start=int(indices[0]),index_stop_exclusive=int(indices[-1])+1,
            rule='numpy.array_split of the ordered full mesh; earlier partitions receive any remainder',
            union='all chunk indices 0 through chunk_count-1 cover the full mesh once',
            weights=('equal KPOINTS weights; spin-merge restores 1/(NX*NY) on the complete mesh'
                     if a.stage == 'spin' else 'equal KPOINTS weights; normalize on the complete mesh'
                     if a.stage == 'wavecar' else 'VASP determines symmetry weights')),
        charge_source=('atomic initialization' if a.stage == 'scf' else 'supplied fresh SCF' if a.charge is None else 'user supplied; source compatibility must be checked'),
        input_sha256={n:sha(a.output_dir/n) for n in names},preparation_sha256=sha(Path(__file__)))
    (a.output_dir/'input_manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Prepared '+str(a.output_dir)+'. Run the matching licensed noncollinear VASP executable in this directory.')


if __name__ == '__main__':main()
