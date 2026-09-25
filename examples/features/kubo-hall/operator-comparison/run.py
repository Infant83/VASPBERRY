#!/usr/bin/env python3
"""Compare native and full PAW velocity on the same completed MoS2 eigenstates."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT/'tools'))
from exported_matrix_kubo import sha256
from wavecar_fukui import Wavecar


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-dir',type=Path,required=True,help='completed opt-in producer run with its run.json')
    p.add_argument('--binary',type=Path,required=True,help='native serial VASPBERRY executable')
    p.add_argument('--pair-band-max',type=int,nargs='+',default=[40,50])
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args()
    if a.output_dir.exists():p.error('output directory exists')
    source=a.run_dir.resolve();out=a.output_dir.resolve();binary=a.binary.resolve()
    wave=Wavecar(source/'WAVECAR',spinor_components=2)
    nk,nb=wave.energies.shape;mesh=int(np.sqrt(nk))
    if mesh*mesh!=nk or nb<=18:p.error('MoS2 comparison requires a full square mesh and empty bands')
    if len(set(a.pair_band_max))!=len(a.pair_band_max):p.error('choose distinct pair cutoffs')
    vbm=float(wave.energies[:,17].max());cbm=float(wave.energies[:,18].min())
    if cbm<=vbm:p.error('the MoS2 example requires a gap above occupied band18')
    out.mkdir(parents=True)
    record=dict(status='RUNNING',mesh=[mesh,mesh],source_nbands=nb,occupied=18,
        pair_band_max=a.pair_band_max,vbm_eV=vbm,cbm_eV=cbm,mu_reference_eV=(vbm+cbm)/2,
        temperatures_K=[0,300],region_radius_inv_A=.35,
        source_wavecar_sha256=sha256(source/'WAVECAR'),commands=[],
        scope='Matched eigenstates and density; two charge-current operators. No material or joint convergence claim.')
    def save():
        (out/'result.json').write_text(json.dumps(record,indent=2)+'\n')
    env=dict(os.environ)
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        env[key]='1'
    def run(name,argv,cwd=None):
        folder=out/'logs'/name;folder.mkdir(parents=True)
        started=time.monotonic()
        with (folder/'stdout.log').open('w') as stdout,(folder/'stderr.log').open('w') as stderr:
            done=subprocess.run(argv,cwd=cwd or ROOT,env=env,stdout=stdout,stderr=stderr)
        record['commands'].append(dict(name=name,argv=argv,returncode=done.returncode,
                                       seconds=time.monotonic()-started));save()
        if done.returncode:raise RuntimeError(name+' failed; inspect its logs')
    cli=[sys.executable,str(ROOT/'tools/vaspberry_kubo.py')]
    try:
        save()
        run('physical-export',cli+['spin-export','--run-dir',str(source),'--output-dir',str(out/'physical')])
        run('physical-pairs',cli+['velocity-pairs','--matrices',str(out/'physical/physical-matrices.npz'),
            '--metadata',str(out/'physical/physical-matrices.json'),'--mesh',str(mesh),str(mesh),
            '--energy-reference','unchanged VASP eigenvalue zero','--output-dir',str(out/'paw-pairs')])
        raw=out/'native';raw.mkdir()
        run('native-export',[str(binary),'-f',str(source/'WAVECAR'),'-s','2','-kubo','2',
                             '-kubo_pairs','PAIRS.csv'],cwd=raw)
        run('native-pairs',cli+['import-pairs','--csv',str(raw/'PAIRS.csv'),'--wavecar',str(source/'WAVECAR'),
            '--spinor-components','2','--spin-multiplicity','1','--mesh',str(mesh),str(mesh),
            '--energy-reference','unchanged VASP eigenvalue zero','--output-dir',str(out/'canonical-pairs')])
        for operator in ('canonical','paw'):
            for cap in a.pair_band_max:
                name=f'{operator}-cap{cap}'
                run(name,cli+['pair-hall','--pairs-dir',str(out/f'{operator}-pairs'),
                    '--mu-min',str(vbm-.2),'--mu-max',str(vbm+.1),'--mu-num','61',
                    '--mu-reference',str((vbm+cbm)/2),'--temperatures','0','300',
                    '--regions',str(HERE.parent/'regions.json'),'--difference','valley:K:Kprime',
                    '--degeneracy-policy','coalesce','--degeneracy-threshold-eV','1e-7',
                    '--pair-band-max',str(cap),'--formats','csv','dat','npz','--output-dir',str(out/name)])
        if sha256(source/'WAVECAR')!=record['source_wavecar_sha256']:
            raise RuntimeError('source WAVECAR changed during comparison')
        record['status']='PASS';save()
    except BaseException as e:
        record.update(status='FAILED',error=type(e).__name__+': '+str(e));save();raise
    print(json.dumps({'status':record['status'],'output':str(out)}))


if __name__=='__main__':main()
