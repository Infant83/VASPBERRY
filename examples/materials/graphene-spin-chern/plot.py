#!/usr/bin/env python3
"""Plot completed native spin-Chern CSVs, optionally beside supplied VASP bands."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'tools'))
from plot_berry_curvature import reciprocal_from_poscar,draw_plaquette_map

def read_csv(path):
    text=path.read_text()
    meta=dict(line[2:].split('=',1) for line in text.splitlines() if line.startswith('# ') and '=' in line)
    rows=list(csv.DictReader(line for line in text.splitlines() if line.strip() and not line.startswith('#')))
    return meta,rows

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',type=Path,default=HERE/'reference/native-12x12')
    p.add_argument('--prefix',default='SPIN',help='Native --output prefix basename')
    p.add_argument('--poscar',type=Path,default=HERE/'inputs/POSCAR')
    band_input=p.add_mutually_exclusive_group()
    band_input.add_argument('--bands-csv',type=Path,help='Optional k_index,q1,q2,q3,E1_eV,... table exported from WAVECAR')
    band_input.add_argument('--path-wavecar',type=Path,help='Optional SOC path WAVECAR; copy its unchanged eigenvalues into bands.csv for plotting')
    p.add_argument('--path-node-indices',type=int,nargs='+')
    p.add_argument('--path-labels',nargs='+')
    p.add_argument('--gap-k-index',type=int,help='Optional one-based band-path point for a micro-eV gap inset')
    p.add_argument('--k-fractional',type=float,nargs=2,default=[1/3,1/3],help='K label in this reciprocal basis; change for another lattice')
    p.add_argument('--title',default='Graphene: projected-spin Chern sectors')
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--formats',nargs='+',choices=['png','pdf','svg'],default=['png','pdf','svg'])
    a=p.parse_args()
    if a.path_wavecar:
        from wavecar_fukui import Wavecar
        wave=Wavecar(a.path_wavecar,spinor_components=2)
        a.output_dir.mkdir(parents=True,exist_ok=True)
        a.bands_csv=a.output_dir/'bands.csv'
        with a.bands_csv.open('w',newline='') as stream:
            writer=csv.writer(stream)
            writer.writerow(['k_index','q1','q2','q3',*[f'E{i+1}_eV' for i in range(wave.header.nbands)]])
            writer.writerows([[i+1,*q,*wave.energies[i]] for i,q in enumerate(wave.kpoints)])
    files={name:a.input_dir/f'{a.prefix}_{name}.csv' for name in ['CHERN','BERRY','SPECTRUM']}
    tables={name:read_csv(path) for name,path in files.items()}
    for name,(meta,rows) in tables.items():
        if meta.get('schema')!='VASPBERRY_SPIN_CHERN_V1' or meta.get('result_status')!='PASS':
            p.error(name+' is not a completed native spin-Chern result')
        if not rows:p.error(name+' contains no numerical rows')
    cm,cr=tables['CHERN'];bm,br=tables['BERRY'];sm,sr=tables['SPECTRUM']
    if cm!=bm or cm!=sm:p.error('Native CHERN/BERRY/SPECTRUM metadata differ; do not mix calculations')
    if len(cr)!=1:p.error('CHERN must contain exactly one native result row')
    result={k:float(v) for k,v in cr[0].items()}
    q=np.array([[float(row[k]) for k in ['k1','k2','k3']] for row in br])
    flux=np.array([[float(row[k]) for k in ['flux_plus_rad','flux_minus_rad']] for row in br])
    vals=np.array([float(row['projected_pauli_eigenvalue']) for row in sr])
    indices=np.array([int(row['k_index']) for row in sr])
    if not all(np.isfinite(x).all() for x in [q,flux,vals,np.array(list(result.values()))]):
        p.error('Nonfinite native output')
    if not np.allclose(flux.sum(axis=0)/(2*np.pi),[result['C_plus'],result['C_minus']],atol=1e-8,rtol=0):
        p.error('Native flux tables do not match the completed native Chern summary')
    reciprocal=reciprocal_from_poscar(a.poscar)
    ncols=2 if a.bands_csv else 3
    fig,axs=plt.subplots(2 if a.bands_csv else 1,ncols,figsize=(10.5,8.5) if a.bands_csv else (12.5,3.8),constrained_layout=True)
    axs=np.asarray(axs).ravel();bound=float(abs(flux).max())
    plot_records=[]
    for axis,j,label in [(axs[0],0,'+'),(axs[1],1,'−')]:
        artist,record=draw_plaquette_map(axis,q,flux[:,j],reciprocal,vmax=bound,k_fractional=a.k_fractional)
        axis.set_title(f"({'ab'[j]}) {label} sector, C{label} = {result['C_plus' if j==0 else 'C_minus']:.6g}")
        fig.colorbar(artist,ax=axis,shrink=.80,label='Plaquette flux (rad)')
        plot_records.append(record)
    spin_ax=axs[-1]
    spin_ax.scatter(indices,vals,s=7,c=np.where(vals>0,'#2166ac','#b2182b'),alpha=.6,rasterized=True)
    spin_ax.axhline(0,color='.4',lw=.7)
    spin_ax.set(xlabel='Full-mesh k-point index',ylabel='Projected Pauli eigenvalue',ylim=(-1.08,1.08),
                title=f"({'d' if a.bands_csv else 'c'}) Projected spin gap")
    spin_ax.text(.04,.50,f"min |λ| = {float(cm['min_abs_spin_eigenvalue']):.6f}\nCspin = {result['C_spin']:.6g}",transform=spin_ax.transAxes,va='center',fontsize=10)
    extra={}
    if a.bands_csv:
        if int(cm['selected_bands'].split(':')[0])!=1:
            p.error('Optional occupied-band panel requires a selected subspace beginning at band 1')
        _,bands=read_csv(a.bands_csv)
        energy_columns=sorted([k for k in bands[0] if k.startswith('E') and k.endswith('_eV')],key=lambda k:int(k[1:-3]))
        energies=np.array([[float(row[k]) for k in energy_columns] for row in bands])
        coords=np.array([[float(row[k]) for k in ['q1','q2','q3']] for row in bands])
        occupied=int(cm['selected_bands'].split(':')[1])
        origin=.5*(energies[:,occupied-1].max()+energies[:,occupied].min())
        distance=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(coords,axis=0)@reciprocal,axis=1))]
        ax=axs[2]
        ax.plot(distance,energies[:,:occupied]-origin,color='#2166ac',lw=.8)
        ax.plot(distance,energies[:,occupied:]-origin,color='.55',lw=.7)
        ax.axhline(0,color='.4',lw=.5)
        ax.set(xlim=(0,distance[-1]),ylim=(-8,5),xlabel='Wave vector',ylabel='E − midgap (eV)',title='(c) VASP band structure')
        if a.path_node_indices or a.path_labels:
            if not a.path_node_indices or not a.path_labels or len(a.path_node_indices)!=len(a.path_labels):p.error('Provide matching path indices and labels')
            nodes=np.array(a.path_node_indices)-1
            if min(nodes)<0 or max(nodes)>=len(distance):p.error('Path node index outside band table')
            ax.set_xticks(distance[nodes],a.path_labels)
        if a.gap_k_index is not None:
            ik=a.gap_k_index-1
            if not 0<=ik<len(energies):p.error('Gap point is outside band table')
            gap=energies[ik,occupied]-energies[ik,occupied-1]
            center=.5*(energies[ik,occupied]+energies[ik,occupied-1])
            inset=ax.inset_axes([.65,.12,.30,.30]);levels=(energies[ik,occupied-2:occupied+2]-center)*1e6
            inset.hlines(levels,[.10,.10,.55,.55],[.45,.45,.90,.90],color=['#2166ac','#2166ac','.4','.4'],lw=2)
            inset.set(xticks=[],ylabel='µeV',title=f'K gap {gap*1e6:.3f} µeV')
            inset.tick_params(labelsize=7);inset.title.set_fontsize(8);inset.yaxis.label.set_size(8)
            extra['displayed_gap_eV']=float(gap)
        extra.update(bands_sha256=sha(a.bands_csv),band_energy_origin_eV=float(origin),band_origin_scope='midgap of supplied path samples')
    fig.suptitle(a.title,fontsize=15)
    a.output_dir.mkdir(parents=True,exist_ok=True)
    outputs=[]
    for suffix in a.formats:
        target=a.output_dir/f'figure.{suffix}';fig.savefig(target,dpi=240);outputs.append(target)
    plt.close(fig)
    record=dict(status='PASS',mode='plot_completed_native_outputs_only',native_input_sha256={k:sha(v) for k,v in files.items()},
                poscar_sha256=sha(a.poscar),native_result=result,maps=plot_records,**extra,
                output_sha256={x.name:sha(x) for x in outputs})
    if a.path_wavecar:record['path_wavecar_sha256']=sha(a.path_wavecar)
    (a.output_dir/'plotting.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))

if __name__=='__main__':main()
