#!/usr/bin/env python3
"""Draw Bi bulk/edge spectra and spin Hall convergence from actual outputs."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.path import Path as PlotPath
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[2]/'tools'))
from plot_berry_panels import periodic_bilinear
from plot_berry_curvature import first_bz_polygon, draw_bz_outline


def load(path):
    with np.load(path,allow_pickle=False) as z:return {k:z[k] for k in z.files}


def save(fig,out,name):
    for fmt in ('png','pdf','svg'):
        fig.savefig(out/(name+'.'+fmt),dpi=220,bbox_inches='tight',facecolor='white')
    plt.close(fig)


def path_positions(q,reciprocal,vertices,ticks):
    result=[]
    for a,b,offset in zip(vertices[:-1],vertices[1:],ticks[:-1]):
        direction=(b-a)@reciprocal;length=np.linalg.norm(direction)
        for i in (-1,0,1):
            for j in (-1,0,1):
                d=(q+[i,j,0]-a)@reciprocal;t=np.dot(d,direction)/length**2
                if -1e-9 <= t <= 1+1e-9 and np.linalg.norm(d-t*direction)<1e-7:
                    x=offset+t*length
                    if not any(abs(x-v)<1e-7 for v in result):result.append(x)
    return result


def spectra(ref,out,meta):
    bulk=load(ref/'bands/bands.npz');edge=load(ref/meta['edge_case']/'edge.npz')
    direct=load(ref/'direct-dft/bands.npz')
    mid=meta['reference_midgap_eV'];occupied=meta['model_occupied']
    reciprocal=2*np.pi*np.linalg.inv(bulk['lattice_A']).T
    fig,(ax,ae)=plt.subplots(1,2,figsize=(9.0,3.8),layout='constrained')
    energy=bulk['energies_eV']-mid;distance=bulk['distance_inv_A'];ticks=bulk['tick_distance_inv_A']
    for b in range(energy.shape[1]):
        ax.plot(distance,energy[:,b],lw=.95,color='#23699A' if b<occupied else '#777777')
    for q,energies in zip(direct['kpoints_fractional'],direct['energies_eV']):
        for x in path_positions(q,reciprocal,bulk['vertices_fractional'],ticks):
            selected=energies-mid
            selected=selected[(selected>=-1.2)&(selected<=1.2)]
            ax.scatter(np.full(len(selected),x),selected,s=15,facecolors='none',edgecolors='#262626',lw=.65,zorder=5)
    ax.set(xlim=(0,ticks[-1]),ylim=(-1.2,1.2),xticks=ticks,
           xticklabels=[r'$\Gamma$','M','K',r'$\Gamma$'],ylabel='Energy − midgap (eV)')
    for x in ticks[1:-1]:ax.axvline(x,color='.85',lw=.6,zorder=0)
    ax.axhline(0,color='.75',lw=.7,zorder=0)
    ax.legend(handles=[Line2D([],[],color='#23699A',lw=1.2,label='VASP-derived Wannier'),
        Line2D([],[],marker='o',ms=4,ls='',mfc='none',mec='.2',label='Direct VASP')],loc='lower left',fontsize=9.5,frameon=False)
    # Sum all states before coloring: invariant under degenerate-state rotations.
    eta=meta.get('edge_display_broadening_eV',.01)
    grid=np.linspace(-1.2,1.2,601)
    ldos=np.empty((len(grid),len(edge['q_parallel'])))
    for k,(e,w) in enumerate(zip(edge['energies_eV']-mid,edge['left_edge_weight'])):
        ldos[:,k]=np.sum(w[None,:]*np.exp(-.5*((grid[:,None]-e[None,:])/eta)**2),axis=1)/(np.sqrt(2*np.pi)*eta)
    edge_meta=json.loads((ref/meta['edge_case']/'edge.json').read_text())
    a=np.linalg.norm(edge['lattice_A'][edge_meta['periodic_axis']])
    parallel=2*np.pi*edge['q_parallel']/a
    vmax=float(ldos.max());vmin=max(vmax*1e-4,1e-5)
    image=ae.pcolormesh(parallel,grid,np.maximum(ldos,vmin),shading='auto',cmap='magma',norm=LogNorm(vmin=vmin,vmax=vmax),rasterized=True)
    ae.axhline(0,color='white',lw=.6,ls='--',alpha=.6)
    ae.set(xlabel=r'$k_\parallel$ (Å$^{-1}$)',ylabel='Energy − midgap (eV)',xlim=(parallel[0],parallel[-1]),ylim=(grid[0],grid[-1]))
    c=fig.colorbar(image,ax=ae,pad=.02);c.set_label(r'Left-edge spectral weight (eV$^{-1}$)',fontsize=10)
    ax.set_title('(a) Bi band structure',loc='left');ae.set_title(f'(b) Edge check of $Z_2=1$, {edge_meta["width_cells"]} cells',loc='left')
    save(fig,out,'bi-bulk-edge')


def response(ref,out,meta):
    data=load(ref/meta['spin_case']/'spin_hall.npz')
    with (ref/'convergence.csv').open() as f:rows=list(csv.DictReader(f))
    if not rows:raise ValueError('Actual convergence data required')
    mesh_rows=[r for r in rows if r['study']=='mesh']
    band_rows=[r for r in rows if r['study']=='bands']
    if len({(r['source_nbands'],r['retained_bands']) for r in mesh_rows})!=1:
        raise ValueError('Mesh curve requires the same source and retained band counts')
    if len({r['mesh'] for r in band_rows})!=1:
        raise ValueError('Band curves require a single fixed mesh')
    keys=[(r['study'],r['mesh'],r['source_nbands'],r['retained_bands']) for r in rows]
    if len(set(keys))!=len(keys):raise ValueError('Duplicate convergence case')
    reciprocal=2*np.pi*np.linalg.inv(data['lattice_A']).T
    polygon=first_bz_polygon(reciprocal)
    xmin,ymin=polygon.min(axis=0);xmax,ymax=polygon.max(axis=0)
    x=np.linspace(xmin,xmax,401);y=np.linspace(ymin,ymax,401);xx,yy=np.meshgrid(x,y)
    xy=np.column_stack((xx.ravel(),yy.ravel()))
    q=np.column_stack((xy@np.linalg.inv(reciprocal[:2,:2]),np.zeros(len(xy))))
    z=periodic_bilinear(data['kpoints_fractional'],data['omega_spin_A2'][:,2,0,1],q)
    z[~PlotPath(polygon).contains_points(xy,radius=1e-10)]=np.nan
    limit=float(np.max(abs(data['omega_spin_A2'][:,2,0,1])))
    fig,axes=plt.subplots(1,3,figsize=(11.5,3.6),layout='constrained',gridspec_kw={'width_ratios':[1.15,1,1]})
    ax,ak,ab=axes
    image=ax.pcolormesh(xx,yy,z.reshape(xx.shape),cmap='RdBu_r',vmin=-limit,vmax=limit,shading='auto',rasterized=True)
    draw_bz_outline(ax,polygon,reciprocal,k_fractional=(2/3,1/3))
    ax.tick_params(top=False,right=False)
    ax.set(xlabel=r'$k_x$ (Å$^{-1}$)',ylabel=r'$k_y$ (Å$^{-1}$)')
    cb=fig.colorbar(image,ax=ax,pad=.025,fraction=.055);cb.set_label(r'$\Omega^z_{xy}$ (Å$^2$)',fontsize=10)
    mesh=[r for r in rows if r['study']=='mesh']
    mesh.sort(key=lambda r:int(r['mesh']))
    ak.plot([int(r['mesh']) for r in mesh],[float(r['sigma_zxy']) for r in mesh],'o-',color='#176B53',ms=5)
    ak.set(xlabel=r'Full $N\times N$ mesh: $N$',ylabel=r'$\sigma^z_{xy}$ [$(\hbar/e)(e^2/h)$]')
    ak.set_xticks([int(r['mesh']) for r in mesh])
    for source in sorted({int(r['source_nbands']) for r in rows if r['study']=='bands'}):
        selected=sorted([r for r in rows if r['study']=='bands' and int(r['source_nbands'])==source],key=lambda r:int(r['retained_bands']))
        ab.plot([int(r['retained_bands']) for r in selected],[float(r['sigma_zxy']) for r in selected],'o-',ms=4,label=f'VASP NBANDS={source}')
    ab.set(xlabel='Retained source bands',ylabel=r'$\sigma^z_{xy}$ [$(\hbar/e)(e^2/h)$]')
    ab.legend(frameon=False,fontsize=8.5)
    for a,title in zip(axes,['(a) Spin Berry curvature','(b) k-mesh refinement','(c) Band-space convergence']):a.set_title(title,loc='left',fontsize=10.5)
    save(fig,out,'bi-spin-hall')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-dir',type=Path,default=HERE/'reference')
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args()
    if a.output_dir.exists():p.error('Output directory exists; choose a new directory')
    meta=json.loads((a.reference_dir/'summary.json').read_text())
    if meta.get('complete') is not True:raise ValueError('Completed actual Bi reference required')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10.5,'axes.labelsize':11,
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42})
    a.output_dir.mkdir(parents=True)
    spectra(a.reference_dir,a.output_dir,meta);response(a.reference_dir,a.output_dir,meta)
    (a.output_dir/'plot.json').write_text(json.dumps(dict(complete=True,source_summary=meta,
        curvature_display='Periodic bilinear display interpolation only; integrals use the original points.',
        edge_display='Gaussian spectral broadening of original strip eigenvalues; all degenerate-state weights summed.'),indent=2)+'\n')


if __name__ == '__main__':main()
