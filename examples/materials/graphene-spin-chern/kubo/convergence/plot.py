#!/usr/bin/env python3
"""Plot completed local radial measurements; distinguish measured data and Dirac diagnostics."""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read(path):
    with path.open() as f: return list(csv.DictReader(f))


def write(path, rows):
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def describe(radial, loops):
    models=[];integrals=[]
    for valley in ['K','Kprime']:
        data=[r for r in radial if r['valley']==valley]
        center=next(r for r in data if float(r['radius_Ainv'])==0)
        mass=float(center['gap_eV'])/2
        radii=sorted(set(float(r['radius_Ainv']) for r in data if float(r['radius_Ainv'])>0))
        gap=np.array([np.mean([float(r['gap_eV']) for r in data if float(r['radius_Ainv'])==q]) for q in radii])
        # Large-r energy slope determines v independently of the curvature.
        q=np.array(radii); select=q>=q.max()/10
        q2=q[select]**2
        v2=float(np.dot(q2, (gap[select]/2)**2-mass**2)/np.dot(q2,q2))
        if v2<=0: raise ValueError('Local dispersion fit is not a massive Dirac cone')
        velocity=np.sqrt(v2);width=mass/velocity
        amplitude=float(center['omega_spin_xy_A2'])*2*width**2
        models.append(dict(valley=valley,gap_at_center_eV=2*mass,energy_slope_eV_A=velocity,q_width_Ainv=width,
                           native_center_omega_spin_A2=float(center['omega_spin_xy_A2']),
                           native_center_over_dispersion_Dirac_curvature=amplitude,
                           model_definition='diagnostic only: sqrt(m^2+v^2*q^2), v fitted to outer measured rings; positive chirality chosen to match the observed sector, not inferred from energy',
                           model_C_local_at_max_radius=.5*(1-width/np.sqrt(width**2+q[-1]**2))))
        omega=np.array([np.mean([float(r['omega_spin_xy_A2']) for r in data if float(r['radius_Ainv'])==radius]) for radius in q])
        # Angular average and finite radial log-trapezoid. No interpolation is
        # presented as a measured point, and the unsampled outer BZ is excluded.
        for stride in [2,1]:
            select=list(range(0,len(q),stride))
            if select[-1]!=len(q)-1:select.append(len(q)-1)
            rr=q[select];oo=omega[select]
            value=.5*float(center['omega_spin_xy_A2'])*rr[0]**2
            for j in range(len(rr)):
                if j: value+=(np.log(rr[j])-np.log(rr[j-1]))*(rr[j]**2*oo[j]+rr[j-1]**2*oo[j-1])/2
                integrals.append(dict(valley=valley,radial_stride=stride,radius_Ainv=rr[j],
                                      C_local_native_log_trapezoid=value,radial_nodes_used=j+1,
                                      scope='Approximate circular local integral; excludes the rest of the BZ and is not a Chern number'))
    return models,integrals


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,default=Path(__file__).resolve().parent/'reference')
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    if args.output_dir.exists():parser.error('Choose a fresh output directory')
    radial=read(args.input_dir/'radial.csv');loops=read(args.input_dir/'loops.csv')
    models,integrals=describe(radial,loops)
    args.output_dir.mkdir(parents=True)
    write(args.output_dir/'model_diagnostics.csv',models);write(args.output_dir/'integration_diagnostics.csv',integrals)
    fig,axs=plt.subplots(2,2,figsize=(11.4,8.1),constrained_layout=True)
    colors=['#235dab','#b74326']
    nonpositive=sum(float(r['omega_spin_xy_A2'])<=0 for r in radial)
    if not all(np.isfinite(float(r['omega_spin_xy_A2'])) and float(r['gap_eV'])>0 for r in radial): raise ValueError('Nonfinite curvature or nonpositive measured gap')
    for color,model in zip(colors,models):
        valley=model['valley'];label="K′" if valley=='Kprime' else 'K'
        data=[r for r in radial if r['valley']==valley and float(r['radius_Ainv'])>0]
        radii=sorted(set(float(r['radius_Ainv']) for r in data));q=np.array(radii)
        gap=np.array([np.mean([float(r['gap_eV']) for r in data if float(r['radius_Ainv'])==radius]) for radius in q])
        omega=np.array([np.mean([float(r['omega_spin_xy_A2']) for r in data if float(r['radius_Ainv'])==radius]) for radius in q])
        axs[0,0].loglog(q,gap*1e6,'o',ms=4,color=color,label=label+' actual VASP')
        dense=np.geomspace(q[0],q[-1],400);m=model['gap_at_center_eV']/2;v=model['energy_slope_eV_A'];q0=model['q_width_Ainv']
        axs[0,0].loglog(dense,2*np.sqrt(m*m+(v*dense)**2)*1e6,'--',lw=1,color=color,label=label+' dispersion fit')
        axs[0,1].plot(q,omega,'o-',ms=4,lw=1,color=color,label=label+' native proxy')
        axs[0,1].plot(dense,q0/(2*(q0*q0+dense*dense)**1.5),'--',lw=1,color=color,label=label+' Dirac shape (+ chirality)')
        for vertices,style in [(4,':'),(8,'-'),(16,'--')]:
            ll=sorted([r for r in loops if r['valley']==valley and int(r['vertices'])==vertices],key=lambda r:float(r['radius_Ainv']))
            if ll:axs[1,0].semilogx([float(r['radius_Ainv']) for r in ll],[float(r['C_local_spin']) for r in ll],style,marker='o',ms=3,color=color,label=f'{label}, {vertices} vertices')
        axs[1,0].semilogx(dense,.5*(1-q0/np.sqrt(q0*q0+dense*dense)),'--',color=color,lw=.8,alpha=.6,label=label+' diagnostic Dirac disk')
        for stride,style in [(2,':'),(1,'-')]:
            ii=[r for r in integrals if r['valley']==valley and r['radial_stride']==stride]
            axs[1,1].semilogx([r['radius_Ainv'] for r in ii],[r['C_local_native_log_trapezoid'] for r in ii],style,marker='o',ms=3,color=color,label=f'{label}, '+('coarse radial subset' if stride==2 else 'all radial nodes'))
    axs[0,1].set_xscale('log')
    axs[0,1].set_yscale('symlog',linthresh=1.) if nonpositive else axs[0,1].set_yscale('log')
    titles=['(a) Actual microscopic gap scale','(b) Native curvature and diagnostic model','(c) Actual WAVECAR polygon flux / 2π','(d) Local native point quadrature']
    ylabels=['Band 9–8 gap (µeV)','Spin-sector curvature (Å²)','Local spin-sector flux / 2π','Local proxy integral / 2π']
    for ax,title,ylabel in zip(axs.ravel(),titles,ylabels):
        ax.set_title(title,fontsize=11);ax.set_xlabel('Cartesian radius from valley (Å⁻¹)');ax.set_ylabel(ylabel)
        ax.grid(alpha=.2);ax.legend(fontsize=7.3,loc='best')
    axs[1,0].axhline(.5,color='.5',lw=.6);axs[1,1].axhline(.5,color='.5',lw=.6)
    fig.suptitle('Graphene intrinsic SOC: local resolution study\nFixed 520 eV / 16-band SCF density; neither local integration is a full-BZ Chern number',fontsize=12)
    for ext in ['png','pdf','svg']:fig.savefig(args.output_dir/('figure.'+ext),dpi=180)
    plt.close(fig)
    (args.output_dir/'plotting.json').write_text(json.dumps(dict(input_dir_label=args.input_dir.name,
      measured_data='Native canonical curvature points; independently computed actual-WAVECAR polygon overlaps',
      model_data='Dashed diagnostic massive-Dirac curves fitted only to measured energy dispersion',
      integration='Finite local polar quadrature; no full-BZ estimate or integer rounding',
      nonpositive_native_samples=nonpositive, native_curvature_scale='symlog' if nonpositive else 'log',
      model_chirality='Positive sign chosen to match the observed sector; energy dispersion alone does not determine chirality',
      models=models),indent=2)+'\n')
    print(json.dumps(models,indent=2))


if __name__=='__main__':main()
