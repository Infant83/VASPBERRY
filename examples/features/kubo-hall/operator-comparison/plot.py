#!/usr/bin/env python3
"""Plot matched charge operators using their unchanged numerical Hall tables."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',type=Path,default=Path(__file__).resolve().parent/'reference')
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();root=a.input_dir
    info=json.loads((root/'result.json').read_text())
    if info['status']!='PASS':raise ValueError('completed comparison required')
    curves={}
    for op in ('canonical','paw'):
        for cap in info['pair_band_max']:
            with (root/f'{op}-cap{cap}/conductivity.csv').open() as f:
                rows=list(csv.DictReader(f))
            curves[op,cap]={region:np.array([[float(r['mu_eV'])-info['vbm_eV'],
                    float(r['delta_sigma_e2_over_h']),float(r['sigma_e2_over_h'])]
                for r in rows if float(r['temperature_K'])==300 and r['region']==region])
                for region in ('K','Kprime','valley','total')}
    first=info['pair_band_max'][0]
    plt.rcParams.update({'font.size':12, 'axes.titlesize':13,
                         'axes.labelsize':12, 'xtick.labelsize':11, 'ytick.labelsize':11})
    fig,ax=plt.subplots(1,2,figsize=(9.8,3.3),layout='constrained')
    colors={'canonical':'#2166ac','paw':'#b2182b'}
    labels={'canonical':'WAVECAR momentum','paw':'PAW full velocity'}
    for op in ('canonical','paw'):
        for region,style in [('K','-'),('Kprime','--')]:
            values=curves[op,first][region]
            ax[0].plot(values[:,0],values[:,1],style,c=colors[op],lw=1.5,
                       label=labels[op]+(' · K' if region=='K' else ' · K′'))
        for i,cap in enumerate(info['pair_band_max']):
            values=curves[op,cap]['valley']
            ax[1].plot(values[:,0],values[:,1],c=colors[op],ls=('-', '--', ':', '-.')[i%4],
                       lw=1.5,label=f'{labels[op]}, M={cap}')
    for axis in ax:
        axis.axhline(0,c='0.7',lw=.6);axis.axvline(0,c='0.7',lw=.6)
        axis.set_xlabel(r'$\mu-E_v$ (eV)');axis.tick_params(direction='in')
        axis.spines[['top','right']].set_visible(False)
        axis.legend(frameon=False,fontsize=11,labelspacing=.25,borderaxespad=.25)
    ax[0].set(title=f'(a) Same states and $M={first}$',
              ylabel=r'Regional $\Delta\sigma_{xy}$ ($e^2/h$)')
    ax[1].set(title='(b) Operator and retained-band comparison',
              ylabel=r'$\Delta\sigma_K-\Delta\sigma_{K^\prime}$ ($e^2/h$)')
    a.output_dir.mkdir(parents=True,exist_ok=False)
    for fmt in ('png','pdf','svg'):
        fig.savefig(a.output_dir/f'mos2-operator-comparison.{fmt}',dpi=220)
    plt.close(fig)
    (a.output_dir/'plot.json').write_text(json.dumps(dict(status='PASS',temperature_K=300,
        source_nbands=info['source_nbands'],pair_band_max=info['pair_band_max'],
        interpretation='Two operators on the same mesh/eigenstates; no k-convergence comparison.'),indent=2)+'\n')


if __name__=='__main__':main()
