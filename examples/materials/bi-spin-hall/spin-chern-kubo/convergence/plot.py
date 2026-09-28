#!/usr/bin/env python3
"""Plot saved Bi spin-sector convergence tables; no electronic calculation."""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--input',type=Path,default=Path(__file__).resolve().parent/'reference',help='Directory containing convergence.csv')
 p.add_argument('--output',type=Path,required=True,help='New figure directory')
 a=p.parse_args()
 if a.output.exists():p.error('output directory exists; choose a new directory')
 source=a.input/'convergence.csv'
 with source.open() as f:rows=list(csv.DictReader(f))
 for row in rows:
  for key in ['mesh','source_nbands','sum_band_max']:row[key]=int(row[key])
  for key in ['C_est_spin','Fukui_C_spin','Gamma_point_contribution']:row[key]=float(row[key])
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
 fig,axes=plt.subplots(2,2,figsize=(10,7.2),layout='constrained')
 for col,(bands,title,color) in enumerate([('1:10','Occupied bands 1–10','#196c93'),('9:10','Selected pair 9–10','#ad4329')]):
  ax=axes[0,col];selected=sorted([r for r in rows if r['series']=='mesh' and r['bands']==bands],key=lambda r:r['mesh'])
  ax.plot([r['mesh'] for r in selected],[r['C_est_spin'] for r in selected],'o-',color=color,label='Canonical proxy, sum bands ≤48')
  ax.axhline(selected[0]['Fukui_C_spin'],color='black',ls='--',lw=1.2,label='Spin Chern number (Fukui)')
  if col==1:ax.plot([r['mesh'] for r in selected],[r['Gamma_point_contribution'] for r in selected],'x:',color='#777777',label='One sampled Γ-point contribution')
  ax.set(xlabel='N in N×N mesh',ylabel='Raw spin-sector integral',title=title,xticks=[6,12,18])
  ax.legend(fontsize=8,loc='best');ax.grid(alpha=.18)
  ax=axes[1,col];selected=sorted([r for r in rows if r['series']=='cutoff' and r['bands']==bands],key=lambda r:r['sum_band_max'])
  good=[r for r in selected if r['sum_band_max']<=48];upper=[r for r in selected if r['sum_band_max']>48]
  ax.plot([r['sum_band_max'] for r in selected],[r['C_est_spin'] for r in selected],'-',color=color,lw=1)
  ax.plot([r['sum_band_max'] for r in good],[r['C_est_spin'] for r in good],'o',color=color,label='Retained ≤48')
  ax.plot([r['sum_band_max'] for r in upper],[r['C_est_spin'] for r in upper],'o',mfc='white',mec=color,label='Includes highest empty states')
  ax.set(xlabel='Response sum cutoff (same 80-band WAVECAR)',ylabel='Raw spin-sector integral',title='Fixed 6×6 source, '+title.lower(),xticks=[16,24,32,48,64,80]);ax.ticklabel_format(axis='y',style='plain',useOffset=False)
  ax.legend(fontsize=8,loc='best');ax.grid(alpha=.18)
 fig.suptitle('Bi: sampling and source-band controls remain distinct',fontsize=14)
 fig.supxlabel('Pseudo-wavefunction canonical derivative proxy; no rounding or full PAW certification.',fontsize=10)
 a.output.mkdir(parents=True)
 for suffix in ['png','pdf','svg']:fig.savefig(a.output/('figure.'+suffix),dpi=190)
 plt.close(fig)
 sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
 (a.output/'plotting.json').write_text(json.dumps({'input_sha256':sha(source),'script_sha256':sha(Path(__file__)),'outputs':{p.name:sha(p) for p in a.output.glob('figure.*')},'scope':'Display of saved tables; straight segments connect actual mesh/cutoff tests. No interpolation model or new physical calculation.'},indent=2)+'\n')
 print(a.output)

if __name__=='__main__':main()
