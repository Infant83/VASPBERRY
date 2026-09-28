#!/usr/bin/env python3
"""Draw completed native spin-sector Kubo proxies beside Chern numbers evaluated with the Fukui method."""
import argparse,csv,hashlib,json
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from matplotlib.ticker import SymmetricalLogLocator
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
sys.path.insert(0,str(ROOT/'tools'))
from plot_berry_curvature import reciprocal_from_poscar,draw_plaquette_map

SCHEMA='VASPBERRY_SPIN_SECTOR_KUBO_V1'

def read(path,schema):
 lines=path.read_text().splitlines()
 meta=dict(x[2:].split('=',1) for x in lines if x.startswith('# ') and '=' in x)
 rows=list(csv.DictReader(x for x in lines if x and not x.startswith('#')))
 if meta.get('schema')!=schema or meta.get('result_status')!='PASS' or not rows:
  raise ValueError('Expected completed '+schema+' data: '+str(path))
 return meta,rows

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def numeric(rows,key):
 a=np.array([float(row[key]) for row in rows])
 if not np.isfinite(a).all():raise ValueError('Nonfinite '+key)
 return a

def same_space(a,b):
 for key in ['selected_bands','spin_axis_cartesian','spin_frame_row_1','spin_frame_row_2','spin_frame_row_3']:
  if a.get(key)!=b.get(key):raise ValueError('Input subspace metadata differ: '+key)

def integral(folder,prefix='SPIN'):
 m,r=read(folder/(prefix+'_KUBO_INTEGRAL.csv'),SCHEMA)
 if len(r)!=1:raise ValueError('Expected one native integral row')
 values={k:float(v) for k,v in r[0].items()}
 if not np.isfinite(list(values.values())).all():raise ValueError('Nonfinite integral')
 return m,values

def fukui(folder,prefix='SPIN'):
 m,r=read(folder/(prefix+'_CHERN.csv'),'VASPBERRY_SPIN_CHERN_V1')
 if len(r)!=1:raise ValueError('Expected one native Chern-number row')
 return m,{k:float(v) for k,v in r[0].items()}

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--input-dir',type=Path,help='Reference layout: occupied-{mesh,path,fukui} and optional pair-* directories')
 for name in ['mesh','path','fukui','pair-mesh','pair-path','pair-fukui']:
  p.add_argument('--'+name+'-dir',type=Path,help='Explicit native output directory; overrides --input-dir layout')
 p.add_argument('--prefix',default='SPIN')
 p.add_argument('--pair-prefix',default='SPIN')
 p.add_argument('--poscar',required=True,type=Path)
 p.add_argument('--title',default='Projected-spin sectors: differential proxy and Chern numbers (Fukui method)')
 p.add_argument('--path-node-indices',nargs='+',type=int,default=[1,17,33,49])
 p.add_argument('--path-labels',nargs='+',default=['Γ','K','M','Γ'])
 p.add_argument('--k-fractional',nargs=2,type=float,default=[1/3,1/3])
 p.add_argument('--linthresh',type=float,default=1.,help='Linear half-width of signed-log curvature scale, Å²')
 p.add_argument('--output-dir',required=True,type=Path)
 a=p.parse_args()
 if not np.isfinite(a.linthresh) or a.linthresh<=0:p.error('--linthresh must be finite and positive')
 files=[]
 folders={}
 for kind in ['occupied-mesh','occupied-path','occupied-fukui','pair-mesh','pair-path','pair-fukui']:
  attr=(kind.replace('occupied-','').replace('-','_'))+'_dir'
  folders[kind]=getattr(a,attr) or (a.input_dir/kind if a.input_dir else None)
 if folders['occupied-fukui'] is None:folders['occupied-fukui']=folders['occupied-mesh']
 if any(folders[x] is None for x in ['occupied-mesh','occupied-path']):p.error('Provide --input-dir or --mesh-dir and --path-dir')
 if a.input_dir is None and folders['pair-path'] is None and folders['pair-mesh'] is None:
  prefix='PAIR' if a.pair_prefix=='SPIN' else a.pair_prefix
  if all((folders['occupied-'+kind]/(prefix+'_KUBO.csv')).is_file() for kind in ['mesh','path']):
   a.pair_prefix=prefix
   for kind in ['mesh','path','fukui']:folders['pair-'+kind]=folders['occupied-'+kind]
 def file_for(kind,suffix):
  prefix=a.pair_prefix if kind.startswith('pair-') else a.prefix
  return folders[kind]/(prefix+'_'+suffix+'.csv')
 def point(folder):
  file=file_for(folder,'KUBO');files.append(file);m,r=read(file,SCHEMA)
  if m.get('derivative')!='CANONICAL_MOMENTUM_DERIVATIVE_PROXY':raise ValueError('Unexpected derivative representation')
  return m,r
 mesh_m,mesh=point('occupied-mesh');path_m,path=point('occupied-path');same_space(mesh_m,path_m)
 for key in ['source_nbands','source_encut_eV','derivative','representation']:
  if mesh_m.get(key)!=path_m.get(key):raise ValueError('Path/mesh source metadata differ: '+key)
 im,iv=integral(folders['occupied-mesh'],a.prefix);fm,fv=fukui(folders['occupied-fukui'],a.prefix);same_space(mesh_m,im);same_space(mesh_m,fm)
 if mesh_m!=im:raise ValueError('Proxy mesh/integral metadata differ')
 area=np.fromstring(im['oriented_cell_area_A-2'],sep=' ')
 expected=np.array([sum(numeric(mesh,'omega_'+sector+'_'+component+'_A2').sum()*area[j] for j,component in enumerate(['yz','zx','xy']))/(2*np.pi) for sector in ['plus','minus','parent']])
 observed=np.array([iv['C_est_plus'],iv['C_est_minus'],iv['C_est_parent']])
 if not np.allclose(expected,observed,atol=1e-8,rtol=1e-10):raise ValueError('Point values do not match raw native mesh integral')
 files += [file_for('occupied-mesh','KUBO_INTEGRAL'),file_for('occupied-fukui','CHERN')]
 bm,berry=read(file_for('occupied-fukui','BERRY'),'VASPBERRY_SPIN_CHERN_V1');files.append(file_for('occupied-fukui','BERRY'))
 if bm!=fm:raise ValueError('FHS flux/summary metadata differ')
 if mesh_m.get('mesh')!=fm.get('mesh'):raise ValueError('Proxy and FHS mesh sizes differ')
 comparison=[dict(subspace='occupied',bands=mesh_m['selected_bands'],C_est_plus=iv['C_est_plus'],C_est_minus=iv['C_est_minus'],C_est_spin=iv['C_est_spin'],C_Fukui_plus=fv['C_plus'],C_Fukui_minus=fv['C_minus'],C_Fukui_spin=fv['C_spin'],fukui_status='PASS')]
 fig,axes=plt.subplots(2,2,figsize=(11,9.2));fig.subplots_adjust(left=.08,right=.96,top=.89,bottom=.18,hspace=.37,wspace=.35)
 def path_panel(ax,meta,rows,label):
  x=numeric(rows,'k_distance_A-1')
  for sector,color in [('plus','#2166ac'),('minus','#b2182b')]:ax.plot(x,numeric(rows,'omega_'+sector+'_xy_A2'),lw=1.2,marker='o',markersize=2.5,label='+' if sector=='plus' else '−',color=color)
  ax.set_yscale('symlog',linthresh=a.linthresh);locator=SymmetricalLogLocator(base=10,linthresh=a.linthresh);locator.set_params(numticks=9);ax.yaxis.set_major_locator(locator);ax.axhline(0,lw=.5,color='.6');ax.set(xlabel='Wave vector',ylabel='Ωxy proxy (Å²; signed log)',title=label+' bands '+meta['selected_bands'])
  indices=np.array(a.path_node_indices)-1
  if len(indices)!=len(a.path_labels) or min(indices)<0 or max(indices)>=len(x):raise ValueError('Invalid path node indices/labels')
  ax.set_xticks(x[indices],a.path_labels);ax.set_xlim(x[0],x[-1]);ax.legend(title='Spin sector',fontsize=8,title_fontsize=8)
 path_panel(axes[0,0],path_m,path,'(a) Occupied path')
 pair_exists=folders['pair-path'] is not None and file_for('pair-path','KUBO').exists()
 if pair_exists:
  pm,pr=point('pair-path');path_panel(axes[0,1],pm,pr,'(b) Selected pair path')
  pim,piv=integral(folders['pair-mesh'],a.pair_prefix);same_space(pm,pim);files.append(file_for('pair-mesh','KUBO_INTEGRAL'))
  result=dict(subspace='pair',bands=pm['selected_bands'],C_est_plus=piv['C_est_plus'],C_est_minus=piv['C_est_minus'],C_est_spin=piv['C_est_spin'])
  if folders['pair-fukui'] is not None and file_for('pair-fukui','CHERN').exists():
   pfm,pfv=fukui(folders['pair-fukui'],a.pair_prefix);same_space(pm,pfm);result.update(C_Fukui_plus=pfv['C_plus'],C_Fukui_minus=pfv['C_minus'],C_Fukui_spin=pfv['C_spin'],fukui_status='PASS');files.append(file_for('pair-fukui','CHERN'))
  else:result.update(C_Fukui_plus=None,C_Fukui_minus=None,C_Fukui_spin=None,fukui_status='NO_COMPLETED_INVARIANT')
  comparison.append(result)
 else:axes[0,1].axis('off')
 reciprocal=reciprocal_from_poscar(a.poscar)
 q=np.column_stack([numeric(mesh,k) for k in ['k1','k2','k3']]);values=numeric(mesh,'omega_plus_xy_A2');bound=max(float(abs(values).max()),a.linthresh)
 artist,map_proxy=draw_plaquette_map(axes[1,0],q,values,reciprocal,norm=SymLogNorm(a.linthresh,vmin=-bound,vmax=bound),k_fractional=a.k_fractional)
 colorbar=fig.colorbar(artist,ax=axes[1,0],shrink=.8,label='Ω+ proxy (Å²; signed log)');locator=SymmetricalLogLocator(base=10,linthresh=a.linthresh);locator.set_params(numticks=9);colorbar.locator=locator;colorbar.update_ticks();axes[1,0].set_title('(c) Occupied full-BZ point samples')
 bq=np.column_stack([numeric(berry,k) for k in ['k1','k2','k3']]);flux=numeric(berry,'flux_plus_rad')
 if abs(flux.sum()/(2*np.pi)-fv['C_plus'])>1e-8:raise ValueError('FHS flux sum and Chern-number summary disagree')
 artist,map_fukui=draw_plaquette_map(axes[1,1],bq,flux,reciprocal,k_fractional=a.k_fractional)
 fig.colorbar(artist,ax=axes[1,1],shrink=.8,label='Native + sector flux (rad)');axes[1,1].set_title('(d) Fukui method: occupied flux')
 table=[['Subspace','Raw C_est,+','Raw C_est,−','C+ (Fukui)','C− (Fukui)']]
 for c in comparison:table.append([c['subspace']+' '+c['bands'],f"{c['C_est_plus']:.7g}",f"{c['C_est_minus']:.7g}",'unassigned' if c['C_Fukui_plus'] is None else f"{c['C_Fukui_plus']:.7g}",'unassigned' if c['C_Fukui_minus'] is None else f"{c['C_Fukui_minus']:.7g}"])
 tx=fig.add_axes([.08,.06,.84,.08]);tx.axis('off');tab=tx.table(cellText=table[1:],colLabels=table[0],loc='center',cellLoc='center');tab.auto_set_font_size(False);tab.set_fontsize(9);tab.scale(1,1.3)
 fig.suptitle(a.title,fontsize=14);fig.text(.5,.93,'Canonical-momentum derivative proxy · pseudo Gram metric · no integer rounding',ha='center',fontsize=10);fig.text(.5,.02,'Lines connect stored points; cell colors are constant. C_est is not a certified topological invariant.',ha='center',fontsize=9)
 a.output_dir.mkdir(parents=True,exist_ok=True)
 for suffix in ['png','pdf','svg']:fig.savefig(a.output_dir/('figure.'+suffix),dpi=240)
 plt.close(fig)
 with (a.output_dir/'comparison.csv').open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(comparison[0]));writer.writeheader();writer.writerows(comparison)
 record=dict(status='PASS',source='completed_native_CSV_only',method='no_Chern_or_proxy_recalculation',comparison=comparison,input_sha256={str(p.relative_to(a.input_dir)) if a.input_dir and p.is_relative_to(a.input_dir) else str(p):sha(p) for p in files},poscar_sha256=sha(a.poscar),plot_script_sha256=sha(Path(__file__)),maps={'proxy':dict(**map_proxy,sample_location='native_k_node; constant centered display cell'),'fukui':map_fukui},outputs_sha256={n:sha(a.output_dir/n) for n in ['figure.png','figure.pdf','figure.svg','comparison.csv']})
 (a.output_dir/'plotting.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
if __name__=='__main__':main()
