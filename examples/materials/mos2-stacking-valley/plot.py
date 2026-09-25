#!/usr/bin/env python3
"""Plot VASP-derived MoS2 stacking bands and circular transition selectivity."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
CASES = ('monolayer', '1h-bilayer', '2h-bilayer', '3r-bilayer')
LABELS = ('Monolayer', '1H bilayer', '2H bilayer', '3R bilayer')
COLORS = ('#23699A', '#AF332E')
NATIVE_INTENSITY_RESOLUTION = 1e-4
NATIVE_ETA_UNCERTAINTY_LIMIT = .005
DISPLAY_PHOTON_RANGE_EV = (1.3,2.2)


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    if not rows: raise ValueError('Empty data table: '+str(path))
    result = {}
    for key in rows[0]:
        result[key] = np.array([np.nan if r[key] == '' else float(r[key]) for r in rows])
    return result


def load_case(directory):
    meta = json.loads((directory/'result.json').read_text())
    if meta.get('status') != 'PASS': raise ValueError('Completed material calculation required: '+str(directory))
    if meta.get('case') != directory.name:
        raise ValueError('The recorded material case does not match its figure label')
    for name in ('bands.csv','native-optical.csv'):
        if sha(directory/name) != meta['output_sha256'][name]:
            raise ValueError('Reference table changed since calculation: '+name)
    bands = read_csv(directory/'bands.csv'); native = read_csv(directory/'native-optical.csv')
    nk, nb = meta['nkpts'], meta['source_nbands']; nocc = meta['nocc']
    if len(bands['energy_eV']) != nk*nb or not np.array_equal(bands['band'], np.tile(np.arange(1,nb+1),nk)):
        raise ValueError('Band table does not preserve the complete source ordering')
    energy = bands['energy_eV'].reshape(nk,nb); distance = bands['path_distance_inv_A'].reshape(nk,nb)[:,0]
    q = np.column_stack([bands[f'k{axis}_fractional'].reshape(nk,nb)[:,0] for axis in 'xyz'])
    vertices = np.arange(5)*((nk-1)//4)
    if (nk-1)%4 or not np.allclose(q[vertices], [[0,0,0],[0,.5,0],[1/3,2/3,0],[0,0,0],[-1/3,-2/3,0]], atol=1e-8):
        raise ValueError('Unexpected reference symmetry path')
    meta['sampled_vbm_eV'] = float(energy[:,nocc-1].max())
    paw = None; paw_meta = None
    if (directory/'paw-optical/optical.json').exists():
        paw_meta = json.loads((directory/'paw-optical/optical.json').read_text())
        if not paw_meta.get('complete'): raise ValueError('Incomplete PAW optics data')
        if (paw_meta['initial_bands_1based'] != meta['optical_initial_bands']
                or paw_meta['final_bands_1based'] != meta['optical_final_bands']):
            raise ValueError('Native and PAW optical windows must match')
        if paw_meta['source_sha256']['WAVECAR'] != meta['source_sha256']['WAVECAR']:
            raise ValueError('Native and PAW optics must use the same original wavefunctions')
        if sha(directory/'paw-optical/spectra.npz') != paw_meta['output_sha256']['spectra.npz']:
            raise ValueError('PAW spectrum archive changed since calculation')
        with np.load(directory/'paw-optical/spectra.npz', allow_pickle=False) as archive:
            paw = {name:archive[name].copy() for name in archive.files}
    return {'meta':meta, 'bands':energy, 'distance':distance, 'q':q, 'vertices':vertices,
            'native':native, 'paw':paw, 'paw_meta':paw_meta}


def eta_values(plus, minus, floor, method):
    """Keep the rounded-channel ratio; only mask its display when unresolved."""
    total = plus+minus
    raw = np.divide(plus-minus,total,out=np.full_like(total,np.nan),where=total>0)
    strength_valid = total > floor
    uncertainty = (np.divide(NATIVE_INTENSITY_RESOLUTION,total,
                            out=np.full_like(total,np.inf),where=total>0)
                   if method == 'native' else np.zeros_like(total))
    valid = strength_valid & (uncertainty <= NATIVE_ETA_UNCERTAINTY_LIMIT)
    return raw, uncertainty, strength_valid, valid


def spectrum(data, k, method):
    table = data[method]; mask = table['k_index'] == k+1
    if method == 'native':
        x = table['photon_energy_eV'][mask]
        plus = table['left_spectrum_au'][mask]; minus = table['right_spectrum_au'][mask]
        floor = max(1e-4, .001*float(np.max(plus+minus)))
    else:
        x = table['photon_eV'][mask]
        plus = table['I_plus_A2_per_eV'][mask]; minus = table['I_minus_A2_per_eV'][mask]
        floor = .001*float(np.max(plus+minus))
    raw,_,_,valid = eta_values(plus,minus,floor,method)
    eta = np.where(valid,raw,np.nan)
    return x,plus,minus,eta


def write_display_diagnostics(cases, out):
    """Archive unmodified ratios and the reproducible figure-only gate."""
    with (out/'native-eta-display.csv').open('w',newline='') as stream:
        writer=csv.writer(stream)
        writer.writerow(['case','k_index','valley','photon_energy_eV','left_spectrum_au',
                         'right_spectrum_au','eta_raw','eta_rounding_absolute_bound',
                         'strength_valid','display_valid','eta_display'])
        for case,data in zip(CASES,cases):
            for valley,k in zip(('K','Kprime'),data['vertices'][[2,4]]):
                x,plus,minus,_=spectrum(data,int(k),'native')
                floor=max(1e-4,.001*float(np.max(plus+minus)))
                raw,bound,strong,valid=eta_values(plus,minus,floor,'native')
                for i,photon in enumerate(x):
                    writer.writerow([case,int(k)+1,valley,photon,plus[i],minus[i],
                                     raw[i] if np.isfinite(raw[i]) else '',
                                     bound[i] if np.isfinite(bound[i]) else '',
                                     int(strong[i]),int(valid[i]),raw[i] if valid[i] else ''])


def save(fig, out, name):
    for suffix in ('png','pdf','svg'):
        fig.savefig(out/f'{name}.{suffix}', dpi=220, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def compact(cases, out):
    fig,axes = plt.subplots(2,4,figsize=(8.8,5.2),layout='constrained',sharey='row')
    for j,(data,label) in enumerate(zip(cases,LABELS)):
        meta=data['meta']; ax=axes[0,j]; energy=data['bands']-meta['sampled_vbm_eV']; nocc=meta['nocc']
        for band in range(energy.shape[1]):
            if energy[:,band].max() >= -1.2 and energy[:,band].min() <= 2.7:
                ax.plot(data['distance'],energy[:,band],color=COLORS[0] if band<nocc else '#777777',lw=.75)
        ax.axhline(0,color='.65',lw=.65)
        ticks=data['distance'][data['vertices']]
        for value in ticks[1:-1]: ax.axvline(value,color='.88',lw=.6,zorder=0)
        ax.set(xticks=ticks,xticklabels=[r'$\Gamma$','M','K',r'$\Gamma$',"K′"],
               xlim=(0,ticks[-1]),ylim=(-1.2,2.7),title=f'({chr(97+j)}) {label}')
        if j==0: ax.set_ylabel('Energy − VBM (eV)')
        axes[1,j].axhline(0,color='.7',lw=.6)
        for k,color in zip(data['vertices'][[2,4]],COLORS):
            x,_,_,eta=spectrum(data,int(k),'native')
            axes[1,j].plot(x,eta,color=color,lw=1.4)
            if data['paw'] is not None:
                x,_,_,eta=spectrum(data,int(k),'paw')
                axes[1,j].plot(x,eta,color=color,lw=1.2,ls='--')
        axes[1,j].set(xlim=DISPLAY_PHOTON_RANGE_EV,ylim=(-1.08,1.08),yticks=[-1,0,1],xlabel='Photon energy (eV)',
                      title=f'({chr(101+j)}) Circular selectivity')
        if j==0: axes[1,j].set_ylabel(r'$\eta=(I_+-I_-)/(I_++I_-)$')
    handles=[Line2D([],[],color=COLORS[0],label='K'),Line2D([],[],color=COLORS[1],label="K′"),
             Line2D([],[],color='.3',lw=1.4,label='WAVECAR momentum')]
    if any(c['paw'] is not None for c in cases):
        handles.append(Line2D([],[],color='.3',lw=1.2,ls='--',label='Optional PAW optical'))
    fig.legend(handles=handles,loc='outside lower center',ncol=len(handles),frameon=False,fontsize=9)
    save(fig,out,'stacking-bands-selectivity')


def strengths(cases,out,method):
    if method=='paw' and not all(c['paw'] is not None for c in cases): return
    fig,axes=plt.subplots(2,4,figsize=(10.6,5.7),layout='constrained',sharex=True)
    for j,(data,label) in enumerate(zip(cases,LABELS)):
        maximum=0.
        for row,k in enumerate(data['vertices'][[2,4]]):
            x,plus,minus,_=spectrum(data,int(k),method); maximum=max(maximum,float(np.max(plus+minus)))
            ax=axes[row,j]
            ax.plot(x,plus,color=COLORS[0],lw=1.4,label=r'$I_+$')
            ax.plot(x,minus,color=COLORS[1],lw=1.3,ls='--',label=r'$I_-$')
            ax.set(title=f'{label}: '+('K' if row==0 else "K′"),xlim=DISPLAY_PHOTON_RANGE_EV)
            if row==1: ax.set_xlabel('Photon energy (eV)')
            if j==0: ax.set_ylabel('Strength (Å²/eV)' if method=='paw' else 'Intensity (arb. units)')
        for row in range(2): axes[row,j].set_ylim(0,max(maximum*1.05,1e-12))
    axes[0,0].legend(frameon=False,fontsize=9)
    save(fig,out,f'stacking-{method}-spectra')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-dir',type=Path,default=HERE/'reference')
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    if args.output_dir.exists(): parser.error('Choose a new output directory')
    cases=[load_case(args.reference_dir/name) for name in CASES]
    args.output_dir.mkdir(parents=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelsize':11,'axes.titlesize':11,
                         'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42})
    compact(cases,args.output_dir)
    for method in ('native','paw'): strengths(cases,args.output_dir,method)
    write_display_diagnostics(cases,args.output_dir)
    record={'status':'PASS','case_order':list(CASES),'band_origin':'Direct VASP path eigenvalues from matching WAVECAR/EIGENVAL',
            'display_photon_range_eV':list(DISPLAY_PHOTON_RANGE_EV),'stored_photon_range_eV':[.5,4.],
            'energy_zero':'Each case own sampled path VBM; not absolute alignment between slabs',
            'helicity':'Native LEFT/RIGHT and PAW plus/minus share (x +/- i y) final-bra/initial-ket convention, +z beam',
            'eta_display_mask':{'relative_intensity_floor_per_k':.001,'native_absolute_intensity_floor_au':.0001,
                'native_intensity_print_resolution_au':NATIVE_INTENSITY_RESOLUTION,
                'native_eta_absolute_uncertainty_max':NATIVE_ETA_UNCERTAINTY_LIMIT,
                'native_eta_uncertainty_bound':'1e-4/(I_plus+I_minus); each nonnegative printed channel has rounding error <=5e-5',
                'raw_data_policy':'Original intensities and ratios retained; no material-specific zeroing',
                'diagnostics':'native-eta-display.csv'},
            'intensity_comparison':'Native photon/path weighting and PAW length-matrix strengths differ; only eta is overlaid',
            'scope':'Independent-particle k-resolved transitions; not PL, absolute absorption or BZ integration',
            'source_summary_sha256':{name:sha(args.reference_dir/name/'result.json') for name in CASES},
            'figure_sha256':{p.name:sha(p) for p in args.output_dir.iterdir() if p.suffix in {'.png','.pdf','.svg'}},
            'display_diagnostics_sha256':sha(args.output_dir/'native-eta-display.csv')}
    (args.output_dir/'plot.json').write_text(json.dumps(record,indent=2)+'\n')
    print(args.output_dir)


if __name__=='__main__': main()
