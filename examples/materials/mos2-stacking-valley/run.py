#!/usr/bin/env python3
"""Evaluate MoS2 path bands and native circular optics, with optional PAW optics."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/'tools'))
from wavecar_fukui import Wavecar
from prepare_vasp import structure


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''): digest.update(block)
    return digest.hexdigest()


def table(path, header, rows):
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(header); writer.writerows(rows)


def check_eigenval(path, wave):
    lines = path.read_text().splitlines(); nelect, nk, nb = map(int, lines[5].split())
    if (nk, nb) != wave.energies.shape: raise ValueError('EIGENVAL/WAVECAR dimensions differ')
    cursor = 6; energies = []; points = []
    for _ in range(nk):
        while not lines[cursor].strip(): cursor += 1
        points.append([float(v) for v in lines[cursor].split()[:3]]); cursor += 1
        energies.append([float(line.split()[1]) for line in lines[cursor:cursor+nb]]); cursor += nb
    if not np.allclose(points, wave.kpoints, atol=5.1e-7, rtol=0):
        raise ValueError('EIGENVAL/WAVECAR coordinates differ')
    error = float(np.max(abs(np.array(energies)-wave.energies)))
    if error > 5.1e-6: raise ValueError('EIGENVAL/WAVECAR eigenvalues differ')
    return nelect, error


def main():
    material = json.loads((HERE/'inputs/provenance.json').read_text())
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=tuple(material['cases']), required=True)
    parser.add_argument('--run-dir', type=Path, required=True, help='Completed ordinary VASP symmetry-path run')
    parser.add_argument('--binary', type=Path, default=ROOT/'build/vaspberry-gfortran')
    parser.add_argument('--paw-optics', action='store_true', help='Also process matching standard VASP WAVEDER')
    parser.add_argument('--final-band', type=int, help='Complete final boundary; otherwise choose at least 2/4 empty states')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(); source = args.run_dir.resolve(); out = args.output_dir.resolve()
    if out.exists(): parser.error('Choose a new output directory')
    if not args.binary.is_file(): parser.error('Build native VASPBERRY with make serial, or supply --binary')
    if not (source/'POTCAR').is_file() or sha(source/'POTCAR') != material['potcar_sha256']:
        parser.error('The completed reference run must retain its matching documented Mo/S POTCAR')
    text = (source/'OUTCAR').read_text(errors='replace')
    if not ('General timing and accounting informations for this job:' in text
            and 'aborting loop because EDIFF is reached' in text):
        parser.error('Completed electronically converged static VASP output required')
    wave = Wavecar(source/'WAVECAR', spinor_components=2); nocc = material['cases'][args.case]['occupied_spinor_bands_expected']
    expected_structure = structure(HERE/'inputs'/args.case/'scf/POSCAR')
    actual_structure = structure(source/'POSCAR')
    displacement = actual_structure[3]-expected_structure[3]; displacement -= np.rint(displacement)
    if (actual_structure[:2] != expected_structure[:2]
            or not np.allclose(actual_structure[2], expected_structure[2], atol=1e-8, rtol=0)
            or not np.allclose(wave.header.lattice, actual_structure[2], atol=1e-8, rtol=0)
            or np.max(abs(displacement)) > 1e-8):
        parser.error('Wavefunction/structure does not match the named fixed-geometry reference')
    nelect, energy_error = check_eigenval(source/'EIGENVAL', wave)
    if wave.header.ispin != 1 or nelect != nocc:
        parser.error('This reference requires the documented neutral two-component SOC system')
    if not (np.allclose(wave.occupations[:, :nocc], 1., atol=1e-5, rtol=0)
            and np.allclose(wave.occupations[:, nocc:], 0., atol=1e-5, rtol=0)):
        parser.error('The occupied/empty partition is inconsistent with this material example')
    energy = wave.energies; nk, nb = energy.shape; threshold = .002
    if np.min(energy[:, nocc]-energy[:, nocc-1]) <= threshold:
        parser.error('Occupied/empty boundary is unresolved at the 2 meV group threshold')
    minimum_empty = 2 if args.case == 'monolayer' else 4
    candidates = [j for j in range(nocc+minimum_empty, nb)
                  if np.min(energy[:, j]-energy[:, j-1]) > threshold]
    final = args.final_band if args.final_band is not None else next(iter(candidates), None)
    if final not in candidates: parser.error('Final boundary must close every 2 meV group and retain a source guard band')
    if (nk-1) % 4 != 0: parser.error('Expected the documented four-segment Gamma-M-K-Gamma-Kprime path')
    node_indices = np.arange(5)*((nk-1)//4)
    nodes = np.array(material['path']['vertices_fractional']); nsegment = (nk-1)//4+1
    expected_points = np.vstack([np.linspace(a,b,nsegment)[:-1] for a,b in zip(nodes[:-1],nodes[1:])]+[nodes[-1:]])
    if not np.allclose(wave.kpoints, expected_points, atol=1e-8, rtol=0):
        parser.error('Path points differ from the documented linear segments')
    distance = np.r_[0., np.cumsum(np.linalg.norm(np.diff(wave.kpoints, axis=0) @ wave.header.reciprocal, axis=1))]
    source_names = ['WAVECAR','EIGENVAL','POSCAR','POTCAR','INCAR','KPOINTS','OUTCAR']+(['WAVEDER'] if args.paw_optics else [])
    hashes = {name:sha(source/name) for name in source_names}
    record = {'status':'RUNNING', 'case':args.case, 'nocc':nocc, 'source_nbands':nb,
              'nkpts':nk, 'optical_initial_bands':[1,nocc], 'optical_final_bands':[nocc+1,final],
              'optical_final_boundary_gap_eV':float(np.min(energy[:,final]-energy[:,final-1])),
              'degenerate_group_threshold_eV':threshold, 'sampled_vbm_eV':float(energy[:,nocc-1].max()),
              'sampled_path_gap_eV':float(energy[:,nocc].min()-energy[:,nocc-1].max()),
              'path_node_indices_zero_based':node_indices.tolist(), 'path_labels':material['path']['labels'],
              'lattice_A':wave.header.lattice.tolist(), 'source_sha256':hashes,
              'eigenval_energy_max_difference_eV':energy_error, 'commands':[],
              'convergence':'Fixed-geometry path demonstration; no full-BZ gap or response-convergence claim.',
              'native_operator':'Pseudo-wavefunction canonical momentum, native photon/path normalization',
              'optical_channel_convention':'Native LEFT corresponds to +; RIGHT to - for (x +/- i y), propagation +z.'}
    out.mkdir(parents=True)
    save = lambda: (out/'result.json').write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
    def execute(name, command, cwd):
        entry = {'name':name, 'argv':[str(v) for v in command], 'status':'RUNNING'}
        record['commands'].append(entry); save(); started = time.monotonic()
        try:
            with (out/(name+'.stdout.log')).open('w') as stdout, (out/(name+'.stderr.log')).open('w') as stderr:
                process = subprocess.run(entry['argv'], cwd=cwd, stdout=stdout, stderr=stderr, timeout=3600,
                    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1'})
        except subprocess.TimeoutExpired:
            entry.update(status='TIMEOUT',elapsed_seconds=time.monotonic()-started);save();raise
        entry.update(returncode=process.returncode, elapsed_seconds=time.monotonic()-started,
                     status='PASS' if process.returncode == 0 else 'FAILED'); save()
        if process.returncode: raise RuntimeError(f'{name} failed; see retained stderr log')
    try:
        table(out/'bands.csv', ['k_index','kx_fractional','ky_fractional','kz_fractional','path_distance_inv_A','band','energy_eV','occupation'],
              ([k+1,*wave.kpoints[k],distance[k],b+1,energy[k,b],wave.occupations[k,b]] for k in range(nk) for b in range(nb)))
        np.savez_compressed(out/'bands.npz', energies_eV=energy, kpoints_fractional=wave.kpoints,
                            distance_inv_A=distance, lattice_A=wave.header.lattice)
        native = out/'native-optical'; native.mkdir()
        execute('native-optical', [args.binary.resolve(),'-f',source/'WAVECAR','-s','2','-kx',nk,'-ky','1',
            '-cd','2','-if',final,'-ien','.5','-fen','4','-nediv','351','-sigma','.05','-theta','0','-phi','0','-o','optical'], native)
        channels = np.array([[np.loadtxt(native/f'CIRC_DICHROISM_W.optical_{channel}_KP{k+1}.dat')
                              for k in range(nk)] for channel in ('LEFT','RIGHT')])
        if channels.shape != (2,nk,351,8) or not np.isfinite(channels).all():
            raise ValueError('Unexpected native spectrum dimensions or nonfinite output')
        if not np.allclose(channels[:,:,:,5:8],wave.kpoints[None,:,None,:],atol=5.1e-7,rtol=0):
            raise ValueError('Native output k points differ from WAVECAR')
        if not np.allclose(channels[:,:,:,3],np.linspace(.5,4.,351)[None,None,:],atol=5.1e-5,rtol=0):
            raise ValueError('Native output photon grid differs from the requested grid')
        photon = channels[0,0,:,3]; plus,minus = channels[:,:,:,4]; total = plus+minus
        if min(float(plus.min()),float(minus.min())) < 0:
            raise ValueError('Native optical channel intensity cannot be negative')
        floor = np.maximum(1e-4, 1e-3*np.max(total, axis=1, keepdims=True)); valid = total > floor
        eta = np.divide(plus-minus,total,out=np.full_like(total,np.nan),where=valid)
        table(out/'native-optical.csv', ['k_index','kx_fractional','ky_fractional','kz_fractional','photon_energy_eV',
              'left_spectrum_au','right_spectrum_au','eta','valid_eta'],
              ([k+1,*wave.kpoints[k],photon[j],plus[k,j],minus[k,j],eta[k,j] if valid[k,j] else '',int(valid[k,j])]
               for k in range(nk) for j in range(len(photon))))
        np.savez_compressed(out/'native-optical.npz', photon_eV=photon, plus_au=plus, minus_au=minus, eta=eta, eta_valid=valid)
        record['native_eta_mask'] = {'relative_floor_per_k':.001,'absolute_floor_au':.0001}
        record['optical_sigma_eV'] = .05
        if args.paw_optics:
            execute('paw-optical', [sys.executable,ROOT/'tools/vaspberry_kubo.py','waveder-optics',
                '--run-dir',source,'--occupied',nocc,'--spinor-components','2',
                '--energy-reference','unchanged VASP eigenvalue zero','--initial','1',nocc,'--final',nocc+1,final,
                '--photon-min','.5','--photon-max','4','--photon-num','351','--sigma-eV','.05',
                '--relative-intensity-floor','.001','--formats','csv','dat','npz','--output-dir',out/'paw-optical'], ROOT)
        if any(sha(source/name) != digest for name,digest in hashes.items()):
            raise ValueError('A source file changed during analysis')
        record['status'] = 'PASS'
    except Exception as exc:
        record.update(status='FAILED', error=f'{type(exc).__name__}: {exc}'); raise
    finally:
        record['output_sha256'] = {p.name:sha(p) for p in out.iterdir() if p.suffix in {'.csv','.npz'}}; save()
    print(f'Completed {args.case}: {out}')


if __name__ == '__main__': main()
