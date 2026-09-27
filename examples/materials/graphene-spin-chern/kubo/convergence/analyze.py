#!/usr/bin/env python3
"""Analyze native local CSVs and independently check actual WAVECAR spin-sector loops."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
sys.path.insert(0, str(REPO/'tools'))
from wavecar_fukui import Wavecar, common_g_indices


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024**2), b''): h.update(b)
    return h.hexdigest()


def read_csv(path):
    with Path(path).open() as f: return list(csv.DictReader(row for row in f if not row.startswith('#')))


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def split(w, k):
    raw = w.coefficients(k, list(range(1,9))).astype(np.complex128)
    c = raw.reshape(8,-1).T
    gram = c.conj().T@c
    values, vectors = np.linalg.eigh(gram)
    if values[0] < 1e-8: raise ValueError('Selected Gram matrix is singular')
    u = c@(vectors/np.sqrt(values))@vectors.conj().T
    ng = raw.shape[-1]
    spin = u.copy(); spin[ng:] *= -1
    evals, eigvec = np.linalg.eigh(u.conj().T@spin)
    if np.min(abs(evals)) < 1e-6 or np.sum(evals>0)!=4: raise ValueError('Projected spin gap/rank guard failed')
    return [(u@eigvec[:,evals>0]).reshape(2,ng,4), (u@eigvec[:,evals<0]).reshape(2,ng,4)]


def loop_data(w, frames, indices, sector):
    product = 1.+0j; smin = 1.
    for i,j in zip(indices, indices[1:]+indices[:1]):
        li,ri = common_g_indices(w.g_vectors(i), w.g_vectors(j), np.zeros(3,dtype=int))
        left = frames[i][sector][:,li,:].reshape(-1,4)
        right = frames[j][sector][:,ri,:].reshape(-1,4)
        overlap = left.conj().T@right
        singular = np.linalg.svd(overlap,compute_uv=False)
        smin = min(smin,float(singular[-1]))
        sign, _ = np.linalg.slogdet(overlap)
        if smin < 1e-8 or abs(sign)<.5: raise ValueError('Independent local overlap guard failed')
        product *= sign
    return -float(np.angle(product)), smin


def analyze(folder, output):
    manifest = json.loads((folder/'input_manifest.json').read_text())
    if manifest.get('stage')!='local_valley_convergence': raise ValueError('Expected matching preparation manifest')
    for name in ['points.csv','KPOINTS','POSCAR','POTCAR','INCAR','CHGCAR']:
        if sha(folder/name)!=manifest['inputs'][name]: raise ValueError('Prepared input changed: '+name)
    outcar = (folder/'OUTCAR').read_text(errors='replace')
    if 'EDIFF is reached' not in outcar or 'General timing and accounting' not in outcar: raise ValueError('VASP did not converge and finish normally')
    native_path = folder/'SPIN_KUBO.csv'
    headers = native_path.read_text().splitlines()
    for required in ['# selected_bands=1:8', '# integration=NONE_ARBITRARY_K_LIST', '# result_status=PASS']:
        if required not in headers: raise ValueError('Expected native occupied local point output: '+required)
    for name, expected in [('spin_axis_cartesian', [0.,0.,1.]), ('spin_frame_row_1', [1.,0.,0.]), ('spin_frame_row_2', [0.,1.,0.]), ('spin_frame_row_3', [0.,0.,1.])]:
        found = [line.split('=',1)[1] for line in headers if line.startswith('# '+name+'=')]
        if len(found)!=1 or not np.allclose([float(x) for x in found[0].split()], expected, atol=1e-12, rtol=0): raise ValueError('Independent validator requires Cartesian z and the standard VASP spin frame')
    points = read_csv(folder/'points.csv'); native = read_csv(native_path)
    w = Wavecar(folder/'WAVECAR')
    if len(points)!=len(native) or len(points)!=w.header.nkpoints or w.header.nbands!=16: raise ValueError('Point/band mismatch')
    if not np.allclose(w.kpoints, [[float(p[k]) for k in ['k1','k2','k3']] for p in points], atol=2e-15, rtol=0): raise ValueError('WAVECAR coordinates do not match preparation')
    frames = [split(w,i) for i in range(len(points))]
    rows = []
    for i,(point, row) in enumerate(zip(points,native)):
        if int(row['k_index'])!=i+1: raise ValueError('Native k order changed')
        for j in range(3):
            if abs(float(row[f'k{j+1}'])-w.kpoints[i,j])>2e-15: raise ValueError('Native/WAVECAR mismatch')
        rows.append({**point,'gap_eV':float(w.energies[i,8]-w.energies[i,7]),
                     'omega_plus_xy_A2':float(row['omega_plus_xy_A2']), 'omega_minus_xy_A2':float(row['omega_minus_xy_A2']),
                     'omega_spin_xy_A2':(float(row['omega_plus_xy_A2'])-float(row['omega_minus_xy_A2']))/2,
                     'spin_mixing_plus_xy_A2':float(row['spin_mixing_plus_xy_A2'])})
    loops = []
    for valley in ['K','Kprime']:
        center = next(i for i,p in enumerate(points) if p['valley']==valley and float(p['radius_Ainv'])==0)
        for radius in manifest['radii_Ainv']:
            ring = [i for i,p in enumerate(points) if p['valley']==valley and float(p['radius_Ainv'])==radius]
            for stride in [2,1]:
                selected = ring[::stride]
                result = dict(valley=valley,radius_Ainv=radius,vertices=len(selected))
                triangles = [[center,a,b] for a,b in zip(selected,selected[1:]+selected[:1])]
                for s,name in enumerate(['plus','minus']):
                    boundary, quality = loop_data(w,frames,selected,s)
                    fans = [loop_data(w,frames,t,s) for t in triangles]
                    if any(abs(f[0]) >= np.pi/2 for f in fans): raise ValueError('A triangle phase reaches pi/2; refine the angular mesh before interpreting the local flux')
                    flux = sum(f[0] for f in fans)
                    phase_error = abs(np.angle(np.exp(1j*(flux-boundary))))
                    if phase_error > 1e-10: raise ValueError('Triangulation/boundary inconsistency')
                    result[f'flux_{name}_rad']=flux
                    result[f'C_local_{name}']=flux/(2*np.pi)
                    result[f'min_link_singular_{name}']=min([quality]+[f[1] for f in fans])
                    result[f'max_triangle_flux_{name}_rad']=max(abs(f[0]) for f in fans)
                result['C_local_spin']=(result['C_local_plus']-result['C_local_minus'])/2
                result['polygon_area_Ainv2']=len(selected)*radius**2*np.sin(2*np.pi/len(selected))/2
                loops.append(result)
    output.mkdir(parents=True, exist_ok=False)
    write_csv(output/'radial.csv',rows); write_csv(output/'loops.csv',loops)
    summary=dict(scope='Actual local points: native canonical-momentum curvature plus independent Gram-normalized projected-spin WAVECAR overlap loops; neither is a full BZ integral',
                 ediff_eV=manifest['ediff_eV'], nelmin=manifest['nelmin'], points=len(rows), angles=manifest['angles'],
                 radius_Ainv=manifest['radii_Ainv'],
                 native_producer='VASPBERRY --task spin-kubo --bands 1:8; source approximation headers retained separately',
                 loop_producer='Independent Python validation from actual WAVECAR; no Hamiltonian model and no native spin-chern task claimed',
                 source_sha256={n:sha(folder/n) for n in ['WAVECAR','OUTCAR','SPIN_KUBO.csv','input_manifest.json','points.csv']},
                 analysis_sha256=sha(Path(__file__)),
                 triangle_phase_guard_rad=float(np.pi/2),
                 sign_branch_scope='Sampled triangle phases must stay below pi/2; this guard does not certify unsampled spatial resolution',
                 minimum_loop_singular=min(min(row['min_link_singular_plus'],row['min_link_singular_minus']) for row in loops),
                 K_gap_eV=rows[0]['gap_eV'], K_omega_plus_A2=rows[0]['omega_plus_xy_A2'])
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,required=True,help='Completed ordinary VASP/native directory containing points.csv')
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    analyze(args.input_dir,args.output_dir)


if __name__=='__main__': main()
