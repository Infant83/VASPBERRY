"""Band dispersion from the same validated real-space Hamiltonian as Hall."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import time

import numpy as np

from berry_data import require
from exported_matrix_kubo import sha256
from wannier_operators import read_operators


def path_grid(vertices, points_per_segment, reciprocal):
    vertices = np.asarray(vertices, dtype=float)
    require(vertices.ndim == 2 and vertices.shape[1] == 3 and len(vertices) >= 2
            and np.isfinite(vertices).all(), 'at least two finite three-dimensional path vertices required')
    require(type(points_per_segment) is int and 2 <= points_per_segment <= 10001,
            'points per segment must be an integer in 2:10001')
    require((len(vertices)-1)*(points_per_segment-1)+1 <= 2_000_000, 'path exceeds two million points')
    lengths = np.linalg.norm(np.diff(vertices, axis=0) @ reciprocal, axis=1)
    require(np.all(lengths > 1e-12), 'successive path vertices must differ')
    q = np.concatenate([np.linspace(a,b,points_per_segment)[:-1]
                        for a,b in zip(vertices[:-1],vertices[1:])] + [vertices[-1:]])
    distance = np.r_[0.,np.cumsum(np.linalg.norm(np.diff(q,axis=0) @ reciprocal,axis=1))]
    return q, distance, np.r_[0.,np.cumsum(lengths)]


def command(args):
    require(not args.output_dir.exists(), 'output exists; choose a new directory')
    require(1 <= args.batch_size <= 128, 'batch size must lie in 1:128')
    require(len(args.vertices) % 3 == 0, 'vertices require triples of fractional coordinates')
    source_hashes = {name:sha256(args.operators/name) for name in ('operators.npz','operators.json')}
    data = read_operators(args.operators)
    require(all(sha256(args.operators/name) == digest for name,digest in source_hashes.items()),
            'operator cache changed while loading')
    reciprocal = 2*np.pi*np.linalg.inv(data.lattice_A).T
    vertices = np.asarray(args.vertices).reshape(-1,3)
    q, distance, ticks = path_grid(vertices,args.points_per_segment,reciprocal)
    require(not args.labels or len(args.labels) == len(vertices), 'one label per vertex required')
    require(bool(args.formats) and len(set(args.formats)) == len(args.formats)
            and set(args.formats) <= {'csv','npz'}, 'choose distinct CSV or NPZ formats')
    energy = np.empty((len(q),data.hamiltonian_eV.shape[-1])); start=time.monotonic()
    for lo in range(0,len(q),args.batch_size):
        phase = np.exp(2j*np.pi*q[lo:lo+args.batch_size] @ data.irvec.T)
        h = (phase @ data.hamiltonian_eV.reshape(len(data.irvec),-1)).reshape(-1,energy.shape[1],energy.shape[1])
        require(np.max(abs(h-h.conj().swapaxes(-1,-2))) <= 1e-8+1e-10*np.max(abs(h)),
                'non-Hermitian Fourier Hamiltonian')
        energy[lo:lo+len(h)] = np.linalg.eigvalsh(h)
    require(np.isfinite(energy).all(), 'nonfinite band energies')
    require(all(sha256(args.operators/name) == digest for name,digest in source_hashes.items()),
            'operator cache changed during band calculation')
    args.output_dir.mkdir(parents=True)
    arrays=dict(kpoints_fractional=q,distance_inv_A=distance,energies_eV=energy,
                lattice_A=data.lattice_A,vertices_fractional=vertices,tick_distance_inv_A=ticks)
    if 'npz' in args.formats: np.savez_compressed(args.output_dir/'bands.npz',**arrays)
    if 'csv' in args.formats:
        with (args.output_dir/'bands.csv').open('w',newline='') as f:
            writer=csv.writer(f);writer.writerow(['k_index','q1','q2','q3','distance_inv_A','band_id','energy_eV'])
            for k in range(len(q)):
                writer.writerows((k+1,*q[k],distance[k],b+1,e) for b,e in enumerate(energy[k]))
    meta=dict(schema='vaspberry.wannier-bands',version=1,complete=True,method='VASPBERRY Hamiltonian Fourier sum and NumPy eigvalsh',
              source_metadata=data.metadata,source_cache_sha256=source_hashes,units={'energy':'eV','distance':'Angstrom^-1'},
              labels=args.labels,output_formats=args.formats,point_count=len(q),band_count=energy.shape[1],
              elapsed_seconds=time.monotonic()-start,implementation_sha256=sha256(Path(__file__)),
              output_sha256={fmt:sha256(args.output_dir/('bands.'+fmt)) for fmt in args.formats})
    (args.output_dir/'bands.json').write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n')
    return meta


def add_command(sub):
    p=sub.add_parser('wannier-bands',help='VASPBERRY dispersion along a fractional path from the same Wannier operators')
    p.add_argument('--operators',type=Path,required=True)
    p.add_argument('--vertices',type=float,nargs='+',required=True)
    p.add_argument('--labels',nargs='*',default=[])
    p.add_argument('--points-per-segment',type=int,default=201)
    p.add_argument('--batch-size',type=int,default=32)
    p.add_argument('--formats',nargs='+',choices=['csv','npz'],default=['npz'])
    p.add_argument('--output-dir',type=Path,required=True)
