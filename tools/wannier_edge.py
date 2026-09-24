"""Finite-strip spectra of a supplied orthonormal Wannier Hamiltonian.

The cut removes hoppings leaving the strip. It models an ideal termination,
without edge relaxation, self-consistency or boundary potentials.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
import time

import numpy as np

from berry_data import require
from exported_matrix_kubo import sha256
from wannier_operators import read_operators, validate_operators


def strip_hamiltonian(data, q_parallel, width, periodic_axis=0, open_axis=1):
    """H[(l,n),(l+R_open,m)] += exp(2 pi i q R_parallel) H_R[n,m]."""
    require(type(width) is int and width >= 2, 'strip width must be an integer >=2')
    require(periodic_axis in (0, 1, 2) and open_axis in (0, 1, 2)
            and periodic_axis != open_axis, 'distinct periodic and open lattice axes required')
    require(np.isfinite(q_parallel), 'finite parallel fractional coordinate required')
    n = data.hamiltonian_eV.shape[-1]
    require(width*n <= 4096, 'strip dimension exceeds 4096; reduce width or model size')
    matrix = np.zeros((width, n, width, n), dtype=np.complex128)
    phase = np.exp(2j*np.pi*q_parallel*data.irvec[:, periodic_axis])
    for shift in np.unique(data.irvec[:, open_axis]):
        shift = int(shift)
        if abs(shift) >= width:
            continue
        mask = data.irvec[:, open_axis] == shift
        block = np.einsum('r,rnm->nm', phase[mask], data.hamiltonian_eV[mask])
        rows = np.arange(max(0, -shift), min(width, width-shift))
        matrix[rows, :, rows+shift, :] += block
    matrix = matrix.reshape(width*n, width*n)
    require(np.max(abs(matrix-matrix.conj().T)) <= 1e-8+1e-10*np.max(abs(matrix)),
            'non-Hermitian strip Hamiltonian; no symmetry repair is applied')
    return matrix


def strip_spectrum(data, q_parallel, width, edge_cells=2, periodic_axis=0, open_axis=1,
                   *, memory_limit_mib=2048., time_limit=3600., progress=None):
    """Eigenenergies and total probability in each selected boundary region.

    Individual weights inside an exact degeneracy depend on the eigenvector
    basis; their sum over that degeneracy, or a broadened spectral function,
    is invariant. No artificial spin polarization or edge localization is made.
    """
    validate_operators(data)
    q = np.asarray(q_parallel, dtype=float)
    require(q.ndim == 1 and 1 <= len(q) <= 10001 and np.isfinite(q).all(),
            'one to 10001 finite parallel coordinates required')
    require(type(width) is int and width >= 2 and type(edge_cells) is int
            and 1 <= edge_cells <= width//2, 'nonoverlapping positive edge-cell regions required')
    n = data.hamiltonian_eV.shape[-1]; dimension = n*width
    require(dimension <= 4096, 'strip dimension exceeds 4096')
    estimate = 16*dimension**2*10 + len(q)*dimension*8*3
    require(np.isfinite(memory_limit_mib) and memory_limit_mib > 0
            and estimate <= memory_limit_mib*2**20, 'strip memory estimate exceeds limit')
    require(np.isfinite(time_limit) and time_limit > 0, 'positive finite time limit required')
    energies = np.empty((len(q), dimension)); left = np.empty_like(energies); right = np.empty_like(energies)
    start = time.monotonic()
    for i, point in enumerate(q):
        h = strip_hamiltonian(data, float(point), width, periodic_axis, open_axis)
        energies[i], u = np.linalg.eigh(h)
        left[i] = np.sum(abs(u[:n*edge_cells])**2, axis=0)
        right[i] = np.sum(abs(u[-n*edge_cells:])**2, axis=0)
        require(np.all(left[i]+right[i] <= 1+1e-10), 'invalid edge probabilities')
        if progress is not None:
            progress(i+1, len(q), time.monotonic()-start)
        require(time.monotonic()-start <= time_limit, 'strip time limit reached between k points')
    require(np.isfinite(energies).all(), 'nonfinite strip energies')
    return dict(energies_eV=energies, left_edge_weight=left, right_edge_weight=right,
                q_parallel=q, estimated_memory_MiB=estimate/2**20,
                elapsed_seconds=time.monotonic()-start)


def command(args):
    require(not args.output_dir.exists(), 'output exists; choose a new directory')
    partial = args.output_dir.with_name(args.output_dir.name+'.partial')
    require(not partial.exists(), 'unfinished output exists; choose a new directory')
    require(args.kpoints >= 2 and np.isfinite([args.q_min, args.q_max]).all()
            and args.q_max > args.q_min, 'increasing finite parallel interval and >=2 k points required')
    require(bool(args.formats) and len(set(args.formats)) == len(args.formats)
            and set(args.formats) <= {'csv','dat','npz'}, 'choose distinct CSV, DAT or NPZ formats')
    hashes = {name:sha256(args.operators/name) for name in ('operators.npz','operators.json')}
    data = read_operators(args.operators)
    require(all(sha256(args.operators/name) == value for name,value in hashes.items()),
            'operator cache changed during loading')
    partial.mkdir(parents=True)
    status = dict(status='RUNNING', completed_points=0, total_points=args.kpoints)
    def record():
        temporary = partial/'run.json.tmp'
        temporary.write_text(json.dumps(status, indent=2, allow_nan=False)+'\n')
        temporary.replace(partial/'run.json')
    def progress(done, total, elapsed):
        status.update(completed_points=done, total_points=total, elapsed_seconds=elapsed); record()
    record()
    try:
        result = strip_spectrum(data, np.linspace(args.q_min,args.q_max,args.kpoints), args.width,
            args.edge_cells,args.periodic_axis,args.open_axis,memory_limit_mib=args.memory_limit_mib,
            time_limit=args.time_limit,progress=progress)
        require(all(sha256(args.operators/name) == value for name,value in hashes.items()),
                'operator cache changed during strip calculation')
        arrays = {key:value for key,value in result.items() if isinstance(value,np.ndarray)}
        arrays['lattice_A'] = data.lattice_A
        if 'npz' in args.formats:
            np.savez_compressed(partial/'edge.npz', **arrays)
        names = ['k_index','q_parallel','strip_band_id','energy_eV','left_edge_weight','right_edge_weight']
        for fmt in ('csv','dat'):
            if fmt not in args.formats: continue
            with (partial/('edge.'+fmt)).open('w', newline='') as f:
                writer = csv.writer(f, delimiter=',' if fmt == 'csv' else ' ')
                writer.writerow(names)
                for k, q in enumerate(result['q_parallel']):
                    writer.writerows((k+1,q,b+1,e,result['left_edge_weight'][k,b],result['right_edge_weight'][k,b])
                                     for b,e in enumerate(result['energies_eV'][k]))
        metadata = dict(schema='vaspberry.wannier-edge', version=1, complete=True,
            producer='VASPBERRY real-space finite strip and NumPy eigh', source_metadata=data.metadata,
            source_cache_sha256=hashes, periodic_axis=args.periodic_axis, open_axis=args.open_axis,
            excluded_fractional_coordinate=0., width_cells=args.width, edge_cells=args.edge_cells,
            point_count=args.kpoints, strip_band_count=result['energies_eV'].shape[1],
            units={'energy':'eV','q_parallel':'fractional reciprocal coordinate','edge_weight':'dimensionless probability'},
            interpretation='Ideal cut of the supplied finite Wannier model; no edge reconstruction, self-consistency or additional boundary potential.',
            degeneracy_scope='Individual eigenvector edge weights may rotate within a degeneracy; sum over the degenerate group or use a spectral function.',
            estimated_memory_MiB=result['estimated_memory_MiB'], memory_limit_mib=args.memory_limit_mib,
            memory_limit_scope='preflight estimate, not an OS memory cap',
            elapsed_seconds=result['elapsed_seconds'], time_limit_seconds=args.time_limit,
            time_limit_scope='checked between k-point eigensolves',
            implementation_sha256=sha256(Path(__file__)),
            output_sha256={fmt:sha256(partial/('edge.'+fmt)) for fmt in args.formats})
        (partial/'edge.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
        status['status']='PASS'; record()
        require(not args.output_dir.exists(), 'output appeared during calculation')
        partial.rename(args.output_dir)
        return metadata
    except BaseException as exc:
        status.update(status='FAILED', error=str(exc)); record(); raise


def add_command(sub):
    p = sub.add_parser('wannier-edge',help='ideal finite-strip edge spectrum from a validated Wannier Hamiltonian')
    p.add_argument('--operators',type=Path,required=True)
    p.add_argument('--width',type=int,required=True,help='number of unit cells along the open lattice axis')
    p.add_argument('--edge-cells',type=int,default=2)
    p.add_argument('--periodic-axis',type=int,choices=[0,1,2],default=0)
    p.add_argument('--open-axis',type=int,choices=[0,1,2],default=1)
    p.add_argument('--q-min',type=float,default=-.5);p.add_argument('--q-max',type=float,default=.5)
    p.add_argument('--kpoints',type=int,default=201)
    p.add_argument('--memory-limit-mib',type=float,default=2048.)
    p.add_argument('--time-limit',type=float,default=3600.)
    p.add_argument('--formats',nargs='+',choices=['csv','dat','npz'],default=['npz'])
    p.add_argument('--output-dir',type=Path,required=True)
