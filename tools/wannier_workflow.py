"""General VASPBERRY workflows for complete Wannier Hamiltonian/position data."""
from __future__ import annotations

from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
import itertools
import json
import os
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np

from berry_data import require, write_hall
from exported_matrix_kubo import sha256
from kubo_pairs import hall_metadata, result_rows
from wannier_hall import occupied_curvature, prepare_fourier
from wavecar_fukui import gauss_reduce_2d_lattice_basis, closest_2d_lattice_residual


def lattice_from_poscar(path):
    lines = Path(path).read_text().splitlines()
    require(len(lines) >= 5, 'truncated POSCAR')
    scale = np.array([float(v) for v in lines[1].split()])
    lattice = np.array([[float(v) for v in line.split()[:3]] for line in lines[2:5]])
    require(lattice.shape == (3, 3) and np.isfinite(lattice).all()
            and abs(np.linalg.det(lattice)) > 1e-12, 'finite nonsingular POSCAR lattice required')
    if len(scale) == 1:
        require(np.isfinite(scale[0]) and scale[0] != 0, 'nonzero POSCAR scale required')
        lattice *= scale[0] if scale[0] > 0 else (-scale[0]/abs(np.linalg.det(lattice)))**(1/3)
    else:
        require(len(scale) == 3 and np.isfinite(scale).all() and np.all(scale > 0),
                'POSCAR requires one nonzero or three positive scales')
        lattice *= scale[None]
    return lattice


def import_command(args):
    from wannier_operators import load_effective_operators, write_operators
    data = load_effective_operators(args.hh, args.aa, lattice_from_poscar(args.poscar),
        spinor_components=args.spinor_components, spin_multiplicity=args.spin_multiplicity,
        energy_reference=args.energy_reference)
    data.metadata['lattice_source'] = {'format':'VASP POSCAR', 'filename':args.poscar.name,
                                     'sha256':sha256(args.poscar)}
    return write_operators(args.output_dir, data)


def quadrature(lattice, mesh, axes=(0, 1), refine=1, radius=None, centers=()):
    """Full periodic 2D cell partition, optionally subdividing selected cells.

    Weights are cell areas divided by total BZ area, never a normalization of
    an arbitrary point cloud. Refinement selects parent-cell centers; children
    replace their parent, including across periodic edges.
    """
    require(len(mesh) == 2 and all(type(v) is int and v >= 2 for v in mesh), 'mesh needs two integers >=2')
    require(len(axes) == 2 and len(set(axes)) == 2 and set(axes) <= {0, 1, 2}, 'two distinct reciprocal axes required')
    require(type(refine) is int and refine >= 1 and refine % 2 == 1, 'positive odd local subdivision required')
    require(refine > 1 or (radius is None and not len(centers)),
            'refinement radius/centers require subdivision greater than one')
    require(refine == 1 or (radius is not None and np.isfinite(radius) and radius > 0 and len(centers)),
            'local refinement needs a positive radius and at least one fractional in-plane center')
    center_array = np.asarray(centers, dtype=float)
    require(not len(centers) or (center_array.shape == (len(centers), 2) and np.isfinite(center_array).all()),
            'refinement centers need two finite fractional coordinates')
    reciprocal = 2*np.pi*np.linalg.inv(lattice).T
    basis = gauss_reduce_2d_lattice_basis(reciprocal[list(axes)])
    points, weights, parents, keys = [], [], [], []
    refined = 0
    require(mesh[0]*mesh[1]*refine**2 <= 100_000_000, 'requested quadrature bound is too large')
    for iy in range(mesh[1]):
        for ix in range(mesh[0]):
            parent = np.array([ix/mesh[0], iy/mesh[1]])
            selected = refine > 1 and any(np.linalg.norm(closest_2d_lattice_residual(
                (parent-center) @ reciprocal[list(axes)], basis)) <= radius for center in center_array)
            child = range(-(refine//2), refine//2+1) if selected else [0]
            mult = refine**2 if selected else 1
            for dy, dx in itertools.product(child, repeat=2):
                key = ((ix*refine+dx) % (mesh[0]*refine), (iy*refine+dy) % (mesh[1]*refine))
                q = np.zeros(3); q[list(axes)] = np.array(key)/(np.array(mesh)*refine)
                points.append(q); keys.append(key); parents.append(iy*mesh[0]+ix)
                weights.append(1/(mesh[0]*mesh[1]*mult))
            refined += int(selected)
            require(len(points) <= 2_000_000, 'quadrature exceeds two million points; split the scientific study')
    q, w, parent = np.asarray(points), np.asarray(weights), np.asarray(parents)
    require(len(set(keys)) == len(keys) and abs(w.sum()-1) < 2e-13, 'invalid periodic quadrature')
    require(np.allclose(np.bincount(parent, weights=w), 1/np.prod(mesh), rtol=0, atol=2e-14),
            'parent-cell weights do not cover the full BZ')
    return q, w, parent, reciprocal, refined


def save_json(path, value):
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def hall_command(args):
    from wannier_operators import read_operators
    require(not args.output_dir.exists(), 'output exists; choose a new directory')
    require(1 <= args.workers <= 8 and 1 <= args.batch_size <= 128, 'workers 1:8 and batch size 1:128 required')
    require(np.isfinite(args.time_limit) and args.time_limit > 0, 'positive finite time limit required')
    require(np.isfinite(args.gap_threshold) and args.gap_threshold > 0, 'positive finite gap threshold required')
    require(bool(args.formats) and len(set(args.formats)) == len(args.formats)
            and set(args.formats) <= {'csv','dat','npz'}, 'choose distinct CSV, DAT or NPZ formats')
    require(args.mu_num >= 1 and np.isfinite([args.mu_min, args.mu_max, args.mu_reference]).all()
            and args.mu_max >= args.mu_min and (args.mu_num > 1 or args.mu_min == args.mu_max), 'invalid mu scan')
    source_hashes = {name:sha256(args.operators/name) for name in ('operators.npz','operators.json')}
    data = read_operators(args.operators)
    require(all(sha256(args.operators/name) == digest for name,digest in source_hashes.items()),
            'operator cache changed while loading')
    require(0 < args.occupied < data.hamiltonian_eV.shape[-1], 'occupied count outside model')
    memory_limit = getattr(args, 'memory_limit_mib', 4096.)
    require(np.isfinite(memory_limit) and memory_limit > 0, 'positive finite memory estimate limit required')
    nb = data.hamiltonian_eV.shape[-1]
    estimated_bytes = (2*(data.hamiltonian_eV.nbytes+data.position_A.nbytes)
                       + 10*data.hamiltonian_eV.nbytes
                       + args.workers*args.batch_size*nb*nb*16*48)
    require(estimated_bytes <= memory_limit*2**20,
            'estimated working memory exceeds limit; reduce workers/batch size or raise the explicit limit')
    q, w, parent, reciprocal, refined = quadrature(data.lattice_A, args.mesh, args.plane_axes,
        args.refine, args.refine_radius, args.refine_center)
    excluded = next(i for i in range(3) if i not in args.plane_axes)
    perpendicular = data.lattice_A[excluded]
    require(all(abs(np.dot(perpendicular, data.lattice_A[i])) <=
                1e-10*np.linalg.norm(perpendicular)*np.linalg.norm(data.lattice_A[i])
                for i in args.plane_axes),
            'sheet Hall requires the remaining real-cell vector perpendicular to the chosen plane')
    mus = np.unique(np.linspace(args.mu_min, args.mu_max, args.mu_num))
    requested = np.r_[mus, args.mu_reference]
    partial = args.output_dir.with_name(args.output_dir.name+'.partial')
    require(not partial.exists(), 'unfinished calculation exists; choose a new output directory')
    partial.mkdir(parents=True)
    start = time.monotonic(); completed = 0
    progress = dict(status='RUNNING', total_points=len(q), completed_points=0)
    save_json(partial/'run.json', progress)
    pool = None
    try:
        fields = prepare_fourier(data)
        terms = np.empty((len(q), 3, 3)); edges = np.empty((len(q), 2))
        max_defect = 0.; minimum_gap = float('inf')
        ranges = iter((lo, min(lo+args.batch_size, len(q))) for lo in range(0, len(q), args.batch_size))
        pool = ThreadPoolExecutor(max_workers=args.workers)
        pending = {}
        def submit():
            bounds = next(ranges, None)
            if bounds is not None:
                lo, hi = bounds
                future = pool.submit(occupied_curvature, data, q[lo:hi], args.occupied,
                                     args.gap_threshold, fields)
                pending[future] = bounds
        for _ in range(args.workers): submit()
        while pending:
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                lo, hi = pending.pop(future); result = future.result()
                terms[lo:hi] = result['omega_terms_A2']
                edges[lo:hi] = result['energies_eV'][:, args.occupied-1:args.occupied+1]
                minimum_gap = min(minimum_gap, result['minimum_gap_eV'])
                max_defect = max(max_defect, result['max_hermiticity_defect'])
                completed += hi-lo
                progress.update(completed_points=completed, elapsed_seconds=time.monotonic()-start)
                save_json(partial/'run.json', progress)
                require(time.monotonic()-start <= args.time_limit, 'wall-time limit reached between batches')
                submit()
        pool.shutdown(); pool = None
        vbm, cbm = float(edges[:, 0].max()), float(edges[:, 1].min())
        require(cbm > vbm and np.all((requested > vbm) & (requested < cbm)),
                'every requested mu must lie strictly inside the sampled global model gap')
        vector = np.cross(*reciprocal[args.plane_axes]); area = float(np.linalg.norm(vector)); normal = vector/area
        mult = data.metadata['spin_multiplicity']
        sigma_terms = -area/(2*np.pi)*mult*np.einsum('k,kja->ja', w, terms)
        sigma = float(sigma_terms.sum(axis=0) @ normal)
        all_mus = np.unique(requested); values = np.zeros((1, 1, len(all_mus), 4))
        values[0, 0, :, 0] = sigma; values[0, 0, :, 1] = args.occupied*mult
        regions = {'total': np.ones(len(q), dtype=bool)}
        metadata = dict(data.metadata, sampling={'kind':'full_2d_cell_partition', 'mesh':args.mesh,
            'plane_axes':args.plane_axes, 'refine':args.refine, 'refine_radius_inv_A':args.refine_radius,
            'refine_centers_fractional':args.refine_center, 'refined_parent_count':refined})
        proxy = SimpleNamespace(metadata=metadata)
        meta = hall_metadata(proxy, area, normal, None, (), regions, args.mu_reference)
        meta.update(method='vaspberry_wannier_full_connection_T0', scope='fixed_occupied_model_bundle',
            producer='VASPBERRY NumPy kernel; no postw90 response input or process',
            formula='Omega=J0+J1+J2; sigma/(e^2/h)=-A_BZ/(2*pi)*g_s*sum_k w_k Omega_normal',
            occupied_model_bands=[1,args.occupied], valence_max_eV=vbm, conduction_min_eV=cbm,
            global_gap_eV=cbm-vbm, minimum_cross_gap_eV=minimum_gap,
            max_hermiticity_defect=max_defect, gap_threshold_eV=args.gap_threshold,
            sigma_terms_e2_over_h=sigma_terms.tolist(), term_order=['J0','J1','J2'],
            axial_component_order=['yz','zx','xy'], weight_sum=float(w.sum()), point_count=len(q),
            integer_rounding_applied=False, electron_count_scope='represented occupied model bands only',
            model_scope='Full response of supplied finite H/position model; source/basis/model convergence is separate.',
            algorithm={'batch_size':args.batch_size,'workers':args.workers,'time_limit_seconds':args.time_limit,
                       'estimated_working_memory_MiB':estimated_bytes/2**20,
                       'memory_estimate_limit_MiB':memory_limit, 'memory_limit_scope':'preflight estimate, not an operating-system RSS cap',
                       'time_limit_scope':'checked between batches; running batches finish on failure',
                       'thread_environment':{k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')}},
            implementation_sha256={name:sha256(Path(__file__).with_name(name))
                for name in ('wannier_hall.py','wannier_operators.py','wannier_workflow.py')},
            source_cache_sha256=source_hashes,
            spatial_scope='Declared two-dimensional periodic model; zero fractional coordinate on the excluded axis, not bulk 3D conductivity.',
            references=['https://doi.org/10.1103/PhysRevB.74.195118', 'https://doi.org/10.1103/PhysRevB.85.014435'])
        np.savez_compressed(partial/'curvature.npz', kpoints_fractional=q, weights=w,
            parent_cell=parent, omega_terms_A2=terms, omega_A2=terms.sum(axis=1),
            band_edges_eV=edges, lattice_A=data.lattice_A, reciprocal_inv_A=reciprocal)
        meta['curvature_npz_sha256'] = sha256(partial/'curvature.npz')
        written = write_hall(partial/'hall', result_rows(values,mus,[0.],all_mus,regions,args.mu_reference),
                             meta, formats=args.formats)
        progress.update(status='PASS', elapsed_seconds=time.monotonic()-start, sigma_e2_over_h=sigma,
                        source_operator_npz_sha256=source_hashes['operators.npz'], source_cache_sha256=source_hashes)
        require(all(sha256(args.operators/name) == digest for name,digest in source_hashes.items()),
                'operator cache changed during calculation')
        save_json(partial/'run.json', progress)
        require(not args.output_dir.exists(), 'output appeared during calculation')
        partial.rename(args.output_dir)
        return written
    except BaseException as exc:
        if pool is not None: pool.shutdown(wait=True, cancel_futures=True)
        progress.update(status='FAILED', elapsed_seconds=time.monotonic()-start, error=str(exc))
        save_json(partial/'run.json', progress)
        raise


def add_commands(sub):
    imp = sub.add_parser('wannier-import', help='full effective-model HH_R/AA_R and VASP lattice -> operator cache')
    imp.add_argument('--hh', type=Path, required=True); imp.add_argument('--aa', type=Path, required=True)
    imp.add_argument('--poscar', type=Path, required=True)
    imp.add_argument('--spinor-components', type=int, choices=[1,2], required=True)
    imp.add_argument('--spin-multiplicity', type=int, choices=[1,2], required=True)
    imp.add_argument('--energy-reference', required=True); imp.add_argument('--output-dir', type=Path, required=True)
    hall = sub.add_parser('wannier-hall', help='VASPBERRY full-connection Wannier occupied-bundle T=0 sheet Hall')
    hall.add_argument('--operators', type=Path, required=True); hall.add_argument('--occupied', type=int, required=True)
    hall.add_argument('--mesh', nargs=2, type=int, required=True); hall.add_argument('--plane-axes', nargs=2, type=int, default=[0,1])
    hall.add_argument('--refine', type=int, default=1); hall.add_argument('--refine-radius', type=float)
    hall.add_argument('--refine-center', nargs=2, type=float, action='append', default=[])
    hall.add_argument('--mu-min', type=float, required=True); hall.add_argument('--mu-max', type=float, required=True)
    hall.add_argument('--mu-num', type=int, default=3); hall.add_argument('--mu-reference', type=float, required=True)
    hall.add_argument('--gap-threshold', type=float, default=1e-8)
    hall.add_argument('--workers', type=int, default=1); hall.add_argument('--batch-size', type=int, default=8)
    hall.add_argument('--time-limit', type=float, default=3600.)
    hall.add_argument('--memory-limit-mib', type=float, default=4096., help='upper bound on the preflight working-memory estimate')
    hall.add_argument('--formats', nargs='+', choices=['csv','dat','npz'], default=['csv','dat','npz'])
    hall.add_argument('--output-dir', type=Path, required=True)


COMMANDS = {'wannier-import':import_command, 'wannier-hall':hall_command}
