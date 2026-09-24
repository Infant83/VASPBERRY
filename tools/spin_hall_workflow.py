"""Checked VASP spin matrices and insulating 2D spin Hall command workflows.

The transport input is produced from VASP eigenstates, not a user-authored
model. Its metadata records the operator, gauge, units and source association.
Numerical checks do not replace validation of the declared operator producer.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
import re
import time
import zipfile

import numpy as np

from berry_data import require
from exported_matrix_kubo import _band_range, sha256
from spin_hall import occupied_spin_curvature, integrate_sheet, CONDUCTANCE_QUANTUM_S
from spin_operators import (extract_wavecar_spin, import_paw_augmentation,
                            write_spin, PHYSICAL_METRIC_ATOL)
from wavecar_fukui import infer_uniform_grid

MATRIX_ARRAYS = ('energies_eV', 'kpoints_fractional', 'weights', 'lattice_A',
                 'spin_pauli', 'overlap', 'velocity_eVA', 'band_indices')
MATRIX_SCHEMA = 'vaspberry.spin-velocity'


def _json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def _hash(value):
    return isinstance(value, str) and re.fullmatch('[0-9a-f]{64}', value) is not None


def read_matrices(path, metadata_path, *, memory_limit_mib=2048., chunk=8):
    """Read complete Cartesian PAW spin and full velocity in one eigenstate gauge.

    A full velocity must contain diagonal and degenerate-block matrix elements;
    i*(E_n-E_m)*connection with missing blocks filled by zero is insufficient.
    No state normalization, Hermitian repair, or missing-block filling occurs.
    """
    path, metadata_path = Path(path), Path(metadata_path)
    hashes = {'matrices': sha256(path), 'metadata': sha256(metadata_path)}
    m = json.loads(metadata_path.read_text())
    require(m.get('schema') == MATRIX_SCHEMA and type(m.get('version')) is int
            and m['version'] == 1 and m.get('complete') is True,
            'complete version1 spin-velocity producer bundle required')
    require(m.get('data_npz_sha256') == hashes['matrices'], 'spin-velocity archive checksum mismatch')
    require(_hash(m.get('source_wavecar_sha256')), 'matching source WAVECAR checksum required')
    require(m.get('producer_status') == 'PASS' and m.get('operator_accuracy') == 'paw_full_velocity'
            and m.get('physical_paw_validated') is True,
            'audited physical PAW spin and full velocity producer required; pseudo-only input is unsupported')
    require(m.get('spin_basis_id') == 'cartesian' and m.get('velocity_basis') == 'cartesian'
            and m.get('gauge') == 'same_eigenvectors_as_source_wavecar',
            'spin and velocity must share the source eigenvectors and Cartesian axes')
    require(m.get('spinor_components') == 2 and m.get('spin_multiplicity') == 1,
            'two-component SOC spinors with multiplicity one required')
    require(m.get('normalization_applied') is False
            and m.get('full_source_band_coverage') is True
            and m.get('diagonal_available') is True and m.get('degenerate_blocks_available') is True,
            'unchanged complete source-band matrices including diagonal and degenerate blocks required')
    require(m.get('units') == {'spin': 'dimensionless Pauli', 'velocity': 'eV angstrom = hbar*v',
                              'energies': 'eV', 'lattice': 'angstrom', 'overlap': 'dimensionless'},
            'explicit Pauli, hbar*v, energy, lattice and overlap units required')
    require(isinstance(m.get('producer'), str) and m['producer'].strip()
            and isinstance(m.get('operator_scope'), str) and m['operator_scope'].strip(),
            'producer and physical operator scope must be recorded')
    source_hashes = m.get('source_sha256')
    require(isinstance(source_hashes, dict) and bool(source_hashes)
            and all(isinstance(k, str) and k and _hash(v) for k,v in source_hashes.items()),
            'audited producer source checksums required')
    require(source_hashes.get('WAVECAR') == m['source_wavecar_sha256'],
            'source WAVECAR association disagrees inside producer metadata')
    nk, nb = m.get('nkpoints'), m.get('source_nbands')
    require(type(nk) is int and nk > 0 and type(nb) is int and nb > 1,
            'positive k-point and band counts required')
    require(type(chunk) is int and chunk > 0, 'positive k-point chunk required')
    # Includes archive arrays, validation temporaries and one spin-current chunk.
    estimate = 16*nb*nb*(20*nk + 36*min(nk,chunk)) + nk*nb*32
    require(np.isfinite(memory_limit_mib) and memory_limit_mib > 0
            and estimate <= memory_limit_mib*2**20, 'spin Hall memory estimate exceeds limit')
    shapes = dict(energies_eV=(nk,nb), kpoints_fractional=(nk,3), weights=(nk,),
        lattice_A=(3,3), spin_pauli=(nk,3,nb,nb), overlap=(nk,nb,nb),
        velocity_eVA=(nk,3,nb,nb), band_indices=(nb,))
    byte_limit = 7*nk*nb*nb*16 + nk*nb*8 + nk*32 + 72 + nb*8 + 65536
    with zipfile.ZipFile(path) as archive:
        entries = archive.infolist()
        require(len(entries) == len(MATRIX_ARRAYS)
                and {e.filename for e in entries} == {k+'.npy' for k in MATRIX_ARRAYS},
                'spin-velocity archive entries disagree')
        require(sum(e.file_size for e in entries) <= byte_limit, 'spin-velocity archive exceeds declared dimensions')
    with np.load(path, allow_pickle=False) as z:
        arrays = {name:z[name] for name in MATRIX_ARRAYS}
    for name, value in arrays.items():
        dtype = np.int64 if name == 'band_indices' else np.complex128 if name in ('spin_pauli','overlap','velocity_eVA') else np.float64
        require(value.shape == shapes[name] and value.dtype == dtype and np.isfinite(value).all(),
                'invalid spin-velocity array '+name)
    require(np.array_equal(arrays['band_indices'], np.arange(1,nb+1)),
            'all source bands in ascending 1-based order are required')
    require(np.all(np.diff(arrays['energies_eV'], axis=1) >= 0), 'ascending source energies required')
    require(abs(np.linalg.det(arrays['lattice_A'])) > 1e-12, 'nonsingular lattice required')
    defect = float(np.max(abs(arrays['overlap']-np.eye(nb))))
    require(defect <= PHYSICAL_METRIC_ATOL, 'PAW overlap fails orthonormality; no normalization is applied')
    herm = {}
    for name in ('spin_pauli','velocity_eVA','overlap'):
        value = arrays[name]
        herm[name] = float(np.max(abs(value-value.conj().swapaxes(-1,-2))))
        require(herm[name] <= 1e-8+1e-10*float(np.max(abs(value))), name+' is not Hermitian')
    # Pauli operators are contractions of norm-one full-space operators.
    bound = 1.
    for start in range(0,nk,chunk):
        o = arrays['overlap'][start:start+chunk,None]
        s = arrays['spin_pauli'][start:start+chunk]
        for sign in (-1,1):
            bound = min(bound, float(np.linalg.eigvalsh(o+sign*s).min()))
    require(bound >= -PHYSICAL_METRIC_ATOL, 'PAW spin violates overlap-relative Pauli bounds')
    require(hashes == {'matrices':sha256(path),'metadata':sha256(metadata_path)},
            'spin-velocity input changed while reading')
    return arrays, m, dict(source_sha256=hashes, estimated_memory_MiB=estimate/2**20,
        overlap_identity_max_abs=defect, hermiticity_max_abs=herm,
        minimum_overlap_plus_minus_spin_eigenvalue=bound)


def spin_matrix_command(args):
    require((args.augmentation is None) == (args.augmentation_metadata is None),
            'PAW augmentation NPZ and metadata must be supplied together')
    data = extract_wavecar_spin(args.wavecar, spin_basis=args.spin_basis, bands=args.bands)
    if args.augmentation is not None:
        data = import_paw_augmentation(data,args.augmentation,args.augmentation_metadata)
    return write_spin(args.output_dir,data)


def spin_hall_command(args):
    require(not args.output_dir.exists(), 'output exists; choose a new directory')
    partial = args.output_dir.with_name(args.output_dir.name+'.partial')
    require(not partial.exists(), 'unfinished output exists; choose a new directory')
    require(bool(args.formats) and len(set(args.formats)) == len(args.formats)
            and set(args.formats) <= {'csv','dat','npz'}, 'choose distinct CSV, DAT or NPZ formats')
    require(len(args.mesh) == 2 and all(type(n) is int and n >= 2 for n in args.mesh), 'two mesh sizes >=2 required')
    axes = args.plane_axes
    require(len(axes) == 2 and len(set(axes)) == 2 and all(a in (0,1,2) for a in axes),
            'two distinct reciprocal plane axes required')
    require(np.isfinite(args.time_limit) and args.time_limit > 0, 'positive finite time limit required')
    require(np.isfinite(args.gap_threshold_eV) and args.gap_threshold_eV > 0,
            'positive finite gap threshold required')
    start_time = time.monotonic()
    arrays, source_meta, check = read_matrices(args.matrices,args.metadata,
        memory_limit_mib=args.memory_limit_mib,chunk=args.k_chunk)
    nk, source_nb = arrays['energies_eV'].shape
    other = next(i for i in range(3) if i not in axes)
    lattice = arrays['lattice_A']; perpendicular = lattice[other]
    require(all(abs(np.dot(perpendicular,lattice[i])) <=
                1e-10*np.linalg.norm(perpendicular)*np.linalg.norm(lattice[i]) for i in axes),
            'sheet spin Hall requires the remaining real-cell vector perpendicular to the chosen plane')
    infer_uniform_grid(arrays['kpoints_fractional'][:,axes+[other]],*args.mesh)
    require(np.allclose(arrays['weights'],1/nk,atol=1e-12,rtol=0),
            'full uniform mesh and uniform normalized weights required; irreducible/path input is unsupported')
    nb = source_nb if args.source_band_limit is None else args.source_band_limit
    require(type(nb) is int and type(args.occupied) is int and 0 < args.occupied < nb <= source_nb,
            '0 < occupied < retained source bands <= source band count required')
    if nb < source_nb:
        require(np.all(arrays['energies_eV'][:,nb]-arrays['energies_eV'][:,nb-1] > args.gap_threshold_eV),
                'source-band truncation would split a degenerate group')
    e = arrays['energies_eV'][:,:nb]
    vbm,cbm = float(e[:,args.occupied-1].max()),float(e[:,args.occupied].min())
    require(cbm-vbm > args.gap_threshold_eV, 'positive sampled global gap required; metallic spin Hall is unsupported')
    mu = .5*(vbm+cbm) if args.mu_eV is None else args.mu_eV
    require(np.isfinite(mu) and vbm < mu < cbm, 'chemical potential must lie strictly inside sampled gap')
    reciprocal = 2*np.pi*np.linalg.inv(arrays['lattice_A']).T
    area_vector = np.cross(*reciprocal[axes]); area = float(np.linalg.norm(area_vector))
    omega = np.empty((nk,3,3,3)); charge = np.empty((nk,3,3))
    partial.mkdir(parents=True)
    status = dict(status='RUNNING',completed_points=0,total_points=nk)
    def record():
        _json(partial/'run.json.tmp',status); (partial/'run.json.tmp').replace(partial/'run.json')
    record()
    try:
        for first in range(0,nk,args.k_chunk):
            stop = min(nk,first+args.k_chunk); sl = slice(first,stop)
            result = occupied_spin_curvature(e[sl],arrays['velocity_eVA'][sl,:,:nb,:nb],args.occupied,
                input_method=source_meta['producer'], spin_pauli=arrays['spin_pauli'][sl,:,:nb,:nb],
                basis_overlap=arrays['overlap'][sl,:nb,:nb], basis_overlap_tolerance=PHYSICAL_METRIC_ATOL,
                gap_threshold_eV=args.gap_threshold_eV,mu_eV=mu)
            omega[sl] = result['omega_spin_A2']; charge[sl] = result['omega_charge_A2']
            status.update(completed_points=stop,elapsed_seconds=time.monotonic()-start_time); record()
            require(status['elapsed_seconds'] <= args.time_limit, 'spin Hall time limit reached between k chunks')
        integrated = integrate_sheet(omega,arrays['weights'],area)
        sigma = integrated['sigma_hbar_over_e_e2_over_h']
        charge_sigma = -area/(2*np.pi)*np.einsum('k,kij->ij',arrays['weights'],charge)
        require(check['source_sha256'] == {'matrices':sha256(args.matrices),'metadata':sha256(args.metadata)},
                'spin-velocity inputs changed during calculation')
        output = dict(kpoints_fractional=arrays['kpoints_fractional'],weights=arrays['weights'],
            lattice_A=arrays['lattice_A'],omega_spin_A2=omega,omega_charge_A2=charge,
            sigma_hbar_over_e_e2_over_h=sigma,sigma_hbar_over_e_S=sigma*CONDUCTANCE_QUANTUM_S,
            charge_sigma_e2_over_h=charge_sigma,mu_eV=np.array(mu))
        if 'npz' in args.formats: np.savez_compressed(partial/'spin_hall.npz',**output)
        for fmt in ('csv','dat'):
            if fmt not in args.formats: continue
            delimiter = ',' if fmt == 'csv' else ' '
            with (partial/('spin_hall.'+fmt)).open('w',newline='') as f:
                writer = csv.writer(f,delimiter=delimiter)
                writer.writerow(['spin','current','electric_field','sigma_hbar_over_e_e2_over_h','sigma_hbar_over_e_S'])
                writer.writerows(('xyz'[s],'xyz'[i],'xyz'[j],sigma[s,i,j],sigma[s,i,j]*CONDUCTANCE_QUANTUM_S)
                                 for s in range(3) for i in range(3) for j in range(3))
            with (partial/('spin_curvature.'+fmt)).open('w',newline='') as f:
                writer = csv.writer(f,delimiter=delimiter)
                writer.writerow(['k_index','kx_fractional','ky_fractional','kz_fractional','weight','spin','current','electric_field','omega_spin_A2'])
                writer.writerows((k+1,*q,arrays['weights'][k],'xyz'[s],'xyz'[i],'xyz'[j],omega[k,s,i,j])
                    for k,q in enumerate(arrays['kpoints_fractional']) for s in range(3) for i in range(3) for j in range(3))
        metadata = dict(schema='vaspberry.spin-hall',version=1,complete=True,
            method='conventional_intrinsic_spin_Kubo_occupied_bundle',source_metadata=source_meta,
            input_validation=check,source_band_count=source_nb,retained_band_count=nb,occupied=args.occupied,
            gap_threshold_eV=args.gap_threshold_eV,
            retained_cutoff_minimum_gap_eV=(float(np.min(arrays['energies_eV'][:,nb]-arrays['energies_eV'][:,nb-1])) if nb < source_nb else None),
            spin_current_method='finite_band_projected_product',
            spin_current_definition='K={P sigma P,P(hbar*v)P}/2; physical J={s,v}/2',
            limitation='The finite-band current omits P sigma Q D P and P D Q sigma P terms; source-band closure and intermediate-band convergence are required. SOC spin is generally not conserved, so spin Hall response is not required to be quantized.',
            spin_formula='sigma/[(hbar/e)*(e^2/h)]=+A_BZ/(4*pi)*sum_k w_k Omega_spin; Omega_spin=-2 Im sum_occ,empty K_nm D_mn/(E_m-E_n)^2',
            charge_formula='sigma_charge/(e^2/h)=-A_BZ/(2*pi)*sum_k w_k Omega_charge',
            tensor_axis_order=['spin','current','electric_field'],cartesian_components=['x','y','z'],
            units={'curvature':'angstrom^2','spin_sheet_conductivity':'(hbar/e)*(e^2/h)','charge_sheet_conductivity':'e^2/h'},
            dimensionality=2,spin_multiplicity=1,temperature_K=0.,mu_eV=mu,
            valence_max_eV=vbm,conduction_min_eV=cbm,global_gap_eV=cbm-vbm,
            minimum_cross_gap_eV=float(np.min(e[:,args.occupied]-e[:,args.occupied-1])),
            metallic_occupations_supported=False,integer_rounding_applied=False,
            sampling={'kind':'uniform_full_2d','mesh':list(args.mesh),'plane_axes':list(axes)},
            bz_area_inv_A2=area,plane_normal_cartesian=(area_vector/area).tolist(),
            tensor_scope='Cartesian matrix response integrated over the declared 2D plane; in-plane current and electric-field contractions give sheet transport.',
            sigma_hbar_over_e_e2_over_h=sigma.tolist(),charge_sigma_e2_over_h=charge_sigma.tolist(),
            memory_limit_MiB=args.memory_limit_mib,memory_limit_scope='preflight estimate, not an OS cap',
            time_limit_seconds=args.time_limit,time_limit_scope='checked between k chunks',k_chunk=args.k_chunk,
            elapsed_seconds=time.monotonic()-start_time,
            implementation_sha256={name:sha256(Path(__file__).with_name(name)) for name in ('spin_hall.py','spin_hall_workflow.py')},
            output_sha256={p.name:sha256(p) for p in sorted(partial.iterdir()) if p.name != 'run.json'})
        _json(partial/'spin_hall.json',metadata)
        status['status']='PASS';record()
        require(not args.output_dir.exists(), 'output appeared during calculation')
        partial.rename(args.output_dir)
        return metadata
    except BaseException as exc:
        status.update(status='FAILED',error=str(exc));record();raise


def add_commands(sub):
    p = sub.add_parser('spin-matrix',help='WAVECAR spinor coefficients -> raw Pauli matrices, optionally with audited PAW augmentation')
    p.add_argument('--wavecar',type=Path,required=True)
    p.add_argument('--spin-basis',required=True,help='explicit spinor axes from source SAXIS/OUTCAR; WAVECAR does not encode them')
    p.add_argument('--bands',type=_band_range,help='1-based inclusive band window; default all source bands')
    p.add_argument('--augmentation',type=Path);p.add_argument('--augmentation-metadata',type=Path)
    p.add_argument('--output-dir',type=Path,required=True)
    p = sub.add_parser('spin-hall',help='audited full PAW spin/velocity -> conventional intrinsic spin Hall of a gapped 2D occupied bundle')
    p.add_argument('--matrices',type=Path,required=True);p.add_argument('--metadata',type=Path,required=True)
    p.add_argument('--mesh',nargs=2,type=int,required=True,metavar=('NX','NY'))
    p.add_argument('--plane-axes',nargs=2,type=int,default=[0,1])
    p.add_argument('--occupied',type=int,required=True,help='complete valence-bundle band count, starting at band1')
    p.add_argument('--source-band-limit',type=int,help='explicit finite-source-band truncation for closure diagnostics')
    p.add_argument('--mu-eV',type=float,help='inside sampled bulk gap; default gap midpoint (unchanged VASP energy zero)')
    p.add_argument('--gap-threshold-eV',type=float,default=1e-8)
    p.add_argument('--k-chunk',type=int,default=8)
    p.add_argument('--memory-limit-mib',type=float,default=2048.)
    p.add_argument('--time-limit',type=float,default=3600.)
    p.add_argument('--formats',nargs='+',choices=['csv','dat','npz'],default=['csv','npz'])
    p.add_argument('--output-dir',type=Path,required=True)


COMMANDS = {'spin-matrix':spin_matrix_command,'spin-hall':spin_hall_command}
