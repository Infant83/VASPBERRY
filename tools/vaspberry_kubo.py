#!/usr/bin/env python3
"""Composable VASP operators -> curvature -> two-dimensional response functions.

Use explicit band windows, spin multiplicity, sampling and source normalization.
No material-specific energy reference, electron count or valley is assumed.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

from berry_data import (CurvatureData, base_metadata, hall_spectrum, read_curvature,
                       require, write_curvature, write_hall)
from exported_matrix_kubo import (ARRAYS as MATRIX_ARRAYS, SCHEMA, UNITS, _band_range,
                                 berry_curvature, read_matrix_bundle, sha256)
from wavecar_fukui import Wavecar
from kubo_hall_workflow import COMMANDS as PAIR_COMMANDS, add_commands as add_pair_commands
from waveder_hall import add_command as add_waveder_command, command as waveder_hall_command
from waveder_optics import add_command as add_optics_command, command as waveder_optics_command
from wannier_workflow import COMMANDS as WANNIER_COMMANDS, add_commands as add_wannier_commands
from wannier_bands import add_command as add_wannier_bands, command as wannier_bands_command
from wannier_edge import add_command as add_wannier_edge, command as wannier_edge_command
from spin_hall_workflow import COMMANDS as SPIN_COMMANDS, add_commands as add_spin_commands
from vasp_spin_export import add_arguments as add_spin_export_arguments, command as spin_export_command
from spin_assembly import add_command as add_spin_merge, command as spin_merge_command
from velocity_pairs import add_command as add_velocity_pairs, command as velocity_pairs_command

__version__ = '1.4.0'


def sampling(args):
    if args.mesh is None:
        return {'kind': 'points'}
    return {'kind': 'uniform_full_2d', 'mesh': args.mesh, 'plane_axes': args.plane_axes}


def provenance(args, paths):
    return {'vaspberry_version': __version__, 'command': vars_as_json(args),
            'source_hashes': {str(p): sha256(p) for p in paths},
            'implementation_sha256': {name: sha256(Path(__file__).with_name(name))
                                      for name in ('vaspberry_kubo.py', 'berry_data.py', 'exported_matrix_kubo.py')}}


def vars_as_json(args):
    return {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}


def matrix_command(args):
    data = read_matrix_bundle(args.matrices, args.metadata, allow_experimental=args.allow_experimental)
    require(not (data.metadata.get('source_nspin') == 2 and args.spin_multiplicity != 1),
            'one explicit collinear spin channel requires multiplicity one; sum channels separately')
    result = berry_curvature(data, interest_band_ids=args.n_bands, intermediate_band_ids=args.m_bands,
                             degeneracy_threshold_eV=args.degeneracy_threshold_eV,
                             degeneracy_policy=args.degeneracy_policy, allow_experimental=args.allow_experimental)
    indices = [int(np.flatnonzero(data.band_ids == n)[0]) for n in result.band_ids]
    meta = base_metadata(source_nbands=data.metadata['source_nbands'], spin_multiplicity=args.spin_multiplicity,
                         method='interband_matrix_kubo', energy_reference=args.energy_reference,
                         sampling=sampling(args), source_operator=data.metadata['operator'],
                         provenance=provenance(args, [args.matrices, args.metadata]))
    meta['kernel_diagnostics'] = result.diagnostics
    meta['source_provenance'] = data.metadata['provenance']
    for key in ('source_nspin', 'spin_channel_1based', 'spinor_components'):
        if key in data.metadata:
            meta[key] = data.metadata[key]
    return write_curvature(args.output_dir, CurvatureData(meta, result.k_ids, result.band_ids,
        result.intermediate_band_ids, data.kpoints_fractional, data.weights, data.energies_eV[:, indices],
        result.omega_A2, result.valid_nondegenerate, result.min_gap_eV, data.lattice_A, data.reciprocal_inv_A))


def import_legacy(args):
    """Explicit one-time migration of old or corrected per-band canonical data."""
    w = Wavecar(args.wavecar, spin=args.spin, spinor_components=args.spinor_components)
    require(not (args.spinor_components == 2 and args.spin_multiplicity != 1),
            'SOC spinors must have multiplicity one')
    require(not (w.header.ispin == 2 and args.spin_multiplicity != 1),
            'one explicit collinear spin channel requires multiplicity one; sum channels separately')
    # Resolve WAVECAR scalar/spinor layout instead of inferring it from ISPIN alone.
    w.coefficients(0, [1])
    with args.csv.open() as f:
        lines = f.readlines()
    comments = [line.strip() for line in lines if line.lstrip().startswith('#')]
    marker = '\n'.join(comments)
    require(not (args.normalization == 'legacy-double' and 'STANDARD_MINUS_TWO_IM' in marker),
            'standard-normalized producer cannot be divided by two again')
    column = 'omega_legacy_A2' if args.normalization == 'legacy-double' else 'omega_z_A2'
    nk, nb = w.energies.shape
    z = np.full((nk, nb), np.nan); gaps = np.full((nk, nb), np.nan)
    seen = set()
    for row in csv.DictReader(line for line in lines if not line.lstrip().startswith('#')):
        if 'spin' in row and int(row['spin']) != args.spin:
            continue
        require(column in row, 'normalization requires explicit CSV column '+column)
        k, b = int(row['k_index'])-1, int(row['band'])-1
        require(0 <= k < nk and 0 <= b < nb and (k, b) not in seen, 'duplicate/out-of-range k or band')
        seen.add((k, b))
        delta = np.array([float(row['k'+a+'_frac']) for a in 'xyz'])-w.kpoints[k]
        delta -= np.rint(delta)
        require(np.max(abs(delta)) < 1e-8, 'CSV and WAVECAR k coordinates disagree')
        require(abs(float(row['energy_eV'])-w.energies[k,b]) < 1e-7, 'CSV and WAVECAR energies disagree')
        z[k,b] = float(row[column]); gaps[k,b] = float(row['min_gap_eV'])
    require(len(seen) == nk*nb and np.isfinite(z).all() and np.isfinite(gaps).all() and np.all(gaps >= 0),
            'complete finite per-band CSV required')
    require(nb >= 2 and np.isfinite(w.energies).all(),
            'at least two finite WAVECAR bands are required for a known isolation gap')
    # The source CSV gap is an audit value, not an authority for state
    # isolation. Adjacent sorted energies suffice for each state's nearest
    # gap and avoid a K*B*B energy-difference allocation.
    order = np.argsort(w.energies, axis=1, kind='stable')
    sorted_energies = np.take_along_axis(w.energies, order, axis=1)
    adjacent_gaps = np.diff(sorted_energies, axis=1)
    sorted_gaps = np.empty_like(sorted_energies)
    sorted_gaps[:, 0] = adjacent_gaps[:, 0]
    sorted_gaps[:, -1] = adjacent_gaps[:, -1]
    if nb > 2:
        sorted_gaps[:, 1:-1] = np.minimum(adjacent_gaps[:, :-1], adjacent_gaps[:, 1:])
    actual_gaps = np.empty_like(sorted_gaps)
    np.put_along_axis(actual_gaps, order, sorted_gaps, axis=1)
    gap_discrepancy = abs(gaps-actual_gaps)
    factor = .5 if args.normalization == 'legacy-double' else 1.
    omega = np.full((nk, nb, 3), np.nan)
    require(np.isfinite(args.degeneracy_threshold_eV) and args.degeneracy_threshold_eV >= 0,
            'finite nonnegative degeneracy threshold required')
    valid = actual_gaps > args.degeneracy_threshold_eV
    omega[:,:,2] = np.where(valid, factor*z, np.nan)
    meta = base_metadata(source_nbands=nb, spin_multiplicity=args.spin_multiplicity,
                         method='canonical_momentum_kubo', energy_reference=args.energy_reference,
                         sampling=sampling(args), available_components=(False, False, True),
                         source_operator={'kind': 'canonical_momentum', 'accuracy_status': 'approximation',
                         'missing_terms': ['PAW augmentation and nonlocal/SOC velocity corrections']},
                         provenance=provenance(args, [args.csv, args.wavecar]))
    meta['normalization_migration'] = {'input_convention': args.normalization, 'factor_applied': factor,
                                      'producer_comments': comments, 'output_already_physical': True}
    meta['spin_channel'] = args.spin
    meta['source_nspin'] = int(w.header.ispin)
    meta['spin_channel_1based'] = args.spin
    meta['spinor_components'] = args.spinor_components
    meta['degeneracy_threshold_eV'] = args.degeneracy_threshold_eV
    meta['gap_validation'] = {
        'source': 'minimum over all other WAVECAR band energies; sorted adjacent gaps then restored band order',
        'max_abs_difference_reported_vs_energies_eV': float(gap_discrepancy.max()),
        'reported_gap_discrepancy_count_above_1e_7_eV': int(np.count_nonzero(gap_discrepancy > 1e-7)),
        'validity_uses': 'independently recomputed WAVECAR energy gaps; CSV values are not trusted for validity',
    }
    return write_curvature(args.output_dir, CurvatureData(meta, np.arange(1,nk+1,dtype=np.int64),
        np.arange(1,nb+1,dtype=np.int64), np.arange(1,nb+1,dtype=np.int64), w.kpoints.copy(),
        np.full(nk,1/nk), w.energies.copy(), omega, valid, actual_gaps, w.header.lattice.copy(), w.header.reciprocal.copy()))


def hall_command(args):
    data = read_curvature(args.curvature)
    require(args.mu_num >= 1 and np.isfinite([args.mu_min,args.mu_max]).all() and args.mu_max >= args.mu_min
            and (args.mu_num > 1 or args.mu_min == args.mu_max), 'valid mu range and positive number required')
    mus = np.unique(np.r_[np.linspace(args.mu_min,args.mu_max,args.mu_num),args.mu_reference])
    spec = json.loads(args.regions.read_text()) if args.regions else None
    differences = []
    for value in args.difference:
        bits = value.split(':')
        require(len(bits) == 3, 'difference syntax is NAME:LEFT:RIGHT')
        differences.append(tuple(bits))
    rows, meta = hall_spectrum(data,mus,args.temperatures,mu_reference=args.mu_reference,
        region_spec=spec,differences=differences,allow_partial_bands=args.allow_partial_bands,
        band_resolved=args.band_resolved,mu_chunk=args.mu_chunk)
    paths = [args.curvature/'curvature.npz',args.curvature/'curvature.json']
    if args.regions: paths.append(args.regions)
    meta['provenance'] = provenance(args,paths)
    return write_hall(args.output_dir,rows,meta,formats=args.formats)


def demo(args):
    """Public two-band QWZ fixture; no WAVECAR or proprietary data needed."""
    require(args.mesh >= 4 and np.isfinite(args.mass) and args.mass not in (-2.,0.,2.),
            'demo needs mesh>=4 and mass away from gap closings -2,0,2')
    require(not args.output_dir.exists(),'output directory exists; choose a new directory')
    n=args.mesh; x,y=np.meshgrid(np.arange(n)/n,np.arange(n)/n,indexing='ij')
    q=np.column_stack((x.ravel(),y.ravel(),np.zeros(n*n))); kx,ky=(2*np.pi*q[:,:2]).T
    pauli=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
    d=np.column_stack((np.sin(kx),np.sin(ky),args.mass+np.cos(kx)+np.cos(ky)))
    e,u=np.linalg.eigh(np.einsum('ka,aij->kij',d,pauli))
    deriv=np.zeros((n*n,3,2,2),complex)
    deriv[:,0]=np.cos(kx)[:,None,None]*pauli[0]-np.sin(kx)[:,None,None]*pauli[2]
    deriv[:,1]=np.cos(ky)[:,None,None]*pauli[1]-np.sin(ky)[:,None,None]*pauli[2]
    D=u.conj().swapaxes(-1,-2)[:,None]@deriv@u[:,None]
    arrays=dict(k_ids=np.arange(1,n*n+1,dtype=np.int64),band_ids=np.array([1,2],dtype=np.int64),
        kpoints_fractional=q,weights=np.full(n*n,1/(n*n)),energies_eV=e,D_eVA=D,
        coverage=np.ones(D.shape,bool),lattice_A=np.eye(3),reciprocal_inv_A=2*np.pi*np.eye(3))
    args.output_dir.mkdir(parents=True)
    np.savez_compressed(args.output_dir/'matrix.npz',**arrays)
    meta=dict(schema=SCHEMA,version=1,complete=True,source_nkpoints=n*n,source_nbands=2,units=UNITS,
        matrix_axes=['k','cartesian','bra_band','ket_band'],matrix_element_convention='<n|D_a|m>',
        cartesian_components=['x','y','z'],reciprocal_convention='2pi',weights_convention='sum_one',
        diagonal_status='present',operator={'kind':'analytic_QWZ_hbar_velocity',
        'definition':'H=sin(kx)sx+sin(ky)sy+(mass+cos(kx)+cos(ky))sz; a=1 Angstrom; D=dH/dk',
        'hermitian':True,'accuracy_status':'validated','included_terms':['complete two-band model'],
        'missing_terms':[],'validation_evidence':['Analytic Hamiltonian derivatives; two-level d-vector oracle tests. Model only.']},
        provenance={'run_id':'QWZ-demo','exporter_revision':__version__,
                    'source_hashes':{'vaspberry_kubo.py':sha256(__file__)}},
        model={'mass':args.mass,'mesh':[n,n],'lower_band_C_A_i':1 if -2<args.mass<0 else -1 if 0<args.mass<2 else 0},
        matrix_npz_sha256=sha256(args.output_dir/'matrix.npz'))
    (args.output_dir/'matrix.json').write_text(json.dumps(meta,indent=2)+'\n')
    return meta


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version',action='version',version='VASPBERRY '+__version__)
    sub=p.add_subparsers(dest='command',required=True)
    m=sub.add_parser('matrix',help='interband matrix bundle -> physical point curvature')
    m.add_argument('--matrices',type=Path,required=True); m.add_argument('--metadata',type=Path,required=True)
    m.add_argument('--n-bands',type=_band_range,required=True); m.add_argument('--m-bands',type=_band_range,required=True)
    m.add_argument('--degeneracy-policy',choices=['error','mask'],default='error')
    m.add_argument('--allow-experimental',action='store_true',help='preserve diagnostic operator label; does not validate physics')
    old=sub.add_parser('import-legacy',help='migrate explicitly declared canonical per-band CSV, checking WAVECAR')
    old.add_argument('--csv',type=Path,required=True); old.add_argument('--wavecar',type=Path,required=True)
    old.add_argument('--normalization',choices=['legacy-double','physical'],required=True,
                     help='legacy-double reads omega_legacy_A2 and halves once; physical reads omega_z_A2 unchanged')
    old.add_argument('--spin',type=int,default=1); old.add_argument('--spinor-components',type=int,choices=[1,2],required=True)
    for cmd in (m,old):
        cmd.add_argument('--mesh',nargs=2,type=int,metavar=('NX','NY'),help='declare and validate full uniform 2D mesh; omitted means points only')
        cmd.add_argument('--plane-axes',nargs=2,type=int,default=[0,1],help='ordered reciprocal axes (zero-based) for oriented 2D plane')
        cmd.add_argument('--spin-multiplicity',type=int,choices=[1,2],required=True,help='1 for SOC or one explicit spin channel; 2 only for scalar degenerate spin')
        cmd.add_argument('--energy-reference',required=True,help='description of unchanged input energy zero; no Fermi shift is inferred')
        cmd.add_argument('--degeneracy-threshold-eV',type=float,required=True)
        cmd.add_argument('--output-dir',type=Path,required=True)
    h=sub.add_parser('hall',help='physical curvature -> 2D sheet Hall, regional and band contributions')
    h.add_argument('--curvature',type=Path,required=True,help='directory containing curvature.npz and curvature.json')
    h.add_argument('--mu-min',type=float,required=True); h.add_argument('--mu-max',type=float,required=True)
    h.add_argument('--mu-num',type=int,default=241); h.add_argument('--mu-reference',type=float,required=True)
    h.add_argument('--temperatures',nargs='+',type=float,default=[0.])
    h.add_argument('--regions',type=Path,help='JSON disjoint named periodic circles or global k-ID sets')
    h.add_argument('--difference',action='append',default=[],metavar='NAME:LEFT:RIGHT',help='signed regional difference; no implicit factor one-half')
    h.add_argument('--allow-partial-bands',action='store_true',help='explicit partial contribution; output cannot claim total occupied response')
    h.add_argument('--band-resolved',action='store_true'); h.add_argument('--mu-chunk',type=int,default=64)
    h.add_argument('--output-dir',type=Path,required=True)
    h.add_argument('--formats',nargs='+',choices=['csv','dat','npz'],default=['csv','npz'])
    d=sub.add_parser('demo',help='generate open analytic QWZ interband matrix fixture')
    d.add_argument('--mesh',type=int,default=32); d.add_argument('--mass',type=float,default=-1.)
    d.add_argument('--output-dir',type=Path,required=True)
    add_pair_commands(sub)
    add_waveder_command(sub)
    add_velocity_pairs(sub)
    add_optics_command(sub)
    add_wannier_commands(sub)
    add_wannier_bands(sub)
    add_wannier_edge(sub)
    add_spin_commands(sub)
    add_spin_merge(sub)
    add_spin_export_arguments(sub.add_parser('spin-export',help='validate instrumented VASP PAW spin/full-velocity output and create portable matrices'))
    return p


def main(argv=None):
    p=parser(); args=p.parse_args(argv)
    try:
        require(not args.output_dir.exists(),'output directory exists; choose a new directory')
        meta={**PAIR_COMMANDS,**WANNIER_COMMANDS,**SPIN_COMMANDS,'velocity-pairs':velocity_pairs_command,'spin-merge':spin_merge_command,'spin-export':spin_export_command,'wannier-bands':wannier_bands_command,'wannier-edge':wannier_edge_command,'waveder-hall':waveder_hall_command,'waveder-optics':waveder_optics_command,'matrix':matrix_command,'import-legacy':import_legacy,'hall':hall_command,'demo':demo}[args.command](args)
    except (ValueError, OSError, KeyError, TypeError) as exc:
        p.error(str(exc))
    print(json.dumps({'output':str(args.output_dir),'schema':meta['schema'],'version':meta['version']}))


if __name__ == '__main__':
    main()
