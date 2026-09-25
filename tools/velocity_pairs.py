"""Audited PAW full velocities -> reusable occupation-weighted charge pairs.

This adapter retains undivided numerators, so filled internal degeneracies
cancel in pair-hall before any energy denominator is evaluated. It currently
accepts the complete spin-velocity bundle from the supported VASP producer;
spin is validated by that reader but is not used for charge Hall response.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from berry_data import require
from exported_matrix_kubo import sha256
from kubo_pairs import PairData, SCHEMA, validate_pairs, write_pairs
from spin_hall_workflow import read_matrices


def velocity_pair_data(arrays, source_metadata, *, sampling, energy_reference):
    """Convert an already validated complete PAW bundle without band truncation.

    Use ``command`` for files: it first applies the full existing producer,
    physical-metric, same-gauge, checksum and operator-scope checks.
    """
    source = source_metadata
    require(source.get('operator_accuracy') == 'paw_full_velocity'
            and source.get('physical_paw_validated') is True
            and source.get('producer_status') == 'PASS'
            and source.get('full_source_band_coverage') is True
            and source.get('diagonal_available') is True
            and source.get('degenerate_blocks_available') is True
            and source.get('normalization_applied') is False,
            'unchanged audited complete physical PAW velocity required')
    require(isinstance(energy_reference, str) and energy_reference.strip(),
            'explicit unchanged energy-reference description required')
    e = arrays['energies_eV']; d = arrays['velocity_eVA']
    require(e.ndim == 2, 'two-dimensional energies required')
    nk, nb = e.shape
    require(d.shape == (nk, 3, nb, nb) and d.dtype == np.complex128
            and np.isfinite(d).all(), 'complete Cartesian complex128 velocities required')
    require(float(np.max(abs(d-d.conj().swapaxes(-1, -2))))
            <= 1e-8 + 1e-10*float(np.max(abs(d))),
            'velocity must be Hermitian; no repair is applied')
    require(source.get('source_nbands') == nb and source.get('nkpoints') == nk,
            'source dimensions disagree')
    n, m = np.triu_indices(nb, 1)
    values = np.stack([-2*np.imag(d[:, a, n, m]*d[:, b, m, n])
                       for a, b in ((1, 2), (2, 0), (0, 1))], axis=-1)
    lattice = arrays['lattice_A'].copy()
    meta = dict(schema=SCHEMA, version=1, complete=True, source_nbands=nb,
        normalization='numerator=-2Im(D_a_nm*D_b_mn)', components=['yz', 'zx', 'xy'],
        units={'energy':'eV', 'numerator':'eV^2 Angstrom^2',
               'lattice':'Angstrom', 'reciprocal':'1/Angstrom'},
        sampling=sampling, spin_multiplicity=1, spin_channel_1based=1,
        source_nspin=1, spinor_components=2, energy_reference=energy_reference,
        source_operator={'kind':'paw_full_velocity',
            'accuracy_status':'audited_within_declared_producer_scope',
            'definition':'D=hbar*v; full PAW, nonlocal and SOC velocity in supported producer',
            'producer':source['producer'], 'scope':source['operator_scope'],
            'normalization_applied':False, 'hermitian_repair_applied':False},
        source_wavecar_sha256=source['source_wavecar_sha256'],
        source_sha256=source['source_sha256'],
        source_producer_metadata=source,
        interpretation='Charge pairs use only velocity. No spin current or band-by-band nondegeneracy is assumed.',
        convergence_status='not_established_by_operator_validation')
    return validate_pairs(PairData(meta, np.arange(1, nk+1), arrays['band_indices'].copy(), n, m,
        arrays['kpoints_fractional'].copy(), arrays['weights'].copy(), e.copy(), values,
        lattice, 2*np.pi*np.linalg.inv(lattice).T))


def command(args):
    require(not args.output_dir.exists(), 'pair output directory exists')
    arrays, meta, check = read_matrices(args.matrices, args.metadata,
                                      memory_limit_mib=args.memory_limit_mib)
    data = velocity_pair_data(arrays, meta,
        sampling={'kind':'uniform_full_2d', 'mesh':args.mesh, 'plane_axes':args.plane_axes},
        energy_reference=args.energy_reference)
    data.metadata['adapter_validation'] = check
    data.metadata['provenance'] = {
        'command':{k:str(v) if isinstance(v, Path) else v for k,v in vars(args).items()},
        'implementation_sha256':{name:sha256(Path(__file__).with_name(name))
                                 for name in ('velocity_pairs.py', 'spin_hall_workflow.py', 'kubo_pairs.py')}}
    return write_pairs(args.output_dir, data)


def add_command(sub):
    p = sub.add_parser('velocity-pairs', help='audited PAW full-velocity bundle -> charge pairs for pair-hall')
    p.add_argument('--matrices', type=Path, required=True, help='physical-matrices.npz from spin-export')
    p.add_argument('--metadata', type=Path, required=True, help='matching physical-matrices.json')
    p.add_argument('--mesh', nargs=2, type=int, required=True, metavar=('NX','NY'))
    p.add_argument('--plane-axes', nargs=2, type=int, default=[0,1])
    p.add_argument('--energy-reference', required=True)
    p.add_argument('--memory-limit-mib', type=float, default=2048.)
    p.add_argument('--output-dir', type=Path, required=True)
