"""Assemble separately audited fixed-charge VASP k chunks without rotating states."""
from __future__ import annotations

import json
from pathlib import Path
import struct

import numpy as np

from berry_data import require
from exported_matrix_kubo import sha256
from spin_hall_workflow import read_matrices
from vasp_spin_export import audit_producer, SCOPE, LIMITATIONS
from waveder_hall import incar_values
from wavecar_fukui import Wavecar, infer_uniform_grid


def assemble_wavecar(paths, waves, order, target):
    """Copy every complete per-k record byte-for-byte, padding the record stride.

    Only global record length and total k count change. The first chunk's
    Fermi header is retained, not interpreted as a full-mesh Fermi level. Per-k eigenvalues,
    occupations, G counts and spinor coefficients remain unchanged.
    """
    nb=waves[0].header.nbands;stride=max(w.header.stride_bytes for w in waves)
    require(all((4+3*nb)*8 <= w.header.stride_bytes
                and np.all(8*w.nplane <= w.header.stride_bytes) for w in waves),
            'source WAVECAR header or coefficient record exceeds its stride')
    require(all(w.header.logical_recl == w.header.stride_bytes for w in waves),
            'assembly supports byte-RECL WAVECAR records only')
    offsets=np.r_[0,np.cumsum([w.header.nkpoints for w in waves])]
    with Path(paths[0]).open('rb') as f:
        h1=bytearray(f.read(waves[0].header.stride_bytes));h2=bytearray(f.read(waves[0].header.stride_bytes))
    struct.pack_into('<d',h1,0,float(stride));struct.pack_into('<d',h2,0,float(len(order)))
    with Path(target).open('xb') as out:
        out.write(h1);out.write(bytes(stride-len(h1)));out.write(h2);out.write(bytes(stride-len(h2)))
        for combined_index in order:
            c=int(np.searchsorted(offsets[1:],combined_index,side='right'));k=int(combined_index-offsets[c])
            w=waves[c];source_stride=w.header.stride_bytes
            with Path(paths[c]).open('rb') as f:
                f.seek((w.record_number(k)-1)*source_stride)
                for _ in range(nb+1):
                    raw=f.read(source_stride)
                    require(len(raw)==source_stride,'truncated source WAVECAR during assembly')
                    out.write(raw);out.write(bytes(stride-source_stride))
    assembled=Wavecar(Path(target),spinor_components=2)
    require(assembled.header.nkpoints==len(order),'assembled WAVECAR k count failed')
    # Verify bytes independently by re-reading the written records, not merely
    # comparing a hash computed while writing the same stream.
    with Path(target).open('rb') as assembled_file:
        for new,combined_index in enumerate(order):
            c=int(np.searchsorted(offsets[1:],combined_index,side='right'));k=int(combined_index-offsets[c])
            w=waves[c];source_stride=w.header.stride_bytes
            with Path(paths[c]).open('rb') as original:
                original.seek((w.record_number(k)-1)*source_stride)
                for band in range(nb+1):
                    assembled_file.seek((assembled.record_number(new)+band-1)*stride)
                    require(assembled_file.read(source_stride)==original.read(source_stride),
                            'assembled per-k record differs from source bytes')
    return assembled


def merge_runs(run_dirs, output_dir, mesh, occupied, *, memory_limit_mib=2048.):
    """Validate same-Hamiltonian chunks, restore one complete periodic mesh."""
    cases=[Path(p).resolve() for p in run_dirs];out=Path(output_dir)
    require(len(cases)>=2 and len(set(cases))==len(cases),'at least two distinct completed run directories required')
    require(not out.exists(),'output directory exists')
    partial=out.with_name(out.name+'.partial');require(not partial.exists(),'unfinished output directory exists')
    require(len(mesh)==2 and all(type(n) is int and n>=2 for n in mesh),'two mesh sizes >=2 required')
    sizes=[(c/'SPIN_VELOCITY.bin').stat().st_size for c in cases]
    estimate=4*sum(sizes)+max(sizes)
    require(np.isfinite(memory_limit_mib) and memory_limit_mib>0
            and estimate<=memory_limit_mib*2**20,'spin assembly memory estimate exceeds limit')
    datasets=[];reports=[];waves=[];conditions=[];fixed_hashes=[]
    for case in cases:
        data,report=audit_producer(case)
        wave=Wavecar(case/'WAVECAR',spinor_components=2)
        hashes={name:sha256(case/name) for name in ('CHGCAR','POTCAR','POSCAR','KPOINTS')}
        run=json.loads((case/'run.json').read_text())
        require(all(run.get('input_sha256',{}).get(name)==value for name,value in hashes.items()),
                'fixed-charge chunk source association failed')
        incar=incar_values((case/'INCAR').read_text())
        require(int(incar.get('ICHARG','0'))==11,'chunk assembly requires ICHARG=11 fixed charge')
        for key in ('SYSTEM','ISTART','NELM','NELMIN','ALGO'):incar.pop(key,None)
        require(np.allclose(data['weights'],1/len(data['weights']),rtol=0,atol=1e-12),
                'every chunk must contain uniformly weighted points')
        if datasets:
            require(all(hashes[n]==fixed_hashes[0][n] for n in ('CHGCAR','POTCAR','POSCAR'))
                    and incar==conditions[0], 'chunks do not share the same fixed Hamiltonian and physical settings')
            require(report['producer_binary_sha256']==reports[0]['producer_binary_sha256'],
                    'chunks use different producer executables')
            require(wave.header.nbands==waves[0].header.nbands
                    and wave.header.encut_ev==waves[0].header.encut_ev
                    and np.array_equal(wave.header.lattice,waves[0].header.lattice),
                    'chunk band coverage, cutoff or lattice disagree')
        datasets.append(data);reports.append(report);waves.append(wave);conditions.append(incar);fixed_hashes.append(hashes)
    q=np.concatenate([d['kpoints_fractional'] for d in datasets]);grid=infer_uniform_grid(q,*mesh)
    order=grid.index.ravel();nk=len(q);nb=waves[0].header.nbands
    require(type(occupied) is int and 0<occupied<nb,'occupied count must be inside source band coverage')
    energies=np.concatenate([d['energies_eV'] for d in datasets])[order]
    occ=np.concatenate([d['occupations'] for d in datasets])[order]
    expected=np.zeros(nb);expected[:occupied]=1.
    require(np.max(abs(occ-expected))<1e-6,'all chunks must share the same insulating occupied bundle')
    vbm=float(energies[:,occupied-1].max());cbm=float(energies[:,occupied].min())
    require(cbm>vbm,'assembled mesh has no positive global gap')
    full=np.concatenate([d['pseudo']+d['delta'] for d in datasets])[order]
    arrays=dict(energies_eV=energies,kpoints_fractional=q[order],weights=np.full(nk,1/nk),
        lattice_A=datasets[0]['lattice_A'],spin_pauli=full[:,1:],overlap=full[:,0],
        velocity_eVA=np.concatenate([d['velocity_eVA'] for d in datasets])[order],
        band_indices=np.arange(1,nb+1,dtype=np.int64))
    partial.mkdir(parents=True)
    status={'status':'RUNNING'}
    def record():(partial/'run.json').write_text(json.dumps(status,indent=2)+'\n')
    record()
    try:
        assembled=assemble_wavecar([c/'WAVECAR' for c in cases],waves,order,partial/'WAVECAR')
        require(np.array_equal(assembled.kpoints,arrays['kpoints_fractional'])
                and np.array_equal(assembled.energies,energies),'assembled WAVECAR and operator gauge association failed')
        np.savez_compressed(partial/'physical-matrices.npz',**arrays)
        source={'WAVECAR':sha256(partial/'WAVECAR')}
        for i,(report,fixed) in enumerate(zip(reports,fixed_hashes)):
            source.update({f'chunk{i:03d}/{name}':digest for name,digest in {**report['source_sha256'],**fixed}.items()})
        meta=dict(schema='vaspberry.spin-velocity',version=1,complete=True,producer_status='PASS',
            data_npz_sha256=sha256(partial/'physical-matrices.npz'),source_wavecar_sha256=source['WAVECAR'],
            source_nbands=nb,nbands=nb,nkpoints=nk,
            units={'spin':'dimensionless Pauli','overlap':'dimensionless','velocity':'eV angstrom = hbar*v','energies':'eV','lattice':'angstrom'},
            spin_basis='Cartesian axes; SAXIS=(0,0,1)',spin_basis_id='cartesian',velocity_basis='cartesian',
            matrix_axes=['k','cartesian','bra','ket'],gauge='same_eigenvectors_as_source_wavecar',
            operator_accuracy='paw_full_velocity',spin_operator_accuracy='paw_augmented',physical_paw_validated=True,
            normalization_applied=False,spinor_components=2,spin_multiplicity=1,full_source_band_coverage=True,
            diagonal_available=True,degenerate_blocks_available=True,
            producer='Audited fixed-charge VASP spin/velocity chunks assembled by VASPBERRY',operator_scope=SCOPE,
            spin_current_construction='finite_band_projected_product',method_limitations=LIMITATIONS,
            source_sha256=source,producer_binary_sha256=reports[0]['producer_binary_sha256'],
            source_energy_reference='Unchanged same-density VASP eigenvalues. The first chunk Fermi header is retained and is not a full-mesh chemical potential; transport selects its own in-gap chemical potential.',
            assembly={'mesh':list(mesh),'source_chunk_names':[c.name for c in cases],
                'source_chunk_point_counts':[w.header.nkpoints for w in waves],
                'per_k_records':'All source bytes independently checked after assembly; only global headers and record padding changed.',
                'weights':'Full uniform mesh: each source point receives1/(NX*NY); no incomplete-grid renormalization.',
                'occupied':occupied,'global_gap_eV':cbm-vbm,'global_midgap_eV':.5*(vbm+cbm),
                'estimated_memory_MiB':estimate/2**20,'memory_limit_MiB':memory_limit_mib},
            implementation_sha256=sha256(Path(__file__)))
        (partial/'physical-matrices.json').write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n')
        _,_,validation=read_matrices(partial/'physical-matrices.npz',partial/'physical-matrices.json',memory_limit_mib=memory_limit_mib)
        (partial/'assembly-audit.json').write_text(json.dumps(dict(status='PASS',chunks=reports,validation=validation),indent=2)+'\n')
        for case,report,fixed in zip(cases,reports,fixed_hashes):
            require(all(sha256(case/name)==digest for name,digest in {**report['source_sha256'],**fixed}.items()),
                    'source chunk changed during assembly')
        status.update(status='PASS',source_chunks=len(cases),nkpoints=nk,nbands=nb);record()
        require(not out.exists(),'output appeared during assembly');partial.rename(out)
        return meta
    except BaseException as exc:
        status.update(status='FAILED',error=str(exc));record();raise


def add_command(sub):
    p=sub.add_parser('spin-merge',help='audit same-density VASP k chunks and assemble full-mesh spin/velocity plus WAVECAR')
    p.add_argument('--run-dirs',nargs='+',type=Path,required=True)
    p.add_argument('--mesh',nargs=2,type=int,required=True)
    p.add_argument('--occupied',type=int,required=True)
    p.add_argument('--memory-limit-mib',type=float,default=2048.)
    p.add_argument('--output-dir',type=Path,required=True)


def command(args):return merge_runs(args.run_dirs,args.output_dir,args.mesh,args.occupied,memory_limit_mib=args.memory_limit_mib)
