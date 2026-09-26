"""CLI operations for reusable native Kubo pairs and charge Hall scans."""
from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np

from berry_data import require, write_hall
from exported_matrix_kubo import sha256
from kubo_pairs import (bundle_hall_spectrum, import_native_pairs, pair_hall_spectrum,
                        read_pairs, write_pairs)


def sampling(args):
    return {'kind': 'uniform_full_2d', 'mesh': args.mesh, 'plane_axes': args.plane_axes}


def scan_options(args):
    require(args.mu_num >= 1 and np.isfinite([args.mu_min, args.mu_max]).all()
            and args.mu_max >= args.mu_min and (args.mu_num > 1 or args.mu_min == args.mu_max),
            'valid mu range and positive number required')
    # The kernel evaluates the reference internally; keep output on the
    # requested scan rather than adding a distant point that stretches plots.
    mus = np.unique(np.linspace(args.mu_min, args.mu_max, args.mu_num))
    spec = json.loads(args.regions.read_text()) if args.regions else None
    differences = []
    for value in args.difference:
        bits = value.split(':')
        require(len(bits) == 3, 'difference syntax is NAME:LEFT:RIGHT')
        differences.append(tuple(bits))
    return mus, dict(mu_reference=args.mu_reference, region_spec=spec, differences=differences)


def pair_options(args):
    return dict(degeneracy_threshold_eV=args.degeneracy_threshold_eV,
                degeneracy_policy=args.degeneracy_policy, allow_partial_bands=args.allow_partial_bands,
                mu_chunk=args.mu_chunk, pair_band_max=args.pair_band_max)


def source_options(args):
    return dict(spin=args.spin, spinor_components=args.spinor_components,
                spin_multiplicity=args.spin_multiplicity, sampling=sampling(args),
                energy_reference=args.energy_reference)


def record_provenance(args, meta):
    meta['provenance'] = {
        'vaspberry_version': '1.4.1',
        'command': {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        'implementation_sha256': {name: sha256(Path(__file__).with_name(name)) for name in
                                 ('kubo_hall_workflow.py', 'kubo_pairs.py', 'berry_data.py', 'vaspberry_transport.py')},
    }
    if getattr(args, 'regions', None):
        meta['provenance']['regions_sha256'] = sha256(args.regions)
    return meta


def import_pairs_command(args):
    data = import_native_pairs(args.csv, args.wavecar, **source_options(args))
    record_provenance(args, data.metadata)
    return write_pairs(args.output_dir, data)


def pair_hall_command(args):
    data = read_pairs(args.pairs_dir)
    mus, options = scan_options(args)
    rows, meta = pair_hall_spectrum(data, mus, args.temperatures, **options, **pair_options(args))
    return write_hall(args.output_dir, rows, record_provenance(args, meta), formats=args.formats)


def bundle_hall_command(args):
    mus, options = scan_options(args)
    require(args.temperatures == [0.], 'fixed-bundle Hall supports only T=0 inside its common global gap')
    rows, meta = bundle_hall_spectrum(args.csv, args.wavecar, mus, occupied=args.occupied,
                                     **source_options(args), **options)
    return write_hall(args.output_dir, rows, record_provenance(args, meta), formats=args.formats)


def wavecar_hall_command(args):
    """One user command, separate native export, reusable cache and integration."""
    mus, options = scan_options(args)
    require(args.mpi_procs >= 1, 'MPI process count must be positive')
    binary = Path(args.binary).expanduser().resolve()
    wavecar = args.wavecar.expanduser().resolve()
    require(binary.is_file() and wavecar.is_file(), 'native binary and WAVECAR must exist')
    out = args.output_dir.resolve()
    require(not out.exists(), 'output directory exists; choose a new directory')
    raw = out/'native'
    argv = [str(binary), '-f', str(wavecar), '-s', str(args.spinor_components),
            '-kubo', '2', '-kubo_pairs', 'PAIRS.csv']
    if args.mpi_procs > 1:
        launcher = shutil.which(args.mpi_launcher)
        require(launcher is not None, 'MPI launcher not found')
        argv = [launcher, '-np', str(args.mpi_procs)] + argv
    raw.mkdir(parents=True)
    manifest = dict(schema='vaspberry.wavecar-hall-run', version=1, complete=False,
                    status='RUNNING', native_argv=argv, native_cwd=str(raw),
                    binary_sha256=sha256(binary), wavecar_sha256=sha256(wavecar))
    def save():
        (out/'workflow.json').write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
    save(); start = time.monotonic()
    try:
        with (raw/'stdout.log').open('w') as stdout, (raw/'stderr.log').open('w') as stderr:
            completed = subprocess.run(argv, cwd=raw, stdout=stdout, stderr=stderr, check=False)
        manifest['native_exit_code'] = completed.returncode
        manifest['native_elapsed_s'] = time.monotonic()-start
        require(completed.returncode == 0, 'native Kubo export failed; inspect native/stdout.log and stderr.log')
        data = import_native_pairs(raw/'PAIRS.csv', wavecar, **source_options(args))
        record_provenance(args, data.metadata)
        write_pairs(out/'pairs', data)
        rows, meta = pair_hall_spectrum(data, mus, args.temperatures, **options, **pair_options(args))
        meta = write_hall(out/'hall', rows, record_provenance(args, meta), formats=args.formats)
        manifest.update(complete=True, status='PASS', output='hall/conductivity.json',
                        pair_cache='pairs', elapsed_s=time.monotonic()-start)
        save()
        return meta
    except BaseException as exc:
        manifest.update(status='FAILED', error=str(exc), elapsed_s=time.monotonic()-start)
        save()
        raise


def add_commands(sub):
    imp = sub.add_parser('import-pairs', help='native PAIRS.csv + WAVECAR -> reusable unordered-pair NPZ cache')
    imp.add_argument('--csv', type=Path, required=True)
    pair = sub.add_parser('pair-hall', help='cached unordered Kubo pairs -> mu/T charge sheet Hall')
    pair.add_argument('--pairs-dir', type=Path, required=True)
    bundle = sub.add_parser('bundle-hall', help='occupied bundle CSV + WAVECAR -> insulating T=0 charge sheet Hall')
    bundle.add_argument('--csv', type=Path, required=True)
    bundle.add_argument('--occupied', type=int, required=True)
    wave = sub.add_parser('wavecar-hall', help='WAVECAR -> native Kubo export -> cache -> charge sheet Hall in one command')
    wave.add_argument('--binary', type=Path, required=True, help='serial or MPI VASPBERRY executable')
    wave.add_argument('--mpi-procs', type=int, default=1, help='native export only; use an MPI build when >1')
    wave.add_argument('--mpi-launcher', default='mpiexec')
    for p in (imp, bundle, wave):
        p.add_argument('--wavecar', type=Path, required=True)
        p.add_argument('--spin', type=int, default=1)
        p.add_argument('--spinor-components', type=int, choices=[1, 2], required=True)
        p.add_argument('--spin-multiplicity', type=int, choices=[1, 2], required=True)
        p.add_argument('--mesh', nargs=2, type=int, required=True, metavar=('NX', 'NY'))
        p.add_argument('--plane-axes', nargs=2, type=int, default=[0, 1])
        p.add_argument('--energy-reference', required=True, help='description of unchanged input energy zero')
    for p in (pair, bundle, wave):
        p.add_argument('--mu-min', type=float, required=True)
        p.add_argument('--mu-max', type=float, required=True)
        p.add_argument('--mu-num', type=int, default=241)
        p.add_argument('--mu-reference', type=float, required=True)
        p.add_argument('--temperatures', nargs='+', type=float, default=[0.])
        p.add_argument('--regions', type=Path, help='optional named periodic circles or full-mesh k-ID sets')
        p.add_argument('--difference', action='append', default=[], metavar='NAME:LEFT:RIGHT')
        p.add_argument('--formats', nargs='+', choices=['csv', 'dat', 'npz'], default=['csv', 'dat', 'npz'],
                       help='independently selectable numerical formats; JSON metadata always included')
    for p in (pair, wave):
        p.add_argument('--degeneracy-threshold-eV', type=float, default=1e-7)
        p.add_argument('--degeneracy-policy', choices=['error', 'coalesce'], default='error',
                       help='coalesce explicitly approximates numerical energy groups by their mean; records shifts')
        p.add_argument('--allow-partial-bands', action='store_true')
        p.add_argument('--mu-chunk', type=int, default=32)
        p.add_argument('--pair-band-max', type=int,
                       help='optional upper pair-band limit for virtual-state convergence on one larger WAVECAR; source bands remain recorded')
    for p in (imp, pair, bundle, wave):
        p.add_argument('--output-dir', type=Path, required=True)


COMMANDS = {'import-pairs': import_pairs_command, 'pair-hall': pair_hall_command,
            'bundle-hall': bundle_hall_command, 'wavecar-hall': wavecar_hall_command}
