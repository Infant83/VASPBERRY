"""CLI operations for reusable native Kubo pairs and charge Hall scans."""
from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

from berry_data import require, write_hall
from exported_matrix_kubo import sha256
from kubo_pairs import (bundle_hall_spectrum, import_native_pairs, pair_hall_spectrum,
                        read_pairs, write_pairs)


WAVECAR_WARNING = ('WAVECAR-only canonical momentum approximation explicitly selected; '
                  'PAW optical matrix elements in WAVEDER are not used. This is not '
                  'equivalent to the standard WAVEDER Kubo route.')


def warn_wavecar_approximation():
    print('WARNING: ' + WAVECAR_WARNING, file=sys.stderr, flush=True)


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
        'vaspberry_version': (Path(__file__).resolve().parents[1]/'VERSION').read_text().strip(),
        'command': {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        'implementation_sha256': {name: sha256(Path(__file__).with_name(name)) for name in
                                 ('kubo_hall_workflow.py', 'kubo_pairs.py', 'native_kubo_csv.py',
                                  'berry_data.py', 'vaspberry_transport.py')},
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
    require(getattr(args, 'kubo_source', 'waveder') == 'wavecar',
            'WAVECAR-only execution requires explicit --kubo-source wavecar; '
            'use kubo-hall for the default WAVEDER protocol')
    warn_wavecar_approximation()
    mus, options = scan_options(args)
    require(args.mpi_procs >= 1, 'MPI process count must be positive')
    require(args.binary is not None, '--binary is required for --kubo-source wavecar')
    binary = Path(args.binary).expanduser().resolve()
    from waveder_hall import resolve_source_paths
    directory, paths = resolve_source_paths(getattr(args, 'input_dir', None),
        {'WAVECAR': args.wavecar} if args.wavecar is not None else None)
    wavecar = paths['WAVECAR']
    require(binary.is_file(), 'native binary does not exist: '+str(binary))
    require(wavecar.is_file(), 'WAVECAR input is missing or not a file: '+str(wavecar))
    out = args.output_dir.resolve()
    require(not out.exists(), 'output directory exists; choose a new directory')
    raw = out/'native'
    argv = [str(binary), '--task', 'kubo-pairs', '--kubo-source', 'wavecar',
            '--wavecar', str(wavecar), '--pairs-csv', 'PAIRS.csv']
    if args.spinor_components is not None:
        argv += ['--spinor', str(args.spinor_components)]
    if args.mpi_procs > 1:
        launcher = shutil.which(args.mpi_launcher)
        require(launcher is not None, 'MPI launcher not found')
        argv = [launcher, '-np', str(args.mpi_procs)] + argv
    raw.mkdir(parents=True)
    manifest = dict(schema='vaspberry.wavecar-hall-run', version=1, complete=False,
                    status='RUNNING', native_argv=argv, native_cwd=str(raw),
                    binary_sha256=sha256(binary), wavecar_sha256=sha256(wavecar),
                    input_directory=str(directory), source_paths={'WAVECAR': str(wavecar)},
                    kubo_source='wavecar', approximation_warning=WAVECAR_WARNING)
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


def kubo_hall_command(args):
    """Choose an explicit physical input route; never fall back on file errors."""
    if args.kubo_source == 'wavecar':
        require(args.run_dir is None and args.occupied is None and args.bands is None
                and all(getattr(args, key, None) is None for key in ('waveder', 'incar', 'outcar')),
                '--run-dir/--occupied/--bands/--waveder/--incar/--outcar apply to the WAVEDER route; '
                'WAVECAR approximation uses input-dir/WAVECAR or --wavecar and the requested mu/T scan')
        return wavecar_hall_command(args)
    from waveder_hall import command

    require(args.binary is None and args.mpi_procs == 1 and args.mpi_launcher == 'mpiexec',
            '--binary and native MPI options apply only to --kubo-source wavecar; '
            'the standard WAVEDER route evaluates optical matrices directly in Python')
    require((args.occupied is None) != (args.bands is None), '--occupied or --bands is required for WAVEDER')
    if args.occupied is not None:
        require(args.temperatures == [0.], 'WAVEDER occupied mode supports only insulating T=0; use --bands for a selected mu/T contribution')
        require(args.occupied > 0, '--occupied must be positive')
    require(args.pair_band_max is None and not args.allow_partial_bands
            and args.degeneracy_policy == 'error' and args.degeneracy_threshold_eV == 1e-7,
            'pair truncation/coalescing and redundant partial acknowledgment are unsupported for WAVEDER; --bands declares selection')
    return command(args)


def add_commands(sub):
    imp = sub.add_parser('import-pairs', help='native PAIRS.csv + WAVECAR -> reusable unordered-pair NPZ cache')
    imp.add_argument('--csv', type=Path, required=True)
    pair = sub.add_parser('pair-hall', help='cached unordered Kubo pairs -> mu/T charge sheet Hall')
    pair.add_argument('--pairs-dir', type=Path, required=True)
    bundle = sub.add_parser('bundle-hall', help='postprocess saved canonical occupied-trace CSV + WAVECAR (approximation)',
        description='Read saved canonical WAVECAR occupied-trace CSVs at insulating T=0. '
                    'PAW WAVEDER trace CSVs require kubo-hall --input-dir with the original optical run '
                    'for same-run validation; no operator relabelling is performed.')
    bundle.add_argument('--csv', type=Path, required=True)
    bundle.add_argument('--occupied', type=int, required=True)
    wave = sub.add_parser('wavecar-hall', help='explicit WAVECAR-only approximation (requires --kubo-source wavecar)',
        description='Legacy canonical momentum approximation. Explicit --kubo-source wavecar is required; '
                    'use kubo-hall for the standard WAVEDER protocol.')
    wave.add_argument('--binary', type=Path, required=True, help='serial or MPI VASPBERRY executable')
    wave.add_argument('--input-dir', type=Path,
                      help='input directory (default: invocation cwd); --wavecar overrides only WAVECAR')
    standard = sub.add_parser('kubo-hall', help='standard WAVEDER Kubo charge Hall; WAVECAR approximation requires opt-in',
        description='WAVEDER is the default Kubo matrix source. The validated standard VASP 5.4.4 '
                    'optical route supports selected --bands mu/T contributions or insulating --occupied T=0 scans. Missing or invalid '
                    'WAVEDER never triggers fallback. --kubo-source wavecar explicitly selects the '
                    'canonical momentum approximation and prints a warning.')
    standard.add_argument('--binary', type=Path, help='native executable; required only for --kubo-source wavecar')
    from waveder_hall import add_input_arguments
    add_input_arguments(standard, include_wavecar=False)
    from band_selection import parse_bands
    selection = standard.add_mutually_exclusive_group()
    selection.add_argument('--occupied', type=int, help='leading occupied bundle; insulating T=0 only')
    selection.add_argument('--bands', type=parse_bands, help='WAVEDER selected IDs/ranges, e.g. 31,33:34; all virtual bands retained')
    for p in (wave, standard):
        p.add_argument('--kubo-source', choices=['waveder', 'wavecar'], default='waveder',
                       help='default: waveder; wavecar explicitly opts into the canonical momentum approximation')
        p.add_argument('--mpi-procs', type=int, default=1, help='WAVECAR native export only; use an MPI build when >1')
        p.add_argument('--mpi-launcher', default='mpiexec')
    for p in (imp, bundle, wave, standard):
        p.add_argument('--wavecar', type=Path, required=p in (imp, bundle),
                       help='WAVECAR only; explicit relative paths resolve from invocation cwd')
        p.add_argument('--spin', type=int, default=1)
        p.add_argument('--spinor-components', type=int, choices=[1, 2],
                       help='optional assertion; default detects the WAVECAR layout')
        p.add_argument('--spin-multiplicity', type=int, choices=[1, 2],
                       help='default: 2 for scalar spin-degenerate states, otherwise 1')
        p.add_argument('--mesh', nargs=2, type=int, required=True, metavar=('NX', 'NY'))
        p.add_argument('--plane-axes', nargs=2, type=int, default=[0, 1])
        p.add_argument('--energy-reference', required=True, help='description of unchanged input energy zero')
    for p in (pair, bundle, wave, standard):
        p.add_argument('--mu-min', type=float, required=True)
        p.add_argument('--mu-max', type=float, required=True)
        p.add_argument('--mu-num', type=int, default=241)
        p.add_argument('--mu-reference', type=float, required=True)
        p.add_argument('--temperatures', nargs='+', type=float, default=[0.])
        p.add_argument('--regions', type=Path, help='optional named periodic circles or full-mesh k-ID sets')
        p.add_argument('--difference', action='append', default=[], metavar='NAME:LEFT:RIGHT')
        p.add_argument('--formats', nargs='+', choices=['csv', 'dat', 'npz'], default=['csv', 'dat', 'npz'],
                       help='independently selectable numerical formats; JSON metadata always included')
    for p in (pair, wave, standard):
        p.add_argument('--degeneracy-threshold-eV', type=float, default=1e-7)
        p.add_argument('--degeneracy-policy', choices=['error', 'coalesce'], default='error',
                       help='coalesce explicitly approximates numerical energy groups by their mean; records shifts')
        p.add_argument('--allow-partial-bands', action='store_true')
        p.add_argument('--mu-chunk', type=int, default=32)
        p.add_argument('--pair-band-max', type=int,
                       help='optional upper pair-band limit for virtual-state convergence on one larger WAVECAR; source bands remain recorded')
    for p in (imp, pair, bundle, wave, standard):
        p.add_argument('--output-dir', type=Path, required=True)


COMMANDS = {'import-pairs': import_pairs_command, 'pair-hall': pair_hall_command,
            'bundle-hall': bundle_hall_command, 'wavecar-hall': wavecar_hall_command,
            'kubo-hall': kubo_hall_command}
