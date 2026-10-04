#!/usr/bin/env python3
"""Execute VASPBERRY and postprocess its results from one settings file.

  python3 tools/vaspberry_post.py check analysis.ini
  python3 tools/vaspberry_post.py run analysis.ini
  python3 tools/vaspberry_post.py plot results/run01

The standard run uses WAVEDER from a completed VASP optical run for selected-band
mu/T contributions or insulating occupied-bundle T=0 charge Hall. Set
[run] kubo_source = wavecar to explicitly select the
WAVECAR-only canonical momentum approximation (a warning is emitted). Missing
or invalid WAVEDER never causes automatic fallback. The plot command reads
saved tables without changing their source operator.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from postprocess_config import load_settings

TOOLS = Path(__file__).resolve().parent
ROOT = TOOLS.parent
SCHEMA = 'vaspberry.postprocess-run'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def source_arguments(run):
    return ['--wavecar', run['wavecar'], '--spin', str(run['spin']),
            '--spinor-components', str(run['spinor_components']),
            '--spin-multiplicity', str(run['spin_multiplicity']),
            '--mesh', *map(str, run['mesh']), '--plane-axes', *map(str, run['plane_axes']),
            '--energy-reference=' + run['energy_reference']]


def scan_arguments(scan):
    return ['--mu-min', str(scan['mu_min']), '--mu-max', str(scan['mu_max']),
            '--mu-num', str(scan['mu_num']), '--mu-reference', str(scan['mu_reference']),
            '--temperatures', *map(str, scan['temperatures'])]


def python_command(script, *arguments):
    return [sys.executable, str(TOOLS / script), *map(str, arguments)]


def read_run(directory):
    root = Path(directory).expanduser().resolve()
    record = json.loads((root / 'run.json').read_text(encoding='utf-8'))
    require(record.get('schema') == SCHEMA and record.get('version') == 1,
            'unrecognized run.json; use a result created by vaspberry_post.py')
    require(record.get('status') == 'PASS' and record.get('complete') is True,
            'the selected previous run did not complete successfully')
    return root, record


def verify_saved(root, record, names):
    recorded = record.get('outputs', {})
    for name in names:
        require(name in recorded and sha256(root / name) == recorded[name],
                'saved output changed or is missing: ' + str(root / name))


def preflight(settings, reuse=None):
    """Check paths and declared source identity before starting native work."""
    from wavecar_fukui import Wavecar, infer_uniform_grid

    run = settings['run']
    output = Path(run['output'])
    require(not output.exists(), 'output directory exists; set [run] output to a new directory')
    source = run.get('kubo_source', 'waveder')
    require(source in ('waveder', 'wavecar'), '[run] kubo_source must be waveder or wavecar')
    require(Path(run['input_dir']).is_dir(), 'input directory does not exist: '+run['input_dir'])
    if source == 'waveder':
        return preflight_waveder(settings, reuse)
    require(settings['hall'].get('occupied') is None and settings['hall'].get('bands') is None,
            '[hall] occupied/bands is only supported for kubo_source = waveder; '
            'the WAVECAR approximation uses the requested mu/T occupations')
    from kubo_hall_workflow import warn_wavecar_approximation
    warn_wavecar_approximation()
    wavecar = Path(run['wavecar'])
    require(wavecar.is_file(), 'WAVECAR does not exist: ' + str(wavecar))
    auto = run['spin_mode'] == 'auto'
    wave = Wavecar(wavecar, spin=run['spin'],
                   spinor_components=None if auto else run['spinor_components'])
    if auto:
        require(wave.header.ispin == 1,
                '[run] this WAVECAR has two collinear channels; select '
                'spin_mode = collinear-up or collinear-down')
    if not auto:
        require(wave.header.ispin == run['expected_nspin'],
                '[run] spin_mode does not match the WAVECAR ISPIN header')
    run['spinor_components'] = wave.spinor_components
    run['spin_multiplicity'] = wave.resolve_spin_multiplicity(
        None if auto else run['spin_multiplicity'])
    run['expected_nspin'] = wave.header.ispin
    run['resolved_spin_mode'] = ('spinor' if wave.spinor_components == 2 else
                                 'scalar-degenerate' if wave.header.ispin == 1 else
                                 'collinear-up' if run['spin'] == 1 else 'collinear-down')
    require(len(wave.kpoints) == run['mesh'][0] * run['mesh'][1],
            '[run] mesh size does not match the number of WAVECAR k points')
    axes = run['plane_axes']
    other = next(i for i in range(3) if i not in axes)
    infer_uniform_grid(wave.kpoints[:, axes + [other]], *run['mesh'])
    wave.coefficients(0, [1])
    projection = settings['projection']
    require(not projection or wave.spinor_components == 2,
            '[projection] requires a two-component spinor WAVECAR')
    inputs = {'wavecar': {'path': str(wavecar), 'sha256': sha256(wavecar)}}
    if projection:
        for name in ('procar', 'outcar'):
            path = Path(projection[name])
            require(path.is_file(), name.upper() + ' does not exist: ' + str(path))
            inputs[name] = {'path': str(path), 'sha256': sha256(path)}
        require(max(projection['bands']) <= wave.energies.shape[1],
                '[projection] bands exceed the WAVECAR band count')
        require(settings['plot']['map_band'] <= wave.energies.shape[1],
                '[plot] map_band exceeds the WAVECAR band count')
    if settings['hall']['pair_band_max'] is not None:
        require(settings['hall']['pair_band_max'] <= wave.energies.shape[1],
                '[hall] pair_band_max exceeds the WAVECAR band count')

    cache = None
    if reuse is None:
        require(run.get('binary') is not None, '[run] binary is required for kubo_source = wavecar')
        binary = Path(run['binary'])
        require(binary.is_file() and os.access(binary, os.X_OK),
                'native executable is missing or not executable: ' + str(binary))
        inputs['binary'] = {'path': str(binary), 'sha256': sha256(binary)}
        if run['mpi_procs'] > 1:
            require(shutil.which(run['mpi_launcher']) is not None,
                    'MPI launcher not found: ' + run['mpi_launcher'])
    else:
        from kubo_pairs import read_pairs

        previous, record = read_run(reuse)
        require(record['inputs']['wavecar']['sha256'] == inputs['wavecar']['sha256'],
                '--reuse requires the same WAVECAR; use a fresh native run for a new VASP state')
        old_run = record['settings']['run']
        # Runs predating source selection used the canonical WAVECAR operator.
        require(old_run.get('kubo_source', 'wavecar') == 'wavecar',
                '--reuse requires a WAVECAR approximation pair cache')
        for key in ('mesh', 'plane_axes', 'spin', 'spinor_components',
                    'spin_multiplicity', 'energy_reference'):
            require(old_run[key] == run[key], '--reuse cannot change source setting: ' + key)
        verify_saved(previous, record, ['pairs/pairs.npz', 'pairs/pairs.json'])
        read_pairs(previous / 'pairs')
        cache = previous / 'pairs'
        inputs['reused_run'] = {'path': str(previous), 'sha256': sha256(previous / 'run.json')}
    return inputs, cache


def preflight_waveder(settings, reuse):
    """Validate and evaluate the optical input before any output directory exists."""
    import numpy as np
    from waveder_hall import waveder_hall_spectrum

    run, hall = settings['run'], settings['hall']
    require(run.get('binary') is None and run['mpi_procs'] == 1 and run['mpi_launcher'] == 'mpiexec',
            '[run] binary and native MPI settings apply only to kubo_source = wavecar; '
            'the standard WAVEDER route evaluates optical matrices directly in Python')
    require(reuse is None, '--reuse is only supported for explicitly selected WAVECAR pair caches')
    require((hall.get('occupied') is None) != (hall.get('bands') is None),
            '[hall] occupied or bands is required for the WAVEDER route')
    if hall.get('occupied') is not None:
        require(hall['temperatures'] == [0.], 'WAVEDER occupied mode supports only insulating T=0; use bands for a selected mu/T contribution')
    require(hall['pair_band_max'] is None, '[hall] pair_band_max is unsupported for WAVEDER; the full source virtual range is retained')
    require(settings['projection'] is None, '[projection] is unsupported for WAVEDER charge Hall; '
            'do not substitute canonical pair projections silently')
    optical = Path(run['input_dir'])
    paths = {name: Path(run[name.lower()]) for name in ('WAVEDER', 'WAVECAR', 'INCAR', 'OUTCAR')}
    rows, meta = waveder_hall_spectrum(optical,
        np.unique(np.linspace(hall['mu_min'], hall['mu_max'], hall['mu_num'])),
        occupied=hall['occupied'], bands=hall.get('bands'), temperatures=hall['temperatures'], spin=run['spin'], spinor_components=run['spinor_components'],
        spin_multiplicity=run['spin_multiplicity'],
        sampling={'kind': 'uniform_full_2d', 'mesh': run['mesh'], 'plane_axes': run['plane_axes']},
        energy_reference=run['energy_reference'], mu_reference=hall['mu_reference'],
        region_spec=settings['regions'], differences=settings['differences'], source_paths=paths)
    source = meta['source_metadata']
    if run['spin_mode'] == 'auto':
        require(source['source_nspin'] == 1, '[run] collinear WAVECAR requires spin_mode = collinear-up or collinear-down')
    else:
        require(source['source_nspin'] == run['expected_nspin'], '[run] spin_mode does not match WAVECAR ISPIN')
    run['spinor_components'] = source['spinor_components']
    run['spin_multiplicity'] = meta['spin_multiplicity']
    run['expected_nspin'] = source['source_nspin']
    run['resolved_spin_mode'] = ('spinor' if run['spinor_components'] == 2 else
                                 'scalar-degenerate' if run['expected_nspin'] == 1 else
                                 'collinear-up' if run['spin'] == 1 else 'collinear-down')
    inputs = {name.lower(): {'path': source['source_paths'][name], 'sha256': checksum}
              for name, checksum in source['source_sha256'].items()}
    return inputs, {'waveder_hall': (rows, meta)}


def settings_summary(settings):
    run, hall = settings['run'], settings['hall']
    print('WAVECAR: ' + run['wavecar'])
    print('Input directory: ' + run['input_dir'])
    print('Kubo matrix source: ' + run.get('kubo_source', 'waveder'))
    if run.get('kubo_source', 'waveder') == 'waveder':
        for name in ('waveder', 'incar', 'outcar'):
            print(name.upper() + ': ' + run[name])
        print('Optical run: ' + run['optical_run_dir'] + ('; selected bands '+str(hall['bands'])
              if hall.get('bands') is not None else '; insulating occupied bundle at T=0'))
        if hall.get('bands') is not None:
            print('Hall scope: selected Fermi-weighted contribution; not automatically total AHC.')
    print('Output:  ' + run['output'])
    mode = run.get('resolved_spin_mode', run['spin_mode'])
    print(f"Source: {mode}; full {run['mesh'][0]} x {run['mesh'][1]} mesh")
    print('Energy zero: ' + run['energy_reference'])
    print(f"Hall: mu {hall['mu_min']:g} to {hall['mu_max']:g} eV ({hall['mu_num']} points), "
          f"reference {hall['mu_reference']:g} eV; T={hall['temperatures']} K")
    if settings['projection']:
        print('PROCAR: selected bands ' + ' '.join(map(str, settings['projection']['bands']))
              + '; virtual-state sum uses all stored pair bands')


def execute_stage(root, record, name, command, cwd=None, native=False):
    cwd = root if cwd is None else Path(cwd)
    logs = root / ('native' if native else 'logs')
    logs.mkdir(exist_ok=True)
    stem = '' if native else name + '.'
    stdout, stderr = logs / (stem + 'stdout.log'), logs / (stem + 'stderr.log')
    stage = dict(name=name, command=command, cwd=str(cwd), status='RUNNING',
                 stdout=str(stdout.relative_to(root)), stderr=str(stderr.relative_to(root)))
    record['stages'].append(stage)
    write_json(root / record['manifest_name'], record)
    print('Running: ' + name, flush=True)
    started = time.monotonic()
    try:
        with stdout.open('w', encoding='utf-8') as out, stderr.open('w', encoding='utf-8') as err:
            result = subprocess.run(command, cwd=cwd, stdout=out, stderr=err, check=False)
        stage['exit_code'] = result.returncode
        require(result.returncode == 0,
                f"{name} failed (exit {result.returncode}); see {stderr} and {stdout}")
        stage['status'] = 'PASS'
    except BaseException as exc:
        stage.update(status='FAILED', error=str(exc))
        raise
    finally:
        stage['elapsed_s'] = time.monotonic() - started
        write_json(root / record['manifest_name'], record)


def run_calculation(settings, reuse=None):
    inputs, cache = preflight(settings, reuse)
    run = settings['run']
    root = Path(run['output'])
    root.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(settings['config_path'], root / 'settings.ini')
    write_json(root / 'regions.json', settings['regions'])
    write_json(root / 'groups.json', settings['groups'])
    record = dict(schema=SCHEMA, version=1, complete=False, status='RUNNING',
                  software_version=(ROOT / 'VERSION').read_text().strip(),
                  settings=settings, inputs=inputs, stages=[], manifest_name='run.json',
                  settings_sha256=sha256(root / 'settings.ini'),
                  implementation_sha256={name: sha256(TOOLS / name) for name in
                                         ('vaspberry_post.py', 'postprocess_config.py')},
                  started_at=datetime.now(timezone.utc).isoformat())
    write_json(root / 'run.json', record)
    started = time.monotonic()
    try:
        if run.get('kubo_source', 'waveder') == 'waveder':
            from berry_data import write_hall

            # preflight evaluated the validated optical matrices. Verify that
            # none of those inputs changed before publishing the cached result.
            require(all(sha256(item['path']) == item['sha256'] for item in inputs.values()),
                    'optical inputs changed after validation')
            rows, metadata = cache['waveder_hall']
            metadata['provenance'] = dict(kubo_source='waveder',
                settings_sha256=record['settings_sha256'],
                implementation_sha256={name: sha256(TOOLS/name) for name in
                                       ('vaspberry_post.py', 'waveder_hall.py', 'waveder_selected.py', 'band_selection.py')})
            write_hall(root/'hall', rows, metadata, formats=['csv', 'dat', 'npz'])
            record['stages'].append(dict(name='waveder-charge-hall', status='PASS', exit_code=0,
                implementation='waveder_hall.waveder_hall_spectrum',
                note='Optical input validated and evaluated before creating the output directory.'))
        elif cache is None:
            raw = root / 'native'
            raw.mkdir()
            command = [run['binary'], '--task', 'kubo-pairs', '--kubo-source', 'wavecar', '--wavecar', run['wavecar'],
                       '--spinor', str(run['spinor_components']), '--pairs-csv', 'PAIRS.csv']
            if run['mpi_procs'] > 1:
                command = [shutil.which(run['mpi_launcher']), '-n', str(run['mpi_procs']), *command]
            execute_stage(root, record, 'native-pairs', command, cwd=raw, native=True)
            execute_stage(root, record, 'import-pairs', python_command('vaspberry_kubo.py',
                'import-pairs', '--csv', raw / 'PAIRS.csv', *source_arguments(run),
                '--output-dir', root / 'pairs'))
        elif run.get('kubo_source', 'waveder') == 'wavecar':
            shutil.copytree(cache, root / 'pairs')
            record['reused_pair_cache'] = str(cache)
            print('Reusing saved pairs; VASPBERRY execution is skipped.', flush=True)
        region_args = ['--regions', str(root / 'regions.json')] if settings['regions']['regions'] else []
        difference_args = [item for row in settings['differences'] for item in ('--difference', ':'.join(row))]
        cap = settings['hall']['pair_band_max']
        cap_args = [] if cap is None else ['--pair-band-max', str(cap)]
        if run.get('kubo_source', 'waveder') == 'wavecar':
            from kubo_hall_workflow import WAVECAR_WARNING
            record['approximation_warning'] = WAVECAR_WARNING
            execute_stage(root, record, 'charge-hall', python_command('vaspberry_kubo.py',
                'pair-hall', '--pairs-dir', root / 'pairs', *scan_arguments(settings['hall']),
                *region_args, *difference_args, *cap_args, '--formats', 'csv', 'dat', 'npz',
                '--output-dir', root / 'hall'))
        projection = settings['projection']
        if projection:
            execute_stage(root, record, 'procar-character', python_command('procar_character.py',
                'project', '--procar', projection['procar'], '--wavecar', run['wavecar'],
                '--outcar', projection['outcar'], '--groups', root / 'groups.json',
                '--axis', *map(str, projection['axis']), '--output-dir', root / 'character'))
            execute_stage(root, record, 'projected-hall', python_command('procar_character.py',
                'hall', '--character-dir', root / 'character', '--pairs-dir', root / 'pairs',
                '--bands', *map(str, projection['bands']), *scan_arguments(projection),
                *region_args, '--output-dir', root / 'character-hall'))
        record['outputs'] = {
            str(path.relative_to(root)): sha256(path)
            for directory in ('pairs', 'hall', 'character', 'character-hall')
            for path in sorted((root / directory).glob('*')) if path.is_file()
        }
        record.update(status='PASS', complete=True)
    except BaseException as exc:
        record.update(status='FAILED', error=str(exc))
        raise
    finally:
        record['elapsed_s'] = time.monotonic() - started
        write_json(root / 'run.json', record)
    if run.get('kubo_source', 'waveder') == 'wavecar':
        print('Saved native pair cache: ' + str(root / 'pairs'))
    print('Saved Hall tables: ' + str(root / 'hall' / 'conductivity.csv'))
    if projection:
        print('Saved character and projected Hall tables: ' + str(root / 'character-hall'))
    print('To draw saved results: python3 tools/vaspberry_post.py plot ' + str(root))
    return record


def plot_results(directory, output_dir=None, *, temperature=None, group=None, band=None, quantity=None):
    source, original = read_run(directory)
    settings = original['settings']
    preferences = dict(settings['plot'])
    verify_saved(source, original, ['hall/conductivity.csv', 'hall/conductivity.json'])
    if settings['projection']:
        verify_saved(source, original, ['character/character.npz', 'character/character.json',
                                        'character-hall/character_hall.npz',
                                        'character-hall/character_hall.json'])
    require(settings['projection'] or (group is None and band is None),
            '--group/--band require saved PROCAR projection results')
    if temperature is not None:
        require(temperature in settings['hall']['temperatures'], 'temperature absent from saved Hall scan')
        preferences['hall_temperatures'] = [temperature]
        if settings['projection']:
            require(temperature in settings['projection']['temperatures'],
                    'temperature absent from saved projected Hall scan')
            preferences['character_temperature'] = temperature
    if group is not None:
        require(group in [row['name'] for row in settings['groups']['groups']], 'unknown saved character group')
        preferences['character_group'] = group
    if band is not None:
        require(band >= 1, 'map band must be positive')
        preferences['map_band'] = band
    if settings['projection']:
        import numpy as np

        with np.load(source / 'character' / 'character.npz', allow_pickle=False) as saved:
            require(preferences['map_band'] <= saved['energies_eV'].shape[1],
                    'map band exceeds the saved character band count')
    if quantity is not None:
        preferences['hall_quantity'] = quantity
        preferences['character_delta'] = quantity == 'delta-sigma'
    root = Path(output_dir).expanduser().resolve() if output_dir else source / 'figures'
    require(not root.exists(), 'figure directory exists; use --output-dir with a new directory')
    require((source / 'hall' / 'conductivity.csv').is_file(), 'saved Hall table is missing')
    root.mkdir(parents=True, exist_ok=False)
    record = dict(schema='vaspberry.postprocess-plots', version=1, complete=False, status='RUNNING',
                  source_run=str(source), source_run_sha256=sha256(source / 'run.json'),
                  preferences=preferences, manifest_name='plots.json', stages=[])
    write_json(root / 'plots.json', record)
    try:
        execute_stage(root, record, 'charge-hall-figures', python_command('plot_hall.py',
            source / 'hall' / 'conductivity.csv', '--regions', *preferences['hall_regions'],
            '--temperatures', *map(str, preferences['hall_temperatures']),
            '--quantity', preferences['hall_quantity'], '--output-dir', root / 'charge-hall'))
        if settings['projection']:
            delta = ['--delta'] if preferences['character_delta'] else []
            execute_stage(root, record, 'character-figures', python_command('procar_character.py',
                'plot', '--character-dir', source / 'character', '--hall-dir', source / 'character-hall',
                '--group', preferences['character_group'], '--band', preferences['map_band'],
                '--temperature', preferences['character_temperature'],
                '--region', preferences['character_region'], *delta,
                '--output-dir', root / 'character'))
        record.update(complete=True, status='PASS')
    except BaseException as exc:
        record.update(status='FAILED', error=str(exc))
        raise
    finally:
        write_json(root / 'plots.json', record)
    print('Saved PNG/PDF/SVG figures: ' + str(root))
    return record


def parser():
    result = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    result.add_argument('--version', action='version', version='VASPBERRY ' + (ROOT / 'VERSION').read_text().strip())
    sub = result.add_subparsers(dest='command', required=True)
    for name, help_text in [('check', 'validate settings and inputs without writing results; WAVEDER also evaluates the requested response'),
                            ('run', 'standard WAVEDER Hall; native WAVECAR approximation only with explicit settings opt-in')]:
        command = sub.add_parser(name, help=help_text, description=help_text)
        command.add_argument('settings', type=Path, help='commented INI settings file')
        command.add_argument('--reuse', type=Path, metavar='PREVIOUS_RUN',
                             help='reuse a pair cache only with [run] kubo_source=wavecar; preserves its approximation')
    plots = sub.add_parser('plot', help='draw saved numerical tables; do not execute VASP or VASPBERRY')
    plots.add_argument('result', type=Path, help='completed run directory')
    plots.add_argument('--output-dir', type=Path, help='new figure directory (default: RESULT/figures)')
    plots.add_argument('--temperature', type=float, help='draw one already calculated temperature (K)')
    plots.add_argument('--group', help='draw another saved atom/orbital group')
    plots.add_argument('--band', type=int, help='band for the character map; the saved Hall band sum is unchanged')
    plots.add_argument('--quantity', choices=['sigma', 'delta-sigma'], help='absolute or reference-subtracted curves')
    return result


def main(argv=None):
    cli = parser()
    args = cli.parse_args(argv)
    try:
        if args.command == 'plot':
            plot_results(args.result, args.output_dir, temperature=args.temperature,
                         group=args.group, band=args.band, quantity=args.quantity)
        else:
            settings = load_settings(args.settings)
            settings_summary(settings)
            if args.command == 'check':
                preflight(settings, args.reuse)
                print('Resolved spin mode: ' + settings['run']['resolved_spin_mode'])
                if settings['run']['kubo_source'] == 'waveder':
                    print('Settings, optical inputs and requested WAVEDER response checked; '
                          'no result files were written. This does not establish material convergence.')
                else:
                    print('Settings, source identity and mesh checked; no result files were written. '
                          'Matrix and response checks continue during run; material convergence is not established.')
            else:
                run_calculation(settings, args.reuse)
    except (ValueError, OSError, KeyError, TypeError) as exc:
        cli.error(str(exc))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
