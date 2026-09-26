#!/usr/bin/env python3
"""Execute VASPBERRY and postprocess its results from one settings file.

  python3 tools/vaspberry_post.py check analysis.ini
  python3 tools/vaspberry_post.py run analysis.ini
  python3 tools/vaspberry_post.py plot results/run01

The run command executes VASPBERRY (with MPI when configured), then uses
Python for numerical Hall integration and optional PROCAR analysis. The plot
command reads saved tables. This front end records the underlying commands
and keeps the numerical data available for independent analysis.
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
    wavecar = Path(run['wavecar'])
    require(wavecar.is_file(), 'WAVECAR does not exist: ' + str(wavecar))
    wave = Wavecar(wavecar, spin=run['spin'], spinor_components=run['spinor_components'])
    require(wave.header.ispin == run['expected_nspin'],
            '[run] spin_mode does not match the WAVECAR ISPIN header')
    require(len(wave.kpoints) == run['mesh'][0] * run['mesh'][1],
            '[run] mesh size does not match the number of WAVECAR k points')
    axes = run['plane_axes']
    other = next(i for i in range(3) if i not in axes)
    infer_uniform_grid(wave.kpoints[:, axes + [other]], *run['mesh'])
    wave.coefficients(0, [1])
    projection = settings['projection']
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
        for key in ('mesh', 'plane_axes', 'spin', 'spinor_components',
                    'spin_multiplicity', 'energy_reference'):
            require(old_run[key] == run[key], '--reuse cannot change source setting: ' + key)
        verify_saved(previous, record, ['pairs/pairs.npz', 'pairs/pairs.json'])
        read_pairs(previous / 'pairs')
        cache = previous / 'pairs'
        inputs['reused_run'] = {'path': str(previous), 'sha256': sha256(previous / 'run.json')}
    return inputs, cache


def settings_summary(settings):
    run, hall = settings['run'], settings['hall']
    print('WAVECAR: ' + run['wavecar'])
    print('Output:  ' + run['output'])
    print(f"Source: {run['spin_mode']}; full {run['mesh'][0]} x {run['mesh'][1]} mesh")
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
        if cache is None:
            raw = root / 'native'
            raw.mkdir()
            command = [run['binary'], '--task', 'kubo-pairs', '--wavecar', run['wavecar'],
                       '--spinor', str(run['spinor_components']), '--pairs-csv', 'PAIRS.csv']
            if run['mpi_procs'] > 1:
                command = [shutil.which(run['mpi_launcher']), '-n', str(run['mpi_procs']), *command]
            execute_stage(root, record, 'native-pairs', command, cwd=raw, native=True)
            execute_stage(root, record, 'import-pairs', python_command('vaspberry_kubo.py',
                'import-pairs', '--csv', raw / 'PAIRS.csv', *source_arguments(run),
                '--output-dir', root / 'pairs'))
        else:
            shutil.copytree(cache, root / 'pairs')
            record['reused_pair_cache'] = str(cache)
            print('Reusing saved pairs; VASPBERRY execution is skipped.', flush=True)
        region_args = ['--regions', str(root / 'regions.json')] if settings['regions']['regions'] else []
        difference_args = [item for row in settings['differences'] for item in ('--difference', ':'.join(row))]
        cap = settings['hall']['pair_band_max']
        cap_args = [] if cap is None else ['--pair-band-max', str(cap)]
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
    for name, help_text in [('check', 'check settings and source paths without running calculations'),
                            ('run', 'execute VASPBERRY, then perform the requested numerical postprocessing')]:
        command = sub.add_parser(name, help=help_text)
        command.add_argument('settings', type=Path, help='commented INI settings file')
        command.add_argument('--reuse', type=Path, metavar='PREVIOUS_RUN',
                             help='reuse its validated pair cache for the same WAVECAR; skip VASPBERRY execution')
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
                print('Settings and source paths checked. Numerical validity is checked during calculation.')
            else:
                run_calculation(settings, args.reuse)
    except (ValueError, OSError, KeyError, TypeError) as exc:
        cli.error(str(exc))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
