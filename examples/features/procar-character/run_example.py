#!/usr/bin/env python3
"""Run public projection CLI stages and check independent analytic values."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[3]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def run(directory):
    out = Path(directory).resolve()
    if out.exists(): raise ValueError('choose a fresh output directory')
    out.mkdir(parents=True)
    log = {'schema': 'vaspberry.procar-example-run', 'version': 1, 'status': 'RUNNING',
           'fixture': 'analytic assigned numerator, not a VASP material result', 'commands': []}
    manifest = out/'run.json'
    def save(): manifest.write_text(json.dumps(log, indent=2, allow_nan=False)+'\n')
    def execute(label, argv):
        command = [sys.executable, *map(str, argv)]
        record = dict(label=label, argv=command, cwd=str(ROOT)); log['commands'].append(record); save()
        start = time.monotonic()
        with (out/(label+'.stdout.log')).open('w') as stdout, (out/(label+'.stderr.log')).open('w') as stderr:
            result = subprocess.run(command, cwd=ROOT, stdout=stdout, stderr=stderr)
        record.update(exit_code=result.returncode, elapsed_s=time.monotonic()-start); save()
        if result.returncode: raise RuntimeError(label+' failed; inspect its stderr/stdout logs')
    try:
        inp, character, pairs, hall = [out/name for name in ('input', 'character', 'pairs', 'hall')]
        execute('fixture', ['examples/features/procar-character/make_fixture.py', '--output-dir', inp])
        execute('project', ['tools/procar_character.py', 'project', '--procar', inp/'PROCAR',
            '--wavecar', inp/'WAVECAR', '--outcar', inp/'OUTCAR', '--groups', inp/'groups.json',
            '--axis', 1, 0, 0, '--output-dir', character])
        execute('pairs', ['tools/vaspberry_kubo.py', 'import-pairs', '--csv', inp/'PAIRS.csv',
            '--wavecar', inp/'WAVECAR', '--spinor-components', 2, '--spin-multiplicity', 1,
            '--mesh', 2, 2, '--energy-reference', 'synthetic analytic fixture zero', '--output-dir', pairs])
        execute('hall', ['tools/procar_character.py', 'hall', '--character-dir', character,
            '--pairs-dir', pairs, '--bands', 1, '--mu-min', -1.2, '--mu-max', .2, '--mu-num', 71,
            '--mu-reference', 0, '--temperatures', 0, 300, '--output-dir', hall])
        execute('plot', ['tools/procar_character.py', 'plot', '--character-dir', character,
            '--hall-dir', hall, '--group', 'lower', '--band', 1, '--temperature', 300,
            '--output-dir', out/'figures'])
        with np.load(hall/'character_hall.npz', allow_pickle=False) as saved:
            a = {name: saved[name] for name in saved.files}
        errors = []
        expected = {('lower', 'charge'): -.4*np.pi, ('lower', 'plus'): -.4*np.pi,
                    ('lower', 'minus'): 0., ('upper', 'charge'): -1.2*np.pi,
                    ('upper', 'plus'): 0., ('upper', 'minus'): -1.2*np.pi,
                    ('$unweighted', 'charge'): -2*np.pi, ('$all_projected', 'charge'): -1.6*np.pi,
                    ('$unprojected_residual', 'charge'): -.4*np.pi}
        values = {}
        for (group, component), target in expected.items():
            mask = ((a['group'] == group) & (a['component'] == component) & (a['region'] == 'total')
                    & (a['temperature_K'] == 0) & np.isclose(a['mu_eV'], 0, atol=1e-14, rtol=0))
            if mask.sum() != 1: raise AssertionError('missing/duplicate analytic evaluation row')
            value = float(a['attribution_e2_over_h'][mask][0])
            errors.append(abs(value-target)); values[group+'/'+component] = value
        if max(errors) > 1e-12: raise AssertionError('analytic Hall values disagree')
        with np.load(character/'character.npz', allow_pickle=False) as saved:
            c = saved['characters']
            joint_error = float(np.max(abs(c[..., 2]+c[..., 3]-c[..., 0])))
            q_range = [float(saved['all_ions_cartesian'][..., 0].min()), float(saved['all_ions_cartesian'][..., 0].max())]
        if joint_error > 1e-14: raise AssertionError('joint spin closure failed')
        checks = dict(status='PASS', analytic_max_error_e2_over_h=max(errors),
                      joint_spin_max_error=joint_error, all_ions_charge_range=q_range,
                      hall_rows=len(a['mu_eV']), expected_T0_mu0=values)
        (out/'verification.json').write_text(json.dumps(checks, indent=2)+'\n')
        log.update(status='PASS', complete=True, checks=checks)
    except BaseException as exc:
        log.update(status='FAILED', complete=False, error=str(exc)); save(); raise
    log['artifact_sha256'] = {str(path.relative_to(out)): sha(path) for path in sorted(out.rglob('*'))
                              if path.is_file() and path != manifest}
    log['implementation_sha256'] = {str(path.relative_to(ROOT)): sha(path) for path in
        [Path(__file__), Path(__file__).with_name('make_fixture.py'), ROOT/'tools/procar_character.py', ROOT/'tools/procar_projection.py']}
    save()
    print(json.dumps(checks, indent=2))
    return checks


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    run(parser.parse_args().output_dir)
