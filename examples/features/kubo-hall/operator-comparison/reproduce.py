#!/usr/bin/env python3
"""Reintegrate the actual matched MoS2 pair caches without running VASP."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-dir', type=Path, default=HERE/'reference')
    p.add_argument('--output-dir', type=Path, required=True)
    a = p.parse_args()
    if a.output_dir.exists():
        p.error('output directory exists; choose a new directory')
    ref, out = a.reference_dir.resolve(), a.output_dir.resolve()
    info = json.loads((ref/'result.json').read_text())
    if info['status'] != 'PASS':
        p.error('completed reference required')
    out.mkdir(parents=True)
    record = dict(info, status='RUNNING', commands=[], numeric_reproduction=[])
    env = dict(os.environ)
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        env[name] = '1'

    def save():
        (out/'result.json').write_text(json.dumps(record, indent=2)+'\n')

    try:
        save()
        for operator in ('canonical', 'paw'):
            for cap in info['pair_band_max']:
                name = f'{operator}-cap{cap}'
                argv = [sys.executable, str(ROOT/'tools/vaspberry_kubo.py'), 'pair-hall',
                    '--pairs-dir', str(ref/f'{operator}-pairs'),
                    '--mu-min', str(info['vbm_eV']-.2), '--mu-max', str(info['vbm_eV']+.1),
                    '--mu-num', '61', '--mu-reference', str(info['mu_reference_eV']),
                    '--temperatures', '0', '300', '--regions', str(HERE.parent/'regions.json'),
                    '--difference', 'valley:K:Kprime', '--degeneracy-policy', 'coalesce',
                    '--degeneracy-threshold-eV', '1e-7', '--pair-band-max', str(cap),
                    '--formats', 'csv', 'dat', 'npz', '--output-dir', str(out/name)]
                with (out/f'{name}.stdout.log').open('w') as stdout, (out/f'{name}.stderr.log').open('w') as stderr:
                    done = subprocess.run(argv, env=env, cwd=ROOT, stdout=stdout, stderr=stderr)
                record['commands'].append(dict(name=name, argv=argv, returncode=done.returncode))
                save()
                if done.returncode:
                    raise RuntimeError(name+' failed; inspect its logs')
                for fmt in ('csv', 'dat', 'npz'):
                    file = f'conductivity.{fmt}'
                    same = (out/name/file).read_bytes() == (ref/name/file).read_bytes()
                    record['numeric_reproduction'].append(dict(file=f'{name}/{file}', byte_identical=same))
                    if not same:
                        raise RuntimeError(f'{name}/{file} differs from the reference')
        record['status'] = 'PASS'
        save()
    except BaseException as exc:
        record.update(status='FAILED', error=f'{type(exc).__name__}: {exc}')
        save()
        raise
    print(json.dumps(dict(status='PASS', output=str(out), byte_identical_files=len(record['numeric_reproduction']))))


if __name__ == '__main__':
    main()
