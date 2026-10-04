#!/usr/bin/env python3
"""Exercise selected WAVEDER CLI/INI/plots with declared synthetic inputs.

Checks software arithmetic and rejection barriers, not material convergence.
Every subprocess command, output and exit status is retained in --output-dir.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--binary', type=Path,
                        help='optional current native binary, enabling geometric checks')
    args = parser.parse_args()
    out = args.output_dir.resolve()
    if out.exists():
        parser.error('output directory exists; use a new directory')
    out.mkdir(parents=True)
    steps = []
    def run(name, command, rejection=None):
        command = [str(x) for x in command]
        result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
        (out/f'{name}.stdout.log').write_text(result.stdout)
        (out/f'{name}.stderr.log').write_text(result.stderr)
        entry = dict(name=name, command=command, cwd=str(ROOT),
                     exit_code=result.returncode, expected_rejection=rejection)
        steps.append(entry)
        (out/'commands.json').write_text(json.dumps(steps, indent=2)+'\n')
        if rejection is None:
            if result.returncode:
                raise RuntimeError(f'{name} failed; see saved logs')
        elif result.returncode == 0 or rejection not in result.stdout+result.stderr:
            raise RuntimeError(f'{name} did not produce expected rejection: {rejection}')

    result = {'status': 'FAIL', 'kind': 'synthetic_software_example',
              'not_material_data': True}
    try:
        fixture = out/'input'
        run('fixture', [sys.executable, Path(__file__).with_name('make_fixture.py'),
                        '--output-dir', fixture])
        post = [sys.executable, ROOT/'tools/vaspberry_post.py']
        run('ini-check', [*post, 'check', fixture/'analysis.ini'])
        run('ini-run', [*post, 'run', fixture/'analysis.ini'])
        hall = fixture/'hall-selected'
        run('ini-plot', [*post, 'plot', hall])
        expected = json.loads((fixture/'expected.json').read_text())
        omega = np.asarray(expected['omega_by_k_band_A2'])[:, 2:4]
        energies = np.asarray(expected['energies_eV'])[2:4]
        with (hall/'hall/conductivity.csv').open() as stream:
            rows = [r for r in csv.DictReader(stream) if r['region'] == 'total']
        residuals = []
        def direct(mu, temperature):
            if temperature == 0:
                occupation = (energies <= mu).astype(float)
            else:
                x = (energies-mu)/(8.617333262145e-5*temperature)
                occupation = np.exp(-np.logaddexp(0., x))
            return -float(np.mean(omega @ occupation))/(2*math.pi)
        for row in rows:
            mu, temperature = float(row['mu_eV']), float(row['temperature_K'])
            sigma = direct(mu, temperature)
            residuals += [abs(float(row['sigma_e2_over_h'])-sigma),
                          abs(float(row['delta_sigma_e2_over_h'])-
                              (sigma-direct(0., temperature)))]
        if len(rows) != 27 or max(residuals) > 1e-12:
            raise AssertionError(f'weighted direct-sum oracle failed: {len(rows)}, {max(residuals)}')
        meta = json.loads((hall/'hall/conductivity.json').read_text())
        assert meta['scope'] == 'selected_band_contribution'
        assert meta['selected_band_ids'] == [3, 4]
        assert meta['intermediate_band_ids'] == [1, 2, 3, 4]
        assert meta['total_AHC_certified'] is False
        assert meta['total_delta_AHC_certified'] is False
        plotmeta = json.loads((hall/'figures/charge-hall/plot.json').read_text())['selected_contribution']
        assert plotmeta['scope'] == 'selected_band_contribution'
        assert plotmeta['selected_band_ids'] == [3, 4]
        run('missing-weighted-pair', [sys.executable, ROOT/'tools/vaspberry_kubo.py',
            'kubo-hall', '--input-dir', fixture/'rectangular', '--bands', '3:4',
            '--mesh', '2', '2', '--energy-reference', 'synthetic fixture energy zero',
            '--mu-min', '0', '--mu-max', '4', '--mu-num', '9', '--mu-reference', '0',
            '--temperatures', '300', '--output-dir', out/'rejected-weighted'],
            'missing rectangular matrix pair')
        result.update(weighted_rows=len(rows), max_weighted_oracle_error=max(residuals))
        if args.binary:
            binary = args.binary.resolve()
            def native(name, source, selection, *extra, rejection=None):
                path = out/f'{name}.csv'
                run(name, [binary, '--task', 'kubo', '--input-dir', fixture/source,
                           '--bands', selection, '--curvature-csv', path, *extra], rejection)
                if rejection:
                    assert not path.exists()
                    return None
                with path.open() as stream:
                    return list(csv.DictReader(line for line in stream if not line.startswith('#')))
            for name, source in [('native-full-trace', 'full'),
                                  ('native-rectangle-trace', 'rectangular')]:
                trace = native(name, source, '3:4', '--mesh', '2,2')
                np.testing.assert_allclose([float(r['omega_z_A2']) for r in trace],
                                           expected['selected_3_4_trace_A2'], atol=1e-14, rtol=0.)
            per = native('native-per-band', 'full', '3:4', '--per-band', '1')
            for row in per:
                np.testing.assert_allclose(float(row['omega_z_A2']),
                    omega[int(row['k_index'])-1, int(row['band'])-3], atol=1e-14, rtol=0.)
            native('missing-geometric-pair', 'rectangular', '3',
                   rejection='required selected-band pair is missing')
            native('split-producer-cluster', 'full', '1',
                   rejection='selected trace splits a full-source transitive producer cluster')
            result['native_checks'] = 'PASS'
            result['binary_sha256'] = hashlib.sha256(binary.read_bytes()).hexdigest()
        else:
            result['native_checks'] = 'NOT_RUN: provide --binary'
        result['status'] = 'PASS'
    except Exception as error:
        result['error'] = str(error)
        raise
    finally:
        (out/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
