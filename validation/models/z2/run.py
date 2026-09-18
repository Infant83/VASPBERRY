#!/usr/bin/env python3
"""Validate/plot stored Bi schema-2 data; optionally rerun WAVECAR postprocessing."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/'examples/Bi_Z2/scripts'))
from plot_nfield import read_field, validate_result, half_sums, integer_style, draw_field


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def validate_rows(path):
    field, metadata = read_field(path)
    nx, ny, z2 = validate_result(field, metadata, path)
    if metadata['schema_version'] != '2':
        raise ValueError('this fixture requires schema version 2')
    with path.open() as f:
        raw = list(csv.DictReader(line for line in f if not line.startswith('#')))
    rows = {int(row['cell_id']): {key: float(value) for key, value in row.items()} for row in raw}
    if len(raw) != len(rows) or set(rows) != set(range(1, nx*ny+1)):
        raise ValueError('duplicate/missing cell IDs')
    values = np.array([list(row.values()) for row in rows.values()])
    if not np.isfinite(values).all():
        raise ValueError('nonfinite stored row')
    for cell, row in rows.items():
        partner_id = int(row['tr_partner'])
        if partner_id not in rows or row['tr_partner'] != partner_id:
            raise ValueError('invalid TR partner')
        partner = rows[partner_id]
        if partner['tr_partner'] != cell or row['plaquette_checks_pass'] != 1:
            raise ValueError('failed row or non-involutive TR map')
        if row['nfield_int'] != round(row['nfield_int']):
            raise ValueError('noninteger n-field')
        if row['pair_nfield_int_sum'] != row['nfield_int']+partner['nfield_int']:
            raise ValueError('inconsistent pair sum')
        q = np.array([row['q1'], row['q2'], row['q3']])
        qp = np.array([partner['q1'], partner['q2'], partner['q3']])
        if abs(q+qp-np.rint(q+qp)).max() > float(metadata['threshold_partner_fractional']):
            raise ValueError('TR partner coordinates are inconsistent')
    flux = np.array([rows[i]['berry_flux_rad'] for i in sorted(rows)])
    recomputed = {
        'total_chern': float(flux.sum()/(2*np.pi)),
        'minimum_link_singular_value': min(row['min_link_singular_value'] for row in rows.values()),
        'max_nfield_integer_residual': max(abs(row['nfield_raw']-row['nfield_int']) for row in rows.values()),
        'max_flux_tr_odd_residual_rad': max(abs(np.angle(np.exp(1j*(row['berry_flux_rad']+rows[int(row['tr_partner'])]['berry_flux_rad'])))) for row in rows.values()),
        'max_abs_flux_rad': float(abs(flux).max()),
    }
    for name, value in recomputed.items():
        if not np.isclose(value, float(metadata[name]), atol=1e-12, rtol=1e-8):
            raise ValueError('stored diagnostic does not match rows: '+name)
    if (abs(recomputed['total_chern']) > float(metadata['threshold_total_chern'])
            or recomputed['minimum_link_singular_value'] < float(metadata['threshold_min_link_singular'])
            or recomputed['max_nfield_integer_residual'] > float(metadata['threshold_nfield_integer'])
            or recomputed['max_flux_tr_odd_residual_rad'] > float(metadata['threshold_flux_tr_odd_rad'])
            or recomputed['max_abs_flux_rad'] >= float(metadata['threshold_max_abs_flux_rad'])):
        raise ValueError('recomputed row diagnostics fail producer thresholds')
    return field, metadata, {'mesh_nx': nx, 'mesh_ny': ny, 'plaquettes': len(rows),
                            'z2': z2, 'half_top_sum': half_sums(field)[0],
                            'half_bottom_sum': half_sums(field)[1], **recomputed}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir', type=Path, required=True, help='fresh result directory')
    p.add_argument('--wavecar', type=Path, help='optional actual tracked Bi WAVECAR; never an LFS pointer')
    p.add_argument('--binary', type=Path, default=ROOT/'build/vaspberry-gfortran', help='current executable, used only with --wavecar')
    p.add_argument('--timeout', type=float, default=300, help='optional postprocessing time limit in seconds')
    args = p.parse_args()
    if args.output_dir.exists():
        p.error('output directory exists; use a fresh path')
    config = json.loads((HERE/'input.json').read_text())
    reference = ROOT/config['reference_field']
    if sha(reference) != config['reference_sha256']:
        p.error('stored reference hash changed; do not silently replace the historical fixture')
    _, _, reference_summary = validate_rows(reference)
    if (reference_summary['z2'] != config['expected_z2'] or
            [reference_summary['half_top_sum'], reference_summary['half_bottom_sum']] != config['expected_half_sums']):
        p.error('reference invariants disagree with fixture input')
    provenance = {'input_sha256': sha(HERE/'input.json'), 'runner_sha256': sha(__file__),
                  'reference_field': config['reference_field'], 'reference_sha256': sha(reference),
                  'plot_helper_sha256': sha(ROOT/'examples/Bi_Z2/scripts/plot_nfield.py')}
    if args.wavecar is not None:
        with args.wavecar.open('rb') as f:
            if f.read(128).startswith(b'version https://git-lfs.github.com/spec/'):
                p.error('WAVECAR is a Git LFS pointer; fetch the actual tracked Bi payload first')
        if args.wavecar.stat().st_size != config['wavecar_bytes'] or sha(args.wavecar) != config['wavecar_sha256']:
            p.error('this optional reproduction requires the exact public Bi WAVECAR fixture')
        if not args.binary.is_file() or not np.isfinite(args.timeout) or args.timeout <= 0:
            p.error('provide a built executable and a finite positive timeout')
        provenance.update(wavecar_sha256=sha(args.wavecar), binary_sha256=sha(args.binary))
    args.output_dir.mkdir(parents=True)
    source = reference
    mode = 'stored_field_validation_and_new_plot'
    if args.wavecar is not None:
        command = [str(args.binary.resolve()), '-f', str(args.wavecar.resolve()), '-o', 'NFIELD',
                   '-z2', '1', '-kx', '12', '-ky', '12', '-s', '2', '-ii', '1', '-if', '10']
        try:
            with (args.output_dir/'fortran.log').open('w') as out:
                run = subprocess.run(command, cwd=args.output_dir, stdout=out, stderr=subprocess.STDOUT,
                                     timeout=args.timeout)
            if run.returncode != 0:
                raise RuntimeError('postprocessing executable failed with exit '+str(run.returncode))
        except (OSError, RuntimeError, subprocess.TimeoutExpired) as error:
            (args.output_dir/'result.json').write_text(json.dumps({'schema_version': 1, 'feature_id': 'z2',
                'status': 'FAILED', 'workflow_mode': 'wavecar_postprocessing',
                'mode': 'wavecar_postprocessing', 'error': str(error), 'provenance': provenance}, indent=2)+'\n')
            p.exit(1, str(error)+'\n')
        source = args.output_dir/'Z2_FIELD.csv'
        mode = 'new_wavecar_postprocessing_compared_to_stored_invariant'
    field, meta, summary = validate_rows(source)
    if (summary['mesh_nx'], summary['mesh_ny']) != tuple(config['mesh']) or any(
            int(meta[key]) != value for key, value in (('band_min', 1), ('band_max', 10),
                                                       ('band_rank', 10), ('spinor_components', 2))):
        raise ValueError('field does not describe the required Bi mesh and occupied subspace')
    if args.wavecar is not None and meta['vaspberry_version'] != (ROOT/'VERSION').read_text().strip():
        raise ValueError('rerun executable version does not match this checkout; rebuild the binary')
    if summary['z2'] != reference_summary['z2']:
        raise ValueError('recomputed Z2 disagrees with the historical reference')
    with (args.output_dir/'summary.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(summary)); writer.writeheader(); writer.writerow(summary)
    cmap, norm, legend = integer_style(field)
    fig, ax = plt.subplots(figsize=(6.4, 6.2))
    draw_field(ax, field, 'Bi reference: stored field revalidated' if args.wavecar is None else 'Bi: new WAVECAR postprocessing',
               cmap, norm, (summary['half_top_sum'], summary['half_bottom_sum']))
    fig.subplots_adjust(left=.14, right=.95, bottom=.19, top=.79)
    fig.legend(handles=legend, loc='lower center', bbox_to_anchor=(.5, .065), ncol=3, frameon=False)
    fig.suptitle('Fukui-Hatsugai Z2 = '+str(summary['z2']), y=.96, fontsize=16)
    fig.text(.5, .015, 'Pointwise n-field is gauge dependent; half-zone parity is the invariant.\n'
             'Default mode does not rerun WAVECAR or DFT.', ha='center', fontsize=9)
    fig.savefig(args.output_dir/'figure.png', dpi=180, bbox_inches='tight'); plt.close(fig)
    result = {'schema_version': 1, 'feature_id': 'z2', 'status': 'PASS', 'workflow_mode': mode, 'mode': mode,
              'software_version': (ROOT/'VERSION').read_text().strip(),
              'field_producer_version': meta['vaspberry_version'], 'summary': summary,
              'numerical_checks': {'status': 'PASS', 'reference_hash': 'PASS', 'schema_and_half_zone_parities': 'PASS',
                                   'row_partner_consistency': 'PASS', 'recomputed_row_diagnostics': 'PASS',
                                   'reference_z2_agreement': 'PASS'},
              'reference_z2_agrees': True, 'source_field_sha256': sha(source), 'provenance': provenance,
              'limitations': ['Default revalidates stored schema-2 rows and regenerates a plot; it does not recalculate wavefunctions or links.',
                              'Optional WAVECAR mode recomputes postprocessing only, not VASP SCF/NSCF.',
                              'Numerical TR reconstruction checks do not independently establish raw-input TR symmetry, gap, PAW completeness or mesh convergence.'],
              'output_sha256': {name: sha(args.output_dir/name) for name in ('summary.csv', 'figure.png')}}
    (args.output_dir/'result.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': 'PASS', 'mode': mode, 'z2': summary['z2']}))


if __name__ == '__main__':
    main()
