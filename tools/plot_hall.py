#!/usr/bin/env python3
"""Plot standardized VASPBERRY Hall tables from CSV, DAT or NPZ."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from berry_data import require
from exported_matrix_kubo import sha256


def read_table(path):
    path = Path(path)
    if path.suffix == '.npz':
        with np.load(path, allow_pickle=False) as data:
            arrays = {key: data[key] for key in data.files}
    else:
        require(path.suffix in ('.csv', '.dat'), 'input must be CSV, DAT or NPZ')
        with path.open() as handle:
            lines = iter(handle)
            first = next(lines)
            if first.startswith('# '):
                first = first[2:]
            rows = list(csv.DictReader([first, *lines], delimiter='\t' if path.suffix == '.dat' else ','))
        require(bool(rows), 'empty Hall table')
        arrays = {key: np.array([row[key] for row in rows], dtype=str if key == 'region' else float)
                  for key in rows[0]}
    required = ('mu_eV', 'mu_minus_reference_eV', 'temperature_K', 'region', 'band_id',
                'sigma_e2_over_h', 'delta_sigma_e2_over_h')
    require(all(k in arrays for k in required), 'Hall table lacks required columns')
    n = len(arrays['region'])
    require(n > 0 and all(arrays[k].shape == (n,) for k in required), 'inconsistent Hall table dimensions')
    require(all(np.isfinite(arrays[k]).all() for k in required if k != 'region'), 'nonfinite Hall values')
    return arrays


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input', type=Path)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--formats', nargs='+', choices=['png', 'pdf', 'svg'], default=['png', 'pdf', 'svg'])
    p.add_argument('--regions', nargs='+', default=['total'])
    p.add_argument('--temperatures', nargs='+', type=float)
    p.add_argument('--quantity', choices=['sigma', 'delta-sigma'], default='sigma')
    p.add_argument('--energy-origin-eV', type=float, help='plot mu minus this energy; default is table mu_reference')
    p.add_argument('--energy-label', default=None, help='x-axis label, e.g. mu - VBM (eV)')
    args = p.parse_args(argv)
    try:
        require(not args.output_dir.exists(), 'output directory exists; choose a new directory')
        require(len(set(args.formats)) == len(args.formats), 'duplicate figure format')
        d = read_table(args.input)
        ts = args.temperatures if args.temperatures is not None else np.unique(d['temperature_K'])
        require(len(ts) and len(set(ts)) == len(ts) and np.isfinite(ts).all(), 'distinct finite temperatures required')
        require(len(set(args.regions)) == len(args.regions), 'duplicate region')
        col = 'sigma_e2_over_h' if args.quantity == 'sigma' else 'delta_sigma_e2_over_h'
        if args.energy_origin_eV is None:
            x = d['mu_minus_reference_eV']; xlabel = r'$\mu-\mu_{\rm ref}$ (eV)'
        else:
            require(np.isfinite(args.energy_origin_eV), 'finite energy origin required')
            x = d['mu_eV']-args.energy_origin_eV; xlabel = r'$\mu-E_0$ (eV)'
        curves = []
        for ri, region in enumerate(args.regions):
            for ti, t in enumerate(ts):
                mask = (d['region'] == region) & (d['temperature_K'] == t) & (d['band_id'] == 0)
                require(mask.any(), f'no total-subspace rows for {region}, T={t:g} K')
                order = np.argsort(x[mask])
                require(len(np.unique(x[mask])) == mask.sum(), 'duplicate energies within curve')
                curves.append((x[mask][order], d[col][mask][order], ri, ti, region, t))
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6.4, 4.2), layout='constrained')
        try:
            styles = ['-', '--', ':', '-.']
            for xs, ys, ri, ti, region, t in curves:
                ax.plot(xs, ys, color=f'C{ri%10}', linestyle=styles[ti%4],
                        label=f'{region}, {t:g} K', linewidth=1.7)
            ax.axhline(0, color='0.7', linewidth=.6, zorder=0)
            ax.set(xlabel=args.energy_label or xlabel,
                   ylabel=(r'$\sigma_{xy}$' if args.quantity == 'sigma' else r'$\Delta\sigma_{xy}$')+r' ($e^2/h$)')
            ax.legend(frameon=False); ax.tick_params(direction='in', top=True, right=True)
            args.output_dir.mkdir(parents=True)
            for fmt in args.formats:
                fig.savefig(args.output_dir/f'hall.{fmt}', dpi=220)
        finally:
            plt.close(fig)
        (args.output_dir/'plot.json').write_text(json.dumps(dict(
            source_sha256=sha256(args.input), command={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
            output_sha256={f'hall.{fmt}': sha256(args.output_dir/f'hall.{fmt}') for fmt in args.formats}), indent=2)+'\n')
    except (ValueError, OSError, KeyError, TypeError, StopIteration) as exc:
        p.error(str(exc))


if __name__ == '__main__':
    main()
