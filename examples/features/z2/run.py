#!/usr/bin/env python3
"""Run current VASPBERRY Z2 on the actual public Bi WAVECAR, validate and plot."""
from __future__ import annotations

import csv
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "examples/features/fukui-chern"))
from bi_common import arguments, execute, write_result
sys.path.insert(0, str(ROOT / "examples/Bi_Z2/scripts"))
from plot_nfield import read_field, validate_result, half_sums, integer_style, draw_field

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
    args = arguments(__doc__)
    out, provenance, edges, gaps = execute(args, "z2", ["-o", "NFIELD", "-z2", "1"])
    field, metadata, summary = validate_rows(out / "Z2_FIELD.csv")
    if (summary["mesh_nx"], summary["mesh_ny"], summary["z2"]) != (12, 12, 1):
        raise ValueError("Bi reference requires a 12 x 12 field and Z2 = 1")
    for key, expected in (("band_min", 1), ("band_max", 10), ("band_rank", 10), ("spinor_components", 2)):
        if int(metadata[key]) != expected:
            raise ValueError("unexpected selected occupied bundle: " + key)
    if [summary["half_top_sum"], summary["half_bottom_sum"]] != [-3, 3]:
        raise ValueError("public input n-field differs from its reference; inspect the field and numerical diagnostics")
    summary.update(gaps)
    cmap, norm, legend = integer_style(field)
    fig, ax = plt.subplots(figsize=(6.4, 6.2))
    draw_field(ax, field, "Bi: actual VASPBERRY calculation from WAVECAR", cmap, norm,
               (summary["half_top_sum"], summary["half_bottom_sum"]))
    fig.subplots_adjust(left=.14, right=.95, bottom=.19, top=.79)
    fig.legend(handles=legend, loc="lower center", bbox_to_anchor=(.5, .065), ncol=3, frameon=False)
    fig.suptitle("Bi 12 × 12 SOC · Fukui-Hatsugai Z2 = 1", y=.96, fontsize=15)
    fig.text(.5, .015, "Pointwise n-field is gauge dependent; the half-zone parity is the invariant.\n"
             "Fresh post-processing of the public VASP 5.4.1 spinor WAVECAR.", ha="center", fontsize=9)
    fig.savefig(out / "figure.png", dpi=180, bbox_inches="tight"); plt.close(fig)
    write_result(out, "z2", provenance, summary, ["NFIELD.dat", "Z2_FIELD.csv"], [
        "Real archived VASP input; reproduces VASPBERRY postprocessing, not an end-to-end new VASP calculation.",
        "The pointwise n-field is gauge and branch dependent; the agreed half-zone parity is the reported invariant.",
        "Numerical checks establish reconstruction consistency, not independent raw-input TR symmetry or mesh convergence.",
        "WAVECAR pseudo-wavefunction overlaps omit PAW augmentation; confirm the gap, symmetry and convergence for another system."])


if __name__ == "__main__":
    main()
