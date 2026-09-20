#!/usr/bin/env python3
"""Run current VASPBERRY Z2 on the actual public Bi WAVECAR, validate and plot."""
from __future__ import annotations

import csv
import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "examples/features/fukui-chern"))
from bi_common import arguments, execute, write_result, wavecar_geometry, sha
sys.path.insert(0, str(ROOT / "examples/Bi_Z2/scripts"))
from plot_nfield import read_field, validate_result, half_sums, integer_style
sys.path.insert(0, str(ROOT / "tools"))
from plot_berry_curvature import draw_plaquette_map, reciprocal_from_poscar

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


def make_figure(source, reciprocal, output):
    """Draw categorical n-field values on the physical reciprocal-cell geometry."""
    field, metadata, summary = validate_rows(Path(source))
    with Path(source).open() as handle:
        rows = list(csv.DictReader(line for line in handle if not line.startswith("#")))
    q = np.array([[float(row[key]) for key in ("q1", "q2", "q3")] for row in rows])
    cartesian = np.array([[float(row[key]) for key in ("kx_A-1", "ky_A-1", "kz_A-1")] for row in rows])
    # Legacy reciprocal-vector conversion differs from exact 2*pi by ~3e-8
    # relative; retain compatibility while rejecting a mismatched lattice.
    if not np.allclose(cartesian, q @ reciprocal, atol=1e-7, rtol=1e-7):
        raise ValueError("Z2 field coordinates disagree with the supplied lattice")
    values = np.array([int(row["nfield_int"]) for row in rows])
    cmap, norm, legend = integer_style(field)
    with plt.rc_context({"font.size": 11, "axes.titlesize": 13, "axes.linewidth": .8}):
        fig, ax = plt.subplots(figsize=(5.2, 4.5), layout="constrained")
        artist, info = draw_plaquette_map(ax, q, values, reciprocal,
                                        mesh=(summary["mesh_nx"], summary["mesh_ny"]), cmap=cmap, norm=norm)
        ax.set_title(r"Bi: $Z_2=" + str(summary["z2"]) + r"$", pad=11)
        colorbar = fig.colorbar(artist, ax=ax, shrink=.72, pad=.035, ticks=sorted(set(values)))
        colorbar.set_label(r"Integer field $n(\mathbf{k})$")
        fig.savefig(output, dpi=300)
        plt.close(fig)
    info.update(quantity="Fukui-Hatsugai integer n-field", source_field_sha256=sha(source),
                plotter_sha256=sha(ROOT / "tools/plot_berry_curvature.py"),
                note="Gauge-dependent integer plaquette field; not a local observable Berry curvature.")
    return info


def plot_existing():
    parser = argparse.ArgumentParser(description="Plot an existing Z2_FIELD.csv in Cartesian reciprocal space")
    parser.add_argument("--plot-only", type=Path, required=True, help="existing reportable Z2_FIELD.csv")
    parser.add_argument("--poscar", type=Path, required=True, help="matching lattice POSCAR")
    parser.add_argument("--figure", type=Path, required=True, help="new PNG/PDF/SVG figure")
    args = parser.parse_args()
    record = args.figure.with_suffix(args.figure.suffix + ".json")
    if args.figure.exists() or record.exists():
        parser.error("figure output exists; choose a new filename")
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    try:
        info = make_figure(args.plot_only, reciprocal_from_poscar(args.poscar), args.figure)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    info.update(workflow_mode="plot_existing_native_output", poscar_sha256=sha(args.poscar),
                figure_sha256=sha(args.figure))
    record.write_text(json.dumps(info, indent=2) + "\n")
    print(args.figure)


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
    _, reciprocal = wavecar_geometry(args.wavecar)
    provenance["figure"] = make_figure(out / "Z2_FIELD.csv", reciprocal, out / "figure.png")
    make_figure(out / "Z2_FIELD.csv", reciprocal, out / "figure.pdf")
    write_result(out, "z2", provenance, summary, ["NFIELD.dat", "Z2_FIELD.csv"], [
        "Real archived VASP input; reproduces VASPBERRY postprocessing, not an end-to-end new VASP calculation.",
        "The pointwise n-field is gauge and branch dependent; the agreed half-zone parity is the reported invariant.",
        "Numerical checks establish reconstruction consistency, not independent raw-input TR symmetry or mesh convergence.",
        "WAVECAR pseudo-wavefunction overlaps omit PAW augmentation; confirm the gap, symmetry and convergence for another system."])


if __name__ == "__main__":
    plot_existing() if "--plot-only" in sys.argv else main()
