#!/usr/bin/env python3
"""Compare existing reportable native Z2 integer fields; plotting only."""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_spec = importlib.util.spec_from_file_location("vaspberry_z2_example", Path(__file__).with_name("run.py"))
_example = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_example)
ROOT, grid, integer_style = _example.ROOT, _example.grid, _example.integer_style
reciprocal_from_poscar, sha, validate_rows = _example.reciprocal_from_poscar, _example.sha, _example.validate_rows


def portable(path):
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return path.name


def load_panel(field_path, poscar, label):
    field, metadata, summary = validate_rows(field_path)
    with field_path.open() as handle:
        rows = list(csv.DictReader(line for line in handle if not line.startswith("#")))
    q = np.array([[float(r[k]) for k in ("q1", "q2", "q3")] for r in rows])
    cart = np.array([[float(r[k]) for k in ("kx_A-1", "ky_A-1", "kz_A-1")] for r in rows])
    if not np.allclose(cart, q @ reciprocal_from_poscar(poscar), atol=1e-7, rtol=1e-7):
        raise ValueError(f"{label}: field and POSCAR lattices disagree")
    xs, ys, values = grid(field)
    nx, ny = summary["mesh_nx"], summary["mesh_ny"]
    xedges, yedges = np.linspace(-.5, .5, nx+1), np.linspace(-.5, .5, ny+1)
    if (not np.allclose(xs, (xedges[:-1]+xedges[1:])/2, atol=1e-6, rtol=0)
            or not np.allclose(ys, (yedges[:-1]+yedges[1:])/2, atol=1e-6, rtol=0)):
        raise ValueError(f"{label}: an unshifted uniform native plaquette mesh is required")
    info = dict(label=label, source_field=portable(field_path), source_field_sha256=sha(field_path),
                poscar=portable(poscar), poscar_sha256=sha(poscar),
                occupied_bands=[int(metadata["band_min"]), int(metadata["band_max"])],
                half_top_parity=abs(summary["half_top_sum"]) % 2,
                half_bottom_parity=abs(summary["half_bottom_sum"]) % 2, **summary)
    return dict(values=values, xedges=xedges, yedges=yedges, info=info)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--field", type=Path, action="append", required=True, help="repeat for each existing Z2_FIELD.csv")
    p.add_argument("--poscar", type=Path, action="append", required=True, help="repeat in field order")
    p.add_argument("--label", action="append", required=True, help="repeat in field order")
    p.add_argument("--output-dir", type=Path, required=True, help="fresh plot-only output directory")
    p.add_argument("--formats", nargs="+", choices=["png", "svg", "pdf"], default=["png", "svg"])
    a = p.parse_args()
    if not 1 <= len(a.field) <= 4 or len(a.field) != len(a.poscar) or len(a.field) != len(a.label):
        p.error("supply one to four matching fields, POSCARs and labels")
    if a.output_dir.exists() or len(a.formats) != len(set(a.formats)):
        p.error("choose a fresh output directory and distinct formats")
    panels = [load_panel(f, s, label) for f, s, label in zip(a.field, a.poscar, a.label)]
    palette = {(i, 0): int(v) for i, v in enumerate(np.unique(np.concatenate([p["values"].ravel() for p in panels])))}
    cmap, norm, legend = integer_style(palette)
    with plt.rc_context({"font.size": 11, "axes.titlesize": 13, "axes.linewidth": .8}):
        fig, axes = plt.subplots(1, len(panels), figsize=(4.4*len(panels), 4.7), squeeze=False)
        for index, (ax, panel) in enumerate(zip(axes[0], panels)):
            info = panel["info"]
            ax.pcolormesh(panel["xedges"], panel["yedges"], panel["values"], cmap=cmap, norm=norm,
                          edgecolors="#d9d9d9", linewidth=.4, shading="flat", antialiased=False)
            ax.axhline(0, color="#222222", linewidth=1.3)
            ax.set(xlim=(-.5, .5), ylim=(-.5, .5), aspect="equal",
                   xlabel=r"Reduced $q_1$", ylabel=r"Reduced $q_2$")
            ax.set_xticks([-.5, -.25, 0, .25, .5]); ax.set_yticks([-.5, -.25, 0, .25, .5])
            ax.set_title(f"({chr(97+index)}) {info['label']}: " + rf"$Z_2={info['z2']}$", pad=10)
            top, bottom = info["half_top_sum"], info["half_bottom_sum"]
            top_text, bottom_text = (f"{v:+d}" if v else "0" for v in (top, bottom))
            ax.text(.5, -.16, f"Upper: {top_text} → {abs(top)%2}     Lower: {bottom_text} → {abs(bottom)%2}",
                    transform=ax.transAxes, ha="center", va="top", fontsize=10.5)
        fig.legend(handles=legend, loc="lower center", bbox_to_anchor=(.5, .01),
                   ncol=len(legend), frameon=False, title=r"Integer field $n(\mathbf{k})$")
        fig.subplots_adjust(left=.08, right=.98, bottom=.27, top=.91, wspace=.28)
        a.output_dir.mkdir(parents=True)
        for fmt in a.formats:
            fig.savefig(a.output_dir/f"figure.{fmt}", dpi=300)
        plt.close(fig)
    metadata = dict(schema="vaspberry.z2-field-comparison-plot", version=1,
                    operation="plot_existing_native_fields", sources=[p["info"] for p in panels],
                    interpolation="none", domain=[[-.5, .5], [-.5, .5]], half_zone_boundary="q2=0",
                    coordinate_definition="k=q1*b1+q2*b2; dimensionless reduced coordinates",
                    interpretation="Gauge- and branch-dependent integer field; agreed half-zone parity is the invariant, not the local tile pattern.",
                    numerical_computation="No wavefunctions, links or invariants recomputed; native result diagnostics and parities checked from CSV rows.",
                    plotter_sha256=sha(Path(__file__)),
                    validator_sha256=sha(Path(__file__).with_name("run.py")),
                    figures={fmt:sha(a.output_dir/f"figure.{fmt}") for fmt in a.formats})
    (a.output_dir/"plotting.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False)+"\n")
    print(json.dumps({"output":str(a.output_dir), "sources":[p["info"] for p in panels]}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
