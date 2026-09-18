#!/usr/bin/env python3
"""Run the public matrix-Kubo CLI and compare with an independent QWZ oracle."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from model_tools import begin, check_curvature, finish, model_curvature, result_provenance, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path, help="new directory; never overwritten")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    feature = Path(__file__).resolve().parent
    config = json.loads((feature / "input.json").read_text())
    manifest = begin(output, feature)
    try:
        data = model_curvature(config, output, manifest)
        checks, _, oracle = check_curvature(data, config)
        n = config["mesh"]
        with (output / "summary.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["k_id", "band_id", "kx_fractional", "ky_fractional",
                             "energy_eV", "omega_z_A2", "oracle_omega_z_A2", "error_A2"])
            # Sixteen reproducible points are sufficient for a compact review table.
            for ix in (0, n//4, n//2, 3*n//4):
                for iy in (0, n//4, n//2, 3*n//4):
                    k = ix*n + iy
                    for b in range(2):
                        value = data["omega_A2"][k,b,2]
                        writer.writerow([int(data["k_ids"][k]), b+1,
                            *data["kpoints_fractional"][k,:2], data["energies_eV"][k,b],
                            value, oracle[k,b], value-oracle[k,b]])
        fig, axes = plt.subplots(1, 2, figsize=(9, 3.9), constrained_layout=True)
        for ax, field, label in zip(axes,
            (data["omega_A2"][:,0,2], abs(data["omega_A2"][:,0,2]-oracle[:,0])),
            (r"Lower-band $\Omega_z$ ($\AA^2$)", r"Absolute oracle error ($\AA^2$)")):
            im = ax.imshow(field.reshape(n,n).T, origin="lower", extent=(0,1,0,1),
                           interpolation="nearest", cmap="viridis", aspect="equal")
            ax.set(xlabel=r"$k_x/(2\pi)$", ylabel=r"$k_y/(2\pi)$", title=label)
            fig.colorbar(im, ax=ax, shrink=.85)
        fig.suptitle(f"QWZ analytic model: matrix Kubo, m = −1, {n} × {n} mesh")
        fig.savefig(output / "figure.png", dpi=170)
        plt.close(fig)
        summary = {"schema_version": 1, "feature_id": "kubo-curvature", "workflow_mode": "analytic_recalculation",
            "model": config, "convention": "A=+i<u|d u>; ordered (kx,ky); a=1 Angstrom",
            "curvature_formula": "Omega_lower=(cos(kx)+cos(ky)+m*cos(kx)*cos(ky))/(2*|d|^3); upper=-lower",
            "units": {"curvature": "Angstrom^2", "energy": "eV", "chern": "dimensionless"},
            "numerical_checks": checks, "provenance": result_provenance(output, manifest),
            "scope": "Analytic two-band implementation/sign check. Finite-grid quadrature is not rounded to an integer. No material convergence claim."}
        write_json(output / "result.json", summary)
        if not checks["passed"]:
            raise RuntimeError("analytic model checks failed; see result.json")
        status = finish(output, feature, manifest, True)
        print(json.dumps(checks, indent=2))
        return status
    except Exception as exc:
        finish(output, feature, manifest, False, exc)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
