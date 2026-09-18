#!/usr/bin/env python3
"""Calculate MoS2 K–Gamma–K' curvature from the supplied VASP WAVECAR."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from workflow import ROOT, begin, fail, finish, run_command, sha256
sys.path.insert(0, str(ROOT / "tools"))
from wavecar_fukui import Wavecar

INPUT_SHA256 = "33f8546512856d6c04ad0a80454b18ac9b60e2af4b98f2f85b376ec49b1a8d9f"
GAP_THRESHOLD_EV = 1e-5


def read_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(line for line in stream if not line.startswith("#")))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--wavecar", type=Path, default=ROOT / "examples/1H-MoS2/KPATH/2.band/WAVECAR")
    parser.add_argument("--binary", type=Path, default=ROOT / "build/vaspberry-gfortran")
    args = parser.parse_args()
    output, wavecar, binary = args.output_dir.resolve(), args.wavecar.resolve(), args.binary.resolve()
    if not binary.is_file():
        parser.error("build the executable first with make serial, or pass --binary")
    feature = Path(__file__).resolve().parent
    record = begin(output, feature, wavecar, INPUT_SHA256)
    record["binary_sha256"] = sha256(binary)
    try:
        native = output / "native"
        native.mkdir()
        run_command([binary, "-f", wavecar, "-s", 2, "-kubo", 2, "-ii", 17, "-if", 18,
                     "-kubo_csv", "KUBO.csv", "-o", "BERRYCURV"], native, output, record, wavecar)
        rows = read_rows(native / "KUBO.csv")
        w = Wavecar(wavecar, spin=1, spinor_components=2)
        q = w.kpoints
        distance = np.r_[0., np.cumsum(np.linalg.norm(np.diff(q @ w.header.reciprocal, axis=0), axis=1))]
        summary = []
        energy_error = 0.
        coordinate_error = 0.
        gap_error = 0.
        for row in rows:
            k, b = int(row["k_index"])-1, int(row["band"])-1
            energy_error = max(energy_error, abs(float(row["energy_eV"])-w.energies[k,b]))
            coordinate_error = max(coordinate_error, float(np.max(abs(
                np.array([float(row[f"k{a}_frac"]) for a in "xyz"])-q[k]))))
            gap = float(np.min(abs(np.delete(w.energies[k], b)-w.energies[k,b])))
            gap_error = max(gap_error, abs(gap-float(row["min_gap_eV"])))
            valid = gap > GAP_THRESHOLD_EV
            summary.append({"k_index": k+1, "band": b+1, "path_inv_A": distance[k],
                "energy_eV": w.energies[k,b], "min_gap_eV": gap,
                "valid_isolated_band": int(valid),
                "omega_z_A2": float(row["omega_z_A2"]) if valid else ""})
        with (output / "summary.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
            writer.writeheader(); writer.writerows(summary)
        by_band = {b: [r for r in summary if r["band"] == b] for b in (17,18)}
        endpoint_residual = max(abs(float(data[0]["omega_z_A2"])+float(data[-1]["omega_z_A2"]))
                                for data in by_band.values())
        valid_counts = {str(b): sum(r["valid_isolated_band"] for r in data) for b,data in by_band.items()}
        reference = feature / "reference/summary.csv"
        reference_error = None
        reference_mask_matches = True
        if reference.exists():
            ref = read_rows(reference)
            reference_mask_matches = len(ref) == len(summary) and all(
                int(a["k_index"]) == b["k_index"] and int(a["band"]) == b["band"]
                and int(a["valid_isolated_band"]) == b["valid_isolated_band"] for a,b in zip(ref,summary))
            reference_error = max((abs(float(a["omega_z_A2"])-float(b["omega_z_A2"]))
                for a,b in zip(ref,summary) if a["omega_z_A2"] and b["omega_z_A2"] != ""), default=0.)
        fig, axes = plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
        for b,data in by_band.items():
            color = {17:"#2563a6",18:"#bb502d"}[b]
            axes[0].plot(distance,[r["energy_eV"] for r in data],color=color,label=f"Band {b}")
            values = np.array([r["omega_z_A2"] if r["valid_isolated_band"] else np.nan for r in data],float)
            axes[1].plot(distance,values,".-",markersize=3,color=color,label=f"Band {b}")
        axes[0].set(ylabel="WAVECAR energy (eV)",title="Valence bands from actual VASP states")
        axes[1].set(ylabel=r"$\Omega_z$ ($\AA^2$)",title="Native Kubo; unresolved states omitted")
        gamma = (distance[23]+distance[24])/2
        for ax in axes:
            ax.set_xticks([distance[0],gamma,distance[-1]],["K",r"$\Gamma$","K′"])
            ax.set_xlabel("K–Γ–K′ path")
            ax.grid(alpha=.2); ax.legend()
            ax.axvline(gamma,color="grey",lw=.7)
        axes[1].text(.5,.02,"Blank near Γ: gap ≤ 10⁻⁵ eV",ha="center",transform=axes[1].transAxes,fontsize=9)
        fig.suptitle("1H-MoS₂ · public VASP WAVECAR · bare-momentum approximation")
        fig.savefig(output / "figure.png",dpi=170); plt.close(fig)
        checks = {"native_rows":len(rows), "source_kpoints":w.header.nkpoints,
            "source_bands":w.header.nbands, "valid_points_per_band":valid_counts,
            "minimum_isolation_gap_eV":min(r["min_gap_eV"] for r in summary),
            "isolation_threshold_eV":GAP_THRESHOLD_EV,
            "energy_agreement_max_eV":energy_error,"coordinate_agreement_max_fractional":coordinate_error,
            "gap_agreement_max_eV":gap_error,"K_plus_Kprime_max_abs_A2":endpoint_residual,
            "reference_comparison_max_abs_A2":reference_error,"reference_validity_mask_matches":reference_mask_matches}
        checks["passed"] = bool(len(rows)==96 and valid_counts=={"17":44,"18":44}
            and energy_error<1e-10 and coordinate_error<1e-10 and gap_error<1e-10
            and endpoint_residual<1e-4 and reference_mask_matches
            and (reference_error is None or reference_error<1e-5))
        finish(output,record,{"feature_id":"kubo-curvature", "material":"1H-MoS2",
            "input_repository_path":"examples/1H-MoS2/KPATH/2.band/WAVECAR",
            "numerical_checks":checks,"units":{"energy":"eV (unchanged WAVECAR zero)","curvature":"Angstrom^2"},
            "sampling":"48-point K-Gamma-Kprime line; duplicate Gamma records 24 and 25",
            "scope":"Actual native Fortran calculation, bands 17:18, intermediate bands 1:32. "
                "Near-degenerate individual bands are masked at 1e-5 eV in the figure and summary; "
                "raw CSV retains diagnostic values. Canonical momentum omits PAW augmentation and nonlocal/SOC velocity terms. "
                "A line cannot define a Brillouin-zone Chern number or Hall conductivity."})
        print(f"PASS: {len(rows)} rows; bands 17/18 each have 44 isolated points; {output / 'figure.png'}")
    except Exception as exc:
        fail(output,record,exc)
        raise


if __name__ == "__main__":
    main()
