#!/usr/bin/env python3
"""Calculate occupied MoS2 K–Gamma–K' trace curvature from the supplied WAVECAR."""
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
from run_fullmesh import validate_bundle
sys.path.insert(0, str(ROOT / "tools"))
from wavecar_fukui import Wavecar
from plot_berry_panels import read_bundle

INPUT_SHA256 = "33f8546512856d6c04ad0a80454b18ac9b60e2af4b98f2f85b376ec49b1a8d9f"
GAP_THRESHOLD_EV = 1e-5
OCCUPIED = 18


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
    record["source_sha256"]["examples/features/kubo-curvature/run_fullmesh.py"] = sha256(feature / "run_fullmesh.py")
    record["source_sha256"]["tools/plot_berry_panels.py"] = sha256(ROOT / "tools/plot_berry_panels.py")
    try:
        w = Wavecar(wavecar, spin=1, spinor_components=2)
        if (w.header.nkpoints != 48 or w.header.nbands != 32
                or np.max(abs(w.occupations[:, :OCCUPIED] - 1)) > 1e-5
                or np.max(abs(w.occupations[:, OCCUPIED:])) > 1e-5):
            raise ValueError("This supplied-path example requires 48 points, 32 spinor bands and 18 occupied bands")
        native = output / "native"
        native.mkdir()
        # An explicit range requests its trace. KUBO.csv is the native default.
        run_command([binary, "--task", "kubo", "--kubo-source", "wavecar", "--wavecar", wavecar,
                     "--bands", "1:18"], native, output, record, wavecar)
        checks = validate_bundle(native / "KUBO.csv", w)
        q, omega, gaps = read_bundle(native / "KUBO.csv", occupied=OCCUPIED, threshold=GAP_THRESHOLD_EV)
        distance = np.r_[0., np.cumsum(np.linalg.norm(np.diff(q @ w.header.reciprocal, axis=0), axis=1))]
        summary = [{"k_index": k+1, "path_inv_A": distance[k], "band_min": 1,
                    "band_max": OCCUPIED, "band_rank": OCCUPIED,
                    "E17_eV": w.energies[k, 16], "E18_eV": w.energies[k, 17],
                    "E19_eV": w.energies[k, 18], "min_external_gap_eV": gaps[k],
                    "omega_z_A2": omega[k]} for k in range(len(q))]
        with (output / "summary.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
            writer.writeheader(); writer.writerows(summary)
        endpoint_residual = float(abs(omega[0] + omega[-1]))
        energy_reference = float(w.energies[:, :OCCUPIED].max())
        fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.8), constrained_layout=True)
        for band, color in [(17, "#2166ac"), (18, "#b35806"), (19, "0.55")]:
            axes[0].plot(distance, w.energies[:, band-1] - energy_reference,
                         color=color, lw=1.5, label=rf"$n={band}$")
        axes[1].plot(distance, omega, color="#2166ac", lw=1.7, label="Occupied bands 1–18")
        axes[0].set(ylabel=r"$E-E_{\mathrm{v}}$ (eV)", title="(a) Bands at the occupied-space boundary")
        axes[1].set(ylabel=r"$\mathrm{Tr}\,\Omega_z$ ($\AA^2$)", title="(b) Occupied-space trace curvature")
        axes[1].axhline(0, color="0.65", lw=.6, zorder=0)
        gamma = (distance[23] + distance[24]) / 2
        for ax in axes:
            ax.set_xticks([distance[0], gamma, distance[-1]], ["K", r"$\Gamma$", "K′"])
            ax.set(xlabel="Wave vector", xlim=(distance[0], distance[-1]))
            ax.tick_params(direction="in", right=True)
            ax.axvline(gamma, color="0.65", lw=.6, zorder=0)
            ax.legend(frameon=False, fontsize=9)
            scale = ax.secondary_xaxis("top")
            scale.set_xticks([distance[0], gamma, distance[-1]],
                             [f"{value:.2f}" for value in (distance[0], gamma, distance[-1])])
            scale.set_xlabel(r"Path distance ($\AA^{-1}$)", fontsize=9)
            scale.tick_params(direction="in", labelsize=8)
        fig.savefig(output / "figure.png", dpi=220)
        fig.savefig(output / "figure.pdf")
        plt.close(fig)
        checks.update(source_kpoints=w.header.nkpoints, source_bands=w.header.nbands,
                      selected_bands=[1, OCCUPIED], isolation_threshold_eV=GAP_THRESHOLD_EV,
                      K_plus_Kprime_max_abs_A2=endpoint_residual,
                      historical_reference_comparison="not_applied: stored reference is a different individual-band quantity")
        checks["passed"] = bool(checks["native_rows"] == 48 and checks["valid_bundle_points"] == 48
                                and endpoint_residual < 1e-3)
        finish(output, record, {"feature_id": "kubo-curvature", "material": "1H-MoS2",
            "input_repository_path": "examples/1H-MoS2/KPATH/2.band/WAVECAR",
            "numerical_checks": checks,
            "units": {"energy": "eV (unchanged WAVECAR zero)", "curvature": "Angstrom^2"},
            "figure_conventions": {"energy_reference_eV": energy_reference,
                "energy_reference_definition": "maximum occupied energy along the supplied path; CSV energies remain unchanged",
                "abscissa": "cumulative Cartesian k distance, reciprocal vectors include 2*pi, Angstrom^-1",
                "symmetry_labels": ["K", "Gamma", "Kprime"],
                "symmetry_path_positions_inv_A": [float(distance[0]), float(gamma), float(distance[-1])]},
            "sampling": "48-point K-Gamma-Kprime line; duplicate Gamma records 24 and 25",
            "scope": "Actual native Fortran occupied1:18 trace; external intermediate bands19:32. "
                "Internal degeneracies are retained in the selected subspace; its external gap must exceed 1e-5 eV. "
                "Historical individual-band17:18 reference files are not replaced or used as a trace benchmark. "
                "Canonical momentum omits PAW augmentation and nonlocal/SOC velocity terms. "
                "A line cannot define a Brillouin-zone Chern number or Hall conductivity."})
        print(f"PASS: 48 occupied1:18 trace rows; {output / 'figure.png'}")
    except Exception as exc:
        fail(output, record, exc)
        raise


if __name__ == "__main__":
    main()
