#!/usr/bin/env python3
"""Run current VASPBERRY on the public Bi WAVECAR and compare Chern/curvature."""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from bi_common import arguments, execute, write_result, wavecar_geometry, sha


COORDINATE_TOLERANCE = 5.1e-7  # Native coordinates have six decimal places.
CURVATURE_TOLERANCE = 5.1e-5   # Native curvature has four decimal places.


def validate_field(raw, kpoints, reciprocal, reference):
    """Check the actual input mesh, unique native coverage and reference map."""
    if (np.shape(raw) != (1296, 7) or not np.isfinite(raw).all()
            or np.shape(kpoints) != (144, 3) or not np.isfinite(kpoints).all()
            or np.shape(reciprocal) != (3, 3) or not np.isfinite(reciprocal).all()):
        raise ValueError("expected finite native 3 x 3 display and actual 144-point input mesh")
    input_indices = np.rint(kpoints[:, :2] * 12).astype(int)
    input_residual = kpoints[:, :2] - input_indices / 12
    if np.max(abs(input_residual)) > 1e-7 or np.max(abs(kpoints[:, 2])) > 1e-7:
        raise ValueError("actual WAVECAR is not the expected Gamma-centered 12 x 12 plane")
    # Native loops are centered half a step below each input vertex.
    expected = {tuple((index - 1) % 12) for index in input_indices}
    if len(expected) != 144:
        raise ValueError("actual WAVECAR has duplicate or missing mesh points")

    def cells(rows):
        if rows.ndim != 2 or rows.shape[1] != 7 or not np.isfinite(rows).all():
            raise ValueError("nonfinite or malformed plaquette table")
        indices = np.rint(rows[:, 4:6] * 12 - .5).astype(int)
        if (np.max(abs(rows[:, 4:6] - (indices + .5) / 12)) > COORDINATE_TOLERANCE
                or np.max(abs(rows[:, 6])) > COORDINATE_TOLERANCE):
            raise ValueError("plaquette coordinates do not match the actual input mesh")
        # Propagate the six-decimal fractional-coordinate rounding into k_cart.
        cart_tolerance = COORDINATE_TOLERANCE * (1 + np.abs(reciprocal).sum(axis=0))
        if np.any(abs(rows[:, :3] - rows[:, 4:] @ reciprocal) > cart_tolerance):
            raise ValueError("Cartesian and fractional plaquette coordinates disagree with WAVECAR")
        return indices

    native_indices = cells(raw)
    if len({tuple(index) for index in native_indices}) != 1296:
        raise ValueError("duplicate plaquette coordinates in native output")
    keys, counts = np.unique(native_indices % 12, axis=0, return_counts=True)
    if {tuple(key) for key in keys} != expected or not np.all(counts == 9):
        raise ValueError("native tiled output does not cover each input plaquette nine times")
    field = raw[(raw[:, 4] >= 0) & (raw[:, 4] < 1) & (raw[:, 5] >= 0) & (raw[:, 5] < 1)]
    if field.shape != (144, 7):
        raise ValueError("expected exactly 144 first-BZ plaquettes")
    field_indices = cells(field)
    if {tuple(index) for index in field_indices} != expected:
        raise ValueError("first-BZ plaquettes have missing or duplicate input coordinates")
    if np.shape(reference) != (144, 7):
        raise ValueError("reference must contain 144 first-BZ plaquettes")
    reference_indices = cells(reference)
    if {tuple(index) for index in reference_indices} != expected:
        raise ValueError("reference plaquettes have missing or duplicate input coordinates")
    order = np.lexsort((field_indices[:, 1], field_indices[:, 0]))
    reference_order = np.lexsort((reference_indices[:, 1], reference_indices[:, 0]))
    maximum_error = float(np.max(abs(field[order, 3] - reference[reference_order, 3])))
    if maximum_error > CURVATURE_TOLERANCE:
        raise ValueError("native curvature map disagrees with the reference at printed precision")
    return field, maximum_error


def main():
    args = arguments(__doc__)
    out, provenance, edges, gaps = execute(args, "fukui-chern", ["-o", "BERRYCURV"])
    raw = np.loadtxt(out / "BERRYCURV.dat")
    kpoints, reciprocal = wavecar_geometry(args.wavecar)
    reference = Path(__file__).resolve().parent / "reference/plaquettes.csv"
    field, map_error = validate_field(raw, kpoints, reciprocal, np.loadtxt(reference, delimiter=",", skiprows=1))
    provenance["reference_plaquettes_sha256"] = sha(reference)
    provenance["coordinate_tolerance_fractional"] = COORDINATE_TOLERANCE
    provenance["curvature_reference_tolerance_A2"] = CURVATURE_TOLERANCE
    value = re.search(r"# Chern Number\s*=\s*([-+0-9.eE]+)", (out / "fortran.log").read_text())
    if not value:
        raise ValueError("no Chern result in native output")
    chern = float(value.group(1))
    if abs(chern) > 1e-5 or np.max(abs(field[:, 3])) > CURVATURE_TOLERANCE:
        raise ValueError("public Bi occupied bundle disagrees with the C = 0 reference")
    summary = {"mesh_nx": 12, "mesh_ny": 12, "plaquettes": 144, "occupied_band_min": 1,
               "occupied_band_max": 10, "chern_from_native_log": chern,
               "max_abs_printed_curvature_A2": float(np.max(abs(field[:, 3]))),
               "max_reference_curvature_difference_A2": map_error, **gaps}
    np.savetxt(out / "plaquettes.csv", field, delimiter=",", header="kx_inv_A,ky_inv_A,kz_inv_A,printed_curvature_A2,q1,q2,q3", comments="", fmt="%.6f")
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.3))
    mesh = axes[0].scatter(field[:, 0], field[:, 1], c=field[:, 3], marker="s", s=42,
                           cmap="RdBu_r", vmin=-1e-4, vmax=1e-4)
    axes[0].set(xlabel=r"$k_x$ ($\AA^{-1}$)", ylabel=r"$k_y$ ($\AA^{-1}$)",
                title="Bi bands 1–10\nC = 0", aspect="equal")
    colorbar = fig.colorbar(mesh, ax=axes[0], pad=.04, format="%.1e")
    colorbar.set_label(r"Native printed curvature ($\AA^2$)")
    colorbar.ax.tick_params(labelsize=8)
    axes[1].plot(edges[:, 0], edges[:, 1], lw=1.1, label="Band 10 (occupied)")
    axes[1].plot(edges[:, 0], edges[:, 2], lw=1.1, label="Band 11 (empty)")
    axes[1].axhspan(edges[:, 1].max(), edges[:, 2].min(), color="#21918c", alpha=.13,
                   label=f"Sampled global gap {gaps['sampled_global_gap_eV']:.6f} eV")
    axes[1].set(xlabel="Full-mesh VASP k-point index", ylabel="Energy (eV; VASP reference)",
                title="Selected bundle is separated on this mesh")
    axes[1].legend(fontsize=8, loc="best")
    fig.suptitle("Actual Bi WAVECAR → VASPBERRY Fukui calculation · 12 × 12 SOC", fontsize=13)
    fig.text(.5, .015, "Native BERRYCURV.dat prints four decimals: every value rounds to zero here.\n"
             "This is a sampled post-processing reference; denser-mesh convergence is not established.", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .10, 1, .94)); fig.savefig(out / "figure.png", dpi=180); plt.close(fig)
    write_result(out, "fukui-chern", provenance, summary, ["BERRYCURV.dat", "plaquettes.csv"], [
        "Real archived VASP input; reproduces VASPBERRY postprocessing, not an end-to-end new VASP calculation.",
        "Bi occupied bundle has C=0; this example does not demonstrate a nonzero-Chern material.",
        "Native curvature output has four decimal places and all values round to zero; no extra precision is inferred.",
        "WAVECAR pseudo-wavefunction overlaps omit PAW augmentation; mesh and gap convergence require separate calculations."])


if __name__ == "__main__":
    main()
