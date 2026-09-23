#!/usr/bin/env python3
"""Plot measured DFT bands, an external full-connection Hall reference and convergence."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter, LogLocator

HERE = Path(__file__).resolve().parent


def records(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def path_positions(q, reciprocal, vertices, ticks):
    """Find every occurrence of an actual periodic k point on the plotted path."""
    found = []
    for start, end, offset in zip(vertices[:-1], vertices[1:], ticks[:-1]):
        direction = (end-start) @ reciprocal
        length = np.linalg.norm(direction)
        for i in (-1, 0, 1):
            for j in (-1, 0, 1):
                displacement = (q + [i, j, 0] - start) @ reciprocal
                fraction = np.dot(displacement, direction) / length**2
                if (-1e-10 <= fraction <= 1+1e-10
                        and np.linalg.norm(displacement-fraction*direction) < 1e-8):
                    coordinate = offset + fraction*length
                    if not any(abs(coordinate-x) < 1e-8 for x in found):
                        found.append(coordinate)
    return found


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference-dir", type=Path, default=HERE / "reference")
    p.add_argument("--output-dir", type=Path, required=True)
    a = p.parse_args()
    ref = a.reference_dir.resolve()
    meta = json.loads((ref / "wannier/metadata.json").read_text())
    bands = records(ref / "wannier/bands.csv")
    direct = records(ref / "direct-dft/sample-bands.csv")
    hall = records(ref / "wannier/conductivity.csv")
    convergence = [r for r in records(ref / "wannier/convergence.csv") if r["plot"] == "true"]
    mid = float(meta["reference_midgap_eV"])
    lattice = np.asarray(meta["lattice_A"])
    reciprocal = 2*np.pi*np.linalg.inv(lattice).T
    vertices = np.array([[0, 0, 0], [.5, 0, 0], [1/3, 1/3, 0], [0, 0, 0]])
    ticks = np.r_[0., np.cumsum(np.linalg.norm(np.diff(vertices, axis=0) @ reciprocal, axis=1))]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11.5,
                         "axes.labelsize": 12, "axes.titlesize": 12,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 3, figsize=(10.0, 4.0), gridspec_kw={"width_ratios": [1.4, 1, 1.1]})
    ax, ah, ac = axes
    inset = ax.inset_axes([.46, .60, .52, .35])
    for band in sorted({int(r["wannier_band"]) for r in bands}):
        selected = [r for r in bands if int(r["wannier_band"]) == band]
        x = np.array([float(r["distance_Ainv"]) for r in selected])
        e = np.array([float(r["energy_eV"])-mid for r in selected])
        if e.max() >= -1.35 and e.min() <= 1.35:
            ax.plot(x, e, color="#2864A0", lw=.85, zorder=1)
        if band in (87, 88):
            mask = x <= .024
            inset.plot(x[mask], 1000*e[mask], color="#2864A0", lw=1.1)
            mask = x >= ticks[-1]-.024
            inset.plot(x[mask]-ticks[-1], 1000*e[mask], color="#2864A0", lw=1.1)
    dots_x, dots_e = [], []
    for row in direct:
        energy = float(row["energy_eV"])-mid
        if not -1.35 <= energy <= 1.35:
            continue
        q = np.array([float(row[f"q{i}"]) for i in (1, 2, 3)])
        positions = path_positions(q, reciprocal, vertices, ticks)
        dots_x.extend(positions); dots_e.extend([energy]*len(positions))
        if int(row["band"]) in (123, 124):
            for position in positions:
                distance = position if position <= .024 else position-ticks[-1]
                if abs(distance) <= .024:
                    inset.plot(distance, 1000*energy, "o", ms=3, mfc="white", mec="#303030", mew=.75)
    ax.scatter(dots_x, dots_e, s=8, facecolors="none", edgecolors="#444444", linewidths=.55, zorder=3)
    ax.set(xlim=(0, ticks[-1]), ylim=(-1.35, 1.35), xticks=ticks,
           xticklabels=[r"$\Gamma$", "M", "K", r"$\Gamma$"], ylabel="Energy − midgap (eV)")
    for value in ticks[1:-1]: ax.axvline(value, color=".85", lw=.6, zorder=0)
    ax.axhline(0, color=".75", lw=.6, zorder=0)
    ax.legend(handles=[Line2D([], [], color="#2864A0", lw=1, label="Wannier90"),
                       Line2D([], [], marker="o", color="none", markeredgecolor=".3",
                              markerfacecolor="white", markersize=3.5, label="VASP")],
              frameon=True, facecolor="white", framealpha=.95, edgecolor="none",
              loc="lower left", fontsize=10.5, handlelength=1.4)
    inset.set(xlim=(-.023, .023), ylim=(-45, 50), xticks=[-.02, 0, .02], yticks=[-40, 0, 40])
    inset.tick_params(labelsize=10.5, pad=1)
    inset.set_xlabel(r"$\Delta k$ (Å$^{-1}$)", fontsize=10.5, labelpad=1)
    inset.set_ylabel("meV", fontsize=10.5, labelpad=1)
    for label in [*inset.get_xticklabels(), *inset.get_yticklabels(),
                  inset.xaxis.label, inset.yaxis.label]:
        label.set_bbox(dict(facecolor="white", edgecolor="none", pad=.4))
    inset.axhline(0, color=".8", lw=.5)
    inset.text(.98, .94, f"{1000*float(meta['dft_sampled_gap_eV']):.1f} meV", transform=inset.transAxes,
               fontsize=10.5, va="top", ha="right",
               bbox=dict(facecolor="white", edgecolor="none", pad=.5))
    mu = np.array([float(r["mu_minus_dft_midgap_eV"])*1000 for r in hall])
    sigma = np.array([float(r["sigma_xy_e2_over_h"]) for r in hall])
    vbm = 1000*(float(meta["dft_sampled_vbm_eV"])-mid)
    cbm = 1000*(float(meta["dft_sampled_cbm_eV"])-mid)
    ah.axvspan(vbm, cbm, color="#E6F1E8", zorder=0, label="DFT gap")
    ah.axhline(1, color=".55", lw=.7, ls="--", zorder=1)
    ah.plot(mu, sigma, "o-", color="#176B53", lw=1.4, ms=3.5, label="postw90 Kubo")
    ah.set(xlabel="Chemical potential\n− midgap (meV)", ylabel=r"$\sigma_{xy}$ ($e^2/h$)")
    ah.text(.97, .07, r"$C_{\mathrm{Fukui}}=-1$", ha="right", transform=ah.transAxes)
    ah.legend(frameon=False, fontsize=10.5, loc="lower left",
              bbox_to_anchor=(0, .17), handlelength=1.1, borderaxespad=.3)
    ah.set_ylim(min(0, float(sigma.min())-.08), max(1.12, float(sigma.max())+.08))
    x = np.array([int(r["point_count"]) for r in convergence])
    values = np.array([float(r["sigma_xy_e2_over_h"]) for r in convergence])
    errors = np.abs(values-1)
    ac.plot(x, errors, "o-", color="#8054A2", lw=1.2, ms=4)
    ac.set(xscale="log", yscale="log", xlabel="Integration points",
           ylabel=r"$|\sigma_{xy}/(e^2/h)-1|$")
    ac.set_xlim(min(2500, .8*x.min()), max(45000, 1.5*x.max()))
    ac.set_ylim(min(2e-5, .5*errors.min()), max(2e-3, 2*errors.max()))
    ac.set_xticks([3000, 10000, 30000], [r"$3\times10^3$", r"$10^4$", r"$3\times10^4$"])
    ac.xaxis.set_minor_formatter(NullFormatter())
    ac.yaxis.set_major_locator(LogLocator(base=10, subs=(1,)))
    ac.yaxis.set_minor_formatter(NullFormatter())
    offsets = [(5, 5), (5, -29), (5, 5), (-5, 6), (-5, -29)]
    for index, (px, py, row) in enumerate(zip(x, errors, convergence)):
        offset = offsets[min(index, len(offsets)-1)]
        ac.annotate(f"{row['base_mesh']} / {row['refinement']}\nR={row['radius_Ainv']}",
                    (px, py), xytext=offset, textcoords="offset points", fontsize=10.5,
                    ha="right" if offset[0] < 0 else "left",
                    bbox=dict(facecolor="white", edgecolor="none", alpha=.94, pad=.5))
    ac.grid(True, which="major", color=".9", lw=.6)
    for panel, label in zip(axes, ["(a) Bands", "(b) Sheet Hall response", "(c) Sampling convergence"]):
        panel.set_title(label, loc="left", pad=10)
    fig.tight_layout(w_pad=1.6)
    a.output_dir.mkdir(parents=True, exist_ok=True)
    for suffix in ["png", "pdf", "svg"]:
        fig.savefig(a.output_dir / f"mnbi2te4-qah.{suffix}", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(a.output_dir)


if __name__ == "__main__":
    main()
