#!/usr/bin/env python3
"""Show a Cartesian Berry-curvature map, a marked path, and matching bands.

Fukui path values are periodic bilinear samples of plaquette averages. Kubo
curves use separately calculated point values on the supplied path. Numerical
energy zeros and raw curvature files are never changed by plotting.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.patches import Polygon
import matplotlib.patheffects as pe
import numpy as np

from plot_berry_curvature import (draw_bz_outline, plaquettes_in_first_bz,
    read_native_curvature, reciprocal_from_poscar, uniform_plaquettes)
from wavecar_fukui import Wavecar


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def export_bands(path_wavecar, destination, occupied=18):
    """Export unchanged VASP energies and path coordinates as a long-form CSV."""
    w = Wavecar(Path(path_wavecar), spinor_components=2)
    if not 0 < occupied < w.header.nbands:
        raise ValueError("occupied boundary must lie inside the stored band window")
    q = w.kpoints
    distance = np.r_[0., np.cumsum(np.linalg.norm(np.diff(q @ w.header.reciprocal, axis=0), axis=1))]
    destination = Path(destination)
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["k_index", "q1", "q2", "q3", "path_inv_A", "band", "energy_eV"])
        for ik, coordinate in enumerate(q):
            for ib, energy in enumerate(w.energies[ik]):
                writer.writerow([ik + 1, *coordinate, distance[ik], ib + 1, energy])
    return dict(source_wavecar_sha256=digest(path_wavecar), nkpoints=len(q),
                nbands=w.header.nbands, occupied_bands=[1, occupied],
                energy_zero="unchanged WAVECAR zero", output_sha256=digest(destination))


def read_bands(path, reciprocal):
    with Path(path).open() as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("empty band table")
    by_key = {(int(r["k_index"]), int(r["band"])): r for r in rows}
    nk = max(key[0] for key in by_key)
    nb = max(key[1] for key in by_key)
    if len(by_key) != len(rows) or set(by_key) != {(i, b) for i in range(1, nk + 1) for b in range(1, nb + 1)}:
        raise ValueError("band table needs unique consecutive 1-based k and band indices")
    q = np.array([[float(by_key[i, 1][f"q{j}"]) for j in (1, 2, 3)] for i in range(1, nk + 1)])
    distance = np.array([float(by_key[i, 1]["path_inv_A"]) for i in range(1, nk + 1)])
    energies = np.array([[float(by_key[i, b]["energy_eV"]) for b in range(1, nb + 1)] for i in range(1, nk + 1)])
    if not all(np.isfinite(a).all() for a in (q, distance, energies)):
        raise ValueError("nonfinite band data")
    for (i, _), row in by_key.items():
        if not np.allclose([float(row[f"q{j}"]) for j in (1, 2, 3)], q[i - 1], atol=1e-10, rtol=0):
            raise ValueError("band coordinates disagree at a shared k point")
        if not np.isfinite(float(row["path_inv_A"])) or abs(float(row["path_inv_A"]) - distance[i - 1]) > 1e-10:
            raise ValueError("band path distances disagree at a shared k point")
    actual = np.r_[0., np.cumsum(np.linalg.norm(np.diff(q @ reciprocal, axis=0), axis=1))]
    if not np.allclose(distance, actual, atol=1e-8, rtol=1e-8):
        raise ValueError("band path distance disagrees with POSCAR")
    if np.any(np.diff(distance) <= 0) or np.max(abs(q[:, 2])) > 1e-8:
        raise ValueError("path must have distinct successive points in the q3=0 plane")
    return q, distance, energies


def native_metadata(path):
    with Path(path).open() as handle:
        lines = handle.readlines()
    metadata = {}
    for line in lines:
        if line.startswith("#") and "=" in line:
            key, value = line[1:].strip().split("=", 1)
            metadata[key.strip()] = value.strip()
    return lines, metadata


def read_kubo(path, band, threshold):
    if not np.isfinite(threshold) or threshold < 0:
        raise ValueError("isolation threshold must be finite and nonnegative")
    lines, metadata = native_metadata(path)
    if metadata.get("normalization") != "STANDARD_MINUS_TWO_IM":
        raise ValueError("Kubo input must use the current physical normalization")
    rows = [r for r in csv.DictReader(line for line in lines if not line.startswith("#")) if int(r["band"]) == band]
    rows.sort(key=lambda r: int(r["k_index"]))
    if not rows or [int(r["k_index"]) for r in rows] != list(range(1, len(rows) + 1)):
        raise ValueError("Kubo input needs exactly one selected-band row per k point")
    q = np.array([[float(r[f"k{axis}_frac"]) for axis in "xyz"] for r in rows])
    energy = np.array([float(r["energy_eV"]) for r in rows])
    gap = np.array([float(r["min_gap_eV"]) for r in rows])
    values = np.array([float(r["omega_z_A2"]) for r in rows])
    if not all(np.isfinite(a).all() for a in (q, energy, gap, values)):
        raise ValueError("nonfinite native Kubo data")
    values[gap <= threshold] = np.nan
    return q, values, energy


def periodic_bilinear(q, values, target):
    """Sample a complete uniform periodic scalar grid; no new physical data."""
    values = np.asarray(values, dtype=float)
    target = np.asarray(target, dtype=float)
    if target.ndim != 2 or target.shape[1] != 3 or not np.isfinite(target).all():
        raise ValueError("interpolation targets must be finite fractional N x 3 coordinates")
    if np.isinf(values).any():
        raise ValueError("infinite curvature cannot be interpolated")
    if np.isnan(values).any():
        # Restore geometry independently; unknown scalar values stay unknown.
        q, indices, mesh = uniform_plaquettes(q, np.arange(len(values), dtype=float))
        values = values[np.rint(indices).astype(int)]
    else:
        q, values, mesh = uniform_plaquettes(q, values)
    axes = [np.unique(q[:, j]) for j in (0, 1)]
    data = np.empty(mesh)
    for point, value in zip(q, values):
        ij = [int(np.argmin(abs(axes[j] - point[j]))) for j in (0, 1)]
        data[tuple(ij)] = value
    u = (np.asarray(target)[:, :2] - [a[0] for a in axes]) * np.array(mesh)
    lower = np.floor(u).astype(int)
    fraction = u - lower
    sampled = np.zeros(len(u))
    valid = np.ones(len(u), dtype=bool)
    for dx, dy in ((0, 0), (0, 1), (1, 0), (1, 1)):
        weight = (fraction[:, 0] if dx else 1 - fraction[:, 0]) * (fraction[:, 1] if dy else 1 - fraction[:, 1])
        neighbour = data[(lower[:, 0] + dx) % mesh[0], (lower[:, 1] + dy) % mesh[1]]
        # Zero-weight unknown neighbours do not invalidate an exact known node.
        valid &= np.isfinite(neighbour) | (weight < 1e-12)
        sampled += weight * np.nan_to_num(neighbour, nan=0.)
    sampled[~valid] = np.nan
    return sampled


def read_bundle(path, occupied, threshold):
    """Read a guarded native bundle trace; never sum ill-conditioned band CSVs."""
    if not np.isfinite(threshold) or threshold < 0:
        raise ValueError("isolation threshold must be finite and nonnegative")
    lines, meta = native_metadata(path)
    if (meta.get("schema") != "VASPBERRY_BARE_MOMENTUM_KUBO_BUNDLE_V1"
            or meta.get("normalization") != "STANDARD_MINUS_TWO_IM"
            or meta.get("operator") != "WAVECAR_BARE_MOMENTUM_NO_PAW_NONLOCAL_VELOCITY"
            or meta.get("berry_connection") != "A_i=i<u|d/dk_i u>"
            or meta.get("intermediate_bands") != "EXTERNAL_TO_SELECTED_BUNDLE_WITHIN_SOURCE_NBANDS"
            or meta.get("result_kind") != "ISOLATED_BUNDLE_TRACE"
            or meta.get("result_status") != "PASS"):
        raise ValueError("a successful current native Kubo bundle CSV is required")
    if tuple(int(meta[k]) for k in ("band_min", "band_max", "band_rank")) != (1, occupied, occupied):
        raise ValueError("Kubo bundle range disagrees with the occupied-group label")
    rows = list(csv.DictReader(line for line in lines if not line.startswith("#")))
    rows.sort(key=lambda r: int(r["k_index"]))
    if not rows or [int(r["k_index"]) for r in rows] != list(range(1, len(rows) + 1)):
        raise ValueError("bundle input needs one selected-spin row per k point")
    if {int(r["spin"]) for r in rows} != {1}:
        raise ValueError("bundle plotter requires the spin-1 channel used by the band table")
    q = np.array([[float(r[f"k{axis}_frac"]) for axis in "xyz"] for r in rows])
    values = np.array([float(r["omega_z_A2"]) for r in rows])
    gap = np.array([float(r["min_external_gap_eV"]) for r in rows])
    if not all(np.isfinite(a).all() for a in (q, values, gap)) or np.any(gap <= threshold):
        raise ValueError("the occupied bundle is not isolated from the excluded states")
    return q, values, gap


def plot_panels(*, method, input_path, poscar_path, bands_csv, output_path,
                band=18, occupied=18, path_input=None, node_indices=(1, 25, 49),
                node_labels=("K", "Gamma", "Kprime"), curve_output=None,
                threshold=1e-5, title="1H-MoS₂", energy_limits=(-3., 3.),
                map_style="cells", display_grid=401):
    if map_style not in ("cells", "smooth"):
        raise ValueError("map style must be cells or smooth")
    if not isinstance(display_grid, int) or not 25 <= display_grid <= 2000:
        raise ValueError("display grid must be an integer between 25 and 2000")
    reciprocal = reciprocal_from_poscar(poscar_path)
    path_q, distance, energies = read_bands(bands_csv, reciprocal)
    nodes = np.array(node_indices, dtype=int) - 1
    if (len(nodes) != len(node_labels) or len(nodes) < 2 or nodes[0] != 0
            or nodes[-1] != len(path_q) - 1 or np.any(np.diff(nodes) <= 0)):
        raise ValueError("path nodes must include the first/last points in increasing order")
    if not 0 < occupied < energies.shape[1] or (method == "kubo" and not 0 < band <= energies.shape[1]):
        raise ValueError("band/occupied boundary outside band table")
    # Each marked segment must be the actual straight path used for the bands.
    for left, right in zip(nodes[:-1], nodes[1:]):
        fractions = (distance[left:right+1] - distance[left]) / (distance[right] - distance[left])
        expected = path_q[left] + fractions[:, None] * (path_q[right] - path_q[left])
        if not np.allclose(path_q[left:right+1], expected, atol=1e-8, rtol=0):
            raise ValueError("band path is not straight between the supplied symmetry nodes")
    if method == "fukui":
        selected = re.search(r"Chern Number for the BANDS\s*:\s*(\d+)\s*-\s*(\d+)", Path(input_path).read_text())
        if selected is None or tuple(map(int, selected.groups())) != (1, occupied):
            raise ValueError("native Fukui band range disagrees with the occupied-group label")
        q, values = read_native_curvature(input_path, reciprocal)
        q, values, mesh = uniform_plaquettes(q, values)
        curve = periodic_bilinear(q, values, path_q)
        curve_kind = "periodic_bilinear_display_cut_of_plaquette_averages"
        map_kind = "native_plaquette_averages"
        quantity = rf"Fukui: occupied bands 1–{occupied}"
        curve_label = "Map cut"
    elif method in ("kubo", "kubo-bundle"):
        if path_input is None:
            raise ValueError("Kubo panels require an actual path KUBO.csv")
        headers = [native_metadata(p)[1] for p in (input_path, path_input)]
        for key in ("intermediate_bands", "operator", "berry_connection", "normalization"):
            if not headers[0].get(key) or headers[0].get(key) != headers[1].get(key):
                raise ValueError("Kubo map/path conventions disagree: " + key)
        if method == "kubo-bundle":
            for key in ("band_min", "band_max", "band_rank", "source_nbands"):
                if not headers[0].get(key) or headers[0].get(key) != headers[1].get(key):
                    raise ValueError("Kubo bundle map/path settings disagree: " + key)
            if int(headers[0]["source_nbands"]) != energies.shape[1]:
                raise ValueError("Kubo bundle intermediate window disagrees with the band table")
            q, values, _ = read_bundle(input_path, occupied, threshold)
            path_native_q, curve, path_gaps = read_bundle(path_input, occupied, threshold)
        else:
            q, values, _ = read_kubo(input_path, band, threshold)
            path_native_q, curve, path_energy = read_kubo(path_input, band, threshold)
        if path_native_q.shape != path_q.shape or not np.allclose(path_native_q, path_q, atol=1e-8, rtol=0):
            raise ValueError("Kubo path coordinates disagree with band table")
        if method == "kubo" and not np.allclose(path_energy, energies[:, band - 1], atol=1e-8, rtol=0):
            raise ValueError("Kubo path energies disagree with band table")
        if method == "kubo-bundle":
            expected_gaps = np.min(abs(energies[:, :occupied, None] - energies[:, None, occupied:]), axis=(1, 2))
            if not np.allclose(path_gaps, expected_gaps, atol=1e-8, rtol=0):
                raise ValueError("Kubo bundle path gaps disagree with band table")
        curve_kind = "native_Kubo_at_actual_path_points"
        map_kind = "pointwise_Kubo_bundle_trace" if method == "kubo-bundle" else "pointwise_Kubo_band_curvature"
        quantity = rf"Kubo: occupied bands 1–{occupied}" if method == "kubo-bundle" else rf"Kubo: band {band}"
        curve_label = "Occupied bundle" if method == "kubo-bundle" else rf"Band {band}"
    else:
        raise ValueError("method must be fukui, kubo or kubo-bundle")
    if not np.isfinite(values).any() or not np.isfinite(curve).any():
        raise ValueError("no reportable curvature values")
    # Geometry uses finite cell IDs; masks are applied only to the colors.
    # Kubo cells visualize point samples, not plaquette-integrated curvature.
    polygons, ids, polygon, mesh = plaquettes_in_first_bz(q, np.arange(len(q), dtype=float), reciprocal)
    colors = values[np.rint(ids).astype(int)]
    bound = float(np.nanmax(abs(values))) or 1e-4
    energy_zero = float(np.max(energies[:, :occupied]))
    label_map = {"Gamma": r"$\Gamma$", "Kprime": "K′"}
    labels = [label_map.get(label, label) for label in node_labels]
    xy = path_q @ reciprocal
    with plt.rc_context({"font.size": 10, "axes.titlesize": 11, "axes.linewidth": .8}):
        fig = plt.figure(figsize=(9., 4.8), layout="constrained")
        gs = fig.add_gridspec(2, 2, width_ratios=(1.25, 1.), hspace=.07, wspace=.10)
        ax_map = fig.add_subplot(gs[:, 0])
        ax_band = fig.add_subplot(gs[0, 1])
        ax_curve = fig.add_subplot(gs[1, 1], sharex=ax_band)
        cmap = plt.get_cmap("RdBu_r").copy()
        cmap.set_bad("#d0d0d0")
        if map_style == "smooth":
            x = np.linspace(polygon[:, 0].min(), polygon[:, 0].max(), display_grid)
            y = np.linspace(polygon[:, 1].min(), polygon[:, 1].max(), display_grid)
            xx, yy = np.meshgrid(x, y)
            cart = np.column_stack((xx.ravel(), yy.ravel(), np.zeros(xx.size)))
            dense = periodic_bilinear(q, values, cart @ np.linalg.inv(reciprocal)).reshape(xx.shape)
            artist = ax_map.pcolormesh(xx, yy, np.ma.masked_invalid(dense), shading="nearest",
                                      cmap=cmap, norm=Normalize(-bound, bound), rasterized=True)
            artist.set_clip_path(Polygon(polygon, closed=True, transform=ax_map.transData))
        else:
            artist = PolyCollection(polygons, array=np.ma.masked_invalid(colors), cmap=cmap,
                                    norm=Normalize(-bound, bound), edgecolors="none", linewidths=0,
                                    antialiaseds=False, rasterized=True)
            ax_map.add_collection(artist)
        draw_bz_outline(ax_map, polygon, reciprocal, symmetry_labels=False)
        # Remove the generic Γ annotation; the chosen path labels are authoritative.
        for text in list(ax_map.texts):
            text.remove()
        ax_map.plot(xy[:, 0], xy[:, 1], color="#272727", lw=1.25, ls="--", zorder=7,
                    path_effects=[pe.Stroke(linewidth=2.7, foreground="white"), pe.Normal()])
        ax_map.scatter(xy[nodes, 0], xy[nodes, 1], s=16, color="#222222", zorder=8)
        for index, label in zip(nodes, labels):
            ax_map.annotate(label, xy[index, :2], xytext=(7, 6), textcoords="offset points", fontsize=11,
                            zorder=9, bbox=dict(facecolor="white", alpha=.85, edgecolor="none", pad=1))
        ax_map.set_title("(a) " + quantity, pad=12)
        bar = fig.colorbar(artist, ax=ax_map, location="bottom", shrink=.76, pad=.035, aspect=27)
        bar.set_label(r"$\Omega_z$ ($\mathrm{\AA}^2$)")
        bar.set_ticks(np.linspace(-bound, bound, 5))
        for b in range(energies.shape[1]):
            color = "#44678a" if b < occupied else "#888888"
            ax_band.plot(distance, energies[:, b] - energy_zero, color=color, lw=.8, alpha=.85)
        if method == "kubo":
            ax_band.plot(distance, energies[:, band - 1] - energy_zero, color="#b35806", lw=1.8,
                         label=rf"$n={band}$")
            ax_band.legend(loc="upper right", frameon=False, fontsize=8)
        ax_band.axhline(0, color=".65", lw=.6, zorder=0)
        ax_band.set(ylabel=r"$E-E_{\mathrm{v}}$ (eV)", ylim=energy_limits, title="(b) Band structure")
        ax_band.tick_params(labelbottom=False)
        ax_curve.plot(distance, curve, color="#b35806" if method.startswith("kubo") else "#2166ac",
                      lw=1.5, marker="o", ms=2., label=curve_label)
        ax_curve.axhline(0, color=".65", lw=.6, zorder=0)
        ax_curve.set(ylabel=r"$\Omega_z$ ($\mathrm{\AA}^2$)", title="(c) Curvature along the marked path",
                     xlabel="Wave vector")
        ax_curve.set_xticks(distance[nodes], labels)
        for ax in (ax_band, ax_curve):
            ax.set_xlim(distance[0], distance[-1])
            ax.tick_params(direction="in", right=True)
            for at in distance[nodes[1:-1]]:
                ax.axvline(at, color=".7", lw=.6, zorder=0)
        fig.suptitle(title, fontsize=13)
        fig.savefig(output_path, dpi=300)
        plt.close(fig)
    if curve_output is not None:
        with Path(curve_output).open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["k_index", "q1", "q2", "q3", "path_inv_A", "omega_z_A2", "valid"])
            for i, (point, at, value) in enumerate(zip(path_q, distance, curve), start=1):
                writer.writerow([i, *point, at, value if np.isfinite(value) else "", int(np.isfinite(value))])
    return dict(method=method, mesh=list(mesh), map_quantity=map_kind, path_quantity=curve_kind,
                axes="Cartesian kx, ky in Angstrom^-1", curvature_units="Angstrom^2",
                map_interpolation="periodic_bilinear" if map_style == "smooth" else "none",
                display_grid=[display_grid, display_grid] if map_style == "smooth" else None,
                path_nodes_1based=list(node_indices), path_labels=list(node_labels),
                path_vertices_fractional=path_q[nodes].tolist(), path_node_distance_inv_A=distance[nodes].tolist(),
                energy_reference_eV=energy_zero, energy_reference="maximum occupied energy on the supplied path",
                occupied_bands=[1, occupied], kubo_band=band if method == "kubo" else None,
                masked_map_points=int(np.isnan(values).sum()), masked_path_points=int(np.isnan(curve).sum()),
                isolation_threshold_eV=threshold if method.startswith("kubo") else None,
                color_limit_A2=bound, source_sha256={str(Path(p).name): digest(p) for p in (input_path, poscar_path, bands_csv)},
                path_input_sha256=digest(path_input) if path_input is not None else None,
                plotter_sha256=digest(__file__))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method", choices=("fukui", "kubo", "kubo-bundle"), required=True)
    parser.add_argument("--input", type=Path, required=True, help="native Fukui BERRYCURV.dat or mesh KUBO.csv")
    parser.add_argument("--poscar", type=Path, required=True)
    bands = parser.add_mutually_exclusive_group(required=True)
    bands.add_argument("--bands-csv", type=Path, help="long-form unchanged VASP band table")
    bands.add_argument("--path-wavecar", type=Path, help="export bands from this SOC path WAVECAR")
    parser.add_argument("--path-input", type=Path, help="native path KUBO.csv for --method kubo")
    parser.add_argument("--band", type=int, default=18)
    parser.add_argument("--occupied", type=int, default=18)
    parser.add_argument("--gap-threshold", type=float, default=1e-5, help="minimum isolated-band separation in eV")
    parser.add_argument("--energy-range", type=float, nargs=2, default=(-3., 3.), metavar=("MIN", "MAX"),
                        help="displayed energies relative to the occupied maximum, in eV")
    parser.add_argument("--path-node-indices", type=int, nargs="+", required=True, help="1-based k indices of symmetry nodes")
    parser.add_argument("--path-labels", nargs="+", required=True, help="e.g. K Gamma Kprime")
    parser.add_argument("--title", default="Berry curvature")
    parser.add_argument("--map-style", choices=("cells", "smooth"), default="cells")
    parser.add_argument("--display-grid", type=int, default=401, help="smooth display pixels per axis; no new k-point calculation")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not np.isfinite(args.energy_range).all() or args.energy_range[0] >= args.energy_range[1]:
        parser.error("--energy-range requires finite increasing limits")
    record = args.output.with_suffix(args.output.suffix + ".json")
    curve = args.output.with_name(args.output.stem + "_path.csv")
    exported = args.output.with_name(args.output.stem + "_bands.csv")
    if any(p.exists() for p in [args.output, record, curve] + ([exported] if args.path_wavecar else [])):
        parser.error("output exists; choose a new output stem")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    try:
        if args.path_wavecar:
            export_bands(args.path_wavecar, exported, args.occupied)
        metadata = plot_panels(method=args.method, input_path=args.input, poscar_path=args.poscar,
            bands_csv=exported if args.path_wavecar else args.bands_csv, output_path=args.output,
            band=args.band, occupied=args.occupied, path_input=args.path_input,
            node_indices=tuple(args.path_node_indices), node_labels=tuple(args.path_labels), curve_output=curve,
            title=args.title, threshold=args.gap_threshold, energy_limits=args.energy_range,
            map_style=args.map_style, display_grid=args.display_grid)
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))
    metadata.update(workflow_mode="plot_existing_native_outputs", figure_sha256=digest(args.output),
                    curve_sha256=digest(curve))
    record.write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
