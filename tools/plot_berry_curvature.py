#!/usr/bin/env python3
"""Plot native Berry curvature in the Cartesian reciprocal Wigner-Seitz cell.

Native Fukui values are drawn as constant plaquette averages. No interpolation
or smoothing is applied. The current interface supports reciprocal planes
parallel to Cartesian xy, with reciprocal vectors including 2*pi.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
import numpy as np


def reciprocal_from_poscar(path):
    lines = Path(path).read_text().splitlines()
    scale = np.array([float(value) for value in lines[1].split()])
    lattice = np.array([[float(value) for value in line.split()[:3]] for line in lines[2:5]])
    if lattice.shape != (3, 3) or not np.isfinite(lattice).all() or abs(np.linalg.det(lattice)) < 1e-12:
        raise ValueError("POSCAR must contain a finite nonsingular lattice")
    if len(scale) == 1:
        value = scale[0]
        if not np.isfinite(value) or value == 0:
            raise ValueError("invalid POSCAR scale")
        lattice *= value if value > 0 else (-value / abs(np.linalg.det(lattice))) ** (1 / 3)
    elif len(scale) == 3 and np.isfinite(scale).all() and np.all(scale > 0):
        lattice *= scale[np.newaxis, :]
    else:
        raise ValueError("POSCAR requires one nonzero or three positive scale values")
    return 2 * np.pi * np.linalg.inv(lattice).T


def xy_basis(reciprocal):
    reciprocal = np.asarray(reciprocal, dtype=float)
    if reciprocal.shape != (3, 3) or not np.isfinite(reciprocal).all():
        raise ValueError("reciprocal vectors must be a finite 3 x 3 array")
    if np.max(abs(reciprocal[:2, 2])) > 1e-10 * max(1., np.max(abs(reciprocal))):
        raise ValueError("this kx/ky plot requires b1 and b2 in the Cartesian xy plane")
    basis = reciprocal[:2, :2]
    if abs(np.linalg.det(basis)) < 1e-12:
        raise ValueError("reciprocal plane is singular")
    return basis


def polygon_area(vertices):
    vertices = np.asarray(vertices)
    if len(vertices) < 3:
        return 0.
    return .5 * abs(np.dot(vertices[:, 0], np.roll(vertices[:, 1], -1))
                    - np.dot(vertices[:, 1], np.roll(vertices[:, 0], -1)))


def _reduce_basis(basis):
    u, v = np.array(basis, dtype=float)
    for _ in range(100):
        if np.dot(v, v) < np.dot(u, u):
            u, v = v, u
        multiple = int(np.rint(np.dot(u, v) / np.dot(u, u)))
        if multiple == 0:
            return np.array([u, v])
        v = v - multiple * u
    raise ValueError("reciprocal basis reduction did not converge")


def _clip_halfplane(vertices, normal, offset, tolerance=1e-12):
    if len(vertices) == 0:
        return np.empty((0, 2))
    clipped = []
    previous = vertices[-1]
    previous_distance = float(np.dot(previous, normal) - offset)
    for current in vertices:
        distance = float(np.dot(current, normal) - offset)
        if (distance <= tolerance) != (previous_distance <= tolerance):
            fraction = previous_distance / (previous_distance - distance)
            clipped.append(previous + fraction * (current - previous))
        if distance <= tolerance:
            clipped.append(current)
        previous, previous_distance = current, distance
    return np.array(clipped).reshape(-1, 2)


def first_bz_polygon(reciprocal):
    """Return the Cartesian first reciprocal Wigner-Seitz polygon, CCW."""
    basis = xy_basis(reciprocal)
    reduced = _reduce_basis(basis)
    neighbours = np.array([i * reduced[0] + j * reduced[1]
                           for i in range(-2, 3) for j in range(-2, 3) if i or j])
    neighbours = neighbours[np.argsort(np.sum(neighbours ** 2, axis=1))]
    radius = 2 * np.linalg.norm(reduced, axis=1).sum()
    polygon = np.array([[-radius, -radius], [radius, -radius], [radius, radius], [-radius, radius]])
    for vector in neighbours:
        polygon = _clip_halfplane(polygon, vector, .5 * np.dot(vector, vector))
    # Remove coincident and collinear vertices introduced by redundant planes.
    cleaned = []
    for point in polygon:
        if not cleaned or np.linalg.norm(point - cleaned[-1]) > 1e-10 * radius:
            cleaned.append(point)
    polygon = np.array(cleaned)
    if len(polygon) > 1 and np.linalg.norm(polygon[0] - polygon[-1]) < 1e-10 * radius:
        polygon = polygon[:-1]
    keep = []
    for index, point in enumerate(polygon):
        left = point - polygon[index - 1]
        right = polygon[(index + 1) % len(polygon)] - point
        if abs(left[0] * right[1] - left[1] * right[0]) > 1e-10 * radius ** 2:
            keep.append(point)
    polygon = np.array(keep)
    if not np.isclose(polygon_area(polygon), abs(np.linalg.det(basis)), rtol=1e-9, atol=1e-12):
        raise ValueError("first-BZ area does not match reciprocal primitive-cell area")
    return polygon


def uniform_plaquettes(q, values, mesh=None):
    """Remove periodic display copies and recover the uniform printed grid."""
    q, values = np.asarray(q, dtype=float), np.asarray(values, dtype=float)
    if q.ndim != 2 or q.shape[1] != 3 or values.shape != (len(q),):
        raise ValueError("expected fractional coordinates N x 3 and N values")
    if not np.isfinite(q).all() or not np.isfinite(values).all():
        raise ValueError("nonfinite plaquette data")
    if np.ptp(q[:, 2]) > 1e-6 or abs(q[0, 2]) > 1e-6:
        raise ValueError("this first-BZ plot expects the q3 = 0 reciprocal plane")
    periodic = np.mod(q[:, :2], 1)
    periodic[np.isclose(periodic, 1, atol=1e-7, rtol=0)] = 0
    groups = {}
    for index, coordinate in enumerate(np.round(periodic, 6)):
        key = tuple(coordinate)
        if key in groups:
            if not np.isclose(values[index], groups[key], rtol=1e-10, atol=1e-10):
                raise ValueError("periodic display copies have inconsistent values")
        else:
            groups[key] = values[index]
    points = np.array(list(groups))
    data = np.array(list(groups.values()))
    axes = [np.unique(points[:, index]) for index in range(2)]
    inferred = tuple(len(axis) for axis in axes)
    if mesh is not None and tuple(mesh) != inferred:
        raise ValueError("declared mesh does not match plaquette coordinates")
    if min(inferred) < 2 or len(points) != inferred[0] * inferred[1]:
        raise ValueError("a complete uniform two-dimensional plaquette mesh is required")
    for index, axis in enumerate(axes):
        count = inferred[index]
        steps = np.diff(np.r_[axis, axis[0] + 1])
        if np.max(abs(steps - 1 / count)) > 2e-6:
            raise ValueError("plaquette mesh is not uniform at native coordinate precision")
        # Restore only rounded coordinates; the curvature values remain untouched.
        phase = float(np.mean(axis - np.arange(count) / count))
        ids = np.argmin(abs(points[:, index, None] - axis[None, :]), axis=1)
        points[:, index] = phase + ids / count
    return np.column_stack((points, np.zeros(len(points)))), data, inferred


def plaquettes_in_first_bz(q, values, reciprocal, mesh=None):
    """Clip periodic native plaquettes to the physical BZ without interpolation."""
    basis = xy_basis(reciprocal)
    polygon = first_bz_polygon(reciprocal)
    q, values, mesh = uniform_plaquettes(q, values, mesh)
    inverse = np.linalg.inv(basis)
    fractional_bounds = polygon @ inverse
    low = np.floor(fractional_bounds.min(axis=0) - 1).astype(int)
    high = np.ceil(fractional_bounds.max(axis=0) + 1).astype(int)
    if int(np.prod(high - low + 1)) * len(q) > 2_000_000:
        raise ValueError("plot tiling is too large; choose a less skew reciprocal primitive basis")
    half = .5 / np.array(mesh)
    offsets = np.array([[-half[0], -half[1]], [half[0], -half[1]],
                        [half[0], half[1]], [-half[0], half[1]]]) @ basis
    if np.linalg.det(basis) < 0:
        offsets = offsets[::-1]
    edges = np.roll(polygon, -1, axis=0) - polygon
    normals = np.column_stack((edges[:, 1], -edges[:, 0]))
    bounds = np.sum(normals * polygon, axis=1)
    polygons, colors = [], []
    minimum, maximum = polygon.min(axis=0), polygon.max(axis=0)
    for i in range(low[0], high[0] + 1):
        for j in range(low[1], high[1] + 1):
            centers = (q[:, :2] + [i, j]) @ basis
            for center, value in zip(centers, values):
                cell = center + offsets
                if np.any(cell.max(axis=0) < minimum - 1e-12) or np.any(cell.min(axis=0) > maximum + 1e-12):
                    continue
                for normal, bound in zip(normals, bounds):
                    cell = _clip_halfplane(cell, normal, bound)
                    if len(cell) < 3:
                        break
                if polygon_area(cell) > polygon_area(polygon) * 1e-13:
                    polygons.append(cell)
                    colors.append(value)
    covered = sum(polygon_area(cell) for cell in polygons)
    if not np.isclose(covered, polygon_area(polygon), rtol=2e-8, atol=1e-12):
        raise ValueError("plotted plaquettes do not cover the first BZ exactly once")
    return polygons, np.array(colors), polygon, mesh


def draw_bz_outline(ax, polygon, reciprocal, symmetry_labels=True, k_fractional=(1 / 3, 2 / 3)):
    closed = np.vstack((polygon, polygon[0]))
    ax.plot(closed[:, 0], closed[:, 1], color="0.15", lw=1.15, zorder=5)
    ax.scatter([0], [0], s=8, color="0.15", zorder=6)
    radius = np.max(np.linalg.norm(polygon, axis=1))
    ax.annotate(r"$\Gamma$", (0, 0), xytext=(5, 5), textcoords="offset points", fontsize=11)
    lengths = np.linalg.norm(np.roll(polygon, -1, axis=0) - polygon, axis=1)
    if symmetry_labels and len(polygon) == 6 and np.allclose(lengths, lengths.mean(), rtol=1e-6):
        # Preserve the declared valley convention rather than naming an arbitrary corner.
        # Default K=(1/3,2/3), K'=-K matches the supplied MoS2 KPOINTS path.
        basis = xy_basis(reciprocal)
        chosen = np.asarray(k_fractional, dtype=float)
        if chosen.shape != (2,) or not np.isfinite(chosen).all():
            raise ValueError("K fractional coordinate must contain two finite numbers")
        for sign, label in [(1, r"$K$"), (-1, r"$K'$" )]:
            differences = (polygon @ np.linalg.inv(basis)) - sign * chosen
            matches = np.max(abs(differences - np.rint(differences)), axis=1) < 1e-7
            if matches.any():
                indices = np.flatnonzero(matches)
                direction = sign * chosen @ basis
                selected = indices[np.argmax(polygon[indices] @ direction)]
                ax.text(*(polygon[selected] * .88), label, ha="center", va="center", fontsize=11,
                        bbox={"facecolor": "white", "alpha": .65, "edgecolor": "none", "pad": .5})
        midpoints = .5 * (polygon + np.roll(polygon, -1, axis=0))
        midpoint = midpoints[np.argmax(midpoints[:, 1])]
        ax.text(*(midpoint * 1.10), r"$M$", ha="center", va="center", fontsize=11)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(polygon[:, 0].min() - .16 * radius, polygon[:, 0].max() + .16 * radius)
    ax.set_ylim(polygon[:, 1].min() - .16 * radius, polygon[:, 1].max() + .16 * radius)
    ax.set_xlabel(r"$k_x$ ($\mathrm{\AA}^{-1}$)")
    ax.set_ylabel(r"$k_y$ ($\mathrm{\AA}^{-1}$)")
    ax.tick_params(direction="in", top=True, right=True)


def draw_plaquette_map(ax, q, values, reciprocal, mesh=None, *, cmap="RdBu_r", norm=None, vmax=None,
                       symmetry_labels=True, k_fractional=(1 / 3, 2 / 3)):
    polygons, colors, polygon, mesh = plaquettes_in_first_bz(q, values, reciprocal, mesh)
    if norm is None:
        bound = float(np.max(abs(colors))) if vmax is None else float(vmax)
        if not np.isfinite(bound) or bound < 0:
            raise ValueError("color limit must be finite and nonnegative")
        if bound == 0:
            bound = 1e-4  # Native curvature's printed resolution, not invented signal.
        norm = Normalize(vmin=-bound, vmax=bound)
    collection = PolyCollection(polygons, array=colors, cmap=cmap, norm=norm,
                                edgecolors="none", linewidths=0, antialiaseds=False, rasterized=True)
    ax.add_collection(collection)
    draw_bz_outline(ax, polygon, reciprocal, symmetry_labels, k_fractional)
    return collection, {"mesh": list(mesh), "first_bz_area_inv_A2": polygon_area(polygon),
                        "rendering": "constant_native_plaquette_values_clipped_to_first_BZ",
                        "interpolation": "none", "axes": "Cartesian kx, ky in Angstrom^-1",
                        "K_fractional": list(k_fractional), "Kprime_fractional": list(-np.asarray(k_fractional)),
                        "color_min": float(norm.vmin), "color_max": float(norm.vmax)}


def read_native_curvature(path, reciprocal):
    with Path(path).open() as handle:
        for line in handle:
            header = line.lstrip().upper()
            if (header.startswith("# KUBO") or "SCHEMA=VASPBERRY_BARE_MOMENTUM_KUBO" in header
                    or (header.startswith("#") and "BERRY CURVATURE" in header and "KUBO" in header)):
                raise ValueError("this plotter expects Fukui plaquette averages, not Kubo point data")
    raw = np.loadtxt(path)
    if raw.ndim != 2 or raw.shape[1] != 7 or not np.isfinite(raw).all():
        raise ValueError("expected native seven-column BERRYCURV.dat")
    basis = xy_basis(reciprocal)
    tolerance = 5.1e-7 * (1 + np.abs(reciprocal).sum(axis=0))
    if np.any(abs(raw[:, :3] - raw[:, 4:] @ reciprocal) > tolerance):
        raise ValueError("BERRYCURV coordinates disagree with the supplied POSCAR lattice")
    # Fukui flux/positive cell area is the component along b1 x b2.
    # Express the displayed scalar as Cartesian Omega_z also for reversed bases.
    return raw[:, 4:], raw[:, 3] * np.sign(np.linalg.det(basis))


def plot_native_curvature(input_path, poscar_path, output_path, *, title="Berry curvature", vmax=None, dpi=300,
                          k_fractional=(1 / 3, 2 / 3)):
    reciprocal = reciprocal_from_poscar(poscar_path)
    q, values = read_native_curvature(input_path, reciprocal)
    with plt.rc_context({"font.size": 11, "axes.titlesize": 13, "axes.linewidth": .8}):
        fig, ax = plt.subplots(figsize=(5.2, 4.5), layout="constrained")
        artist, metadata = draw_plaquette_map(ax, q, values, reciprocal, vmax=vmax, k_fractional=k_fractional)
        metadata["native_to_cartesian_z_sign"] = int(np.sign(np.linalg.det(xy_basis(reciprocal))))
        ax.set_title(title, pad=11)
        colorbar = fig.colorbar(artist, ax=ax, shrink=.72, pad=.035,
                               ticks=np.linspace(artist.norm.vmin, artist.norm.vmax, 5))
        colorbar.set_label(r"$\Omega_z$ ($\mathrm{\AA}^{2}$)")
        colorbar.formatter.set_powerlimits((-3, 3))
        colorbar.update_ticks()
        fig.savefig(output_path, dpi=dpi)
        plt.close(fig)
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="native BERRYCURV.dat")
    parser.add_argument("--poscar", type=Path, required=True, help="matching POSCAR defining reciprocal lattice")
    parser.add_argument("--output", type=Path, required=True, help="new figure filename, PNG/PDF/SVG")
    parser.add_argument("--title", default="Berry curvature")
    parser.add_argument("--vmax", type=float, help="symmetric color limit in Angstrom squared")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--k-fractional", default="0.3333333333333333,0.6666666666666666",
                        help="two fractional coordinates of K; K-prime is -K (default matches MoS2 KPOINTS)")
    args = parser.parse_args()
    record_path = args.output.with_suffix(args.output.suffix + ".json")
    if args.output.exists() or record_path.exists():
        parser.error("output exists; choose a new figure filename")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    try:
        k_fractional = tuple(float(value) for value in args.k_fractional.split(","))
        if len(k_fractional) != 2 or not np.isfinite(k_fractional).all():
            raise ValueError("--k-fractional requires two finite comma-separated coordinates")
        metadata = plot_native_curvature(args.input, args.poscar, args.output, title=args.title,
                                         vmax=args.vmax, dpi=args.dpi, k_fractional=k_fractional)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    digest = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
    metadata.update(schema_version=1, workflow_mode="plot_existing_native_output", input=str(args.input),
                    poscar=str(args.poscar), input_sha256=digest(args.input), poscar_sha256=digest(args.poscar),
                    plotter_sha256=digest(__file__), figure_sha256=digest(args.output))
    record_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
