#!/usr/bin/env python3
"""Plot a reportable VASPBERRY Fukui-Hatsugai integer n-field."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
import numpy as np


KEY_DIGITS = 6


def centered(value: float) -> float:
    """Map a reduced coordinate to [-0.5, 0.5)."""
    return (value + 0.5) % 1.0 - 0.5


def metadata(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            body = line[1:].strip()
            if "=" in body:
                key, value = body.split("=", 1)
                values[key.strip()] = value.strip()
    return values


def read_field(
    path: Path,
) -> tuple[dict[tuple[float, float], int], dict[str, str]]:
    field: dict[tuple[float, float], int] = {}
    with path.open(encoding="utf-8") as handle:
        rows = (line for line in handle if not line.startswith("#"))
        reader = csv.DictReader(rows)
        required = {"q1", "q2", "nfield_int"}
        if (
            reader.fieldnames is None
            or not required.issubset(reader.fieldnames)
        ):
            raise ValueError(f"{path} is not a VASPBERRY Z2 field CSV")
        for row_number, row in enumerate(reader, start=2):
            key = (
                round(centered(float(row["q1"])), KEY_DIGITS),
                round(centered(float(row["q2"])), KEY_DIGITS),
            )
            if key in field:
                raise ValueError(
                    f"{path}: duplicate plaquette coordinate at CSV row "
                    f"{row_number}"
                )
            field[key] = int(row["nfield_int"])
    if not field:
        raise ValueError(f"{path} contains no plaquettes")
    return field, metadata(path)


def reported_half_sums(values: dict[str, str]) -> tuple[int, int] | None:
    """Read the half-zone sums from Z2_FIELD.csv metadata."""
    try:
        return (
            int(values["half_top_nfield_sum"]),
            int(values["half_bottom_nfield_sum"]),
        )
    except (KeyError, ValueError):
        return None


def grid(
    field: dict[tuple[float, float], int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xs = np.array(sorted({key[0] for key in field}), dtype=float)
    ys = np.array(sorted({key[1] for key in field}), dtype=float)
    if len(field) != len(xs) * len(ys):
        raise ValueError("n-field does not form a complete rectangular mesh")
    values = np.empty((len(ys), len(xs)), dtype=int)
    for iy, y in enumerate(ys):
        for ix, x in enumerate(xs):
            values[iy, ix] = field[
                (round(x, KEY_DIGITS), round(y, KEY_DIGITS))
            ]
    return xs, ys, values


def edges(points: np.ndarray) -> np.ndarray:
    if len(points) < 2:
        raise ValueError("at least two mesh points are required")
    step = float(np.median(np.diff(points)))
    return np.r_[points - step / 2.0, points[-1] + step / 2.0]


def half_sums(field: dict[tuple[float, float], int]) -> tuple[int, int]:
    top = sum(value for (_, q2), value in field.items() if q2 > 0.0)
    bottom = sum(value for (_, q2), value in field.items() if q2 < 0.0)
    return top, bottom


def parity(value: int) -> int:
    return abs(value) % 2


def validate_result(
    field: dict[tuple[float, float], int],
    values: dict[str, str],
    path: Path,
) -> tuple[int, int, int]:
    """Validate the current result contract and return ``(nx, ny, z2)``."""
    if values.get("schema") != "VASPBERRY_Z2_FIELD":
        raise ValueError(f"{path}: missing VASPBERRY_Z2_FIELD schema")
    try:
        schema_version = int(values["schema_version"])
        nx = int(values["nkx"])
        ny = int(values["nky"])
        z2 = int(values["z2_invariant"])
        top_parity = int(values["half_top_z2_parity"])
        bottom_parity = int(values["half_bottom_z2_parity"])
    except (KeyError, ValueError) as error:
        raise ValueError(f"{path}: incomplete Z2 result metadata") from error

    if schema_version < 2:
        raise ValueError(f"{path}: schema version 2 or newer is required")
    if values.get("result_kind") != "FUKUI_HATSUGAI_NFIELD_Z2":
        raise ValueError(f"{path}: unsupported Z2 result kind")
    if values.get("result_status") != "PASS":
        raise ValueError(f"{path}: result_status must be PASS")
    if values.get("reportable_invariant") != "1":
        raise ValueError(f"{path}: the Z2 invariant is not reportable")
    if values.get("half_bz_parity_consistent") != "1":
        raise ValueError(f"{path}: half-zone parities are inconsistent")
    if values.get("numerical_self_consistency_checks_pass") != "1":
        raise ValueError(f"{path}: numerical field checks did not pass")
    if z2 not in (0, 1) or top_parity != z2 or bottom_parity != z2:
        raise ValueError(f"{path}: inconsistent Z2 and half-zone parity metadata")
    if nx < 4 or ny < 4 or nx % 2 != 0 or ny % 2 != 0:
        raise ValueError(
            f"{path}: nkx and nky must be even integers greater than or "
            "equal to 4"
        )
    if len(field) != nx * ny:
        raise ValueError(
            f"{path}: expected {nx} x {ny}={nx * ny} plaquettes, "
            f"found {len(field)}"
        )

    xs, ys, _ = grid(field)
    if len(xs) != nx or len(ys) != ny:
        raise ValueError(
            f"{path}: metadata and rectangular mesh dimensions differ"
        )
    reported = reported_half_sums(values)
    calculated = half_sums(field)
    if reported is None or reported != calculated:
        raise ValueError(f"{path}: half-zone sums do not reproduce the metadata")
    if (
        parity(calculated[0]) != top_parity
        or parity(calculated[1]) != bottom_parity
    ):
        raise ValueError(f"{path}: n-field sums do not reproduce Z2 parities")
    return nx, ny, z2


def integer_style(
    field: dict[tuple[float, float], int],
) -> tuple[ListedColormap, BoundaryNorm, list[Patch]]:
    """Build a discrete, zero-centered color scale without clipping |n| > 1."""
    max_abs = max(1, max(abs(value) for value in field.values()))
    integers = np.arange(-max_abs, max_abs + 1, dtype=int)
    colors = plt.get_cmap("RdBu_r")(
        np.linspace(0.08, 0.92, len(integers))
    )
    colors[max_abs, :3] = 1.0
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(
        np.arange(-max_abs - 0.5, max_abs + 1.5, 1.0), cmap.N,
    )
    legend = [
        Patch(
            facecolor=colors[index], edgecolor="#999999",
            label=f"n = {value:+d}" if value else "n = 0",
        )
        for index, value in enumerate(integers)
    ]
    return cmap, norm, legend


def draw_field(
    ax: plt.Axes,
    field: dict[tuple[float, float], int],
    title: str,
    cmap: ListedColormap,
    norm: BoundaryNorm,
    reported_sums: tuple[int, int] | None = None,
) -> None:
    xs, ys, values = grid(field)
    ax.pcolormesh(
        edges(xs), edges(ys), values, cmap=cmap, norm=norm,
        edgecolors="#d9d9d9", linewidth=0.45, shading="flat",
    )
    ax.axhline(0.0, color="#222222", linewidth=1.7)
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(-0.5, 0.5)
    ax.set_aspect("equal")
    ax.set_xlabel(r"reduced $q_1$")
    ax.set_ylabel(r"reduced $q_2$")
    ax.set_xticks((-0.5, -0.25, 0.0, 0.25, 0.5))
    ax.set_yticks((-0.5, -0.25, 0.0, 0.25, 0.5))
    top, bottom = reported_sums if reported_sums is not None else half_sums(field)
    ax.set_title(
        f"{title}\n"
        f"top {top:+d} -> {parity(top)}, bottom {bottom:+d} -> {parity(bottom)}",
        fontsize=11,
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot a VASPBERRY Fukui-Hatsugai integer n-field.",
    )
    parser.add_argument("field_csv", type=Path, help="Z2_FIELD.csv input")
    parser.add_argument(
        "--output", type=Path, default=Path("Z2_nfield.pdf"),
        help="PDF or PNG output path (default: Z2_nfield.pdf)",
    )
    args = parser.parse_args()

    field, meta = read_field(args.field_csv)
    nx, ny, z2 = validate_result(field, meta, args.field_csv)
    cmap, norm, legend = integer_style(field)
    fig, ax = plt.subplots(
        1, 1, figsize=(6.2, 5.8), constrained_layout=False,
    )
    draw_field(
        ax, field, "Fundamental n-field", cmap, norm,
        reported_half_sums(meta),
    )

    fig.subplots_adjust(left=0.14, right=0.96, bottom=0.21, top=0.76)
    fig.legend(
        handles=legend, loc="lower center", bbox_to_anchor=(0.5, 0.075),
        ncol=3, frameon=False,
    )
    fig.suptitle(
        "Fukui-Hatsugai integer n-field\n"
        f"{nx} x {ny} fundamental mesh; Z2 = {z2}; result: PASS",
        fontsize=14, y=0.97,
    )
    fig.text(
        0.5, 0.025,
        "Half-zone sums reproduce the output metadata.\n"
        "Pointwise n-field values are gauge/branch dependent; "
        "Z2 is the half-zone sum modulo 2.",
        ha="center", va="bottom", fontsize=9,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=240, bbox_inches="tight")
    if args.output.suffix.lower() == ".pdf":
        fig.savefig(args.output.with_suffix(".png"), dpi=240, bbox_inches="tight")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
