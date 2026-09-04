#!/usr/bin/env python3
"""Plot a fundamental Fukui-Hatsugai integer n-field and optional legacy delta."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch, Rectangle
import numpy as np


COLORS = ("#2166ac", "#ffffff", "#b2182b")
CMAP = ListedColormap(COLORS)
NORM = BoundaryNorm((-1.5, -0.5, 0.5, 1.5), CMAP.N)
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


def read_current(path: Path) -> tuple[dict[tuple[float, float], int], dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        rows = (line for line in handle if not line.startswith("#"))
        reader = csv.DictReader(rows)
        required = {"q1", "q2", "nfield_int"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} is not a VASPBERRY Z2 field CSV")
        field = {
            (
                round(centered(float(row["q1"])), KEY_DIGITS),
                round(centered(float(row["q2"])), KEY_DIGITS),
            ): int(row["nfield_int"])
            for row in reader
        }
    if not field:
        raise ValueError(f"{path} contains no plaquettes")
    return field, metadata(path)


def read_legacy(path: Path) -> dict[tuple[float, float], int]:
    field: dict[tuple[float, float], int] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            columns = line.split()
            if len(columns) < 7:
                continue
            q1, q2 = float(columns[4]), float(columns[5])
            if -0.5 < q1 < 0.5 and -0.5 < q2 < 0.5:
                key = (round(q1, KEY_DIGITS), round(q2, KEY_DIGITS))
                field[key] = int(round(float(columns[3])))
    if not field:
        raise ValueError(f"{path} has no fundamental-zone n-field rows")
    return field


def legacy_reported_half_sums(path: Path) -> tuple[int, int] | None:
    """Read the half-zone sums printed by the archived implementation."""
    values: dict[str, int] = {}
    pattern = re.compile(
        r"^#\s*Z2 Invariant \((top|bottom)\)\s*=\s*([+-]?\d+)\s*$",
        re.IGNORECASE,
    )
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            match = pattern.match(line)
            if match:
                values[match.group(1).lower()] = int(match.group(2))
    if set(values) == {"top", "bottom"}:
        return values["top"], values["bottom"]
    return None


def current_reported_half_sums(values: dict[str, str]) -> tuple[int, int] | None:
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


def draw_field(
    ax: plt.Axes,
    field: dict[tuple[float, float], int],
    title: str,
    reported_sums: tuple[int, int] | None = None,
) -> None:
    xs, ys, values = grid(field)
    ax.pcolormesh(
        edges(xs), edges(ys), values, cmap=CMAP, norm=NORM,
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


def changed_cells(
    old: dict[tuple[float, float], int],
    new: dict[tuple[float, float], int],
) -> Iterable[tuple[float, float]]:
    if set(old) != set(new):
        missing = len(set(old) ^ set(new))
        raise ValueError(f"legacy/current grids differ at {missing} coordinates")
    return (key for key in new if old[key] != new[key])


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot a VASPBERRY Fukui-Hatsugai integer n-field.",
    )
    parser.add_argument("field_csv", type=Path, help="Z2_FIELD.csv input")
    parser.add_argument(
        "--legacy", type=Path,
        help="optional pre-correction NFIELD.dat to compare",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("Z2_nfield.pdf"),
        help="PDF or PNG output path (default: Z2_nfield.pdf)",
    )
    args = parser.parse_args()

    current, meta = read_current(args.field_csv)
    if len(current) != 144:
        raise ValueError(
            f"expected 144 plaquettes for this example, found {len(current)}"
        )

    legacy = read_legacy(args.legacy) if args.legacy else None
    if legacy is None:
        fig, axes = plt.subplots(
            1, 1, figsize=(6.0, 5.6), constrained_layout=False,
        )
        axes_list = [axes]
    else:
        fig, axes = plt.subplots(
            1, 2, figsize=(11.2, 5.2), constrained_layout=False,
        )
        axes_list = list(axes)
        draw_field(
            axes_list[0], legacy, "Incomplete pre-v1.1.1 field",
            legacy_reported_half_sums(args.legacy),
        )

    draw_field(
        axes_list[-1], current, "Corrected 12 x 12 field",
        current_reported_half_sums(meta),
    )

    if legacy is not None:
        xs, _, _ = grid(current)
        cell = float(np.median(np.diff(xs)))
        changed = list(changed_cells(legacy, current))
        for q1, q2 in changed:
            axes_list[-1].add_patch(
                Rectangle(
                    (q1 - cell / 2.0, q2 - cell / 2.0), cell, cell,
                    fill=False, edgecolor="#111111", linewidth=2.0,
                )
            )
        axes_list[-1].text(
            0.02, 0.02, f"outlined changes: {len(changed)}/144",
            transform=axes_list[-1].transAxes, ha="left", va="bottom",
            fontsize=9,
            bbox={
                "facecolor": "white", "alpha": 0.85,
                "edgecolor": "none", "pad": 2.0,
            },
        )

    legend = [
        Patch(facecolor=COLORS[0], edgecolor="#999999", label="n = -1"),
        Patch(facecolor=COLORS[1], edgecolor="#999999", label="n = 0"),
        Patch(facecolor=COLORS[2], edgecolor="#999999", label="n = +1"),
    ]
    fig.subplots_adjust(left=0.075, right=0.98, bottom=0.21, top=0.76, wspace=0.3)
    fig.legend(
        handles=legend, loc="lower center", bbox_to_anchor=(0.5, 0.075),
        ncol=3, frameon=False,
    )
    status = meta.get("result_status", "UNKNOWN")
    fig.suptitle(
        "Bi Fukui-Hatsugai integer n-field\n"
        f"fundamental Gamma-centered mesh; corrected field checks: {status}",
        fontsize=14, y=0.97,
    )
    fig.text(
        0.5, 0.025,
        "Panel sums reproduce output metadata. The pointwise n-field is "
        "gauge/branch dependent; Z2 is the half-zone sum modulo 2.",
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
