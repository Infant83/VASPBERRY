#!/usr/bin/env python3
"""Plot selected channels from a standardized conductivity.csv; no new integration."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


LABELS = {
    "sigma_e2_over_h": r"Sheet Hall response ($e^2/h$)",
    "delta_sigma_e2_over_h": r"Change in sheet Hall response ($e^2/h$)",
    "sigma_S": "Sheet Hall conductance (S)",
    "electrons_per_cell": "Occupation sum (electrons / cell)",
    "delta_electrons_per_cell": "Occupation change (electrons / cell)",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="standard conductivity.csv")
    parser.add_argument("--output", type=Path, required=True, help="new PNG/PDF/SVG file")
    parser.add_argument("--quantity", choices=LABELS, default="sigma_e2_over_h")
    parser.add_argument("--regions", nargs="+", default=["total"])
    parser.add_argument("--bands", nargs="+", type=int, default=[0], help="0 means represented-band sum")
    parser.add_argument("--temperatures", nargs="+", type=float, help="default: all temperatures")
    parser.add_argument("--relative-mu", action="store_true", help="plot mu minus the declared reference")
    parser.add_argument("--title", default="Charge Hall post-processing")
    parser.add_argument("--dpi", type=int, default=180)
    args = parser.parse_args()
    try:
        if args.output.exists():
            raise ValueError("output exists; choose a fresh filename")
        if args.dpi <= 0:
            raise ValueError("dpi must be positive")
        xkey = "mu_minus_reference_eV" if args.relative_mu else "mu_eV"
        with args.input.open(newline="") as handle:
            reader = csv.DictReader(handle)
            required = {xkey, args.quantity, "temperature_K", "region", "band_id"}
            if not required <= set(reader.fieldnames or []):
                raise ValueError("CSV lacks required columns: " + ", ".join(sorted(required-set(reader.fieldnames or []))))
            rows = list(reader)
        channels = {}
        for row in rows:
            region, band, temperature = row["region"], int(row["band_id"]), float(row["temperature_K"])
            if region not in args.regions or band not in args.bands:
                continue
            if args.temperatures is not None and not any(math.isclose(temperature, t, abs_tol=1e-9, rel_tol=0)
                                                       for t in args.temperatures):
                continue
            x, y = float(row[xkey]), float(row[args.quantity])
            if not all(math.isfinite(v) for v in (x, y, temperature)):
                raise ValueError("selected CSV channel contains a nonfinite value")
            channels.setdefault((region, band, temperature), []).append((x, y))
        if not channels:
            raise ValueError("no rows match the requested channels")
        for region in args.regions:
            if not any(key[0] == region for key in channels):
                raise ValueError("requested region has no matching rows: " + region)
        for band in args.bands:
            if not any(key[1] == band for key in channels):
                raise ValueError("requested band has no matching rows: " + str(band))
        if args.temperatures is not None:
            for temperature in args.temperatures:
                if not any(math.isclose(key[2], temperature, abs_tol=1e-9, rel_tol=0) for key in channels):
                    raise ValueError("requested temperature has no matching rows: " + str(temperature))
        metadata_path = args.input.with_suffix(".json")
        metadata = json.loads(metadata_path.read_text()) if metadata_path.exists() else {}
        notes = []
        if metadata.get("scope") == "partial_band_contribution":
            notes.append("partial band contribution")
        operator = metadata.get("source_curvature_metadata", {}).get("source_operator", {})
        if operator.get("accuracy_status") in ("experimental", "approximation"):
            notes.append("operator: " + operator["accuracy_status"])
        fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
        for (region, band, temperature), points in sorted(channels.items()):
            points.sort()
            if len({x for x, _ in points}) != len(points):
                raise ValueError("duplicate chemical potentials in one selected channel")
            label = f"{region}; {'band sum' if band == 0 else 'band '+str(band)}; {temperature:g} K"
            ax.plot(*zip(*points), label=label.replace("$", r"\$"), linewidth=1.7)
        ax.axhline(0, color="0.65", linewidth=0.7, zorder=0)
        ax.set_xlabel(r"$\mu-\mu_{\rm ref}$ (eV)" if args.relative_mu else r"$\mu$ (eV, input reference)")
        ax.set_ylabel(LABELS[args.quantity])
        ax.set_title(args.title + ("\n" + "; ".join(notes) if notes else ""))
        ax.grid(alpha=0.18)
        ax.legend(fontsize=8)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.output, dpi=args.dpi)
        plt.close(fig)
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.error(str(error))
    print(args.output)


if __name__ == "__main__":
    main()
