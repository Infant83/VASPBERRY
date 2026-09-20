#!/usr/bin/env python3
"""MoS2 Fukui curvature: calculate from a full WAVECAR mesh or replot its archive."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from plot_berry_curvature import (plot_native_curvature, read_native_curvature,
                                 reciprocal_from_poscar, uniform_plaquettes)
from wavecar_fukui import Wavecar, infer_uniform_grid


def sha(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--wavecar", type=Path, help="Actual full 12x12 SOC MoS2 WAVECAR; a line input is rejected")
    mode.add_argument("--archive-reference", action="store_true", help="Replot the supplied historical output without running VASPBERRY")
    parser.add_argument("--poscar", type=Path, default=ROOT / "examples/1H-MoS2/POSCAR")
    parser.add_argument("--binary", type=Path, default=ROOT / "build/vaspberry-gfortran")
    parser.add_argument("--output-dir", type=Path, required=True, help="New directory; existing results are preserved")
    args = parser.parse_args()
    poscar = args.poscar.resolve()
    reciprocal = reciprocal_from_poscar(poscar)
    if args.output_dir.exists():
        parser.error("output directory exists; choose a new path")
    wavecar = None
    if args.wavecar is not None:
        wavecar = Wavecar(args.wavecar.resolve(), spinor_components=2)
        infer_uniform_grid(wavecar.kpoints, 12, 12)
        if wavecar.header.nbands <= 18:
            parser.error("the tutorial needs occupied bands 1:18 plus at least one empty state")
        if not np.allclose(wavecar.header.reciprocal, reciprocal, rtol=0, atol=1e-8):
            parser.error("WAVECAR and POSCAR lattices disagree")
        if not (np.allclose(wavecar.occupations[:, :18], 1, atol=1e-6) and
                np.allclose(wavecar.occupations[:, 18:], 0, atol=1e-6)):
            parser.error("the tutorial requires 18 occupied SOC bands throughout the mesh")
        if not args.binary.is_file():
            parser.error("build VASPBERRY with make serial, or pass --binary")
    output = args.output_dir.resolve()
    output.mkdir(parents=True)
    report = {"schema_version": 1, "feature_id": "fukui-berry-curvature", "status": "RUNNING",
              "material": "1H-MoS2", "mesh": [12, 12], "occupied_bands": [1, 18],
              "execution_mode": "calculate_and_postprocess" if wavecar else "archived_native_output_replot",
              "workflow_mode": "vasp_wavecar_calculation" if wavecar else "archived_output_visualization",
              "commands": [], "provenance": {"poscar_sha256": sha(poscar), "analysis_script_sha256": sha(__file__),
                  "plotter_sha256": sha(ROOT / "tools/plot_berry_curvature.py"), "producer_revision": "unknown"},
              "units": {"k_cartesian": "Angstrom^-1", "berry_curvature": "Angstrom^2", "plaquette_flux": "radian"}}
    start = time.monotonic()
    try:
        native = output / "BERRYCURV.dat"
        if wavecar is not None:
            command = [str(args.binary.resolve()), "-f", str(wavecar.path), "-s", "2", "-kx", "12", "-ky", "12",
                       "-ii", "1", "-if", "18", "-o", "BERRYCURV"]
            with (output / "vaspberry.stdout.log").open("w") as stdout, (output / "vaspberry.stderr.log").open("w") as stderr:
                run = subprocess.run(command, cwd=output, stdout=stdout, stderr=stderr, timeout=180)
            report["commands"].append({"argv": command, "exit_code": run.returncode, "cwd": "."})
            report["provenance"].update(wavecar_sha256=sha(wavecar.path), source_sha256=sha(ROOT / "vaspberry.f"),
                                        binary_sha256=sha(args.binary), source_bands=wavecar.header.nbands,
                                        producer_revision="source-sha256:" + sha(ROOT / "vaspberry.f"))
            if run.returncode:
                raise RuntimeError("VASPBERRY failed; inspect vaspberry.stderr.log")
            sampled_gap = float(np.min(wavecar.energies[:, 18]-wavecar.energies[:, 17]))
            report["sampled_occupied_boundary_gap_eV"] = sampled_gap
            report["sampled_global_gap_eV"] = float(wavecar.energies[:, 18].min()-wavecar.energies[:, 17].max())
            if sampled_gap <= 1e-6:
                raise ValueError("occupied subspace is not separated from band 19 on the sampled mesh")
            if report["sampled_global_gap_eV"] <= 1e-6:
                raise ValueError("the sampled bands overlap in energy; this tutorial requires an insulator")
        else:
            archive = ROOT / "examples/1H-MoS2/BERRYCURV.dat"
            shutil.copyfile(archive, native)
            report["provenance"].update(input_repository_path="examples/1H-MoS2/BERRYCURV.dat",
                                        archived_output_sha256=sha(archive), source_bands=26,
                                        input_scope="Historical VASP-based native output; matching full-mesh WAVECAR is not supplied")
        raw_q, raw_values = read_native_curvature(native, reciprocal)
        q, values, mesh = uniform_plaquettes(raw_q, raw_values, mesh=(12, 12))
        area = float(np.linalg.norm(np.cross(reciprocal[0], reciprocal[1])) / np.prod(mesh))
        total_chern = float(values.sum()*area/(2*np.pi))
        tr_residual = 0.
        for coordinate, value in zip(q, values):
            distance = q[:, :2] + coordinate[:2]
            distance -= np.rint(distance)
            matches = np.flatnonzero(np.max(abs(distance), axis=1) < 2e-6)
            if len(matches) != 1:
                raise ValueError("missing or ambiguous time-reversed plaquette center")
            tr_residual = max(tr_residual, abs(value+values[matches[0]]))
        rows = []
        for coordinate, value in zip(q, values):
            centered = coordinate.copy()
            centered[:2] -= np.floor(centered[:2]+.5)
            rows.append([*centered, *(centered @ reciprocal), value, value*area])
        rows.sort(key=lambda row: (row[0], row[1]))
        with (output / "summary.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["q1", "q2", "q3", "kx_inv_A", "ky_inv_A", "kz_inv_A", "omega_z_A2", "flux_rad"])
            writer.writerows(rows)
        checks = {"chern_absolute_tolerance": 1e-4, "time_reversal_curvature_tolerance_A2": 1e-3,
                  "chern_near_zero": abs(total_chern) < 1e-4,
                  "time_reversal_odd": bool(tr_residual < 1e-3),
                  "nonzero_local_curvature": float(np.max(abs(values))) > 1,
                  "reference_comparison": "not_applicable"}
        reference = HERE / "reference"
        if wavecar is not None and (reference / "result.json").is_file():
            saved = json.loads((reference / "result.json").read_text())
            if saved["provenance"].get("wavecar_sha256") == report["provenance"]["wavecar_sha256"]:
                expected = np.loadtxt(reference / "summary.csv", delimiter=",", skiprows=1)
                actual = np.asarray(rows)
                np.testing.assert_allclose(actual[:, :6], expected[:, :6], rtol=0, atol=1e-7)
                error = float(np.max(abs(actual[:, 6]-expected[:, 6])))
                checks.update(reference_comparison="PASS" if error <= 1.1e-4 else "FAIL",
                              reference_max_curvature_difference_A2=error,
                              reference_absolute_tolerance_A2=1.1e-4)
        report["numerical_checks"] = checks
        if not (checks["chern_near_zero"] and checks["time_reversal_odd"] and checks["nonzero_local_curvature"]
                and checks["reference_comparison"] != "FAIL"):
            raise AssertionError("MoS2 reference sanity check failed; inspect result.json")
        title = "1H-MoS₂"
        plot = plot_native_curvature(native, poscar, output / "figure.png", title=title)
        plot_native_curvature(native, poscar, output / "figure.pdf", title=title)
        report.update(status="PASS", native_display_rows=len(raw_values), unique_plaquettes=len(values),
                      plaquette_area_inv_A2=area, total_chern_from_printed_values=total_chern,
                      minimum_curvature_A2=float(values.min()), maximum_curvature_A2=float(values.max()),
                      max_time_reversal_odd_residual_A2=float(tr_residual), plot=plot,
                      historical_unit_label_note="Native historical A^-2 label is incorrect; flux/deltaS has units A^2.",
                      interpretation="Opposite-sign K/Kprime valley curvature can coexist with vanishing total Chern number.",
                      output_sha256={name: sha(output / name) for name in ["BERRYCURV.dat", "summary.csv", "figure.png", "figure.pdf"]})
    except Exception as error:
        report.update(status="FAILED", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        report["elapsed_seconds"] = time.monotonic()-start
        (output / "result.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps({key: report[key] for key in ["status", "execution_mode", "unique_plaquettes", "total_chern_from_printed_values"]}, indent=2))


if __name__ == "__main__":
    main()
