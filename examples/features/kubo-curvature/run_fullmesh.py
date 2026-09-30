#!/usr/bin/env python3
"""Reproduce the MoS2 Kubo BZ map, symmetry-path curvature and band panel.

Both WAVECARs must use the documented SOC MoS2 settings, 26 bands, and the
same lattice. The path uses 49 explicit points K--Gamma--Kprime, with Gamma
included once. This tutorial validates the sampling instead of pinning an
input checksum, so regenerated VASP files can be used.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
from wavecar_fukui import Wavecar

GAP_THRESHOLD_EV = 1e-5


def sha256(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def validate_inputs(mesh, path):
    for label, wavecar, nk in (("mesh", mesh, 144), ("path", path, 49)):
        h = wavecar.header
        if h.ispin != 1 or h.nbands != 26 or h.nkpoints != nk:
            raise ValueError(f"{label} WAVECAR requires ISPIN=1, NBANDS=26 and {nk} k points")
        if not np.isclose(h.encut_ev, 400., atol=1e-8, rtol=0):
            raise ValueError(f"{label} WAVECAR must use the documented 400 eV cutoff")
        if np.max(abs(wavecar.kpoints[:, 2])) > 1e-8:
            raise ValueError(f"{label} WAVECAR must sample the q3=0 plane")
        if (np.max(abs(wavecar.occupations[:, :18] - 1)) > 1e-5
                or np.max(abs(wavecar.occupations[:, 18:])) > 1e-5):
            raise ValueError(f"{label} WAVECAR must have exactly 18 occupied SOC bands")
        if np.min(wavecar.energies[:, 18] - wavecar.energies[:, 17]) <= GAP_THRESHOLD_EV:
            raise ValueError(f"{label} occupied manifold is not isolated from band 19")
    if not np.allclose(mesh.header.lattice, path.header.lattice, atol=1e-10, rtol=0):
        raise ValueError("mesh and path WAVECAR lattice vectors differ")
    mesh_ids = np.rint(mesh.kpoints[:, :2] * 12).astype(int)
    if np.max(abs(mesh.kpoints[:, :2] * 12 - mesh_ids)) > 1e-7:
        raise ValueError("mesh WAVECAR must use a Gamma-centered 12 x 12 grid")
    if len(np.unique(mesh_ids % 12, axis=0)) != 144:
        raise ValueError("mesh WAVECAR does not contain each 12 x 12 grid point once")
    endpoint = np.array([1 / 3, 2 / 3, 0.])
    expected = np.linspace(endpoint, -endpoint, 49)
    if not np.allclose(path.kpoints, expected, atol=1e-8, rtol=0):
        raise ValueError("path WAVECAR must follow 49 explicit K-Gamma-Kprime points, Gamma at index 25")
    return {
        "mesh_kpoints": 144, "path_kpoints": 49, "source_bands": 26,
        "occupied_bands": "1:18", "selected_bands": [1, 18],
        "intermediate_bands": "19:26 (internal selected-space pairs cancel)", "mesh": [12, 12],
        "path_node_indices_one_based": [1, 25, 49],
        "minimum_occupied_direct_gap_mesh_eV": float(np.min(mesh.energies[:, 18] - mesh.energies[:, 17])),
        "minimum_occupied_direct_gap_path_eV": float(np.min(path.energies[:, 18] - path.energies[:, 17])),
    }


def validate_bundle(csv_path, wavecar):
    from plot_berry_panels import read_bundle
    q, omega, gaps = read_bundle(csv_path, occupied=18, threshold=GAP_THRESHOLD_EV)
    if q.shape != wavecar.kpoints.shape or not np.allclose(q, wavecar.kpoints, atol=1e-10, rtol=0):
        raise ValueError("bundle CSV coordinates disagree with WAVECAR")
    external_gaps = np.min(abs(wavecar.energies[:, :18, None] - wavecar.energies[:, None, 18:]), axis=(1, 2))
    if not np.allclose(gaps, external_gaps, atol=1e-10, rtol=0):
        raise ValueError("bundle CSV external gaps disagree with WAVECAR")
    return dict(native_rows=len(q), valid_bundle_points=len(q),
                minimum_external_gap_eV=float(gaps.min()),
                gap_agreement_max_eV=float(np.max(abs(gaps - external_gaps))),
                coordinate_agreement_max_fractional=float(np.max(abs(q - wavecar.kpoints))),
                omega_range_A2=[float(omega.min()), float(omega.max())])


def run_native(binary, wavecar, name, output, record):
    directory = output / name
    directory.mkdir()
    command = [str(binary), "--task", "kubo", "--wavecar", str(wavecar),
               "--bands", "1:18"]
    stage = {"argv": command, "cwd": str(directory), "timeout_s": 180}
    record["commands"].append(stage)
    write_json(output / "provenance.json", record)
    start = time.monotonic()
    with (directory / "stdout.log").open("w") as stdout, (directory / "stderr.log").open("w") as stderr:
        try:
            proc = subprocess.run(command, cwd=directory, stdout=stdout, stderr=stderr, timeout=180)
            stage["exit_code"] = proc.returncode
        except subprocess.TimeoutExpired:
            stage["exit_code"] = None
            stage["status"] = "FAIL_TIMEOUT"
            raise
        finally:
            stage["elapsed_s"] = time.monotonic() - start
            write_json(output / "provenance.json", record)
    if proc.returncode:
        raise RuntimeError(f"{name} native calculation failed; see its stderr.log")
    stage["status"] = "PASS"
    return directory / "KUBO.csv"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wavecar", required=True, type=Path, help="Gamma-centered 12 x 12 SOC WAVECAR")
    parser.add_argument("--path-wavecar", required=True, type=Path, help="matched 49-point K-Gamma-Kprime SOC WAVECAR")
    parser.add_argument("--binary", type=Path, default=ROOT / "build/vaspberry-gfortran")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--map-style", choices=("smooth", "cells"), default="smooth")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    record = {"status": "RUNNING", "feature_id": "kubo-curvature-fullmesh",
              "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "commands": [], "inputs": {}}
    write_json(output / "provenance.json", record)
    try:
        from plot_berry_panels import export_bands, plot_panels
        wavecar, path_wavecar, binary = args.wavecar.resolve(), args.path_wavecar.resolve(), args.binary.resolve()
        if not binary.is_file():
            raise ValueError("build the executable first with make serial, or pass --binary")
        for name, source in (("wavecar", wavecar), ("path_wavecar", path_wavecar)):
            record["inputs"][name] = {"path": str(source), "bytes": source.stat().st_size, "sha256": sha256(source)}
        record["binary_sha256"] = sha256(binary)
        record["version"] = (ROOT / "VERSION").read_text().strip()
        record["source_sha256"] = {str(source.relative_to(ROOT)): sha256(source) for source in (
            Path(__file__).resolve(), ROOT / "vaspberry.f", ROOT / "tools/wavecar_fukui.py",
            ROOT / "tools/plot_berry_panels.py", ROOT / "tools/plot_berry_curvature.py")}
        mesh = Wavecar(wavecar, spin=1, spinor_components=2)
        path = Wavecar(path_wavecar, spin=1, spinor_components=2)
        settings = validate_inputs(mesh, path)
        settings["quantity"] = "occupied_band_trace"
        poscar = output / "POSCAR.lattice"
        poscar.write_text("MoS2 reciprocal geometry from WAVECAR\n1.0\n" +
            "\n".join(" ".join(f"{value:.17g}" for value in vector) for vector in mesh.header.lattice) +
            "\nX\n1\nDirect\n0 0 0\n")
        mesh_csv = run_native(binary, wavecar, "native-map", output, record)
        path_csv = run_native(binary, path_wavecar, "native-path", output, record)
        numerical = {"mesh": validate_bundle(mesh_csv, mesh), "path": validate_bundle(path_csv, path)}
        bands_csv = output / "bands.csv"
        band_metadata = export_bands(path_wavecar, bands_csv, occupied=18)
        common = dict(method="kubo-bundle", input_path=mesh_csv, poscar_path=poscar,
                      bands_csv=bands_csv, band=18, path_input=path_csv,
                      node_indices=(1, 25, 49), node_labels=("K", "Gamma", "Kprime"), map_style=args.map_style)
        plot_metadata = plot_panels(**common, output_path=output / "figure.png",
                                    curve_output=output / "path_curvature.csv")
        plot_panels(**common, output_path=output / "figure.pdf")
        output_hashes = {str(p.relative_to(output)): sha256(p) for p in sorted(output.rglob("*"))
                         if p.is_file() and p.name not in ("result.json", "provenance.json")}
        result = {"status": "PASS", "schema_version": 1, "feature_id": "kubo-curvature-fullmesh",
                  "workflow_mode": "vasp_wavecar_calculation", "material": "1H-MoS2",
                  "settings": settings, "numerical_checks": {**numerical, "passed": True},
                  "band_export": band_metadata, "figure_conventions": plot_metadata,
                  "units": {"energy": "eV (unchanged WAVECAR zero)", "curvature": "Angstrom^2"},
                  "scope": "Native occupied1:18 Kubo trace; internal pairs excluded before division, external states19:26. " +
                           "Canonical momentum omits PAW augmentation and nonlocal/SOC "
                           "velocity terms. The 12 x 12 mesh and 26-band sum are tutorial settings, not a "
                           "material-convergence study.", "output_sha256": output_hashes}
        record.update(status="PASS", output_sha256=output_hashes,
                      finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        write_json(output / "provenance.json", record)
        result["provenance"] = record
        write_json(output / "result.json", result)
        print(f"PASS: map and path native Kubo outputs; {output / 'figure.png'}")
    except Exception as exc:
        record.update(status="FAIL", error=str(exc),
                      finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        write_json(output / "provenance.json", record)
        raise


if __name__ == "__main__":
    main()
