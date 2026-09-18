#!/usr/bin/env python3
"""Verify real/imaginary Gamma-point wavefunction output from production Fortran."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
FIXTURE = {"kind": "gamma_scalar_three_plane_waves", "cell_side_A": 2*np.pi,
           "cutoff_eV": 5, "bands": 2, "spinor": False, "k_fractional": [0, 0, 0]}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixture(directory):
    """Self-authored scalar coefficients; no VASP executable or material input."""
    recl = 128
    records = bytearray(5 * recl)
    def put(index, values, dtype):
        data = np.asarray(values, dtype=dtype).tobytes()
        assert len(data) <= recl
        records[index * recl:index * recl + len(data)] = data
    lattice = 2*np.pi*np.eye(3)
    put(0, [recl, 1, 45200], "<f8")
    put(1, [1, 2, 5, *lattice.ravel()], "<f8")
    put(2, [7, 0, 0, 0, 0, 0, 2, 2, 0, 0], "<f8")
    valence = np.zeros(7, dtype=np.complex64)
    conduction = valence.copy()
    valence[[0, 1, 3]] = 1 / np.sqrt(3)
    w = np.exp(2j*np.pi/3)
    conduction[[0, 1, 3]] = np.array([1, w, w.conjugate()]) / np.sqrt(3)
    put(3, valence, "<c8")
    put(4, conduction, "<c8")
    (directory / "WAVECAR.synthetic").write_bytes(records)
    (directory / "POSCAR").write_text("Synthetic visualization cell; no material calculation\n1.0\n" +
        "\n".join(" ".join(f"{x:.17g}" for x in row) for row in lattice) + "\nH\n1\nDirect\n0 0 0\n")
    (directory / "EIGENVAL").write_text("1 1 1 1\n0 0 0 0 0\n0\nSynthetic\nVisualization only\n2 1 2\n\n0 0 0 1\n1 0 2\n2 2 0\n")
    return {p.name: sha(p) for p in (directory / "WAVECAR.synthetic", directory / "POSCAR", directory / "EIGENVAL")}


def read_grid(path, grid):
    lines = path.read_text().splitlines()
    matches = [i for i, line in enumerate(lines) if line.split() == [str(n) for n in grid]]
    if len(matches) != 1:
        raise AssertionError(f"expected one scalar grid in {path.name}")
    values = np.fromstring(" ".join(lines[matches[0]+1:]), sep=" ")
    if values.size != np.prod(grid) or not np.all(np.isfinite(values)):
        raise AssertionError(f"invalid values in {path.name}")
    return values.reshape(tuple(reversed(grid)))  # x varies fastest in the file.


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True, help="New directory; existing paths are rejected")
    parser.add_argument("--binary", type=Path, default=REPO / "build/vaspberry-gfortran")
    parser.add_argument("--input", type=Path, default=HERE / "input.json", help="Example input JSON; fixed Gamma fixture, configurable grid")
    args = parser.parse_args()
    config = json.loads(args.input.read_text())
    if (set(config) != {"schema_version", "feature_id", "workflow_mode", "fixture", "grid"}
            or config["schema_version"] != 1 or config["feature_id"] != "wavefunction"
            or config["workflow_mode"] != "synthetic-production-fortran" or config["fixture"] != FIXTURE):
        parser.error("unsupported configuration: preserve the documented Gamma/scalar fixture and schema")
    grid = config["grid"]
    if (not isinstance(grid, list) or len(grid) != 3 or any(type(n) is not int or not 5 <= n <= 64 for n in grid)
            or grid[0] != grid[1]):
        parser.error("grid must contain three integers from 5 to 64, with nx=ny for the diagonal summary")
    binary = args.binary.resolve()
    if not binary.is_file():
        parser.error("build the production binary with 'make serial', or pass --binary")
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    report = {"status": "running", "example": "wavefunction", "commands": [],
              "schema_version": 1, "feature_id": "wavefunction", "workflow_mode": config["workflow_mode"],
              "input": config, "config_sha256": sha(args.input),
              "binary_sha256": sha(binary), "script_sha256": sha(Path(__file__)),
              "source_sha256": sha(REPO / "vaspberry.f"),
              "scope": "Synthetic scalar Gamma-point wavefunction; raw output is amplitude, not charge density."}
    started = time.monotonic()
    try:
        report["input_sha256"] = fixture(out)
        command = [str(binary), "-f", "WAVECAR.synthetic", "-s", "1", "-kx", "1", "-ky", "1",
                   "-wf", "1", "-k", "1", "-ng", ",".join(map(str, grid)), "-im", "1"]
        with (out / "stdout.log").open("w") as stdout, (out / "stderr.log").open("w") as stderr:
            result = subprocess.run(command, cwd=out, stdout=stdout, stderr=stderr, timeout=30,
                                    env={**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"})
        report["commands"].append({"cwd_relative_to_output": ".", "argv": command, "returncode": result.returncode})
        if result.returncode:
            raise RuntimeError("production Fortran failed; inspect stderr.log")
        real = read_grid(out / "PARCHG-W-K001-E001-SPIN1", grid)
        imag = read_grid(out / "PARCHG-W-K001-E001-IM-SPIN1", grid)
        volume = (2*np.pi)**3
        actual_u = (real + 1j*imag) / np.sqrt(volume)
        x = np.arange(grid[0]) / grid[0]
        y = np.arange(grid[1]) / grid[1]
        expected_u = (1 + np.exp(2j*np.pi*x)[None, :] + np.exp(2j*np.pi*y)[:, None]) / np.sqrt(3)
        error = float(np.max(np.abs(actual_u - expected_u[None, :, :])))
        if error > 2e-6:
            raise AssertionError(f"wavefunction differs from the analytic oracle by {error}")
        mean_norm = float(np.mean(np.abs(actual_u)**2))
        if abs(mean_norm - 1) > 2e-6:
            raise AssertionError("dimensionless u has incorrect average norm")
        with (out / "summary.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["x_equals_y_fractional", "fortran_Re_u", "fortran_Im_u", "analytic_Re_u", "analytic_Im_u"])
            for i, coordinate in enumerate(x):
                a, e = actual_u[0, i, i], expected_u[i, i]
                writer.writerow([coordinate, a.real, a.imag, e.real, e.imag])
        fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.6), constrained_layout=True)
        for ax, values, title, cmap in zip(axes,
                [actual_u[0].real, actual_u[0].imag, np.abs(actual_u[0])**2],
                ["Re u(x,y)", "Im u(x,y)", "|u(x,y)|²"], ["RdBu_r", "RdBu_r", "viridis"]):
            artist = ax.imshow(values, origin="lower", extent=(0, 1, 0, 1), interpolation="nearest", cmap=cmap)
            ax.set(title=title, xlabel="x / a", ylabel="y / a")
            fig.colorbar(artist, ax=ax, shrink=.78)
        fig.suptitle("Production Fortran: synthetic Gamma state, z = 0; u = raw / √V")
        fig.savefig(out / "figure.png", dpi=170)
        plt.close(fig)
        for name, digest in report["input_sha256"].items():
            if sha(out / name) != digest:
                raise AssertionError(f"input changed: {name}")
        report.update(status="PASS", max_absolute_error=error, tolerance=2e-6, grid=list(grid),
                      cell_volume_A3=volume, mean_abs_u_squared=mean_norm,
                      analytic_u="(1+exp(2*pi*i*x)+exp(2*pi*i*y))/sqrt(3)",
                      raw_output="sqrt(V)*u; normalized physical psi=raw/V", k_fractional=[0, 0, 0],
                      output_sha256={p.name: sha(p) for p in (out / "summary.csv", out / "figure.png")})
    except Exception as exc:
        report.update(status="FAILED", error=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        report["wall_seconds"] = time.monotonic() - started
        (out / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("status", "max_absolute_error", "mean_abs_u_squared")}, indent=2))


if __name__ == "__main__":
    main()
