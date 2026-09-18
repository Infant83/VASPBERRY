#!/usr/bin/env python3
"""Generate a public synthetic WAVECAR and verify actual optical selectivity."""
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


def fixture(path):
    """Gamma, a=2*pi A, cutoff=5 eV: G=0,+x,-x,+y,-y,+z,-z."""
    recl = 128
    records = bytearray(5 * recl)
    def put(index, values, dtype):
        data = np.asarray(values, dtype=dtype).tobytes()
        assert len(data) <= recl
        records[index * recl:index * recl + len(data)] = data
    put(0, [recl, 1, 45200], "<f8")
    put(1, [1, 2, 5, *(2*np.pi*np.eye(3)).ravel()], "<f8")
    put(2, [7, 0, 0, 0, 0, 0, 2, 2, 0, 0], "<f8")
    valence = np.zeros(7, dtype=np.complex64)
    conduction = valence.copy()
    valence[[0, 1, 3]] = 1 / np.sqrt(3)
    w = np.exp(2j*np.pi/3)
    conduction[[0, 1, 3]] = np.array([1, w, w.conjugate()]) / np.sqrt(3)
    put(3, valence, "<c8")
    put(4, conduction, "<c8")
    path.write_bytes(records)
    return {"bytes": len(records), "sha256": sha(path), "k_fractional": [0, 0, 0],
            "record_bytes": recl, "bands": 2, "plane_waves": 7,
            "coefficient_norms": [float(np.vdot(c, c).real) for c in (valence, conduction)],
            "interband_overlap_absolute": float(abs(np.vdot(valence, conduction)))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True, help="New directory; existing paths are rejected")
    parser.add_argument("--binary", type=Path, default=REPO / "build/vaspberry-gfortran")
    parser.add_argument("--input", type=Path, default=HERE / "input.json", help="Example input JSON; fixed fixture, configurable theta samples")
    args = parser.parse_args()
    config = json.loads(args.input.read_text())
    if (set(config) != {"schema_version", "feature_id", "workflow_mode", "fixture", "theta_degrees", "phi_degrees"}
            or config["schema_version"] != 1 or config["feature_id"] != "circular-dichroism"
            or config["workflow_mode"] != "synthetic-production-fortran" or config["fixture"] != FIXTURE
            or config["phi_degrees"] != 0):
        parser.error("unsupported configuration: preserve the documented Gamma/scalar fixture, schema, and phi=0")
    angles = config["theta_degrees"]
    if (not isinstance(angles, list) or not 1 <= len(angles) <= 181
            or any(type(t) not in (int, float) or not np.isfinite(t) or not 0 <= t <= 180 for t in angles)
            or len(set(angles)) != len(angles)):
        parser.error("theta_degrees must contain 1..181 distinct finite angles between 0 and 180")
    binary = args.binary.resolve()
    if not binary.is_file():
        parser.error("build the production binary with 'make serial', or pass --binary")
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    report = {"status": "running", "example": "circular-dichroism", "commands": [],
              "schema_version": 1, "feature_id": "circular-dichroism", "workflow_mode": config["workflow_mode"],
              "input": config, "input_sha256": sha(args.input),
              "binary_sha256": sha(binary), "script_sha256": sha(Path(__file__)),
              "source_sha256": sha(REPO / "vaspberry.f"), "scope": "Synthetic Gamma-point bare-momentum optical selectivity; no material or absolute rate."}
    started = time.monotonic()
    try:
        report["fixture"] = fixture(out / "WAVECAR.synthetic")
        rows = []
        for index, theta in enumerate(angles):
            case = out / f"theta-{index:03d}"
            case.mkdir()
            command = [str(binary), "-f", "../WAVECAR.synthetic", "-s", "1",
                       "-kx", "1", "-ky", "1", "-ii", "1", "-if", "2",
                       "-cd", "1", "-theta", str(theta), "-phi", "0", "-kp", "1", "-o", "selectivity"]
            with (case / "stdout.log").open("w") as stdout, (case / "stderr.log").open("w") as stderr:
                result = subprocess.run(command, cwd=case, stdout=stdout, stderr=stderr,
                                        timeout=30, env={**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"})
            report["commands"].append({"cwd_relative_to_output": case.name, "argv": command, "returncode": result.returncode})
            if result.returncode:
                raise RuntimeError(f"production Fortran failed in {case.name}; inspect stderr.log")
            data = np.atleast_2d(np.loadtxt(case / "selectivity.dat"))
            if data.shape != (1, 7) or not np.all(np.isfinite(data)):
                raise AssertionError(f"unexpected result dimensions or nonfinite data: {data.shape}")
            dat_value = float(data[0, 3])
            printed = [line for line in (case / "stdout.log").read_text().splitlines()
                       if line.startswith("# IK, K(reci), SELECTIVITY(n) :")]
            if len(printed) != 1:
                raise AssertionError("expected one Gamma-point selectivity in stdout")
            actual = float(printed[0].split()[-1])
            cosine = np.cos(np.deg2rad(theta))
            expected = float(np.sqrt(3)*cosine/(1+cosine*cosine))
            rows.append([theta, actual, expected, actual - expected, dat_value, dat_value - expected])
        rows = np.asarray(rows)
        error = float(np.max(np.abs(rows[:, 3])))
        dat_error = float(np.max(np.abs(rows[:, 5])))
        if error > 1e-6:  # stdout has six decimal places; .dat has four.
            raise AssertionError(f"selectivity differs from the analytic oracle by {error}")
        if dat_error > 5.1e-5:
            raise AssertionError(f".dat value exceeds its four-decimal rounding tolerance: {dat_error}")
        with (out / "summary.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["theta_degrees", "fortran_stdout_selectivity", "analytic_selectivity", "stdout_error",
                             "fortran_dat_selectivity", "dat_error"])
            writer.writerows(rows)
        fig, ax = plt.subplots(figsize=(7, 4.3), constrained_layout=True)
        angles = np.linspace(0, 180, 361)
        c = np.cos(np.deg2rad(angles))
        ax.plot(angles, np.sqrt(3)*c/(1+c*c), color="#334155", label="Analytic model")
        ax.scatter(rows[:, 0], rows[:, 1], color="#dc6b32", s=30, zorder=3, label="Production Fortran")
        ax.axhline(0, color="#cbd5e1", linewidth=.8)
        ax.set(xlabel="Incident polar angle θ (degrees), φ = 0", ylabel="Circular optical selectivity η",
               title="Synthetic Gamma-point transition", xlim=(0, 180), ylim=(-1, 1))
        ax.legend(frameon=False)
        fig.savefig(out / "figure.png", dpi=170)
        plt.close(fig)
        if sha(out / "WAVECAR.synthetic") != report["fixture"]["sha256"]:
            raise AssertionError("input fixture changed")
        report.update(status="PASS", max_absolute_error=error, tolerance=1e-6, points=len(rows),
                      dat_max_absolute_error=dat_error, dat_tolerance=5.1e-5,
                      analytic_formula="sqrt(3)*cos(theta)/(1+cos(theta)^2)",
                      output_sha256={p.name: sha(p) for p in (out / "summary.csv", out / "figure.png")})
    except Exception as exc:
        report.update(status="FAILED", error=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        report["wall_seconds"] = time.monotonic() - started
        (out / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("status", "max_absolute_error", "points")}, indent=2))


if __name__ == "__main__":
    main()
