#!/usr/bin/env python3
"""Run VASPBERRY on an actual MoS2 Gamma-point spinor state and plot its density."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PUBLIC_WAVECAR_SHA256 = "33f8546512856d6c04ad0a80454b18ac9b60e2af4b98f2f85b376ec49b1a8d9f"
sys.path.insert(0, str(REPO / "tools"))
from wavecar_fukui import Wavecar


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_grid(path, grid):
    """The SOC output contains two scalar grids, x fastest, one per spinor component."""
    lines = path.read_text().splitlines()
    starts = [i for i, line in enumerate(lines) if line.split() == [str(n) for n in grid]]
    if len(starts) != 2:
        raise AssertionError(f"expected two spinor grids in {path.name}")
    blocks = []
    for block, start in enumerate(starts):
        end = starts[block+1] if block+1 < len(starts) else len(lines)
        values = np.fromstring(" ".join(lines[start+1:end]), sep=" ")
        if values.size != np.prod(grid) or not np.isfinite(values).all():
            raise AssertionError(f"invalid scalar grid in {path.name}")
        blocks.append(values.reshape(tuple(reversed(grid))))
    return np.asarray(blocks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, default=REPO / "examples/1H-MoS2/KPATH/2.band")
    parser.add_argument("--wavecar", type=Path, help="Override WAVECAR; POSCAR and EIGENVAL still come from --input-dir")
    parser.add_argument("--binary", type=Path, default=REPO / "build/vaspberry-gfortran")
    parser.add_argument("--band", type=int, default=18)
    parser.add_argument("--k-index", type=int, default=24, help="One-based Gamma index in WAVECAR")
    parser.add_argument("--postprocess-only", action="store_true", help="Read an existing raw calculation directory; do not run Fortran")
    args = parser.parse_args()
    source = (args.wavecar or args.input_dir / "WAVECAR").resolve()
    material = "1H-MoS2" if sha(source) == PUBLIC_WAVECAR_SHA256 else "User-supplied VASP data"
    wavecar = Wavecar(source, spinor_components=2)
    if not 1 <= args.k_index <= wavecar.header.nkpoints or not 1 <= args.band <= wavecar.header.nbands:
        parser.error("band or k index is outside WAVECAR")
    ik = args.k_index-1
    if not np.allclose(wavecar.kpoints[ik], 0, atol=1e-10):
        parser.error("select a Gamma point (k=0) for this tutorial")
    inputs = {"WAVECAR": source, **{name: (args.input_dir / name).resolve() for name in ("POSCAR", "EIGENVAL")}}
    if not all(p.is_file() for p in inputs.values()):
        parser.error("WAVECAR, POSCAR and EIGENVAL from the same VASP calculation are required")
    poscar = inputs["POSCAR"].read_text().splitlines()
    lattice = np.asarray([[float(x) for x in line.split()] for line in poscar[2:5]]) * float(poscar[1])
    if not np.allclose(lattice, wavecar.header.lattice, rtol=0, atol=1e-8):
        parser.error("POSCAR lattice does not match WAVECAR")
    if not (np.allclose(lattice[:2, 2], 0, atol=1e-8) and np.allclose(lattice[2, :2], 0, atol=1e-8) and lattice[2, 2] > 0):
        parser.error("this z-profile plotter needs an xy-plane slab with a positive z-directed third cell vector")
    grid = [24, 24, 64]
    out = args.output_dir.resolve()
    if args.postprocess_only:
        if not out.is_dir() or (out / "result.json").exists():
            parser.error("postprocessing needs an existing raw directory without result.json")
    else:
        if not args.binary.resolve().is_file():
            parser.error("run 'make serial' or supply --binary")
        out.mkdir(parents=True, exist_ok=False)
    report = {"schema_version": 1, "feature_id": "wavefunction", "status": "RUNNING",
              "workflow_mode": "vasp_wavecar_calculation", "commands": [],
              "execution_mode": "postprocess_existing_outputs" if args.postprocess_only else "calculate_and_postprocess",
              "provenance": {"material": material, "input_paths": {name: str(p) for name, p in inputs.items()},
                             "input_sha256": {name: sha(p) for name, p in inputs.items()},
                             "script_sha256": sha(Path(__file__)),
                             "local_analysis_source_sha256": {"run.py": sha(Path(__file__)),
                                                              "tools/wavecar_fukui.py": sha(REPO / "tools/wavecar_fukui.py")},
                             "analysis_checkout_version": (REPO / "VERSION").read_text().strip(),
                             "producer_revision": "unknown",
                             "input_origin": "Existing public VASP output; no new VASP calculation"},
              "parameters": {"spinor_components": 2, "band": args.band, "k_index": args.k_index, "grid": grid},
              "scope": "Stored Gamma-point pseudo-wavefunction spinor, with both components retained; no PAW augmentation density."}
    start = time.monotonic()
    try:
        if not args.postprocess_only:
            report["provenance"].update(source_sha256=sha(REPO / "vaspberry.f"),
                                        source_sha256_scope="local Fortran source at calculation time",
                                        producer_revision="source-sha256:" + sha(REPO / "vaspberry.f"),
                                        vaspberry_version=(REPO / "VERSION").read_text().strip())
            (out / "WAVECAR").symlink_to(source)
            for name in ("POSCAR", "EIGENVAL"):
                shutil.copyfile(inputs[name], out / name)
            command = [str(args.binary.resolve()), "-f", "WAVECAR", "-s", "2", "-kx", str(wavecar.header.nkpoints),
                       "-ky", "1", "-wf", str(args.band), "-k", str(args.k_index), "-ng", ",".join(map(str, grid)), "-im", "1"]
            with (out / "stdout.log").open("w") as stdout, (out / "stderr.log").open("w") as stderr:
                run = subprocess.run(command, cwd=out, stdout=stdout, stderr=stderr, timeout=120,
                                     env={**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"})
            report["commands"].append({"argv": command, "cwd_relative_to_output": ".", "returncode": run.returncode})
            report["provenance"]["binary_sha256"] = sha(args.binary.resolve())
            if run.returncode:
                raise RuntimeError("VASPBERRY failed; see stderr.log")
        basename = f"PARCHG-W-K{args.k_index:03d}-E{args.band:03d}"
        real = read_grid(out / (basename + "-SPIN1"), grid)
        imag = read_grid(out / (basename + "-IM-SPIN1"), grid)
        volume = abs(float(np.linalg.det(lattice)))
        psi = (real + 1j*imag) / volume
        density = np.sum(abs(psi)**2, axis=0)
        coeff = wavecar.coefficients(ik, [args.band])[0].astype(complex)
        g = wavecar.g_vectors(ik)
        if any(np.ptp(g[:, axis]) >= grid[axis] for axis in range(3)):
            raise ValueError("tutorial grid is too small for this WAVECAR cutoff; adapt -ng and the postprocessor")
        fft_grid = np.zeros((2, *grid[::-1]), complex)
        for spinor in range(2):
            fft_grid[spinor, g[:, 2] % grid[2], g[:, 1] % grid[1], g[:, 0] % grid[0]] = coeff[spinor]
        expected = np.fft.ifftn(fft_grid, axes=(-3, -2, -1)) * np.prod(grid) / np.sqrt(volume)
        error = float(np.max(abs(psi - expected)))
        coefficient_norm = float(np.sum(abs(coeff)**2))
        density_integral = float(np.mean(density)*volume)
        if error > 1e-6 or abs(density_integral-coefficient_norm) > 1e-6:
            raise AssertionError("raw Fortran grids disagree with WAVECAR coefficients beyond output precision")
        with (out / "summary.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["z_fractional", "z_A", "mean_spinor1_density_inv_A3", "mean_spinor2_density_inv_A3", "mean_total_density_inv_A3"])
            means = np.mean(abs(psi)**2, axis=(2, 3))
            for iz in range(grid[2]):
                z = iz / grid[2]
                writer.writerow([z, z*lattice[2, 2], means[0, iz], means[1, iz], means[:, iz].sum()])
        reference = HERE / "reference/summary.csv"
        reference_info = HERE / "reference/result.json"
        comparison = {"status": "NOT_AVAILABLE"}
        if reference.is_file() and reference_info.is_file():
            saved = json.loads(reference_info.read_text())
            if (report["provenance"]["input_sha256"] == saved["provenance"]["input_sha256"] and
                    report["parameters"] == saved["parameters"]):
                actual = np.loadtxt(out / "summary.csv", delimiter=",", skiprows=1)
                expected_rows = np.loadtxt(reference, delimiter=",", skiprows=1)
                np.testing.assert_allclose(actual, expected_rows, rtol=1e-6, atol=1e-9)
                comparison = {"status": "PASS", "relative_tolerance": 1e-6, "absolute_tolerance": 1e-9,
                              "reference": "reference/summary.csv"}
            else:
                comparison = {"status": "NOT_APPLICABLE", "reason": "different input/state; no MoS2 numerical match asserted"}
        fig, axes = plt.subplots(1, 2, figsize=(9.6, 4), constrained_layout=True)
        projection = density.mean(axis=0)*lattice[2, 2]
        artist = axes[0].imshow(projection, origin="lower", extent=(0, 1, 0, 1), interpolation="nearest", cmap="viridis")
        axes[0].set(xlabel="Fraction along a₁", ylabel="Fraction along a₂", title="Density integrated along z")
        fig.colorbar(artist, ax=axes[0], label="Pseudo-density (Å⁻²)", shrink=.82)
        z = np.arange(grid[2])*lattice[2, 2]/grid[2]
        axes[1].plot(z, density.mean(axis=(1, 2)), color="#265d92", linewidth=2)
        axes[1].set(xlabel="z (Å)", ylabel="Plane-averaged pseudo-density (Å⁻³)", title="Both spinor components included")
        axes[1].grid(alpha=.2)
        fig.suptitle(f"{material.replace('MoS2', 'MoS₂')} · VASPBERRY -wf {args.band} -k {args.k_index} · Γ")
        fig.savefig(out / "figure.png", dpi=170)
        plt.close(fig)
        report.update(status="PASS", grid=grid, k_fractional=wavecar.kpoints[ik].tolist(),
                      state_energy_eV=float(wavecar.energies[ik, args.band-1]), cell_volume_A3=volume,
                      integrated_pseudo_density=density_integral, coefficient_norm=coefficient_norm,
                      fft_comparison_max_error_inv_A_three_halves=error, fft_absolute_tolerance=1e-6,
                      raw_output_convention="raw = sqrt(V) * sum_G c_G exp(i G.r); psi_pseudo = raw / V",
                      reference_comparison=comparison,
                      output_sha256={p.name: sha(p) for p in (out / "summary.csv", out / "figure.png")})
        if {name: sha(p) for name, p in inputs.items()} != report["provenance"]["input_sha256"]:
            raise AssertionError("VASP inputs changed during calculation")
    except Exception as exc:
        report.update(status="FAILED", error=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        report["wall_seconds"] = time.monotonic() - start
        (out / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("status", "integrated_pseudo_density", "coefficient_norm")}, indent=2))


if __name__ == "__main__":
    main()
