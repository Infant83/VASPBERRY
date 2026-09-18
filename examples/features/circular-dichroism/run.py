#!/usr/bin/env python3
"""Run VASPBERRY optical spectra on the supplied VASP MoS2 WAVECAR."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, default=REPO / "examples/1H-MoS2/KPATH/2.band")
    parser.add_argument("--wavecar", type=Path, help="Override WAVECAR; the documented band window must remain appropriate")
    parser.add_argument("--binary", type=Path, default=REPO / "build/vaspberry-gfortran")
    parser.add_argument("--postprocess-only", action="store_true", help="Read an existing raw calculation directory; do not run Fortran")
    args = parser.parse_args()
    source = (args.wavecar or args.input_dir / "WAVECAR").resolve()
    input_hash = sha(source)
    material = "1H-MoS2" if input_hash == PUBLIC_WAVECAR_SHA256 else "User-supplied VASP data"
    wavecar = Wavecar(source, spinor_components=2)
    if wavecar.header.nbands < 20 or wavecar.header.ispin != 1:
        parser.error("this tutorial requires an ISPIN=1 spinor WAVECAR with at least 20 bands")
    if not (np.allclose(wavecar.occupations[:, :18], 1, atol=1e-6) and
            np.allclose(wavecar.occupations[:, 18:], 0, atol=1e-6)):
        parser.error("this tutorial's 18 occupied bands do not match the WAVECAR; adapt the documented CLI to your system")
    out = args.output_dir.resolve()
    if args.postprocess_only:
        if not out.is_dir() or (out / "result.json").exists():
            parser.error("postprocessing needs an existing raw directory without result.json")
    else:
        if not args.binary.resolve().is_file():
            parser.error("run 'make serial' or supply --binary")
        out.mkdir(parents=True, exist_ok=False)
    report = {"schema_version": 1, "feature_id": "circular-dichroism", "status": "RUNNING",
              "workflow_mode": "vasp_wavecar_calculation", "commands": [],
              "execution_mode": "postprocess_existing_outputs" if args.postprocess_only else "calculate_and_postprocess",
              "provenance": {"material": material, "input_path": str(source), "input_sha256": {"WAVECAR": input_hash},
                             "script_sha256": sha(Path(__file__)),
                             "local_analysis_source_sha256": {"run.py": sha(Path(__file__)),
                                                              "tools/wavecar_fukui.py": sha(REPO / "tools/wavecar_fukui.py")},
                             "analysis_checkout_version": (REPO / "VERSION").read_text().strip(),
                             "producer_revision": "unknown",
                             "input_origin": "Existing public VASP output; no new VASP calculation"},
              "parameters": {"spinor_components": 2, "occupied_bands": [1, 18], "empty_bands": [19, 20],
                             "photon_energy_eV": [1, 3], "energy_points": 201, "sigma_eV": 0.05,
                             "theta_degrees": 0, "phi_degrees": 0},
              "scope": "Bare-momentum k-resolved optical spectra on a line; arbitrary units, not a full-BZ absorption integral."}
    start = time.monotonic()
    try:
        if not args.postprocess_only:
            report["provenance"].update(source_sha256=sha(REPO / "vaspberry.f"),
                                        source_sha256_scope="local Fortran source at calculation time",
                                        producer_revision="source-sha256:" + sha(REPO / "vaspberry.f"),
                                        vaspberry_version=(REPO / "VERSION").read_text().strip())
            (out / "WAVECAR").symlink_to(source)
            command = [str(args.binary.resolve()), "-f", "WAVECAR", "-s", "2", "-kx", str(wavecar.header.nkpoints),
                       "-ky", "1", "-cd", "2", "-if", "20", "-ien", "1", "-fen", "3", "-nediv", "201",
                       "-sigma", "0.05", "-theta", "0", "-phi", "0", "-o", "optical"]
            with (out / "stdout.log").open("w") as stdout, (out / "stderr.log").open("w") as stderr:
                run = subprocess.run(command, cwd=out, stdout=stdout, stderr=stderr, timeout=120,
                                     env={**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"})
            report["commands"].append({"argv": command, "cwd_relative_to_output": ".", "returncode": run.returncode})
            report["provenance"]["binary_sha256"] = sha(args.binary.resolve())
            if run.returncode:
                raise RuntimeError("VASPBERRY failed; see stderr.log")
        channels = []
        sum_errors = []
        for channel in ("LEFT", "RIGHT"):
            records = [np.loadtxt(out / f"CIRC_DICHROISM_W.optical_{channel}_KP{i+1}.dat")
                       for i in range(wavecar.header.nkpoints)]
            data = np.asarray(records)
            if data.shape != (wavecar.header.nkpoints, 201, 8) or not np.isfinite(data).all():
                raise AssertionError("unexpected optical output dimensions or nonfinite data")
            if not np.allclose(data[:, :, 5:8], wavecar.kpoints[:, None, :], atol=5.1e-7):
                raise AssertionError("output k points do not match WAVECAR")
            energy = data[0, :, 3]
            if not np.allclose(data[:, :, 3], np.linspace(1, 3, 201)[None, :], atol=5.1e-5):
                raise AssertionError("unexpected photon-energy grid")
            spectrum = data[:, :, 4]
            if np.min(spectrum) < 0:
                raise AssertionError("negative optical intensity")
            total = np.loadtxt(out / f"CIRC_DICHROISM_W.optical_{channel}.dat")
            error = float(np.max(abs(total[:, 1] - spectrum.sum(axis=0))))
            if error > (wavecar.header.nkpoints + 1) * 5.1e-5:
                raise AssertionError("total file and per-k files disagree beyond output rounding")
            sum_errors.append(error)
            channels.append(spectrum)
        left, right = channels
        intensity = left + right
        threshold = max(1e-4, 1e-3 * float(intensity.max()))
        visible = intensity > threshold
        selectivity = np.zeros_like(intensity)
        np.divide(left - right, intensity, out=selectivity, where=visible)
        distances = np.r_[0, np.cumsum(np.linalg.norm(np.diff(wavecar.kpoints, axis=0) @ wavecar.header.reciprocal, axis=1))]
        rows = []
        for ik in range(wavecar.header.nkpoints):
            for ie, photon in enumerate(energy):
                rows.append([ik+1, *wavecar.kpoints[ik], distances[ik], photon, left[ik, ie], right[ik, ie],
                             selectivity[ik, ie] if visible[ik, ie] else "", int(visible[ik, ie])])
        with (out / "summary.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["k_index", "kx_fractional", "ky_fractional", "kz_fractional", "path_distance_inv_A",
                             "photon_energy_eV", "left_spectrum_au", "right_spectrum_au", "selectivity", "selectivity_is_resolved"])
            writer.writerows(rows)
        reference = HERE / "reference/summary.csv"
        reference_info = HERE / "reference/result.json"
        comparison = {"status": "NOT_AVAILABLE"}
        if reference.is_file() and reference_info.is_file():
            expected_hash = json.loads(reference_info.read_text())["provenance"]["input_sha256"]["WAVECAR"]
            if sha(source) == expected_hash:
                expected = np.genfromtxt(reference, delimiter=",", skip_header=1)
                actual = np.genfromtxt(out / "summary.csv", delimiter=",", skip_header=1)
                np.testing.assert_allclose(actual, expected, rtol=0, atol=1.1e-4, equal_nan=True)
                comparison = {"status": "PASS", "absolute_tolerance": 1.1e-4, "reference": "reference/summary.csv"}
            else:
                comparison = {"status": "NOT_APPLICABLE", "reason": "different WAVECAR; no MoS2 numerical match asserted"}
        fig, axes = plt.subplots(1, 3, figsize=(11.2, 4), constrained_layout=True)
        extent = (.5, wavecar.header.nkpoints+.5, float(energy[0]), float(energy[-1]))
        for ax, values, title, cmap in zip(axes, [left, right, np.ma.array(selectivity, mask=~visible)],
                                         ["Left channel", "Right channel", "(Left − Right) / (Left + Right)"],
                                         ["magma", "magma", "RdBu_r"]):
            limits = {"vmin": -1, "vmax": 1} if ax is axes[2] else {"vmin": 0, "vmax": float(intensity.max())}
            artist = ax.imshow(values.T, origin="lower", aspect="auto", interpolation="nearest", extent=extent, cmap=cmap, **limits)
            ax.set(title=title, xlabel="WAVECAR k-point index", ylabel="Photon energy (eV)")
            if input_hash == PUBLIC_WAVECAR_SHA256:
                ax.set_xticks([1, 24.5, 48], ["K (1)", "Γ (24/25)", "K′ (48)"])
            fig.colorbar(artist, ax=ax, shrink=.8, label="Selectivity" if ax is axes[2] else "Intensity (a.u.)")
        fig.suptitle(f"{material.replace('MoS2', 'MoS₂')} · VASP WAVECAR → VASPBERRY -cd 2 · normal incidence")
        fig.savefig(out / "figure.png", dpi=170)
        plt.close(fig)
        peak = int(np.argmax(intensity[0]))
        report.update(status="PASS", points=wavecar.header.nkpoints, spectrum_rows=len(rows),
                      reference_comparison=comparison, total_file_rounding_max_errors=sum_errors,
                      selectivity_intensity_threshold_au=threshold, weak_intensity_policy="blank CSV value and masked figure pixel",
                      first_k_peak_energy_eV=float(energy[peak]), first_k_peak_left=float(left[0, peak]),
                      first_k_peak_right=float(right[0, peak]), first_k_peak_selectivity=float(selectivity[0, peak]),
                      last_k_peak_selectivity=float(selectivity[-1, np.argmax(intensity[-1])]),
                      output_sha256={p.name: sha(p) for p in (out / "summary.csv", out / "figure.png")})
        if sha(source) != report["provenance"]["input_sha256"]["WAVECAR"]:
            raise AssertionError("WAVECAR changed during calculation")
    except Exception as exc:
        report.update(status="FAILED", error=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        report["wall_seconds"] = time.monotonic() - start
        (out / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("status", "spectrum_rows", "first_k_peak_selectivity", "last_k_peak_selectivity")}, indent=2))


if __name__ == "__main__":
    main()
