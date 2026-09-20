"""Shared input/provenance checks for the two real Bi WAVECAR tutorials."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
WAVECAR_SHA256 = "a8d81854f2efc561938e478dde1be29a17ccc95d5d37325c122b9d9e82fa0838"
WAVECAR_BYTES = 200421600
WAVECAR_URL = "https://media.githubusercontent.com/media/Infant83/VASPBERRY/a692e21482d24767859d02b7885dd70b234c6a2c/examples/Bi_Z2/WAVECAR"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def arguments(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--output-dir", required=True, type=Path, help="new directory for native VASPBERRY outputs, comparison and figure")
    parser.add_argument("--wavecar", type=Path, default=ROOT / "examples/Bi_Z2/WAVECAR", help="actual public Bi WAVECAR payload (the tutorial verifies its checksum)")
    parser.add_argument("--binary", type=Path, default=ROOT / "build/vaspberry-gfortran", help="current serial executable; build with make serial")
    parser.add_argument("--timeout", type=float, default=300, help="maximum runtime in seconds")
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error("output directory exists; choose a new path")
    if not args.wavecar.is_file():
        parser.error("missing Bi WAVECAR; run git lfs pull --include='examples/Bi_Z2/WAVECAR' or use the download command in README.md")
    with args.wavecar.open("rb") as handle:
        if handle.read(128).startswith(b"version https://git-lfs.github.com/spec/"):
            parser.error("WAVECAR is a Git LFS pointer, not VASP data; run git lfs pull --include='examples/Bi_Z2/WAVECAR' or download the public payload using README.md")
    if args.wavecar.stat().st_size != WAVECAR_BYTES or sha(args.wavecar) != WAVECAR_SHA256:
        parser.error("this reference tutorial requires the exact public Bi WAVECAR (SHA-256 " + WAVECAR_SHA256 + "); use the README commands with your own parameters for another system")
    if not args.binary.is_file():
        parser.error("missing executable; run make serial or provide --binary")
    if not np.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("timeout must be finite and positive")
    return args


def band_edges(wavecar):
    """Read energies from WAVECAR and cross-check the archived VASP EIGENVAL."""
    with wavecar.open("rb") as handle:
        recl, nspin, rtag = np.fromfile(handle, dtype="<f8", count=3)
        handle.seek(int(recl))
        header = np.fromfile(handle, dtype="<f8", count=12)
        nk, nb = map(int, header[:2])
        if (int(nspin), int(rtag), nk, nb) != (1, 45200, 144, 18):
            raise ValueError("unexpected WAVECAR header")
        energies = []
        for ik in range(nk):
            handle.seek(int(recl) * (2 + ik * (nb + 1)))
            values = np.fromfile(handle, dtype="<f8", count=4 + 3 * nb)
            energies.append(values[4::3])
    energies = np.array(energies)
    eigenval = ROOT / "examples/Bi_Z2/archive-2016-run/EIGENVAL"
    lines = eigenval.read_text().splitlines()
    if list(map(int, lines[5].split())) != [10, 144, 18]:
        raise ValueError("unexpected EIGENVAL header")
    archived = np.array([[float(lines[8 + ik * 20 + ib].split()[1]) for ib in range(nb)] for ik in range(nk)])
    if not np.allclose(energies, archived, atol=5.1e-7, rtol=0):
        raise ValueError("WAVECAR energies disagree with archived EIGENVAL")
    edges = np.column_stack((np.arange(1, nk + 1), energies[:, 9], energies[:, 10], energies[:, 10] - energies[:, 9]))
    checks = {"minimum_direct_gap_eV": float(edges[:, 3].min()),
              "sampled_global_gap_eV": float(edges[:, 2].min() - edges[:, 1].max())}
    if checks["minimum_direct_gap_eV"] <= 0 or checks["sampled_global_gap_eV"] <= 0:
        raise ValueError("selected occupied bundle is not separated on this mesh")
    return edges, checks


def wavecar_geometry(wavecar):
    """Read the actual mesh and lattice used by the native calculation."""
    with wavecar.open("rb") as handle:
        recl, nspin, rtag = np.fromfile(handle, dtype="<f8", count=3)
        handle.seek(int(recl))
        header = np.fromfile(handle, dtype="<f8", count=12)
        nk, nb = map(int, header[:2])
        if (int(nspin), int(rtag), nk, nb) != (1, 45200, 144, 18):
            raise ValueError("unexpected WAVECAR header")
        lattice = header[3:].reshape(3, 3)
        kpoints = []
        for ik in range(nk):
            handle.seek(int(recl) * (2 + ik * (nb + 1)))
            kpoints.append(np.fromfile(handle, dtype="<f8", count=4)[1:])
    return np.array(kpoints), 2 * np.pi * np.linalg.inv(lattice).T


def execute(args, feature_id, options):
    edges, gaps = band_edges(args.wavecar)
    out = args.output_dir.resolve()
    out.mkdir(parents=True)
    # Use a short, portable native command; never copy the 200 MB input into results.
    (out / "WAVECAR").symlink_to(args.wavecar.resolve())
    argv = [str(args.binary.resolve()), "-f", "WAVECAR", *options,
            "-kx", "12", "-ky", "12", "-s", "2", "-ii", "1", "-if", "10"]
    portable = ["build/vaspberry-gfortran", *argv[1:]]
    provenance = {"wavecar": "examples/Bi_Z2/WAVECAR", "wavecar_url": WAVECAR_URL,
                  "wavecar_bytes": WAVECAR_BYTES, "wavecar_sha256": WAVECAR_SHA256,
                  "eigenval": "examples/Bi_Z2/archive-2016-run/EIGENVAL",
                  "eigenval_sha256": sha(ROOT / "examples/Bi_Z2/archive-2016-run/EIGENVAL"),
                  "source_sha256": sha(ROOT / "vaspberry.f"), "binary_sha256": sha(args.binary),
                  "input_kind": "real_archived_VASP_5.4.1_SOC_WAVECAR",
                  "command": portable, "runner_sha256": sha(ROOT / "examples/features" / feature_id / "run.py"),
                  "helper_sha256": sha(__file__)}
    (out / "command.json").write_text(json.dumps(portable, indent=2) + "\n")
    start = time.monotonic()
    try:
        with (out / "fortran.log").open("w") as handle:
            completed = subprocess.run(argv, cwd=out, stdout=handle, stderr=subprocess.STDOUT, timeout=args.timeout)
        if completed.returncode:
            raise RuntimeError("VASPBERRY exited with status " + str(completed.returncode))
        log = (out / "fortran.log").read_text()
        if "VASPBERRY (Ver " + (ROOT / "VERSION").read_text().strip() + ")" not in log:
            raise RuntimeError("executable version does not match this checkout; rebuild")
    except Exception as error:
        (out / "result.json").write_text(json.dumps({"schema_version": 1, "feature_id": feature_id,
            "status": "FAILED", "workflow_mode": "vasp_wavecar_calculation", "error": str(error),
            "provenance": provenance}, indent=2) + "\n")
        raise
    finally:
        (out / "WAVECAR").unlink()
    provenance["elapsed_seconds"] = time.monotonic() - start
    np.savetxt(out / "band_edges.csv", edges, delimiter=",", header="k_index,band_10_eV,band_11_eV,direct_gap_eV", comments="", fmt=["%d", "%.12g", "%.12g", "%.12g"])
    return out, provenance, edges, gaps


def write_result(out, feature_id, provenance, summary, outputs, limitations):
    np.savetxt(out / "summary.csv", [[summary[key] for key in summary]], delimiter=",", header=",".join(summary), comments="", fmt="%.12g")
    outputs = [*outputs, "band_edges.csv", "summary.csv", "figure.png", "figure.pdf", "fortran.log", "command.json"]
    result = {"schema_version": 1, "feature_id": feature_id, "status": "PASS",
              "workflow_mode": "vasp_wavecar_calculation", "software_version": (ROOT / "VERSION").read_text().strip(),
              "summary": summary, "numerical_checks": {"status": "PASS", "actual_wavecar_sha256": "PASS",
              "wavecar_energies_match_archived_eigenval": "PASS", "positive_sampled_gap": "PASS"},
              "provenance": provenance, "limitations": limitations,
              "output_sha256": {name: sha(out / name) for name in outputs}}
    (out / "result.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": "PASS", "feature_id": feature_id, "summary": summary}))
