#!/usr/bin/env python3
"""Compare real MoS2 native serial/MPI Kubo exports; keep reproducible receipts.

The pair-to-bundle comparison checks the normalization and external-transition
sum independently of the native bundle loop. It is not a PAW accuracy test.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
INPUT_SHA256 = "33f8546512856d6c04ad0a80454b18ac9b60e2af4b98f2f85b376ec49b1a8d9f"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_export(path, kind, version):
    metadata, lines = {}, []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            key, separator, value = line[1:].strip().partition("=")
            if separator:
                metadata[key] = value
        elif line.strip():
            lines.append(line)
    expected = {
        "schema": f"VASPBERRY_BARE_MOMENTUM_KUBO_{kind}_V1",
        "normalization": "STANDARD_MINUS_TWO_IM",
        "result_status": "PASS",
        "vaspberry_version": version,
    }
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError(f"{path}: expected {key}={value}")
    rows = [{key: float(value) for key, value in row.items()}
            for row in csv.DictReader(lines)]
    if not rows or not all(math.isfinite(value) for row in rows for value in row.values()):
        raise ValueError(f"{path}: empty or non-finite export")
    return rows


def compare_rows(left, right, *, label, atol=2e-10, rtol=2e-11):
    if len(left) != len(right) or not left:
        raise ValueError(f"{label}: row count differs or is empty")
    maximum = 0.0
    for row_index, (a, b) in enumerate(zip(left, right), 1):
        if a.keys() != b.keys():
            raise ValueError(f"{label}: columns differ")
        for key in a:
            if not (math.isfinite(a[key]) and math.isfinite(b[key])):
                raise ValueError(f"{label}: non-finite {key} at row {row_index}")
            delta = abs(a[key] - b[key])
            maximum = max(maximum, delta)
            if not math.isclose(a[key], b[key], abs_tol=atol, rel_tol=rtol):
                raise ValueError(f"{label}: {key} differs at row {row_index}: {a[key]} != {b[key]}")
    return {"rows": len(left), "max_absolute_difference": maximum, "atol": atol, "rtol": rtol}


def pair_bundle_check(bundle, pairs):
    expected_keys = {(1, k) for k in range(1, 49)}
    bundle_keys = [(r["spin"], r["k_index"]) for r in bundle]
    if len(bundle) != 48 or set(bundle_keys) != expected_keys:
        raise ValueError("expected 48 distinct MoS2 bundle k points")
    expected_pairs = {(n, m) for n in range(1, 33) for m in range(n + 1, 33)}
    seen = {key: set() for key in expected_keys}
    contributions = {key: [] for key in expected_keys}
    gaps = {key: [] for key in expected_keys}
    coordinates = {(r["spin"], r["k_index"]): tuple(r[f"k{axis}_frac"] for axis in "xyz")
                   for r in bundle}
    for row in pairs:
        key = (row["spin"], row["k_index"])
        pair = (row["n_band"], row["m_band"])
        if key not in seen or pair not in expected_pairs or pair in seen[key]:
            raise ValueError("invalid or duplicate MoS2 unordered pair")
        seen[key].add(pair)
        if coordinates[key] != tuple(row[f"k{axis}_frac"] for axis in "xyz"):
            raise ValueError("pair/bundle coordinates differ")
        gap = abs(row["energy_n_eV"] - row["energy_m_eV"])
        if not math.isclose(gap, row["gap_eV"], abs_tol=1e-12, rel_tol=1e-12):
            raise ValueError("pair energy/gap mismatch")
        if pair[0] <= 18 < pair[1]:
            if gap <= 1e-5:
                raise ValueError("MoS2 occupied bundle has a closed external gap")
            contributions[key].append(row["numerator_xy_eV2_A2"] / gap**2)
            gaps[key].append(gap)
    if any(values != expected_pairs for values in seen.values()):
        raise ValueError("incomplete MoS2 unordered pair coverage")
    reconstructed = []
    for row in bundle:
        key = row["spin"], row["k_index"]
        reconstructed.append({**row, "omega_z_A2": math.fsum(contributions[key]),
                              "min_external_gap_eV": min(gaps[key])})
    return compare_rows(bundle, reconstructed, label="external-pair sum versus bundle")


def run_recorded(command, directory):
    directory.mkdir()
    receipt = {"command": command, "cwd": str(directory), "status": "RUNNING"}
    record_path = directory / "command.json"
    record_path.write_text(json.dumps(receipt, indent=2) + "\n")
    started = time.monotonic()
    try:
        with (directory / "stdout.log").open("wb") as out, (directory / "stderr.log").open("wb") as err:
            result = subprocess.run(command, cwd=directory, stdout=out, stderr=err, timeout=180)
        receipt["exit_code"] = result.returncode
        if result.returncode:
            raise RuntimeError(f"native command failed ({result.returncode}); see {directory}")
        receipt["status"] = "PASS"
    except Exception as exc:
        receipt.update(status="FAIL", error=str(exc))
        raise
    finally:
        receipt["elapsed_seconds"] = time.monotonic() - started
        record_path.write_text(json.dumps(receipt, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--serial", type=Path, required=True)
    parser.add_argument("--mpi", type=Path, required=True)
    parser.add_argument("--mpiexec", default="mpiexec.hydra")
    parser.add_argument("--mpi-flag", action="append", default=[])
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    serial, mpi = args.serial.resolve(), args.mpi.resolve()
    launcher = shutil.which(args.mpiexec)
    if not launcher:
        parser.error(f"MPI launcher not found: {args.mpiexec}")
    wavecar = ROOT / "examples/1H-MoS2/KPATH/2.band/WAVECAR"
    if sha256(wavecar) != INPUT_SHA256:
        raise ValueError("public MoS2 WAVECAR checksum mismatch")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    summary = {"status": "RUNNING", "input_sha256": INPUT_SHA256,
               "serial_sha256": sha256(serial), "mpi_sha256": sha256(mpi),
               "mpi_launcher": launcher, "ranks": 2, "checks": {}}
    try:
        tables = {}
        version = (ROOT / "VERSION").read_text().strip()
        for mode, prefix in (("serial", [str(serial)]),
                             ("mpi", [launcher, *args.mpi_flag, "-n", "2", str(mpi)])):
            tables[mode] = {}
            for task, kind, extra in (("kubo", "BUNDLE", ["--bands", "1:18", "--bundle", "1"]),
                                      ("kubo-pairs", "PAIRS", [])):
                directory = output / f"{mode}-{task}"
                destination = directory / f"{kind}.csv"
                option = "--curvature-csv" if task == "kubo" else "--pairs-csv"
                command = prefix + ["--task", task, "--wavecar", str(wavecar), "--spinor", "2",
                                    *extra, option, str(destination)]
                run_recorded(command, directory)
                tables[mode][kind] = read_export(destination, kind, version)
            summary["checks"][f"{mode}_pair_to_bundle"] = pair_bundle_check(
                tables[mode]["BUNDLE"], tables[mode]["PAIRS"])
        for kind in ("BUNDLE", "PAIRS"):
            summary["checks"][f"serial_vs_mpi_{kind.lower()}"] = compare_rows(
                tables["serial"][kind], tables["mpi"][kind], label=f"serial/MPI {kind}")
        summary["status"] = "PASS"
    except Exception as exc:
        summary.update(status="FAIL", error=str(exc))
        raise
    finally:
        summary["outputs"] = {str(path.relative_to(output)): sha256(path)
                              for path in sorted(output.rglob("*")) if path.is_file()}
        (output / "validation.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["checks"], indent=2))


if __name__ == "__main__":
    main()
