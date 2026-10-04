#!/usr/bin/env python3
"""Compare real MoS2 native serial/MPI Kubo exports; keep reproducible receipts.

The pair-to-bundle comparison checks the normalization and external-transition
sum independently of the native bundle loop. It is not a PAW accuracy test.
"""
from __future__ import annotations

import argparse
import cmath
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from native_kubo_csv import (BAND_SCHEMA, TRACE_SCHEMA, metadata, read_curvature_csv,
                             require_terminal_pass, positive_integer)
INPUT_SHA256 = "33f8546512856d6c04ad0a80454b18ac9b60e2af4b98f2f85b376ec49b1a8d9f"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_export(path, kind, version):
    all_lines = Path(path).read_text().splitlines()
    meta = metadata(all_lines)
    require_terminal_pass(all_lines)
    schema = {"BUNDLE": TRACE_SCHEMA, "FULLTRACE": TRACE_SCHEMA, "BAND": BAND_SCHEMA,
              "PAIRS": "VASPBERRY_BARE_MOMENTUM_KUBO_PAIRS_V1"}[kind]
    expected = {
        "schema": schema,
        "normalization": "STANDARD_MINUS_TWO_IM",
        "result_status": "PASS",
        "vaspberry_version": version,
    }
    for key, value in expected.items():
        if meta.get(key) != value:
            raise ValueError(f"{path}: expected {key}={value}")
    if kind in ("BUNDLE", "FULLTRACE", "BAND"):
        return read_curvature_csv(path, kind="band" if kind == "BAND" else "trace")[1]
    rows = [{key: float(value) for key, value in row.items()}
            for row in csv.DictReader(line for line in all_lines if line.strip() and not line.startswith("#"))]
    if not rows or not all(math.isfinite(value) for row in rows for value in row.values()):
        raise ValueError(f"{path}: empty or non-finite export")
    nk, ns, nb, count = [positive_integer(meta.get(key), key) for key in
                         ("source_nkpoints", "source_nspin", "source_nbands", "expected_rows")]
    if ns not in (1, 2) or nb < 2 or count != nk*ns*nb*(nb-1)//2 or len(rows) != count:
        raise ValueError(f"{path}: incomplete native pair row count")
    seen = set()
    for row in rows:
        key = tuple(row[name] for name in ("spin", "k_index", "n_band", "m_band"))
        if (any(value != int(value) for value in key) or not 1 <= key[0] <= ns
                or not 1 <= key[1] <= nk or not 1 <= key[2] < key[3] <= nb or key in seen):
            raise ValueError(f"{path}: invalid or duplicate native pair coverage")
        seen.add(key)
    return rows


def compare_rows(left, right, *, label, atol=2e-10, rtol=2e-11):
    if len(left) != len(right) or not left:
        raise ValueError(f"{label}: row count differs or is empty")
    maximum = 0.0
    for row_index, (a, b) in enumerate(zip(left, right), 1):
        if a.keys() != b.keys():
            raise ValueError(f"{label}: columns differ")
        for key in a:
            if a[key] is None and b[key] is None:
                # The full retained-space trace has no external gap (native NA).
                continue
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


def run_recorded(command, directory, *, expected_success=True, diagnostic=None):
    directory.mkdir()
    receipt = {"command": command, "cwd": str(directory), "status": "RUNNING",
               "expected_success": expected_success}
    record_path = directory / "command.json"
    record_path.write_text(json.dumps(receipt, indent=2) + "\n")
    started = time.monotonic()
    try:
        with (directory / "stdout.log").open("wb") as out, (directory / "stderr.log").open("wb") as err:
            result = subprocess.run(command, cwd=directory, stdout=out, stderr=err, timeout=180)
        receipt["exit_code"] = result.returncode
        if (result.returncode == 0) != expected_success:
            raise RuntimeError(f"native command failed ({result.returncode}); see {directory}")
        if not expected_success and result.returncode != 1:
            raise ValueError(f"expected controlled exit 1, got {result.returncode}; see {directory}")
        if diagnostic:
            logs = "\n".join((directory / name).read_text() for name in ("stdout.log", "stderr.log"))
            # List-directed Fortran output may wrap a diagnostic across lines.
            if " ".join(diagnostic.lower().split()) not in " ".join(logs.lower().split()):
                raise ValueError(f"expected diagnostic {diagnostic!r}; see {directory}")
        if not expected_success:
            for path in directory.glob("*.csv*"):
                if path.read_text().strip().endswith("# result_status=PASS"):
                    raise ValueError(f"failed native command left a completed CSV: {path}")
        receipt["status"] = "PASS" if expected_success else "EXPECTED_REJECTION"
    except Exception as exc:
        receipt.update(status="FAIL", error=str(exc))
        raise
    finally:
        receipt["elapsed_seconds"] = time.monotonic() - started
        record_path.write_text(json.dumps(receipt, indent=2) + "\n")


def write_synthetic_wavecar(path, *, word4=False, corruption=None):
    """Tiny independent scalar two-spin fixture; no material/research input."""
    stride, nk, nb, ns = 1024, 2, 3, 2
    data = bytearray(stride * (2 + ns * nk * (nb + 1)))
    struct.pack_into("<3d", data, 0, stride // 4 if word4 else stride, ns, 45200)
    struct.pack_into("<12d", data, stride, nk, nb, 120.,
                     1., 0., 0., 0., 1., 0., 0., 0., 1.)
    for spin in range(ns):
        for k, point in enumerate(((.25, .25, 0.), (-.25, -.25, 0.))):
            record = 2 + (spin * nk + k) * (nb + 1)
            values = [3., *point]
            for band, energy in enumerate((-2., -1., 2.)):
                values.extend((energy + .1*spin, 0., 1. if band < 2 else 0.))
            struct.pack_into("<13d", data, record * stride, *values)
            for band in range(nb):
                for pw in range(3):
                    coefficient = cmath.exp(2j * math.pi * band * pw / 3) / math.sqrt(3)
                    struct.pack_into("<2f", data, (record + 1 + band) * stride + pw * 8,
                                     coefficient.real, coefficient.imag)
    if corruption == "truncated":
        del data[-stride:]
    elif corruption in ("nan", "inf"):
        struct.pack_into("<f", data, len(data)-stride, float(corruption))
    elif corruption == "dimensions":
        struct.pack_into("<d", data, stride, float("nan"))
    elif corruption is not None:
        raise ValueError("unknown synthetic corruption")
    path.write_bytes(data)


def contract_checks(prefix, mode, output, public_wavecar, fixtures, version):
    """Exercise the release interface on each compiler/MPI build with bounded inputs."""
    checks, tables = {}, {}
    for label, flags, diagnostic in (
            ("per-band-gap", ["--task", "kubo", "--kubo-source", "wavecar", "--bands", "17:18", "--per-band", "1"], "not isolated"),
            ("single-gap", ["--task", "kubo", "--kubo-source", "wavecar", "--bands", "18"], "not isolated"),
            ("removed-modern", ["--bundle", "1"], "removed"),
            ("removed-legacy", ["-kubo_bundle", "1"], "removed")):
        directory = output / f"{mode}-{label}"
        run_recorded(prefix + ["--wavecar", str(public_wavecar), *flags], directory,
                     expected_success=False, diagnostic=diagnostic)
        if list(directory.glob("*.csv")) or list(directory.glob("*.dat")):
            raise ValueError(f"selection rejection must precede output: {directory}")
        checks[label] = "EXPECTED_REJECTION_BEFORE_OUTPUT"
    for layout in ("byte", "word4"):
        tables[layout] = {}
        for label, kind, flags, filename in (
                ("trace", "BUNDLE", ["--task", "kubo", "--kubo-source", "wavecar", "--bands", "1:2"], "KUBO.csv"),
                ("full-trace", "FULLTRACE", ["--task", "kubo", "--kubo-source", "wavecar", "--bands", "1:3"], "KUBO.csv"),
                ("per-band", "BAND", ["--task", "kubo", "--kubo-source", "wavecar", "--bands", "1:2", "--per-band", "1",
                                       "--curvature-csv", "BAND.csv"], "BAND.csv"),
                ("pairs", "PAIRS", ["--task", "kubo-pairs", "--kubo-source", "wavecar", "--pairs-csv", "PAIRS.csv"], "PAIRS.csv")):
            directory = output / f"{mode}-{layout}-{label}"
            run_recorded(prefix + ["--wavecar", str(fixtures[layout]), "--spinor", "1", *flags], directory)
            tables[layout][kind] = read_export(directory / filename, kind, version)
    # Directory INQUIRE semantics differ between GNU and Intel. Exercise the
    # real source-directory route on every serial/MPI compiler, including spaces,
    # and ensure explicit file overrides never waive a bad directory.
    source_directory = output / f"{mode} input directory"
    source_directory.mkdir()
    shutil.copyfile(fixtures["byte"], source_directory / "WAVECAR")
    directory = output / f"{mode}-input-directory"
    run_recorded(prefix + ["--input-dir", str(source_directory), "--spinor", "1",
                          "--task", "kubo", "--kubo-source", "wavecar", "--bands", "1:2"], directory)
    selected = read_export(directory / "KUBO.csv", "BUNDLE", version)
    checks["input_directory_with_spaces"] = compare_rows(
        selected, tables["byte"]["BUNDLE"], label="input directory versus explicit WAVECAR")
    for label, bad_directory in (("missing", output / f"{mode}-missing-input"),
                                  ("regular-file", fixtures["byte"])):
        directory = output / f"{mode}-input-directory-{label}"
        run_recorded(prefix + ["--input-dir", str(bad_directory), "--wavecar", str(fixtures["byte"]),
                              "--spinor", "1", "--task", "kubo", "--kubo-source", "wavecar", "--bands", "1:2"],
                     directory, expected_success=False, diagnostic="must be an existing directory")
        if list(directory.glob("*.csv")):
            raise ValueError(f"directory rejection must precede output: {directory}")
        checks[f"input_directory_{label}"] = "EXPECTED_REJECTION_BEFORE_OUTPUT"
    for kind in ("BUNDLE", "FULLTRACE", "BAND", "PAIRS"):
        checks[f"byte_vs_word4_{kind.lower()}"] = compare_rows(tables["byte"][kind], tables["word4"][kind], label=kind)
    for corruption in ("truncated", "nan", "inf", "dimensions"):
        for label, flags in (("trace", []), ("per-band", ["--per-band", "1", "--curvature-csv", "BAND.csv"])):
            directory = output / f"{mode}-{corruption}-{label}"
            run_recorded(prefix + ["--wavecar", str(fixtures[corruption]), "--spinor", "1",
                                  "--task", "kubo", "--kubo-source", "wavecar", "--bands", "1:2", *flags], directory,
                         expected_success=False)
            checks[f"{corruption}-{label}"] = "EXPECTED_REJECTION_WITHOUT_COMPLETED_CSV"
    return checks


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
        fixture_dir = output / "synthetic-inputs"
        fixture_dir.mkdir()
        fixtures = {name: fixture_dir / name for name in ("byte", "word4", "truncated", "nan", "inf", "dimensions")}
        for name, path in fixtures.items():
            write_synthetic_wavecar(path, word4=name == "word4",
                                    corruption=None if name in ("byte", "word4") else name)
        for mode, prefix in (("serial", [str(serial)]),
                             ("mpi", [launcher, *args.mpi_flag, "-n", "2", str(mpi)])):
            tables[mode] = {}
            for task, kind, extra in (("kubo", "BUNDLE", ["--bands", "1:18"]),
                                      ("kubo-pairs", "PAIRS", [])):
                directory = output / f"{mode}-{task}"
                destination = directory / ("KUBO.csv" if task == "kubo" else "PAIRS.csv")
                command = prefix + ["--task", task, "--kubo-source", "wavecar", "--wavecar", str(wavecar), "--spinor", "2",
                                    *extra] + ([] if task == "kubo" else ["--pairs-csv", str(destination)])
                run_recorded(command, directory)
                tables[mode][kind] = read_export(destination, kind, version)
            summary["checks"][f"{mode}_pair_to_bundle"] = pair_bundle_check(
                tables[mode]["BUNDLE"], tables[mode]["PAIRS"])
            summary["checks"][f"{mode}_interface_and_layout"] = contract_checks(
                prefix, mode, output, wavecar, fixtures, version)
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
