#!/usr/bin/env python3
"""Validate supplied serial/MPI executables on synthetic native spin fixtures.

Needs NumPy, never invokes a compiler, and never substitutes another executable.
The QWZ fixture has known opposite spin-Chern integers. The finite-band spin-Kubo
fixture is checked using tiny overlap loops of displaced raw wavefunctions,
independently of the Fortran sector-derivative formula. Neither is a material
calculation or a validation of physical PAW spin-Hall response.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
# Importing these helpers does not run their unittest classes or compile code.
sys.path.insert(0, str(ROOT / "tests"))
from mpi_validation_launcher import mpi_validation_command
from test_fortran_spin_chern import fixture, table
from test_fortran_spin_kubo import response_fixture

COMPONENTS = ("yz", "zx", "xy")
SECTORS = ("plus", "minus", "parent")


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require_current_producer(metadata):
    version = (ROOT / "VERSION").read_text().strip()
    if metadata.get("producer") != f"VASPBERRY {version}":
        raise AssertionError(f"binary output producer does not match VERSION {version}")


def executable(value):
    candidate = Path(value).expanduser()
    resolved = candidate.resolve() if candidate.is_file() else shutil.which(value)
    if not resolved or not os.access(resolved, os.X_OK):
        raise ValueError(f"executable not found or not executable: {value}")
    return str(Path(resolved).resolve())


def invoke(command, case, expected_error=None, timeout=120):
    """Always retain an invocation receipt, including timeout and launch failure."""
    receipt = {"command": command, "cwd": str(case), "expected_error": expected_error,
               "started_utc": datetime.now(timezone.utc).isoformat(), "status": "RUNNING"}
    dump(case / "command.json", receipt)
    started = time.monotonic()
    try:
        with (case / "stdout.log").open("w") as out, (case / "stderr.log").open("w") as err:
            result = subprocess.run(command, cwd=case, stdout=out, stderr=err,
                                    check=False, timeout=timeout)
        receipt["exit_code"] = result.returncode
        stderr = (case / "stderr.log").read_text(errors="replace")
        if expected_error:
            if result.returncode != 1 or expected_error not in stderr:
                raise AssertionError(f"{case.name}: expected rejection {expected_error!r}")
            if list(case.glob("SPIN*.csv")):
                raise AssertionError(f"{case.name}: rejected input produced result files")
        elif result.returncode != 0:
            raise AssertionError(f"{case.name}: exit {result.returncode}; see stdout/stderr.log")
        receipt["status"] = "EXPECTED_REJECTION" if expected_error else "PASS"
    except Exception as error:
        receipt.update(status="FAIL", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        receipt["elapsed_seconds"] = time.monotonic() - started
        dump(case / "command.json", receipt)


def response_case(path, *, split=True, zero_omitted=False):
    response_fixture(path)
    data = bytearray((path / "WAVECAR").read_bytes())
    recl = int(struct.unpack_from("d", data)[0])
    for k in range(16):
        record = 2 + 5*k
        if split:
            struct.pack_into("d", data, record*recl + 8*(4+3*3), 2.)
        if zero_omitted:
            nc = int(struct.unpack_from("d", data, record*recl)[0])
            start = (record+4)*recl
            data[start:start+8*nc] = bytes(8*nc)
    (path / "WAVECAR").write_bytes(data)


def sector_frames(c):
    """Orthonormalize by QR and diagonalize spin, without tangent formulae."""
    q, _ = np.linalg.qr(c)
    ng = len(q)//2
    spin_q = np.concatenate((q[:ng], -q[ng:]))
    values, vectors = np.linalg.eigh(q.conj().T @ spin_q)
    if abs(values).min() < .1:
        raise AssertionError("reference projected-spin gap unexpectedly small")
    return (q @ vectors[:, values > 0], q @ vectors[:, values < 0], q)


def loop_curvature(c, dc, step):
    """Berry curvature from overlap phases, with no eigenvector derivatives."""
    result = np.empty((3, 3))
    for j, (a, b) in enumerate(((1, 2), (2, 0), (0, 1))):
        corners = [sector_frames(c + step*(x*dc[a] + y*dc[b]))
                   for x, y in ((-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5))]
        for sector in range(3):
            product = 1.+0j
            for t in range(4):
                overlap = corners[t][sector].conj().T @ corners[(t+1) % 4][sector]
                product *= np.linalg.slogdet(overlap)[0]
            result[sector, j] = -np.angle(product)/step**2
    return result


def numpy_reference(path, sum_max):
    """Read the exact complex64 fixture, build canonical tangents, then loops."""
    data = path.read_bytes()
    recl = int(struct.unpack_from("d", data)[0])
    nk, nb, cutoff = struct.unpack_from("3d", data, recl)
    if (int(nk), int(nb), cutoff) != (16, 4, 100.):
        raise AssertionError("reference expects the bounded four-band fixture")
    sequence = list(range(9)) + list(range(-8, 0))
    all_g = np.array([(x, y, z) for z in sequence for y in sequence for x in sequence])
    factor = (6.582119569e-16 * 2.99792458e8 * 1e10)**2 / .510998950e6
    reference, refinement = [], []
    for k in range(int(nk)):
        record = 2 + 5*k
        header = np.frombuffer(data, dtype=np.float64, count=16, offset=record*recl)
        count, point = int(header[0]), header[1:4]
        energies = header[4:].reshape(4, 3)[:, 0]
        g = all_g[np.sum((all_g+point)**2, axis=1)/.262465831 < cutoff]
        if 2*len(g) != count:
            raise AssertionError("reference plane-wave count mismatch")
        states = np.array([np.frombuffer(data, dtype=np.complex64, count=count,
                                        offset=(record+n+1)*recl)
                           for n in range(sum_max)], dtype=complex).T
        parent, outside = states[:, :2], states[:, 2:]
        denominator = energies[:2][None, :] - energies[2:sum_max, None]
        dc = np.array([outside @ ((outside.conj().T @
                                  (factor*np.tile((g+point)[:, a], 2)[:, None]*parent))
                                 / denominator) for a in range(3)])
        coarse = loop_curvature(parent, dc, 4e-4)
        fine = loop_curvature(parent, dc, 2e-4)
        reference.append((4*fine-coarse)/3)
        refinement.append(float(abs(fine-coarse).max()))
    return np.array(reference), max(refinement)


def curves(case):
    meta, rows = table(case / "SPIN_KUBO.csv")
    require_current_producer(meta)
    if meta.get("result_status") != "PASS" or len(rows) != 16:
        raise AssertionError(f"{case.name}: incomplete spin-Kubo table")
    omega = np.array([[[float(row[f"omega_{sector}_{component}_A2"])
                        for component in COMPONENTS] for sector in SECTORS] for row in rows])
    base = np.array([[[float(row[f"parent_projected_{sector}_{component}_A2"])
                       for component in COMPONENTS] for sector in SECTORS[:2]] for row in rows])
    mixing = np.array([[[float(row[f"spin_mixing_{sector}_{component}_A2"])
                         for component in COMPONENTS] for sector in SECTORS[:2]] for row in rows])
    np.testing.assert_allclose(omega[:, :2], base+mixing, rtol=2e-11, atol=2e-11)
    np.testing.assert_allclose(omega[:, :2].sum(axis=1), omega[:, 2], rtol=2e-11, atol=2e-11)
    np.testing.assert_allclose(mixing.sum(axis=1), 0, atol=2e-11)
    if not np.isfinite(omega).all() or abs(mixing).max() <= 1e-6:
        raise AssertionError("fixture must expose a finite, nonzero sector derivative")
    integral_meta, integral_rows = table(case / "SPIN_KUBO_INTEGRAL.csv")
    if integral_meta.get("result_status") != "PASS" or len(integral_rows) != 1:
        raise AssertionError("missing completed mesh integral")
    observed = [float(integral_rows[0][f"C_est_{s}"]) for s in SECTORS]
    np.testing.assert_allclose(observed, omega[:, :, 2].sum(axis=0)/(16*2*np.pi),
                               rtol=2e-11, atol=2e-11)
    return meta, omega


def equal_outputs(left, right, names):
    for name in names:
        # Each k-point is owned by one MPI rank; only disjoint zero-padded
        # arrays are reduced. The same build's serial/MPI CSVs must agree.
        if (left / name).read_bytes() != (right / name).read_bytes():
            raise AssertionError(f"{left.name} versus {right.name}: {name} differs")


def validate(args, output):
    serial, mpi, launcher = map(executable, (args.serial, args.mpi, args.mpiexec))
    report = {"scope": "synthetic pseudo-wavefunction kernel validation; not material/PAW validation",
              "version": (ROOT / "VERSION").read_text().strip(),
              "executables": {p: sha256(Path(p)) for p in (serial, mpi, launcher)},
              "mpi_ranks": 2, "numpy_version": np.__version__, "checks": [],
              "mpi_signal_policy": "ignore_sigpipe" if getattr(args, "mpi_ignore_sigpipe", False) else "default"}
    dump(output / "environment.json", {
        "python": sys.version, "platform": sys.platform,
        "runtime": {k: os.environ[k] for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS", "I_MPI_FABRICS", "I_MPI_PIN",
                    "OMPI_ALLOW_RUN_AS_ROOT", "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM") if k in os.environ},
        **report})

    def run_case(name, task, generate, flags=(), expected_error=None):
        paths = []
        for mode, binary in (("serial", serial), ("mpi", mpi)):
            case = output / f"{name}-{mode}"
            generate(case)
            command = [binary, "--task", task, "--bands", "1:2", *flags]
            if task == "spin-kubo":
                command += ["--kubo-source", "wavecar"]
            if mode == "mpi":
                command = mpi_validation_command([launcher, *args.mpi_arg, "-n", "2", *command],
                                                 ignore_sigpipe=getattr(args, "mpi_ignore_sigpipe", False))
            invoke(command, case, expected_error, args.timeout)
            paths.append(case)
        return paths

    chern = run_case("chern", "spin-chern", lambda p: fixture(p, n=8), ["--mesh", "8,8"])
    equal_outputs(*chern, ["SPIN_CHERN.csv", "SPIN_BERRY.csv", "SPIN_SPECTRUM.csv"])
    meta, rows = table(chern[0] / "SPIN_CHERN.csv")
    require_current_producer(meta)
    if meta.get("result_status") != "PASS" or len(rows) != 1:
        raise AssertionError("incomplete spin-Chern result")
    np.testing.assert_allclose([float(rows[0][k]) for k in
                               ("C_plus", "C_minus", "C_spin", "C_charge", "C_parent")],
                               [1, -1, 1, 0, 0], atol=2e-10, rtol=0)
    if len(table(chern[0] / "SPIN_BERRY.csv")[1]) != 64:
        raise AssertionError("spin-Chern plaquette count mismatch")
    if len(table(chern[0] / "SPIN_SPECTRUM.csv")[1]) != 128:
        raise AssertionError("spin-Chern spectrum count mismatch")
    report["checks"].append("known C_plus=1, C_minus=-1, C_spin=1; serial/MPI identical")

    cases, values = {}, {}
    names = ["SPIN_KUBO.csv", "SPIN_KUBO_SPECTRUM.csv", "SPIN_KUBO_INTEGRAL.csv"]
    for name, flags in (("default", []), ("explicit", ["--sum-bands", "4"]),
                        ("truncated", ["--sum-bands", "3"]), ("zero-omitted", [])):
        paths = run_case(name, "spin-kubo",
                         lambda p: response_case(p, zero_omitted=name == "zero-omitted"),
                         ["--mesh", "4,4", *flags])
        equal_outputs(*paths, names)
        meta, values[name] = curves(paths[0])
        sum_max = 3 if name == "truncated" else 4
        if (meta.get("source_nbands"), meta.get("sum_band_max"),
                meta.get("summed_external_bands")) != ("4", str(sum_max), str(sum_max-2)):
            raise AssertionError(f"{name}: incorrect sum-band metadata")
        if sum_max == 3:
            np.testing.assert_allclose(float(meta["sum_boundary_gap_eV"]), 1., rtol=0, atol=1e-12)
        if len(table(paths[0] / "SPIN_KUBO_SPECTRUM.csv")[1]) != 32:
            raise AssertionError("spin-Kubo spectrum count mismatch")
        if name in ("default", "truncated"):
            reference, refinement = numpy_reference(paths[0] / "WAVECAR", sum_max)
            dump(paths[0] / "numpy-reference.json", {
                "method": "Richardson-extrapolated 4e-4 and 2e-4 overlap loops",
                "components": COMPONENTS, "sectors": SECTORS, "omega_A2": reference.tolist(),
                "max_step_refinement_A2": refinement,
                "max_absolute_error_A2": float(abs(values[name]-reference).max()),
                "comparison_atol_A2": 2e-7, "comparison_rtol": 2e-5})
            np.testing.assert_allclose(values[name], reference, rtol=2e-5, atol=2e-7)
        cases[name] = paths[0]
    equal_outputs(cases["default"], cases["explicit"], names)
    np.testing.assert_allclose(values["truncated"], values["zero-omitted"], rtol=2e-11, atol=2e-11)
    if abs(values["truncated"]-values["default"]).max() <= 1e-6:
        raise AssertionError("truncating nonzero external couplings must change the result")
    report["checks"].extend([
        "all spin-Kubo components: NumPy overlap-loop reference for full/truncated sums",
        "sector sum and projected-parent/mixing closure; explicit mesh integral",
        "default equals explicit NBANDS; truncation equals zero omitted couplings",
        "full/default/truncated/zero-omitted serial/MPI outputs identical"])
    run_case("split-degeneracy", "spin-kubo", lambda p: response_case(p, split=False),
             ["--sum-bands", "3"], "cuts an unresolved external degeneracy")
    report["checks"].append("serial/MPI reject a sum cutoff through a degenerate multiplet without output")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--serial", required=True, help="prebuilt serial executable")
    parser.add_argument("--mpi", required=True, help="prebuilt MPI executable")
    parser.add_argument("--mpiexec", default="mpiexec", help="MPI launcher executable")
    parser.add_argument("--mpi-arg", action="append", default=[], help="repeatable launcher argument; use --mpi-arg=--flag")
    parser.add_argument("--mpi-ignore-sigpipe", action="store_true",
                        help="Intel Hydra validation: ignore SIGPIPE at launcher exec; keep real exit status")
    parser.add_argument("--output-dir", required=True, type=Path, help="new directory; existing paths are refused")
    parser.add_argument("--timeout", default=120., type=float, help="maximum seconds per native invocation")
    args = parser.parse_args(argv)
    if not np.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("--timeout must be finite and positive")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    report = {"status": "RUNNING", "started_utc": datetime.now(timezone.utc).isoformat()}
    dump(output / "summary.json", report)
    try:
        report.update(validate(args, output), status="PASS")
    except Exception as error:
        report.update(status="FAIL", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        report["finished_utc"] = datetime.now(timezone.utc).isoformat()
        dump(output / "summary.json", report)
        dump(output / "checksums.json", {str(p.relative_to(output)): {"bytes": p.stat().st_size,
            "sha256": sha256(p)} for p in sorted(output.rglob("*")) if p.is_file()})
    print(f"PASS: native spin-Chern and spin-Kubo validation; receipts in {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
