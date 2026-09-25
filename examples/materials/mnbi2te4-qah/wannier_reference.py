#!/usr/bin/env python3
"""Prepare, run and collect the full-connection Wannier90 Hall reference."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import datetime
import hashlib
import importlib.util
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import threading
import time

import numpy as np

HERE = Path(__file__).resolve().parent
PACKAGE = HERE / "inputs/wannier/operators"
PRODUCER = "Wannier90 3.1.0 postw90 with effective-model dimension initialization fix"
OPERATOR = "External full J0+J1+J2 Wannier connection; exact H and Hermitian position matrices"
INPUTS = ("wannier90.win", "kpoint.dat", "wannier90_HH_R.dat", "wannier90_AA_R.dat")
OUTPUTS = ("wannier90.wpout", "wannier90-ahc-fermiscan.dat")
COLUMNS = ("mu_eV", "mu_minus_dft_midgap_eV", "sigma_yz_S_per_cm",
           "sigma_zx_S_per_cm", "sigma_xy_S_per_cm", "sigma_xy_e2_over_h")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def fixed_quadrature(lattice, base, refine, radius):
    """Replace selected torus cells by odd, centered 2D submeshes."""
    require(base >= 4 and refine >= 1 and refine % 2 == 1,
            "base must be at least 4 and refine must be a positive odd integer")
    require(np.isfinite(radius) and radius > 0, "radius must be positive and finite")
    reciprocal = 2 * np.pi * np.linalg.inv(lattice).T
    images = np.array([(i, j, 0) for i, j in itertools.product([-1, 0, 1], repeat=2)])
    points, weights, parents, keys = [], [], [], []
    refined = 0
    denominator = base * refine
    for iy in range(base):
        for ix in range(base):
            center = np.array([ix / base, iy / base, 0.])
            distance = np.linalg.norm((center[None, :] + images) @ reciprocal, axis=1).min()
            selected = distance <= radius
            children = range(-(refine // 2), refine // 2 + 1) if selected else [0]
            multiplicity = refine**2 if selected else 1
            first = len(weights)
            for dy, dx in itertools.product(children, repeat=2):
                key = ((ix * refine + dx) % denominator, (iy * refine + dy) % denominator)
                keys.append(key)
                points.append([key[0] / denominator, key[1] / denominator, 0.])
                weights.append(1. / (base**2 * multiplicity))
                parents.append(iy * base + ix)
            require(abs(sum(weights[first:]) - 1 / base**2) < 1e-15, "cell weight mismatch")
            refined += int(selected)
    q, w = np.asarray(points), np.asarray(weights)
    require(len(keys) == len(set(keys)) and np.isfinite(q).all(), "invalid or repeated k points")
    require(np.all(q[:, 2] == 0) and np.all(w > 0) and abs(w.sum() - 1) < 2e-14,
            "the complete two-dimensional weights must sum to one")
    return q, w, np.asarray(parents), reciprocal, refined


def prepare(args):
    output = args.output_dir.resolve()
    require(not output.exists(), "choose a new output directory")
    require(1 <= args.chunks <= 8, "chunks must be between 1 and 8")
    package = json.loads((PACKAGE / "export.json").read_text())
    reference_path = HERE / "reference/direct-dft/metadata.json"
    reference = json.loads(reference_path.read_text())
    reader_patch = Path("inputs/wannier/toolchain/effective-reader-modern31.patch")
    reader_patch_sha256 = sha(HERE / reader_patch)
    lattice = np.asarray(package["lattice_A"], dtype=float)
    require(lattice.shape == (3, 3) and np.isfinite(lattice).all(), "invalid source lattice")
    require(np.array_equal(lattice, reference["lattice_A"]), "operator and energy-reference lattices differ")
    require(package["status"] == "VALIDATED_FULL_OPERATOR_INPUT"
            and package["source_use_ws_distance"] is True
            and package["source_transl_inv"] is False, "unsupported operator package")
    vbm, cbm = reference["vbm_eV"], reference["cbm_eV"]
    require(np.isfinite([vbm, cbm]).all() and cbm > vbm, "reference must have a positive sampled gap")
    midgap = (vbm + cbm) / 2
    mu = np.linspace(vbm + .05 * (cbm - vbm), cbm - .05 * (cbm - vbm), 3)
    q, w, parents, reciprocal, refined = fixed_quadrature(lattice, args.base, args.refine, args.radius)
    require(args.chunks <= len(q), "more partitions than k points")
    expected = {f"wannier90_{key}_R.dat": package[key]["sha256"] for key in ("HH", "AA")}
    if args.operators_dir is not None:
        for name, digest in expected.items():
            require(sha(args.operators_dir / name) == digest, f"operator mismatch: {name}")
    # A retained staging directory distinguishes failed preparation from a
    # complete, runnable directory. Nothing overwrites an existing calculation.
    staging = output.with_name(output.name + ".preparing")
    require(not staging.exists(), f"unfinished preparation exists: {staging}")
    staging.mkdir(parents=True)
    save(staging / "preparation.json", dict(status="PREPARING", started_utc=now()))
    try:
        operators = staging / "operators"
        if args.operators_dir is None:
            spec = importlib.util.spec_from_file_location("restore_wannier_operators", PACKAGE / "restore.py")
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            module.restore(PACKAGE, operators)
        else:
            operators.mkdir()
            for name in expected:
                shutil.copyfile(args.operators_dir / name, operators / name)
        for name, digest in expected.items():
            require(sha(operators / name) == digest, f"restored operator mismatch: {name}")
        cell = "\n".join(" ".join(f"{x:.15g}" for x in row) for row in lattice)
        win = (f"num_wann = {package['num_wann']}\neffective_model = true\nspinors = true\n"
               "use_ws_distance = false\ntransl_inv = false\niprint = 2\n"
               f"begin unit_cell_cart\nAng\n{cell}\nend unit_cell_cart\n"
               "berry = true\nberry_task = ahc\nberry_curv_unit = ang2\n"
               f"berry_kmesh = {args.base} {args.base} 1\nberry_curv_adpt_kmesh = 1\n"
               "berry_curv_adpt_kmesh_thresh = 1e100\nwanint_kpoint_file = true\n"
               f"fermi_energy_min = {mu[0]:.16f}\nfermi_energy_max = {mu[-1]:.16f}\n"
               f"fermi_energy_step = {(mu[-1]-mu[0])/2:.16f}\n")
        np.savez_compressed(staging / "quadrature.npz", kpoints_fractional=q, weights=w,
                            parent_cell=parents, lattice_A=lattice, reciprocal_Ainv=reciprocal)
        parts = []
        for index in range(args.chunks):
            part = staging / f"part{index:02d}"
            part.mkdir()
            for name in expected:
                try:
                    os.link(operators / name, part / name)
                except OSError:
                    shutil.copyfile(operators / name, part / name)
            (part / "wannier90.win").write_text(win)
            indices = np.arange(index, len(q), args.chunks)
            with (part / "kpoint.dat").open("w") as handle:
                handle.write(str(len(indices)) + "\n")
                for i in indices:
                    handle.write(" ".join(f"{x:.17g}" for x in (*q[i], w[i])) + "\n")
            hashes = dict(expected, **{name: sha(part / name) for name in ("wannier90.win", "kpoint.dat")})
            parts.append(dict(directory=part.name, point_count=len(indices),
                              weight_sum=float(w[indices].sum()), input_sha256=hashes))
        manifest = dict(schema="VASPBERRY_EXTERNAL_WANNIER_REFERENCE_V1", status="PREPARED",
                        producer=PRODUCER, operator=OPERATOR, base_mesh=[args.base, args.base, 1],
                        refined_submesh=[args.refine, args.refine, 1],
                        fixed_refinement_radius_Ainv=args.radius, refined_parent_count=refined,
                        local_equivalent_uniform_mesh=[args.base*args.refine]*2+[1],
                        point_count=len(q), weight_sum=float(w.sum()), fixed_kz=0,
                        requested_mu_eV=mu.tolist(), temperature_K=0,
                        mu_scope="three actual samples spanning the central 90% of the sampled DFT gap",
                        reference_midgap_eV=midgap, sampled_dft_vbm_eV=vbm, sampled_dft_cbm_eV=cbm,
                        energy_reference=reference["energy_reference"], lattice_A=lattice.tolist(),
                        original_use_ws_distance=True, effective_use_ws_distance=False,
                        partial_sums_are_not_individually_normalized=True,
                        all_chemical_potentials_share_identical_quadrature=True,
                        operator_package_sha256=sha(PACKAGE / "export.json"),
                        reader_patch=dict(source_relative_path=str(reader_patch), sha256=reader_patch_sha256),
                        energy_reference_metadata_sha256=sha(reference_path),
                        source_operator_sha256=expected, quadrature_sha256=sha(staging / "quadrature.npz"),
                        parts=parts, interpretation="Full finite-model response. Assess mesh and model convergence separately.")
        save(staging / "quadrature.json", manifest)
        save(staging / "preparation.json", dict(status="PREPARED", finished_utc=now()))
        require(not output.exists(), "output appeared during preparation; staging retained")
        staging.rename(output)
    except Exception as error:
        save(staging / "preparation.json", dict(status="FAILED_PREPARATION", error=str(error), finished_utc=now()))
        raise
    return dict(status="PREPARED", output_dir=str(output), points=len(q), partitions=args.chunks,
                chemical_potentials_eV=mu.tolist())


def read_plan(output):
    plan = json.loads((output / "quadrature.json").read_text())
    require(plan.get("schema") == "VASPBERRY_EXTERNAL_WANNIER_REFERENCE_V1", "unsupported preparation schema")
    require(1 <= len(plan["parts"]) <= 8, "invalid partition count")
    for index, part in enumerate(plan["parts"]):
        require(part["directory"] == f"part{index:02d}", "invalid partition order or path")
    require(sha(output / "quadrature.npz") == plan["quadrature_sha256"], "quadrature changed after preparation")
    return plan


def check_inputs(directory, expected):
    require(set(expected) == set(INPUTS), "incomplete input record")
    for name, digest in expected.items():
        require(sha(directory / name) == digest, f"input changed: {directory.name}/{name}")


def validate_response(directory, requested_mu):
    text = (directory / "wannier90.wpout").read_text(errors="replace")
    require("All done: postw90 exiting" in text and not re.search(r"(?im)^\s*(Error:|Exiting\.\.\.)", text),
            f"postw90 did not complete normally: {directory.name}")
    require(re.search(r"Release:\s*3\.1\.0\b", text), "expected the documented Wannier90 3.1.0 reader")
    require(not (directory / "wannier90.werr").exists(), "postw90 left an error file")
    response = np.loadtxt(directory / "wannier90-ahc-fermiscan.dat", ndmin=2)
    requested_mu = np.asarray(requested_mu)
    require(response.shape == (len(requested_mu), 4) and np.isfinite(response).all(), "invalid AHC response table")
    require(np.allclose(response[:, 0], requested_mu, rtol=0, atol=5.1e-7), "response chemical potentials differ from prepared values")
    terms = []
    for label in ("J0", "J1", "J2"):
        rows = re.findall(r"(?m)^\s*" + label + r" term\s*:\s*([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)", text)
        values = np.asarray(rows, dtype=float)
        require(values.shape == response[:, 1:].shape and np.isfinite(values).all(), f"missing full {label} contribution")
        terms.append(values)
    terms = np.asarray(terms)
    require(np.max(np.abs(terms.sum(axis=0) - response[:, 1:])) < 1.6e-4,
            "J0+J1+J2 does not match the full response within printed precision")
    return response, terms


def run_part(output, part, binary, binary_hash, timeout, requested_mu, stop):
    directory = output / part["directory"]
    record = dict(status="RUNNING", started_utc=now(), command=[str(binary), "wannier90"],
                  binary_sha256=binary_hash, input_sha256=part["input_sha256"],
                  timeout_seconds=timeout, threads_per_process=1)
    record_path = directory / "run.json"
    save(record_path, record)
    start = time.monotonic()
    proc = None
    try:
        require(not stop.is_set(), "another partition failed or the run was interrupted")
        check_inputs(directory, part["input_sha256"])
        env = dict(os.environ)
        for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
            env[key] = "1"
        with (directory / "stdout.log").open("x") as stdout, (directory / "stderr.log").open("x") as stderr:
            proc = subprocess.Popen(record["command"], cwd=directory, env=env,
                                    stdout=stdout, stderr=stderr, start_new_session=True)
            record["pid"] = proc.pid
            save(record_path, record)
            while proc.poll() is None:
                if stop.is_set() or time.monotonic() - start > timeout:
                    raise RuntimeError("execution interrupted" if stop.is_set() else "partition exceeded its time limit")
                time.sleep(.2)
        record["returncode"] = proc.returncode
        require(proc.returncode == 0, f"postw90 returned {proc.returncode}")
        validate_response(directory, requested_mu)
        check_inputs(directory, part["input_sha256"])
        require(sha(binary) == binary_hash, "executable changed during calculation")
        record["status"] = "FINISHED"
    except Exception as error:
        if proc is not None and proc.poll() is None:
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
        record.update(status="FAILED", error=str(error), returncode=None if proc is None else proc.returncode)
        stop.set()
    record.update(finished_utc=now(), wall_seconds=time.monotonic()-start,
                  output_sha256={name:sha(directory/name) for name in (*OUTPUTS, "stdout.log", "stderr.log", "wannier90.werr")
                                 if (directory/name).is_file()})
    save(record_path, record)
    return record


def run(args):
    output = args.output_dir.resolve()
    require(1 <= args.workers <= 8 and np.isfinite(args.timeout_seconds) and args.timeout_seconds > 0,
            "workers must be 1–8 and timeout must be positive")
    binary = args.postw90.expanduser().resolve()
    require(binary.is_file() and os.access(binary, os.X_OK), "postw90 must be an executable file")
    plan = read_plan(output)
    require(not (output / "execution.json").exists(), "this prepared calculation has already been run; prepare a new directory")
    for part in plan["parts"]:
        directory = output / part["directory"]
        check_inputs(directory, part["input_sha256"])
        require(set(path.name for path in directory.iterdir()) == set(INPUTS),
                f"stale or unexpected files in {directory.name}; prepare a new directory")
    binary_hash = sha(binary)
    execution = dict(status="RUNNING", producer=PRODUCER, started_utc=now(),
                     binary_sha256=binary_hash, workers=min(args.workers, len(plan["parts"])),
                     threads_per_process=1)
    save(output / "execution.json", execution)
    stop = threading.Event()
    results = []
    failure = None
    with ThreadPoolExecutor(max_workers=execution["workers"]) as pool:
        futures = [pool.submit(run_part, output, part, binary, binary_hash, args.timeout_seconds,
                               plan["requested_mu_eV"], stop) for part in plan["parts"]]
        try:
            for future in as_completed(futures):
                results.append(future.result())
        except BaseException as error:
            stop.set()
            failure = f"{type(error).__name__}: {error}"
            results = []
            for future in futures:
                try:
                    results.append(future.result())
                except Exception as child_error:
                    results.append(dict(status="FAILED", error=str(child_error)))
    execution.update(status="FINISHED" if failure is None and len(results) == len(plan["parts"]) and all(r["status"] == "FINISHED" for r in results) else "FAILED",
                     finished_utc=now(), completed_partitions=sum(r["status"] == "FINISHED" for r in results))
    if failure is not None:
        execution["error"] = failure
    save(output / "execution.json", execution)
    require(execution["status"] == "FINISHED", "calculation failed; inspect retained execution and partition logs")
    return execution


def collect(args):
    output = args.output_dir.resolve()
    destination = args.result_dir.resolve() if args.result_dir else output / "results"
    require(not destination.exists(), "choose a new result directory")
    plan = read_plan(output)
    execution = json.loads((output / "execution.json").read_text())
    require(execution["status"] == "FINISHED", "execution did not finish successfully")
    with np.load(output / "quadrature.npz", allow_pickle=False) as grid:
        q, weights, lattice = grid["kpoints_fractional"], grid["weights"], grid["lattice_A"]
    expected_q, expected_w, _, _, _ = fixed_quadrature(lattice, plan["base_mesh"][0],
                                                       plan["refined_submesh"][0], plan["fixed_refinement_radius_Ainv"])
    require(np.array_equal(q, expected_q) and np.array_equal(weights, expected_w), "prepared geometry or weights are inconsistent")
    total = np.zeros((len(plan["requested_mu_eV"]), 3))
    terms = np.zeros((3, *total.shape))
    printed_mu = None
    timings = []
    for index, part in enumerate(plan["parts"]):
        directory = output / part["directory"]
        record = json.loads((directory / "run.json").read_text())
        require(record["status"] == "FINISHED" and record["returncode"] == 0, "a partition failed")
        require(record["binary_sha256"] == execution["binary_sha256"], "partitions used different executables")
        require(record["input_sha256"] == part["input_sha256"], "partition input records differ")
        check_inputs(directory, part["input_sha256"])
        for name in OUTPUTS:
            require(sha(directory/name) == record["output_sha256"][name], f"output changed after execution: {name}")
        data = np.loadtxt(directory / "kpoint.dat", skiprows=1, ndmin=2)
        indices = np.arange(index, len(q), len(plan["parts"]))
        require(data.shape == (len(indices), 4) and np.array_equal(data[:, :3], q[indices])
                and np.array_equal(data[:, 3], weights[indices]), "partition does not match the complete fixed quadrature")
        require(int((directory/"kpoint.dat").read_text().splitlines()[0]) == len(indices), "wrong partition point count")
        response, decomposition = validate_response(directory, plan["requested_mu_eV"])
        if printed_mu is None:
            printed_mu = response[:, 0]
        require(np.array_equal(printed_mu, response[:, 0]), "partitions used different chemical potentials")
        total += response[:, 1:]
        terms += decomposition
        timings.append(dict(part=index, wall_seconds=record["wall_seconds"]))
    height = abs(np.linalg.det(lattice)) / np.linalg.norm(np.cross(lattice[0], lattice[1]))
    g0 = (1.602176634e-19)**2 / 6.62607015e-34
    factor = height * 1e-8 / g0
    rows = np.column_stack([printed_mu, printed_mu-plan["reference_midgap_eV"], total, total[:, 2]*factor])
    metadata = dict(schema="VASPBERRY_EXTERNAL_WANNIER_AHC_V1", status="INTEGRATED", producer=PRODUCER,
                    operator=OPERATOR, temperature_K=0, quadrature=plan,
                    source_binary_sha256=execution["binary_sha256"],
                    actual_output_mu_eV=printed_mu.tolist(), requested_mu_eV=plan["requested_mu_eV"],
                    all_components_S_per_cm=total.tolist(), sigma_xy_e2_over_h=(total[:, 2]*factor).tolist(),
                    decomposition_J0_J1_J2_S_per_cm=terms.tolist(),
                    decomposition_print_precision_S_per_cm=1e-4,
                    simulation_cell_normal_height_A=height, conductance_quantum_S=g0,
                    sheet_conversion_factor=factor, conversion_constants="exact modern SI e and h; raw producer S/cm retained",
                    chern_hall_convention="sigma_xy/(e^2/h)=-C", timings=timings,
                    no_rounding_or_symmetry_average_applied=True,
                    interpretation="External full-connection finite-model reference. Neither gap flatness nor representation parity proves integration or material convergence.")
    destination.mkdir(parents=True, exist_ok=False)
    with (destination / "conductivity.csv").open("x", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(COLUMNS)
        writer.writerows(rows)
    np.savetxt(destination / "conductivity.dat", rows,
               header="External full J0+J1+J2 postw90 reference\n" + " ".join(COLUMNS))
    np.savez_compressed(destination / "conductivity.npz", **{name:rows[:, i] for i, name in enumerate(COLUMNS)},
                        producer=np.asarray(PRODUCER), operator=np.asarray(OPERATOR),
                        metadata_json=np.asarray(json.dumps(metadata, allow_nan=False)))
    save(destination / "integration.json", metadata)
    return dict(status="INTEGRATED", result_dir=str(destination), sigma_xy_e2_over_h=metadata["sigma_xy_e2_over_h"])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("prepare", help="restore full operators and prepare fixed 2D partitions")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--operators-dir", type=Path, help="already restored exact HH_R and AA_R directory")
    p.add_argument("--base", type=int, default=80)
    p.add_argument("--refine", type=int, default=9)
    p.add_argument("--radius", type=float, default=.18, help="periodic Cartesian refinement radius in inverse angstrom")
    p.add_argument("--chunks", type=int, default=8)
    p.set_defaults(function=prepare)
    p = commands.add_parser("run", help="run serial postw90 partitions with one thread per process")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--postw90", type=Path, required=True, help="3.1.0 executable with the documented effective-reader fix")
    p.add_argument("--workers", type=int, default=1, help="concurrent serial processes, 1–8 (default: 1)")
    p.add_argument("--timeout-seconds", type=float, default=3600, help="maximum wall time per partition (default: 3600)")
    p.set_defaults(function=run)
    p = commands.add_parser("collect", help="validate all partitions and write external CSV, DAT and NPZ results")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--result-dir", type=Path, help="new destination; default: OUTPUT/results")
    p.set_defaults(function=collect)
    args = parser.parse_args(argv)
    try:
        result = args.function(args)
    except (ValueError, OSError, KeyError, json.JSONDecodeError) as error:
        parser.exit(1, f"{error}\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
