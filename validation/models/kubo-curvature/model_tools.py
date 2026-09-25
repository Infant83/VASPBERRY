"""Small example-only QWZ oracle and CLI/provenance helpers.

The Hall feature imports this sibling helper. No production kernel is imported
for the d-vector oracle or the independent Hall quadrature.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

import matplotlib
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BASE_COMMIT = "49fd65e521b7eca40abc4ecc2abcfff62130fcb9"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def source_hashes(feature_dir):
    paths = [ROOT / "VERSION", *(ROOT / "tools" / name for name in (
        "vaspberry_kubo.py", "berry_data.py", "exported_matrix_kubo.py",
        "vaspberry_transport.py", "wavecar_fukui.py"))]
    paths += list(feature_dir.glob("*.py")) + list(feature_dir.glob("*.json"))
    paths.append(Path(__file__).resolve())
    return {p.relative_to(ROOT).as_posix(): sha256(p) for p in sorted(set(paths))}


def begin(output, feature_dir):
    output.mkdir(parents=True, exist_ok=False)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                            capture_output=True, text=True, check=False)
    manifest = {
        "feature": feature_dir.name,
        "model_scope": "public analytic two-band model; not a material benchmark",
        "reference_base_commit": BASE_COMMIT,
        "source_commit": commit.stdout.strip() if commit.returncode == 0 else None,
        "source_commit_lookup_exit_code": commit.returncode,
        "version": (ROOT / "VERSION").read_text().strip(),
        "python_version": sys.version.split()[0],
        "numpy_version": np.__version__,
        "matplotlib_version": matplotlib.__version__,
        "source_sha256_before": source_hashes(feature_dir),
        "path_convention": "repository-relative sources; generated files relative to output directory",
        "commands": [],
        "status": "running",
    }
    write_json(output / "run.json", manifest)
    return manifest


def run_cli(arguments, output, manifest):
    command = [sys.executable, str(ROOT / "tools/vaspberry_kubo.py"), *map(str, arguments)]
    stage = f"{len(manifest['commands']) + 1:02d}-{arguments[0]}"
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=False)
    (output / f"{stage}.stdout").write_text(result.stdout)
    (output / f"{stage}.stderr").write_text(result.stderr)
    portable = ["python3", "tools/vaspberry_kubo.py"] + [
        str(arg).replace(str(output), "<output-dir>").replace(str(ROOT) + "/", "")
        for arg in arguments]
    manifest["commands"].append({"argv": portable, "exit_code": result.returncode,
        "stdout": f"{stage}.stdout", "stderr": f"{stage}.stderr"})
    write_json(output / "run.json", manifest)
    if result.returncode:
        raise RuntimeError(f"{stage} failed with exit {result.returncode}; see saved stderr")


def model_curvature(config, output, manifest):
    mesh = config["mesh"]
    run_cli(["demo", "--mesh", mesh, "--mass", config["mass"],
             "--output-dir", output / "model"], output, manifest)
    run_cli(["matrix", "--matrices", output / "model/matrix.npz",
             "--metadata", output / "model/matrix.json", "--n-bands", "1:2",
             "--m-bands", "1:2", "--mesh", mesh, mesh, "--plane-axes", 0, 1,
             "--spin-multiplicity", config["spin_multiplicity"],
             "--energy-reference", config["energy_reference"],
             "--degeneracy-threshold-eV", config["degeneracy_threshold_eV"],
             "--output-dir", output / "curvature"], output, manifest)
    with np.load(output / "curvature/curvature.npz", allow_pickle=False) as source:
        return {name: source[name] for name in source.files}


def oracle(q, mass):
    """Independent closed form for a=1 Angstrom, energy coefficient 1 eV.

    A=+i<u|d u>, ordered (kx,ky): lower Omega=+d.(dx d cross dy d)/(2|d|^3).
    """
    x, y = (2 * np.pi * q[:, :2]).T
    norm = np.sqrt(np.sin(x)**2 + np.sin(y)**2 + (mass + np.cos(x) + np.cos(y))**2)
    lower = (np.cos(x) + np.cos(y) + mass * np.cos(x) * np.cos(y)) / (2 * norm**3)
    return np.column_stack((-norm, norm)), np.column_stack((lower, -lower))


def check_curvature(data, config):
    energies, omega = oracle(data["kpoints_fractional"], config["mass"])
    actual = data["omega_A2"][:, :, 2]
    point_error = float(np.max(abs(actual - omega)))
    energy_error = float(np.max(abs(data["energies_eV"] - energies)))
    pair_error = float(np.max(abs(actual.sum(axis=1))))
    chern = 2 * np.pi * np.einsum("k,kb->b", data["weights"], actual)
    expected = np.array([1., -1.])
    if config["mass"] != -1 or config["spin_multiplicity"] != 1:
        raise ValueError("These fixed reference checks require mass=-1 and spin_multiplicity=1")
    checks = {
        "pointwise_oracle_max_abs_error_A2": point_error,
        "analytic_energy_max_abs_error_eV": energy_error,
        "band_pair_sum_max_abs_A2": pair_error,
        "sampled_band_chern": chern.tolist(),
        "continuum_expected_band_chern": expected.tolist(),
        "sampled_chern_max_abs_error": float(np.max(abs(chern - expected))),
        "minimum_isolation_gap_eV": float(data["min_gap_eV"].min()),
        "all_states_valid_nondegenerate": bool(data["valid_nondegenerate"].all()),
        "pointwise_absolute_tolerance_A2": config["pointwise_absolute_tolerance_A2"],
        "chern_absolute_tolerance": config["chern_absolute_tolerance"],
    }
    checks["passed"] = bool(point_error < config["pointwise_absolute_tolerance_A2"]
        and energy_error < 1e-12 and pair_error < 1e-12
        and checks["sampled_chern_max_abs_error"] < config["chern_absolute_tolerance"]
        and checks["minimum_isolation_gap_eV"] > 0 and checks["all_states_valid_nondegenerate"])
    return checks, energies, omega


def finish(output, feature_dir, manifest, success, error=None):
    manifest["source_sha256_after"] = source_hashes(feature_dir)
    unchanged = manifest["source_sha256_before"] == manifest["source_sha256_after"]
    manifest["sources_unchanged_during_run"] = unchanged
    manifest["status"] = "passed" if success and unchanged else "failed"
    manifest["exit_code"] = 0 if success and unchanged else 1
    if error:
        manifest["error"] = str(error)
    if (output / "result.json").exists():
        result = json.loads((output / "result.json").read_text())
        result["status"] = "PASS" if success and unchanged else "FAIL"
        result["output_sha256"] = {name: sha256(output/name) for name in ("summary.csv", "figure.png")}
        write_json(output / "result.json", result)
    manifest["output_sha256"] = {p.relative_to(output).as_posix(): sha256(p)
        for p in sorted(output.rglob("*")) if p.is_file() and p.name != "run.json"}
    write_json(output / "run.json", manifest)
    return manifest["exit_code"]


def result_provenance(output, manifest):
    return {key: manifest[key] for key in ("reference_base_commit", "source_commit", "version",
        "python_version", "numpy_version", "matplotlib_version", "source_sha256_before", "commands", "path_convention")} | {
        "generated_data_sha256": {p.relative_to(output).as_posix(): sha256(p)
            for p in sorted(output.rglob("*")) if p.is_file() and p.suffix in (".npz", ".json")
            and p.name not in ("run.json", "result.json")}}
