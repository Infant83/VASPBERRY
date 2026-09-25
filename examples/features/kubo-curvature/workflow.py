"""Logging helpers shared by the real-WAVECAR Kubo and Hall tutorials."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]


def sha256(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def begin(output, feature, wavecar, expected_sha256):
    if not wavecar.is_file() or wavecar.stat().st_size < 1024:
        raise ValueError("A real WAVECAR payload is required; see the tutorial input download instructions.")
    digest = sha256(wavecar)
    if digest != expected_sha256:
        raise ValueError("This reference tutorial needs its documented WAVECAR (SHA256 mismatch). "
                         "Use the README commands with appropriate parameters for another system.")
    output.mkdir(parents=True, exist_ok=False)
    sources = [feature / "run.py", Path(__file__).resolve(), ROOT / "vaspberry.f",
               ROOT / "tools/wavecar_fukui.py", ROOT / "VERSION"]
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                            capture_output=True, text=True, check=False)
    record = {"feature_id": feature.name, "version": (ROOT / "VERSION").read_text().strip(),
              "source_commit": commit.stdout.strip() if commit.returncode == 0 else None,
              "python_version": sys.version.split()[0],
              "input": {"name": "WAVECAR", "sha256": digest, "bytes": wavecar.stat().st_size},
              "source_sha256": {p.relative_to(ROOT).as_posix(): sha256(p) for p in sources},
              "commands": [], "status": "RUNNING"}
    write_json(output / "provenance.json", record)
    return record


def run_command(argv, cwd, output, record, wavecar):
    stage = f"{len(record['commands']) + 1:02d}"
    started = time.monotonic()
    proc = subprocess.run(list(map(str, argv)), cwd=cwd, capture_output=True, text=True, check=False)
    (output / f"{stage}.stdout.log").write_text(proc.stdout)
    (output / f"{stage}.stderr.log").write_text(proc.stderr)
    portable = [str(arg).replace(str(wavecar), "<WAVECAR>")
                .replace(str(output), "<output-dir>")
                .replace(str(ROOT) + "/", "") for arg in argv]
    if portable[0] == sys.executable:
        portable[0] = "python3"
    record["commands"].append({"argv": portable,
        "cwd": str(cwd).replace(str(output), "<output-dir>").replace(str(ROOT), "<repository>"),
        "exit_code": proc.returncode, "elapsed_s": time.monotonic() - started,
        "stdout": f"{stage}.stdout.log", "stderr": f"{stage}.stderr.log"})
    write_json(output / "provenance.json", record)
    if proc.returncode:
        raise RuntimeError(f"Calculation exited with {proc.returncode}; see {stage}.stderr.log")


def finish(output, record, result):
    result["status"] = "PASS" if result["numerical_checks"]["passed"] else "FAIL"
    record["status"] = result["status"]
    result["schema_version"] = 1
    result["workflow_mode"] = "vasp_wavecar_calculation"
    result["provenance"] = record.copy()
    result["output_sha256"] = {p.relative_to(output).as_posix(): sha256(p)
        for p in sorted(output.rglob("*")) if p.is_file()
        and p.name not in ("result.json", "provenance.json")}
    write_json(output / "result.json", result)
    record["output_sha256"] = result["output_sha256"]
    write_json(output / "provenance.json", record)
    if result["status"] != "PASS":
        raise RuntimeError("Reference checks failed; see result.json")


def fail(output, record, error):
    record["status"] = "FAIL"
    record["error"] = str(error)
    write_json(output / "provenance.json", record)
