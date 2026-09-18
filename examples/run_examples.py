#!/usr/bin/env python3
"""List or run the public feature examples, retaining failures and their logs."""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import time
from datetime import datetime, timezone


EXAMPLES = Path(__file__).resolve().parent
ROOT = EXAMPLES.parent


def main() -> int:
    catalog = json.loads((EXAMPLES / "catalog.json").read_text())
    features = {item["id"]: item for item in catalog["features"]}
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("features", nargs="*", help="Feature IDs shown by --list")
    parser.add_argument("--list", action="store_true", help="List inputs and default workflow modes")
    parser.add_argument("--all", action="store_true", help="Run every default workflow, including stored Bi validation")
    parser.add_argument("--output-dir", type=Path, help="New directory for all generated files and logs")
    parser.add_argument("--binary", type=Path, default=ROOT / "build/vaspberry-gfortran",
                        help="Current Fortran binary for optical and wavefunction examples")
    args = parser.parse_args()
    if args.all and args.features:
        parser.error("use feature IDs or --all, not both")
    unknown = sorted(set(args.features) - features.keys())
    if unknown:
        parser.error("unknown feature IDs: " + ", ".join(unknown))
    if args.list or not (args.all or args.features):
        for item in features.values():
            print(f"{item['id']:20} {item['default_mode']:28} {item['title']}")
        print("\nSee examples/README.md for inputs, reference results and full recalculation.")
        return 0
    selected = list(features) if args.all else list(dict.fromkeys(args.features))
    if args.output_dir is None:
        parser.error("--output-dir is required when running examples")
    if args.output_dir.exists():
        parser.error("output directory already exists; choose a new path")
    for module in ("numpy", "matplotlib"):
        if importlib.util.find_spec(module) is None:
            parser.error("install requirements-transport.txt in this Python environment")
    binary = args.binary.resolve()
    if any(features[key]["requires_fortran"] for key in selected) and not binary.is_file():
        parser.error(f"missing Fortran binary: {binary}; run make serial or pass --binary")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    report = {
        "schema_version": 1,
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "features": [],
        "status": "RUNNING",
    }
    report_path = output / "run_manifest.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    failed = False
    for key in selected:
        item = features[key]
        command = [sys.executable, str(EXAMPLES / item["runner"]),
                   "--output-dir", str(output / key)]
        if item["requires_fortran"]:
            command.extend(["--binary", str(binary)])
        start = time.monotonic()
        entry = {"id": key, "default_mode": item["default_mode"], "command": command}
        with (output / f"{key}.stdout.log").open("w") as stdout, \
                (output / f"{key}.stderr.log").open("w") as stderr:
            try:
                result = subprocess.run(command, cwd=ROOT, stdout=stdout, stderr=stderr, check=False)
                entry["exit_code"] = result.returncode
                missing = [name for name in item["required_outputs"]
                           if not (output / key / name).is_file()]
                entry["missing_outputs"] = missing
                entry["status"] = "PASS" if result.returncode == 0 and not missing else "FAIL"
                if entry["status"] == "PASS":
                    try:
                        payload = json.loads((output / key / "result.json").read_text())
                        if (not isinstance(payload, dict) or payload.get("schema_version") != 1
                                or payload.get("feature_id") != key or payload.get("status") != "PASS"):
                            raise ValueError("result.json must report schema_version=1, the selected feature_id and status=PASS")
                    except (OSError, ValueError) as exc:
                        entry.update(status="FAIL", error=str(exc))
            except OSError as exc:
                entry.update(status="FAIL", exit_code=None, error=str(exc))
        entry["elapsed_seconds"] = time.monotonic() - start
        entry["stdout"] = f"{key}.stdout.log"
        entry["stderr"] = f"{key}.stderr.log"
        failed |= entry["status"] != "PASS"
        report["features"].append(entry)
        report_path.write_text(json.dumps(report, indent=2) + "\n")
        print(f"{key}: {entry['status']} ({entry['elapsed_seconds']:.2f} s)", flush=True)
    report["status"] = "FAIL" if failed else "PASS"
    report["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Run manifest: {report_path}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
