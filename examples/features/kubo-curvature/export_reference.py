#!/usr/bin/env python3
"""Package compact reference files from a successful actual-WAVECAR tutorial run.

Runtime manifests retain every generated file. A distributed reference has its
own hash set for only the files copied here; original full-run hashes remain
under an explicitly named provenance field and are not reference-file paths.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from workflow import sha256, write_json


def export_reference(run_dir, destination):
    run_dir, destination = Path(run_dir), Path(destination)
    result = json.loads((run_dir / "result.json").read_text())
    if result.get("status") != "PASS" or result.get("workflow_mode") != "vasp_wavecar_calculation":
        raise ValueError("only a successful actual-WAVECAR calculation can become a reference")
    feature = result["feature_id"]
    mapping = {name: name for name in ("summary.csv", "figure.png")}
    if feature == "kubo-curvature":
        mapping["native/KUBO.csv"] = "KUBO.csv"
    elif feature == "hall-valley":
        for name in ("transport_full_t0.csv", "transport_full_t0_diagnostics.json", "fukui_occupied.csv"):
            mapping["calculation/" + name] = name
    else:
        raise ValueError("this exporter handles kubo-curvature and hall-valley only")
    full_hashes = result["output_sha256"].copy()
    for source in mapping:
        if full_hashes.get(source) != sha256(run_dir / source):
            raise ValueError("runtime manifest does not match the payload: " + source)
    destination.mkdir(parents=True, exist_ok=False)
    for source, target in mapping.items():
        shutil.copyfile(run_dir / source, destination / target)
    payload_hashes = {target: sha256(destination / target) for target in mapping.values()}
    provenance = json.loads((run_dir / "provenance.json").read_text())
    provenance["original_full_run_output_sha256"] = provenance.pop("output_sha256")
    provenance["distributed_payload_sha256"] = payload_hashes
    provenance["reference_packaging"] = {
        "path_scope": "result.json output_sha256 and provenance.json distributed_payload_sha256 "
                      "both use reference-directory-relative paths to shipped files; "
                      "original_full_run_output_sha256 uses original-run-relative paths to the "
                      "complete calculation record, including files not shipped in this reference",
        "original_run_path_to_distributed_path": mapping,
        "exporter": "examples/features/kubo-curvature/export_reference.py",
        "exporter_sha256": sha256(__file__),
    }
    write_json(destination / "provenance.json", provenance)
    result["original_full_run_output_sha256"] = result.pop("output_sha256")
    result["output_sha256"] = {**payload_hashes, "provenance.json": sha256(destination / "provenance.json")}
    result["reference_packaging"] = provenance["reference_packaging"]
    write_json(destination / "result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True, help="new reference directory")
    args = parser.parse_args()
    export_reference(args.run_dir, args.output_dir)


if __name__ == "__main__":
    main()
