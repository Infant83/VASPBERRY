#!/usr/bin/env python3
"""Prepare the MoS2 K–Gamma–K' NSCF path using the full-mesh VASP inputs."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil

HERE = Path(__file__).resolve().parent


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def incar_settings(path):
    settings = {}
    for line in path.read_text().splitlines():
        for entry in re.split(r";", re.split(r"[!#]", line, maxsplit=1)[0]):
            if "=" in entry:
                key, value = entry.split("=", 1)
                settings[key.strip().upper()] = value.strip().upper()
    return settings


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh-dir", type=Path, required=True,
                        help="Completed 12x12 MoS2 VASP directory containing INCAR, POSCAR, CHGCAR and POTCAR")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="New VASP working directory for the 49-point path")
    args = parser.parse_args()
    source, output = args.mesh_dir.resolve(), args.output_dir.resolve()
    names = ("INCAR", "POSCAR", "CHGCAR", "POTCAR")
    missing = [name for name in names if not (source / name).is_file()]
    if missing:
        parser.error("mesh directory is missing: " + ", ".join(missing))
    if output.exists():
        parser.error("output directory must be new")
    settings = incar_settings(source / "INCAR")
    required_numbers = {"NBANDS": 26, "ENCUT": 400, "ICHARG": 11,
                        "ISTART": 0, "ISYM": -1, "NSW": 0}
    for key, expected in required_numbers.items():
        try:
            value = float(settings.get(key, "nan"))
        except ValueError:
            value = float("nan")
        if value != expected:
            parser.error(f"this reference path requires {key}={expected} in the full-mesh INCAR")
    for key in ("LSORBIT", "LWAVE"):
        if settings.get(key, "").strip(".") not in ("T", "TRUE"):
            parser.error(f"this reference path requires {key}=.TRUE.")
    if sha(source / "POSCAR") != sha(HERE / "inputs/POSCAR"):
        parser.error("the full-mesh POSCAR must match the supplied MoS2 reference structure")
    output.mkdir(parents=True, exist_ok=False)
    for name in names:
        shutil.copyfile(source / name, output / name)
    shutil.copyfile(HERE / "inputs/path/KPOINTS", output / "KPOINTS")
    manifest = {
        "schema_version": 1, "status": "PREPARED", "material": "1H-MoS2",
        "method": "SOC ICHARG=11 NSCF along K-Gamma-Kprime using unchanged full-mesh inputs",
        "path_node_indices_1based": [1, 25, 49],
        "path_nodes_fractional": [[1 / 3, 2 / 3, 0], [0, 0, 0], [-1 / 3, -2 / 3, 0]],
        "path_kpoints": 49, "bands": 26,
        "input_sha256": {name: sha(output / name) for name in (*names, "KPOINTS")},
        "preparation_script_sha256": sha(__file__),
    }
    (output / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Prepared {output}. Run your licensed VASP noncollinear executable in this directory.")


if __name__ == "__main__":
    main()
