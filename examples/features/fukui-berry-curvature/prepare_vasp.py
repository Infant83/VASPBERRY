#!/usr/bin/env python3
"""Prepare the MoS2 full-mesh NSCF calculation from public charge-density data."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
POTCAR_SHA256 = "b8f5bccb188e7a8b1082445abe307e62f27d3cd9f9f7b128325dacee3e9192a0"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--potcar", type=Path, required=True, help="Licensed concatenated Mo+S PAW-PBE file matching the material specification")
    parser.add_argument("--output-dir", type=Path, required=True, help="New VASP working directory")
    args = parser.parse_args()
    if not args.potcar.is_file() or sha(args.potcar) != POTCAR_SHA256:
        parser.error("POTCAR must match examples/1H-MoS2/PSEUDOPOTENTIAL.md for this reference calculation")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    for name in ("INCAR", "POSCAR", "KPOINTS"):
        shutil.copyfile(HERE / "inputs" / name, output / name)
    shutil.copyfile(args.potcar, output / "POTCAR")
    charge = ROOT / "examples/1H-MoS2/KPATH/1.scf/CHGCAR.gz"
    with gzip.open(charge, "rb") as source, (output / "CHGCAR").open("wb") as target:
        shutil.copyfileobj(source, target)
    manifest = {"schema_version": 1, "status": "PREPARED", "material": "1H-MoS2",
                "method": "SOC ICHARG=11 NSCF on a full12x12 mesh using the supplied public SCF charge density",
                "input_sha256": {name: sha(output / name) for name in ("INCAR", "POSCAR", "KPOINTS", "CHGCAR", "POTCAR")},
                "source_chgcar_gz_sha256": sha(charge), "preparation_script_sha256": sha(__file__)}
    (output / "input_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(f"Prepared {output}. Run your licensed VASP noncollinear executable in this directory.")


if __name__ == "__main__":
    main()
