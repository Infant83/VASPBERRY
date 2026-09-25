#!/usr/bin/env python3
"""Prepare an actual MnBi2Te4 SOC optical mesh from the supplied SCF density."""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil

HERE = Path(__file__).resolve().parent


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--potcar", type=Path, required=True,
                        help="licensed concatenated Mn+Bi+Te PAW-PBE datasets matching PSEUDOPOTENTIAL.md")
    parser.add_argument("--mesh", type=int, default=6,
                        help="complete Gamma-centered N x N x 1 mesh")
    parser.add_argument("--nbands", type=int, default=192,
                        help="stored SOC bands; the occupied bundle contains 123 bands")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.mesh < 2 or args.nbands <= 123:
        parser.error("mesh must be at least 2 and NBANDS must exceed 123 occupied spinor bands")
    meta = json.loads((HERE / "inputs/provenance.json").read_text())
    if not args.potcar.is_file() or sha(args.potcar) != meta["potcar_sha256"]:
        parser.error("POTCAR must match the Mn/Bi/Te datasets in PSEUDOPOTENTIAL.md")
    target = args.output_dir.resolve()
    if target.exists():
        parser.error("output directory exists; choose a new directory")
    incar = (HERE / "inputs/INCAR.nscf").read_text()
    incar, count = re.subn(r"(?m)^NBANDS\s*=.*$", f"NBANDS = {args.nbands}", incar)
    if count != 1:
        raise ValueError("reference INCAR must specify NBANDS exactly once")
    target.mkdir(parents=True)
    (target / "INCAR").write_text(incar)
    shutil.copyfile(HERE / "inputs/POSCAR", target / "POSCAR")
    shutil.copyfile(args.potcar, target / "POTCAR")
    (target / "KPOINTS").write_text(
        f"MnBi2Te4 complete Gamma-centered optical mesh\n0\nGamma\n{args.mesh} {args.mesh} 1\n0 0 0\n")
    with gzip.open(HERE / "inputs/scf/CHGCAR.gz", "rb") as source, (target / "CHGCAR").open("wb") as out:
        shutil.copyfileobj(source, out)
    if sha(target / "CHGCAR") != meta["charge_sha256"]:
        raise ValueError("supplied SCF charge density failed its integrity check")
    manifest = {
        "status": "PREPARED", "material": "MnBi2Te4 three septuple layers",
        "geometry": meta["geometry_status"], "mesh": [args.mesh, args.mesh, 1],
        "nbands": args.nbands, "occupied_bands": [1, 123],
        "method": "SOC fixed-charge NSCF, ICHARG=11, ISYM=-1; standard longitudinal PAW optical WAVEDER",
        "input_sha256": {name: sha(target / name) for name in ["INCAR", "POSCAR", "KPOINTS", "CHGCAR", "POTCAR"]},
    }
    (target / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Prepared {target}. Run your licensed noncollinear VASP executable in this directory.")


if __name__ == "__main__":
    main()
