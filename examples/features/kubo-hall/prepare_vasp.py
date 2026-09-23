#!/usr/bin/env python3
"""Prepare a MoS2 SOC mesh for Kubo-Hall sampling and band-window checks."""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import shutil

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = HERE.parent / "fukui-berry-curvature/inputs"
POTCAR_SHA256 = "b8f5bccb188e7a8b1082445abe307e62f27d3cd9f9f7b128325dacee3e9192a0"


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--potcar", type=Path, required=True,
                        help="licensed concatenated Mo+S PAW-PBE file matching the example specification")
    parser.add_argument("--mesh", type=int, default=24, help="full Gamma-centered N x N x 1 mesh")
    parser.add_argument("--nbands", type=int, default=60, help="stored SOC bands, including empty intermediate states")
    parser.add_argument("--ediff", type=float, default=1e-8, help="electronic convergence tolerance in eV")
    parser.add_argument("--nelmin", type=int, default=16, help="minimum iterations; occupied-energy convergence alone does not validate empty states")
    parser.add_argument("--output-dir", type=Path, required=True, help="new VASP working directory")
    args = parser.parse_args()
    if not math.isfinite(args.ediff) or args.ediff <= 0 or args.nelmin < 1:
        parser.error("ediff must be finite and positive; nelmin must be at least one")
    source_incar = (SOURCE/"INCAR").read_text()
    nelm_match = re.search(r"(?m)^NELM\s*=\s*(\d+)", source_incar)
    if nelm_match is None:
        parser.error("reference INCAR must explicitly specify NELM")
    nelm = int(nelm_match.group(1))
    if args.nelmin > nelm:
        parser.error(f"nelmin must not exceed reference NELM ({nelm})")
    if args.mesh < 2 or args.nbands < 20:
        parser.error("this MoS2 example requires mesh >= 2 and NBANDS >= 20")
    if not args.potcar.is_file() or sha(args.potcar) != POTCAR_SHA256:
        parser.error("POTCAR must match examples/1H-MoS2/PSEUDOPOTENTIAL.md")
    output = args.output_dir.resolve()
    if output.exists():
        parser.error("output directory exists; choose a new directory")
    incar, count = re.subn(r"(?m)^NBANDS\s*=.*$", f"NBANDS = {args.nbands}",
                           source_incar)
    if count != 1:
        parser.error("reference INCAR must contain exactly one NBANDS setting")
    incar = re.sub(r"(?m)^EDIFF\s*=.*$", f"EDIFF = {args.ediff:g}", incar)
    incar = re.sub(r"(?m)^NELMIN\s*=.*$", "", incar).rstrip() + f"\nNELMIN = {args.nelmin}\n"
    output.mkdir(parents=True)
    (output/"INCAR").write_text(incar)
    shutil.copyfile(SOURCE/"POSCAR", output/"POSCAR")
    shutil.copyfile(args.potcar, output/"POTCAR")
    (output/"KPOINTS").write_text(
        f"MoS2 full Gamma mesh for transport convergence\n0\nGamma\n{args.mesh} {args.mesh} 1\n0 0 0\n")
    charge = ROOT/"examples/1H-MoS2/KPATH/1.scf/CHGCAR.gz"
    with gzip.open(charge,"rb") as source, (output/"CHGCAR").open("wb") as target:
        shutil.copyfileobj(source,target)
    record = {
        "status":"PREPARED", "material":"1H-MoS2", "mesh":[args.mesh,args.mesh,1],
        "nbands":args.nbands, "ediff_eV":args.ediff, "nelmin":args.nelmin, "occupied_bands":[1,18], "encut_eV":400,
        "method":"SOC ICHARG=11 NSCF using the supplied public SCF charge density; ISYM=-1",
        "input_sha256":{name:sha(output/name) for name in ("INCAR","POSCAR","KPOINTS","CHGCAR","POTCAR")},
        "preparation_script_sha256":sha(__file__), "source_charge_gz_sha256":sha(charge),
    }
    (output/"input_manifest.json").write_text(json.dumps(record,indent=2)+"\n")
    print(f"Prepared {output}. Run your licensed VASP noncollinear executable here.")


if __name__ == "__main__":
    main()
