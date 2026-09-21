#!/usr/bin/env python3
"""Prepare Cartesian 9 x 9 VASP sampling patches around MoS2 K and K'."""
from pathlib import Path
import argparse
import json
import shutil
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools"))
from plot_berry_curvature import reciprocal_from_poscar


def valley_points(reciprocal):
    offsets = np.linspace(-0.12, 0.12, 9)
    delta = np.array([[x, y, 0.] for x in offsets for y in offsets])
    dq = delta @ np.linalg.inv(reciprocal)
    centers = np.array([[1/3, 2/3, 0.], [-1/3, -2/3, 0.]])
    return np.concatenate([center + dq for center in centers])


def prepare(input_dir, output_dir):
    required = ("INCAR", "POSCAR", "POTCAR", "CHGCAR")
    for name in required:
        if not (input_dir / name).is_file():
            raise ValueError(f"missing VASP input: {input_dir / name}")
    points = valley_points(reciprocal_from_poscar(input_dir / "POSCAR"))
    output_dir.mkdir(parents=True, exist_ok=False)
    for name in required:
        shutil.copyfile(input_dir / name, output_dir / name)
    (output_dir / "KPOINTS").write_text(
        "MoS2 K and Kprime Cartesian valley patches; 9x9 each\n162\nReciprocal\n"
        + "".join(" ".join(f"{q:.16f}" for q in point) + " 1\n" for point in points))
    (output_dir / "sampling.json").write_text(json.dumps({
        "valleys": ["K", "Kprime"], "centers_fractional": [[1/3, 2/3, 0], [-1/3, -2/3, 0]],
        "offset_axes": "Cartesian delta kx and delta ky, Angstrom^-1",
        "offset_min": -0.12, "offset_max": 0.12, "points_per_axis": 9,
        "ordering": "valley, delta kx, delta ky (fastest)",
        "input_source": "User-provided converged MoS2 full-mesh NSCF setup; only KPOINTS replaced"
    }, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True,
                        help="full-mesh MoS2 NSCF directory with INCAR/POSCAR/POTCAR/CHGCAR")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        prepare(args.input_dir.resolve(), args.output_dir.resolve())
    except (ValueError, OSError) as error:
        parser.error(str(error))
    print(args.output_dir)


if __name__ == "__main__":
    main()
