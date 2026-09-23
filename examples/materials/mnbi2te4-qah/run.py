#!/usr/bin/env python3
"""Run direct Fukui and standard-WAVEDER Hall checks on a completed MnBi2Te4 mesh."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from wavecar_fukui import Wavecar, infer_uniform_grid


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, nargs="+", required=True,
                        help="completed full-mesh VASP optical run, or fixed-charge chunks covering that mesh")
    parser.add_argument("--wavecar", type=Path,
                        help="complete mesh WAVECAR for Fukui; required for several source chunks")
    parser.add_argument("--mesh", type=int, default=6)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    sources = [path.resolve() for path in args.run_dir]
    if len(sources) > 1 and args.wavecar is None:
        parser.error("several optical chunks require a complete --wavecar for the Fukui calculation")
    wavecar_path = args.wavecar.resolve() if args.wavecar else sources[0] / "WAVECAR"
    output = args.output_dir.resolve()
    if args.mesh < 2 or output.exists():
        parser.error("mesh must be at least 2; choose a new output directory")
    wave = Wavecar(wavecar_path)
    if wave.header.ispin != 1 or wave.spinor_components != 2 or wave.header.nbands <= 123:
        parser.error("this material example requires SOC spinors with more than 123 stored bands")
    infer_uniform_grid(wave.kpoints, args.mesh, args.mesh)
    covered = []
    for source in sources:
        part = Wavecar(source / "WAVECAR")
        if (part.header.nbands != wave.header.nbands or part.header.ispin != wave.header.ispin
                or part.header.encut_ev != wave.header.encut_ev
                or not np.array_equal(part.header.lattice, wave.header.lattice)):
            parser.error("the Fukui and optical source wavefunctions must use the same lattice, cutoff, spin and band count")
        for ik, point in enumerate(part.kpoints):
            delta = wave.kpoints - point
            delta -= np.rint(delta)
            matches = np.flatnonzero(np.max(np.abs(delta), axis=1) < 1e-12)
            if len(matches) != 1:
                parser.error("optical source k points do not match the complete Fukui mesh")
            target = int(matches[0])
            if not np.allclose(part.energies[ik], wave.energies[target], atol=1e-8, rtol=0):
                parser.error("Fukui and optical source energies differ")
            if (source / "WAVECAR").resolve() != wavecar_path:
                if (part.nplane[ik] != wave.nplane[target]
                        or not np.allclose(point, wave.kpoints[target], atol=1e-12, rtol=0)
                        or not np.array_equal(part.g_vectors(ik), wave.g_vectors(target))):
                    parser.error("the complete Fukui WAVECAR must preserve each source plane-wave basis")
                for first in range(1, wave.header.nbands + 1, 16):
                    bands = range(first, min(first + 16, wave.header.nbands + 1))
                    if not np.array_equal(part.coefficients(ik, bands), wave.coefficients(target, bands)):
                        parser.error("the complete Fukui WAVECAR must preserve all optical-source coefficients exactly")
            covered.append(target)
    if sorted(covered) != list(range(wave.header.nkpoints)):
        parser.error("optical source runs must cover the complete mesh exactly once")
    vbm = float(wave.energies[:, 122].max())
    cbm = float(wave.energies[:, 123].min())
    if cbm <= vbm:
        parser.error("the supplied calculation has no sampled global gap above band 123")
    output.mkdir(parents=True)
    midgap = (vbm + cbm) / 2
    commands = [
        [sys.executable, str(ROOT / "tools/wavecar_fukui.py"), str(wavecar_path),
         "--nx", str(args.mesh), "--ny", str(args.mesh), "--spinor-components", "2",
         "--map", "occupied=1:123", "--map", "deep_s=1:36", "--energy-band", "124",
         "--output-dir", str(output / "fukui")],
        [sys.executable, str(ROOT / "tools/vaspberry_kubo.py"), "waveder-hall",
         "--run-dir", *map(str, sources), "--occupied", "123", "--spinor-components", "2",
         "--spin-multiplicity", "1", "--mesh", str(args.mesh), str(args.mesh),
         "--energy-reference", "unchanged VASP eigenvalue energy zero",
         "--mu-min", str(vbm + .1 * (cbm - vbm)), "--mu-max", str(cbm - .1 * (cbm - vbm)),
         "--mu-num", "31", "--mu-reference", str(midgap),
         "--formats", "csv", "dat", "npz", "--output-dir", str(output / "waveder")],
    ]
    for name, command in zip(["fukui", "waveder"], commands):
        with (output / f"{name}.stdout.log").open("w") as out, (output / f"{name}.stderr.log").open("w") as err:
            result = subprocess.run(command, stdout=out, stderr=err)
        if result.returncode:
            raise SystemExit(f"{name} failed; inspect its retained stderr log in {output}")
    diagnostics = json.loads((output / "fukui/diagnostics.json").read_text())
    record = {
        "status": "COMPLETED", "material": "MnBi2Te4 three septuple layers",
        "mesh": [args.mesh, args.mesh, 1], "stored_bands": wave.header.nbands,
        "occupied_bands": [1, 123], "sampled_vbm_eV": vbm, "sampled_cbm_eV": cbm,
        "sampled_indirect_gap_eV": cbm - vbm, "reference_midgap_eV": midgap,
        "fukui_maps": diagnostics["maps"],
        "source_state_check": "identical source file or complete plane-wave basis and coefficient equality",
        "waveder_status": "finite-grid Kubo result; assess band and k-mesh convergence independently",
        "commands": commands,
    }
    (output / "result.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"Completed direct Fukui and standard-WAVEDER calculations in {output}.")
    print("An integer Fukui invariant does not establish convergence of the pointwise Kubo integral.")


if __name__ == "__main__":
    main()
