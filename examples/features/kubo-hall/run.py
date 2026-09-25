#!/usr/bin/env python3
"""Run the general WAVECAR Hall command for the supplied MoS2 mesh."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from wavecar_fukui import Wavecar


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--wavecar", type=Path, required=True)
    p.add_argument("--binary", type=Path, default=ROOT / "build/vaspberry-gfortran")
    p.add_argument("--mesh", type=int, help="square mesh size; inferred from WAVECAR when omitted")
    p.add_argument("--mpi-procs", type=int, default=1)
    p.add_argument("--pair-band-max", type=int, help="optional intermediate-band cap within the stored source bands")
    p.add_argument("--output-dir", type=Path, required=True)
    a = p.parse_args()
    if a.output_dir.exists():
        p.error("output directory exists; choose a new directory")
    w = Wavecar(a.wavecar, spinor_components=2)
    if a.mesh is None:
        a.mesh = math.isqrt(w.header.nkpoints)
    if a.mesh < 2 or w.header.nkpoints != a.mesh ** 2 or w.header.nbands < 20:
        p.error("the MoS2 example requires a complete N x N mesh and at least 20 bands")
    cells = (w.kpoints[:, :2] - w.kpoints[0, :2]) * a.mesh
    rounded = np.rint(cells)
    if (np.max(abs(cells-rounded)) > 1e-7
            or len(np.unique(rounded.astype(int) % a.mesh, axis=0)) != a.mesh**2
            or np.max(abs(w.kpoints[:, 2]-np.rint(w.kpoints[:, 2]))) > 1e-8):
        p.error("WAVECAR must contain each point of a regular full square mesh in the q3=0 plane")
    ev = float(np.max(w.energies[:, :18]))
    ec = float(np.min(w.energies[:, 18:]))
    if ec <= ev:
        p.error("bands 1–18 must be separated from the empty bands by a global gap")
    ref = (ev + ec) / 2
    out = a.output_dir.resolve()
    out.mkdir(parents=True)
    cli = ROOT / "tools/vaspberry_kubo.py"
    command = [sys.executable, str(cli), "wavecar-hall", "--binary", str(a.binary.resolve()),
        "--wavecar", str(a.wavecar.resolve()), "--spinor-components", "2", "--spin-multiplicity", "1",
        "--mesh", str(a.mesh), str(a.mesh), "--energy-reference", "unchanged VASP eigenvalue zero",
        "--mu-min", str(ev - .20), "--mu-max", str(ev + .10), "--mu-num", "61",
        "--mu-reference", str(ref), "--temperatures", "0", "300", "--regions", str(HERE / "regions.json"),
        "--difference", "valley:K:Kprime", "--degeneracy-policy", "coalesce",
        "--degeneracy-threshold-eV", "1e-7", "--formats", "csv", "dat", "npz",
        "--mpi-procs", str(a.mpi_procs), "--output-dir", str(out / "calculation")]
    if a.pair_band_max is not None:
        command.extend(["--pair-band-max", str(a.pair_band_max)])
    record = dict(schema_version=1, feature_id="kubo-hall", workflow_mode="vasp_wavecar_calculation", status="RUNNING", mesh=[a.mesh, a.mesh], nbands=w.header.nbands,
        pair_band_max=a.pair_band_max or w.header.nbands,
        vbm_eV=ev, cbm_eV=ec, gap_eV=ec-ev, mu_reference_eV=ref,
        mu_minus_vbm_range_eV=[-.20, .10], temperatures_K=[0, 300],
        source_wavecar_sha256=sha(a.wavecar), native_binary_sha256=sha(a.binary),
        source_sha256={str(f.relative_to(ROOT)): sha(f) for f in
            [Path(__file__), HERE / "regions.json", cli, ROOT / "tools/kubo_pairs.py", ROOT / "tools/berry_data.py", ROOT / "tools/plot_hall.py"]},
        command=["python", "tools/vaspberry_kubo.py", *command[2:]],
        interpretation="rigid-band intrinsic regional charge Hall response; valley = K - Kprime without a factor of one half")
    # Public metadata identifies external inputs without publishing workstation paths.
    for i, token in enumerate(record["command"]):
        if token in ("--binary", "--wavecar", "--regions", "--output-dir"):
            record["command"][i+1] = {"--binary": "build/vaspberry-gfortran", "--wavecar": "INPUT/WAVECAR",
                "--regions": "examples/features/kubo-hall/regions.json", "--output-dir": "OUTPUT/calculation"}[token]
    def save():
        for name in ("result.json", "results.json"):
            (out / name).write_text(json.dumps(record, indent=2) + "\n")
    save()
    env = dict(os.environ)
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        env[key] = "1"
    start = time.monotonic()
    try:
        with (out / "stdout.log").open("w") as stdout, (out / "stderr.log").open("w") as stderr:
            result = subprocess.run(command, stdout=stdout, stderr=stderr, env=env)
        record.update(returncode=result.returncode, wall_seconds=time.monotonic()-start)
        if result.returncode:
            raise RuntimeError("wavecar-hall failed; inspect stdout.log and stderr.log")
        with (out / "calculation/hall/conductivity.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        keys = {(float(r["temperature_K"]), r["region"], round(float(r["mu_eV"])-ev, 8)): r for r in rows}
        offsets = np.round(np.linspace(-.2, .1, 61), 8)
        def value(t, region, offset, field):
            return float(keys[t, region, offset][field])
        checks = {"row_count": len(rows),
            "max_abs_total_sigma_e2h": max(abs(float(r["sigma_e2_over_h"])) for r in rows if r["region"] == "total"),
            "max_abs_total_delta_sigma_e2h": max(abs(float(r["delta_sigma_e2_over_h"])) for r in rows if r["region"] == "total")}
        checks["regional_sum_max_error_e2h"] = max(abs(sum(value(t,r,x,"sigma_e2_over_h") for r in ("K","Kprime","rest"))
            -value(t,"total",x,"sigma_e2_over_h")) for t in (0.,300.) for x in offsets)
        checks["valley_difference_max_error_e2h"] = max(abs(value(t,"K",x,"delta_sigma_e2_over_h")-value(t,"Kprime",x,"delta_sigma_e2_over_h")
            -value(t,"valley",x,"delta_sigma_e2_over_h")) for t in (0.,300.) for x in offsets)
        checks["probe_delta_valley_e2h"] = {str(t): {str(x): value(t,"valley",x,"delta_sigma_e2_over_h")
            for x in (-.2,-.15,-.1,-.05,0.,.05)} for t in (0.,300.)}
        if len(rows) != 610 or max(checks[k] for k in ("regional_sum_max_error_e2h", "valley_difference_max_error_e2h")) > 1e-10:
            raise RuntimeError("row-count or region algebra sanity check failed")
        plot_command = [sys.executable, str(ROOT / "tools/plot_hall.py"),
            str(out / "calculation/hall/conductivity.csv"), "--regions", "K", "Kprime", "total", "valley",
            "--temperatures", "300", "--quantity", "delta-sigma", "--energy-origin-eV", str(ev),
            "--energy-label", "μ − Ev (eV)", "--formats", "png", "pdf", "svg", "--output-dir", str(out / "plot")]
        with (out / "plot.stdout.log").open("w") as stdout, (out / "plot.stderr.log").open("w") as stderr:
            plotted = subprocess.run(plot_command, stdout=stdout, stderr=stderr, env=env)
        if plotted.returncode:
            raise RuntimeError("plot_hall failed; inspect plot.stderr.log")
        record["figure_sha256"] = {str(f.relative_to(out)): sha(f) for f in (out / "plot").glob("hall.*")}
        record.update(status="PASS", checks=checks, wall_seconds=time.monotonic()-start,
            output_sha256={str(f.relative_to(out)): sha(f) for f in sorted((out / "calculation/hall").glob("*")) if f.is_file()})
        save()
        print(f"PASS: {a.mesh} x {a.mesh}, {w.header.nbands} bands; outputs in {out}")
    except Exception as exc:
        record.update(status="FAILED", error=str(exc), wall_seconds=time.monotonic()-start)
        save()
        raise


if __name__ == "__main__":
    main()
