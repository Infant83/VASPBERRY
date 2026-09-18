#!/usr/bin/env python3
"""Calculate the insulating Bi charge-Hall plateau from its actual VASP WAVECAR."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "kubo-curvature"))
from workflow import ROOT, begin, fail, finish, run_command

INPUT_SHA256 = "a8d81854f2efc561938e478dde1be29a17ccc95d5d37325c122b9d9e82fa0838"


def read_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--wavecar", type=Path, default=ROOT / "examples/Bi_Z2/WAVECAR")
    args = parser.parse_args()
    output, wavecar = args.output_dir.resolve(), args.wavecar.resolve()
    feature = Path(__file__).resolve().parent
    record = begin(output,feature,wavecar,INPUT_SHA256)
    try:
        calculation = output / "calculation"
        run_command([sys.executable,ROOT / "tools/wavecar_fukui.py",wavecar,
            "--nx",12,"--ny",12,"--spinor-components",2,"--energy-band",11,
            "--map","occupied=1:10","--transport-full-t0",10,
            "--mu-min",-1.3,"--mu-max",-.9,"--mu-num",41,
            "--valley-k","0.6666666666666666,0.3333333333333333,0",
            "--valley-kp","0.3333333333333333,0.6666666666666666,0",
            "--output-dir",calculation],ROOT,output,record,wavecar)
        rows = read_rows(calculation / "transport_full_t0.csv")
        diagnostics = json.loads((calculation / "transport_full_t0_diagnostics.json").read_text())
        map_rows = read_rows(calculation / "fukui_occupied.csv")
        summary = [{key:row[key] for key in ("mu_eV","quality","continuous_scan_quality",
            "occupied_band_count_min","occupied_band_count_max","chern_total","sigma_total_e2_over_h")}
            for row in rows]
        with (output / "summary.csv").open("w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(summary[0]));writer.writeheader();writer.writerows(summary)
        mu=np.array([float(row["mu_eV"]) for row in rows])
        sigma=np.array([float(row["sigma_total_e2_over_h"]) for row in rows])
        flux=np.empty((12,12))
        for row in map_rows:
            flux[int(row["ix"]),int(row["iy"])]=float(row["phi_rad"])
        map_chern=float(flux.sum()/(2*np.pi))
        sum_rule_error=max(abs(float(row["sigma_total_e2_over_h"])-sum(float(row[name]) for name in
            ("sigma_K_e2_over_h","sigma_Kp_e2_over_h","sigma_outside_e2_over_h"))) for row in rows)
        reference_error=None
        reference=feature / "reference/summary.csv"
        if reference.exists():
            ref=read_rows(reference)
            if len(ref)!=len(rows) or any(abs(float(a["mu_eV"])-float(b["mu_eV"]))>1e-12 for a,b in zip(ref,rows)):
                raise RuntimeError("Reference chemical-potential grid mismatch")
            reference_error=max(abs(float(a["sigma_total_e2_over_h"])-float(b["sigma_total_e2_over_h"])) for a,b in zip(ref,rows))
        fig,axes=plt.subplots(1,2,figsize=(10,4.1),constrained_layout=True)
        axes[0].axhline(0,color="grey",linestyle="--",label="TRS charge-Hall expectation")
        axes[0].plot(mu,sigma,color="#2166ac",lw=2,label="Actual WAVECAR calculation")
        axes[0].set(xlabel="Chemical potential μ (WAVECAR eV)",ylabel=r"Sheet $\sigma_{xy}$ ($e^2/h$)",
                    title="Bi · occupied bands 1:10 · T = 0 K",ylim=(-.05,.05))
        axes[0].grid(alpha=.2);axes[0].legend(fontsize=8,loc="upper right")
        axes[0].text(.5,.1,f"max |σ| = {np.max(abs(sigma)):.2e} e²/h\nAll 41 points pass the guards",
                     ha="center",transform=axes[0].transAxes,fontsize=9)
        image=axes[1].imshow(flux.T*1e6,origin="lower",extent=(0,1,0,1),cmap="RdBu_r",
                            vmin=-max(abs(flux.min()),abs(flux.max()))*1e6,
                            vmax=max(abs(flux.min()),abs(flux.max()))*1e6,interpolation="nearest")
        axes[1].set(xlabel="Fractional q₁",ylabel="Fractional q₂",title="Occupied-subspace plaquette flux residual")
        fig.colorbar(image,ax=axes[1],label="Flux (10⁻⁶ rad)",shrink=.85)
        fig.suptitle("Bi bilayer · public 12 × 12 VASP WAVECAR · charge response")
        fig.savefig(output / "figure.png",dpi=170);plt.close(fig)
        checks={"mu_count":len(rows),"sampled_global_gap_eV":diagnostics["indirect_gap_above_max_bundle_ev"],
            "valence_max_eV":diagnostics["max_bundle_energy_max_ev"],"conduction_min_eV":diagnostics["sentinel_min_ev"],
            "minimum_direct_occupied_subspace_gap_eV":diagnostics["max_bundle_quality"]["min_gap_to_sentinel_ev"],
            "all_discrete_and_continuous_mu_validated":diagnostics["validated"],
            "all_rows_have_10_occupied_bands":all(int(r["occupied_band_count_min"])==10 and int(r["occupied_band_count_max"])==10 for r in rows),
            "max_abs_sigma_e2_over_h":float(np.max(abs(sigma))),"map_chern":map_chern,
            "map_vs_transport_chern_abs_error":abs(map_chern-diagnostics["max_bundle_chern"]),
            "geometric_partition_sum_error_e2_over_h":sum_rule_error,
            "minimum_plane_wave_coverage":diagnostics["minimum_plane_wave_coverage"],
            "minimum_occupied_subspace_link_sv":diagnostics["max_bundle_quality"]["min_link_singular_value"],
            "max_occupied_flux_abs_rad":float(np.max(abs(flux))),"reference_comparison_max_abs_e2_over_h":reference_error}
        checks["passed"]=bool(len(rows)==41 and diagnostics["validated"] and checks["all_rows_have_10_occupied_bands"]
            and all(row["quality"]=="PASS" and row["continuous_scan_quality"]=="PASS" for row in rows)
            and checks["max_abs_sigma_e2_over_h"]<1e-10 and abs(map_chern)<1e-10
            and checks["map_vs_transport_chern_abs_error"]<1e-10 and sum_rule_error<1e-10
            and (reference_error is None or reference_error<1e-10))
        finish(output,record,{"feature_id":"hall-valley","material":"Bi bilayer",
            "input_repository_path":"examples/Bi_Z2/WAVECAR (Git LFS payload)",
            "numerical_checks":checks,"units":{"energy":"eV (unchanged WAVECAR zero)","sigma":"e^2/h (2D sheet response)","flux":"radian"},
            "method":"Guarded cumulative occupied-subspace Fukui transport; MAX_BAND=10, sentinel=11, T=0",
            "scope":"Actual VASP WAVECAR calculation of the insulating total charge response. "
                "Bi's near-degenerate Kramers partners require occupied subspaces; this is not individual-band point-Kubo integration. "
                "The required K/Kprime centers define a geometric partition only, not a validated physical valley observable. "
                "Regional values are retained for reproducibility; no finite valley or spin Hall effect is claimed. "
                "The 12x12 grid is a regression reference, not a material-convergence study."})
        print(f"PASS: all 41 μ points; max |σ|={np.max(abs(sigma)):.3g} e²/h; {output / 'figure.png'}")
    except Exception as exc:
        fail(output,record,exc)
        raise


if __name__ == "__main__":
    main()
