#!/usr/bin/env python3
"""Run fixed-state metallic Hall scans and illustrative regional decomposition."""
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
from model_tools import begin, check_curvature, finish, model_curvature, result_provenance, run_cli, write_json

KB_EV_K = 8.617333262145e-5
E2_OVER_H_S = 1.602176634e-19**2 / 6.62607015e-34


def occupations(energies, mu, temperature):
    if temperature == 0:
        return (energies <= mu).astype(float)
    return 1 / (1 + np.exp(np.clip((energies-mu)/(KB_EV_K*temperature), -700, 700)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path, help="new directory; never overwritten")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    feature = Path(__file__).resolve().parent
    config = json.loads((feature / "input.json").read_text())
    regions = json.loads((feature / "regions.json").read_text())
    manifest = begin(output, feature)
    try:
        data = model_curvature(config, output, manifest)
        curvature_checks, energies, omega = check_curvature(data, config)
        run_cli(["hall", "--curvature", output / "curvature", "--mu-min", config["mu_min_eV"],
            "--mu-max", config["mu_max_eV"], "--mu-num", config["mu_num"],
            "--mu-reference", config["mu_reference_eV"], "--temperatures", *config["temperatures_K"],
            "--regions", feature / "regions.json", "--difference", "A_minus_B:patch_A:patch_B",
            "--band-resolved", "--output-dir", output / "hall"], output, manifest)
        with (output / "hall/conductivity.csv").open(newline="") as stream:
            reader = csv.DictReader(stream)
            columns = reader.fieldnames
            rows = [{key: value if key == "region" else int(value) if key == "band_id" else float(value)
                     for key, value in row.items()} for row in reader]
        meta = json.loads((output / "hall/conductivity.json").read_text())
        # Independent square-cell minimum-image masks: the model reciprocal
        # lattice is 2*pi*I. This does not use the production region helper.
        masks = {"total": np.ones(len(energies))}
        for region in regions["regions"]:
            delta = data["kpoints_fractional"] - region["center_fractional"]
            delta -= np.rint(delta)
            masks[region["name"]] = (np.linalg.norm(2*np.pi*delta, axis=1) <= region["radius_inv_A"]).astype(float)
        masks["rest"] = 1-masks["patch_A"]-masks["patch_B"]
        masks["A_minus_B"] = masks["patch_A"]-masks["patch_B"]
        lookup = {(r["temperature_K"],r["mu_eV"],r["region"],r["band_id"]): r for r in rows}
        oracle_errors = {name: 0. for name in ("sigma_e2_over_h", "delta_sigma_e2_over_h",
                         "electrons_per_cell", "delta_electrons_per_cell", "sigma_S")}
        sum_errors = {"region_partition": 0., "band_sum": 0., "signed_difference": 0.}
        quantities = ("sigma_e2_over_h", "delta_sigma_e2_over_h", "electrons_per_cell", "delta_electrons_per_cell")
        for row in rows:
            temperature, mu, band = row["temperature_K"], row["mu_eV"], row["band_id"]
            occ = occupations(energies, mu, temperature)
            delta_occ = occ - occupations(energies, config["mu_reference_eV"], temperature)
            weights = masks[row["region"]][:,None] * np.ones_like(energies) / len(energies)
            if band:
                weights[:, 2-band] = 0  # band 1 keeps column 0; band 2 keeps column 1.
            expected = dict(sigma_e2_over_h=float(-2*np.pi*np.sum(weights*occ*omega)),
                delta_sigma_e2_over_h=float(-2*np.pi*np.sum(weights*delta_occ*omega)),
                electrons_per_cell=float(np.sum(weights*occ)),
                delta_electrons_per_cell=float(np.sum(weights*delta_occ)))
            expected["sigma_S"] = expected["sigma_e2_over_h"] * E2_OVER_H_S
            for quantity, value in expected.items():
                oracle_errors[quantity] = max(oracle_errors[quantity], abs(row[quantity]-value))
            if row["region"] == "total":
                for quantity in quantities:
                    regional_sum = sum(lookup[(temperature,mu,part,band)][quantity]
                                       for part in ("patch_A","patch_B","rest"))
                    sum_errors["region_partition"] = max(sum_errors["region_partition"], abs(row[quantity]-regional_sum))
            if row["region"] == "A_minus_B":
                for quantity in quantities:
                    contrast = lookup[(temperature,mu,"patch_A",band)][quantity]-lookup[(temperature,mu,"patch_B",band)][quantity]
                    sum_errors["signed_difference"] = max(sum_errors["signed_difference"], abs(row[quantity]-contrast))
            if band == 0:
                for quantity in quantities:
                    band_sum = sum(lookup[(temperature,mu,row["region"],b)][quantity] for b in (1,2))
                    sum_errors["band_sum"] = max(sum_errors["band_sum"], abs(row[quantity]-band_sum))
        gap = lookup[(0.,0.,"total",0)]
        metallic = [r for r in rows if r["temperature_K"] == 0 and r["region"] == "total"
                    and r["band_id"] == 0 and 0 < r["electrons_per_cell"] < 1]
        checks = {"curvature": curvature_checks, "independent_oracle_max_abs_errors": oracle_errors,
            "sum_rule_max_abs_errors": sum_errors, "metallic_zero_T_sample_count": len(metallic),
            "zero_T_mu_zero_sigma_e2_over_h": gap["sigma_e2_over_h"],
            "zero_T_mu_zero_electrons_per_cell": gap["electrons_per_cell"],
            "expected_gap_sigma_e2_over_h": -1., "expected_gap_electrons_per_cell": 1.,
            "highest_exported_band_max_occupation": meta["highest_band_max_occupation"],
            "integration_scope": meta["scope"], "region_points": meta["region_points"],
            "hall_absolute_tolerance_e2_over_h": config["hall_absolute_tolerance_e2_over_h"],
            "occupation_absolute_tolerance": config["occupation_absolute_tolerance"]}
        checks["passed"] = bool(curvature_checks["passed"] and len(metallic) > 0
            and max(oracle_errors.values()) < config["hall_absolute_tolerance_e2_over_h"]
            and max(sum_errors.values()) < config["occupation_absolute_tolerance"]
            and abs(gap["sigma_e2_over_h"]+1) < config["chern_absolute_tolerance"]
            and abs(gap["electrons_per_cell"]-1) < config["occupation_absolute_tolerance"]
            and meta["scope"] == "all_represented_bands")
        mus = sorted({r["mu_eV"] for r in rows})
        selected = {mus[i] for i in np.linspace(0,len(mus)-1,9,dtype=int)} | {0.}
        with (output / "summary.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=columns)
            writer.writeheader()
            writer.writerows(r for r in rows if r["band_id"] == 0 and r["mu_eV"] in selected)
        fig, axes = plt.subplots(1,2,figsize=(10,4.1),constrained_layout=True)
        colors = {"total":"black", "patch_A":"#1368aa", "patch_B":"#bd5b00", "rest":"#4d8e51", "A_minus_B":"#8a4f9e"}
        for name, color in colors.items():
            selected_rows = sorted((r for r in rows if r["region"] == name and r["band_id"] == 0
                                    and r["temperature_K"] == 0), key=lambda r:r["mu_eV"])
            axes[0].plot([r["mu_eV"] for r in selected_rows], [r["sigma_e2_over_h"] for r in selected_rows],
                         label=name, color=color, linewidth=1.6)
        for temperature, style in ((0.,"-"),(300.,"--")):
            selected_rows = sorted((r for r in rows if r["region"] == "total" and r["band_id"] == 0
                                    and r["temperature_K"] == temperature), key=lambda r:r["mu_eV"])
            axes[1].plot([r["mu_eV"] for r in selected_rows], [r["electrons_per_cell"] for r in selected_rows],
                         style, label=f"{temperature:g} K", linewidth=1.7)
        axes[0].set(ylabel=r"Sheet $\sigma_{xy}$ ($e^2/h$)", title="Illustrative regions, T = 0 K")
        axes[1].set(ylabel="Electrons per cell", title="Fixed states: metal → gap")
        for ax in axes:
            ax.axvspan(-1,config["mu_max_eV"],color="grey",alpha=.1,label="Gap range")
            ax.axvline(0,color="grey",linewidth=.7)
            ax.set(xlabel=r"Chemical potential $\mu$ (eV)", xlim=(config["mu_min_eV"],config["mu_max_eV"]))
            ax.grid(alpha=.2)
            ax.legend(fontsize=8)
        fig.suptitle("QWZ analytic model: arbitrary patches, not physical valleys")
        fig.savefig(output / "figure.png",dpi=170)
        plt.close(fig)
        write_json(output / "result.json", {"schema_version":1, "feature_id":"hall-valley",
            "workflow_mode":"analytic_recalculation", "model":config, "regions":regions,
            "numerical_checks":checks, "provenance":result_provenance(output,manifest),
            "formula":"sigma/(e^2/h)=-2*pi*mean_k sum_n f_n Omega_n; a=1 Angstrom, spin multiplicity 1",
            "units":{"sigma":"e^2/h (sheet response)","sigma_S":"siemens","energy":"eV", "curvature":"Angstrom^2"},
            "scope":"Fixed-state analytic-model hole metal and insulating gap. patch_A/patch_B are arbitrary masks, not uniquely defined valleys; A_minus_B is left minus right without a factor of one-half. No material or metallic mesh-convergence claim."})
        if not checks["passed"]:
            raise RuntimeError("Hall oracle or sum-rule checks failed; see result.json")
        status = finish(output, feature, manifest, True)
        print(json.dumps(checks,indent=2))
        return status
    except Exception as exc:
        finish(output, feature, manifest, False, exc)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
