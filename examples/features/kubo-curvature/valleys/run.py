#!/usr/bin/env python3
"""Calculate and plot isolated MoS2 valence-band curvature near K and K'."""
from __future__ import annotations
import argparse
import csv
import datetime
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools"))
from wavecar_fukui import Wavecar
from prepare import valley_points


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bilinear(data, x, y):
    """Bounded bilinear interpolation inside a complete Cartesian valley patch."""
    data = np.asarray(data, float)
    x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
    u = (x + .12) / .03
    v = (y + .12) / .03
    if np.any(u < -1e-10) or np.any(u > 8+1e-10) or np.any(v < -1e-10) or np.any(v > 8+1e-10):
        raise ValueError("interpolation cannot leave the sampled valley patch")
    i = np.clip(np.floor(u).astype(int), 0, 7)
    j = np.clip(np.floor(v).astype(int), 0, 7)
    a, b = np.clip(u-i, 0, 1), np.clip(v-j, 0, 1)
    return ((1-a)*(1-b)*data[i,j] + a*(1-b)*data[i+1,j]
            + (1-a)*b*data[i,j+1] + a*b*data[i+1,j+1])


def plot(data, output):
    offsets = np.linspace(-.12, .12, 9)
    if len(data) != 162:
        raise ValueError("expected 162 valley samples")
    for k, row in enumerate(data):
        expected_label = "K" if k < 81 else "Kprime"
        expected_xy = [offsets[(k % 81)//9], offsets[k % 9]]
        if row["k_index"] != k+1 or row["valley"] != expected_label:
            raise ValueError("sample ordering differs from the prepared valley mesh")
        if not np.allclose([row["delta_kx_inv_A"], row["delta_ky_inv_A"]], expected_xy, atol=1e-12, rtol=0):
            raise ValueError("sample coordinates differ from the prepared Cartesian mesh")
        if not np.isfinite([row[key] for key in ("E17_eV", "E18_eV", "E19_eV", "min_gap18_eV", "omega_z_A2")]).all():
            raise ValueError("nonfinite valley data")
        if row["min_gap18_eV"] <= .01 or not row["E17_eV"] < row["E18_eV"] < row["E19_eV"]:
            raise ValueError("band 18 must be isolated throughout both patches")
    fine = np.linspace(-.12, .12, 401)
    xx, yy = np.meshgrid(fine, fine)
    omega = np.array([r["omega_z_A2"] for r in data]).reshape(2,9,9)
    energies = np.array([[r["E17_eV"],r["E18_eV"]] for r in data]).reshape(2,9,9,2)
    bound = float(np.max(np.abs(omega)))
    zero = float(np.max(energies[:,:,:,1]))
    with plt.rc_context({"font.size":10, "axes.titlesize":11, "axes.linewidth":.8}):
        fig, axes = plt.subplots(2,2,figsize=(8.2,7.1),layout="constrained")
        for valley,label in enumerate(("K", "K′")):
            ax=axes[0,valley]
            artist=ax.pcolormesh(fine,fine,bilinear(omega[valley],xx,yy),
                shading="auto",cmap="RdBu_r",vmin=-bound,vmax=bound,rasterized=True)
            ax.axhline(0,color="#202020",ls="--",lw=1.)
            ax.scatter([0],[0],s=22,facecolors="white",edgecolors="#202020",zorder=3)
            ax.set(xlabel=r"$\delta k_x$ ($\AA^{-1}$)",ylabel=r"$\delta k_y$ ($\AA^{-1}$)",
                title=f"({'ab'[valley]}) Band 18 near {label}",aspect="equal",
                xlim=(-.12,.12),ylim=(-.12,.12))
            ax.set_xticks([-.1,0,.1]);ax.set_yticks([-.1,0,.1])
        bar=fig.colorbar(artist,ax=axes[0,:],shrink=.85,aspect=25,pad=.025)
        bar.set_label(r"$\Omega_{18,z}$ ($\AA^2$)")
        for valley,label in enumerate(("K", "K′")):
            for band,color in ((0,"#2166ac"),(1,"#b35806")):
                axes[1,0].plot(fine,bilinear(energies[valley,:,:,band],fine,0)-zero,
                    color=color,ls="-" if valley==0 else "--",lw=1.5,
                    label=f"{label}, band {band+17}")
            color=("#b35806","#2166ac")[valley]
            axes[1,1].plot(fine,bilinear(omega[valley],fine,0),color=color,lw=1.7,label=label)
            axes[1,1].plot(offsets,omega[valley,:,4],"o",ms=3,color=color)
        axes[1,0].set(ylabel=r"$E-E_{\rm v}$ (eV)",title="(c) Valence bands along the dashed cuts")
        axes[1,1].set(ylabel=r"$\Omega_{18,z}$ ($\AA^2$)",title="(d) Curvature along the dashed cuts")
        for ax in axes[1]:
            ax.set(xlabel=r"$\delta k_x$ ($\AA^{-1}$), $\delta k_y=0$",xlim=(-.12,.12))
            ax.axvline(0,color=".65",lw=.6,zorder=0)
            ax.legend(frameon=False,fontsize=8,ncol=2 if ax is axes[1,0] else 1)
        for ax in axes.ravel():ax.tick_params(direction="in",right=True,top=True)
        fig.suptitle("1H-MoS₂: isolated upper valence band near K and K′",fontsize=12)
        for extension in ("png","pdf"):fig.savefig(output/f"figure.{extension}",dpi=250)
        plt.close(fig)
    with (output/"line.csv").open("w",newline="") as handle:
        writer=csv.writer(handle)
        writer.writerow(["valley","delta_kx_inv_A","omega18_z_A2","E17_eV","E18_eV"])
        for valley,label in enumerate(("K","Kprime")):
            for x,omega_value,e17,e18 in zip(fine,bilinear(omega[valley],fine,0),
                    bilinear(energies[valley,:,:,0],fine,0),bilinear(energies[valley,:,:,1],fine,0)):
                writer.writerow([label,x,omega_value,e17,e18])
    return {"map_and_cut_interpolation":"NumPy bilinear; same interpolant for map and line",
        "display_grid_per_patch":[401,401],"raw_grid_per_patch":[9,9],
        "energy_reference_eV":zero,"energy_reference":"maximum E18 on the two sampled patches",
        "color_bound_A2":bound,"native_omega18_range_A2":[float(omega.min()),float(omega.max())]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wavecar",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--binary",type=Path,default=ROOT/"build/vaspberry-gfortran")
    args=parser.parse_args()
    source=args.wavecar.resolve();binary=args.binary.resolve();output=args.output_dir.resolve()
    if not source.is_file() or not binary.is_file():parser.error("WAVECAR and executable must exist")
    if output.exists():parser.error("output directory exists; choose a new directory")
    output.mkdir(parents=True)
    record={"status":"RUNNING","created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "workflow":"VASP WAVECAR -> native VASPBERRY Kubo -> band-18 valley maps",
        "input_source":"MoS2 SOC valley NSCF prepared from public full-mesh example",
        "source_wavecar_sha256":sha(source),"binary_sha256":sha(binary)}
    def save():(output/"result.json").write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
    save()
    try:
        w=Wavecar(source,spinor_components=2)
        expected=valley_points(w.header.reciprocal)
        if w.energies.shape != (162,26):raise ValueError("this example requires 162 k points and 26 bands")
        if not np.allclose(w.kpoints,expected,atol=1e-12,rtol=0):raise ValueError("WAVECAR coordinates do not match the valley sampling")
        if abs(w.header.encut_ev-400)>1e-10:raise ValueError("reference setup requires ENCUT=400 eV")
        native=output/"native";native.mkdir()
        command=[str(binary),"-f",str(source),"-s","2","-kubo","2","-ii","17","-if","18",
                 "-kubo_csv","KUBO.csv","-o","BERRYCURV"]
        record["command"]=["vaspberry", "-f", "WAVECAR", *command[3:]]
        start=time.monotonic()
        with (output/"native.stdout.log").open("w") as stdout,(output/"native.stderr.log").open("w") as stderr:
            process=subprocess.run(command,cwd=native,stdout=stdout,stderr=stderr,timeout=300)
        record.update(native_returncode=process.returncode,native_seconds=time.monotonic()-start)
        if process.returncode:raise RuntimeError("native Kubo calculation failed; see native.stderr.log")
        shutil.copyfile(native/"KUBO.csv",output/"KUBO.csv")
        with (output/"KUBO.csv").open() as stream:
            lines=stream.readlines()
        metadata=dict(line[1:].strip().split("=",1) for line in lines if line.startswith("#") and "=" in line)
        if metadata.get("normalization")!="STANDARD_MINUS_TWO_IM" or metadata.get("intermediate_bands")!="1:26":
            raise ValueError("native Kubo normalization or intermediate-band range differs from the reference")
        rows=list(csv.DictReader(line for line in lines if not line.startswith("#")))
        indexed={(int(r["k_index"]),int(r["band"])):r for r in rows}
        if len(indexed)!=324 or len(rows)!=324:raise ValueError("unexpected or duplicate native Kubo rows")
        if set(indexed)!={(k,b) for k in range(1,163) for b in (17,18)}:raise ValueError("incomplete native output")
        summary=[];max_gap_error=0.;max_energy_error=0.;max_coordinate_error=0.
        offsets=np.linspace(-.12,.12,9)
        with (output/"bands.csv").open("w",newline="") as handle:
            writer=csv.writer(handle);writer.writerow(["k_index","q1","q2","q3","band","energy_eV"])
            for k in range(162):
                for b in range(26):writer.writerow([k+1,*w.kpoints[k],b+1,w.energies[k,b]])
        for k in range(162):
            row=indexed[k+1,18]
            gap=float(np.min(np.abs(np.delete(w.energies[k],17)-w.energies[k,17])))
            max_gap_error=max(max_gap_error,abs(gap-float(row["min_gap_eV"])))
            max_energy_error=max(max_energy_error,abs(w.energies[k,17]-float(row["energy_eV"])))
            max_coordinate_error=max(max_coordinate_error,float(np.max(np.abs(w.kpoints[k]-[float(row[f"k{a}_frac"]) for a in "xyz"]))))
            omega=float(row["omega_z_A2"])
            if not np.isfinite(omega) or gap<=.01:raise ValueError("band 18 is not reliably isolated throughout the patch")
            summary.append({"k_index":k+1,"valley":"K" if k<81 else "Kprime",
                "delta_kx_inv_A":offsets[(k%81)//9],"delta_ky_inv_A":offsets[k%9],
                "E17_eV":w.energies[k,16],"E18_eV":w.energies[k,17],"E19_eV":w.energies[k,18],
                "min_gap18_eV":gap,"omega_z_A2":omega})
        if max_gap_error>1e-10 or max_energy_error>1e-10 or max_coordinate_error>1e-10:
            raise ValueError("native CSV disagrees with independent WAVECAR metadata")
        with (output/"summary.csv").open("w",newline="") as handle:
            writer=csv.DictWriter(handle,fieldnames=list(summary[0]));writer.writeheader();writer.writerows(summary)
        omega=np.array([r["omega_z_A2"] for r in summary]).reshape(2,9,9)
        residual=float(np.max(np.abs(omega[0]+omega[1,::-1,::-1])))
        if residual>1e-3:raise ValueError("time-reversed valley curvature residual exceeds 1e-3 Angstrom^2")
        record.update(status="PASS",material="1H-MoS2",target_band=18,other_exported_band=17,
            intermediate_bands=[1,26],sampling="9 x 9 Cartesian patches at each of K and Kprime",
            numerical_checks={"passed":True,"minimum_band18_isolation_gap_eV":min(r["min_gap18_eV"] for r in summary),
                "minimum_required_gap_eV":.01,"isolated_points":162,"gap_agreement_max_eV":max_gap_error,
                "energy_agreement_max_eV":max_energy_error,"coordinate_agreement_max_fractional":max_coordinate_error,
                "time_reversed_valley_residual_max_A2":residual,
                "K_Kprime_omega18_A2":[float(omega[0,4,4]),float(omega[1,4,4])]},
            figure_conventions=plot(summary,output),
            scope="Locally isolated single-band curvature near the valleys; band 18 is not isolated across the full BZ. "
                "Canonical-momentum Kubo with 26 intermediate bands omits PAW augmentation and nonlocal/SOC velocity terms. "
                "Interpolation changes display only; these patches do not determine a full-BZ Chern number or conductivity.")
        record["output_sha256"]={name:sha(output/name) for name in ("KUBO.csv","bands.csv","summary.csv","line.csv","figure.png","figure.pdf")}
        save();print(json.dumps(record["numerical_checks"],indent=2))
    except Exception as error:
        record.update(status="FAILED",error=str(error));save();raise


if __name__ == "__main__":main()
