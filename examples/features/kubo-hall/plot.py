#!/usr/bin/env python3
"""Plot the MoS2 Hall example and measured k-mesh / band-window changes."""
from __future__ import annotations
import argparse
import hashlib
import json
import logging
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle, Patch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from plot_hall import read_table
from plot_berry_curvature import reciprocal_from_poscar, first_bz_polygon, draw_bz_outline
from plot_berry_panels import read_bands


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    meta = json.loads((path / "results.json").read_text())
    table = path / "calculation/hall/conductivity.csv"
    if not table.exists():
        table = path / "conductivity.csv"
    if meta["status"] != "PASS":
        raise ValueError("only validated completed calculations may be plotted")
    return meta, read_table(table), table


def curve(case, region="valley", temperature=300):
    meta, data, _ = case
    mask = (data["region"] == region) & (data["temperature_K"] == temperature)
    x = data["mu_eV"][mask] - meta["vbm_eV"]
    order = np.argsort(x)
    return x[order], data["delta_sigma_e2_over_h"][mask][order]


def cap(case):
    return case[0].get("pair_band_max", case[0]["nbands"])


def comparison(cases):
    records = []
    for parameter in ("mesh", "pair_band_max"):
        groups = {}
        for case in cases:
            m = case[0]
            # Mesh comparisons keep both stored NBANDS and retained cap fixed.
            # Cutoff comparisons additionally require exactly the same WAVECAR.
            key = (m["nbands"], cap(case)) if parameter == "mesh" else (
                m["nbands"], tuple(m["mesh"]), m["source_wavecar_sha256"])
            groups.setdefault(key, []).append(case)
        for fixed, group in groups.items():
            group.sort(key=lambda c: c[0]["mesh"] if parameter == "mesh" else cap(c))
            for a,b in zip(group[:-1], group[1:]):
                coarse = a[0]["mesh"] if parameter == "mesh" else cap(a)
                fine = b[0]["mesh"] if parameter == "mesh" else cap(b)
                if coarse == fine:
                    continue
                xa,ya = curve(a); xb,yb = curve(b)
                if not np.allclose(xa,xb,rtol=0,atol=1e-10):
                    raise ValueError("convergence curves must use the same mu-VBM grid")
                records.append(dict(parameter=parameter, source_nbands=fixed[0],
                    fixed_value=fixed[1],
                    coarse=a[0]["mesh"] if parameter == "mesh" else cap(a),
                    fine=b[0]["mesh"] if parameter == "mesh" else cap(b), temperature_K=300,
                    max_abs_difference_e2h=float(np.max(abs(ya-yb))),
                    relative_l2_difference=float(np.linalg.norm(ya-yb)/np.linalg.norm(yb)),
                    criterion_relative_l2=.01,
                    below_one_percent=bool(np.linalg.norm(ya-yb)/np.linalg.norm(yb)<.01)))
    return records


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--case", action="append", type=Path, required=True,
                   help="repeat for completed run.py outputs or compact reference case folders")
    p.add_argument("--primary", type=Path, required=True)
    p.add_argument("--mesh-study-source-bands", type=int, default=60)
    p.add_argument("--cutoff-study-source-bands", type=int, default=96)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--bands", type=Path, default=HERE.parent/"fukui-berry-curvature/reference/path/bands.csv")
    p.add_argument("--poscar", type=Path, default=HERE.parent/"fukui-berry-curvature/inputs/POSCAR")
    a=p.parse_args()
    if a.output_dir.exists(): p.error("output directory exists; choose a new directory")
    primary=load(a.primary); cases=[load(path) for path in a.case]
    case_keys=[(tuple(c[0]["mesh"]),c[0]["nbands"],cap(c)) for c in cases]
    if len(set(case_keys)) != len(case_keys):
        p.error("each mesh, source-band count and pair cutoff must occur once")
    kcases=sorted([c for c in cases if cap(c)==40 and c[0]["nbands"]==a.mesh_study_source_bands],key=lambda c:c[0]["mesh"])
    bcases=sorted([c for c in cases if c[0]["mesh"]==[12,12] and c[0]["nbands"]==a.cutoff_study_source_bands],key=lambda c:cap(c))
    if len({c[0]["source_wavecar_sha256"] for c in bcases}) > 1:
        p.error("cutoff comparison must use one identical source WAVECAR")
    if len(kcases) < 2 or len(bcases) < 2:
        p.error("study figures require at least two meshes at one source band count and two cutoffs from one WAVECAR")
    comparisons=comparison(cases)
    reciprocal=reciprocal_from_poscar(a.poscar); polygon=first_bz_polygon(reciprocal)
    q,distance,energies=read_bands(a.bands,reciprocal)
    logging.getLogger("fontTools.ttLib.tables._h_e_a_d").setLevel(logging.ERROR)
    plt.rcParams.update({"font.size":10,"axes.titlesize":11,"axes.labelsize":11,"legend.fontsize":8.5,
        "pdf.fonttype":42,"svg.fonttype":"none"})
    fig,axs=plt.subplots(2,2,figsize=(10.4,8.0),layout="constrained")
    ax=axs[0,0]
    patch=Polygon(polygon,facecolor="#edf0f3",edgecolor="none"); ax.add_patch(patch)
    regions=json.loads((HERE/"regions.json").read_text())["regions"]
    colors={"K":"#0077bb", "Kprime":"#ee7733", "total":"#111111", "valley":"#aa3377", "rest":"#777777"}
    for region in regions:
        for i in range(-2,3):
            for j in range(-2,3):
                center=(np.array(region["center_fractional"])+[i,j,0])@reciprocal
                circle=Circle(center[:2],region["radius_inv_A"],facecolor=colors[region["name"]],alpha=.65,edgecolor="none")
                circle.set_clip_path(patch); ax.add_patch(circle)
    path=q@reciprocal
    ax.plot(path[:,0],path[:,1],"k--",lw=1.0)
    draw_bz_outline(ax,polygon,reciprocal)
    ax.legend(handles=[Patch(facecolor=colors["K"],label="$K$ region"),Patch(facecolor=colors["Kprime"],label="$K'$ region"),Patch(facecolor="#edf0f3",label="rest")],loc="lower center",ncol=3)
    ax.set_title("(a) Periodic valley regions, radius 0.35 Å⁻¹",loc="left")
    ax=axs[0,1]
    # This independently calculated path is a reference for the common fixed-charge setup.
    path_ev=float(energies[:,:18].max())
    for ib in range(14,min(22,energies.shape[1])):
        ax.plot(distance,energies[:,ib]-path_ev,color=colors["K"] if ib<18 else "0.5",lw=1.1)
    ax.axhspan(-.20,.10,color=colors["valley"],alpha=.15,label="chemical-potential scan")
    ax.axhline(0,color="0.5",lw=.6)
    for x in distance[[0,24,48]]: ax.axvline(x,color="0.8",lw=.5)
    ax.set(xlim=(distance[0],distance[-1]),ylim=(-.65,2.05),xticks=distance[[0,24,48]],xticklabels=["$K$",r"$\Gamma$","$K'$"],ylabel="$E-E_v$ (eV)")
    ax.legend(loc="upper center",frameon=False)
    ax.set_title("(b) Bands along the marked path",loc="left")
    ax=axs[1,0]
    for region,label,style in [("K","$K$","-"),("Kprime","$K'$","-"),("total","total","--"),("valley","$K-K'$","-")]:
        x,y=curve(primary,region);ax.plot(x,y,color=colors[region],ls=style,lw=1.65,label=label)
    ax.axhline(0,color="0.75",lw=.5);ax.legend(frameon=False,ncol=2)
    ax.set_title(f"(c) Regional response, 300 K; {primary[0]['mesh'][0]} × {primary[0]['mesh'][1]}, M = {cap(primary)}",loc="left")
    ax.set(xlabel=r"$\mu-E_v$ (eV)",ylabel=r"$\Delta\sigma_{xy}$ ($e^2/h$)",xlim=(-.2,.1))
    ax=axs[1,1]
    for i,c in enumerate(kcases):
        x,y=curve(c);n=c[0]["mesh"][0]
        ax.plot(x,y,lw=1.65,label=f"{n} × {n}",ls=[":","--","-"][i%3])
    ax.axhline(0,color="0.75",lw=.5);ax.legend(frameon=False,title=f"M = 40; {a.mesh_study_source_bands} stored bands")
    ax.set(xlabel=r"$\mu-E_v$ (eV)",ylabel=r"$\Delta\sigma_{xy}^{K}-\Delta\sigma_{xy}^{K'}$ ($e^2/h$)",xlim=(-.2,.1))
    ax.set_title("(d) k-mesh refinement, 300 K",loc="left")
    for ax in axs.flat:ax.tick_params(direction="in",top=True,right=True)
    a.output_dir.mkdir(parents=True)
    for suffix in ("png","pdf","svg"):fig.savefig(a.output_dir/f"mos2-hall.{suffix}",dpi=220)
    plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10.0,4.0),layout="constrained")
    for c in bcases:
        x,y=curve(c);axs[0].plot(x,y,label=f"M = {cap(c)}",lw=1.6)
    axs[0].set_title(f"(a) Pair cutoff: 12 × 12, source {a.cutoff_study_source_bands}, 300 K",loc="left")
    for t,style in [(0,":"),(300,"-")]:
        x,y=curve(primary,temperature=t);axs[1].plot(x,y,ls=style,label=f"{t} K",lw=1.6)
    axs[1].set_title("(b) Temperature and finite-mesh steps",loc="left")
    for ax in axs:
        ax.set(xlabel=r"$\mu-E_v$ (eV)",ylabel=r"$\Delta\sigma_{xy}^{K}-\Delta\sigma_{xy}^{K'}$ ($e^2/h$)",xlim=(-.2,.1));ax.legend(frameon=False);ax.tick_params(direction="in",top=True,right=True)
    for suffix in ("png","pdf","svg"):fig.savefig(a.output_dir/f"mos2-hall-checks.{suffix}",dpi=220)
    plt.close(fig)
    result=dict(status="PASS",temperature_K=300,alignment="each mesh's own VBM; common mu-VBM grid",
        comparison=comparisons,source_table_sha256={f"{c[0]['mesh'][0]}x{c[0]['mesh'][1]}-source{c[0]['nbands']}-cap{cap(c)}":digest(c[2]) for c in cases},
        band_path_sha256=digest(a.bands),structure_sha256=digest(a.poscar),plot_script_sha256=digest(__file__),
        source_sha256={str(f.relative_to(ROOT)):digest(f) for f in
            [Path(__file__).resolve(),ROOT/"tools/plot_hall.py",ROOT/"tools/plot_berry_curvature.py",ROOT/"tools/plot_berry_panels.py"]},
        output_sha256={f.name:digest(f) for f in a.output_dir.glob("*") if f.is_file()})
    (a.output_dir/"plot.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":main()
