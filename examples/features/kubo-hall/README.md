# MoS₂: intrinsic charge and regional valley Hall response

This example starts from an actual VASP spinor `WAVECAR`. VASPBERRY exports
interband matrix-element pairs and integrates their occupation-weighted Kubo
response while scanning the chemical potential. The calculation uses the
same monolayer 1H-MoS₂ structure and fixed SCF density as the
[Berry-curvature example](../fukui-berry-curvature/README.md).

![MoS2 regional Hall response and k-mesh study](reference/figures/mos2-hall.png)

**Figure.** (a) Cartesian first Brillouin zone and periodic disks of radius
0.35 Å⁻¹ around K and K′; the remaining area is `rest`. (b) The independently
calculated 26-band K–Γ–K′ path for the same structure, 400 eV cutoff and SCF
density; this band plot is an energy reference, not the integration mesh.
(c) Regional charge Hall changes at 300 K.
(d) k-mesh refinement with 60 stored bands and fixed pair cutoff `M=40`. Each calculation uses its own
valence-band maximum (E_v) to align the horizontal axis.

Here **Δσ(μ) = σ(μ) − σ(μref)**, with μref at the middle of the global band gap.
The valley difference is **ΔσK − ΔσK′**, with no factor of one half. Time-reversal
symmetry makes the total charge Hall response vanish; the opposite regional
responses remain visible. These spatially partitioned charge responses depend
on the specified regions and are not a calculation of a separately defined
conserved valley-current operator.

## 1. Prepare and run VASP

From the repository root, use Python with NumPy and Matplotlib, and a built
native VASPBERRY executable. The [build instructions](../../../README.md)
cover serial and MPI builds.

```bash
python3 examples/features/kubo-hall/prepare_vasp.py \
  --potcar /path/to/licensed/Mo-S/POTCAR \
  --mesh 24 --nbands 60 --output-dir results/mos2-24-b60-vasp

cd results/mos2-24-b60-vasp
vasp_ncl > stdout.log 2> stderr.log
cd ../..
```

The preparation script writes ordinary `INCAR`, `POSCAR`, `KPOINTS`, `CHGCAR`
and licensed `POTCAR` files. It reuses the public
[SCF charge density](../../1H-MoS2/KPATH/1.scf/CHGCAR.gz),
[structure](../fukui-berry-curvature/inputs/POSCAR) and
[NSCF settings](../fukui-berry-curvature/inputs/INCAR), setting `NBANDS`, the complete Γ-centered k mesh, `EDIFF=1e-8` and
`NELMIN=16`. Extra stored empty bands and explicit iteration checks are needed
because occupied-energy convergence alone can leave the highest empty states
inaccurate. This command starts from scratch; the reference conditions also record
the warm restarts used for the published calculations. The matching
[PAW dataset specification](../../1H-MoS2/PSEUDOPOTENTIAL.md) identifies the
required Mo/S datasets. `POTCAR` and newly generated `WAVECAR` files are not
distributed here.

Confirm that VASP reached `EDIFF` and ended normally before using its output.
The setup is SOC, `ICHARG=11`, `ISYM=-1`, 400 eV, with 18 occupied spinor bands.
Use a complete two-dimensional mesh; the supplied line-path `WAVECAR` cannot
replace it.

## 2. Calculate the Hall curves

```bash
python3 examples/features/kubo-hall/run.py \
  --wavecar results/mos2-24-b60-vasp/WAVECAR \
  --binary build/vaspberry-gfortran \
  --mesh 24 --pair-band-max 40 --output-dir results/mos2-24-cap40-hall
```

This MoS₂-specific driver identifies the gap above band18 and calls the
general `wavecar-hall` command. For another material, set its own spin convention,
chemical-potential range and regions through the general CLI; see
[applying the workflow to your system](../../APPLY_TO_YOUR_SYSTEM.md). It scans 61
chemical potentials from (E_v-0.20) to (E_v+0.10) eV at 0 and 300 K.
The [region file](regions.json) defines the two disks; a `rest` region closes
the partition. JSON describes this optional integration geometry; the
electronic-structure input remains the VASP `WAVECAR`.

Outputs include:

| Output | Use |
|---|---|
| `calculation/native/PAIRS.csv` | Native matrix-element pair export |
| `calculation/pairs/` | Validated reusable pair cache |
| `calculation/hall/conductivity.csv`, `.dat`, `.npz` | Equivalent numerical Hall tables |
| `calculation/hall/conductivity.json` | Units, sign convention, regions, approximation and degeneracy diagnostics |
| `result.json` (also `results.json`) | Band edges, run status, timings and numerical checks |
| `plot/hall.png`, `.pdf`, `.svg` | The 300 K regional curves aligned to the VBM |

The tables contain **absolute σ and Δσ**, `total`, `K`, `Kprime`, `rest` and
`valley`, plus carrier counts. Conductivity is a two-dimensional sheet response
in (e^2/h) and siemens; no arbitrary slab thickness is applied. The regional
sum and the unhalved valley difference are checked automatically.

To use MPI for native pair export, pass `--binary build/vaspberry-mpi
--mpi-procs 4`. Python integration then uses the exported cache. Additional
chemical-potential scans can reuse that cache with the general `pair-hall`
command; see `python3 tools/vaspberry_kubo.py pair-hall --help`.

## 3. Plot and inspect the data

```bash
python3 tools/plot_hall.py \
  results/mos2-24-cap40-hall/calculation/hall/conductivity.csv \
  --regions total K Kprime valley --temperatures 300 \
  --quantity delta-sigma --output-dir results/mos2-hall-curves
```

The general plotter accepts CSV, DAT or NPZ and writes PNG, PDF and SVG.
Its default horizontal axis is μ − μref. For an axis relative to the VBM,
pass the `vbm_eV` value from `results.json` as `--energy-origin-eV`.

The published composite figures can be regenerated directly from their
compact numerical references using [plot.py](plot.py). The command and measured
convergence values are in [reference/README.md](reference/README.md).

## 4. Check k sampling and the intermediate-band window

`NBANDS` is the number of states stored by VASP. `--pair-band-max M` restricts
the pair sum to bands 1–M while retaining the actual source-band count in the
metadata. Storing additional states above M helps converge the states used in
the sum. The example uses **60 stored bands with M=40** for mesh refinement;
its cutoff study uses several values of M within one 96-band source.

The k mesh controls integration resolution. The reference compares the 300 K
valley-change curves on a common μ − Ev grid, reporting the maximum absolute
change and relative L2 norm of each refinement. A 1% relative-L2 threshold
is used to describe these comparisons. The 0 K finite-mesh steps are retained.
The final source-band counts, retained cutoffs and measured residuals are
listed in [reference/README.md](reference/README.md); increasing a cutoff to
the highest stored band does not by itself establish convergence.

![Pair-cutoff and temperature checks](reference/figures/mos2-hall-checks.png)

An input sanity check was decisive here: a loose occupied-energy stopping
criterion left the highest empty states at k and −k with different energies.
Using a converged lower pair window inside a larger WAVECAR restored the
expected total-charge cancellation. The complete source export is kept in
the cache, and the retained window is always reported explicitly. No response
is forced to zero or averaged with its time-reversed counterpart. VASP
[documents the slower convergence of the highest iteratively calculated states](https://vasp.at/wiki/NBANDS).

This is a rigid-band intrinsic response: changing μ changes occupations, not
the self-consistent potential. The present native operator uses canonical
momentum without PAW, nonlocal or SOC velocity corrections. Increasing the
mesh or `NBANDS` tests numerical convergence of that operator and does not
remove those operator approximations. Numerical energy groups within
10⁻⁷ eV are explicitly coalesced for the pair integration; unchanged input
energies and the resulting energy/occupation shifts are recorded. A 10⁻⁶ eV
sensitivity check accompanies the reference.
