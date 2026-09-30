# MoS₂: Kubo curvature of the occupied valence-band bundle

The main figure shows the trace Berry curvature of **occupied bands 1–18**
across the Brillouin zone and along **K–Γ–K′**. It uses the same occupied
space and VASP inputs as the [Fukui method example](../fukui-berry-curvature/).

![MoS2 occupied-bundle Kubo map, bands and symmetry-path curve](reference/bundle/figure.png)

**Figure.** (a) Occupied-bundle curvature from a full 12×12 mesh, displayed
with NumPy periodic bilinear interpolation in Cartesian reciprocal space.
The dashed line marks the path used in both right-hand panels. (b) Actual
VASP bands at 49 path points, with occupied states in blue and empty states
in gray; energies are relative to the valence-band maximum. (c) Bundle Kubo
curvature calculated directly at those path points.

## Method and inputs

For an isolated bundle $`\mathcal V`$, VASPBERRY evaluates the Kubo formula

```math
\Omega^{\mathcal V}_{xy}=-2\,\mathrm{Im}
\sum_{n\in\mathcal V}\sum_{m\notin\mathcal V}
\frac{D^x_{nm}D^y_{mn}}{(E_n-E_m)^2},\qquad D^a=\hbar v_a.
```

This is the sum of the individual-band curvatures wherever those bands are
resolved. Internal pair contributions cancel, so excluding them before
calculating denominators gives a stable bundle result even at internal
band degeneracies. The selected bundle must be separated from excluded
states. See [Wang et al., Eq. (11) and Sec. III D](https://doi.org/10.1103/PhysRevB.74.195118).

The native WAVECAR implementation uses canonical momentum, $`v_a=p_a/m_e`$,
and exports curvature in Å². PAW augmentation and nonlocal/SOC velocity
corrections are not included. The finite stored empty-band window must be
converged for quantitative work.

| Setting | BZ map and matching path |
|---|---|
| VASP input | Two SOC WAVECARs: full mesh and band path |
| Electronic structure | Same public SCF density, structure, potentials, 400 eV cutoff and 26 bands |
| BZ sampling | Full Γ-centered 12×12×1 mesh |
| Path sampling | 49 K–Γ–K′ points; Γ occurs once at index 25 |
| Reciprocal vertices | K = (1/3, 2/3, 0), Γ = (0, 0, 0), K′ = −K |
| Selected bundle | Occupied bands 1–18 |
| External intermediate states | Empty bands 19–26 |
| Minimum bundle-to-external gap | Approximately 1.674 eV |

Use the [full-mesh preparation](../fukui-berry-curvature/inputs/README.md)
and [matching-path preparation](../fukui-berry-curvature/inputs/path/README.md)
to generate `results/mos2-fullmesh-vasp/WAVECAR` and
`results/mos2-path-vasp/WAVECAR` with your licensed VASP and PAW datasets.
Only the k-point list changes between these calculations.

## 1. Calculate the mesh and path with native Fortran

After completing the two VASP preparations above, run from the repository root:

```bash
make serial
repo_dir="$PWD"
mkdir -p results/mos2-kubo-mesh results/mos2-kubo-path
(
  cd results/mos2-kubo-mesh
  "$repo_dir/build/vaspberry" \
    --wavecar "$repo_dir/results/mos2-fullmesh-vasp/WAVECAR" \
    --task kubo --bands 1:18 \
    --curvature-csv KUBO.csv > vaspberry.log
)
(
  cd results/mos2-kubo-path
  "$repo_dir/build/vaspberry" \
    --wavecar "$repo_dir/results/mos2-path-vasp/WAVECAR" \
    --task kubo --bands 1:18 \
    --curvature-csv KUBO.csv > vaspberry.log
)
```

Both calculations read the VASP wavefunctions directly. No Wannierization
or fitted Hamiltonian is required. `--task kubo` uses every stored point and
creates no new sampling. A range such as `--bands 1:18` automatically writes
one trace-curvature row per k point and spin, including the minimum gap to
excluded states. `KUBO.csv` is the default output for a multi-band trace. It rejects
external gaps ≤10⁻⁵ eV before creating the CSV; internal valence-band
degeneracies are allowed. The native outputs are
`results/mos2-kubo-mesh/KUBO.csv` and `results/mos2-kubo-path/KUBO.csv`.
See the [output format](../../../docs/OUTPUT_FORMAT.md) for their columns.

For a single isolated band, use `--bands 18`. To request separate curves
for multiple bands, add `--per-band 1`; every selected band must remain
separated from all other stored bands by more than 10⁻⁵ eV. The
[valley example](valleys/) provides such a region for bands 17 and 18.

For MPI, build with `make mpi` and replace the executable with
`mpiexec -n 2 "$repo_dir/build/vaspberry-mpi"` in each command.

## 2. Plot the native results

```bash
python3 -m pip install -r requirements-transport.txt
python3 tools/plot_berry_panels.py \
  --method kubo-bundle --input results/mos2-kubo-mesh/KUBO.csv \
  --path-input results/mos2-kubo-path/KUBO.csv \
  --path-wavecar results/mos2-path-vasp/WAVECAR \
  --poscar examples/features/fukui-berry-curvature/inputs/POSCAR \
  --occupied 18 --path-node-indices 1 25 49 --path-labels K Gamma Kprime \
  --map-style smooth --display-grid 401 \
  --title '1H-MoS2' --output results/mos2-kubo-panels/figure.png
```

Here `--method kubo-bundle` selects the plotting reader for a trace CSV; it
is not a native calculation option. This Python step reads the two native
curvature tables and the stored VASP band energies; it does not recalculate curvature. It writes the figure,
plotted path values and an unchanged-energy band table. Use a PDF output
suffix for a vector figure.

## Redraw the supplied results

```bash
python3 tools/plot_berry_panels.py \
  --method kubo-bundle \
  --input examples/features/kubo-curvature/reference/bundle/KUBO_mesh.csv \
  --path-input examples/features/kubo-curvature/reference/bundle/KUBO_path.csv \
  --bands-csv examples/features/kubo-curvature/reference/bundle/bands.csv \
  --poscar examples/features/fukui-berry-curvature/inputs/POSCAR \
  --occupied 18 --path-node-indices 1 25 49 --path-labels K Gamma Kprime \
  --map-style smooth --display-grid 401 \
  --title '1H-MoS2' --output results/mos2-kubo-replot.png
```

Use a PDF suffix for a vector figure. `--path-wavecar` can replace
`--bands-csv` for new data. The general plotter accepts other structures,
occupied counts, path vertices and labels; its bundle band table currently
uses spin channel 1. The tutorial runner validates this specific MoS₂ setup.

`--map-style cells` displays the original samples as cells. The smooth
option uses **NumPy periodic bilinear interpolation**, then clips to the
physical first BZ. The 401×401 display grid adds no calculated k points;
raw native CSV values and any physical integrals remain unchanged.
The Kubo path curve uses actual path calculations, independently of the
map interpolation. No smoothing is applied across invalid single-band data.

## Reference results

| Quantity | Occupied bands 1–18 |
|---|---:|
| Ωz(K), full mesh | −13.17809 Å² |
| Ωz(K′), full mesh | +13.17800 Å² |
| Ωz(K), matching path | −13.17813 Å² |
| Ωz(K′), matching path | +13.17815 Å² |
| Valid bundle points | 144/144 on the mesh; 49/49 on the path |

[Panel PNG](reference/bundle/figure.png) · [PDF](reference/bundle/figure.pdf) ·
[Full-mesh Kubo CSV](reference/bundle/KUBO_mesh.csv) ·
[Path Kubo CSV](reference/bundle/KUBO_path.csv) ·
[Band energies](reference/bundle/bands.csv) ·
[Plotted path curve](reference/bundle/path_curvature.csv) ·
[Calculation record](reference/bundle/result.json)

Bands and curvature share cumulative Cartesian path distance; each K–Γ
segment is approximately 1.320704 Å⁻¹. The Fukui method and Kubo formula cover the same
occupied space. Their finite-resolution values can differ: the Fukui method uses
plaquette averages, whereas native Kubo uses point samples, a finite
empty-band window and the canonical-momentum approximation.

## Optional reproduction helper

The following command combines the two native calculations, input checks
and plotting for this specific MoS₂ reference:

```bash
python3 examples/features/kubo-curvature/run_fullmesh.py \
  --wavecar results/mos2-fullmesh-vasp/WAVECAR \
  --path-wavecar results/mos2-path-vasp/WAVECAR \
  --output-dir results/mos2-kubo-checked
```

The helper selects `--bands 1:18` and uses `--map-style smooth` by default.
There is no separate bundle mode. Use the native commands above when selecting
another material, band range or sampling.

## A meaningful single-band example

The [local K/K′ valley example](valleys/) provides two actual 9×9 VASP patches
where band 18 is separated from every other band by at least 0.130 eV.
It contains input preparation, native calculation, smooth maps, line cuts,
band structure and numerical reference files. This region supports an
individual-band interpretation; Γ does not because its partners touch there.

The earlier whole-zone band-18 [cell map](reference/map-path/figure.png)
and [raw results](reference/map-path/result.json) remain as historical data.
Those old calculations retained unresolved individual-band rows and masked
them during plotting. Current native calculations reject a selected single
band or a `--per-band 1` request when any band gap is ≤10⁻⁵ eV; they do not
enlarge the selected band range. Use the local valley example for separate
band curves, or the isolated occupied range for a full path or mesh trace.

## Supplied 32-band path example

The original downloadable [48-point WAVECAR](../../1H-MoS2/KPATH/2.band/WAVECAR)
is ready for a direct calculation without a new VASP run. Use the complete
occupied bundle so that its internal degeneracy at Γ is handled correctly:

```bash
mkdir -p results/mos2-path
build/vaspberry --task kubo \
  --wavecar examples/1H-MoS2/KPATH/2.band/WAVECAR \
  --bands 1:18 \
  --curvature-csv results/mos2-path/KUBO.csv
```

`KUBO.csv` has one occupied-bundle row per stored k point. This 32-band input
has 48 points, including Γ twice; it is distinct from the 26-band, 49-point
map/path calculation above. Inspect the native CSV directly. The composite
plotter above requires a full mesh and a matching path, so it is not a
standalone plotting command for this line file.

To calculate and draw this supplied path in one step:

```bash
python3 examples/features/kubo-curvature/run.py \
  --binary build/vaspberry --output-dir results/mos2-kubo-supplied-path
```

The helper checks all 48 trace rows against the WAVECAR coordinates and the
occupied-to-empty gap. It plots bands 17–19 for energy context and the
occupied-space trace curvature, including the two stored Γ points.

### Historical band-resolved path reference

[Original path figure](reference/figure.png) · [PDF](reference/figure.pdf) ·
[Summary](reference/summary.csv) · [Native Kubo CSV](reference/KUBO.csv).
The archived calculation contains bands 17 and 18 separately: each has 44
isolated points and four unresolved points on the supplied 48-point path.
Its band-18 valley values are about ±6.67 Å² with intermediate bands 1–32.
The saved raw CSV and execution provenance remain unchanged; they are not
outputs of the current helper and are not a benchmark for its different,
occupied-space trace. Current `--per-band 1` correctly rejects this whole
path at the unresolved points. See the valley tutorial for valid individual
band calculations.

The supplied 32-band path and the matched 26-band map/path above are separate
VASP inputs. A band path alone cannot determine a Brillouin-zone Chern number
or Hall conductivity; see the [Bi Hall example](../hall-valley/) for a full
occupied-space integral.
