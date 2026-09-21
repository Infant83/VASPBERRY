# MoS₂: Kubo curvature of the occupied valence-band bundle

The main figure shows the trace Berry curvature of **occupied bands 1–18**
across the Brillouin zone and along **K–Γ–K′**. It uses the same occupied
space and VASP inputs as the [Fukui example](../fukui-berry-curvature/).

![MoS2 occupied-bundle Kubo map, bands and symmetry-path curve](reference/bundle/figure.png)

**Figure.** (a) Occupied-bundle curvature from a full 12×12 mesh, displayed
with NumPy periodic bilinear interpolation in Cartesian reciprocal space.
The dashed line marks the path used in both right-hand panels. (b) Actual
VASP bands at 49 path points, with occupied states in blue and empty states
in gray; energies are relative to the valence-band maximum. (c) Bundle Kubo
curvature calculated directly at those path points.

## Method and inputs

For an isolated bundle $\mathcal V$, VASPBERRY evaluates

$$
\Omega^{\mathcal V}_{xy}=-2\,\mathrm{Im}
\sum_{n\in\mathcal V}\sum_{m\notin\mathcal V}
\frac{D^x_{nm}D^y_{mn}}{(E_n-E_m)^2},\qquad D^a=\hbar v_a.
$$

This is the sum of the individual-band curvatures wherever those bands are
resolved. Internal pair contributions cancel, so excluding them before
calculating denominators gives a stable bundle result even at internal
band degeneracies. The selected bundle must be separated from excluded
states. See [Wang et al., Eq. (11) and Sec. III D](https://doi.org/10.1103/PhysRevB.74.195118).

The native WAVECAR implementation uses canonical momentum, $v_a=p_a/m_e$,
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

## Calculate and plot

From the repository root:

```bash
make serial
python3 -m pip install -r requirements-transport.txt
python3 examples/features/kubo-curvature/run_fullmesh.py \
  --wavecar results/mos2-fullmesh-vasp/WAVECAR \
  --path-wavecar results/mos2-path-vasp/WAVECAR \
  --output-dir results/mos2-kubo-panels
```

The runner defaults to `--mode bundle --map-style smooth`. It executes
VASPBERRY, validates coordinates and external gaps against both WAVECARs,
and exports unchanged VASP energies, native CSV results, PNG/PDF figures
and execution records. The prepared-input workflow takes about four seconds
on the reference machine.

The native command in each calculation directory is:

```bash
/path/to/vaspberry-gfortran \
  -f /path/to/WAVECAR -s 2 -kubo 2 -kubo_bundle 1 -ii 1 -if 18 \
  -kubo_csv KUBO.csv > vaspberry.log
```

Run it separately for the mesh and path WAVECARs. `-kubo 2` uses their actual
stored points; it creates no new sampling. `-kubo_bundle 1` writes one
trace-curvature row per k point and spin, with its minimum external gap.
It rejects external gaps ≤10⁻⁵ eV before creating the CSV. Internal
valence-band degeneracies are allowed. This option writes a standalone
bundle CSV; see the [output format](../../../docs/OUTPUT_FORMAT.md).
For MPI, use `make mpi` and prefix the MPI executable with `mpiexec -n 2`.

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
segment is approximately 1.320704 Å⁻¹. Fukui and Kubo now cover the same
occupied space. Their finite-resolution values can differ: Fukui uses
plaquette averages, whereas native Kubo uses point samples, a finite
empty-band window and the canonical-momentum approximation.

## A meaningful single-band example

The [local K/K′ valley example](valleys/) provides two actual 9×9 VASP patches
where band 18 is separated from every other band by at least 0.130 eV.
It contains input preparation, native calculation, smooth maps, line cuts,
band structure and numerical reference files. This region supports an
individual-band interpretation; Γ does not because its partners touch there.

The earlier whole-zone band-18 [cell map](reference/map-path/figure.png)
and [raw results](reference/map-path/result.json) remain available. To
reproduce that version, add `--mode band --map-style cells` to the runner
above. It masks unresolved individual bands and is not the occupied-bundle
quantity used in the main figure.

## Supplied 32-band path example

The original downloadable [48-point WAVECAR](../../1H-MoS2/KPATH/2.band/WAVECAR)
and its 32-band Kubo calculation remain available:

```bash
python3 examples/features/kubo-curvature/run.py \
  --wavecar examples/1H-MoS2/KPATH/2.band/WAVECAR \
  --output-dir results/mos2-kubo-supplied-path
```

[Original path figure](reference/figure.png) · [PDF](reference/figure.pdf) ·
[Summary](reference/summary.csv) · [Native Kubo CSV](reference/KUBO.csv).
It gives band-18 valley values about ±6.67 Å² with intermediate bands 1–32.
That stored result is separate from the matched 26-band map/path shown above.
A band path alone cannot determine a Brillouin-zone Chern number or Hall
conductivity; see the [Bi Hall example](../hall-valley/) for a full occupied-space integral.
