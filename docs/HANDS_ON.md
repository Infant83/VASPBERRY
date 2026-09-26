# Hands-on: native calculation, saved data, and reusable plots

The main workflow is **VASP → WAVECAR → VASPBERRY Fortran with Intel MPI → numerical files → analysis and plots**. Native curvature is already a numerical result that you can plot in your preferred program. A charge-Hall scan additionally integrates Kubo pairs with occupations.

For routine Hall and PROCAR analysis, start with the [single-file postprocessing guide](POSTPROCESSING.md) and [public Bi example](../examples/features/simple-postprocess/README.md): `vaspberry_post.py run analysis.ini` calculates the requested tables and `vaspberry_post.py plot results/run01` draws them. The explicit stages below explain the native outputs and how to use them independently.

Use VASPBERRY 1.5.0 for this walkthrough. Start in the repository root in Bash on a Linux host with Intel oneAPI Fortran, oneMKL and Intel MPI available. The following uses four MPI ranks; select the rank count allowed by your cluster allocation. Choose a fresh output directory when repeating a calculation.

```bash
source /opt/intel/oneapi/setvars.sh
make ifx-mpi
repo_dir="$PWD"
vb_bin="$repo_dir/build/vaspberry-ifx-mpi"
mpiexec -n 4 "$vb_bin" --help
```

Use your site's oneAPI setup path or modules when they differ. A retained Intel Classic installation can use `make ifort-mpi` and `vb_bin="$repo_dir/build/vaspberry-ifort-mpi"`. [Build details and GNU alternatives](BUILD.md) include the byte-RECL input requirement and compiler validation status. [Native commands](NATIVE_COMMANDS.md) explain arguments; the [report-to-example map](../examples/REPORT_REPRODUCTION.md) links material inputs and figures.

| Stage | Reads → writes | Purpose |
|---|---|---|
| Native `--task kubo` | WAVECAR → `KUBO.csv` | Calculate point curvature Ωxy in Å²; ready for a curvature plot |
| Native `--task kubo-pairs` | WAVECAR → `PAIRS.csv` | Calculate raw interband numerators N in eV² Å²; reusable transport input |
| Python `import-pairs` | `PAIRS.csv` + same WAVECAR → `pairs.npz/json` | Validate, organize and cache those numbers; no Hall integration |
| Python `pair-hall` | pair cache + μ/T/regions → `conductivity.csv/dat/npz/json` | Calculate the occupation-weighted Kubo transport integral |
| `plot_hall.py` or another plotter | completed conductivity table → figures | Plot existing results; no new curvature or transport calculation |

## 1. Calculate a real supplied WAVECAR with Fortran

This first calculation uses the ordinary SOC MoS₂ band-path WAVECAR already in the repository. It evaluates occupied-bundle point curvature, with bands 1–18 and the stored empty bands as intermediate states.

```bash
mkdir -p results/first-kubo
(
  cd results/first-kubo
  mpiexec -n 4 "$vb_bin" --task kubo \
    --wavecar "$repo_dir/examples/1H-MoS2/KPATH/2.band/WAVECAR" \
    --spinor 2 --bands 1:18 --bundle 1 --curvature-csv KUBO.csv \
    > vaspberry.log 2> vaspberry.err
)
```

| Argument | Choice made here |
|---|---|
| `--task kubo` | Calculate point curvature at the k points already in WAVECAR |
| `--wavecar` | Read this VASP wavefunction file |
| `--spinor 2` | Use its two-component SOC states; no extra spin factor of two |
| `--bands 1:18 --bundle 1` | Trace the isolated group as a whole, allowing internal degeneracies |
| `--curvature-csv KUBO.csv` | Save the numerical table separately from the console log |

The input WAVECAR supplies eigenvalues, plane-wave coefficients, k points and the lattice. `results/first-kubo/KUBO.csv` contains 48 rows, with `spin,k_index,kx_frac,ky_frac,kz_frac,omega_z_A2,min_external_gap_eV`. Here `omega_z_A2` is the already calculated occupied-bundle Ωxy in Å²; the gap is in eV. Comment lines beginning `#` record the operator, normalization, band selection and PASS status. Keep the CSV and native log.

This input is a line path: it supports a curvature-versus-path-sample plot. A full-zone map or Hall integral needs the separately supplied full-mesh recipe.

## 2. Plot the saved CSV with your own tool

**Python is optional for plotting this native result.** In Origin, a spreadsheet or another CSV plotter: import comma-separated values, skip `#` metadata lines, use the next line as column names, select `spin=1`, then plot `k_index` on x and `omega_z_A2` on y. That gives the curvature along the ordered input samples. Keep the metadata beside any exported table.

The equivalent Python/Matplotlib example follows. It only reads the finished CSV and saves `curvature.png`, `.pdf` and `.svg`; no VASPBERRY module or Fortran execution is involved. Use your usual Python 3.10+ environment with Matplotlib. The [dependency file](../requirements-transport.txt) lists supported versions.

```bash
python3 - <<'PY'
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

source = Path('results/first-kubo/KUBO.csv')
with source.open() as stream:
    rows = list(csv.DictReader(line for line in stream if not line.startswith('#')))
rows = [r for r in rows if int(r['spin']) == 1]
fig, ax = plt.subplots()
ax.plot([int(r['k_index']) for r in rows],
        [float(r['omega_z_A2']) for r in rows])
ax.set(xlabel='Ordered k-point index', ylabel='Bundle curvature (Å²)')
fig.tight_layout()
for extension in ('png', 'pdf', 'svg'):
    fig.savefig(source.parent / ('curvature.' + extension))
PY
```

Other tools can read the same CSV after skipping `#` metadata lines. Preserve those lines alongside exported data. The [output specification](OUTPUT_FORMAT.md) describes the columns, axes and units. For a physical Cartesian BZ map, matching band path and symmetry labels, use the [full-mesh/path MoS₂ tutorial](../examples/features/kubo-curvature/); its inputs and plot command reproduce report Figure 2. Display interpolation changes the drawing only.

## 3. Reintegrate saved Kubo data without repeating Fortran

The public MoS₂ cache below contains native **canonical-momentum** pair numerators on a 12×12 mesh, with 60 stored bands. It was saved during the optional matched-operator study; reusing this canonical cache requires neither that study's PAW producer nor VASP. The small scan demonstrates the general integration interface; it is not the full mesh-convergence study.

```bash
python3 tools/vaspberry_kubo.py pair-hall \
  --pairs-dir examples/features/kubo-hall/operator-comparison/reference/canonical-pairs \
  --pair-band-max 40 \
  --mu-min -1.47 --mu-max -1.17 --mu-num 11 --mu-reference -0.4381 \
  --temperatures 0 300 \
  --regions examples/features/kubo-hall/regions.json --difference valley:K:Kprime \
  --degeneracy-policy coalesce --formats csv dat npz \
  --output-dir results/saved-hall

python3 tools/plot_hall.py results/saved-hall/conductivity.csv \
  --regions total K Kprime valley --temperatures 300 --quantity delta-sigma \
  --output-dir results/saved-hall-plot
```

`pair-hall` performs the **transport calculation**: it applies Fermi occupations, energy denominators, k weights and the BZ/region integral. It writes σ in e²/h and siemens, Δσ relative to μref, and represented electron counts. `plot_hall.py` then reads those finished rows and writes `results/saved-hall-plot/hall.png`, `.pdf` and `.svg`: the selected 300 K regional Δσ curves versus μ−μref. The `.json` sidecar records units, input hashes and integration choices. The explicit `coalesce` policy approximates tiny numerical energy splittings by a group mean and records the shifts.

To study another temperature or μ range, rerun only `pair-hall` into a new directory, then plot it. To change the visual selection, rerun only `plot_hall.py`. For new region shapes use a supported periodic-circle or k-ID JSON definition; K/K′ coordinates depend on the actual reciprocal basis. The [transport guide](KUBO_TRANSPORT.md) explains this scope and the difference between regional charge response and a valley-current operator.

## 4. Generate that reusable data for your material

First prepare a converged VASP calculation and a full uniform 2D SOC WAVECAR, following [input selection](../examples/APPLY_TO_YOUR_SYSTEM.md). A line path and a symmetry-reduced mesh are unsuitable for this BZ integral. The three stages below are general; the input path and `NX`, `NY` must describe your actual calculation.

```bash
wavecar=/absolute/path/to/completed-full-mesh/WAVECAR
NX=24
NY=24
mkdir -p results/my-native
(
  cd results/my-native
  mpiexec -n 4 "$vb_bin" --task kubo-pairs --wavecar "$wavecar" \
    --spinor 2 --pairs-csv PAIRS.csv > vaspberry.log 2> vaspberry.err
)

python3 tools/vaspberry_kubo.py import-pairs \
  --csv results/my-native/PAIRS.csv --wavecar "$wavecar" \
  --spinor-components 2 --spin-multiplicity 1 --mesh "$NX" "$NY" \
  --energy-reference 'unchanged VASP eigenvalue zero' \
  --output-dir results/my-pairs
```

`PAIRS.csv` has one row for each k/spin and unordered band pair n<m. It contains the two energies, their gap, k coordinates and `numerator_yz_eV2_A2`, `numerator_zx_eV2_A2`, `numerator_xy_eV2_A2`. These are **Nab = −2 Im(Da,nm Db,mn)** in eV² Å², not yet curvature or conductivity. Native pair export does not divide by the squared energy gap or apply occupations.

`import-pairs` checks that every row matches WAVECAR and that the declared mesh is complete. It writes the same numerical information as compact arrays in `pairs.npz`, with checked geometry, units and provenance in `pairs.json`. It does **not** calculate Hall conductivity. Keep both cache files together.

Then use the `pair-hall` and plot commands from step 3 with `--pairs-dir results/my-pairs` and **your** μ range/reference, regions and retained-band cutoff. Pair export uses every stored band and k point, so it needs neither `--bands` nor `--mesh`. `--pair-band-max` selects the retained virtual-state window during integration. The [complete MoS₂ tutorial](../examples/features/kubo-hall/) provides a specified full-mesh VASP input, density, reference energies and convergence cases.

### Optional one-command export and integration

After preparing the tutorial's actual 24×24, 60-band MoS₂ WAVECAR,
`wavecar-hall` can run the same native MPI export, cache validation and Hall
integration together:

```bash
python3 tools/vaspberry_kubo.py wavecar-hall \
  --wavecar results/mos2-24-b60-vasp/WAVECAR --binary "$vb_bin" --mpi-procs 4 \
  --spinor-components 2 --spin-multiplicity 1 --mesh 24 24 \
  --energy-reference 'unchanged VASP eigenvalue zero' --pair-band-max 40 \
  --mu-min -1.47487388 --mu-max -1.17487388 --mu-num 61 \
  --mu-reference -0.43809870 --temperatures 0 300 \
  --regions examples/features/kubo-hall/regions.json --difference valley:K:Kprime \
  --degeneracy-policy coalesce --formats csv dat npz \
  --output-dir results/mos2-hall-workflow
```

The new directory contains `native/PAIRS.csv` and native logs, `pairs/` caches,
`hall/conductivity.*` and `workflow.json`. This is an execution shortcut;
the Fortran calculation and Python integration retain their distinct roles.
Plot `results/mos2-hall-workflow/hall/conductivity.csv` with the command in
step 3. Only the native stage uses MPI; the Python integration is a separate
NumPy calculation. Reuse `pairs/` for subsequent μ/T scans.

## 5. Add layer, orbital and spin character

Use matching SOC `PROCAR`, `WAVECAR` and `OUTCAR` from one final VASP run. The public [PROCAR example](../examples/features/procar-character/) gives complete commands and a small runnable analytic fixture. Its three stages have distinct jobs:

| Command | Inputs → outputs | Result or figure |
|---|---|---|
| `procar_character.py project` | VASP files + atom/orbital groups + spin axis → `character.csv/npz/json`, `projection_diagnostics.csv` | Raw group charge, Pauli weights and joint spin projections of each state |
| `procar_character.py hall` | character cache + native pair cache + isolated bands + μ/T/regions → `character_hall.csv/npz/json`, `selected_curvature.npz` | Calculate selected-band, character-weighted charge-Hall contributions and their reference changes |
| `procar_character.py plot` | saved character/Hall files + selected group/band/T → `character.*`, `character_hall.*` | k-space character maps and Hall curves in PNG/PDF/SVG |

The atom IDs and Cartesian spin axis are user inputs; no material-specific layer names are built into the tool. The `hall` stage is a numerical attribution calculation. The `plot` stage only visualizes its results. These charge contributions explain state character; conventional spin-current response is a separate observable.

## 6. Choose what can be inferred from the outputs

- Native Fukui/Z₂ use wavefunction overlaps. Kubo point curvature and occupation-weighted transport have different discretizations and convergence checks.
- Native Kubo uses pseudo-wavefunction canonical momentum. It is a useful baseline for charge and regional analysis; full material velocity may include missing PAW, nonlocal, SOC and Hubbard-U terms.
- Matching SOC `PROCAR`, `WAVECAR` and `OUTCAR` support atom/layer/orbital and chosen-axis spin character. The [PROCAR hands-on example](../examples/features/procar-character/) uses `tools/procar_character.py project`, `hall` and `plot` to save those weights, attribute selected isolated bands' charge Hall response, and make reusable figures. Define layer groups by actual atom IDs. This is character-weighted charge attribution; physical orbital-, layer- and spin-current operators are separate observables.
- The optional [spin Hall route](SPIN_HALL.md) needs full spin/velocity matrices and supports a gapped 2D occupied group at T=0. The supplied instrumented VASP producer excludes Hubbard U; do not apply its Bi instructions unchanged to a DFT+U magnetic material.
- Converge k sampling, source `NBANDS`, retained pair window and region choice separately. A smooth μ curve or picture does not refine the original k mesh.

For report figures, use the linked feature/material recipes and their compact reference outputs. Those example-specific assembly scripts arrange known panels; general native commands and `tools/` interfaces accept user data within their documented contracts. Optional Wannier examples are collected in the report's supporting appendix.
