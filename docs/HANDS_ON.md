# Hands-on: native calculation, saved data, and reusable plots

The main workflow is **VASP → WAVECAR → VASPBERRY Fortran → numerical files → analysis and plots**. Charge-Hall transport adds an occupation-weighted Kubo integration after the Fortran matrix-element calculation. Wannier is an optional supporting route.

Use VASPBERRY 1.4.0 or later for the commands below. Start at the repository root with GNU Fortran, BLAS/LAPACK and Python 3.10+ installed. Commands create new outputs; choose a new result directory when repeating a calculation.

```bash
make serial
python3 -m pip install -r requirements-transport.txt
repo_dir="$PWD"
./build/vaspberry --help
```

[Native command reference](NATIVE_COMMANDS.md) explains tasks, arguments and filenames. [Report-to-example map](../examples/REPORT_REPRODUCTION.md) links material-specific inputs and expected figures.

## 1. Calculate a real supplied WAVECAR with Fortran

This first calculation uses the ordinary SOC MoS₂ band-path WAVECAR already in the repository. No VASP execution or Wannier input is needed. It evaluates occupied-bundle point curvature, with bands 1–18 and the stored empty bands as intermediate states.

```bash
mkdir -p results/first-kubo
(
  cd results/first-kubo
  "$repo_dir/build/vaspberry" --task kubo \
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

Keep the native CSV and log. Native legacy side products also stay in this result directory. Inspect the header's normalization, operator and gap information. The bundle must be isolated from excluded states. This input is a line path: it cannot supply a full-zone Hall integral or reproduce the report's separately generated full-mesh figure by itself.

## 2. Plot the saved CSV with your own tool

The following is ordinary Python/Matplotlib; it does not invoke the Fortran code or recalculate a matrix element. The example plots the ordered samples so no lattice geometry or symmetry labels are invented.

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

`pair-hall` calculates Fermi occupations and the BZ integral. `plot_hall.py` reads the finished table and writes PNG/PDF/SVG. `conductivity.csv`, `.dat`, `.npz` and `.json` retain absolute σ, Δσ relative to μref, carrier counts, regions, units and provenance. The explicit `coalesce` policy approximates tiny numerical energy splittings by a group mean and records the shifts; it is not a universal way to remove physical near-crossings.

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
  "$repo_dir/build/vaspberry" --task kubo-pairs --wavecar "$wavecar" \
    --spinor 2 --pairs-csv PAIRS.csv > vaspberry.log 2> vaspberry.err
)

python3 tools/vaspberry_kubo.py import-pairs \
  --csv results/my-native/PAIRS.csv --wavecar "$wavecar" \
  --spinor-components 2 --spin-multiplicity 1 --mesh "$NX" "$NY" \
  --energy-reference 'unchanged VASP eigenvalue zero' \
  --output-dir results/my-pairs
```

Then use the `pair-hall` and plot commands from step 3 with `--pairs-dir results/my-pairs` and **your** μ range/reference, regions and retained-band cutoff. The numeric MoS₂ choices in step 3 are not material-independent defaults. Pair export uses every stored band and all stored k points, so it needs neither `--bands` nor `--mesh`. The importer checks mesh coverage against WAVECAR; integration selects the virtual-state cutoff with `--pair-band-max`.

For MPI, run `make mpi` and replace the native executable by `mpiexec -n 4 "$repo_dir/build/vaspberry-mpi"`. Keep the remaining arguments. Start with the [complete MoS₂ VASP → pairs → Hall tutorial](../examples/features/kubo-hall/) for a specified material, including its PAW datasets, SCF density, preparation commands and reference energies.

## 5. Choose what can be inferred from the outputs

- Native Fukui/Z₂ use wavefunction overlaps. Kubo point curvature and occupation-weighted transport have different discretizations and convergence checks.
- Native Kubo uses pseudo-wavefunction canonical momentum. It is a useful baseline for charge and regional analysis; full material velocity may include missing PAW, nonlocal, SOC and Hubbard-U terms.
- Matching SOC `PROCAR`, `WAVECAR` and `OUTCAR` support atom/layer/orbital and chosen-axis spin character. The [PROCAR hands-on example](../examples/features/procar-character/) uses `tools/procar_character.py project`, `hall` and `plot` to save those weights, attribute selected isolated bands' charge Hall response, and make reusable figures. Define layer groups by actual atom IDs. This is character-weighted charge attribution; physical orbital-, layer- and spin-current operators are separate observables.
- The optional [spin Hall route](SPIN_HALL.md) needs full spin/velocity matrices and supports a gapped 2D occupied group at T=0. The supplied instrumented VASP producer excludes Hubbard U; do not apply its Bi instructions unchanged to a DFT+U magnetic material.
- Converge k sampling, source `NBANDS`, retained pair window and region choice separately. A smooth μ curve or picture does not refine the original k mesh.

For report figures, use the linked feature/material recipes and their compact reference outputs. Those example-specific assembly scripts arrange known panels; general native commands and `tools/` interfaces accept user data within their documented contracts. Optional Wannier examples are collected in the report's supporting appendix.
