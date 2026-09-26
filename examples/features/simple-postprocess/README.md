# One settings file: VASPBERRY Bi pairs, Hall table and figures

This example uses the public **Bi bilayer SOC WAVECAR on a complete 12×12
mesh**. One INI file specifies the VASPBERRY executable, input, output and Hall
scan. `run` launches VASPBERRY to calculate interband pairs, then performs
numerical postprocessing to integrate them. `plot` draws the saved table.

At T=0, the supplied μ range lies inside the sampled gap between occupied
bands 1–10 and empty bands 11–18. The ideal time-reversal-symmetric total
charge-Hall response is zero; this fixed input and finite pair sum give a
small residual, documented below. This is an installation/workflow
demonstration with a real VASP input, not a convergence study.
Its Kubo pair integral is distinct from the occupied-subspace Fukui result
shown in the [insulating Hall example](../hall-valley/README.md).

## Choose how far to go

Start with steps 1–3 for the first figure. The remaining steps each change
one thing; no group, region or plot section is required for the first run.

| Input file | What changes | Result directory |
|---|---|---|
| [`bi.ini`](bi.ini) | First run: VASPBERRY pair export and a total Hall scan | `results/simple-bi` |
| [`bi-rescan.ini`](bi-rescan.ini) | A different μ scan, reusing the first run | `results/simple-bi-rescan` |
| [`bi-regions.ini`](bi-regions.ini) | A named subset of two k points, reusing the first run | `results/simple-bi-regions` |

Read the [beginner guide](../../../docs/POSTPROCESSING.md) when adapting this
example to your material. Keep the [settings reference](../../../docs/POSTPROCESSING_REFERENCE.md)
for looking up a key rather than reading all options before starting.

## 1. Build and obtain the input

Run from the repository root. Load Intel Fortran, Intel MPI and oneMKL using
your actual installation or site modules; a common Linux setup is:

```bash
source /opt/intel/oneapi/setvars.sh
make ifx-mpi
make check-ifx-mpi
python3 examples/fetch_inputs.py bi --output-dir results/simple-bi-input
```

The fetcher verifies the public 200,421,600-byte WAVECAR against its stored
SHA-256. `results/simple-bi-input` must be new. The VASP input is already
calculated, so this example needs no VASP executable or licensed POTCAR.
Use your usual Python environment; plotting uses Matplotlib. Compiler and
MPI alternatives are in the [build guide](../../../docs/BUILD.md).

## 2. Check and calculate

Inspect [`bi.ini`](bi.ini). Its paths are relative to the INI directory;
the `../../../` prefixes point back to the repository root. The calculation
settings are:

| Setting | Value for this input |
|---|---|
| Source | SOC spinors; complete 12×12 mesh; 18 stored bands |
| μ scan | −1.3 to −0.9 eV, 9 points |
| Reference μ | −1.1 eV |
| Temperature | 0 K |
| VASPBERRY execution | Intel MPI, 4 ranks |

```bash
python3 tools/vaspberry_post.py check examples/features/simple-postprocess/bi.ini
python3 tools/vaspberry_post.py run examples/features/simple-postprocess/bi.ini
```

`run` executes the compiled VASPBERRY program through MPI, then runs the
Python Hall integration. Only VASPBERRY runs under MPI. Do not prepend
`mpiexec` to the Python command. Choose a rank count allowed by your machine
or scheduler allocation. For GNU MPI, build with `make mpi`, then change
`binary` to `../../../build/vaspberry-mpi` and `mpi_launcher` to your matching
GNU MPI launcher in the same INI. For GNU serial, use `make serial`,
`binary = ../../../build/vaspberry`, and `mpi_procs = 1`.

<details>
<summary>Optional: run VASPBERRY directly for the raw pair data</summary>

For the full Hall table and figures, continue with the INI commands above and
below. You can also run VASPBERRY directly to obtain the raw pair data,
without the Python front end. After step 1, the equivalent pair export is:

```bash
mkdir results/direct-bi
mpiexec.hydra -n 4 build/vaspberry-ifx-mpi \
  --task kubo-pairs --wavecar results/simple-bi-input/WAVECAR \
  --spinor 2 --pairs-csv results/direct-bi/PAIRS.csv
```

`results/direct-bi` must be new. VASPBERRY reads `WAVECAR` and writes energies,
k coordinates and interband matrix-element numerators to `PAIRS.csv`. That
file is the input for numerical Hall integration, not a conductivity table.
This direct command does not create an INI-run result directory, so do not
pass `results/direct-bi` to the front end's `plot` or `--reuse` options. The
[manual Hall workflow](../kubo-hall/README.md) shows the individual import and
integration commands when you need them.

</details>

## 3. Draw the saved data

```bash
python3 tools/vaspberry_post.py plot results/simple-bi
```

Outputs under `results/simple-bi/` are:

| File | Meaning |
|---|---|
| `native/PAIRS.csv` | Written by VASPBERRY: 22,032 unordered band-pair rows (144 k points × 18×17/2 pairs). Each stores energies, k coordinates and three interband numerators in eV² Å². |
| `pairs/pairs.npz`, `pairs/pairs.json` | Checked reusable pair arrays and their units, lattice and source identity. |
| `hall/conductivity.csv`, `.dat`, `.npz` | Completed charge-Hall scan, with absolute σ, Δσ relative to −1.1 eV and represented carrier counts. There are 18 rows: 9 μ values for each of `total` and `rest`. With no explicit regions, both cover the full mesh. |
| `figures/charge-hall/hall.png`, `.pdf`, `.svg` | Total sheet σ in e²/h versus μ−reference, at T=0. |
| `settings.ini`, `run.json`, `logs/` | Saved settings, actual stage commands/status and logs. |

The JSON sidecars are generated records. The editable user input is the
INI file. The CSV/DAT tables can also be plotted using your own tools.
`plot` does not repeat the wavefunction calculation or the Hall integral.

## 4. Reuse the pairs for a second scan

[`bi-rescan.ini`](bi-rescan.ini) keeps the same source and changes the scan
to 13 points between −1.25 and −0.95 eV, with a new output directory:

```bash
python3 tools/vaspberry_post.py check examples/features/simple-postprocess/bi-rescan.ini --reuse results/simple-bi
python3 tools/vaspberry_post.py run examples/features/simple-postprocess/bi-rescan.ini --reuse results/simple-bi
python3 tools/vaspberry_post.py plot results/simple-bi-rescan
```

This skips VASPBERRY execution, copies the validated pair cache, and
calculates the new occupation-weighted table. Keep the same WAVECAR available
for the source-identity check. The VASPBERRY executable and MPI launcher are
not used with `--reuse`.
For a later fresh calculation, set them to match the compiler/MPI environment.

## 5. Name a k-space subset

[`bi-regions.ini`](bi-regions.ini) adds just these sections to the first
example and changes `[run] output` to `../../../results/simple-bi-regions`:

```ini
[region sampled_k]
k_ids = 1 2

[plot]
hall_regions = total sampled_k rest
```

`sampled_k` is an arbitrary label. `k_ids` selects the **1-based k-point
indices** in the WAVECAR/VASPBERRY output. Here it selects the first two of the
144 saved points to demonstrate the interface. It is not a definition of a
physical valley or M point. A symmetry-point label such as `M` also needs
actual coordinates and a radius, as explained in the
[region guide](../../../docs/POSTPROCESSING.md#add-a-k-space-region).

```bash
python3 tools/vaspberry_post.py check examples/features/simple-postprocess/bi-regions.ini --reuse results/simple-bi
python3 tools/vaspberry_post.py run examples/features/simple-postprocess/bi-regions.ini --reuse results/simple-bi
python3 tools/vaspberry_post.py plot results/simple-bi-regions
```

This does not rerun VASPBERRY. `results/simple-bi-regions/hall/conductivity.csv`
has 27 rows: 9 μ values for each of `total`, `sampled_k` and `rest`. Each subset
keeps its original full-BZ integration weights; `rest` is the complement of
the selected points, so `sampled_k + rest = total` at each μ and temperature.
The three curves are drawn together in
`results/simple-bi-regions/figures/charge-hall/hall.png` (also PDF/SVG).
They illustrate partitioning; their magnitudes are not a valley-convergence result.

## 6. Redraw a completed table

For a new plot of Δσ rather than absolute σ, use:

```bash
python3 tools/vaspberry_post.py plot results/simple-bi --temperature 0 --quantity delta-sigma --output-dir results/simple-bi-redraw
```

This reads the existing Hall table without using WAVECAR, starting MPI or
repeating integration. All numerical and figure output directories must be
new. Plot choices stored in an edited INI do not modify an earlier run;
`plot` reads that run's saved settings, with supported command-line overrides.

## Scope and the next calculation

The sampled gap is −1.348963526 to −0.838919163 eV, in the unchanged VASP
energy zero. Both scans stay inside it with ten occupied states at every
k point. The gap and sampling come from this fixed public input; they do
not establish material convergence. The reference run (GNU Fortran/OpenMPI,
two ranks) gives σ ≈ −6.27668023×10⁻⁶ e²/h at every μ, Δσ = 0 and ten
represented electrons per cell. The rescan reproduces the same values at its
13 μ points. This residual is a regression reference for the fixed input,
operator and stored-band sum; its cause and convergence are not established
by this example. It should not be interpreted as an anomalous Hall signal.

The VASPBERRY pair export uses canonical momentum of the stored
pseudo-wavefunctions.
The input includes near-degenerate Kramers partners, so this example uses the
total occupied response and does not request individual-band PROCAR attribution.
The near-degenerate occupied pairs have equal T=0 occupations and cancel
before division in the pair integral.

For your own material, replace the input, mesh, spin mode, energy range and
reference. A full periodic mesh is required; the supplied MoS₂ band-path
WAVECAR is not a Hall-integration input. The
[beginner guide](../../../docs/POSTPROCESSING.md#use-your-own-vasp-calculation)
shows which settings to change first. Layer/orbital/spin character additionally
needs matching SOC `PROCAR` and `OUTCAR` files; they are not supplied as a
matched projection dataset for this Bi walkthrough. Follow the
[PROCAR tutorial](../procar-character/) for those inputs and its separate,
explicitly analytic fixture. The [technical report](../../../docs/TECHNICAL_REPORT.md)
and [report reproduction table](../../REPORT_REPRODUCTION.md) connect the
feature-specific physical examples to reference figures.
