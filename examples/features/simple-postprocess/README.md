# One settings file: native Bi pairs, Hall table and figures

This example uses the public **Bi bilayer SOC WAVECAR on a complete 12×12
mesh**. One INI file specifies the native executable, input, output and Hall
scan. Fortran calculates interband pairs; numerical postprocessing integrates
them; a separate command draws the saved table.

At T=0, the supplied μ range lies inside the sampled gap between occupied
bands 1–10 and empty bands 11–18. The ideal time-reversal-symmetric total
charge-Hall response is zero; this fixed input and finite pair sum give a
small residual, documented below. This is an installation/workflow
demonstration with a real VASP input, not a convergence study.
Its Kubo pair integral is distinct from the occupied-subspace Fukui result
shown in the [insulating Hall example](../hall-valley/README.md).

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
| Native execution | Intel MPI, 4 ranks |

```bash
python3 tools/vaspberry_post.py check examples/features/simple-postprocess/bi.ini
python3 tools/vaspberry_post.py run examples/features/simple-postprocess/bi.ini
```

The front end launches only the Fortran stage with MPI. Do not prepend
`mpiexec` to the Python command. Choose a rank count allowed by your machine
or scheduler allocation. For GNU MPI, build with `make mpi`, then change
`binary` to `../../../build/vaspberry-mpi` and `mpi_launcher` to your matching
GNU MPI launcher in the same INI. For GNU serial, use `make serial`,
`binary = ../../../build/vaspberry`, and `mpi_procs = 1`.

The native stage is equivalent to running `--task kubo-pairs` on the named
WAVECAR. You do not need to assemble separate import and integration
commands or write a JSON input file.

## 3. Draw the saved data

```bash
python3 tools/vaspberry_post.py plot results/simple-bi
```

Outputs under `results/simple-bi/` are:

| File | Meaning |
|---|---|
| `native/PAIRS.csv` | 22,032 unordered band-pair rows: 144 k points × 18×17/2 pairs. Each stores energies, k coordinates and three interband numerators in eV² Å². |
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

This skips native Fortran export, copies the validated pair cache, and
calculates the new occupation-weighted table. Keep the same WAVECAR available
for the source-identity check. If you selected GNU for the first run, make
the corresponding executable settings consistent in `bi-rescan.ini` too.
For another plot destination, use:

```bash
python3 tools/vaspberry_post.py plot results/simple-bi --temperature 0 --quantity delta-sigma --output-dir results/simple-bi-redraw
```

This redraws Δσ from the existing table without using WAVECAR or starting
MPI. All numerical and figure output directories must be new.

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

The native operator is canonical momentum of the stored pseudo-wavefunctions.
The input includes near-degenerate Kramers partners, so this example uses the
total occupied response and does not request individual-band PROCAR attribution.
The near-degenerate occupied pairs have equal T=0 occupations and cancel
before division in the pair integral.

For your own material, replace the input, mesh, spin mode, energy range and
reference. A full periodic mesh is required; the supplied MoS₂ band-path
WAVECAR is not a Hall-integration input. The
[settings guide](../../../docs/POSTPROCESSING.md) shows optional atom/layer/
orbital/spin groups and k-space regions in the same INI, explains every
output stage, and separates native calculation, occupation integration and
plotting.
