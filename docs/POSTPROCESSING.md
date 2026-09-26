# Postprocessing: from a VASP file to a figure

Start with **one settings file and three commands**. Add layer/spin or
k-space regions only when your question needs them. You do not need to
write JSON files or assemble the individual Python stages.

| Command | What it does | What you get |
|---|---|---|
| `check analysis.ini` | Check settings, input files and mesh before running. | A summary or an actionable error; no calculation directory. |
| `run analysis.ini` | Execute VASPBERRY, with MPI when configured, then calculate the requested Hall tables. | Numerical data in the INI's `output` directory. |
| `plot results/run01` | Read a completed result and draw it. | PNG, PDF and SVG figures; no new wavefunction calculation. |

Prefix each command with `python3 tools/vaspberry_post.py`. `run` launches
the compiled **VASPBERRY executable** named by `[run] binary`,
then applies the requested numerical postprocessing to its output. `run`
also performs the checks, so the separate `check` is useful when preparing
an input. Numerical validity checks continue during `run`; this preliminary
check does not establish convergence.

**Choose a starting point:** [try the public example](#start-with-the-public-bi-example),
[execute VASPBERRY directly](#execute-vaspberry-directly-for-a-curvature-table),
[use your own data](#use-your-own-vasp-calculation),
[add layer/spin analysis](#add-layer-and-spin-character), or
[add a region](#add-a-k-space-region). Look up individual options in the
[settings reference](POSTPROCESSING_REFERENCE.md).

## Start with the public Bi example

Run commands from the repository directory containing `Makefile` and
`tools/`. Load Intel Fortran, Intel MPI and oneMKL with your site's modules
or, for a common Linux installation, the setup script below. Use your usual
Python environment; the supplied figures use Matplotlib.

```bash
source /opt/intel/oneapi/setvars.sh
make ifx-mpi
make check-ifx-mpi
python3 examples/fetch_inputs.py bi --output-dir results/simple-bi-input
python3 tools/vaspberry_post.py check examples/features/simple-postprocess/bi.ini
python3 tools/vaspberry_post.py run examples/features/simple-postprocess/bi.ini
python3 tools/vaspberry_post.py plot results/simple-bi
```

Replace the setup path with your installed oneAPI path or site modules.
The supplied [`bi.ini`](../examples/features/simple-postprocess/bi.ini) already
specifies this public input, the full 12×12 mesh, the energy scan and four
Intel MPI ranks. Use an allocation permitting those ranks. **Run Python
once**: the helper launches VASPBERRY through MPI with the configured
`mpi_procs` and `mpi_launcher`.
GNU/serial substitutions are in the [example's build instructions](../examples/features/simple-postprocess/README.md#2-check-and-calculate)
and the [build guide](BUILD.md). Use new result directories when repeating a run.

Open **`results/simple-bi/figures/charge-hall/hall.png`**. It shows sheet Hall
conductivity versus chemical potential relative to the reference. For this
particular insulating Bi input, expect a nearly flat, near-zero result
(about −6.28×10⁻⁶ e²/h). That is the documented finite-input residual, not
a missing plot or an established anomalous Hall signal. The
[complete example](../examples/features/simple-postprocess/README.md) lists
its reference values, inputs, outputs and limitations.

### What was calculated, and which file should I use?

```text
WAVECAR → VASPBERRY execution (MPI when configured) → native/PAIRS.csv
        → Python occupation/energy-denominator/BZ integration → hall/conductivity.csv
        → plotting → figures/charge-hall/hall.png (also PDF and SVG)
```

VASPBERRY calculates and exports the reusable interband matrix data.
Python then performs the numerical Hall integral. Plotting is a separate
operation on the completed table.
Neither step reruns VASP.

| Your purpose | File inside the result directory | Contents |
|---|---|---|
| Inspect the original VASPBERRY output | `native/PAIRS.csv` | k points, paired bands and energies, and interband numerators in eV² Å²; occupations and energy denominators have not yet been applied. |
| Plot or analyze Hall response in any program | `hall/conductivity.csv` or `.dat` | μ, temperature, region, sheet σ and reference-subtracted Δσ in e²/h, conductivity in siemens, and carrier counts. |
| Reuse the expensive matrix calculation | `pairs/pairs.npz` **and** `pairs/pairs.json` | Validated pair arrays and their source, units and geometry. Keep both. |
| Understand a failed stage or recover its command | `run.json`, `logs/`, `native/*.log` | Settings, source identity, stage commands, status and messages. |

Use Origin, gnuplot or your own plotting program with the tables if preferred.
The [output specification](OUTPUT_FORMAT.md) defines the columns. JSON files
in a result are generated records; your editable input is the INI.
A fresh `run` executes VASPBERRY and writes `native/PAIRS.csv`. A reused
run copies the pair cache and does not execute VASPBERRY again.

## Execute VASPBERRY directly for a curvature table

For a Berry-curvature curve, a direct VASPBERRY command is enough. Load
the Intel environment and compile with `make ifx-mpi` as described in the
[build guide](BUILD.md). Then run this from the repository root in an
allocation allowing four ranks:

```bash
repo_dir="$PWD"
mkdir -p results
mkdir results/mos2-direct-kubo && (
  cd results/mos2-direct-kubo
  mpiexec.hydra -n 4 "$repo_dir/build/vaspberry-ifx-mpi" \
    --task kubo \
    --wavecar "$repo_dir/examples/1H-MoS2/KPATH/2.band/WAVECAR" \
    --spinor 2 --bands 1:18 --bundle 1 \
    --curvature-csv KUBO.csv > vaspberry.log
)
```

This executes **VASPBERRY through MPI** without a Python helper. It reads
the supplied SOC MoS₂ WAVECAR and creates
`results/mos2-direct-kubo/KUBO.csv`, plus the execution log. The CSV has one
occupied-bundle curvature row for bands 1–18 at each of the 48 stored path
points: k-point index, fractional k coordinates, `omega_z_A2` in Å² and
the minimum energy gap to excluded bands in eV.

Import the CSV into Origin, gnuplot or another plotting program, skipping
lines beginning with `#`. Use `k_index` as x and `omega_z_A2` as y to inspect the
path curve. Python is not required for this calculation or for using its
CSV. The [supplied-path example](../examples/features/kubo-curvature/README.md#supplied-32-band-path-example)
explains this input and connects it to the full-mesh/map examples; the
[output specification](OUTPUT_FORMAT.md#native-kubo-bundle-csv) defines all columns.

A chemical-potential/temperature Hall scan requires additional occupation
weighting and integration over a complete 2D mesh; this path CSV does not
provide that integral. The INI `run` command above executes VASPBERRY's
`--task kubo-pairs` calculation and then performs those numerical steps.
If you prefer to execute each stage yourself, use the explicit VASPBERRY
and postprocessing commands in the [hands-on guide](HANDS_ON.md).

## Use your own VASP calculation

Create `analysis.ini` in the **repository root**. Use your own complete,
uniform two-dimensional mesh WAVECAR; a symmetry-reduced mesh or band path
cannot provide this Hall integral. Here is an illustrative SOC calculation:

```ini
[run]
wavecar = /path/to/your/vasp-calculation/WAVECAR
binary = build/vaspberry-ifx-mpi
output = results/run01
mesh = 24 24
spin_mode = soc
energy_reference = unchanged VASP eigenvalue zero
mpi_procs = 4
mpi_launcher = mpiexec.hydra

[hall]
# Example only: reference energy is 5.2 eV in this WAVECAR's energy zero.
# Scan from 0.2 eV below to 0.2 eV above it, at 101 points.
mu = 5.0 5.4 101
reference = 5.2
temperatures = 0 100 300
```

Set `wavecar` and `mesh` from your input; choose `mu`, `reference` and
`temperatures` for your research question. Set the executable/launcher once
for your installation. `spin_mode = soc` is for two-component spinors;
other spin conventions are in the [reference](POSTPROCESSING_REFERENCE.md).

**Energy numbers must use the WAVECAR eigenvalue zero.** In this example the
plot spans −0.2 to +0.2 eV because its x coordinate is `mu - reference`.
The tool does not read a Fermi energy automatically. The text
`energy_reference` documents the convention; it does not shift energies.
`reference` also defines Δσ(μ,T) = σ(μ,T) − σ(reference,T).

```bash
python3 tools/vaspberry_post.py check analysis.ini
python3 tools/vaspberry_post.py run analysis.ini
python3 tools/vaspberry_post.py plot results/run01
```

Paths **inside the INI** are relative to its directory. Paths on the command
line, such as `results/run01` or `--reuse results/run01`, are relative to your
shell's current directory. Moving an INI may require changing its paths.
The commands above assume it stays in the repository root.

For ordinary total Hall curves, **this is the complete input**. No group,
region, projection or plot section is required. All stored bands enter the
pair sum by default; mesh and band convergence remain part of the material
study. For point-curvature maps or symmetry-path curves calculated directly
by VASPBERRY, follow the [Kubo example](../examples/features/kubo-curvature/README.md).

## Add layer and spin character

Use this extension with your own SOC calculation's matching `WAVECAR`,
`PROCAR` and `OUTCAR`. PROCAR must contain the supported noncollinear
`LORBIT=11` charge/spin blocks. Add the following to the same INI **before
running**, replacing the illustrative atom and band IDs:

```ini
[projection]
bands = 1
axis = 0 0 1

[group layer1]
ions = 1

[group layer2]
ions = 2

[plot]
character_group = layer1
```

`layer1` and `layer2` are names you choose. **`ions` selects the actual
1-based POSCAR/PROCAR atom numbers**; for several atoms use, for example,
`ions = 1 3 5`. Naming a group does not discover a layer. Optional
`orbitals` selects actual PROCAR orbital labels. The axis is Cartesian +z
here; the OUTCAR supplies the source spin-frame rotation. PROCAR and OUTCAR
default to the WAVECAR directory.

`bands` selects the bands included in the character-weighted Hall sum.
Each selected band must be isolated from every other stored band; the public Bi input's
Kramers partners are unsuitable for this individual-band example. To learn
what the projection data mean without your own files, use the
[small analytic PROCAR example](../examples/features/procar-character/README.md).
It checks assigned analytic data; it is separate from the real Bi run.

Run the same `check`, `run`, `plot` commands. In addition to total Hall data,
you get `character/character.csv` (state/group charge and spin weights),
`character-hall/character_hall.csv` (selected-band, group/spin-weighted
charge-Hall contributions), and `figures/character/character.png` plus
`character_hall.png` (also PDF/SVG). These figures show where the states have
a given character and how that character contributes to charge Hall.
The [reference](POSTPROCESSING_REFERENCE.md) explains overlapping groups,
projection residuals and the distinction from a spin-current conductivity.

## Add a k-space region

A region selects k points for an additional Hall curve. Its name is also
yours to choose: `K`, `M`, `m` and `pocket1` are labels. The tool uses the
coordinates or IDs you supply, not the label, to choose points.

For example, **if the desired M point is (1/2, 0, 0) in your input's
reciprocal basis**, add:

```ini
[region M]
center = 1/2 0 0
radius = 0.10

[plot]
hall_regions = total M
```

`center` uses fractional reciprocal coordinates; `radius` is in Å⁻¹. This
selects a periodic disk **around** M. Verify the coordinate in your cell and
choose a radius containing sampled k points. Distinct regions must not
overlap. `total` covers the full mesh and `rest` covers points outside all
your regions; these two names are reserved.

If you already have `[plot]`, add the new key there; do not create a second
`[plot]` section. Names are case-sensitive and start with an English letter;
then letters, digits, `_`, `.` and `-` are allowed. When renaming a group or
region, update the corresponding plot and difference references too.

Try the [public Bi region exercise](../examples/features/simple-postprocess/README.md)
with [`bi-regions.ini`](../examples/features/simple-postprocess/bi-regions.ini).
It selects known k-point IDs and reuses the first run. The
[MoS₂ regional Hall tutorial](../examples/features/kubo-hall/README.md)
shows a physically specified K/K′ partition and its regional contrast.

## Change a calculation or redraw a figure

| What changed? | What to do |
|---|---|
| Only the figure: saved temperature, σ versus Δσ, or saved group/band to display | Use `plot` with an override and a new figure directory. |
| μ range, temperature grid, region, or group definition/name | Edit the INI, choose a new `output`, then `run ... --reuse PREVIOUS_RESULT`; plot that new result. |
| VASP electronic structure, geometry or WAVECAR | Execute VASPBERRY again with a fresh `run` and new output directory. |

For a new numerical scan, copy `analysis.ini` to `analysis-next.ini` in the
same directory, change the desired settings and set `output = results/run02`:

```bash
python3 tools/vaspberry_post.py check analysis-next.ini --reuse results/run01
python3 tools/vaspberry_post.py run analysis-next.ini --reuse results/run01
python3 tools/vaspberry_post.py plot results/run02
```

`--reuse` skips VASPBERRY execution, verifies the same WAVECAR/source
conventions and recalculates the postprocessing. Keep WAVECAR available;
projection also needs matching PROCAR/OUTCAR. Follow the
[ready-to-run Bi rescan](../examples/features/simple-postprocess/README.md#4-reuse-the-pairs-for-a-second-scan)
for a complete example.

To draw an already calculated temperature and Δσ in a different folder:

```bash
python3 tools/vaspberry_post.py plot results/run01 --temperature 300 --quantity delta-sigma --output-dir results/run01-redraw
```

For an already calculated projection group, add `--group layer2` to a plot
command. **Plotting uses the saved result's settings**, so editing your
original INI does not change an existing result. `plot` needs only the
completed result files; it requires no VASP inputs, VASPBERRY executable
or MPI launcher. Preserve the saved result and choose a new figure
directory for another redraw.

## Find a setting or solve a problem

- [Settings reference](POSTPROCESSING_REFERENCE.md): every INI key, defaults,
  units, naming rules, input requirements and common errors.
- [Public Bi example](../examples/features/simple-postprocess/README.md): first
  run, expected table/figure, rescan and region exercise.
- [PROCAR example](../examples/features/procar-character/README.md): analytic
  projection check and applying atom/orbital/spin groups to your own data.
- [Output formats](OUTPUT_FORMAT.md): exact columns for your plotting tools.
- [Technical report](TECHNICAL_REPORT.md): scientific definitions and scope;
  [example/report map](../examples/REPORT_REPRODUCTION.md) connects its figures
  to the relevant example and data.

Keep explicit band/mesh/operator convergence choices when adapting the
example. Advanced degeneracy controls and individual pipeline stages remain
in the [Kubo reference](KUBO_TRANSPORT.md).
