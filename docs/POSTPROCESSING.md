# Postprocessing from one settings file

Keep input paths, the chemical-potential range, temperatures and optional
atom groups in one commented INI file. VASPBERRY still performs the
wavefunction matrix calculation in Fortran; the Python front end runs that
executable and the numerical postprocessing stages in order.

```bash
python3 tools/vaspberry_post.py check analysis.ini
python3 tools/vaspberry_post.py run analysis.ini
python3 tools/vaspberry_post.py plot results/run01
```

`check` checks the settings and source files without starting a calculation.
`run` exports native pairs, creates a reusable cache, and calculates Hall
tables. `plot` reads completed tables and draws figures. Numerical validity
checks also run during the calculation; passing `check` alone is not a
convergence result.

For a complete public example, start with the
[Bi settings and commands](../examples/features/simple-postprocess/README.md).
It uses an actual full-mesh VASP WAVECAR and no manually written JSON.

## Set up the native executable once

Load your site's Intel Fortran, Intel MPI and oneMKL environment. For a common
Linux installation, the commands are:

```bash
source /opt/intel/oneapi/setvars.sh
make ifx-mpi
make check-ifx-mpi
```

Use the actual installation path or equivalent cluster modules. The
[build guide](BUILD.md) also covers retained Intel Classic and GNU builds.
Run the Python command once, without `mpiexec` in front of it. The
`mpi_procs` setting launches the **Fortran executable** on that many ranks
using the matching MPI launcher. Run within an allocation permitting those
ranks. Keep the compiler/library environment loaded.

Use your usual Python environment. The supplied plotting commands use
Matplotlib; the native Fortran calculation does not require Python.

## Describe your calculation

Put `analysis.ini` beside your own full-mesh `WAVECAR`. Here is the structure
of a minimal SOC study; replace the executable path, mesh and energy values
with those for your calculation:

```ini
[run]
wavecar = WAVECAR
binary = /path/to/VASPBERRY/build/vaspberry-ifx-mpi
output = results/run01
mesh = 24 24
spin_mode = soc
energy_reference = unchanged VASP eigenvalue zero
mpi_procs = 4
mpi_launcher = mpiexec.hydra

[hall]
# Chemical potentials in eV: MIN MAX NUMBER_OF_POINTS.
mu = -0.5 0.5 101
reference = 0.0
# Temperatures in kelvin, separated by spaces.
temperatures = 0 100 300
```

Paths are relative to the **INI file's directory**. Absolute paths and `~`
are accepted; shell variables and shell commands are not expanded. If you
move a settings file, update its relative paths. Section names and keys are
case-sensitive; lists use spaces rather than commas. Use `#` or `;` for
comments. Each `output` must name a new directory.

`mesh` describes all points stored in WAVECAR. A symmetry-reduced mesh or a
band path cannot replace the complete uniform two-dimensional mesh needed
for Hall integration. `spin_mode = soc` selects two-component spinors and
counts each stored spinor state once. The other explicit choices are
`scalar-degenerate`, `collinear-up`, and `collinear-down`; select the mode
matching the actual VASP input. A collinear channel is one channel's
response, not the sum of both spin channels.

`mu` and `reference` use the **same energy zero as the WAVECAR eigenvalues**.
The `energy_reference` text records that convention; it does not shift the
energies. `reference` defines Δσ(μ) = σ(μ) − σ(reference). The figures use
μ−reference on their horizontal axis. Choose these values from your own
band energies rather than copying a different material's example.

All stored bands enter the pair sum by default. If a convergence study
requires a smaller retained window, add `pair_band_max = 40` under `[hall]`,
for example. This retains bands 1–40 in the total charge-Hall integration;
it does not change what native Fortran exports. The source needs at least
that many stored bands. The retained window and number of stored bands
remain distinct quantities to converge.

## What is saved

With `output = results/run01`, the result directory contains:

| File | Contents and use |
|---|---|
| `native/PAIRS.csv` | Raw Fortran output: each k point/spin and unordered pair n<m, the two energies, fractional k coordinates, gap, and three antisymmetric matrix products in eV² Å². These numerators have no occupation factors or energy denominators yet. |
| `native/stdout.log`, `native/stderr.log` | The native calculation's messages. |
| `pairs/pairs.npz`, `pairs/pairs.json` | Validated pair arrays and their geometry, units and source identity; keep both files together for rescans. |
| `hall/conductivity.csv`, `.dat`, `.npz` | The completed numerical table: μ, T, region, sheet σ in e²/h and siemens, reference-subtracted Δσ in e²/h, and represented electron/carrier counts. |
| `hall/conductivity.json` | Units, operator, band-window and integration diagnostics for that table. |
| `settings.ini`, `run.json`, `logs/` | A settings snapshot, resolved paths/source hashes/stage status, and postprocessing logs. Edit your original INI for the next study. |
| `groups.json`, `regions.json` | Automatically generated definitions from optional sections of your INI. You do not need to create or edit these JSON files. |

The import stage checks and stores the native values. The Hall stage performs
a calculation: it divides by squared energy differences, applies Fermi
occupations, and integrates over the Brillouin zone or selected regions.
It does not rerun VASP or change the self-consistent electronic structure.

After `plot`, the default `figures/charge-hall/` directory contains
`hall.png`, `hall.pdf`, and `hall.svg`. These plot sheet σ or Δσ against
chemical potential at the saved temperatures. You can also load the CSV in
Origin, gnuplot, a spreadsheet, or your own Python plotting program. The
[output specification](OUTPUT_FORMAT.md) defines individual columns and
conventions.

For a directly calculated curvature map or band-path curve, native
`--task kubo` already writes `KUBO.csv` in Å². Those commands and plotting
examples are in the [hands-on guide](HANDS_ON.md); the Hall INI workflow
above uses the more general pair output to scan occupations.

## Change the scan without repeating Fortran

Keep the same WAVECAR and source settings. Change `[hall]` and the output
name in your INI, then point `--reuse` to the completed previous run:

```bash
python3 tools/vaspberry_post.py check analysis-next.ini --reuse results/run01
python3 tools/vaspberry_post.py run analysis-next.ini --reuse results/run01
python3 tools/vaspberry_post.py plot results/run02
```

Here `analysis-next.ini` must specify `output = .../run02` using its own
relative path convention. `--reuse` reads a completed result from this
front end, verifies the same WAVECAR and mesh/spin/energy conventions, and
copies the validated pair cache into the new result. It skips native pair
export and recalculates the requested occupation integral. Keep the WAVECAR
available for this identity check. A changed Hamiltonian, geometry or VASP
calculation requires a new native run.

To draw one saved temperature and the reference-subtracted response into
another folder:

```bash
python3 tools/vaspberry_post.py plot results/run01 --temperature 300 --quantity delta-sigma --output-dir results/run01-figures
```

The plot command does not require the original WAVECAR, native executable
or MPI launcher. It uses the completed tables and, by default, the plot
choices saved when that run was created. `--temperature` selects a temperature
already calculated; `--quantity` chooses `sigma` or `delta-sigma`. Editing
the original INI does not change an existing result. The result and figure
directories are preserved; choose new names when rerunning.

## Optional atom, layer, orbital and spin character

For your own SOC calculation, keep matching `PROCAR`, `OUTCAR` and
`WAVECAR` from the same final static calculation. Add `[projection]` and
named groups to the same INI. The following illustrates the syntax for
atoms 1 and 2 and isolated band 1; replace the IDs with those of your actual
structure and the bands you intend to analyze:

```ini
[projection]
bands = 1
axis = 0 0 1

[group lower]
ions = 1

[group upper]
ions = 2
```

`PROCAR` and `OUTCAR` default to the WAVECAR directory. Set `procar` and
`outcar` under `[projection]` if their paths differ. `axis` is a Cartesian
unit vector. Atom and band IDs are 1-based; layer names are your definitions
from the actual structure. Add an optional `orbitals = dxy dyz dz2 dxz x2-y2`
line to a group only when those labels occur in its PROCAR header. Projection
uses the Hall scan by default; its own `mu`, `reference`, and `temperatures`
can be specified under `[projection]`. Its virtual-state sum uses all stored
pair bands, regardless of `[hall] pair_band_max`; that setting limits only
the total charge-Hall table.

This adds `character/character.csv` and `.npz` for raw charge/spin character,
and `character-hall/character_hall.csv` and `.npz` for the selected-band,
character-weighted charge-Hall contributions. Projection diagnostics and
JSON metadata accompany both. The joint spin weights use
`plus=(q+m_axis)/2` and `minus=(q-m_axis)/2`; the raw local projections are
not normalized to sum to one.

With projection enabled, `plot` also creates state-character maps and
selected-group Hall curves under `figures/character/`. Set default plot
choices before `run`, for example:

```ini
[plot]
character_group = lower
map_band = 1
character_temperature = 0
character_region = total
character_delta = false
```

You can redraw another saved group's character and Hall curves without
reintegrating:

```bash
python3 tools/vaspberry_post.py plot results/run01 --group upper --band 1 --temperature 0 --output-dir results/run01-upper
```

`--band` selects the state-character map; it does not change the bands
included in the saved projected Hall sum. Changing that sum requires a new
`run` with `--reuse` and a different `[projection] bands` selection.

Selected bands must be isolated from every other stored band; unresolved
degeneracies are rejected. The public Bi example has Kramers-degenerate
states and deliberately omits this band-resolved projection. These are
character-weighted **charge-Hall contributions**; independent spin-, layer-
and orbital-current responses are different observables. The
[projection guide](../examples/features/procar-character/README.md) supplies
an explicitly analytic fixture and the detailed matching checks.

## Optional regional curves in the same file

Define periodic disks with fractional reciprocal coordinates and radii in
Å⁻¹. This example is the **MoS₂ reciprocal basis** used by the
[regional Hall tutorial](../examples/features/kubo-hall/README.md); it is not
a transferable valley definition for the Bi demonstration or another cell:

```ini
[region K]
center = 1/3 2/3 0
radius = 0.35

[region Kprime]
center = -1/3 -2/3 0
radius = 0.35

[differences]
valley = K Kprime

[plot]
hall_regions = total K Kprime valley
hall_temperatures = 300
hall_quantity = delta-sigma
```

Include 300 K in your `[hall] temperatures` for this plot selection.
If your INI already has a `[plot]` section, merge these keys into that section;
do not add a second `[plot]` header.
Alternatively, a region can use `k_ids = 1 2 3` to select known 1-based
k-point IDs. Regions must form a valid disjoint partition; `rest` covers the
remaining points. The `valley` difference above is K−K′ without a factor of
one half, and is a regional charge response defined by your partition.
Differences apply to the total Hall table; projected Hall keeps the named
regions themselves.

This compact front end uses the existing strict numerical defaults. For an
explicit numerical-degeneracy treatment or other advanced controls, use the
individual [Kubo commands](KUBO_TRANSPORT.md). Source-band, k-mesh and
operator accuracy still need to be assessed for the material under study.
