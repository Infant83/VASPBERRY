# VASPBERRY command reference

**1.6.6 input policy:** the standard charge-Kubo route uses WAVEDER. Examples
that retain the WAVECAR canonical-momentum approximation explicitly select
`--kubo-source wavecar` (or INI `kubo_source = wavecar`) to reproduce their
existing results. See the [standard protocol](WAVEDER_KUBO_PROTOCOL.md)
for PAW optical selected bands, required-pair checks and occupied compatibility.

The VASPBERRY executable reads VASP WAVECAR and writes numerical results for
analysis and plotting. These examples use **Intel oneAPI Fortran and Intel MPI**;
the commands below execute VASPBERRY directly, without a Python script. The
[hands-on guide](HANDS_ON.md) follows the saved files through transport and plots.

This guide follows VASPBERRY 1.6.6. Its Kubo range selection replaces the
native `--bundle`/`-kubo_bundle` flags; see the
[migration guide](MIGRATION.md#native-kubo-band-selection) when updating from 1.6.2.

## Learn one task or option at a time

VASPBERRY includes topic help in the native Fortran executable.
After [building](BUILD.md), use `--help task` to choose a calculation,
`--help kubo` to read its inputs and outputs, and `--help bands` to check an
option before executing the example below:

```bash
mpiexec -n 4 ./build/vaspberry-ifx-mpi --help task
mpiexec -n 4 ./build/vaspberry-ifx-mpi --help kubo
mpiexec -n 4 ./build/vaspberry-ifx-mpi --help bands
```

These commands need neither Python nor material input files. Replace `kubo`
with another calculation, such as `spin-chern` or `spin-kubo`, to read that
task's requirements. Replace `bands` with an option name, such as `mesh`.
Run help as a separate command, without calculation arguments. `--help --bands`
also selects the `bands` topic; `-h` accepts the same topics as `--help`.

| Help request | Shows |
|---|---|
| `--help` or `-h` | Short overview and where to go next |
| `--help task` or `--help tasks` | Available calculations |
| `--help kubo` | One task's purpose, inputs, outputs and example |
| `--help options` | Available option topics |
| `--help bands` | One option's values and usage |
| `--help all` or `--help legacy` | Complete flag reference, including historical options |

## A runnable first command

In Bash on a Linux oneAPI host, start at the repository root and use the
included MoS₂ path WAVECAR. Replace the setup path with your site's setup or
module environment when needed:

```bash
source /opt/intel/oneapi/setvars.sh
make ifx-mpi
VB_ROOT="$PWD"
VB_BIN="$VB_ROOT/build/vaspberry-ifx-mpi"
mkdir results-native-path-01
cd results-native-path-01
mpiexec -n 4 "$VB_BIN" --task kubo --kubo-source wavecar \
  --wavecar "$VB_ROOT/examples/1H-MoS2/KPATH/2.band/WAVECAR" \
  --bands 1:18 --curvature-csv KUBO.csv
cd "$VB_ROOT"
```

This writes 48 occupied-bundle curvature rows to `KUBO.csv`. It uses every
stored empty band as an intermediate state; it does not integrate a path.
The [curvature example](../examples/features/kubo-curvature/) supplies the
matched full-mesh calculation and figure commands. To repeat this command,
choose another fresh directory.

For Intel Classic already installed at the site, use `make ifort-mpi` and
`build/vaspberry-ifort-mpi`. For serial Intel execution use `make ifx` and
`build/vaspberry-ifx` without `mpiexec`. [Build details](BUILD.md) retain GNU
alternatives and distinguish compiler smoke tests from numerical validation.

`mpiexec -n 4` launches four MPI ranks; the remaining arguments select the
scientific task and its files. All ranks read the same WAVECAR, and the native
program assembles the output. Use the same Intel MPI environment at build and
run time and a rank count consistent with the scheduler allocation.

## Choose the result

All commands accept `--wavecar PATH` and detect the spinor components
automatically. The mesh sizes and band ranges below are example choices;
use the values appropriate to your source calculation.

The Chern-number tasks use the Fukui–Hatsugai–Suzuki (FHS) link-variable
method, also called the Fukui method. The Z₂ task uses the separate
Fukui–Hatsugai (FH) n-field method.
Kubo-formula tasks evaluate Berry curvature or response from interband matrix
elements. Task names below are unchanged command identifiers.

| Task | Additional arguments | Numerical result and required sampling |
|---|---|---|
| `chern` | `--mesh 12,12 --bands 1:18` | Occupied-subspace Berry flux from the Fukui method in `BERRYCURV.dat`; Chern number in its header. Full periodic 2D mesh. |
| `spin-chern` | `--mesh 12,12 --bands 1:8 --spin-axis z` | Both projected-spin sector Chern numbers and flux/spectrum CSVs. Full periodic mesh and matching OUTCAR; [spin Chern number guide](SPIN_CHERN.md). Available in 1.6.1. |
| `spin-kubo` | `--kubo-source wavecar --bands 1:10` | Projected-spin sector Kubo-formula Berry-curvature proxy at every source k point; native CSVs, arbitrary path allowed. Explicit `--mesh NX,NY` adds a raw approximation integral. [Spin-sector Kubo-formula guide](SPIN_KUBO.md); available in 1.6.1. |
| `z2` | `--mesh 12,12 --bands 1:18` | `Z2_FIELD.csv` and `NFIELD.dat` after PASS. Full even Γ-centered mesh, SOC, time-reversal symmetry, occupied rank even and gap open; [full requirements](Z2_FUKUI_HATSUGAI.md). |
| `kubo` | `--curvature-csv KUBO.csv` | Inferred insulating occupied trace; complete same-run optical files. |
| `kubo` | `--bands 31,33:34 --curvature-csv SELECTED.csv` | Selected geometric WAVEDER trace Ωxy in Å²; all source intermediate bands and required pairs checked. |
| `kubo` | `--bands 31:32 --per-band 1 --curvature-csv BANDS.csv` | Separate selected WAVEDER bands when individually resolvable; no occupations applied. |
| `kubo` | `--kubo-source wavecar --bands 18 --curvature-csv BAND18.csv` | Individual-band Ωxy, energy and minimum gap, plus legacy DAT companions. Inspect isolation before interpretation. |
| `kubo` | `--kubo-source wavecar --bands 18:19 --per-band 1 --curvature-csv BANDS18_19.csv` | Separate curvature rows for bands 18 and 19. Each must be isolated from every other source band by more than 1e−5 eV. |
| `kubo-pairs` | `--kubo-source wavecar --pairs-csv PAIRS.csv` | Every source-band pair's undivided momentum numerator, energies and k coordinates; all yz/zx/xy components. No band selector. |
| `kubo-integral` | `--mesh 12,12 --bands 1:18 --curvature-csv KUBO.csv` | Same bundle CSV plus a finite-grid integral in stdout. Requires the full mesh; this diagnostic is not an occupation-weighted transport scan. |
| `optical` | `--mesh 12,12 --bands 18:19 --theta 0 --phi 0` | Selected transition 18→19 circular selectivity in `BERRYCURV.dat` (legacy filename). |
| `spectrum` | `--mesh 12,12 --bands 1:20 -ien 1 -fen 3 -nediv 201 -sigma 0.05` | Broadened left/right transition spectra in `OPT_TRANS_RATE_LEFT.dat` and `OPT_TRANS_RATE_RIGHT.dat`. Energies/broadening in eV. [Optical workflow](../examples/features/circular-dichroism/). |
| `wavefunction` | `--wavefunction-band 18 --kpoint 24 --real-grid 24,24,64 --imaginary 1` | CHGCAR-style `PARCHG-W-K024-E018-SPIN1` and its `-IM-SPIN1` companion. Index 24 is Γ in the included 48-point path only. Matching POSCAR/EIGENVAL must be in the working directory. [Wavefunction workflow](../examples/features/wavefunction/). |
| `velocity` | `--bands 18` | Canonical-momentum velocity expectation in `VEL_EXPT.dat`, x/y components in m/s; a single-band diagnostic. |

`kubo-line` is a synonym of `kubo`; both preserve the supplied k-point order.
For the default WAVEDER source, `--bands N`, `--bands FIRST:LAST` and
`--bands N,FIRST:LAST` select a geometric target space; the intermediate sum
still covers all source NBANDS. Omit the selector for the validated insulating
occupied-input route. Selected traces must include complete producer-degenerate
groups (2 meV transitive clusters), and every noncanceling pair must be
available in at least one orientation. `--per-band 1` rejects an unresolved
individual band. Selected mesh integrals have unit geometric weights and no
total-Hall label; Python `kubo-hall --bands` supplies explicit μ/T weighting.
For the explicit `--kubo-source wavecar` approximation, `--bands N` selects one band;
a multi-band range such as `--bands 1:18` selects the trace of that entire
subspace. Internal degeneracies are allowed in this trace, but its gap to
every excluded source band must exceed 1e−5 eV. Add `--per-band 1` only when
you need separate bands: every selected band must then satisfy the same gap
threshold against every other source band, including other selected bands.
If a multi-band trace has no explicit `--curvature-csv`, it writes `KUBO.csv`.
Optical band endpoints
instead select the initial/final transition; `spectrum` uses the selected band
window and source occupied boundary. `velocity` uses the first selected band,
so specify a single band.

Prefer one `--bands` selector. For compatibility, a later `--bands` replaces
an earlier selection, and legacy `-ii`, `-if` or `-is` after a scalar/range
selection retain their historical endpoint/single-band overrides. These
overrides print a warning; inspect the effective band IDs. Legacy selectors
cannot modify an earlier comma-separated list; a later modern `--bands`
can replace earlier legacy settings. A legacy `-ii FIRST -if LAST` pair also
remains available by itself. Modern standard WAVEDER Kubo tasks reject controls
from unrelated tasks, including `--kpoint`, real-space grid/imaginary controls,
optical angles and spectral broadening/energy-grid controls. They evaluate
every k point stored in the supplied source; no unrelated flag selects a
single point or changes the Kubo response.

## A small set of common arguments

| Argument | Meaning |
|---|---|
| `--task NAME` | One calculation. Explicit task names reject conflicting legacy task flags. Default is `chern`, which uses the Fukui method. |
| `--input-dir DIR` | Directory for default input filenames; omitted means the invocation working directory. Does not change output paths. |
| `--wavecar PATH` | Override only WAVECAR, whose default is `input-dir/WAVECAR`. Never changes the directory used for other files. Quote paths containing spaces. |
| `--kubo-source waveder\|wavecar` | Default `waveder`: same-run PAW optical input. Explicit `wavecar` opts into the warned canonical approximation. |
| `--waveder PATH` / `--incar PATH` | Override only the named standard-Kubo input; defaults are `input-dir/WAVEDER` and `input-dir/INCAR`. |
| `--outcar PATH` | Override the standard-Kubo producer record or spin-frame metadata; default `input-dir/OUTCAR`. |
| `--sum-bands N` | `spin-kubo` only: use source bands 1:N outside the selected group in the intermediate sum; default all stored bands. This keeps the VASP eigenstates fixed for sum-convergence checks. |
| `--spin-axis z` | Cartesian analysis axis for both spin-sector tasks: `x`, `y`, `z` or a comma-separated unit vector, e.g. `0,0,1`. Default z; OUTCAR defines the source spin frame. |
| `--energy-gap-tol VALUE` / `--spin-gap-tol VALUE` | Both spin-sector tasks: default 1e-8 eV for the selected energy subspace and 1e-6 for the dimensionless projected-Pauli distance from zero. |
| `--spinor auto` / `--spinor 1` / `--spinor 2` | Default `auto`: detect the component count from WAVECAR. Explicit `1` or `2` checks the file layout; it is not a spin-degeneracy factor. |
| `--mesh NX,NY` | Mesh dimensions for mesh-based algorithms. For standard WAVEDER Kubo, this explicitly requests full-grid validation and integration; legacy `-kx`/`-ky` alone only set dimensions. It does not generate k points, interpolate, or convert a path into a mesh. |
| `--bands FIRST:LAST` / `--bands N` | Inclusive one-based range or singleton. WAVEDER Kubo also accepts a comma-separated list such as `31,33:34`; it selects a geometric trace unless `--per-band 1`. |
| `--per-band 1` | Native `kubo`, `kubo-line` and `kubo-integral` only: calculate separate isolated bands rather than the multi-band trace. |
| `--curvature-csv PATH` / `--pairs-csv PATH` | Exact CSV output filename. Parent directory must already exist. |
| `--output LABEL` | Output naming label, not a directory. For `spin-chern`, default `SPIN` gives `SPIN_CHERN.csv`, `SPIN_BERRY.csv`, `SPIN_SPECTRUM.csv`; it does not replace explicit Kubo CSV paths. |

All explicit relative CLI paths are relative to the invocation working
directory, including per-file overrides. They are not rebased onto
`--input-dir`. Output filenames and prefixes keep their existing meaning.
WAVEDER/INCAR are required by standard charge Kubo, not by every native task.
The native directory option affects only omitted WAVECAR/WAVEDER/INCAR/OUTCAR
paths. It does not relocate unrelated legacy POSCAR/EIGENVAL files, input
lists, or outputs; those keep their existing working-directory rules.

For `spin-kubo`, the same default prefix gives `SPIN_KUBO.csv` and
`SPIN_KUBO_SPECTRUM.csv`; an explicitly supplied full `--mesh` additionally
gives `SPIN_KUBO_INTEGRAL.csv`. No integer invariant is assigned by this task.

Pass option values as separate arguments: use `--mesh 12,12`, not
`--mesh=12,12` or `--mesh 12 12`. Lists contain no spaces. Modern names and
legacy flags can be mixed where needed; optical energy-grid controls above
retain their short names. There is no need to learn every historical flag.

Two components identify a noncollinear spinor, not whether VASP enabled SOC.
The spin-axis convention still comes from the matching OUTCAR.

The current native reader detects byte and legacy four-byte-word RECL
layouts automatically for supported WAVECAR coefficients (`RTAG=45200`).
It validates the header and record layout without rewriting the input.
Keep the supplied Intel `-assume byterecl` build flag: it sets the reader's
internal units, not a restriction on the detected file layout. Retain the
original VASP producer/compiler provenance.

## Where outputs go

Use a fresh working directory per native calculation. Select the source
directory explicitly with `--input-dir`; its default remains the invocation
working directory. Individual input overrides do not move the output files.
Legacy DAT and wavefunction outputs can replace existing files. Kubo
and spin Chern number CSV exporters refuse an existing filename; they never append another
calculation to it. Current native curvature and pair exports must end with
`result_status=PASS`. Trace CSV V2 and individual-band CSV V3 begin as
INCOMPLETE and record their expected row coverage; a header or partial file
is not a completed calculation. Curvature and pair exporters use a `.partial` file
and publish the final filename only on successful completion. Existing
final or partial filenames are refused. Spin Chern number and Z₂ validation status
must also be checked before use.

Legacy labels are retained for compatibility: `--output sample` produces
`BERRYCURV.sample.dat` for Berry flux from the Fukui method, `sample.dat` for optical selectivity,
`CIRC_DICHROISM_W.sample_LEFT/RIGHT.dat` for spectra, and
`VEL_EXPT.sample.dat` for velocity. Individual-band Kubo adds `.EIG-N`, and
its sum has the base label. Scalar `ISPIN=2` DAT outputs add `.UP`/`.DN`.
Wavefunction filenames are selected by k and band; `--output` does not rename
them. `--curvature-csv PATH` chooses an exact filename for either charge-Kubo
source. If omitted, standard WAVEDER uses `KUBO_WAVEDER.csv` for occupied,
selected-trace and separate-band results. The explicit WAVECAR approximation
defaults to `KUBO.csv` for a multi-band trace; its separate-band CSV needs
an explicit path. Pair export requires its explicit CSV path.

The native program currently stores most paths in 256-character fields.
Prefer short working paths. The CSV interfaces are the preferred precision
and interchange route for Kubo results; legacy DAT tables have fewer printed
digits. [Output schemas](OUTPUT_FORMAT.md) define units, columns, index bases,
operator assumptions, and the NPZ/JSON caches used by postprocessing.

## Which file can be plotted immediately?

| Native output | Contents | Next use |
|---|---|---|
| Native trace CSV: `KUBO_WAVEDER.csv` by default for standard WAVEDER, `KUBO.csv` for the explicit WAVECAR approximation | `spin,k_index,kx_frac,ky_frac,kz_frac,omega_z_A2,min_external_gap_eV` | Plot the computed Ωxy in Å² against k index/path position, or use the lattice for a BZ map; retain the operator metadata |
| Individual-band curvature CSV | k/band IDs, fractional k, band energy, gap and `omega_z_A2` | Plot energy and curvature for the selected isolated band |
| `PAIRS.csv` | k/spin and n/m IDs, two energies, gap and three `numerator_*_eV2_A2` columns | Integrate occupations and squared energy denominators to obtain Hall response |

Native CSV files start with `#` metadata lines and then the named column
header. To plot bundle curvature in Origin, a spreadsheet or another CSV tool,
skip the comments, select `spin=1`, and use `k_index` and `omega_z_A2` as x/y.
This does not need Python. Preserve the metadata when exporting a table.

The pair numerators are `Nab=-2 Im(Da,nm Db,mn)` in eV² Å² for n<m. They are
not complex velocity matrices, curvature or conductivity: the squared gap
division and occupation weighting have not yet been applied. The
[output specification](OUTPUT_FORMAT.md#native-kubo-pair-export-and-reusable-cache)
lists their exact columns and units.

## Charge transport and reusable figures

The standard charge-Hall workflow reads the optical matrix elements directly:

```text
Same-run full-mesh WAVEDER + WAVECAR + INCAR + OUTCAR
  -> Python kubo-hall: validated optical connections, occupations and BZ integration
  -> conductivity.csv / .dat / .npz + .json
  -> plot_hall.py or your preferred plotting program
```

Use `kubo-hall --bands SELECTOR` for an occupation-weighted selected
contribution, or `--occupied N` for the insulating T=0 response. These
selectors are mutually exclusive. The [standard protocol](WAVEDER_KUBO_PROTOCOL.md)
gives complete commands and required inputs. The INI `run` command selects
the same backend by default. A selected sum over the `total` k-space region
still denotes the selected bands' contribution, not certified total AHC.

The explicitly requested WAVECAR approximation retains a reusable pair cache:

```text
VASP full-mesh WAVECAR
  -> VASPBERRY --task kubo-pairs --kubo-source wavecar --pairs-csv PAIRS.csv
  -> Python import-pairs: validation/cache -> pairs.npz + pairs.json
  -> Python pair-hall: transport integration -> conductivity.csv / .dat / .npz + .json
  -> plot_hall.py or your preferred plotting program
```

Pair export does not use source occupations or divide by energy gaps;
k-dependent metallic/smeared occupations therefore need no manual `-ne`
override. Postprocessing supplies μ, temperature, spin multiplicity and
regions. Its complete, uniform 2D mesh and degeneracy checks still apply.
Once this canonical cache exists, varying μ, temperature, regions or the
retained pair-band window does not require another VASPBERRY execution.
`import-pairs` checks and reorganizes the native output. `pair-hall` performs
the occupation-weighted transport integral. `plot_hall.py` only reads the
finished conductivity table and draws figures. The optional `wavecar-hall`
command executes VASPBERRY, then imports its pair data and integrates in one call;
it records the same intermediate files.
See the [three-stage commands](KUBO_TRANSPORT.md#native-pairs-to-charge-hall)
and [actual MoS₂ Hall example](../examples/features/kubo-hall/).

For atom/layer/orbital and chosen-axis spin character, retain the matching
SOC `PROCAR` and `OUTCAR` as well. The separate
[PROCAR workflow](../examples/features/procar-character/) accepts user-defined
groups, combines their projections with the saved canonical pair cache for
selected isolated bands, and writes charge-attribution tables and plots.
The saved VASPBERRY pair data can be reused for these projections.

Native Kubo-formula calculations from WAVECAR use canonical momentum of
pseudo-wavefunctions. This approximation does not supply missing
PAW/nonlocal/SOC velocity terms. These outputs support
the documented intrinsic charge-response approximation, not a complete
longitudinal conductivity, relaxation-time transport or spin-Hall calculation.
For additional physical operators and spin Hall, use the separately documented
[operator routes](OPERATOR_ROUTES.md). Wannier is an optional model route.

Regenerate historical native `-vel 1` outputs: the previous velocity conversion
omitted the factor c² needed with electron mass expressed in eV/c², and its
extrema used an invalid k index. Rounded legacy zero values cannot be repaired
afterward. The corrected diagnostic prints scientific-notation m/s values and
still represents only canonical momentum, not the full material velocity.

The task arguments are identical for Intel serial/MPI and GNU builds; choose
the executable and launcher shown in the [build guide](BUILD.md).
