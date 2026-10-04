# Postprocessing settings reference

Use the [step-by-step guide](POSTPROCESSING.md) first. This page lists the
accepted settings for `tools/vaspberry_post.py`; optional settings can simply
be omitted. Start from the [selected WAVEDER software example](../examples/features/waveder-selected/)
or the canonical [Bi example](../examples/features/simple-postprocess/README.md)
and its [initial scan](../examples/features/simple-postprocess/bi.ini) or
[rescan](../examples/features/simple-postprocess/bi-rescan.ini).

## File syntax and names

- Sections and keys are case-sensitive. Use `key = value`, space-separated
  lists, and `#` or `;` comments. Inline comments need preceding whitespace.
- Explicit relative paths in the INI are relative to that INI's directory,
  including `input_dir` and individual file overrides. Absolute paths and
  `~` work; shell variables, commands and quoted shell syntax are not expanded.
  Write a path containing spaces without adding quotes.
- Omitted `[run] input_dir` means the invocation working directory, not the
  INI directory or the parent of a `wavecar` override. Input-directory selection
  never changes how output paths are resolved.
- Each section and key may appear once. Unknown keys, empty values and
  `[DEFAULT]` are rejected. Merge additions into an existing `[plot]` section.
- Group, region and difference names begin with an ASCII letter, followed by
  letters, digits, `_`, `.` or `-`. Names such as `layer1`, `lower`, `M` and
  `m` are labels, and `M` and `m` are different names. Use the same spelling
  wherever a plot or difference refers to a label.
- Atom, band and k-point IDs start at **1**. `plane_axes` is the exception:
  reciprocal-basis axis indices start at **0**.

## Required settings

Both `[run]` and `[hall]` are required. From 1.6.5 the default is
`[run] kubo_source = waveder`. Same-run optical input is mandatory for that
route; missing files or unsupported calculations stop. Set `kubo_source =
wavecar` explicitly to run the warned canonical-momentum approximation.
The [protocol](WAVEDER_KUBO_PROTOCOL.md) supplies a standard INI example.
WAVEDER requires exactly one `[hall] bands` or `occupied` selector. The
selected route rejects missing nonzero-weight pairs and unequal weights
inside producer-erased groups. A finite-temperature selection remains a
contribution; it does not become the total AHC.
`--reuse`, pair cutoffs and PROCAR attribution belong to the explicit
WAVECAR route and are rejected with WAVEDER.

| Section / key | Value and meaning |
|---|---|
| `[run] input_dir` | Directory for default source filenames. Optional; defaults to the invocation working directory. Explicit relative values are relative to the INI. |
| `[run] wavecar` | Optional override of only WAVECAR; default `input_dir/WAVECAR`. Its k points must form a complete uniform 2D mesh, not an irreducible mesh or path. |
| `[run] waveder`, `incar`, `outcar` | Optional per-file optical overrides, resolved relative to the INI when relative. Each changes only its named file; defaults are those filenames in `input_dir`. WAVEDER route only. |
| `[run] kubo_source` | `waveder` (default) or explicit `wavecar`. The latter prints an approximation warning. |
| `[run] optical_run_dir` | Legacy WAVEDER-only alias for `input_dir`. Mutually exclusive with `input_dir`; rejected with the WAVECAR approximation. No directory is inferred from WAVECAR. |
| `[run] binary` | Compiled VASPBERRY path, required only for the explicit WAVECAR route, including its `--reuse` configuration. Rejected by the Python WAVEDER backend; omit it for that route. |
| `[run] output` | New result directory. An existing directory is rejected, including an empty one. |
| `[run] mesh` | `NX NY`, two integers at least 2. Their product equals the source k-point count; coordinates must match that uniform grid. |
| `[run] energy_reference` | Descriptive text recording the WAVECAR energy convention, e.g. `unchanged VASP eigenvalue zero`. It does not shift any energy. |
| `[hall] mu` | `MIN MAX N`: N equally spaced chemical potentials in eV, including endpoints. Use `MIN < MAX` and integer `N >= 2`, or `MIN = MAX` and `N = 1` for a single Hall-response point. |
| `[hall] reference` | Reference chemical potential in eV for Δσ and the plot's horizontal origin. Uses the same energy zero as `mu` and WAVECAR; need not be a sampled `mu` value. |
| `[hall] temperatures` | Distinct nonnegative temperatures in K. WAVEDER `occupied` mode requires `0`; selected `bands` mode checks all requested weighted pairs, including finite-temperature tails. |
| `[hall] occupied` | WAVEDER insulating compatibility selector; μ/reference strictly inside the global gap, T=0. Mutually exclusive with `bands`; rejected with WAVECAR source. |
| `[hall] bands` | WAVEDER target selector, e.g. `31`, `31:32` or `31,33:34`; mutually exclusive with `occupied`. All source intermediate bands remain included. Gives a selected contribution, not certified total AHC; rejected with WAVECAR source. |

`[run] spin_mode` is optional and defaults to `auto`.

| `spin_mode` | Source and counting |
|---|---|
| `auto` | Detect components from WAVECAR: two-component states count once; scalar `ISPIN=1` states have multiplicity 2. For `ISPIN=2`, choose a channel below. |
| `spinor` | Assert two-component `ISPIN=1` wavefunctions, counted once. |
| `soc` | Compatibility alias for `spinor`; it does not establish that SOC was enabled. |
| `scalar-degenerate` | Scalar states, WAVECAR `ISPIN=1`; spin multiplicity 2. |
| `collinear-up` | First scalar channel of an `ISPIN=2` WAVECAR; multiplicity 1. |
| `collinear-down` | Second scalar channel of an `ISPIN=2` WAVECAR; multiplicity 1. |

A collinear run gives one channel's response. It does not sum both channels.

### Energy zero, reference and scan

The front end uses the eigenvalues stored in WAVECAR without shifting them.
It does not read a Fermi energy from OUTCAR to set the energy zero. For
example, if your chosen reference is 5.2 eV on that eigenvalue scale and you
want ±0.5 eV around it, write:

```ini
[hall]
mu = 4.7 5.7 101
reference = 5.2
temperatures = 0 300
```

The CSV retains `mu_eV` and `mu_minus_reference_eV`; figures use the latter
on the horizontal axis. `sigma` shows σ(μ,T), while `delta-sigma` shows
σ(μ,T) − σ(reference,T), evaluated using occupation differences. Changing
`reference` requires a new numerical run; the explicit WAVECAR route can
reuse its pair cache. The standard WAVEDER route rereads its optical inputs.
Changing `energy_reference` text alone does not perform an energy conversion.

## Common optional settings

### Executing the explicit WAVECAR approximation with MPI

| `[run]` key | Default | Meaning |
|---|---|---|
| `mpi_procs` | `1` | Positive integer. Values greater than 1 execute VASPBERRY with `MPI_LAUNCHER -n N`; 1 runs it directly. |
| `mpi_launcher` | `mpiexec` | One launcher executable name on `PATH`, or a path. No launcher options in this value. For Intel MPI, e.g. `mpiexec.hydra`. |

Run the Python helper once. Its `run` command executes the VASPBERRY
binary selected in `[run] binary`; `mpi_procs` applies to that VASPBERRY
calculation. The subsequent Python postprocessing runs once. Use a matching
VASPBERRY MPI build and launcher within an allocation allowing those ranks.
See [build and runtime requirements](BUILD.md).

### Atom, layer, orbital and spin character

`[projection]` is optional. It requires detected two-component wavefunctions
and at least one `[group NAME]`. Groups also require `[projection]`.

| `[projection]` key | Default | Meaning |
|---|---|---|
| `bands` | Required | Unique 1-based band IDs, e.g. `12 13`. These select the bands included with their occupations in the character-weighted Hall sum; the virtual-state sum still uses all stored pair bands. |
| `axis` | Required | Three Cartesian components of a unit spin-analysis vector, e.g. `0 0 1`. It must already have unit norm; no automatic normalization. |
| `procar` | `input_dir/PROCAR` | Matching noncollinear projection file. |
| `outcar` | `input_dir/OUTCAR` | Matching completed static-run output, including the SAXIS-to-Cartesian transformation. |
| `mu` | `[hall] mu` | Same `MIN MAX N` syntax, but requires `MIN < MAX` and `N >= 2`. |
| `reference` | `[hall] reference` | Reference μ in eV for the character-weighted response. |
| `temperatures` | `[hall] temperatures` | Distinct nonnegative temperatures in K. |

| `[group NAME]` key | Default | Meaning |
|---|---|---|
| `ions` | Required | Unique 1-based atom IDs in the POSCAR/PROCAR ordering, e.g. `1 3 5`. The label does not select atoms automatically. |
| `orbitals` | All orbitals in each selected ion's printed total | Optional exact orbital labels from the PROCAR header, e.g. `dxy dyz dz2 dxz x2-y2`. |

For example, `[group layer1]` with `ions = 1 3` and `[group layer2]`
with `ions = 2 4` names two atom groups. Neither `layer1` nor `upper` has
special meaning; choose IDs from your structure. Groups may overlap, but
overlapping groups do not provide an additive partition.

Use PROCAR, WAVECAR and OUTCAR from the **same unmodified final static run**.
The supported PROCAR reader handles the noncollinear `LORBIT=11` layout
with charge and three magnetization blocks. OUTCAR must report
`LNONCOLLINEAR=T`, `ISPIN=1`, `NSW=0`, `LORBIT=11`, completed execution and
one unambiguous spin rotation. The numerical stage checks state ordering,
coordinates, energies, occupations, dimensions and lattice; passing `check`
alone does not run all those matching checks.

Each selected band must be isolated from every other stored band at every
k point: gaps at or below 10⁻⁵ eV are rejected for character-weighted Hall.
Raw character weights are not normalized to unity. The `charge`,
`pauli_axis`, `plus` and `minus` components describe projected character;
`plus/minus = (charge ± pauli_axis)/2`. Their weighted Hall contributions
have charge-Hall units, not independent spin- or layer-current units.
See the [analytic PROCAR example](../examples/features/procar-character/README.md)
and [character output columns](OUTPUT_FORMAT.md#procar-character-and-charge-hall-attribution).

### Reciprocal-space regions

Regions are optional. `total` is the full mesh and `rest` is the complement
of all user regions; both are automatic and reserved names. With no user
regions, `rest` coincides with `total`.

| `[region NAME]` key | Default | Meaning |
|---|---|---|
| `center` | Required for a disk | Three reciprocal fractional coordinates. Fractions such as `1/2` work here. The center must lie in the sampled 2D plane, modulo periodicity. |
| `radius` | Required for a disk | Positive Cartesian reciprocal distance in Å⁻¹. Periodic images are restricted to the sampled plane. |
| `k_ids` | Alternative to a disk | Unique 1-based source k-point IDs. Use this key alone, without `center` or `radius`. |

For example, `[region M]` with `center = 1/2 0 0` and `radius = 0.10`
selects the disk about those coordinates. Use that center only if it is the
M point in **your reciprocal basis**. Naming a region `M` does not locate
M automatically, interpolate missing points or request a band-path plot.

Every user region must contain sampled k points. Distinct regions cannot
share sampled points; two disks additionally cannot overlap or touch
geometrically, even if a coarse mesh samples no shared points. These tests
run during integration. See the [regional Hall example](../examples/features/kubo-hall/README.md)
for a material-specific K/K′ partition.

| `[differences]` entry | Meaning |
|---|---|
| `NAME = LEFT RIGHT` | Adds the **LEFT − RIGHT** charge-Hall channel, without a factor of one half. Both names must already exist: user regions, `total`, `rest`, or an earlier difference in the same section. The new name cannot duplicate an existing region/difference. |

Differences apply to the total Hall table. Character-weighted Hall supports
the regions themselves, not difference channels.

### Figure defaults

`[plot]` is optional. These preferences are saved at `run` time; they do not
change which numerical data are computed. `plot` reads the saved preferences.

| `[plot]` key | Default | Meaning |
|---|---|---|
| `hall_regions` | `total` | Space-separated existing region or difference names. |
| `hall_temperatures` | All `[hall] temperatures` | Saved temperatures to draw; no new temperatures are calculated. |
| `hall_quantity` | `sigma` | `sigma` or `delta-sigma`. |
| `character_group` | First group in the file | Saved group name; requires `[projection]`. |
| `map_band` | First `[projection] bands` entry | 1-based band for the raw character map. May be any stored band, independently of the selected Hall band sum. |
| `character_temperature` | First `[projection] temperatures` entry | Saved temperature for character-weighted Hall curves. |
| `character_region` | `total` | A user region, `total` or `rest`, not a difference channel. Selects the Hall curve; the raw character map still covers the full stored mesh. |
| `character_delta` | `false` | `true` for reference-subtracted character-weighted Hall curves; `false` for absolute curves. |

The last five keys require `[projection]`. Boolean values accept `true` or
`false`, regardless of case. Character maps do not depend on temperature or
the Hall reference; the accompanying Hall curves do.

## Advanced settings and numerical scope

| Section / key | Default | Meaning |
|---|---|---|
| `[run] plane_axes` | `0 1` | Two distinct reciprocal-basis indices from `0 1 2`. `mesh` sizes follow this order; the remaining fractional coordinate is fixed. The ordered cross product sets the oriented integration normal. |
| `[hall] pair_band_max` | All stored bands | WAVECAR approximation only; rejected for WAVEDER. Integer from 2 through the source band count. Both endpoints of retained pairs lie in bands `1..pair_band_max`. Does not reduce the VASPBERRY pair export or the virtual-state sum in character-weighted Hall. |

The compact front end uses the strict numerical defaults of the underlying
tools. In the explicit WAVECAR pair-Hall route, pairs with gaps at or below 10⁻⁷ eV must have equal
occupations at every requested μ/T and the reference. A retained-band
cutoff must not split an unresolved group, and its highest band must be
unoccupied within the numerical tolerance (10⁻⁸). These checks do not
replace convergence of the source bands, pair window, k mesh and operator
approximation. Controls such as explicit degeneracy coalescing or partial
band acceptance are not INI keys; see [advanced Kubo-formula transport commands](KUBO_TRANSPORT.md).

## Commands, saved results and reuse

These commands use the prefix `python3 tools/vaspberry_post.py`. For direct
VASPBERRY commands that produce a curvature table without this helper, see
[direct execution](POSTPROCESSING.md#execute-vaspberry-directly-for-a-curvature-table)
and the [hands-on guide](HANDS_ON.md).

| Command | Input dependencies and action |
|---|---|
| `check analysis.ini` | Checks settings, a new output path and WAVECAR/header/grid. WAVEDER validates same-run optical files and evaluates the requested μ/T spectrum to check weighted pairs and producer clusters, without saving the spectrum; it may take comparable time to `run`. The explicit WAVECAR route checks its executable/launcher and any projection inputs without executing the native pair calculation. Creates no result directory. |
| `run analysis.ini` | Performs the checks and computes Hall tables. WAVEDER contracts optical connections directly; the explicit WAVECAR route executes VASPBERRY to write `native/PAIRS.csv` (with MPI when configured), then imports the pairs and any character attribution. Saves settings, source/operator scope, logs and table checksums. |
| `check analysis-next.ini --reuse results/run01` | Requires a completed previous front-end run, its intact pair cache and the same WAVECAR. Checks compatibility without calculating. Executable and launcher availability are not needed. |
| `run analysis-next.ini --reuse results/run01` | Copies the verified pair cache into a new output directory and recalculates the requested postprocessing. Requires source PROCAR/OUTCAR again if projection is enabled. |
| `plot results/run01` | Uses a completed `run.json` and its intact saved Hall tables; with projection, also the saved character and character-Hall data. Requires no source WAVECAR, PROCAR, OUTCAR, VASPBERRY executable or MPI launcher. |

`--reuse` is available only for the explicit WAVECAR approximation. It requires
an identical WAVECAR hash and unchanged mesh, plane axes,
spin convention and `energy_reference` text. Paths may change if the file
content is identical. Scan values, groups, regions and projection selections
can change; these are recalculated. Reuse skips VASPBERRY execution and
recomputes the numerical occupation integral. The original result remains
unchanged.

For plotting, optional command-line overrides are:

| Option | Effect |
|---|---|
| `--output-dir PATH` | New figure directory; default `RESULT/figures`. Existing directories are rejected. |
| `--temperature T` | One saved temperature in K. With projection, T must be present in both total and projected scans. |
| `--quantity sigma` / `--quantity delta-sigma` | Sets both total and character-weighted Hall curves to absolute/reference-subtracted values. |
| `--group NAME` | Another saved character group; requires projection results. |
| `--band N` | Another saved band's raw character map; the Hall band sum is unchanged. |

Command-line paths are relative to the current working directory, unlike
paths inside an INI. Editing the original INI, including renaming a group or
region, does **not** change a completed result. Use a new run with `--reuse`
for changed numerical definitions or saved labels. Preserve saved tables
and metadata together; the front end checks their recorded hashes.

For filenames, data columns and units, use the
[output specification](OUTPUT_FORMAT.md). The JSON/NPZ files written by the
front end are saved metadata/caches, not additional settings to author.

## Errors and fixes

| Message or symptom | What to check |
|---|---|
| `unknown key`, `unknown section`, or duplicate section | Use the spelling in this reference; combine all plot keys under one `[plot]`. Lists use spaces, not commas. |
| `output directory exists` / `figure directory exists` | Choose a new `[run] output` or `--output-dir`; keep the previous result. |
| WAVECAR, executable or MPI launcher missing | Resolve INI paths from the INI's directory. Build the executable and load its compiler/library/MPI environment; see [BUILD.md](BUILD.md). |
| Mesh mismatch, incomplete/off-grid points | Use the full uniform 2D WAVECAR and its correct `mesh`/`plane_axes`. A band path or irreducible mesh cannot be integrated here. |
| `region ... has no sampled points` | Check reciprocal basis, sampled plane and radius/IDs. A name such as M does not add a k point. |
| Regions overlap or touch | Define disjoint regions; disks must also be geometrically separated under periodicity. |
| `highest ... band occupied` / cutoff splits a group | Increase the retained/source band window as appropriate, or reconsider a scan extending beyond that window; converge the virtual-state sum. |
| Unequal occupation in an unresolved pair | Review the source gaps and μ/T/reference. See [Kubo-formula numerical policy](KUBO_TRANSPORT.md); no INI key silently bypasses the check. |
| PROCAR/OUTCAR mismatch or unresolved individual-state gap | Use files from one matching static run and isolated selected bands. Degenerate-band character-weighted Hall cannot be obtained by relabeling groups. |
| Saved temperature/group absent | Select a value already calculated. Otherwise create a new run, using `--reuse` when source-compatible. |
| Saved output changed or missing / reuse source mismatch | Restore the original complete result, or calculate a new one. Keep caches with their metadata and the matching WAVECAR for reuse. |

A failed numerical stage is recorded in `run.json`; its message points to
the corresponding `native/` or `logs/` files. A failed result is not accepted
as a reusable or plottable completed run.
