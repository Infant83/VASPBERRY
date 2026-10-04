# Migrating to the current source

VASPBERRY 1.6.5 changes the default Kubo input as described first below.
Version 1.6.3 simplified native Kubo band selection and strengthened
curvature-output validation as specified below. Use `--help task`, then `--help NAME` for a task or
option; the complete flag reference is `--help all` (also `--help legacy`).
The [current release](https://github.com/Infant83/VASPBERRY/releases/tag/v1.6.5),
[native reference](NATIVE_COMMANDS.md) and [hands-on tutorial](HANDS_ON.md)
provide the current setup and commands. Earlier migration requirements follow.

<a id="waveder-default-kubo-input-in-164"></a>
<a id="waveder-default-kubo-input-in-170"></a>

## WAVEDER-default Kubo input in 1.6.5

Version 1.6.5 supersedes the unpublished 1.6.4 candidate. These are the
changes from published 1.6.3; no separate 1.6.4 release is required.

`--input-dir DIR` now selects the default input directory in native commands
and Python `kubo-hall`; omit it to use the invocation working directory.
Each per-file option overrides only that file. `--wavecar elsewhere/WAVECAR`
never changes where WAVEDER/INCAR/OUTCAR are found. Explicit relative CLI
paths and all output paths remain relative to the invocation working directory.
INI `input_dir` works similarly, but explicit relative INI paths remain
relative to the INI file. The standard INI no longer requires `wavecar`.

`--task kubo` and `kubo-integral` now select `--kubo-source waveder`.
Keep same-run WAVEDER, WAVECAR, INCAR and OUTCAR. A missing file or unsupported
producer/filling stops with an explanation. There is no automatic change of
operator. Omit native `--bands` to retain the complete insulating occupied route.
An explicit single/range/list selector, such as `31,33:34`, now selects a
geometric target space; `--per-band 1` requests separate resolvable bands.
Even an explicit `1:N` matching the occupied set is labeled geometric: only
the omitted-band occupied route retains the native physical total-Hall
metadata. Pair export and spin-sector Kubo remain separate routes.

To reproduce an old canonical-momentum command, add `--kubo-source wavecar`,
including for `kubo-pairs`, `spin-kubo` and legacy `-kubo` selectors. The code
prints a warning and marks the approximation in its output. Missing PAW or
nonlocal terms are not restored by increasing the mesh or band count.

For Python, prefer `kubo-hall` (WAVEDER by default). Historical
`wavecar-hall` also requires `--kubo-source wavecar`. The INI front end
defaults to `[run] kubo_source = waveder`, with exactly one `[hall] bands`
(selected μ/T contribution) or `[hall] occupied` (insulating T=0). Required
pairs are checked across the full source virtual space. Old
canonical pair/metallic/projected recipes must explicitly set `kubo_source = wavecar`.
See [commands and supported scope](WAVEDER_KUBO_PROTOCOL.md). Historical
reference data keep their original operator and version metadata.

<a id="native-kubo-band-selection"></a>

## Historical canonical Kubo band selection in 1.6.3

Version 1.6.3 removes `--bundle` and `-kubo_bundle`. Both names now
fail with migration guidance instead of selecting a mode. The table below
describes canonical tasks `kubo`, `kubo-line` and `kubo-integral`; add
`--kubo-source wavecar` when reproducing these commands with 1.6.5:

| Requested result | Current arguments | Output |
|---|---|---|
| One isolated band | `--bands 18 --curvature-csv BAND18.csv` | Individual-band CSV and legacy DAT |
| Trace of bands 1–18 | `--bands 1:18` | `KUBO.csv`, one row per k/spin |
| Separate bands 18 and 19 | `--bands 18:19 --per-band 1 --curvature-csv BANDS.csv` | Individual-band CSV and legacy DAT |

Use `--curvature-csv PATH` to name the trace file explicitly. Trace CSV V2
and individual-band CSV V3 keep the numerical columns and add completion
and coverage metadata; see [output formats](OUTPUT_FORMAT.md#native-kubo-bundle-csv).
The physical bundle concept and Python `bundle-hall` integration remain valid.
Older schemas are accepted as legacy-unverified inputs, with their original
metadata preserved; rerun when the new completion guarantee is required.

A named Kubo task applies this selection rule even with legacy `-ii`/`-if`
endpoints. A purely legacy `-kubo` command retains individual-band output.
Single/per-band modes now reject any selected-to-other-source-band gap at
or below 1e-5 eV before writing results. The trace checks only selected-to-
excluded gaps and permits internal degeneracies. Pair, spin-sector and
Chern-number tasks retain their existing semantics.

These are 1.6.3 interface changes; fixed 1.6.2 binaries retain their original behavior. Use the [tagged 1.6.2 guide](https://github.com/Infant83/VASPBERRY/blob/v1.6.2/docs/NATIVE_COMMANDS.md)
when reproducing that release's commands.

<a id="migrating-to-the-140-source"></a>

## Migrating from releases before 1.4.0

Version 1.4.0 is available as a versioned source release. Preserve the exact producer
version/commit and original files when migrating a calculation. The existing
Chern-number interface based on the Fukui–Hatsugai–Suzuki (FHS) method, also
called the Fukui method, remains available. The Z₂ interface retains the
separate Fukui–Hatsugai (FH) n-field method.

## Version 1.4.0 commands, projections and velocity fixes

Use the [native reference](https://github.com/Infant83/VASPBERRY/blob/v1.4.0/docs/NATIVE_COMMANDS.md)
and [hands-on tutorial](https://github.com/Infant83/VASPBERRY/blob/v1.4.0/docs/HANDS_ON.md)
for the readable commands included in the fixed 1.4.0 tag; record the exact
source commit. Original short flags remain supported. The serial executable
is `build/vaspberry`; `build/vaspberry-gfortran` remains a compatibility name.

The [PROCAR character tutorial](../examples/features/procar-character/) uses
matching wavefunctions, explicit groups and the actual spin frame. Its
projected-Hall command consumes already normalized native pair data. It does
not apply the historical factor-of-two correction. Existing private scripts
that divide legacy curvature by two must not be applied unchanged to these
outputs. Historical reference files keep their original producer metadata.

`--task kubo-pairs` now permits k-dependent stored occupations, including
metallic/smeared WAVECAR inputs. It exports all source-band pair numerators,
which contain no occupations. The later `pair-hall` command evaluates them at
the requested chemical potential and temperature. There is no need to force
an artificial `-ne` value just to export pairs. Existing checks in fixed-subspace
modes remain in place; this is not permission to treat a metal as an insulator.

Native `--task velocity` / `-vel 1` had an incorrect unit conversion from the
electron rest energy and could print all values as zero at its fixed decimal
precision. Its extrema code also used a k index after the loop, outside the
allocated range. The corrected routine reports the bare canonical-momentum
expectation in **m/s**, uses valid extrema indices, and prints scientific notation.
Regenerate old `VEL_EXPT*.dat` from WAVECAR; do not infer physical zero velocity
from those old rounded values. It remains a pseudo-wavefunction momentum
diagnostic, not the full PAW/nonlocal/SOC/U group velocity. Separate Kubo-formula
Berry-curvature/pair kernels are unchanged by the velocity correction, and no new
factor should be applied to their outputs.

## Legacy Kubo factor of two

The older Fortran Kubo path formed circular matrix elements `Px+iPy` and
`Px-iPy` without the `1/sqrt(2)` normalization. Their squared difference is

```math
|P_x+iP_y|^2-|P_x-iP_y|^2=4\mathrm{Im}(P_xP_y^*).
```

Consequently, its negated numerator contained `-4 Im` instead of the standard
`-2 Im` for `A=i<u|grad_k u>`. The corrected 1.3.0 Kubo path uses the standard
factor. At otherwise identical inputs and numerical settings, the affected
legacy curvature magnitude is halved. This is a normalization correction,
not a guarantee of physical velocity completeness or mesh convergence.

| Data source | Normalization action |
|---|---|
| Confirmed affected older Fortran Kubo output | Import with the explicit historical doubled convention |
| Corrected 1.3.0 Fortran Kubo output | Import as standard normalization |
| Generic interband-matrix Kubo result using `-2 Im` | Already normalized; do not divide by two again |
| Plaquette flux from the Fukui method / FH Z₂ n-field | Unaffected by this Kubo correction |
| Unknown or modified producer | Establish its formula from source/provenance before importing |

The `import-legacy` command requires an explicit source convention. It does
not infer a factor from a plausible Chern number, a filename, or a rounded
output. Keep the original data and the generated normalization metadata so a
later reader can tell which transformation was applied. Historical example
files remain historical and must not be silently overwritten or relabeled as
new calculations.

### Corrected Fortran output

The legacy optional `-kubo_csv PATH` writes full-precision per-band rows with spin,
k index, band, fractional coordinates, energy, minimum gap and `omega_z_A2`.
Its comments identify `STANDARD_MINUS_TWO_IM`. Request all bands present in
the WAVECAR when preparing the current importer's input; a selected n-band
subset is not a complete import file. Use `-kubo 1` for a full 2D mesh and
`-kubo 2` only for point/path output. For example, after a full-grid calculation
has produced `kubo.csv`, a **scalar, spin-degenerate** 12×12 input can be imported
as follows:

```bash
python3 tools/vaspberry_kubo.py import-legacy \
  --csv kubo.csv --wavecar WAVECAR --normalization physical \
  --spin 1 --spinor-components 1 --spin-multiplicity 2 \
  --mesh 12 12 --plane-axes 0 1 \
  --energy-reference "Unchanged WAVECAR energy zero" \
  --degeneracy-threshold-eV 1e-6 \
  --output-dir imported-kubo
```

Replace the mesh and gap threshold with the actual calculation's choices.
The current importer detects both values when they are omitted. Explicit
`--spinor-components` checks the file; `--spin-multiplicity` sets the counting.
For explicit collinear spin channels, select `--spin` and use multiplicity 1;
each imported channel remains a charge-response contribution. The importer
checks coordinates and energies against the supplied WAVECAR.

### Historical doubled output

To prepare an older result, preserve its original files and make a CSV with
these columns:

```text
k_index,band,kx_frac,ky_frac,kz_frac,energy_eV,min_gap_eV,omega_legacy_A2
```

An optional `spin` column selects channels. There must be exactly one row per
source k point and band for the selected channel; the CSV must contain actual
matching energies and gap information. Old plot tables without these fields
are not self-contained transport input, and the importer does not guess the
missing values. Use `--normalization legacy-double` to apply 0.5 exactly once.
For an already normalized file the column is `omega_z_A2` and the flag is
`--normalization physical`. A standard-normalization producer marker prevents
accidentally importing a new Fortran CSV as historical doubled data.

The normalized output retains `normalization_migration`, including
`input_convention`, `factor_applied` and producer comments. Small-gap entries
are marked invalid; the Hall command rejects invalid band curvature. This is
not a replacement for the producer's degeneracy or physical-operator treatment.

The correction does not supply nonlocal, PAW, SOC or Hubbard-projector terms
missing from a canonical-momentum approximation. It also does not change a
line-mode calculation into a two-dimensional integration mesh. See the
[Kubo method guide](KUBO_TRANSPORT.md).

## Choose the correct output kind

- Keep Berry flux from the Fukui method in radians with its cell and vertex-energy interpretation.
  Its old output is not an input to the point-curvature importer.
- Use the standardized point-curvature NPZ/JSON pair for the new `hall`
  command. Arrays and metadata travel together; the JSON binds the NPZ bytes
  by checksum.
- Carry the actual energy reference and occupations. A selected-band Chern number
  sum is not a replacement for a chemical-potential-dependent metal response.
- Preserve invalid or experimental status. Converting storage formats does
  not turn a failed matrix check into a validated physical calculation.

## Version and release records

`VERSION`, CLI version output and current citation metadata identify 1.6.5.
Historical changelog entries, the 1.2.0 release notes and archived reference
data retain their original versions. The CFF schema version is independent
of the software version.

The 2018 DOI `10.5281/zenodo.1402593` identifies VASPBERRY V1.0. It is not a
DOI for version 1.6.5. Use the immutable
[`v1.6.5` release](https://github.com/Infant83/VASPBERRY/releases/tag/v1.6.5)
for version-pinned source, or `master` for current source. Build binaries locally.
Cite the software version, exact commit and method references appropriate to
the calculation; see the [version policy](RELEASING.md).
