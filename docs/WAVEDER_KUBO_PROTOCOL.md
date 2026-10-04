# Standard Kubo protocol: WAVEDER input

From **1.6.5**, charge Kubo calculations use `--kubo-source waveder` by
default. A missing WAVEDER or unsupported request stops. To choose the
previous WAVECAR canonical-momentum approximation, pass
`--kubo-source wavecar` explicitly; the program warns and records that
choice. A WAVEDER file does not change an explicitly selected WAVECAR
operator. These optical files are additional inputs for standard charge
Kubo; unrelated WAVECAR-based topology tasks do not require WAVEDER.

## 1. Prepare and preserve one optical calculation

Converge the VASP electronic density, structure, magnetic state, cutoff and
PAW setup first. Generate the required eigenstates on the actual full BZ
mesh, or a path when only point curvature is required. The supported
producer is VASP **5.4.4** with the standard static longitudinal branch:

```text
LOPTICS = .TRUE.
LPEAD   = .FALSE.
LNABLA  = .FALSE.
NSW     = 0
ISYM    = -1
LREAL   = .FALSE.
```

Keep its completed `WAVEDER`, `WAVECAR`, `INCAR` and `OUTCAR` together, with
the structure, PAW identities, density, source/binary identity and exact
command. The adapter checks the implemented producer contract; other
versions, hybrid/meta-GGA and unsupported optical branches are rejected.
File acceptance is not a claim that every possible SOC/+U response term is
complete. Converge the chosen operator and electronic model independently.

## 2. Calculate point curvature with the native executable

For the complete insulating occupied space inferred from the source:

```bash
make serial
mkdir -p results/optical-curvature
build/vaspberry --task kubo --input-dir path/to/optics \
  --curvature-csv results/optical-curvature/KUBO.csv
```

`--input-dir DIR` selects the input directory. If omitted, it defaults to
the invocation working directory. The default files are `DIR/WAVECAR`,
`DIR/WAVEDER`, `DIR/INCAR` and `DIR/OUTCAR`. Override one with `--wavecar PATH`,
`--waveder PATH`, `--incar PATH` or `--outcar PATH`; **each override changes
only that file**. In particular, `--wavecar elsewhere/WAVECAR` does not move
the search for any other file. All files must still belong to the same
completed source run.

Relative CLI directory, per-file and output paths are relative to the
invocation working directory. `--input-dir` does not change the working
directory or prepend itself to explicit paths. For example, launching from
`/work/analysis` with `--input-dir ../optics --wavecar saved/WAVECAR` uses
`/work/analysis/saved/WAVECAR` and `/work/optics/{WAVEDER,INCAR,OUTCAR}`;
`--curvature-csv result.csv` still writes `/work/analysis/result.csv`.

Directory forms such as `job1`, `job1/optics`, `./job1`, `../job1` and
`/absolute/path/to/job1` are accepted. Bare `~/job1` is expanded by the shell.
For a quoted home path in the native command use `"$HOME/job1"`; native
Fortran does not expand a literal `"~/job1"`. Python CLI/INI paths expand
`~`. Quote paths containing spaces. These rules apply to the source directory
without changing the explicit-file or output-path rules above.

Omitting `--bands` retains the complete insulating occupied-space route.
An explicit selector instead requests a geometric band or bundle:

```bash
build/vaspberry --task kubo --input-dir path/to/optics \
  --bands 31,33:34 --curvature-csv results/optical-curvature/SELECTED.csv

build/vaspberry --task kubo --input-dir path/to/optics \
  --bands 31:32 --per-band 1 --curvature-csv results/optical-curvature/BANDS.csv
```

Band indices are one-based; comma-separated singletons and inclusive ranges
form the target set. The virtual index still runs over **all source NBANDS**.
Selecting `31:32` does not restrict the intermediate sum to those two bands.
The default is the trace of the selected space; `--per-band 1` asks for each
band separately and rejects an unresolved individual state. A trace may
contain a complete producer-degenerate group because equal-weight internal
pairs cancel. Splitting such a group is rejected. The checked group includes
transitive connections within the producer's **2 meV** threshold.

The default output name is `KUBO_WAVEDER.csv`. It contains source k points,
the Cartesian `omega_z_A2` trace curvature in Å² and the minimum external
gap, with selected-band, operator/source/completion metadata. It does not
write canonical per-band DAT companions. A path gives only point curvature.
For integration, supply the actual complete uniform mesh:

```bash
build/vaspberry --task kubo-integral --input-dir path/to/optics \
  --bands 31:32 --mesh 12,12 --curvature-csv results/optical-curvature/MESH.csv
```

Never substitute a path or irreducible mesh. An explicit selected-space
integral uses unit geometric band weights, without physical spin multiplicity,
and has no total-Hall label. Even if the explicit set matches the occupied
bands, it remains a selected geometric result. The omitted-band, validated
occupied-input route retains its physical T=0 Hall integral. Neither finite-grid
integral is rounded to a Fukui Chern invariant.

The current native WAVEDER contraction runs on MPI rank 0; other ranks
wait, so adding MPI ranks does not accelerate this route. Its reader accepts
at most 32 million complex matrix elements (about 256 MB for the stored
complex64 matrix) and rejects larger inputs before allocation. This is an
implementation memory bound, not a physical convergence limit. For larger
meshes, generate matching fixed-charge optical chunks and use the Python
Hall route below to validate and combine the complete mesh. Do not discard
source bands or k points merely to bypass the native bound.

## 3. Apply occupations and calculate a selected charge contribution

Python `kubo-hall` and `waveder-hall` accept either `--bands SELECTOR` or
`--occupied N`, never both. For selected bands, the chemical potential and
temperature determine their weights on unchanged source eigenstates. For
example, with mesh and energies chosen from the actual source:

```bash
python3 tools/vaspberry_kubo.py kubo-hall \
  --input-dir path/to/optics --bands 31:32 --mesh "$NX" "$NY" \
  --energy-reference 'unchanged VASP eigenvalue zero' \
  --mu-min "$MU_LO" --mu-max "$MU_HI" --mu-num 3 \
  --mu-reference "$MU_REF" --temperatures 0 300 \
  --output-dir results/optical-selected-hall
```

For this WAVEDER CLI command, `conductivity.csv`, `conductivity.json` and
requested companion formats are written directly in `--output-dir`.
For example, plot its selected contribution with:

```bash
python3 tools/plot_hall.py results/optical-selected-hall/conductivity.csv \
  --output-dir results/optical-selected-hall-plot
```

The INI front end shown below instead writes tables under its `output/hall/`
subdirectory and its `plot` subcommand locates them automatically.

The Python `kubo-hall` command follows the same input-directory and per-file
override rules as the native command: omitted `--input-dir` means the
invocation working directory, and `--wavecar` changes only WAVECAR. It also
accepts `--waveder`, `--incar` and `--outcar`. Explicit relative paths remain
relative to the invocation working directory; output paths are unchanged.

The compatibility `--run-dir DIR [DIR ...]` alternative names one optical
run or matching chunks; do not combine it with `--input-dir`. One run permits
per-file overrides. Several runs require conventional filenames in each
directory and reject individual overrides.

Do not pass `--binary` or MPI launch controls
to this Python WAVEDER backend; those options are for the explicit WAVECAR
route. For collinear ISPIN=2 input, Python `kubo-hall` selects one channel
with `--spin` (default 1): calculate both and sum for the two-channel
contribution. A selected sum still is not certified total AHC. The native
WAVEDER mesh integral sums both stored spin channels, using unit geometric
weights for an explicit selection and physical multiplicity only for its
omitted-band occupied route.

The dedicated `waveder-hall` command remains available, including validated
multi-run mesh assembly. `--occupied N` keeps the legacy insulating T=0
contract: μ and reference must lie inside the common gap. `--bands` allows
μ/T scans only when every required pair is available and its weight can be
resolved. A small omitted thermal weight is not silently treated as zero. The [full input contract](KUBO_TRANSPORT.md#standard-waveder-insulating-paw-hall-response)
explains its source-identity checks. Numerical tables use sheet conductance
in e²/h from reciprocal area: no guessed slab thickness or extra simulation
cell-height multiplication. A selected subset gives that subset's charge
contribution; finite temperature does not turn it into the total AHC. Record
the selection and scope beside every σ or Δσ plot. An insulating full-occupied
T=0 gap scan is flat because occupations are unchanged; this alone does not
establish convergence or quantization.

The single-file postprocessor uses the same default:

```ini
[run]
kubo_source = waveder
input_dir = /absolute/path/to/optics
output = /absolute/path/to/new-results
mesh = 12 12
energy_reference = unchanged VASP eigenvalue zero

[hall]
bands = 31:32
mu = 1.10 1.20 3
reference = 1.15
temperatures = 0 300
```

These band indices, mesh and energies are illustrative, not a material
reference. Replace them with the values in your source run. Explicit relative
INI paths, including `input_dir` and per-file overrides, are relative to the
INI file. An omitted `input_dir` instead means the invocation working
directory; it is not inferred from the INI location or a `wavecar` override.
The optional `[run] wavecar`, `waveder`, `incar` and `outcar` keys override
only their respective default file. The legacy `optical_run_dir` key remains
a WAVEDER-only alias for `input_dir`; do not combine both keys.
Run `vaspberry_post.py check analysis.ini`, then
`run analysis.ini` and `plot /absolute/path/to/new-results`. A compiled
binary is required only for the explicit WAVECAR route. A WAVEDER INI must
omit `[run] binary`, use `mpi_procs = 1` if supplied, and omit any custom
MPI launcher. An explicit `optical_run_dir` is rejected with the WAVECAR
approximation. For collinear input,
choose INI `spin_mode = collinear-up` or `collinear-down` explicitly and sum
the two separately computed channel contributions. Remove `[hall] occupied`
and `bands` when choosing the WAVECAR route; these selectors are WAVEDER-only. The
WAVEDER `[hall] bands` selector is mutually exclusive with `occupied`; it
uses the same single/range/list syntax as the CLI.

## 4. Why a rectangular WAVEDER can suffice

The file stores `NBANDS` bra rows by `NDBANDS` derivative-ket columns for
each k, spin and Cartesian direction. `NBANDS=80` does not require the
producer to output 80 derivative columns. The audited 5.4.4 branch can use
`NDBANDS=min(2*LAST_FILLED_OPTICS,NBANDS)`: an occupied bound of 30 can give
**80×60**, with every stored bra 1–80 connected to ket columns 1–60.
`NDBANDS` is not the intermediate-state summation cutoff.

For the complete occupied trace, occupied–occupied pair terms cancel.
The empty-bra/occupied-ket block is sufficient:

```math
\Omega^{\mathrm{occ}}_{xy}=-2\,\mathrm{Im}
\sum_{n\in\mathrm{occ}}\sum_{m\in\mathrm{empty}}
C_{mn,x}^{*}C_{mn,y}.
```

For a general selected response, define `g_n = s_n f(E_n,μ,T)`, where
`s_n` is one for selected target bands and zero otherwise. The pair sum uses
`g_n-g_m`. Native geometric traces instead use `g_n=s_n`, without Fermi
weights. Equal-weight internal pairs cancel; every nonzero-weight pair must
be represented by at least one stored orientation. A valid reverse
orientation supplies the same pair with the appropriate sign, without
constructing a square matrix.

For an 80×60 file, a target band ≤60 can connect to all 80 source bras.
A missing 61–80 internal pair can still matter for a different selection or
μ/T scan; if its weights differ, the calculation rejects it. The code does
not zero-fill that block, restrict the virtual sum to the selected bands, or
assume high conduction weights vanish at finite temperature.

The producer zeroes couplings within its 2 meV degeneracy treatment. A
complete group with equal response weights needs no internal pair and is
allowed. Unequal weights within that erased group are unresolved and are
rejected, including at finite temperature; temperature does not recover the
lost matrix elements. Required numerical zeros outside such groups remain
valid zeros. These coverage and weight checks determine applicability for
each request; a rectangular shape alone neither grants nor denies support.

The PAW projector, overlap and augmentation terms enter **inside VASP's
optical matrix construction, before writing WAVEDER**. VASPBERRY consumes
those values without adding PAW corrections again. The connection already
contains the energy denominator; the formula above does not divide by the
gap a second time. Retain the original complex precision and producer cutoff
in the provenance.

## 5. Explicit WAVECAR approximation and unsupported requests

| Request | Standard WAVEDER route | Explicit WAVECAR route |
|---|---|---|
| Complete gapped occupied trace | Omit native `--bands`, or Python `--occupied N` for T=0 Hall | Canonical-momentum approximation |
| Selected single band, range or list | Geometric trace; optional native `--per-band 1` if resolvable | Available with canonical isolation rules |
| Selected μ/T contribution | Python `--bands`; all nonzero-weight pairs and erased-group weights checked | Approximate pair calculation |
| Native all-pair cache export | Not provided by the WAVEDER route | `kubo-pairs --kubo-source wavecar` |
| Projected-spin Kubo proxy | Unsupported | `spin-kubo --kubo-source wavecar`; not conventional spin Hall |
| Conventional spin-current response | Requires separately validated spin/current operators | Not supplied by ordinary WAVECAR momentum |

Reproduce a saved canonical example with:

```bash
build/vaspberry --task kubo --kubo-source wavecar \
  --wavecar examples/1H-MoS2/KPATH/2.band/WAVECAR \
  --bands 1:18 --curvature-csv results/canonical.csv

build/vaspberry --task kubo-pairs --kubo-source wavecar \
  --wavecar path/to/WAVECAR --pairs-csv results/PAIRS.csv
```

The fallback uses pseudo-wavefunction canonical momentum, omitting the
longitudinal PAW augmentation/nonlocal velocity terms in the standard optical
route. Denser k meshes and more intermediate bands converge that approximation;
they cannot add its missing operator terms. A missing or invalid WAVEDER never
authorizes fallback. Resolve the input issue or select this approximation
explicitly. For requests whose required matrix pairs or spin/current
operators are unavailable, use a separately validated input route in
[OPERATOR_ROUTES.md](OPERATOR_ROUTES.md).

## 6. Reusable example and convergence

The [selected-band example](../examples/features/waveder-selected/) provides
a clearly synthetic format/arithmetic fixture plus reusable source commands.
It distinguishes native geometric traces from occupation-weighted contributions
and retains expected rejection cases. It is not a material benchmark.

Compare k meshes at fixed electronic model and operator. Separately increase
accurately computed source NBANDS and test the virtual-state sum. Preserve
source occupations, energies, operator identity, included bands, thresholds,
μ/T, reciprocal area and all input/output hashes. Compare operators only on
matching source states and sampling; convergence in one route does not prove
agreement with another operator.

The [MnBi₂Te₄ example](../examples/materials/mnbi2te4-qah/) provides actual
standard optical inputs and the unrounded 6×6 Hall diagnostic. Its large
coarse-grid value is retained to demonstrate the need for mesh convergence;
this protocol update does not relabel it as a converged material prediction.
Historical MoS₂ and Bi WAVECAR references remain explicitly canonical-momentum
results with their original producer provenance.
