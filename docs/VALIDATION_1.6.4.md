# VASPBERRY 1.6.4 validation scope

**Historical candidate record.** Version 1.6.4 was prepared and validated
locally but never published. It is superseded by [1.6.5](VALIDATION_1.6.5.md).
The observations and counts below retain their original candidate scope.

The implementation checks below completed on 2026-10-02 for the
WAVEDER-default source update, before its version labels were corrected
to the maintainer-selected 1.6.4. They ran on macOS with GNU Fortran and
Open MPI; their original run paths and evidence remain unchanged. The
version-correction checks are recorded separately. This is a prepared
source version; these results do not assert a published tag or remote
CI status.

| Check | Result |
|---|---|
| Full `python -m unittest discover -s tests -v` | **909 passed**, no skipped tests |
| `make check-gnu BUILD_DIR=build-waveder170` | Serial/MPI builds, help and two-rank runtime checks passed |
| Native WAVEDER executable checks | 12 passed; independent rectangular matrix contraction, source validation, scope barriers and MPI equality |
| Native/Python actual standard optical input | Six MnBi₂Te₄ k points, 192 stored bands, 123 occupied; maximum curvature difference 5.55e-16 Å² |
| Generic default `kubo-hall` actual input | Six optical chunks forming a 6×6 mesh; conductivity CSV byte-identical to the earlier validated WAVEDER calculation |
| Explicit WAVECAR example | Actual MoS₂ 48-point path completed with warning and preserved canonical operator metadata |
| Plot/reader integration | PAW trace accepted with its own operator; mixed PAW/canonical panels rejected; 32 plotting checks passed |
| Native include dependency | Rebuild verified after WAVEDER include content changes, including unchanged timestamps |
| Technical report PDF | 34 pages, all 16 figure captions, source/figure checksums and visual rendering checked |

The optical CSV SHA256 is
`4a87786ae59e5a8326df1e5190f18f908a86c8406ba09406677d507429457710`.
The actual optical comparison is a software/operator reproduction. Its
coarse 6×6 material result remains unconverged. The source contains Hubbard U;
this check does not newly certify complete SOC/+U velocity terms.

The first full test pass exposed two historical parser fixtures that still
assumed implicit WAVECAR Kubo. They were migrated to explicit source selection,
and the final full suite passed. Separate retained failures include a real
OUTCAR plane-wave-report parsing case, corrected before the final native test,
and a plotting test's expected metadata label. Independent review also added
ordered final-OUTCAR-table validation and rejected unused cross-route options.

Detailed local commands, stdout/stderr, exit status, source/input hashes and
failed attempts are retained in the research workspace at
`runs/kubo-waveder-default-20261002/`, with separate native, Python, documentation
and integration records. The public test suite contains the independent
synthetic fixtures and barriers used here.

Intel/Intel-MPI builds, Linux and remote GitHub CI were not executed for this
update. Existing GNU legacy MPI argument-mismatch warnings remain. No new VASP
calculation, SH.LEE patch, licensed source redistribution or material
convergence claim is part of this version change. The new native WAVEDER
contraction runs on rank 0 and bounds its matrix at 32 million complex
elements; the Python optical route and chunked input retain their own limits.

## Version correction validation

The maintainer subsequently fixed the intended number at **1.6.4**. Only
version strings, release-document names/links and corresponding test expectations
changed; calculation logic and saved scientific results were preserved.
For this final version, 14 version/compatibility tests and 65 source/help/Kubo
workflow tests passed. The default `build/vaspberry` and `build/vaspberry-mpi`
were rebuilt; `make check-gnu` passed, including serial/MPI help and two-rank
runtime checks. The report PDF was regenerated and its two changed pages
visually checked; the other 32 rendered pages matched the preceding PDF.

These correction checks are retained separately in
`runs/version-correction-164-20261002/`. The initial version test ran before
the concurrent changelog edit completed and failed its changelog assertion;
the completed-tree rerun passed all 14 tests. The earlier 909-test full-suite
record remains preserved as implementation evidence, rather than relabeled
as a fresh full-suite run of the corrected version.

## Final code review and release preparation

On 2026-10-02, the corrected 1.6.4 tree passed a fresh full suite of **945
tests, with no skips**. Independent review found two executable callers
that still assumed implicit WAVECAR input: the native-spin CI driver and
the full-velocity/canonical comparison example. Both now explicitly select
`--kubo-source wavecar`. Their pre-fix failures were retained; their actual
post-fix executions passed. A compiled regression covers the example caller.
No additional blocking defect was found within the reviewed native/Python
WAVEDER paths; physical operator completeness remains outside this verdict.

The actual MoS₂ serial/two-rank MPI validation and the complete synthetic
spin-Chern/spin-Kubo CI driver both passed. The MoS₂ pair table contained
23,808 rows with zero serial/MPI numerical difference. These are local GNU
checks; the historical name of `run_intel_mpi_validation.py` does not mean
Intel compilers were executed locally.

The prepared manual-only 1.6.4 publisher passed 35 publication-contract tests,
workflow lint and an independent offline metadata/release-body review. It
requires all thirteen remote validation jobs on the exact publication commit
and preserves the verified v1.6.3 commit. Publication was not dispatched.
The source-candidate archive and clean-install receipts are retained with
`runs/followup-164-20261002/`; the candidate is an uncommitted-tree snapshot,
not a claim that remote CI or public release validation has completed.


## Explicit input directory follow-up

On 2026-10-04, the prepared 1.6.4 interface changed to `--input-dir DIR`,
defaulting to the invocation directory. Explicit file paths change only that
file; WAVECAR no longer selects the auxiliary input directory. Native, Python
and INI routing, help, examples and the technical report use this rule.

The completed tree passed a fresh full suite of **966 tests with no skips**
and GNU serial/MPI builds, help and runtime checks. Routing regressions cover
cwd defaults, independent overrides, option order, missing inputs, directory
validation and Python chunk compatibility. Actual six-point MnBi2Te4 optical
input retained the previous native values, with identical serial/MPI results
and a maximum native/Python difference of 5.55e-16 A². This is input routing
and implementation verification, not a new material-convergence calculation.

Review replaced a GNU-specific working-directory call with POSIX C binding
to retain the existing GNU/Intel portability approach. Intel execution remains
unverified locally. Initial path-display assertions assumed macOS `/var` and
`/private/var` were different; those test expectations were corrected. New
Python tests avoid the Python 3.11-only `contextlib.chdir` helper so they remain
compatible with the Python 3.10 CI configuration.

The report PDF has 34 pages and all 16 figure captions. The input-directory
page was visually checked; the other 33 rendered pages match the preceding
reviewed PDF exactly. Detailed receipts and failed attempts are retained in
`runs/kubo-input-directory-20261004/` in the research workspace. No remote CI,
commit, tag or publication is implied by these local checks.

## Selected WAVEDER bands and occupation-weighted response

The subsequent 2026-10-04 implementation removes the occupied-only restriction
for explicit band selections. Native single/range/list selectors produce
geometric traces or isolated-band rows. Python CLI and INI paths evaluate
selected Fermi-weighted μ/T contributions with all source NBANDS retained as
virtual states. Required-pair coverage, full-source producer clusters and
equal-weight cancellation determine whether a request is supported.

The final full suite passed **998 tests with no skips**. GNU serial/MPI
builds, help and runtime checks passed. Independent exhaustive checks covered
all 127 nonempty subsets of a seven-band Python fixture and 180 compiled-native
CLI cases spanning rectangular storage and scalar, spinor and collinear
inputs. The 3,207 accepted Python contractions and 134 accepted native outputs
agreed with independent direct sums; 984 Python and 46 native missing-coverage
cases rejected as expected. Exact degenerate-group rotations, unequal
occupations and nonzero finite-temperature tails were checked separately.

On actual standard 80×60 optical data, public native n31, n32 and 31:32 trace
outputs agreed with direct columns using all 80 virtual states to within
5.33e-15 Å². On existing standard MnBi₂Te₄ 6×6 data, bands 123:124 with M192
gave 36 μ/T points at 0, 100 and 300 K within 7.47e-14 e²/h of independent
`sum_n f_n Omega_n` integration. Unified and standalone CSVs were identical.
These checks do not certify material convergence, full physical SOC/+U
velocity completeness or total AHC from a selected contribution.

The [synthetic executable example](../examples/features/waveder-selected/README.md)
checks INI validation, calculation, plot labels, native trace/per-band output
and required rejections. All synthetic inputs are explicitly labeled.
Implementation hashes and version now accompany the direct WAVEDER CLI
outputs; this provenance follow-up passed its focused 13-test suite.

Review found and corrected a repeated-selector compatibility regression and
semicolon-containing SYSTEM title handling. The initial full suite's one
selector failure and all failed harness attempts remain distinct from final
PASS records under `runs/kubo-selected-bundle-implementation-20261004/` in the
research workspace. Version 1.6.4 remains a locally prepared release candidate.
