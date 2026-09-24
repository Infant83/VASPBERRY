# Local validation of the 1.3.0 development candidate

Initial checks on 2026-09-19; scientific examples and transport checks updated through 2026-09-24. This is a local source validation, not a hosted release
or a completed material-convergence study.

## Reproducible public checks

```bash
python3 -m pip install -r requirements-transport.txt
python3 -m unittest discover -s tests
make check-gnu
```

The initial Kubo/Hall validation run passed **224 tests**, including compiled
Fortran physics/parser checks. Compiler-dependent tests explicitly skip if
their compiler is unavailable; a skip is not compiler validation. GNU serial,
GNU/OpenMPI build/help/runtime and Intel Classic serial build checks passed
locally on that source. Hosted CI and Intel MPI remain separate checks.

Coverage includes existing Fukui/Z2/transport tests; independent two-band
curvature and QWZ sign/phase oracles; Hermiticity and missing vertices; full-mesh
geometry and oriented planes; occupations, partial bands and spin multiplicity;
schema types and hashes; one-time normalization migration; imported energy-gap
checks; regional/band sum rules; and small doping responses on large baselines.

The [method guide](KUBO_TRANSPORT.md) includes runnable public matrix → curvature
→ Hall commands. Those commands and the standalone
[CSV plotting example](../examples/kubo/README.md) passed local smoke checks.

## Actual VASP tutorials

The user-facing examples now start from the supplied **real MoS₂ and Bi VASP
WAVECAR files**. All six were recalculated with the current production
Fortran or Python WAVECAR path on 2026-09-19, with original outputs, figures,
checksums and numerical comparisons retained in each `reference/` directory.
The combined local run passed in approximately 76 seconds; timings depend on
hardware. The public Bi download was independently retrieved and verified
against its pinned 200,421,600-byte SHA-256 payload.

```bash
make serial
python3 examples/fetch_inputs.py bi --output-dir results/inputs/bi
python3 examples/run_examples.py \
  fukui-chern z2 hall-valley kubo-curvature circular-dichroism wavefunction \
  --bi-wavecar results/inputs/bi/WAVECAR --output-dir results/feature-examples
python3 -m unittest discover -s tests -v
```

| Actual input and calculation | Checked result |
|---|---|
| Bi 12×12 occupied bands 1:10, native Fukui | C = 0; sampled direct/global gaps 0.592448500 / 0.510044362 eV |
| Bi 12×12, native Z₂ | Z₂ = 1; two half-zone integer sums −3 and +3 |
| Bi occupied-subspace Fukui transport, T=0 | All 41 chemical potentials within the gap pass; max charge Hall residual 1.17e−16 e²/h |
| MoS₂ actual 48-point path, native Kubo | 96 band/point rows; 88 meet the stated 1e−5 eV isolation threshold; K/K′ opposite signs |
| MoS₂ actual path, native `-cd 2` | 9,648 k/energy rows; independent momentum/spectral sum agrees within 4.991e−5 a.u., within text-output rounding |
| MoS₂ Γ band 18, native `-wf` | Both spinor components; independent Fourier amplitude error 3.42e−7 Å^(−3/2), integrated pseudo-density 0.8258141586 versus coefficient norm 0.8258141556 |

The full local Python suite passed **238 tests** with no skips, including an
actual MoS₂ native Kubo batch calculation, immutable reference hashes and
archive contents, input-pointer rejection, malformed Fukui mesh/coordinate rejection, retained failures, and compiled
Fortran regressions. See the individual tutorials for comparison tolerances.
The CI `feature-examples` job downloads the actual Bi input, recalculates all
six examples and retains outputs/logs. Local success alone does not establish
hosted CI success.

An actual two-rank OpenMPI MoS₂ Kubo run also passed. Its 96-row high-precision
CSV is byte-identical to the serial reference, including the 88 valid / 8
unresolved-state mask. This check covers the native Kubo path; Python tutorial
wrappers and Python transport do not gain MPI support from it.

The [feature index](../examples/README.md) links the exact inputs, production
commands, reference outputs and [own-material guide](../examples/APPLY_TO_YOUR_SYSTEM.md).
Existing `examples/1H-MoS2/` and `examples/Bi_Z2/` inputs and historical outputs
remain unchanged. Bi is a zero charge-Hall sanity check, not a nonzero-Chern
or valley-Hall material demonstration. The MoS₂ path cannot supply a BZ
integral. Native canonical momentum does not include all PAW/nonlocal/SOC
velocity terms.

## Cartesian figures and full-zone MoS₂ example

The 2026-09-20 update adds a seventh real-material tutorial: occupied-subspace
Fukui curvature on a full 12×12 MoS₂ mesh. VASP 5.4.4 regenerated the 26-band
SOC WAVECAR from the public charge density and matching licensed potentials.
The fixed-charge run converged in eight electronic iterations (about 144 s,
peak resident memory 602 MB on the measured machine). Its minimum direct and
global gaps are approximately 1.674 eV.

Native VASPBERRY gives curvature extrema −12.3125 and +12.3124 Å²; integration
of the rounded text map gives C ≈ −5.0e−7. An independent Python overlap
calculation gives C ≈ 5.65e−16, with minimum link singular value 0.714 and
minimum plane-wave coverage 0.978. The historical MoS₂ map is retained as a
separate reference and is not claimed to be identical to the regenerated NSCF
calculation.

The final local suite passed **252 tests** with no skips. The first-BZ plotter
preserves native plaquette values and Cartesian lattice geometry. Geometry tests cover hexagonal and square cells, skew and rotated
bases, reversed orientation, area preservation, malformed meshes and lattice
mismatches. Existing Bi native outputs and the six earlier numerical reference
tables remain unchanged; the figure updates affect presentation only.

The seventh tutorial requires a locally generated VASP mesh. Public CI runs
the six tutorials with downloadable wavefunctions and tests the seventh
workflow's discovery/input checks. The [scientific report](TECHNICAL_REPORT.md)
and [MoS₂ tutorial](../examples/features/fukui-berry-curvature/) show the
resulting figures and the complete VASP preparation commands.

## Occupied-bundle Kubo and local single-band examples

The 2026-09-21 follow-up adds `-kubo_bundle 1`, which excludes internal
selected-band transitions before evaluating the Kubo sum. The final local
suite passes **287 tests**, including scalar/spinor bundle sums, internal
degeneracies, external-gap rejection, legacy normalization, periodic masked
interpolation and native-value preservation. GNU serial and MPI builds and
the help/runtime checks pass; a two-rank actual MoS₂ bundle CSV is identical
to serial. Existing individual-band Kubo CSV and all three DAT files remain
byte-identical when rerun on the prior 144-point input.

The occupied MoS₂ bundle has 144 valid mesh points and 49 valid path points,
with minimum external gap 1.673550 eV. At 18 shared mesh/path samples,
an independent complex128 WAVECAR momentum contraction agrees within
3.56e−14 Å² when using the native constants convention. The main figure
shows all occupied bands 1–18; the K/K′ values are approximately ∓13.178 Å².

A separate VASP calculation supplies two 9×9 valley patches for band 18.
Every one of the 162 points is isolated, with minimum separation 0.130203 eV.
The VASP run took 161 s and 653 MB peak resident memory. Independent EIGENVAL
gap checks, native rerun and figure regeneration pass. The new smooth maps
use NumPy bilinear interpolation for display only. The original numerical
references and un-smoothed figures remain available.

## Native Kubo pairs and occupation-weighted transport

The 2026-09-22 transport update passed the full **331-test** local suite,
including compiled native pair checks, independent pair-occupation formulas,
exact-degenerate-subspace rotations, finite-temperature integration, regional
partitions, orientation and explicit spin multiplicity. Native serial/MPI
help and two-rank runtime checks pass. Failure tests distinguish a native
zero exit code from the required completed-export footer.

On actual 12×12 MoS₂ data, the native pair export reproduces occupied-bundle
curvature within 2.67e−15 Å² for 26 bands; separate 40/60-band comparisons agree
within 4e−15 Å². An independent complex128 wavefunction contraction at five
representative points verifies all three numerator components for 40 and 60
bands. Existing per-band CSV/DAT and bundle CSV outputs remain byte-identical.
Serial and two-rank Hall tables agree, as do CSV, DAT and NPZ values. NPZ-only
output and subsequent cached rescans are tested; plots read all three formats.

Large native CSV files are imported in bounded row chunks with complete
coordinate, energy, index and row-coverage checks. On the measured 26-band
example, the same WAVECAR-to-Hall workflow took 3.22 s after this change versus
13.66 s before it; timings on the shared machine are illustrative. Cached
arrays and common-grid Hall values remain exactly equal.

The real-material study also distinguishes numerical pipeline checks from
input convergence. Occupied MoS₂ energies can be stable while the highest
empty states still violate time-reversal energy equality. Therefore a VASP
`EDIFF` completion alone is not used as evidence for an accurate virtual-state
sum. The [transport example](../examples/features/kubo-hall/) reports the
measured k/band-window dependence and source-state checks.

The standard WAVEDER adapter adds 15 checks covering occupied-bundle signs,
all Cartesian components, internal occupied-space rotations, selection-rule
zeros, producer restrictions, collinear filling, formats and guarded chunk
combination. Genuine chunks reproduce the same full-grid fixture; mixed
Hamiltonians, duplicate points and incomplete meshes are rejected. An actual
standard VASP MnBi₂Te₄ Γ/K/M optical probe agrees with an independent optical
contraction, and its three-point sampling is correctly rejected as a full-BZ
Hall input. That probe alone establishes no material Chern number or plateau.

The eighth feature, MoS₂ Kubo Hall, also passes an actual catalog-driven run
with all required numerical and figure outputs. Its baseline reference and
the seven catalog/orchestration checks pass.

## Actual MoS₂ transport refinement

The completed 2026-09-24 study separates mesh refinement, the retained pair
window and accuracy of the stored VASP states. At 300 K, on the common
chemical-potential grid, the K−K′ Hall-change curve differs by 12.890% for
12×12 → 24×24 and 1.482% for 24×24 → 36×36. Both comparisons use 60 stored
bands and pairs within bands 1–40. The last mesh increment does not meet the
stated 1% criterion.

On one 12×12, 96-band source, increasing the pair cutoff from 60 to 80 changes
that curve by 1.263%; 80 → 90 changes it by 0.202%. These are separate
refinement tests, not evidence of joint mesh/window convergence. At the same
12×12 mesh and cutoff 40, replacing the 60-band source with the 96-band source
changes the curve by at most 3.664e−9 e²/h. The final 36×36 total charge-Hall
residual is below 2.24e−9 e²/h, without imposed time-reversal averaging.
Increasing the explicitly requested energy-grouping tolerance from 1e−7 to
1e−6 eV changes the Hall difference by less than 2.2e−13 e²/h in the checked
final cases.

All eight refinement cases have exactly matching CSV, DAT and NPZ fields
(240 comparisons). The public references retain the VASP input conditions,
compressed original outputs, calculation metadata and two figures. Source
and output integrity checks, figure regeneration and relative-link checks
pass. The original 12×12, 26-band catalog reference remains separate.

The largest calculation used six fixed-charge 216-point jobs for the 36×36,
60-band mesh, followed by a byte-preserving WAVECAR assembly. The final warm
stages took 161 worker-minutes in total and at most 1.27 GiB resident memory
per job on the measured host; earlier seed/checkpoint stages are additional.
These measurements describe this run, rather than a general runtime promise.
See the [reference tables and conditions](../examples/features/kubo-hall/reference/).

## Actual magnetic Chern and optical checks

The three-septuple-layer MnBi₂Te₄ example uses a 21-atom film with alternating
out-of-plane Mn moments. Six completed, fixed-charge VASP optical runs form
the full 6×6 mesh, with 192 stored SOC bands and 123 occupied bands. Independent
checks compare every assembled coefficient record with its original source.
The calculation has a sampled direct/global gap of 17.101 meV.

The Python WAVECAR Fukui routine gives occupied C = −1, with maximum
plaquette phase 2.060 rad and minimum link singular value 0.305. Independently,
the PAW overlap matrices give C = −1. The standard WAVEDER adapter and an
independent optical contraction agree within 2.85e−14 e²/h on the full mesh.
Its coarse Hall result is +163.109 e²/h: the sharp Γ curvature is unresolved.
A constant T=0 response while the chemical potential remains inside the gap
only demonstrates unchanged occupations; it is not evidence of quantization.

The public preparation/runner audit passes source-state association,
Hamiltonian consistency, charge-density integrity and original-output checks.
It validates the actual VASP inputs and guarded workflow. The
[material example](../examples/materials/mnbi2te4-qah/) reports the separate
dense-integration reference and its validation scope.

The full-connection Wannier reference retains both Hamiltonian and position
matrices. At the seven near-Γ comparison points, its band-edge differences
from direct VASP are at most 1.26 meV and its curvature differences at most
1.44%. Across all 12 comparison points, the largest band-edge difference is
18.14 meV. These compare different finite subspaces: the 36 omitted deep
bands have zero total Chern number but can contribute local curvature.
The model uses 300 localization iterations and did not reach its requested
spread tolerance; it is not presented as a converged Wannier optimization.

An independent Hamiltonian-only scan of 160×160 points finds the minimum
direct and global gaps at Γ, both 17.101 meV. All three chemical potentials
spanning the central 90% of the DFT gap remain inside this sampled model gap.
This energy check took 271.8 s with batched Fourier transforms; it is separate
from the full-connection Hall calculation.

For portable reproduction, exact full operators and pair-dependent
translations are exported without dropping nonzero elements. The documented
Wannier90 3.1.0 effective-model reader needs one dimension-initialization fix;
its response formulas are unchanged. Independent operator, derivative and
curl comparisons agree within 2.2e−14. Local band/curvature values, a weighted
50-point test, and the entire 3,088-point integral reproduce the unmodified
original producer at its printed precision, including each J0/J1/J2 term.
The latter gives 1.0001934839919 e²/h before further quadrature refinement.
This is an operator-equivalence check, not a mesh-convergence estimate.

Nine independently compressed parts reproduce the complete numerical
operators. Public restoration checks all parts and complete operators before
creating an output directory; damaged input and existing output are rejected.
The patches, source version, GPL license and regeneration instructions are
provided alongside the numerical inputs. Licensed VASP code is not included.

Ten additional public workflow tests pass using synthetic process fixtures.
They cover incomplete Fortran termination, nonzero exits, missing connection
terms, wrong producer versions, stale or modified files, worker limits and
CSV/DAT/NPZ agreement. These test the external runner's safeguards; they are
separate from the 331-test full-suite run and the actual material calculations.
The public prepare/run/collect commands also pass a genuine 4×4, two-part
execution with the supplied full operators and patched reader. Execution and
collection took 23.6 s and 2.0 s, respectively, using one scientific worker.
Original exit statuses, completion messages, partition sums, units and all
three output formats agree. This deliberately coarse execution check is not
used as a quantization or convergence reference.

The final full-connection integral uses 27,600 weighted two-dimensional
points: an 80×80 base mesh with 9×9 subdivisions of cells within 0.18 Å⁻¹
of Γ. All eight producers finish normally with zero exit status. The result
is **1.0003849765685957 e²/h** at each of the three gap chemical potentials.
The preceding 60×60/refinement-11 calculation over the same region gives
1.0005919239860739 e²/h; their difference, **2.06947e−4 e²/h**, meets the
predeclared 1e−3 criterion for this finite model's quadrature.

The controls are not monotonic: expanding the refined radius from 0.12 to
0.18 Å⁻¹ at fixed 60×60/refinement-11 changes the result by 5.42162e−4 e²/h,
whereas increasing the inner subdivision from 7 to 11 at radius 0.12 Å⁻¹
changes it by only 5.03784e−6 e²/h. All completed controls are retained.
The final deviation from +1 is 3.84977e−4 e²/h; no symmetry average or integer
rounding is imposed. This verifies numerical integration at the stated
settings, separately from convergence of the VASP/Wannier model itself.

## Separate developer checks and historical wavefunction fix

The earlier QWZ, analytic optical, synthetic wavefunction and stored-field
checks are preserved under [validation/models](../validation/models/), outside
the user tutorials. All six still pass after the move (about 8.4 seconds
locally). They check signs, formulas, sum rules and regressions; they are not
actual VASP material demonstrations.

The synthetic wavefunction fixture previously exposed a phase-array extent
mismatch (`npmax` versus active `ncnt`) and an uninitialized POSCAR-header loop
flag. Both complete Fortran sources pass a regression compiled with
`-fcheck=all -finit-integer=1`, checking real/imaginary grids and headers.
The actual MoS₂ tutorial independently checks the retained amplitude convention.
No production numerical kernel changed during the real-input tutorial update.

## Historical normalization

The [migration guide](MIGRATION.md) and public import tests distinguish legacy
doubled Kubo data from already corrected data. Reproducibility targets the
corrected normalization; importing a historical doubled file applies the
documented conversion exactly once. Public material comparisons in this
record use the named reproducible examples above.

## Scope

These checks validate formulas, data contracts, migration and the checked
execution paths. They do not establish k-mesh/intermediate-band convergence
for an arbitrary material or the physical completeness of an exporter.
Experimental matrix inputs retain their labels, and the bare-momentum path
retains its missing-velocity-term limitations.

The WAVECAR Hall path computes intrinsic charge sheet response on uniform full
2D meshes. Invalid individual-band curvature is rejected. The full-connection
Wannier path described below additionally supports local cell refinement.
Spin/layer currents and general 3D bulk response need separate implementations.
Existing Fukui subspace workflows remain available for their documented scope.

## VASPBERRY full-connection Wannier backend

The complete local suite now contains **372 tests**, all executed and passing
with no skips (76.221 s). The GNU Fortran/Open MPI regression tests ran in this
suite. Production source hashes were unchanged before and after the run.

The 31 new tests comprise 12 effective-operator parser/cache checks, seven
independent physics checks and 12 command/workflow checks. Physical oracles
include analytical two-band curvature, Wilson plaquettes of explicit physical
wavefunctions with a nonzero basis connection, k-dependent integer orbital
gauge changes, exact internal degeneracies, Cartesian axis signs and length
scaling. Workflow tests cover actual CLI import, normalized output, whole-zone
cell partition, threaded/serial equality, source mutation, no-clobber, failed
partial output, memory preflight and band-path serialization.

The actual 138-orbital effective HH_R/AA_R input imports within 8.9e−16 of the
independently expanded source arrays. VASPBERRY's own 12-point energies and all
three curvature components agree with unmodified postw90 within each original
printed value's rounding interval. A 50-point weighted comparison also checks
all three J0/J1/J2 components. S/cm comparisons retain the reference producer's
CODATA-2006 conversion; it differs from the modern SI conductance convention
by a relative 3.668e−9. This conversion difference is not a kernel discrepancy.

The public `wannier-import` command was executed on the actual archived inputs,
and `wannier-bands` computed 3,601 path points with all 138 model bands.
The original postw90 output remains a separate independent reference.

The final public-command 27,600-point VASPBERRY calculation gives
**1.0003849781827807 e²/h** at all three gap chemical potentials. Its complete
quadrature coordinates, weights and parent cells are exactly equal to the
independent reference. The preceding wider-region 21,720-point result is
1.0005919297781898 e²/h, a final change of 2.06951595409e−4 e²/h.
The independently accumulated totals and all J0/J1/J2 components agree within
the original postw90 output rounding intervals. Modern/legacy physical constants
are accounted for before comparing the three components. The final command
completed in 576.64 s using four single-BLAS-thread workers while another
four-worker calculation ran; observed peak process RSS was 2.365 GB.

All five production conditions pass independent audits of command exit status,
operator and output hashes, complete integration partitions, all three curvature
components, J0/J1/J2 contributions and common CSV/DAT/NPZ table fields. The
small-region results are 1.0001934913592478 (3,088 points),
1.0000548041196191 (6,528 points), and 1.0000497725580064 (10,920 points) e²/h.
The wider-region controls above are retained because the convergence is not
monotonic with refinement radius. After accounting for the physical-constant
convention, the maximum scalar difference from postw90 is 6.73e−9 e²/h, within
the producer's printed precision. The complete five-condition queue took about
12.3 minutes; observed peak RSS remained below 2.37 GB per calculation.
