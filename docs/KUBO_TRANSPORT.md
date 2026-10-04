<a id="kubo-curvature-and-two-dimensional-charge-hall-transport"></a>

# Kubo-formula Berry curvature and two-dimensional charge Hall response

Actual VASP-based tutorials and reference results are provided for
[native MoS₂ Kubo-formula Berry curvature](../examples/features/kubo-curvature/) and
[MoS₂ Kubo-formula charge-Hall scans](../examples/features/kubo-hall/), and
[Bi occupied-subspace Hall response](../examples/features/hall-valley/). Their
methods and required sampling are explicitly different: the Bi example uses
Berry flux from the Fukui–Hatsugai–Suzuki (FHS) link-variable method, also
called the Fukui method, whereas the MoS₂ examples use Kubo-formula point
Berry curvature.

From 1.6.5, the standard charge Kubo source is same-run WAVEDER plus
WAVECAR/INCAR/OUTCAR. Read the [standard protocol](WAVEDER_KUBO_PROTOCOL.md)
first: native `--task kubo --bands SELECTOR` writes selected geometric
curvature, and `kubo-hall --bands SELECTOR` integrates the selected μ/T
charge contribution when all required pairs are valid. Omitting native
`--bands`, or choosing Python `--occupied N`, retains occupied T=0 behavior. Missing
or unsupported WAVEDER input stops rather than changing the operator.

The explicit `--kubo-source wavecar` approximation reads WAVECAR and
exports point curvature or interband pair numerators. The supplied `tools/vaspberry_kubo.py` tool
postprocesses those results: it validates a reusable cache, applies occupations
and integrates the two-dimensional intrinsic charge Hall response. Its
additional matrix import and analytic-check commands are advanced interfaces.
Python 3.10+ and `requirements-transport.txt` are required for this stage.

Use `python tools/vaspberry_kubo.py --help` and each subcommand's `--help` for
the available options. The [output specification](OUTPUT_FORMAT.md) describes
the interchange files; [migration](MIGRATION.md) covers older normalization.
The [operator input guide](OPERATOR_ROUTES.md) compares the ordinary WAVECAR
approximation with the standard optical and optional full-matrix routes, including
which routes require VASP source instrumentation.

## Start from actual VASP output

Use the [MoS₂ tutorial](../examples/features/kubo-curvature/) to execute
VASPBERRY for a Kubo-formula Berry-curvature calculation on real VASP WAVECARs generated with the linked
preparation recipes. It includes the actual
commands, high-precision CSV files, occupied-bundle maps and a single-band
valley example. The full mesh and band path are calculated separately; a path
is not an integration mesh. The native operator is
canonical momentum; PAW/nonlocal/SOC velocity corrections are not supplied by
that approximation.

The [Bi charge-Hall tutorial](../examples/features/hall-valley/) uses a different
production route: occupied-subspace Hall response in the insulating gap,
using plaquette flux calculated by the Fukui method from a full-mesh WAVECAR.
Bi's unresolved Kramers pairs prevent treating its individual
bands as isolated-band input for pointwise Kubo-formula Berry curvature. Its zero charge Hall plateau is an actual
material sanity check, not a nonzero valley-Hall demonstration.

<a id="wavecar-to-charge-hall-in-one-command"></a>

## Explicit WAVECAR approximation: charge Hall in one command

The standard entry is `kubo-hall` with WAVEDER. To reproduce the historical
canonical-momentum approximation, `wavecar-hall --kubo-source wavecar`
combines the native export, cache
validation and numerical Hall integration. **Python launches the
VASPBERRY executable**, then reads its saved interband pairs. It performs no
VASP calculation. Plotting remains a separate command so the same numerical
results can be displayed in several ways.

WAVECAR-based commands detect the component count and spin multiplicity.
Use `--spin 1` or `--spin 2` to choose a collinear channel; explicit
`--spinor-components` values check the detected layout.

First choose the executable and native rank count. From the repository root,
with the Python dependencies installed, a serial GNU build uses:

```bash
make serial
binary=build/vaspberry
ranks=1
```

For Intel MPI, load the site's oneAPI compiler, Intel MPI and oneMKL environment
as described in [BUILD.md](BUILD.md#intel-oneapi-build), then use:

```bash
make ifx-mpi
binary=build/vaspberry-ifx-mpi
ranks=4
```

`mpiexec` must belong to the MPI implementation linked into the executable.
The GNU MPI alternative is `make mpi`, `binary=build/vaspberry-mpi`, and the
desired `ranks`; use its matching launcher. These build choices do not change
the Python scan or plotting commands. The launcher is used only when
`--mpi-procs` exceeds one; do not launch the Python script with `mpiexec`.

After preparing the actual MoS₂ input described below, run:

```bash
python3 tools/vaspberry_kubo.py wavecar-hall --kubo-source wavecar \
  --binary "$binary" --mpi-procs "$ranks" --mpi-launcher mpiexec \
  --wavecar results/mos2-24-b60-vasp/WAVECAR \
  --mesh 24 24 \
  --energy-reference 'unchanged VASP eigenvalue zero' --pair-band-max 40 \
  --mu-min -1.47487388 --mu-max -1.17487388 --mu-num 121 \
  --mu-reference -0.43809870 --temperatures 0 300 \
  --regions examples/features/kubo-hall/regions.json --difference valley:K:Kprime \
  --degeneracy-policy coalesce --degeneracy-threshold-eV 1e-7 \
  --output-dir results/mos2-hall-wrapped

python3 tools/plot_hall.py results/mos2-hall-wrapped/hall/conductivity.csv \
  --regions total K Kprime valley --temperatures 300 --quantity delta-sigma \
  --energy-origin-eV -1.27487388 --energy-label 'mu - Ev (eV)' \
  --output-dir results/mos2-hall-wrapped-plot
```

Use a new `--output-dir` that does not already exist. The energy range, reference,
regions and explicit numerical-degeneracy approximation above belong to this
MoS₂ example; replace them for another material.

| File under `results/mos2-hall-wrapped/` | What it contains | Next use |
|---|---|---|
| `native/PAIRS.csv` | VASPBERRY k/band-pair energies, coordinates and undivided Cartesian matrix-product numerators | Import or inspect the raw electronic response data |
| `native/stdout.log`, `native/stderr.log` | Native calculation diagnostics | Check that the export completed |
| `pairs/pairs.npz`, `pairs/pairs.json` | Validated full source pair cache, energies, mesh, lattice and provenance | Reuse for a different μ/T/region/window scan or PROCAR attribution |
| `hall/conductivity.csv`, `.dat`, `.npz` | Equivalent tables of absolute σ, reference-subtracted Δσ, carrier counts, μ, T and region | Plot in VASPBERRY, NumPy or another analysis tool |
| `hall/conductivity.json` | Units, operator, retained/source bands, region definitions and diagnostics | Interpret the numerical table |
| `workflow.json` | Native command, hashes, completion status and failure information | Inspect the complete run |

The separate plot directory contains `hall.png`, `hall.pdf`, `hall.svg` and
`plot.json`. This example draws the **300 K Δσ versus μ−Ev** curves for total,
K, K′ and K−K′. It does not draw a band structure or a k-resolved curvature map;
those use the [native curvature examples](../examples/features/kubo-curvature/).

For the next scan, reuse the pair cache directly. For example, retain the
same μ range and bands while calculating 100 K only:

```bash
python3 tools/vaspberry_kubo.py pair-hall \
  --pairs-dir results/mos2-hall-wrapped/pairs --pair-band-max 40 \
  --mu-min -1.47487388 --mu-max -1.17487388 --mu-num 121 \
  --mu-reference -0.43809870 --temperatures 100 \
  --regions examples/features/kubo-hall/regions.json --difference valley:K:Kprime \
  --degeneracy-policy coalesce --degeneracy-threshold-eV 1e-7 \
  --output-dir results/mos2-hall-100K
```

This command reads `pairs/`; it does not execute VASPBERRY or reread WAVECAR. Its
`conductivity.*` files are directly under `results/mos2-hall-100K/`; pass that
CSV, DAT or NPZ to the same plotter. For the explicit native-first procedure,
see [Native pairs to charge Hall](#native-pairs-to-charge-hall) below.

## Native pairs to charge Hall

This explicit WAVECAR approximation has three stages: execute VASPBERRY to export pair
data, perform numerical postprocessing, and plot the results. First generate the 24×24, 60-band MoS₂ WAVECAR with the
[material tutorial](../examples/features/kubo-hall/#1-prepare-and-run-vasp).
Run from the repository root with a fresh result directory:

```bash
mkdir -p results/mos2-hall/native
build/vaspberry --task kubo-pairs --kubo-source wavecar \
  --wavecar results/mos2-24-b60-vasp/WAVECAR \
  --pairs-csv results/mos2-hall/native/PAIRS.csv

python3 tools/vaspberry_kubo.py import-pairs \
  --csv results/mos2-hall/native/PAIRS.csv \
  --wavecar results/mos2-24-b60-vasp/WAVECAR \
  --mesh 24 24 \
  --energy-reference 'unchanged VASP eigenvalue zero' \
  --output-dir results/mos2-hall/pairs

python3 tools/vaspberry_kubo.py pair-hall \
  --pairs-dir results/mos2-hall/pairs --pair-band-max 40 \
  --mu-min -1.47487388 --mu-max -1.17487388 --mu-num 121 \
  --mu-reference -0.43809870 --temperatures 0 300 \
  --degeneracy-policy coalesce --degeneracy-threshold-eV 1e-7 \
  --formats csv dat npz --output-dir results/mos2-hall/hall

python3 tools/plot_hall.py results/mos2-hall/hall/conductivity.npz \
  --quantity delta-sigma --formats png pdf svg \
  --output-dir results/mos2-hall-plot
```

These energies belong to this MoS₂ reference. Choose the range and reference
from your own bands. The [complete example](../examples/features/kubo-hall/)
adds periodic K/K′ regions, their difference, band plots and convergence
comparisons. The default degeneracy policy is `error`; the explicit `coalesce`
choice above is explained below.

The pair export includes every stored band pair and all three Cartesian
components. Omit band selectors in pair mode. The importer checks coordinates,
energies, lattice and complete selected-spin pair coverage against WAVECAR.
For a collinear two-channel calculation, import and integrate each spin with
multiplicity one, then sum the two charge responses.

`pairs/pairs.npz` and `pairs.json` form a reusable cache. To change temperature,
chemical potential or region definitions, repeat only `pair-hall` with that
cache and a new output directory. The postprocessor writes
`conductivity.csv`, `.dat`, `.npz` and `.json`; `--formats npz` alone is also
valid. JSON records units, operator, occupations, regions and diagnostics.
CSV and DAT are readable tables; NPZ stores typed arrays. The plotter reads
any of the three numerical formats without a conversion step.

The numerical tools do not encode MoS₂, its occupied count, its valley
coordinates or its energy zero. Those choices belong to the input data and
the example commands. Keep the same cache when changing the following:

| Question | Change only these `pair-hall` inputs |
|---|---|
| Response versus chemical potential or temperature | `--mu-min`, `--mu-max`, `--mu-num`, `--mu-reference`, `--temperatures` |
| Different valleys, pockets or sampled patches | `--regions`; optionally `--difference NAME:LEFT:RIGHT` |
| Retained pair-band window on the same eigenstates (both endpoints) | `--pair-band-max`, keeping the complete occupied window and unresolved degenerate groups; this tests intermediate-state convergence only while the occupied contribution stays complete |
| Smaller temporary memory use | `--mu-chunk`; this does not change the physical model |

Changing the Hamiltonian, structure, density, k mesh or available source
bands requires new VASP states and a new export. Changing the plot style,
axis range or selected curves requires only the completed conductivity
table. It cannot produce an uncalculated temperature or region; rescan the
pair cache for those. See the [independent plotting recipe](OUTPUT_FORMAT.md#read-and-plot-with-other-tools)
for reading numerical results without importing any VASPBERRY modules.

For MPI, replace the native command with:

```bash
mpiexec -n 4 build/vaspberry-mpi --task kubo-pairs --kubo-source wavecar \
  --wavecar results/mos2-24-b60-vasp/WAVECAR \
  --pairs-csv results/mos2-hall/native/PAIRS.csv
```

Native k-point work is distributed across ranks. Python integration then uses
NumPy and bounded chemical-potential chunks. Coefficient caching is capped at
64 MiB per native rank, with a direct-read fallback. The complete pair cache
still scales as the number of k points times the square of `NBANDS`.

### Intermediate-state convergence

Generate a sufficiently large, well-converged WAVECAR once and rescan its
cache with `--pair-band-max M`. This restricts both ends of each pair to bands
1:M while retaining the source `NBANDS` in the metadata. A reduced window is
labeled `truncated_pair_space`; using all stored bands is `all_source_bands`.
The cutoff must leave the occupied window complete and must not split an
unresolved degenerate group. Comparing M values on identical eigenstates
separates truncation of the virtual-state sum from changes in the VASP solver.
Check the highest empty eigenstates themselves: convergence of the occupied
total energy alone does not guarantee their accuracy. For nonmagnetic MoS₂,
the equality of energies at k and −k provides a useful additional check.

### Layer, orbital and spin character from PROCAR

The [PROCAR example](../examples/features/procar-character/) adds atom/layer,
orbital and chosen-axis spin attribution to the same native pair workflow.
It uses the matching SOC `PROCAR`, `WAVECAR` and `OUTCAR`, user-defined atom
groups and `tools/procar_character.py project`, `hall` and `plot`.
The public analytic fixture makes the complete command sequence runnable
without VASP; it is a software check, not a material benchmark.

For selected isolated bands, the character tool weights canonical charge
curvature with the measured state projections and scans occupations. Keep
the unweighted selected-band response beside those channels. A selection
containing only conduction bands describes that contribution, not the total
filled-valence plus conduction response. Arbitrary overlapping atom/orbital
groups do not form a partition. Read the normalization and spin-axis metadata
when comparing their weights.

This is a practical analysis of which states contribute to a charge response.
It does not construct a spin-, layer- or orbital-current operator. SOC spin
projections also do not create independent up/down eigenvalue channels.
Exactly degenerate states cannot generally be assigned a unique individual
projection-weighted curvature; use the documented isolation checks rather
than interpreting an arbitrary eigensolver basis as a physical decomposition.

### Occupations and degeneracies

For an unordered pair, the exporter stores
$`N^{ab}_{nm}=-2\mathrm{Im}(D^a_{nm}D^b_{mn})`$. Integration uses

```math
\frac{\sigma_{xy}}{e^2/h}=-\frac{S_{\rm BZ}}{2\pi}g_s
\sum_k w_k\sum_{n\lt m}(f_n-f_m)
\frac{N^{xy}_{nm}}{(E_n-E_m)^2}.
```

Equal-occupation pairs cancel before division. Thus exact internal degeneracies
of a filled band bundle need no individual-band curvature. A small denominator
with unequal occupations rejects the calculation by default. This is the
intrinsic clean-limit formula; it does not introduce lifetime broadening.

`--degeneracy-policy coalesce` explicitly approximates numerical energy groups
whose total spread is within the supplied tolerance by their mean energy.
Fermi occupations and denominators use those means; the original energies are
retained. Output records the largest energy and occupation changes. At T=0, a
chemical potential cutting the original unresolved spread is rejected. Check
that the result is insensitive to a tolerance small compared with physical
band splittings and thermal energies; this option does not resolve a real
band crossing or establish transport in a disordered metal.

For an existing canonical WAVECAR occupied-bundle CSV, `bundle-hall` provides
a simpler adapter. It rejects the new PAW trace schema: use `kubo-hall` with
the original optical run files for same-run validation. For the canonical CSV,
use `--csv`, `--wavecar`, `--occupied`, mesh/spin options and a μ range strictly
inside its common global insulating gap. It accepts only T=0 and a sampling
plane parallel to Cartesian xy, because the native bundle CSV stores only
the xy curvature. Use `pair-hall` for changes in occupation across band edges
or at finite temperature.

Both workflows preserve the native bare-momentum approximation. Increasing
k sampling or `NBANDS` cannot supply missing PAW/nonlocal/SOC velocity terms.
The total charge response of nonmagnetic MoS₂ should vanish by time-reversal
symmetry. Regional K/K′ contributions can be nonzero, but their cancellation
alone does not validate the physical magnitude of the approximate operator.
A quantitative QAH plateau additionally requires an appropriate full velocity
operator and convergence against an independent topological calculation.

## Spin-sector curvature

Use native [`--task spin-kubo`](SPIN_KUBO.md), available in 1.6.1, to resolve the
positive/negative projected-spin sectors of an isolated SOC band group.
The native executable writes point curvature and its spin-mixing correction
directly on arbitrary k paths. An explicit full `--mesh` also produces raw
approximation integrals. Use the separate [`spin-chern`](SPIN_CHERN.md)
task for geometric sector Chern numbers. The [graphene](../examples/materials/graphene-spin-chern/kubo/)
and [Bi](../examples/materials/bi-spin-hall/spin-chern-kubo/) examples show both
results and their numerical limits. This fixed-rank geometric calculation
is separate from occupation-weighted charge Hall and conventional spin-current
Hall transport.

## Native curvature of a band bundle

This section describes the explicit `--kubo-source wavecar` approximation.
For native `--task kubo`, `kubo-line` or `kubo-integral` with that option,
select a multi-band range to calculate its trace curvature from WAVECAR. For the eighteen occupied SOC bands in the
MoS₂ example:

```bash
build/vaspberry --task kubo --kubo-source wavecar --wavecar WAVECAR \
  --bands 1:18 --curvature-csv KUBO_BUNDLE.csv
```

Choose `--bands FIRST:LAST` for the physical subspace in your own material.
The same options work with the MPI binary. `--task kubo` and `kubo-line`
preserve the input point order without a Brillouin-zone integral.
`--task kubo-integral` also prints the finite-mesh trace integral; it requires
the full uniform mesh and matching `--mesh NX,NY`.

For a selected group $`I`$, let $`D_{\alpha,nm}`$ denote the matrix
element of $`\partial_{k_\alpha}H`$. Its eigenstate representation is

```math
\Omega_{I,xy}(k)=-2\mathrm{Im}
\sum_{n\in I}\sum_{m\notin I}
\frac{D_{x,nm}(k)D_{y,mn}(k)}{(E_n-E_m)^2}.
```

When individual bands are isolated, this equals their summed curvature.
Internal pairs contribute opposite terms and cancel. The bundle routine
removes those pairs before division, so internal degeneracies do not produce
undefined individual-band terms. For the occupied manifold, this is the
occupied-to-empty formulation of [Wang et al., Eq. (11), Sec. III D and
Appendix B](https://doi.org/10.1103/PhysRevB.74.195118).

Before calculating, the native routine checks every spin channel and sampled
k point. Any gap from the selected group to an excluded source band at or
below **1e-5 eV** rejects the run before an output CSV is created. This checks
the sampled points and available source states; k and `NBANDS` convergence
remain necessary. The native operator retains the existing bare-momentum
approximation.

The `--curvature-csv` path must be new; without it, a multi-band trace writes
`KUBO.csv` in the working directory. The trace writes a
[bundle CSV](OUTPUT_FORMAT.md#native-kubo-bundle-csv) with the selected range,
source band count and external gaps; it produces no legacy `.dat` files.
Trace CSV V2 starts as INCOMPLETE and records the expected source spin/k
coverage. Only the final PASS marks successful completion. The readers
reject truncated, duplicate or nonfinite results. Historical V1 inputs
remain readable as legacy-unverified data, with their original provenance.
The bundle has unit occupation throughout its selected subspace.
It cannot replace band energies and occupation weights in an arbitrary
chemical-potential or temperature scan.

Use `--bands N` for one isolated band. To request several separate band
curvatures, use, for example, `--bands 18:19 --per-band 1`. Each selected band
must be separated from every other stored band by more than 1e-5 eV at every
sampled k point and spin channel. These calculations write legacy DAT files
and, when requested with `--curvature-csv`, an individual-band CSV containing
band IDs and energies. They reject unresolved degeneracies before producing
results; use the trace of an externally isolated group when internal band
degeneracies are present.

The removed native `--bundle` and `-kubo_bundle` flags give an error with
migration guidance. A named Kubo task uses the new range semantics even with
legacy `-ii`/`-if` endpoints. A command using only legacy `-kubo` and band
flags retains individual-band output with the same isolation check. The
Python `bundle-hall` command remains a separate integration of an already
computed canonical native trace CSV; it does not select the native calculation
mode or accept the new PAW optical trace as a replacement for its source files.

The remainder of this guide describes the generic matrix and point-curvature
interfaces. They require the actual exported operator/curvature data specified
below. JSON metadata records what those arrays mean; it does not replace the
VASP wavefunctions or create a material calculation.

## Applying the workflow to other data

Write a [generic matrix bundle](OUTPUT_FORMAT.md#generic-interband-matrix-input)
with your actual Cartesian operator and provenance, then use `matrix` with
explicit n and m windows. Ranges such as `1:40` and lists such as `1,3,5:8`
use global one-based band IDs. `--energy-reference` describes the existing
energy zero; it does not shift energies or guess a Fermi level. Select
`--spin-multiplicity` from the represented physical states.

`--mesh NX NY` declares a full uniform 2D mesh, which is checked against the
coordinates and weights. Without it, output is marked `points` and cannot be
integrated by `hall`. `--plane-axes 0 1` selects an ordered reciprocal-axis
pair; other pairs are permitted when the data contain the required curvature
components. A line path and an irreducible k mesh are not full integration
meshes.

For regional contributions, save a JSON region specification as described in
[the format guide](OUTPUT_FORMAT.md#user-defined-reciprocal-space-regions).
Then add, for example (choose the chemical-potential values for your own
energy zero and available band window):

```bash
python3 tools/vaspberry_kubo.py hall \
  --curvature results/material-curvature \
  --mu-min -0.5 --mu-max 0.5 --mu-num 101 --mu-reference 0 \
  --temperatures 0 300 --regions regions.json \
  --difference A_minus_B:region_A:region_B \
  --output-dir results/material-regions
```

The difference is exactly left minus right; there is no implicit factor of
one-half. Region names and locations are supplied by the user. `rest` records
the part of the sampled plane outside the chosen regions, and `total` retains
the whole-plane result. A spin/layer label on a region or band does not turn
this charge response into a spin/layer-current calculation.

At T=0, a state at `E==mu` is fully occupied. At finite temperature the code
uses Fermi occupations, processing chemical potentials in chunks controlled
by `--mu-chunk` (default 64 for `hall`, 32 for `pair-hall` and `wavecar-hall`).
Changes from `--mu-reference` are integrated
from occupation differences. The zero-temperature path accumulates states
around the reference occupation boundary, and the finite-temperature path
uses `f(mu,T)-f(mu_reference,T)`. This avoids subtracting two large filled-band
responses to obtain a small doping change. If bands are omitted or the highest available
band is occupied over the requested scan, the default command refuses to
claim a complete occupation window. `--allow-partial-bands` explicitly reports
only the represented contribution; it does not supply a missing baseline.

## Observable, sign and units

The connection convention is

```math
A_{n,\alpha}=i\langle u_n|\partial_{k_\alpha}u_n\rangle.
```

For an isolated band and a Hermitian physical vertex
$`D_{\alpha,nm}=\langle u_n|\partial_{k_\alpha}H|u_m\rangle`$,

```math
\Omega_{n,\alpha\beta}(k)=
-2\mathrm{Im}\sum_{m\ne n}
\frac{D_{\alpha,nm}(k)D_{\beta,mn}(k)}{(E_n-E_m)^2}.
```

Cartesian k is measured in Å⁻¹, D in eV·Å, energies in eV and curvature in
Å². The three output components are $`(\Omega_{yz},\Omega_{zx},\Omega_{xy})`$,
equivalently the Cartesian curvature vector. Reciprocal vectors include
$`2\pi`$. If the supplied operator is velocity, convert with $`D=\hbar v`$
before writing the matrix contract. A bare momentum matrix requires the
appropriate physical conversion and is not generally the full velocity of a
nonlocal, PAW or SOC Hamiltonian.

For the plane spanned by the ordered reciprocal vectors $`b_1,b_2`$, define
$`S_{\rm BZ}=|b_1\times b_2|`$ and $`\hat n=(b_1\times b_2)/S_{\rm BZ}`$.
Normalized integration weights satisfy $`\sum_k w_k=1`$. The sheet response is

```math
C_{\rm occ}(\mu,T)\simeq\frac{S_{\rm BZ}}{2\pi}
g_s\sum_{kn}w_k f(E_{nk}-\mu,T)\,\Omega_n(k)\cdot\hat n,
\qquad \frac{\sigma_{xy}}{e^2/h}=-C_{\rm occ}.
```

Here `xy` denotes the declared oriented plane. It is the literal Cartesian
xy response only when the normal points along +z. No effective layer
thickness is assumed. The command does not calculate 3D bulk conductivity,
spin Hall conductivity, extrinsic scattering contributions, or finite-frequency
optical conductivity. The separate [spin Hall workflow](SPIN_HALL.md) uses
its own spin-current operator, input contract and tensor output.

The explicit multiplicity $`g_s`$ is 1 for SOC spinors or one spin channel,
and 2 only for a declared scalar spin-degenerate calculation. It counts
represented physical states; it is separate from the corrected circular-momentum
normalization and must not be added to an already enumerated spinor spectrum.

## Point curvature and plaquette flux

The Kubo formula yields Berry curvature at each k point. The Fukui method
yields the flux through a mesh plaquette in radians from a loop of link
variables. Dividing that flux by cell area is an area-averaged curvature
estimate, not a measurement of point curvature at the cell center.
Their finite-mesh occupation approximations differ. Keep the existing
`wavecar_fukui.py` workflows and their [valley-transport guide](VALLEY_TRANSPORT.md)
when working with geometric plaquette data.

Do not round a finite-mesh Kubo-formula curvature integral to make it integer.
A lattice Chern number computed with the Fukui method also needs a sufficiently
resolved isolated band or bundle before it can be identified with the continuum
invariant.

## Band windows and degeneracies

The selected output bands n and intermediate bands m are independent inputs.
For a single-band curvature, the complete expression includes **all other
bands**, below and above the band of interest. M truncation must be converged
separately from k sampling and from the accuracy of the computed eigenstates.
The number of states in an exported window is not automatically the source
calculation's full `NBANDS`.

The matrix calculator checks the requested gap threshold against **all band
energies present in the export**, independently of the chosen m sum, and
rejects small-gap pairs by default. Its explicit masking mode preserves undefined band
curvatures as invalid entries instead of adding a denominator broadening or
silently dropping pairs. `min_gap_eV` therefore does not increase merely
because a nearby band was removed from the m window. Diagnostics record
`isolation_check_scope` and `all_source_band_energies_available`; even a full
source-energy check does not prove isolation from states absent from the source
or between sampled k points. Internally
degenerate groups require an appropriate separated-subspace calculation; this
single-band matrix command is not a general non-Abelian bundle algorithm.
The native multi-band trace route above explicitly supports an internally
degenerate, externally separated selected group.

An empty isolated band's unit-occupation Chern number is a valid geometric question.
It is not the occupied Hall response of that material. In metals, electron and
hole occupations must enter the integration. Partial band-window contributions
must be identified as such; omitted occupied bands cannot silently be assigned
zero response. Region-restricted integrals are not generally integer invariants.

## Operator provenance and numerical checks

The generic input declares the operator's definition, included and missing
terms, accuracy label and supporting evidence. Hermiticity, coverage, units,
array dimensions and source hashes are checked without silently repairing
matrices. An `experimental` operator requires an explicit opt-in and remains
diagnostic in the output. A producer's `validated` label is preserved; the
postprocessor does not independently establish that physical claim.

For PAW/nonorthogonal representations, overlap and basis-connection terms
matter. Differentiating a generalized eigenproblem involves
$`\partial H-E_n\partial S`$; that expression alone is not automatically the
Hermitian physical vertex expected here. Nonlocal potentials, SOC and Hubbard
projectors can contribute to velocity. Validate the actual exporter against
the declared Hamiltonian and basis.

The optional `vasp_optical_export.py` adapter reads its explicitly versioned
experimental optical stream. It is not a claim that an unmodified VASP
installation writes that stream, and no licensed VASP source or exporter build
is supplied by this repository. Use the generic matrix contract for other
producers.

The native adapter requires a successful producer record that binds
`BERRY_CONNECTION.bin` to its SHA256 through `outputs` or `output_sha256`.
A checksum generated only at post-processing time cannot replace that producer
record. The adapter preserves the experimental operator status and raw data;
matching provenance does not certify physical completeness.

## Reproducible calculation record

Preserve the source matrices or legacy files, their metadata and checksums,
the exact software commit, CLI invocation, generated NPZ/JSON pair and Hall
outputs. Record the original energy reference, source band count, chosen n/m
windows, mesh/weights, region definitions, temperature and chemical-potential
range. Repeat with denser k sampling and larger accurately computed m windows.
File-format validation and physical convergence answer different questions.

<a id="standard-waveder-insulating-paw-hall-response"></a>

## Standard WAVEDER: selected and occupied PAW responses

`kubo-hall` defaults to WAVEDER; the dedicated `waveder-hall` command remains
available. Both read the standard `WAVEDER`, `WAVECAR`, `INCAR` and `OUTCAR`
from one completed VASP optical run. Choose mutually exclusive `--bands
SELECTOR` for a selected μ/T contribution or `--occupied N` for the original
insulating T=0 response. A selector accepts a single band, an inclusive range
or a comma-separated combination, such as `31,33:34`. All source bands remain
in the virtual sum. The selected output is labeled a contribution and does
not certify total AHC, including at finite temperature.

The initial supported producer is **VASP 5.4.4** with explicit
`LOPTICS=.TRUE.`, `LPEAD=.FALSE.` and `LNABLA=.FALSE.`, a static `NSW=0` run,
`ISYM=-1`, and `LREAL=.FALSE.`. Hybrid, meta-GGA, PEAD and other optical
branches are rejected. The producer must report exactly `DEG_THRESHOLD=0.002`
eV; altered producer thresholds are unsupported. Transitive producer-degenerate
groups may be retained when their response weights are equal. Unequal weights
require erased internal couplings and are rejected; this applies to both a
split geometric selection and μ/T weighting of a complete group. No rounding
of thermal tails authorizes a missing pair.

For the occupied compatibility mode, use the actual occupied count, full
mesh and energies from your run. Every μ/reference must stay in the global
gap and temperatures must be zero. For
example, after defining `N_OCC`, `NX`, `NY`, `MU_LO`, `MU_HI` and `MU_REF`
for a spinor calculation:

```bash
python3 tools/vaspberry_kubo.py waveder-hall \
  --input-dir path/to/optics-run --occupied "$N_OCC" \
  --mesh "$NX" "$NY" \
  --energy-reference 'unchanged VASP eigenvalue zero' \
  --mu-min "$MU_LO" --mu-max "$MU_HI" --mu-reference "$MU_REF" \
  --formats csv dat npz --output-dir results/paw-insulating-hall
```

The equivalent standalone command is `python3 tools/waveder_hall.py` with
the same arguments. These Hall commands default `--input-dir` to the
invocation working directory and read the four conventional filenames there.
The `--wavecar`, `--waveder`, `--incar` and `--outcar` options override only
the named file. Explicit relative CLI paths are relative to the invocation
working directory, not the selected input directory. Output paths do not move. Output uses the standard Hall tables and can be plotted
with `tools/plot_hall.py`. The insulating scan gives a constant response because
the T=0 occupations do not change within the gap; constancy alone does not
establish quantization. `delta_sigma_e2_over_h` is zero there. Collinear spin channels are
selected with `--spin` and summed separately, with multiplicity one.

For a mesh calculated as several VASP jobs, pass their genuine run directories
together, for example `--run-dir optics-part01 optics-part02 optics-part03`.
This compatibility alternative cannot be combined with `--input-dir`.
Individual file overrides are allowed for one run only; multiple chunks
require conventional filenames in each directory.
Every job must use `ICHARG=11`, `LCHARG=.FALSE.`, the same retained `CHGCAR`,
`POTCAR` and `POSCAR`, and identical INCAR parameters apart from `SYSTEM`.
The adapter checks each completed run and integrates only after their union
forms the declared full mesh without duplicates. It combines the calculated
curvatures in memory; it does not construct replacement WAVEDER or OUTCAR
files. For region definitions using k-point IDs, IDs follow run-directory
order followed by each run's original k-point order.

For a selected contribution, use for example:

```bash
python3 tools/vaspberry_kubo.py kubo-hall \
  --input-dir path/to/optics-run --bands 31:32 --mesh "$NX" "$NY" \
  --energy-reference 'unchanged VASP eigenvalue zero' \
  --mu-min "$MU_LO" --mu-max "$MU_HI" --mu-num 101 \
  --mu-reference "$MU_REF" --temperatures 0 300 \
  --output-dir results/paw-selected-contribution
```

These source bands and energies must match the intended selection. `--bands`
does not need `--allow-partial-bands`; that redundant flag, pair-band truncation
and coalescing controls are rejected on this route. Regions and differences
remain available. The [runnable selected example](../examples/features/waveder-selected/)
uses declared synthetic input to demonstrate the full command flow.

The stored band axes are `NBANDS` bras by `NDBANDS` derivative kets; they
need not be square. Let `g_n=s_n f_n`, with `s_n=1` for selected target bands
and zero otherwise. Every pair with nonzero `g_n-g_m` is required, using
either available orientation with its proper sign. The intermediate index
always spans all source NBANDS. Equal-weight internal terms cancel. For the
complete occupied T=0 trace this reduces to the empty-bra/occupied-ket block.
Missing required high-band pairs cause rejection; the adapter never pads
zeros, constructs a square replacement, or restricts virtual bands to `S`.

WAVEDER stores optical matrix elements in Å with complex64 precision. The
adapter accumulates the pair contractions in complex128, without an
extra energy denominator, factor one-half, or integer rounding. Valid zero
matrix elements remain zero. This is the longitudinal PAW optical operator
including projector and augmentation terms described by
[Gajdoš et al., Eqs. (29)–(30)](https://doi.org/10.1103/PhysRevB.73.045112).
Those terms enter VASP's matrix construction before WAVEDER is written;
the Hall adapter does not add the same PAW correction again after reading it.

The adapter checks dimensions, final OUTCAR eigenvalues and occupations,
lattice, k points, electron count and mesh against WAVECAR. Keep the four
files from the same run: WAVEDER itself has no energies or coordinates, so
their association cannot be proven from its header. Source records are
generated automatically. Finite-smearing source occupations are allowed;
selected-mode response occupations come from the requested μ/T on those
fixed states. The `--occupied` compatibility mode remains explicitly T=0.
Selected metadata records `scope=selected_band_contribution`,
`selected_band_ids`, `intermediate_band_ids`, `rectangular_coverage` and
`producer_cluster_checks`; `total_AHC_certified` and
`total_delta_AHC_certified` remain false. A subset's finite-temperature curve
must not be presented as the total response.

These checks establish the file, optical-branch and filling contract, not
complete physical-velocity validation for every Hamiltonian. In particular,
the current WAVEDER validator does not audit `LDAU`. A readable SOC/+U file
or a reproducible integral does not independently certify all response terms.
The separate full-velocity/spin producer's exclusion of Hubbard U must also
be respected; it is not an automatically validated substitute for such a run.

Converge both the k mesh and accurately computed empty states. The
[MnBi₂Te₄ example](../examples/materials/mnbi2te4-qah/) compares the Chern number
computed with the Fukui method against results from actual unmodified VASP
optical runs. It separates the coarse
WAVEDER integral from VASPBERRY's dense full-connection Wannier calculation; a
successful file/producer check does not establish integration convergence.

## Full-connection Wannier bands and Hall response

This is an optional supporting calculation on externally prepared operators,
separate from the main WAVECAR/VASPBERRY workflow. It is retained to reproduce
the technical report’s model checks; the native tutorials require none of it.

For dense integration of a VASP-derived Wannier model, use `wannier-import`,
`wannier-bands` and `wannier-hall`. VASPBERRY evaluates Hamiltonian and position
Fourier sums, occupied-bundle J0/J1/J2 and the weighted 2D Hall integral itself.
The separate `postw90` calculation is an independent verification. Both full
operator matrices are required; a conventional Hamiltonian-only `_hr.dat`
does not provide the missing basis connection.

See [the general guide](WANNIER_TRANSPORT.md) for input contracts, commands,
sampling and resource controls, and [the actual MnBi₂Te₄ tutorial](../examples/materials/mnbi2te4-qah/NATIVE_WANNIER.md)
for reproduction. This backend currently supports fixed insulating bundles
at T=0. It complements the existing WAVECAR μ/T workflow and retains explicit
model-convergence and operator provenance.

## References

- [Xiao, Chang and Niu, Rev. Mod. Phys. 82, 1959 (2010), Eq.1.13](https://doi.org/10.1103/RevModPhys.82.1959): spectral curvature formula.
- [Wang et al., Phys. Rev. B 74, 195118 (2006)](https://doi.org/10.1103/PhysRevB.74.195118): occupied-to-empty curvature, cancellation of occupied pairs and the subspace trace. See also the [2007 erratum](https://doi.org/10.1103/PhysRevB.76.169902).
- [Wannier90 Berry module documentation](https://wannier90.readthedocs.io/en/latest/user_guide/postw90/berry/): connection, velocity and intrinsic Hall conventions.
- [Fukui, Hatsugai and Suzuki, JPSJ 74, 1674 (2005)](https://doi.org/10.1143/JPSJ.74.1674): geometric lattice Chern-number calculation and continuum limit.
- [Gajdoš et al., Phys. Rev. B 73, 045112 (2006)](https://doi.org/10.1103/PhysRevB.73.045112): PAW optical matrix elements and generalized overlap terms.

## Developer analytic check (separate from VASP tutorials)

The [developer fixture catalog](../validation/models/) retains these formula
checks. Use the real-input tutorials above to learn VASPBERRY from VASP files.

From the repository root, use fresh output directories:

```bash
python3 -m pip install -r requirements-transport.txt
python3 tools/vaspberry_kubo.py demo \
  --mesh 32 --mass -1 --output-dir results/qwz-matrix

python3 tools/vaspberry_kubo.py matrix \
  --matrices results/qwz-matrix/matrix.npz \
  --metadata results/qwz-matrix/matrix.json \
  --n-bands 1:2 --m-bands 1:2 \
  --mesh 32 32 --plane-axes 0 1 \
  --spin-multiplicity 1 --energy-reference "QWZ model zero (eV)" \
  --degeneracy-threshold-eV 1e-8 \
  --output-dir results/qwz-curvature

python3 tools/vaspberry_kubo.py hall \
  --curvature results/qwz-curvature \
  --mu-min -0.5 --mu-max 0.5 --mu-num 101 --mu-reference 0 \
  --temperatures 0 300 --band-resolved \
  --output-dir results/qwz-hall
```

`demo` creates the two-band Qi–Wu–Zhang model's `matrix.npz` and `matrix.json`.
It is a complete model fixture, not a simulated material. With this mass and
orientation the lower-band continuum Chern number is +1, so the zero-temperature
occupied sheet response at μ=0 approaches −1 in units of e²/h. Increase the
mesh to inspect convergence; the tool does not round the result. The upper
band is retained to make both the spectral sum and occupation window explicit.

The `matrix` command writes `curvature.npz`, `curvature.json` and
`curvature.csv`. `hall` writes `conductivity.npz`, `conductivity.json` and
`conductivity.csv`. In Hall output `band_id=0` means the sum of all represented
bands; `--band-resolved` adds rows for each band. The μ reference is included
in the output even when it is not one of the requested evenly spaced points.

Plot that standardized CSV with the standalone public example:

```bash
python3 examples/kubo/plot_hall.py \
  --input results/qwz-hall/conductivity.csv \
  --output results/qwz-hall/hall.png
```

See [plotting options](../examples/kubo/) for different regions, bands,
temperatures, quantities and figure formats. The plotter displays stored
values and does not make a convergence claim.
