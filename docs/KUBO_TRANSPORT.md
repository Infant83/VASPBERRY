# Kubo curvature and two-dimensional charge Hall transport

Actual VASP-based tutorials and reference results are provided for
[native MoS₂ Kubo curvature](../examples/features/kubo-curvature/) and
[Bi occupied-subspace Fukui Hall](../examples/features/hall-valley/). Their
methods and required sampling are explicitly different.

The `tools/vaspberry_kubo.py` command converts declared interband matrices or
legacy Kubo data into a common point-curvature format and integrates a
two-dimensional intrinsic charge Hall response. It also provides a public
analytic model demonstration. Python 3.10+ and the dependencies in
`requirements-transport.txt` are required.

Use `python tools/vaspberry_kubo.py --help` and each subcommand's `--help` for
the available options. The [output specification](OUTPUT_FORMAT.md) describes
the interchange files; [migration](MIGRATION.md) covers older normalization.

## Start from actual VASP output

Use the [MoS₂ tutorial](../examples/features/kubo-curvature/) to run the native
Fortran Kubo calculation on the supplied real WAVECAR. It includes the actual
command, raw high-precision CSV, plotted reference and unresolved-band mask.
The band-path input is not a full integration mesh. The native operator is
canonical momentum; PAW/nonlocal/SOC velocity corrections are not supplied by
that approximation.

The [Bi charge-Hall tutorial](../examples/features/hall-valley/) uses a different
production route: full-mesh WAVECAR occupied-subspace Fukui transport in the
insulating gap. Bi's unresolved Kramers pairs prevent treating its individual
bands as isolated point-Kubo input. Its zero charge Hall plateau is an actual
material sanity check, not a nonzero valley-Hall demonstration.

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
by `--mu-chunk` (default 64). Changes from `--mu-reference` are integrated
from occupation differences. The zero-temperature path accumulates states
around the reference occupation boundary, and the finite-temperature path
uses `f(mu,T)-f(mu_reference,T)`. This avoids subtracting two large filled-band
responses to obtain a small doping change. If bands are omitted or the highest available
band is occupied over the requested scan, the default command refuses to
claim a complete occupation window. `--allow-partial-bands` explicitly reports
only the represented contribution; it does not supply a missing baseline.

## Observable, sign and units

The connection convention is

\[
A_{n,\alpha}=i\langle u_n|\partial_{k_\alpha}u_n\rangle.
\]

For an isolated band and a Hermitian physical vertex
\(D_{\alpha,nm}=\langle u_n|\partial_{k_\alpha}H|u_m\rangle\),

\[
\Omega_{n,\alpha\beta}(k)=
-2\operatorname{Im}\sum_{m\ne n}
\frac{D_{\alpha,nm}(k)D_{\beta,mn}(k)}{(E_n-E_m)^2}.
\]

Cartesian k is measured in Å⁻¹, D in eV·Å, energies in eV and curvature in
Å². The three output components are \((\Omega_{yz},\Omega_{zx},\Omega_{xy})\),
equivalently the Cartesian curvature vector. Reciprocal vectors include
\(2\pi\). If the supplied operator is velocity, convert with \(D=\hbar v\)
before writing the matrix contract. A bare momentum matrix requires the
appropriate physical conversion and is not generally the full velocity of a
nonlocal, PAW or SOC Hamiltonian.

For the plane spanned by the ordered reciprocal vectors \(b_1,b_2\), define
\(S_{\rm BZ}=|b_1\times b_2|\) and \(\hat n=(b_1\times b_2)/S_{\rm BZ}\).
Normalized integration weights satisfy \(\sum_k w_k=1\). The sheet response is

\[
C_{\rm occ}(\mu,T)\simeq\frac{S_{\rm BZ}}{2\pi}
g_s\sum_{kn}w_k f(E_{nk}-\mu,T)\,\Omega_n(k)\cdot\hat n,
\qquad \frac{\sigma_{xy}}{e^2/h}=-C_{\rm occ}.
\]

Here `xy` denotes the declared oriented plane. It is the literal Cartesian
xy response only when the normal points along +z. No effective layer
thickness is assumed. The command does not calculate 3D bulk conductivity,
spin Hall conductivity, extrinsic scattering contributions, or finite-frequency
optical conductivity.

The explicit multiplicity \(g_s\) is 1 for SOC spinors or one spin channel,
and 2 only for a declared scalar spin-degenerate calculation. It counts
represented physical states; it is separate from the corrected circular-momentum
normalization and must not be added to an already enumerated spinor spectrum.

## Point curvature and plaquette flux

Kubo yields a value at each k point. A Fukui loop yields the flux through a
mesh plaquette in radians. Dividing that flux by cell area is an area-averaged
curvature estimate, not a measurement of point curvature at the cell center.
Their finite-mesh occupation approximations differ. Keep the existing
`wavecar_fukui.py` workflows and their [valley-transport guide](VALLEY_TRANSPORT.md)
when working with geometric plaquette data.

Do not round a finite-mesh Kubo integral to make it integer. A lattice Fukui
integer also needs a sufficiently resolved isolated band or bundle before it
can be identified with the continuum invariant.

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

An empty isolated band's unit-occupation Chern is a valid geometric question.
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
\(\partial H-E_n\partial S\); that expression alone is not automatically the
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

## References

- [Xiao, Chang and Niu, Rev. Mod. Phys. 82, 1959 (2010), Eq.1.13](https://doi.org/10.1103/RevModPhys.82.1959): spectral curvature formula.
- [Wannier90 Berry module documentation](https://wannier90.readthedocs.io/en/latest/user_guide/postw90/berry/): connection, velocity and intrinsic Hall conventions.
- [Fukui, Hatsugai and Suzuki, JPSJ 74, 1674 (2005)](https://doi.org/10.1143/JPSJ.74.1674): geometric lattice Chern calculation and continuum limit.
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
orientation the lower-band continuum Chern is +1, so the zero-temperature
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
