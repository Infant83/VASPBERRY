# Changelog

All notable changes to VASPBERRY are recorded here.

## [Unreleased]

### Added

- Readable native Fortran task and input aliases, including `--task chern`,
  `z2`, `kubo` and `kubo-pairs`, with argument validation. Existing short
  options and numerical methods remain available.
- Full-program serial/MPI Fukui regression for curvature, Chern number and
  extrema metadata, exercised by CI alongside the native command tests.
- Canonical serial executable `build/vaspberry`, with the previous
  `build/vaspberry-gfortran` name retained for existing scripts.

### Fixed

- Populate Fukui MPI output-header curvature extrema and their k points from
  the gathered mesh, matching the printed numerical data.
- Avoid a zero divisor in the wavefunction progress display for real-space
  grids with fewer than five z points. The reconstruction formula is unchanged.

### Documentation

- Lead README, examples and the technical report with direct VASP WAVECAR →
  Fortran calculation → saved results → postprocessing/plotting workflows.
  Separate optional PAW/spin operator extensions and Wannier supporting checks.
- Show native commands before optional example runners, and distinguish
  Python occupation-weighted Hall integration from plotting.

## [1.3.0] - 2026-09-25

### Added

- Actual ordinary-VASP MoS₂ monolayer and 1H/2H/3R bilayer examples, with
  documented fixed structures, band energies, complete-group circular
  strengths, optional standard PAW optics and reproducible comparison figures.
- Matched MoS₂ charge-Hall calculations from the same VASP eigenstates using
  canonical momentum and optional full PAW velocity. Public pair caches
  reproduce the chemical-potential, temperature and band-cutoff comparisons.

- Standard-VASP `waveder-optics` for PAW circular transition strengths and
  k-resolved spectra, with complete initial/final degenerate groups, explicit
  polarization conventions and selectable CSV/DAT/NPZ output.
- `velocity-pairs` converts an audited full PAW velocity bundle to the same
  occupation-pair format used by the ordinary WAVECAR charge-Hall workflow.
  Operator provenance and all original eigenstates remain explicit.

- Conventional insulating 2D `spin-hall`, full complex `spin-matrix`, and
  audited PAW `spin-export` commands. Spin and full physical velocity share
  one VASP eigenstate basis, including diagonal/degenerate velocity blocks.
  Spin-current finite-band closure, sheet units and PAW augmentation are
  explicit; k sampling and source bands require separate convergence.
- `wannier-edge` for ideal finite-strip bands and boundary probabilities,
  with memory/time estimates, no-overwrite output and failure records.

- VASPBERRY-owned `wannier-import`, `wannier-bands` and `wannier-hall` commands
  for complete real-space Hamiltonian/position inputs. The NumPy kernel retains
  full J0/J1/J2, uses only cross-gap denominators and handles internal occupied
  degeneracies. Fixed insulating T=0 sheet response supports uniform or locally
  refined full-zone quadrature, shared-memory batches, resource estimates,
  source-integrity checks and the common CSV/DAT/NPZ Hall output.
- Reorganize the MnBi₂Te₄ example around VASP-derived Wannier bands and Hall
  response evaluated from those operators, with postw90 retained as an
  independent comparison. Distinguish the 123-band VASP Fukui bundle from the
  87 occupied model states and retain the unconverged direct optical diagnostic.

- Native `-kubo_pairs` export of unordered interband numerators, with serial/MPI
  parity and bounded coefficient caching. Occupation-weighted pair integration
  cancels equal-occupation internal transitions before denominator evaluation.
- `wavecar-hall`, `import-pairs`, `pair-hall` and insulating `bundle-hall`
  commands, reusable NPZ data, independently selectable CSV/DAT/NPZ tables and
  a common PNG/PDF/SVG plotting tool. Numerical degeneracy coalescing is explicit
  and records the energy/occupation approximation; the default rejects unresolved
  unequal occupations.
- An explicit `--pair-band-max` cutoff separates the Kubo intermediate-state
  window from the stored VASP `NBANDS`, retaining source dimensions and checking
  occupations and unresolved degeneracies at the cutoff.
- `waveder-hall` for zero-temperature insulating occupied-bundle response from
  standard VASP 5.4.4 longitudinal PAW optical files, with source consistency,
  operator, electron-count and global-gap validation.
- Actual MoS₂ Hall scans and separate mesh/pair-window studies, with original
  VASP conditions, numerical references in CSV/DAT/NPZ and PNG/PDF/SVG figures.
  A 21-atom MnBi₂Te₄ film adds magnetic Chern, optical integration and independent
  Wannier-reference examples, with preparation files and a public SCF density.
- Native selected-bundle Kubo trace via `-kubo_bundle 1`, with analytic
  exclusion of internal transitions, external-gap checks and a versioned
  `-kubo_csv` output. Serial and MPI support the same interface; the default
  individual-band calculation remains unchanged.
- A scientific technical report and an explicit MoS2 full-zone Fukui curvature
  tutorial, recalculated from a new 12x12 SOC VASP mesh using the public charge
  density. The VASP preparation, native output, numerical table and figures
  are supplied; the large mesh WAVECAR is generated by the preparation recipe.
- Cartesian first-Brillouin-zone curvature maps and reduced-coordinate
  n-field plots, plus publication-style band-path, optical, Hall and
  real-space figures in PNG and PDF. Existing numerical results are preserved.
- Real VASP-based feature tutorials: Bi occupied-subspace Fukui/Chern, Z2 and
  gap Hall; MoS2 occupied-bundle and valley-resolved Kubo curvature, optical
  response and Gamma wavefunctions. Each
  gives actual input files, explicit production commands, freshly calculated
  original outputs, numerical/figure references and steps for another system.
  Public Bi input can be fetched with size/checksum verification. Model and
  synthetic checks are retained separately under validation/models/.
- General matrix-to-curvature, explicit legacy-normalization import, and
  two-dimensional intrinsic charge Hall commands in `tools/vaspberry_kubo.py`.
- A versioned NPZ plus JSON format for pointwise Berry curvature, energies,
  reciprocal coordinates, integration weights, selected bands and provenance.
- User-defined reciprocal-space regions for decomposing the sheet Hall
  response, with a public analytic two-band demonstration that needs no VASP
  files or proprietary data.
- A standalone Matplotlib example that plots selected regions, bands,
  temperatures and observables from the standardized Hall CSV.
- Method, output-format and migration guides documenting sign, units,
  band-sum truncation, physical-operator scope and the distinction between
  point curvature and Fukui plaquette flux.

### Fixed

- Stabilize the optional full-velocity producer's optical consistency check
  near numerical degeneracies. Apply the recorded 2 meV optical exclusion
  only in that branch, preserve undivided velocity/spin matrices exactly,
  and retain the strict comparison tolerance and legacy optical-only policy.
- Clarify VASP electronic-structure provenance throughout the technical report
  and band figure labels; present the Bi ideal edge as corroboration of its
  bulk Z₂ index. Separate ordinary and optional operator preparation guides.

- Match the wavefunction phase-array assignment to its active plane-wave count
  and initialize the POSCAR-header parsing loop in both Fortran sources. The
  synthetic Gamma example now runs with array bounds checks and verifies the
  original real/imaginary output against a known state.
- Corrected the legacy Fortran Kubo circular-momentum normalization: the
  curvature now uses the standard `-2 Im` factor for
  `A = i<u|d_k u>`. Earlier unnormalized circular components produced twice
  this curvature. The migration guide explains how to identify and import
  old outputs without silently rescaling already normalized results.
- Preserve custom `-o` output names by avoiding overlapping Fortran string
  reads and writes. Corrected the serial job-displacement array bounds.
- Check individual-band isolation against all exported energies even when
  the intermediate sum is truncated, and recompute imported gaps from the
  matching WAVECAR. Invalid individual-band data cannot enter Hall integration.
- Check full-mesh geometry, spin multiplicity, reciprocal-region overlap and
  producer checksums at the portable data boundary.

### Performance and numerical stability

- Validate native pair CSV in bounded row chunks and stream Hall integration
  over k points and chemical-potential blocks. Reuse the imported pair cache
  for additional scans without repeating native matrix-element calculations.
- Evaluate only selected-to-excluded pairs for native bundle curvature and
  construct plane-wave momentum coordinates once per k point. Coefficient
  storage remains proportional to the plane-wave count.
- Use sorted cumulative sums at zero temperature and share finite-temperature
  occupation arrays across regions and bands in bounded chemical-potential
  chunks. Band-resolved work avoids a full state array for each band mask.
- Integrate doping differences from the reference occupation directly, so
  subtraction of a large filled-band baseline does not erase a small signal.
- The default GNU serial and MPI builds now use the same current source.
  The historical reduced serial source remains explicitly available.

### Compatibility and scientific scope

- Existing Fukui command lines and plaquette outputs remain available.
  They are not relabeled as pointwise Kubo curvature.
- The normalization fix changes the legacy Kubo magnitude; it does not add
  missing PAW, nonlocal, SOC or Hubbard-potential velocity terms. A supplied
  matrix must represent the declared Hamiltonian and physical operator.
- Validation of array shapes, units and numerical properties is distinct
  from material-specific k-mesh, intermediate-band and basis convergence.
  Integer Chern numbers are not enforced by rounding finite-mesh integrals.
- The charge Hall commands compute a two-dimensional intrinsic **charge**
  sheet response. The separate `spin-hall` command adds conventional spin
  sheet response for insulating T=0 bundles; no 3D bulk conversion is inferred.
- Source release `v1.3.0` fixes this version's code, documentation and examples
  under an immutable Git tag and GitHub Release. The default branch is `master`.
  Build executables locally; no new DOI or prebuilt binary is assigned.

## [1.2.0] - 2026-09-04

### Added

- Promoted `-z2 1` to the documented two-dimensional Fukui-Hatsugai lattice
  n-field Z2 interface. `-h` now states the full even Gamma-centered mesh,
  SOC spinor, occupied-band, output, parity, and mesh-convergence contract.
- Added reviewed two-stage Bi input templates: an SOC SCF calculation that
  writes `CHGCAR`, followed by an `ICHARG=11` fixed-charge calculation that
  writes the full-mesh spinor `WAVECAR`.
- Added a plotting helper for reportable schema-2 `Z2_FIELD.csv` outputs and a
  rendered 12 x 12 Bi n-field. The checked reference has top/bottom sums
  -3/+3 and matching parity 1.
- Added the schema-2 Bi 12 x 12 reference CSV used by the guide, plotter, and
  example-layout regression tests.
- Added GNU serial/OpenMPI Makefile targets, a two-rank MPI datatype/status
  smoke test, compiler-specific build documentation, and reviewed manual
  recipes for Intel `ifx` and legacy `ifort`.
- Added a roadmap that explicitly defers an importable Python library port
  until it has golden numerical comparisons with the Fortran implementations.

### Changed

- PASS field results now record
  `result_kind=FUKUI_HATSUGAI_NFIELD_Z2`,
  `reportable_invariant=1`, a numeric `z2_invariant`, and matching top/bottom
  parity fields. Rejected results remain non-reportable. The Z2 field schema
  is version 2; old ownership markers remain recognized for safe cleanup.
- Replaced the nonstandard `MPI_REAL8` datatype in the Z2 reductions with
  `MPI_DOUBLE_PRECISION`, and made MPI help exit through `MPI_FINALIZE`.
- Made the legacy extended-BZ parser pass an explicit length-one character to
  `ICHAR`, avoiding an Intel ifx 2025.0 front-end failure while retaining the
  parser's behavior.
- Moved the Bi and 1H-MoS2 data under `examples/`, separated the Bi 2016 raw
  archive from reviewed new-run templates, and separated full-BZ MoS2 data
  from its K-Gamma-K' line workflow.
- Removed the previous Python Wilson-loop CLI, its detailed guide, and its
  dedicated tests from the active tree. An independent method remains an
  optional literature cross-check and is not part of `-z2`.
- Removed historical-field comparison from the regular n-field plotting CLI
  and active Z2 guide. Earlier result details remain confined to this
  changelog and the explicitly named legacy archive.

### Data and portability scope

- Removed all four tracked MoS2 `POTCAR` copies from the current tree and
  added non-proprietary provenance plus a CI guard against future tracking.
  Historical Git objects require a separately authorized history rewrite to
  purge completely.
- GNU serial and OpenMPI builds are configured for CI on Ubuntu 22.04 and
  24.04. Intel `ifx` 2025.0 and retained `ifort` 2021.10 serial builds are
  compile/help-tested in CI; Intel MPI numerical validation remains manual.
- The bundled Bi `WAVECAR` reproduces VASPBERRY post-processing. The complete
  preceding SCF provenance for its archived `CHGCAR` is unavailable, so it is
  not described as an end-to-end DFT reproduction package.

## [1.1.1] - 2026-08-31

### Fixed

- Corrected spinor time-reversal reconstruction at nonzero TRIMs with
  `G_target = round(-k_source-k_target)-G_source`. At a represented TRIM,
  this is `G_target=-G_source-2*TRIM`; the old shift-free rule was valid
  only at Gamma.
- Replaced coefficient-position assumptions with an explicit reciprocal-G
  bijection and hard checks for missing, duplicate, noninteger, or
  norm-changing mappings.
- Applied the spin-1/2 operation
  `Theta(C_up,C_down)=(-conjg(C_down),conjg(C_up))` through a shared helper
  and retained the preceding odd state as the source of an even TRIM Kramers
  partner.
- Promoted the real and imaginary parts of single-precision WAVECAR
  coefficients before norm squaring, avoiding false failures from
  single-precision `abs` rounding.
- Replaced the rank-dependent `abs(det S)>1e-14` gate with a `ZGESVD`
  minimum-singular-value measurement. Determinant phases now come from a
  separate `ZGETRF` copy, including pivot parity, without determinant
  multiplication.
- In the archived Bi 12 x 12 regression, these corrections changed five of
  144 fundamental plaquettes and replaced the inconsistent -1/0 half-zone
  sums with -3/+3.
- Corrected the sign and units in the legacy `NFIELD.dat` explanation.
- Preserved a custom legacy Z2 `-o` base without the general
  `BERRYCURV.` prefix, enlarged the serial CLI value buffer, and added
  executable-level 71-character-pass/72-character-reject coverage.

### Added

- `Z2_FIELD.csv` on the fundamental mesh for PASS results, and
  `Z2_FIELD.invalid.csv` for INVALID or INCOMPLETE runs. Same-directory
  temporary output and C/POSIX `rename()` prevent partial publication.
  Reserved basenames and POSIX `realpath()` checks reject direct, relative,
  absolute, and symbolic-link input aliases. Schema/marker preflight refuses
  to delete unowned files. After successful preflight, stale PASS plus legacy
  NFIELD products are removed before WAVECAR processing. Legacy output is
  closed and atomically published first; the INCOMPLETE sentinel is removed
  next; `Z2_FIELD.csv` is the final PASS commit marker. If that final rename
  fails, no regular PASS or sentinel remains and the staged temporary CSV is
  retained with a nonzero exit. Z2 `-o` bases longer than 71 characters are
  rejected before cleanup to prevent path truncation.
- Version, grid, band range, spinor rank, overlap backend, numerical
  thresholds, minimum link singular values, per-plaquette status, and an
  explicit check-scope statement in the field CSV.
- Production-linked Fortran regression tests for Γ/M1/M2/M3 reciprocal
  mapping, `Theta^2=-1`, an integrated M1 `get_z2_state`
  coefficient reconstruction, complex*8 norm accumulation, LU pivot phase,
  minimum singular values, canonical-path aliases, owned-output guards,
  legacy/final atomic commit failures, and output-name length boundaries.
  Both production objects are linked into and exercised by helper drivers;
  MPI communication ordering is source-checked.
- Public 1H-MoS₂ regression coverage that verifies all nine periodic copies
  of each of the 144 fundamental plaquettes before checking the wrapped
  time-reversal-odd Berry curvature and zero occupied-subspace Chern number.

### Scientific scope

- The fixed reciprocal-G shift is an independently verified code-level
  defect. Quantitative attribution of a previously observed MoS₂ Z2-field
  asymmetry requires a new corrected run on the corresponding full-mesh
  WAVECAR; that private input is not present in CI.
- The Fortran B-minus and even-TRIM states are constructed with time reversal.
  Its flux oddness and zero-Chern checks therefore test reconstruction
  self-consistency, not the raw WAVECAR's physical time-reversal symmetry.
  The CSV records `input_trs_independently_verified=0`.
- The direct overlap backend remains
  `WAVECAR_PSEUDO_NO_PAW_AUGMENTATION`. PAW augmentation can change the
  complex finite-neighbor overlap matrices, including phases and conditioning,
  but it is not the source of the omitted reciprocal-G shift. With a correct
  TR construction and paired mesh, its consistent omission is not by itself
  an exact TR-covariance breaker, although it can amplify coarse-mesh or
  branch-cut failures.
- The pointwise Fukui n-field remains gauge- and logarithm-branch-dependent.
  The physical relation is the wrapped flux condition
  `wrap[flux(-k)+flux(k)]=0` modulo `2*pi`.
- Thresholds are conservative diagnostic policies and can yield false
  rejection on a coarse mesh. At the time of version 1.1.1, the guarded
  Wilson-loop workflow and mesh convergence were required before reporting
  an invariant; version 1.2.0 supersedes that reporting policy for the
  validated Fukui-Hatsugai n-field result.

### Distribution status

- Version 1.1.1 identifies the reviewed source-tree candidate. It does not
  create a Git tag, GitHub Release archive, or Zenodo record.
- The archived DOI `10.5281/zenodo.1402593` remains specific to VASPBERRY
  V1.0 (2018).

## [1.1.0] - 2026-08-31

### Added

- Read-only Python tools for high-precision Fukui maps and guarded,
  chemical-potential-dependent Hall transport directly from a full-mesh
  `WAVECAR`.
- Cumulative occupied-subspace transport through a chosen valence-band
  maximum, including separate K and K' contributions, strict CSV/JSON
  diagnostics, and an optional Cartesian first-Brillouin-zone plot.
- A class-AII Wilson-loop/Wannier-charge-centre workflow that reports a
  Z2 invariant only after its time-reversal, gap, link, Chern, Kramers,
  and fixed-grid convergence checks pass.
- Regression tests for analytic topological models, source-level failure
  paths, output safety, release metadata, and the transport formulas.
- GitHub Actions checks for the Python suite and serial/MPI GNU Fortran
  compilation.

### Changed

- Corrected the legacy Fortran Z2 control-flow, accumulator initialisation,
  periodic indexing, neighbour validation, and fatal-error handling.
- Labelled the Fortran half-zone result as a **legacy Z2 candidate**. It is
  not a reportable invariant until the guarded Wilson-loop workflow passes.
- Corrected the spin-resolved Kubo accumulator lifetime while preserving the
  existing default Fukui calculation and output format.
- Added collision checks and stale-output cleanup to prevent a rejected
  calculation from leaving apparently valid transport or Z2 products.
- Updated all program and tool version strings to 1.1.0.

### Scientific scope

- A candidate `Z2=0` or `Z2=1` from an unconverged mesh remains diagnostic;
  it is emitted as `z2: null` in strict JSON output.
- Numerical time-reversal residuals do not by themselves establish the
  material's magnetic state or physical symmetry.
- The first-BZ option changes the map representation only. Fukui integration
  continues on the original periodic reciprocal-coordinate torus.
- The archived Zenodo DOI `10.5281/zenodo.1402593` identifies VASPBERRY V1.0
  (2018). It is not assigned to version 1.1.0.

### Transport interfaces

- `tools/wavecar_fukui.py --transport-t0` writes guarded three-manifold
  `transport_t0.csv` data and `transport_t0_diagnostics.json`.
- `--transport-full-t0 MAX_BAND` writes the cumulative valence-window
  `transport_full_t0.csv` and its diagnostics JSON. K and K' contributions are
  reported separately from their sum and contrast.
- `--plot` derives PNG views from validated numerical products.
  `--plot-domain first-bz` changes only the k-resolved display coordinates;
  Fukui integration remains on the original periodic reciprocal-coordinate
  torus.
- `--allow-invalid-transport` produces explicitly unvalidated diagnostics; it
  does not turn a rejected point into a validated conductivity.

### Z2 interface and result status

Run the guarded class-AII calculation on a uniform full reciprocal mesh with,
for example:

```text
python tools/wavecar_z2.py WAVECAR --nx 12 --ny 12 \
  --occupied-bands 18 --axis both --output-dir z2_result --plot
```

- `z2_wilson_wcc.csv` contains the WCC data,
  `z2_diagnostics.json` records every guard and the validated result, and
  `z2_wilson_wcc.png` is the optional plot.
- Exit 0 means the enabled guards pass and `z2` is reportable. Exit 2 records
  `z2: null` and a diagnostic `candidate_z2`; exit 1 is a hard error and
  removes planned products.
- Time reversal constrains the WCC spectrum as a set at pump coordinates `q`
  and `-q`, with Kramers pairing at the invariant endpoints. Sorted branches
  may exchange partners.
- POS/GAP/MOVE-style checks are fixed-grid convergence proxies, not an
  adaptive Z2Pack calculation.

Method references: Fukui and Hatsugai,
<https://doi.org/10.1143/JPSJ.76.053702>; Yu *et al.*,
<https://doi.org/10.1103/PhysRevB.84.075119>; and Gresch *et al.*,
<https://doi.org/10.1103/PhysRevB.95.075146>.

### Distribution status

- Version 1.1.0 identifies the updated source tree. This update does not create
  a Git tag, GitHub Release archive, or Zenodo record.
- A formal packaged release remains subject to maintainer review of the
  repository licence and tracked VASP input artefacts, including any licensed
  `POTCAR` data.

## [1.0] - 2018-08-23

- Initial archived release: Berry curvature and Chern number by the Fukui
  method, circular dichroism, and real-space wavefunction output.
- Archive: <https://doi.org/10.5281/zenodo.1402593>.
