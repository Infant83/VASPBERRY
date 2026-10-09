# Changelog

All notable changes to VASPBERRY are recorded here.

## [Unreleased]

## [1.6.6] - 2026-10-04

This patch retains the 1.6.5 WAVEDER input protocol and packages the
CLI-consistency fixes and report reproduction corrections below. Historical
material reference data retain their original producer versions.

### Help and workflow clarity — updated 2026-10-09

- Organize every native task and option page into readable plain-text sections,
  with copyable commands, required inputs, output quantities and units,
  defaults, validity conditions and related guides. Keep lines within 80 columns.
- Correct legacy option names, Chern output labels, optical spectrum quantities,
  periodic-repeat dimensions and historical defaults. Qualify canonical Kubo
  compatibility rules and CSV schemas by source; retain existing help aliases.
- State required wavefunction POSCAR/EIGENVAL inputs, spin task prerequisites,
  singleton band ranges and task-specific file replacement rules.
- Group Python CLI options, expose useful defaults/units and saved artifact
  types, and carry scientific scope into direct subcommand help. Clarify
  WAVEDER direct calculation versus native WAVECAR execution, conditional
  selections, plot input requirements and retained legacy replacement behavior.
- Keep numerical kernels, selection grammar, defaults and output schemas intact;
  extend help checks to cover presentation and parser contracts.

### Report reproduction

- Recalculate the technical-report figures and numerical controls from
  retained inputs, and record the dated reproduction scope separately from
  electronic-structure regeneration and physical convergence.
- Correct source-specific Bi conditioning diagnostics in four reference
  rows, retaining the original records and unchanged Chern/Kubo values.
- Clarify rounded report values and the Hall example's pair-window help;
  wrap Kubo plot operator labels to prevent clipping.

### Command-line clarity and validation

- Reject wavefunction, optical and folding controls that modern standard
  WAVEDER Kubo tasks do not use. Warn when a later band selector overrides
  an earlier one, preserving the supported historical precedence.
- Keep native `--mesh` integration distinct from legacy mesh-size settings.
  Reject nondefault `--mu-chunk` on the WAVEDER route, where it has no effect.
- Match CLI and INI chemical-potential scans: equal endpoints require one
  point; a multi-point scan requires increasing endpoints.
- Explain that INI `check` evaluates the requested WAVEDER response without
  writing results, show every resolved optical input, and distinguish this
  validation from material convergence.
- Correct standard/canonical output filenames, rescan and cache-reuse
  guidance, protocol navigation and the two-endpoint pair-band cutoff label.

### CI launcher validation

- Stabilize Intel MPI CI launcher signal handling while retaining controlled
  rejection status, expected diagnostics and output guards.
- Make a best-effort flush of diagnostic streams before fatal termination;
  retain the existing MPI_ABORT call and error code without adding collectives
  or changing numerical kernels. Linux MPICH archive validation explicitly
  retains per-process diagnostic logs for its rejection checks.

## [1.6.5] - 2026-10-04

This release incorporates and supersedes the locally prepared 1.6.4 candidate,
which was never published. Saved scientific reference data keep their
original producer versions.

### Release documentation

- Align the README, native/Python commands, selected-band example and
  technical report with the 1.6.5 release. Correct the current-source
  installation links and migration anchors.
- Distinguish direct WAVEDER Hall integration from the separate native
  WAVECAR pair-export workflow throughout the introductory guides.
- Retain the original 1.6.4 validation record as candidate history and
  record new release checks separately in `docs/VALIDATION_1.6.5.md`.

### Explicit input directory

- Add `--input-dir DIR` to native and Python Hall commands, with the
  invocation working directory as the default. Per-file overrides change
  only that file; `--wavecar` never relocates auxiliary inputs.
- Keep explicit relative CLI paths and output paths relative to the launch
  directory. INI `input_dir` defaults to the launch directory, while explicit
  relative INI paths stay relative to their configuration file.
- Retain mutually exclusive `--run-dir` chunk compatibility and the
  WAVEDER-only `optical_run_dir` INI alias. Standard INIs may omit `wavecar`.
- Validate native input directories consistently with GNU and Intel
  compilers. Preserve the input location when CI relocates an example INI.

### WAVEDER standard Kubo protocol

- Default supported Kubo charge calculations to standard VASP WAVEDER
  optical matrices. Require explicit `--kubo-source wavecar` for the
  canonical-momentum approximation and print an operator warning.
- Reject missing, invalid or unsupported WAVEDER input without silently
  changing the operator. Keep source choice and the actual operator in
  output metadata and calculation records.
- Support explicit single/range/list WAVEDER target bands, geometric native
  traces or resolvable per-band output, and Python μ/T selected contributions.
  Retain the full source intermediate-band space. Keep inferred native occupied
  and Python `--occupied N` insulating T=0 compatibility routes.
- Check every required weighted pair in either stored orientation. Permit
  complete producer-degenerate groups only with equal response weights; reject
  erased or missing pairs when needed. Never pad a missing block, truncate
  virtual bands to the target selection, or apply a second optical denominator.
- Label selected geometric integrals and charge contributions separately from
  total physical Hall, including finite-temperature subset calculations.
- Make the supported standard route explicit in native/Python help,
  examples, migration instructions and the technical report. Historical
  canonical results keep their original provenance and gain explicit
  approximation selections in reproduction commands.
- Update the native-spin CI driver and full-velocity/canonical comparison
  example to select the intended WAVECAR approximation explicitly; verify
  these executable callers as part of release preparation.
- Preserve explicit coverage/weight limits for selected metallic and
  finite-temperature responses, and separate spin-current operator contracts.
  Operator applicability and k-mesh/band convergence remain independent of
  source-file validation.

The maintainer selected the 1.6.5 patch-number increment for this release.
The changed default relative to published 1.6.3 means implicit WAVECAR Kubo
commands require the documented opt-in. Earlier release tags remain fixed.

## [1.6.3] - 2026-09-30

### Native Kubo band selection

- Remove native `--bundle` and `-kubo_bundle`; rejected uses explain how to
  migrate. Named `kubo`, `kubo-line` and `kubo-integral` tasks calculate one
  selected band or the trace of a multi-band range by default. Add
  `--per-band 1` for separate band curvatures.
- Default multi-band trace output to `KUBO.csv` when no explicit curvature
  path is supplied. Keep its numerical columns and the separate Python
  `bundle-hall` integration workflow.
- Reject unresolved individual-band degeneracies before producing output:
  every selected band must be separated from every other stored band by
  more than 1e-5 eV. Multi-band traces retain their external-gap check and
  permit internal degeneracies. Pure legacy `-kubo` invocations retain
  individual-band selection with the same stricter gap check.
- Update native walkthroughs, output descriptions and the technical report;
  preserve original result metadata and historical release instructions.

### WAVECAR compatibility and result validity

- Detect byte and legacy four-byte-word RECL layouts automatically, without
  modifying the WAVECAR. Reject malformed or ambiguous record layouts.
- Reject nonfinite Kubo inputs/results and incomplete coefficient reads.
  Write trace CSV schema V2 and individual-band CSV schema V3 with source
  dimensions, expected row counts and an initial INCOMPLETE status. Write
  terminal PASS only after every selected spin/k/band result is complete.
- Stage curvature and pair CSVs as `.partial` files and publish their final
  names only after successful output completion; preserve existing files.
- Validate finite values and complete row coverage in the Python readers.
  Retain older trace V1 and individual-band V2 files as legacy-unverified
  inputs; their original data and producer metadata are not rewritten.
- Use the maintainer-selected version 1.6.3 for this update; removed flags
  and changed multi-band defaults require the documented command migration.

## [1.6.2] - 2026-09-29

### Native command help

- Make `--help` and `-h` show a short overview. Add task and option help,
  including `--help task`, `--help kubo`, `--help bands` and
  `--help spin-chern`; retain the complete flag reference through
  `--help all` (also `--help legacy`).
- Add an Intel MPI walkthrough from task selection to option help and a
  native calculation. Help requires neither Python nor material input files.

### Maintenance

- Keep the completed v1.6.1 publisher manual so subsequent validation changes
  do not restart publication.

### Documentation and terminology

- Distinguish Chern numbers and spin Chern numbers from the
  Fukui-Hatsugai-Suzuki (FHS, also called Fukui) method used to calculate them, and describe
  Berry-curvature/response calculations in terms of the Kubo formula.
- Identify the separate Z₂ n-field construction as the Fukui-Hatsugai (FH)
  method. Update guides, examples, figure labels, help and the technical report.
  These editorial updates leave numerical algorithms, results and data formats
  unchanged.
- Add author affiliations, a present-address note and visible links to guides,
  calculation commands and data in the technical report. Its scientific content
  remains the 1.6.1 method baseline; saved results keep their original producers.

## [1.6.1] - 2026-09-27

### Spin Chern numbers

- Add native `--task spin-chern` for two-component WAVECARs. The selected
  isolated band subspace is split by the full projected-spin operator;
  positive and negative sectors use geometric determinant links. Report both
  sector Chern numbers, their sum and their half-difference.
- Use the complete pseudo-wavefunction Gram metric and the Cartesian spin
  frame from matching `OUTCAR`. Check energy isolation, the projected-spin
  gap, sector ranks and overlap conditioning before writing final results.
- Default native inputs to `WAVECAR` and `OUTCAR` in the execution directory;
  retain explicit file paths and independent energy/spin-gap tolerances.
- Add an intrinsic-SOC graphene example, numerical outputs, optional plotting
  and a technical-report discussion. This task is separate from
  conventional spin-current Hall response and PROCAR character weighting.

### Spin-sector Kubo curvature

- Add native `--task spin-kubo` for arbitrary k lists and paths. Compute all
  three Cartesian curvature components with the projected-spin derivative
  and consistent pseudo Gram normalization.
- Add `spin-kubo --sum-bands N` to test intermediate sums on unchanged VASP
  eigenstates, recording source and summed band counts separately and rejecting
  a cutoff through an unresolved external degeneracy.
- An explicit `--mesh NX,NY` adds a validated full-mesh approximation integral;
  retain its raw values and never round them to a Chern number. Report the
  canonical-momentum derivative-proxy approximation in every output.
- Save native point and spin-spectrum CSVs; keep Python optional for plotting.
  Reject unresolved energy/spin gaps, invalid metrics, changing ranks and
  nonfinite results. Preserve existing Kubo and Fukui task behavior.
- Add analytic-model geometry tests and real graphene/Bi examples comparing
  path curves, selected degenerate pairs and full-BZ Fukui invariants.
- Add controlled Bi mesh/source/retained-band comparisons and actual local
  graphene SOC sampling, with separate native curvature, geometric polygon
  flux, electronic-convergence checks and clearly labelled diagnostic fits.

### Compatibility and examples

- Retain existing native tasks, automatic spinor detection and postprocessing
  commands. The new tasks use matching WAVECAR/OUTCAR inputs and write separate
  CSV schemas; historical reference outputs retain their original metadata.
- Provide Intel Fortran/MPI commands, GNU alternatives, input preparation and
  saved-table figure replay. Native Fortran computes the new spin-sector
  results; the plotting scripts read completed numerical outputs.
- Keep software checks, geometric invariants and material convergence distinct.
  The graphene local-sampling controls do not establish a converged full-BZ
  Kubo integral, and projected-spin sectors are not conventional spin-current
  Hall conductivity or PROCAR-weighted charge attribution.

## [1.5.1] - 2026-09-27

### Automatic WAVECAR layout

- VASPBERRY and the WAVECAR-based Python tools detect scalar versus
  two-component wavefunctions from coefficient counts and the full plane-wave
  basis. All k points and stored spin channels must agree.
- Existing `-s 1|2`, `--spinor 1|2` and Python component arguments remain
  optional assertions. A conflicting value fails instead of being ignored.
- The INI `spin_mode` now defaults to `auto`. Ordinary scalar/spinor occupation
  multiplicities follow the detected source; collinear INI runs still require
  an explicit up/down channel choice. `soc` remains a compatibility alias for
  two-component storage, not evidence that SOC was enabled.
- Native logs report spinor components without inferring `LSORBIT`. Unsupported
  or inconsistent coefficient layouts fail before numerical work.
- Direct commands, hands-on examples and settings files omit redundant spinor
  arguments. Numerical formulas and saved-result schemas are unchanged.
- Accept the reviewed 1.5.1 Z₂ producer explicitly; unreviewed versions remain
  rejected. Publication requires the same thirteen exact-source validation jobs.

### Additional fixes and documentation

- Usage guides and CLI help explicitly name VASPBERRY execution, distinguish
  it from Python integration/plotting, and show direct Intel MPI commands for
  curvature CSVs and reusable pair data. Numerical behavior is unchanged.

- The hands-on guide now explicitly separates strict INI workflows from
  advanced coalesced calculations and distinguishes their saved-result formats.

- Postprocessing documentation now separates a first-run guide from the complete
  INI reference, with linked Bi rescan/region exercises, user-defined group and
  region names, and explicit numerical-output versus plotting instructions.

- Native MPI help now prints once from rank zero while all ranks finalize,
  preventing interleaved help text and intermittent multi-rank help checks.

- CI retries the checksum-pinned public Bi input fetch up to three times
  and retains each attempt's download/verification record, including failures
  before the numerical examples start. Persistent failures still fail CI.

## [1.5.0] - 2026-09-26

### Added

- A short `vaspberry_post.py check/run/plot` interface using one commented
  INI file for native execution, μ/T scans, optional PROCAR atom/orbital/spin
  groups and reciprocal-space regions. No manually written JSON is needed.
- Reuse of a completed run's validated pair cache for the same WAVECAR;
  saved-table plotting without original VASP inputs or the native executable.
- Source/mesh/spin checks, no-overwrite outputs, settings snapshots and
  command/status records. Existing numerical engines and tables are reused.
- A public full-mesh Bi example covering native pairs, insulating Hall,
  cache rescans and figures, plus focused parser/workflow regression tests.

### Documentation

- Lead repeated postprocessing with short run/plot commands; explain the
  editable settings, native raw data, numerical integration and drawing.
- Clarify initial Kubo contributions by Sun-Woo Kim and ongoing subsequent
  development by Hyun-Jung Kim in the README contributors section.
- Retain existing commands, numerical schemas and the scientific report.
  Fortran kernels and integration/projection formulas are unchanged.

## [1.4.2] - 2026-09-26

### Fixed

- Honor environment-selected compilers, MPI wrappers and launchers in Make;
  retain `gfortran` as the default instead of GNU Make's built-in `f77`.
- Rebuild GNU executables when their source choice, compiler command or build
  flags change, avoiding silent reuse of a previously configured binary.
- Remove incomplete GNU binaries after failed compilation and reject symlink
  build directories before writing or cleaning generated files.

### Build and documentation

- Add `make help`, Intel serial help checks and standard extra compiler/link
  flag controls while retaining byte-RECL, fixed-form and LP64 defaults.
- Document native-only prerequisites separately from Python postprocessing:
  compiler and numerical-library development packages, MPI SDK/runtime,
  shell/Make, matching launcher and runtime library paths.
- Add fresh-source-archive numerical installation checks for GNU/MPICH on
  Linux and GNU/OpenMPI/OpenBLAS on macOS ARM64, alongside the existing
  Linux GNU/OpenMPI and Intel ifx/ifort/Intel MPI checks.
- Clarify Intel Mac existing-toolchain support, current Homebrew limitations,
  WSL2 recipe scope, cluster modules, release archives and troubleshooting.
- Numerical kernels, output schemas and scientific reference data are
  unchanged. Rebuild after updating to use the corrected build behavior.
- Document the observed Ubuntu24.04 MPICH 4.2.0-5build3 PMIx/Hydra
  launcher incompatibility; validate MPICH on Ubuntu22.04 and retain
  Open MPI checks on Ubuntu24.04.

## [1.4.1] - 2026-09-26

### Build and validation

- Add Intel MPI runtime/help checks for `mpiifx` and `mpiifort`, retaining
  byte-based WAVECAR records and sequential LP64 oneMKL defaults.
- Validate actual MoS₂ curvature and pair export in Intel serial and
  two-rank MPI builds. Release publication requires these jobs in addition
  to the existing full CI and actual Bi serial/MPI checks on the same commit.

### Documentation

- Lead installation and native commands with Intel oneAPI + MPI; retain
  Intel Classic and GNU alternatives and matching-launcher instructions.
- Specify inputs, filenames, numerical columns/units and possible figures
  for native curvature, reusable pairs, Hall integration and PROCAR analysis.
- Explain why Python cache import, occupation-weighted integration and
  plotting are distinct, and document the existing `wavecar-hall` shortcut
  with saved intermediates and cache reuse.
- Preserve the v1.4.0 scientific report, examples and reference data.
  This compatible patch changes build checks and instructions, not the
  numerical kernels or output schemas; no scientific-result migration is needed.

## [1.4.0] - 2026-09-25

### Added

- General SOC PROCAR character postprocessing with named atom/layer/orbital
  groups, an explicit spin axis and the actual OUTCAR spin-frame rotation.
  Matching native pair caches support selected-band projected charge-Hall
  scans and reusable numerical outputs and plots. An analytic public fixture
  gives step-by-step commands; conventional spin-current Hall is a separate
  observable.
- Readable native Fortran task and input aliases, including `--task chern`,
  `z2`, `kubo` and `kubo-pairs`, with argument validation. Existing short
  options and numerical methods remain available.
- Full-program serial/MPI Fukui regression for curvature, Chern number and
  extrema metadata, exercised by CI alongside the native command tests.
- Canonical serial executable `build/vaspberry`, with the previous
  `build/vaspberry-gfortran` name retained for existing scripts.

### Fixed

- Export Kubo pairs from metallic or smeared states without requiring a fixed
  rounded occupation count across k points. Pair numerators do not use source
  occupations; fixed-subspace checks in other native modes are retained.
- Correct native canonical-momentum velocity output to m/s, fix extrema
  indexing after the k loop, and use scientific notation. Older velocity
  files may contain rounded zeros or invalid extrema and must be regenerated;
  this correction does not change the separate Kubo kernels.

- Reject incompatible occupied-bundle and spinor metadata before plotting
  native Z₂ fields. Numerical topology routines are unchanged.
- Populate Fukui MPI output-header curvature extrema and their k points from
  the gathered mesh, matching the printed numerical data.
- Avoid a zero divisor in the wavefunction progress display for real-space
  grids with fewer than five z points. The reconstruction formula is unchanged.

### Documentation

- Add a native command/output reference and a hands-on calculation → saved
  data → independent plotting tutorial, with explicit reusable Kubo scans.
- Prepare ordinary Bi WAVECAR inputs without the optional PAW spin producer;
  clarify which report examples need VASP regeneration versus saved results.
- Document generic postprocessing and PROCAR attribution, with explicit
  spin/layer/orbital-current scope and optional-producer limits for DFT+U.

- Pair direct-WAVECAR MoS₂ Z₂ = 0 and Bi Z₂ = 1 results in the technical
  report and feature examples, with native n-field data and reproduction commands.
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

Historical interface: the following Wilson-loop command belongs to this
earlier 1.1.0 release. `tools/wavecar_z2.py` is retired from the current
source. For version1.4.0 use native `--task z2` and the
[current Fukui–Hatsugai guide](docs/Z2_FUKUI_HATSUGAI.md); the old command and
output names below are retained as release history.

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
