# Changelog

All notable changes to VASPBERRY are recorded here.

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
- Corrected the sign and units in the legacy `NFIELD.dat` explanation.

### Added

- `Z2_FIELD.csv` on the fundamental mesh for PASS results, and
  `Z2_FIELD.invalid.csv` for INVALID or INCOMPLETE runs. Same-directory
  temporary output and C/POSIX `rename()` prevent partial publication;
  stale PASS files are removed before WAVECAR processing.
- Version, grid, band range, spinor rank, overlap backend, numerical
  thresholds, minimum link singular values, per-plaquette status, and an
  explicit check-scope statement in the field CSV.
- Production-linked Fortran regression tests for Γ/M1/M2/M3 reciprocal
  mapping, `Theta^2=-1`, complex*8 norm accumulation, LU pivot phase, and
  minimum singular values. Serial and MPI sources are both linked and run.
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
  `WAVECAR_PSEUDO_NO_PAW_AUGMENTATION`. PAW augmentation can affect
  quantitative finite-neighbor overlaps, but it is not the source of the
  omitted reciprocal-G shift.
- The pointwise Fukui n-field remains gauge- and logarithm-branch-dependent.
  The physical relation is the wrapped flux condition
  `wrap[flux(-k)+flux(k)]=0` modulo `2*pi`.
- Thresholds are conservative diagnostic policies and can yield false
  rejection on a coarse mesh. The guarded Wilson-loop workflow and mesh
  convergence remain required before reporting an invariant.

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
