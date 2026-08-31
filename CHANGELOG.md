# Changelog

All notable changes to VASPBERRY are recorded here.

## [1.1.1] - 2026-08-31

### Fixed

- Corrected spinor time-reversal reconstruction at nonzero TRIMs by mapping
  reciprocal-G labels with
  `G_target = round(-k_source-k_target)-G_source`. The previous shift-free
  `-G` rule was valid at Gamma but not at the three folded M points.
- Replaced coefficient-position assumptions with an explicit source-to-target
  G-vector bijection and hard checks for missing, duplicate, noninteger, or
  norm-changing mappings.
- Separated the physical Berry flux from the gauge- and logarithm-branch-
  dependent Fukui integer n-field, using double-precision pi and guarded link
  determinants.
- Restricted MPI receive-buffer use to rank zero and reduced the new flux
  array explicitly before validation.

### Added

- `Z2_FIELD.csv` on the fundamental mesh, with reciprocal-cell centers,
  time-reversal partners, Berry flux and curvature, n-field values, and
  numerical residuals.
- Hard validation for partner-map bijection/involution, odd Berry flux under
  time reversal, integer n-field residual, phase margin, zero total Chern
  number, balanced half zones, and matching half-zone parity.
- Synthetic Gamma/M1/M2/M3 regression tests that distinguish the exact
  folded-TRIM mapping from the former shift-free reconstruction.

### Scientific scope

- The direct WAVECAR overlap backend is labelled
  `WAVECAR_PSEUDO_NO_PAW_AUGMENTATION`. Missing PAW augmentation is a
  quantitative limitation, but it is not the identified cause of the
  M-point-only asymmetry fixed here.
- Pointwise n-field symmetry is diagnostic rather than a hard physical gate;
  the physical time-reversal test is
  `berry_flux(-k) = -berry_flux(k)`.
- The Fortran result remains a legacy Z2 candidate. A full-mesh WAVECAR,
  convergence checks, and the guarded Wilson-loop workflow are still required
  before reporting an invariant.

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
