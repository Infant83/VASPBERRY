# Independent numerical diagnostics

These records check actual VASP-derived inputs and VASPBERRY outputs. They
are case-specific validation snapshots; the [material guide](../../) states
the current interpretation and convergence limits.

- `actual-6x6-independent-audit.json` and `actual-12x12-independent-audit.json`:
  direct decoding of the 48-band raw PAW operator streams and independent
  summation of all 27 spin-conductivity and 9 charge-Hall tensor components.
- `actual-6x6-b64-independent-audit.json`: the same check for the 64-band
  source. Its highest states remain an eigensolver-convergence diagnostic.
- `actual-three-case-independent-audit.json`: comparison of those source
  calculations and retained-band windows. Implementation agreement does
  not establish k-mesh or current-product convergence.
- `actual-6x6-b48-b64-cutoff-diagnostic.json`: measured doublet splittings,
  retained-boundary gaps and changes between source band counts.
- `occupied-tr-subspace-12x12.json`: independent time-reversal closure of
  the 10 occupied raw-pseudo WAVECAR states, with explicit spin and reciprocal
  basis mapping. The native n-field table comes from the identical WAVECAR.

- `actual-edge-width-independent-audit.json`: width 20/40 spectra, finite-width
  Γ splitting, sampled common bulk gap and odd edge-crossing checks. Boundary
  localization uses eigenvalues of the edge projector in each crossing
  subspace, so it is independent of degenerate-eigenvector choices.

- `actual-6x6-b80-retained-independent-audit.json` and
  `actual-12x12-b64-retained-independent-audit.json`: independent full-tensor
  checks in explicitly retained spaces of padded VASP calculations. The
  inaccurate highest source states are excluded.
- `padded-source-independent-summary.json`: finite-band and same-cutoff
  eigensolver stability controls, with retained doublet and boundary gaps.

The spin-current oracle uses the same exported physical operators with an
independent contraction, so this verifies calculation and data association,
not an independent electronic-structure approximation. The time-reversal
check compares raw-pseudo coefficient subspaces; it is not a PAW-metric
fidelity. The physical spin Hall response includes its separately validated
PAW spin and full-velocity corrections. No source matrices or wavefunctions
are repaired, and no result is rounded to an integer.
