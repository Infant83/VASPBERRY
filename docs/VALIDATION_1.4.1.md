# Validation of VASPBERRY 1.4.1

The immutable release's source SHA and successful hosted checks are linked
in its [publication record](https://github.com/Infant83/VASPBERRY/releases/tag/v1.4.1).

This patch changes build checks and documentation. Fortran numerical source
changes are limited to version labels; Python calculation modules also only
change producer-version labels. Output schemas, scientific references and
the v1.4.0 technical-report PDF remain unchanged.

## Checks required before publication

- Full CI: Python 3.10 and 3.12 regression suites, feature examples,
  Fortran guards, GNU 11/13 serial/MPI portability and Intel serial builds.
- Actual Bi: independently regenerated serial and two-rank MPI n-fields,
  compared with the retained v1.2.0 golden result, Z₂ = 1.
- Intel MPI: separate ifx and ifort jobs build default Intel MPI/oneMKL
  targets, check the two-rank runtime and native help, then run actual
  MoS₂ occupied-bundle curvature and all-band pair export in serial and
  MPI. Saved outputs are compared numerically and against independent or
  stored curvature references. Compiler versions and artifacts are retained
  by [the workflow](../.github/workflows/intel-mpi-validation.yml).

All eleven jobs must succeed on the publication commit. The publisher
rejects missing, skipped or failed jobs and checks source SHA, workflow,
run attempt, unchanged previous tag and current default-branch head.
Historical material-convergence studies retain their original provenance.

## Reproduction and interpretation

The guides separate native point curvature from raw pair numerators,
cache import from Hall integration, and integration from plotting. The
existing `wavecar-hall` shortcut preserves all raw/cached/final files;
manual `pair-hall` scans use the same cached input. No Python plotting
step is needed to obtain native Fortran data.

These are software and workflow checks, not new material convergence.
The direct WAVECAR route retains the canonical-momentum approximation;
PROCAR gives state character and selected-band charge-Hall attribution.
K-mesh, band-window and physical-operator convergence must be assessed
for each material. See [v1.4.0 validation](VALIDATION_1.4.0.md) for the
retained scientific example evidence and limits.
