# VASPBERRY 1.6.5 validation scope

Version 1.6.5 supersedes the unpublished 1.6.4 candidate. This record
separates prior implementation evidence, checks executed for the release,
and the exact-commit remote checks required before publication. Existing
material reference data are not regenerated or relabeled by a version change.

## Prior implementation evidence

The [1.6.4 candidate record](VALIDATION_1.6.4.md) retains the original run
history. Its final selected-WAVEDER implementation passed **998 tests** on
macOS with GNU Fortran/Open MPI, an independent Python/native review, and a
clean source-archive build and example replay. These are candidate-baseline
results, not a claim that the later 1.6.5 commit already passed remote CI.

The independent arithmetic checks covered all selections of small full and
rectangular matrices, source-cluster and missing-pair rejection barriers,
and standard VASP 5.4.4 optical data. On the actual full 6×6 MnBi₂Te₄ source,
36 selected μ/T rows matched a direct band-column sum within
7.5×10⁻¹⁴ e²/h. This verifies contraction and integration on those inputs;
it does not establish mesh or source-band convergence of that material.

## Release checks

The release preparation rechecks the current README, CLI help, native/Python
commands, INI workflow, selected-band example, source comments and technical
report against the implemented operator contract. The executable synthetic
example records each command, output, exit status and expected rejection,
and compares its 27 μ/T rows with an independent full-matrix expression.

On 2026-10-04, the rebuilt 1.6.5 GNU executable and Python tools passed the
selected-band example: **27 μ/T rows**, maximum direct-sum discrepancy
**1.39×10⁻¹⁷ e²/h**, native full/rectangular trace agreement, per-band output,
and three expected rejection cases. Nine further public help/calculation
commands passed, including disjoint `--bands 1:2,4`, six weighted μ/T points,
and native mesh integration. Its plot records the selected contribution and
band IDs in both the visible title and checksum-bound metadata.

The release-wide local suite passed **1,002 tests with no skips** on macOS.
GNU serial/MPI build and runtime checks passed. The real standard optical
6×6 selected μ/T calculation was repeated and matched the independent
direct-column reference within 7.47×10⁻¹⁴ e²/h; the existing occupied-space
serial/MPI results were unchanged. The 34-page PDF was rendered again:
four changed pages passed visual inspection, and the other thirty page
renders matched the previously reviewed report.

The first remote candidate exposed a compiler-dependent directory check
in the Intel MPI runs and an input-directory override missing from the
relocated-INI example harness. Both were corrected before publication.
The final native directory probe passed 25 focused CLI tests, and the
expanded serial/MPI validation driver passed locally with GNU on actual
MoS₂ data. The relocated public Bi example also passed calculation, cached
rescan and plot generation. Final Intel results are required from the
exact-commit remote checks below.

Exact-commit remote receipts are recorded in the publication body after
the required checks finish. No unexecuted compiler or workflow job is
counted as passed here.

## Required publication checks

The manual `publish-v1.6.5.yml` workflow requires all thirteen jobs from the
four validation workflows on its exact source commit: full tests, actual
Bi validation, Intel MPI variants, and clean Linux/macOS source-archive
installation. It verifies the previous immutable v1.6.3 commit and refuses
a conflicting tag or an advanced default branch. The published release
body records the checked source commit and successful workflow receipts.
A fresh version-pinned source archive is checked again after publication.

## Scientific scope

- Native explicit `--bands` selects geometric band/bundle curvature; its
  mesh integral uses unit weights and is labeled as a selected result.
- Python `--bands` integrates the selected Fermi-weighted μ/T contribution.
  All source NBANDS remain in the virtual sum. A selected result is not
  automatically the total AHC or total ΔAHC.
- The complete insulating occupied T=0 route remains available. Required
  pairs must exist and survive the producer's degeneracy treatment;
  missing information is rejected rather than padded or silently omitted.
- Software/source validation does not certify complete SOC/+U physical
  velocity terms, a new material's mesh/NBANDS convergence, or equality of
  distinct operators. Existing coarse optical and canonical examples keep
  their original scientific limits.

The [standard protocol](WAVEDER_KUBO_PROTOCOL.md) and
[reproducible selected-band example](../examples/features/waveder-selected/)
define the supported procedure.
