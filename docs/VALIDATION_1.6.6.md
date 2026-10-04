# VASPBERRY 1.6.6 validation scope

Version 1.6.6 packages the command-validation and report corrections made
after published 1.6.5. This record distinguishes the completed numerical
reproduction from exact-commit release checks. Updating the software version
does not relabel historical material producers or establish convergence.

## Numerical report reproduction on 2026-10-04

The VASPBERRY and associated analysis calculations underlying Figures 1–16
and Section 3.6 were repeated from retained electronic-state, operator and
fixed Wannier-model inputs. The five numerical result tables were
recalculated, and the task/output and material/settings tables were checked
against the executed workflows and input records. The run used the 1.6.5
source with the subsequent CLI/documentation corrections now in 1.6.6.

- All numerical values in 42 MoS₂ CSV tables (189,333 rows) and 63 native
  topology/spin tables matched the references exactly. Fresh native pairs
  covered every MoS₂ mesh/cutoff case; all 34 Bi convergence rows were checked.
- Actual four-stacking WAVEDER spectra, raw wavefunction reconstruction,
  independent topology diagnostics and public example workflows were rerun.
  Synthetic four-state, PROCAR and selected-WAVEDER checks are recorded
  separately from material data.
- Bi physical-matrix charge/spin integrations, model bands and strip states
  were repeated. Original operator streams were independently decoded.
- All five MnBi₂Te₄ fixed-model quadratures were rerun. One independent
  postw90 quadrature agreed within 7.37×10⁻⁹ e²/h; other external controls
  retain historical provenance. The actual complete-mesh native Fukui
  check gave occupied Chern −1 and deep-space Chern 0. The coarse direct
  optical response was reproduced without implying convergence.
- Eight supplementary Bi conditioning cells were corrected to the matching
  48-band source. Principal Chern/Kubo values were unchanged. Two report
  rounding descriptions and a clipped operator label were corrected.

The [report Appendix C](TECHNICAL_REPORT.md#appendix-c-numerical-reproduction-check-of-2026-10-04)
and [reproduction guide](../examples/REPORT_REPRODUCTION.md) give the scope,
input availability and commands. No fresh VASP SCF/NSCF states, producer run
or Wannier localization were generated. The published unconverged cases
remain unconverged; reproducing them is not a convergence claim.

## Release source checks

The preceding CLI-consistency audit passed 1,013 local tests. Those results
belong to that audit source; they do not replace validation of the final
versioned release commit. Report-specific checks additionally exercised the
four-state model, revised plot labels, all eight catalogue features, INI
calculation/rescan/region workflows and selected-WAVEDER/PROCAR examples.

Final release preparation rebuilds the versioned executable and reruns the
appropriate local tests and examples. Exact-commit remote checks and the
publication receipt are recorded in the published release body; an
unexecuted compiler or workflow job is never counted as passed here.

## Publication requirements

The manual `publish-v1.6.6.yml` workflow requires thirteen successful jobs
from four workflows on the exact source commit: full tests, actual Bi,
both Intel MPI compiler variants, and clean Linux/macOS archive installs.
It rechecks the latest attempts, preserves the immutable v1.6.5 commit
`511aac3c16d787ac37979a4b49b37dfa5d92e37a`, and refuses a conflicting tag,
an advanced default branch or modification of an existing release.
A fresh version-pinned source archive is checked after publication.

## Scientific boundaries

The [WAVEDER protocol](WAVEDER_KUBO_PROTOCOL.md) and
[selected-band example](../examples/features/waveder-selected/) retain
separate geometric-trace, selected Fermi-weighted contribution and complete
occupied insulating total-Hall contracts. Missing required pairs or
unsupported producer-degenerate weights are rejected. Source validation
neither certifies complete SOC/+U velocity terms nor proves k-mesh/NBANDS
convergence or equality between distinct operators.
