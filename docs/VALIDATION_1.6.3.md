# VASPBERRY 1.6.3 validation scope

This release changes native Kubo selection and output-completion checks,
adds automatic byte/word RECL detection, and updates readers and public
examples. The [release notes](releases/v1.6.3.md) give a direct Intel MPI
calculation and the command migration.

## Exact-release evidence

This document defines required validation; it does not declare an untested
candidate successful. The [release body](https://github.com/Infant83/VASPBERRY/releases/tag/v1.6.3)
records the exact source commit and successful run links only after the
publisher accepts all required checks. Local development results and tests
on earlier commits do not substitute for those final records.

The release requires **thirteen successful jobs on the same source commit**:

| Workflow | Required scope |
|---|---|
| [CI](../.github/workflows/ci.yml), eight jobs | Python 3.10/3.12 regression suites; public feature examples; actual-source Z₂ checks; GNU 11/13 serial/MPI and Intel ifx/ifort serial portability |
| [Bi Z₂](../.github/workflows/bi-z2-validation.yml), one job | Actual public Bi WAVECAR, serial and two-rank MPI fields compared with the retained reference |
| [Intel MPI](../.github/workflows/intel-mpi-validation.yml), two jobs | ifx/ifort help/runtime, actual MoS₂ serial/MPI curvature and pair data, layout and new CLI checks, native spin validation |
| [Archive installation](../.github/workflows/installation-validation.yml), two jobs | Fresh source archives on macOS ARM64 GNU/OpenMPI and Linux GNU/MPICH; native build/runtime and numerical checks without `.git` |

The [publisher](../.github/workflows/publish-v1.6.3.yml) verifies exact workflow,
job, source and attempt identities. Missing, failed or skipped required jobs,
an advanced default branch or a conflicting tag prevent publication. Earlier
tags remain fixed. After publication, check the actual source archive against
the release commit, including expanded Git LFS payloads, then build and replay
representative public commands from that archive.

## Native and reader checks

- **Band selection:** singleton and `N:N` selections, multi-band trace,
  explicit per-band output, named-task aliases and pure legacy selection.
  Removed flags must give actionable errors. Invalid single/per-band gaps
  must fail before creating a completed result.
- **Input layout:** byte and four-byte-word RECL representations must recover
  the same physical input and numerical result. Malformed, truncated or
  ambiguous layouts must be rejected; detection must not modify the source.
- **Numerical validity:** reject nonfinite coefficients and results, incomplete
  reads and failed output writes. A later spin/channel failure must not leave
  a final CSV marked PASS. Partial-file staging and final publication must
  preserve existing files and reject failed flush/close operations.
- **CSV completion:** trace V2 and band V3 require finite values, complete
  unique row coverage consistent with source dimensions, and terminal PASS.
  Readers must reject truncated, duplicated or corrupted output. Older
  schemas remain explicitly legacy-unverified rather than acquiring a new
  completion guarantee.
- **Independent numerical checks:** compare serial/MPI occupied MoS₂ trace
  and pair exports, and reconstruct the trace from the external-pair sum.
  The shared [Intel/archive harness](../tests/run_intel_mpi_validation.py)
  also exercises the new default-output and input/selection behavior.

Native help requires no material input or Python and must terminate cleanly
on all ranks. Source archives must contain all shared native includes.
Full Python suites also cover postprocessing and plotting interfaces; help
and file-format checks alone are not numerical validation.

## Scientific scope and retained references

The [technical report](TECHNICAL_REPORT.md) ([PDF](TECHNICAL_REPORT.pdf)) and
[report-to-example map](../examples/REPORT_REPRODUCTION.md) identify each
public dataset, its regeneration requirements and its numerical limits.
The supplied 48-point path helper now computes occupied bands 1–18. Its
historical band-17/18 references represent a different quantity and remain
unchanged; the guide does not claim the new trace reproduces those curves.

Canonical momentum remains the WAVECAR Kubo operator. The release does not
add missing PAW/nonlocal/SOC velocity terms or establish new material
convergence. Geometric spin-sector results, conventional spin-current Hall
response and PROCAR-weighted charge attribution remain distinct quantities.
Original tables, figures and producer labels are preserved. New outputs
identify version 1.6.3; record the exact commit, executable, compiler, inputs,
sampling and band selection alongside each calculation.
