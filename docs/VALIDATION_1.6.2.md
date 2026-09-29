# VASPBERRY 1.6.2 validation scope

This patch adds native topic help and improves documentation. Calculation
commands, numerical methods and output formats retain the 1.6.1 scientific
baseline. The [release notes](releases/v1.6.2.md) show the Intel MPI help and
calculation workflow.

## Exact-release evidence

This document specifies the required checks; it does not declare an untested
candidate successful. The publisher records the exact source commit and links
to successful validation runs in the
[v1.6.2 release body](https://github.com/Infant83/VASPBERRY/releases/tag/v1.6.2)
only after all required checks pass. Use those commit-specific records when
citing release validation. A local development check or a successful run on
another commit does not replace them.

The release requires **thirteen successful jobs on the same source commit**:

| Workflow | Required scope |
|---|---|
| [CI](../.github/workflows/ci.yml), eight jobs | Python 3.10/3.12 regression suites; public feature examples; actual-source Fortran Z₂ checks; GNU 11/13 serial/MPI portability; Intel ifx/ifort serial portability |
| [Bi n-field Z₂](../.github/workflows/bi-z2-validation.yml), one job | Actual Bi WAVECAR, serial and two-rank MPI fields against the retained material reference |
| [Intel MPI](../.github/workflows/intel-mpi-validation.yml), two jobs | ifx and ifort builds, help/runtime checks, actual MoS₂ serial/MPI curvature and pair exports, and native spin calculations |
| [Archive installation](../.github/workflows/installation-validation.yml), two jobs | Clean source archives on macOS ARM64 GNU/OpenMPI and Linux GNU/MPICH, including help and native numerical checks |

The [publisher](../.github/workflows/publish-v1.6.2.yml) checks workflow and job
identities, source SHA and run attempt. It rejects missing, failed or skipped
required jobs, a conflicting tag, and an advanced default branch. Existing
release tags remain fixed. Publication is followed by a fresh archive check
against the same source, including Git LFS payload identity where applicable.

## Native help checks

The [topic-help tests](../tests/test_native_help.py) exercise the short
overview, task and option pages, aliases, malformed requests and calls from
an empty directory. Help must terminate without opening WAVECAR or creating
calculation output. MPI help prints once from rank zero and all ranks finish.
The current and historical serial sources use the shared help include while
advertising only the calculations each source supports.

The [compiler help check](../tests/check_fortran_help.sh) checks the short
overview and complete flag reference on the built executables. Calculation regressions continue to run separately;
printing a correct help page is not evidence that a numerical kernel works.
The archive jobs verify that `vaspberry_help.inc` and the existing spin includes
are available when compiling without a Git checkout.

## Scientific scope and retained references

The [1.6.1 validation scope](VALIDATION_1.6.1.md) and
[technical report](TECHNICAL_REPORT.md) ([PDF](TECHNICAL_REPORT.pdf)) describe
the numerical methods, material examples and unresolved convergence questions.
This release retains the existing report Markdown and PDF, including their
1.6.1 method-baseline label. New command help does not alter its calculations
or saved results.

Ordinary WAVECAR Kubo calculations retain the canonical-momentum approximation.
Spin-sector Chern numbers, spin-sector Kubo curvature, conventional spin-current
response and PROCAR-weighted charge attribution remain distinct quantities.
This usability patch makes no new material-convergence claim.

Reference tables, figures and their input/output checksums keep their original
producer labels. Version 1.6.2 metadata identifies newly executed outputs; it
does not relabel older results as new calculations. Record the executable
version, exact source, compiler, inputs, sampling and band selections when
reproducing a calculation.
