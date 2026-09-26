# Validation of VASPBERRY 1.4.2

The exact source commit and hosted validation URLs are recorded in the
[release publication](https://github.com/Infant83/VASPBERRY/releases/tag/v1.4.2).
This patch changes build configuration handling and installation instructions;
numerical source changes are limited to version labels.

## Required release checks

The existing eight CI jobs cover Python 3.10/3.12, feature examples,
Fortran guards, GNU11/13 with Open MPI and Intel serial compiler portability.
The separate actual Bi job compares serial/MPI n-fields and Z₂=1.
Two Intel MPI jobs use real ifx/ifort, Intel MPI and sequential LP64 oneMKL
for runtime/help and actual MoS₂ serial/MPI curvature/pair comparisons.

The new [installation workflow](../.github/workflows/installation-validation.yml)
adds two jobs. Each extracts a fresh archive of the exact source commit
and builds without the checkout's generated files or Git metadata:

- Ubuntu22.04, GNU Fortran + MPICH + LP64 system BLAS/LAPACK.
- macOS15 ARM64, Homebrew GNU Fortran + Open MPI + LP64 OpenBLAS.

These jobs check serial help, two-rank communication/help, a known-value
complex ZGESVD call through the LP64 LAPACK interface, and actual MoS₂
bundle curvature (48 rows) and all stored band pairs (23,808 rows).
Serial/MPI outputs and an independently assembled external-pair sum are
compared within the recorded tolerances. Archive, source, input and output
identities and compiler/library versions are retained in job artifacts.
All thirteen jobs must pass on the release commit before the tag is created.

## Build-behavior checks and platform limits

Focused regression tests cover environment and command-line overrides,
configuration-dependent rebuilds, unchanged incremental builds and parallel
Make behavior, preserving the protected build directory contract.
Historical numerical references and the scientific report remain unchanged.

Intel Linux, GNU Linux and macOS ARM64 have hosted checks as listed above.
The preparation environment also permits a clean archive check on an existing
Intel Mac GNU/OpenMPI stack; that does not establish a fresh current Homebrew
installation on Intel Macs. Native Windows is unsupported; WSL2 is a Linux
installation recipe, not a separately executed hosted test. Other compiler,
MPI or OS combinations need site validation. See the [build guide](BUILD.md).

An initial Ubuntu24.04 MPICH archive job failed the independent two-rank
communication check: the distribution's MPICH 4.2.0-5build3 package launched
two singleton ranks. The [upstream diagnosis](https://github.com/pmodels/mpich/issues/7064#issuecomment-2301026290)
identifies incompatible PMIx/Hydra packaging. The MPICH job therefore uses
Ubuntu22.04; Ubuntu24.04 remains covered by the Open MPI jobs. The build guide
records the affected package and the required runtime check.

Software installation checks do not replace material-specific mesh, source
band and physical-operator convergence.
