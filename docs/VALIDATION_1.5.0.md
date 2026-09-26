# Validation of VASPBERRY 1.5.0

The [release publication](https://github.com/Infant83/VASPBERRY/releases/tag/v1.5.0)
records the exact source commit and hosted validation URLs. This release
adds a settings parser and orchestration front end; it reuses the established
Fortran, pair integration, PROCAR projection and plotting implementations.

## New interface checks

The parser tests cover unknown/duplicate settings, paths, spin conventions,
finite scans, regions, groups and plot selections. Workflow tests exercise
input rejection before native execution, failed-stage records, cache reuse
and plotting from saved results. Commands are passed as argument lists;
settings do not execute shell expressions.

The [public Bi example](../examples/features/simple-postprocess/README.md)
uses the real 12×12 SOC input with 18 stored bands. Its T=0 scan stays inside
the common band-10/band-11 gap. The output is the canonical-momentum Kubo
charge response, followed by a cached rescan and PNG/PDF/SVG plots. The fixed
input gives a small residual σ ≈ −6.27668×10⁻⁶ e²/h and Δσ = 0 across these
in-gap scans; CI uses that value as a regression reference. Time reversal
requires zero charge Hall in the converged physical calculation. This
workflow check does not establish material convergence or a nonzero effect.

During preparation, a separate matched 14×14, 40-band SOC dataset reproduced
the prior explicit commands exactly in five numerical NPZ products. Reusing
the cache at 100 K matched the corresponding saved scan rows exactly, with
no native export. Plotting copied numerical results produced all nine figures
with original VASP paths unavailable. That dataset is private and is not
distributed as a public fixture; the public example and automated fixtures
provide the reproducible release checks.

## Required release checks

All thirteen jobs must pass on the publication commit: the eight CI jobs
(Python 3.10/3.12, feature examples, Fortran guards, GNU11/13/Open MPI and
Intel serial), actual Bi serial/MPI Z₂, Intel ifx/ifort with Intel MPI and
LP64 oneMKL, and clean archive installations on Ubuntu22.04 GNU/MPICH and
macOS ARM64 GNU/OpenMPI/OpenBLAS. The feature-example job also runs the
new public settings workflow. Artifacts retain inputs, native results,
commands and validation receipts.

The [build guide](BUILD.md) distinguishes executed platforms from other
installation recipes. The retained [1.4.2 installation record](VALIDATION_1.4.2.md)
explains the environment coverage and MPICH packaging limitation. Native
source changes in this release are version labels only. Material-specific
k-mesh, band-window and operator convergence remain separate from software
and interface validation.
