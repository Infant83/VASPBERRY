# VASPBERRY 1.6.0 validation scope

The release publisher requires thirteen successful validation jobs on the
exact source commit: Python 3.10/3.12, GNU and Intel serial builds, native Z₂
regressions, public feature examples, actual Bi serial/MPI comparison, Intel
MPI with ifx and ifort, and clean Linux/macOS source-archive installations.
The GitHub release body records those runs and the immutable source commit.

Automatic-layout regressions cover scalar, two-component and collinear
WAVECARs, omitted options, explicit assertions and conflicts, inconsistent
k points/channels, unsupported layouts, MPI behavior, inferred occupation
multiplicity, and INI/PROCAR/cache-reuse workflows. Old explicit commands
remain part of the compatibility checks.

Native MoS₂ curvature and Bi pair exports are compared with the existing
numerical references. The public Bi INI workflow checks end-to-end export,
integration and saved-table plots with automatic spinor detection.

These checks validate software behavior for supported full-complex WAVECARs.
They do not establish material convergence or infer SOC, magnetization direction,
or PROCAR spin-frame settings from the coefficient layout. Historical scientific
reference outputs retain their original inputs and provenance.
