# VASPBERRY 1.5.1 validation scope

The patch-numbered release retains the automatic-spinor functionality published
as 1.6.0. Its publisher requires thirteen successful validation jobs on the new
exact source commit: Python 3.10/3.12, GNU and Intel serial builds, native Z₂
regressions, public feature examples, actual Bi serial/MPI comparison, Intel
MPI with ifx and ifort, and clean Linux/macOS source-archive installations.
The release body records those runs and the immutable source commit; both
published v1.5.0 and v1.6.0 tags remain fixed.

Automatic-layout checks cover scalar, two-component and collinear WAVECARs,
omitted options, explicit assertions/conflicts, inconsistent k points/channels,
unsupported layouts, MPI behavior and INI/PROCAR/cache reuse. Kubo Hall checks
cover inferred occupation multiplicity while existing one-channel conventions
remain unchanged. The Z₂ comparator accepts 1.5.1 explicitly and still rejects
unreviewed producer versions.

Native MoS₂ curvature and Bi pair exports retain the numerical reference checks.
The public Bi INI workflow covers export, integration, cache rescans and
saved-table plots with automatic detection.

These checks validate supported software behavior, not material convergence.
The coefficient layout does not establish SOC or the PROCAR spin frame.
Historical scientific results and the retained technical-report edition keep
their original inputs and provenance.
