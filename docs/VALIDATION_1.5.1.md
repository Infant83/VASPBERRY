# VASPBERRY 1.5.1 validation scope

Release 1.5.1 was validated by thirteen successful jobs on its exact
[source commit](https://github.com/Infant83/VASPBERRY/commit/a9b2aa15f9b1d89d7da1fc10889a62bbf026059a):
Python 3.10/3.12, GNU and Intel serial builds, native Z₂ regressions, public
feature examples, actual Bi serial/MPI comparison, Intel MPI with ifx and
ifort, and clean Linux/macOS source-archive installations. Its release body
records those runs. The 1.5.1 tag and source remain unchanged after withdrawal
of the duplicate 1.6.0 release/tag; the one-time publishers are retired.

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
