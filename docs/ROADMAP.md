# Roadmap

## Python library interface

A packaged Python port of the Z2 implementation is a future development
direction. The intended result is an importable package with a stable array-level
API, for example a `vaspberry.z2` module, while retaining the Fortran command
line as the reference implementation during migration.

Before publishing that API, the work should be split into independently tested
layers: WAVECAR metadata and coefficient I/O, reciprocal-G/time-reversal
mapping, overlap and n-field kernels, validation policy, and plotting/output.
Each layer needs golden comparisons against the serial and MPI Fortran results,
including invalid-input and near-threshold cases. Packaging, semantic versioning,
typed return objects, and a documented exception model should follow only after
the numerical equivalence tests pass.

No Python library compatibility or release date is promised by the present
source tree. The existing Python scripts remain command-line post-processing
tools until this separation and equivalence work is complete.

Version 1.3.0 adds general Kubo/charge-Hall command-line tools and a versioned
point-curvature file format. Their internal modules are not a separately
packaged stable Python library, and they do not replace the Fortran Z2 kernel.
See [Kubo transport](KUBO_TRANSPORT.md) and [output format](OUTPUT_FORMAT.md).
