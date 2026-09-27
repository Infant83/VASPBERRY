# Versions and source releases

## Choose the source for a calculation

`master` is the repository's default branch and contains the latest merged
source. Clone it with `git clone https://github.com/Infant83/VASPBERRY.git`.
Update a clean checkout with `git pull --ff-only`.

An active release tag, such as `v1.5.1`, stays fixed. Use the tag when reproducing a
published calculation, and record the version plus `git rev-parse HEAD` in
the calculation record. `git describe --tags --always --dirty` also identifies
changes after a release. A development branch is not the default download.

[GitHub Releases](https://github.com/Infant83/VASPBERRY/releases) provide
versioned source archives and notes. Executables are built locally; large
Git LFS inputs may need the [input-fetch procedure](../examples/INPUTS.md).

## Version policy

- Patch releases, for example 1.5.1, contain compatible fixes, including input
  detection and simpler use of existing commands. Numerical bug fixes must
  describe changed results and any migration needed.
- Minor releases, for example 1.4.0, add compatible features and commands.
- Major releases change public command or data contracts incompatibly.

Agree on the version with the maintainer before publishing. Do not create an
additional release merely to renumber the same update. Version 1.6.0 was
withdrawn and must not be reused; a future 1.6 minor release must use 1.6.1
or later.

`VERSION` is the authoritative software version. CLI/Fortran version labels,
`CITATION.cff`, release notes and the changelog must agree. A release records
its date and a fixed tag; development work remains under an Unreleased
changelog entry until publication. The CFF format version is separate from
the software version. A historical DOI must not be reused as a new version DOI.

The automatic-spinor update uses **1.5.1**. At the maintainer's explicit
request, the duplicate 1.6.0 release and tag were withdrawn. Its
[source commit](https://github.com/Infant83/VASPBERRY/commit/eb0b468aa0674f5ec579e75a014fe4591ec4cdfe)
and [withdrawal record](releases/v1.6.0.md) preserve provenance.

Ordinarily, published tags remain fixed: correct numerical results in a new
release without moving the old tag. An explicit maintainer-requested withdrawal
is an exception, recorded with the reason and original commit; do not reuse its
version number or recreate its release/tag. The active 1.5.1 release and tag
are unchanged. Software validation and material convergence remain separate.

## Publication checks

For a release, run the relevant unit, compiler/MPI and actual-input checks,
review the source and reference changes, and verify the documentation links.
Publish only the checked commit. Record the CI runs and input/operator scope
with the release; test a fresh version-pinned checkout after publication.

The one-time `publish-v1.3.0.yml` workflow waits for the release commit's CI
and verifies that its scientific files match the independently tested PR26
merge. That merge's actual Bi serial/MPI validation is retained as source
validation. Only the explicitly listed release documentation, metadata,
metadata tests and publication files may differ. The workflow refuses a tag
that points elsewhere and never edits an existing published release.

Future releases need their own target and validation review. The 1.3.0
workflow must not be used to label changed scientific code as already tested.

The `publish-v1.4.0.yml` workflow validates versioned metadata and waits for
both the full CI and actual Bi serial/MPI workflow on its own exact commit.
It requires successful jobs, checks that master has not advanced and that
the previous tag remains unchanged, and refuses to move a conflicting tag
or edit an existing release. Its publication artifact and release body record
the source commit and successful checks. Follow publication with a fresh
version-pinned checkout and the public hands-on/example commands.

The `publish-v1.4.1.yml` workflow additionally requires both actual
Intel MPI jobs (`Intel MPI (ifx)` and `Intel MPI (ifort)`) on its exact commit.
All eleven validation jobs must pass before publication. The v1.4.0 tag and
its publication workflow remain unchanged. The v1.4.0 technical report is
the retained method edition; current build and command instructions are in
the v1.4.1 guides.

The `publish-v1.4.2.yml` workflow retains all eleven v1.4.1 validation jobs
and adds two clean source-archive installation jobs: Linux GNU/MPICH and
macOS ARM64 GNU/OpenMPI/OpenBLAS. All thirteen jobs must pass on the exact
publication commit. Previous tags and publication scripts remain immutable.

The `publish-v1.5.0.yml` workflow requires the same thirteen exact-source
jobs, with the new parser/workflow tests and public settings example in CI.
It also requires the postprocessing command, guide and example files, and
preserves the immutable v1.4.2 tag. The scientific report remains the
retained method edition; new command instructions are in the current guides.

The former 1.6.0 and 1.5.1 publishers each required thirteen successful jobs
on their exact release commit, including Bi, both Intel MPI builds and
Linux/macOS archive installations. They have been retired so an old workflow
cannot recreate the withdrawn tag or require it for future publication.
Current downloads use 1.5.1; its original source and validation record remain
unchanged. The withdrawal creates no new software version.
