# Versions and source releases

## Choose the source for a calculation

`master` is the repository's default branch and contains the latest merged
source. Clone it with `git clone https://github.com/Infant83/VASPBERRY.git`.
Update a clean checkout with `git pull --ff-only`.

A release tag, such as `v1.4.0`, stays fixed. Use the tag when reproducing a
published calculation, and record the version plus `git rev-parse HEAD` in
the calculation record. `git describe --tags --always --dirty` also identifies
changes after a release. A development branch is not the default download.

[GitHub Releases](https://github.com/Infant83/VASPBERRY/releases) provide
versioned source archives and notes. Executables are built locally; large
Git LFS inputs may need the [input-fetch procedure](../examples/INPUTS.md).

## Version policy

- Patch releases, for example 1.3.1, contain compatible fixes. Numerical bug
  fixes must describe changed results and any migration needed.
- Minor releases, for example 1.4.0, add compatible features and commands.
- Major releases change public command or data contracts incompatibly.

`VERSION` is the authoritative software version. CLI/Fortran version labels,
`CITATION.cff`, release notes and the changelog must agree. A release records
its date and a fixed tag; development work remains under an Unreleased
changelog entry until publication. The CFF format version is separate from
the software version. A historical DOI must not be reused as a new version DOI.

Never move an existing release tag to fix a result. Publish the correction
under a new version, preserve the earlier reference files, and explain the
change. Software validation and material-specific convergence are separate;
release notes retain the scientific scope and limitations.

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
