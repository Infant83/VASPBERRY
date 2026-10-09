# Install VASPBERRY from its canonical release

The plugin contains workflow instructions and helpers. The engine remains at
https://github.com/Infant83/VASPBERRY; do not copy it into the plugin.

## Resolve a new installation or requested update

Preserve an explicitly requested version or supplied checkout. When installing
without a version request:

1. Read `https://api.github.com/repos/Infant83/VASPBERRY/releases/latest`, or
   inspect `https://github.com/Infant83/VASPBERRY/releases/latest` with available
   web tools. Select its published, non-draft, non-prerelease `tag_name`.
2. Resolve the observed tag to its exact commit, for example through
   `https://api.github.com/repos/Infant83/VASPBERRY/commits/<tag_name>` with a
   URL-encoded tag. Record release URL, tag, full commit and lookup time. If
   lookup fails, say that the latest release could not be confirmed. Do not
   silently substitute the default branch or call a remembered version latest.
3. Announce the selected release. Clone that tag into a new versioned directory
   with `GIT_LFS_SKIP_SMUDGE=1 git clone --branch <observed-tag> --depth 1
   https://github.com/Infant83/VASPBERRY.git VASPBERRY-<version>-<commit-prefix>`, substituting
   verified values and checking that the destination does not already exist.
   Include the observed commit prefix so a same-version amendment does not
   collide with an earlier checkout. Use a separate new environment as well.
4. Compare checkout HEAD with the resolved commit and check VERSION/CITATION.cff.
   Stop on an identity mismatch. Read that release's README, build guide,
   changelog and CLI help; run its applicable validation and record actual
   source identity in the result. The old plugin baseline is not proof that
   a newly published version has been tested by this plugin.

For an update request, repeat this lookup and install in a separate versioned
directory/environment. Preserve earlier runs and modified checkouts. If an
existing checkout is suitable and no update was requested, use it and report
its version. This workflow runs when the agent is invoked; there is no separate
update daemon or second engine repository to maintain.

If the helper rejects a future metadata/interface change, inspect the release's
own metadata and supported CLI; do not invent a citation or bypass input checks.

## Last validated release and reproducibility

The plugin was checked against the 2026-10-09 help amendment of v1.6.6,
commit `2e7067b8b8b6d3ab444eef5bc025f4ce5720ca4d`. The upstream v1.6.6 tag moved
from the original `de84ed3337b900b36990f0230e5d497ff989061f` when that amendment
was published. Reproduce a recorded analysis by its exact commit, not just the
version or a mutable tag. For a new installation use the release lookup above;
the following commands reproduce this plugin's checked source:

```bash
git init VASPBERRY-1.6.6-2e7067b
git -C VASPBERRY-1.6.6-2e7067b remote add origin https://github.com/Infant83/VASPBERRY.git
git -C VASPBERRY-1.6.6-2e7067b fetch --depth 1 origin 2e7067b8b8b6d3ab444eef5bc025f4ce5720ca4d
GIT_LFS_SKIP_SMUDGE=1 git -C VASPBERRY-1.6.6-2e7067b checkout --detach FETCH_HEAD
git -C VASPBERRY-1.6.6-2e7067b rev-parse HEAD
```

Expected source commit: `2e7067b8b8b6d3ab444eef5bc025f4ce5720ca4d`. Resolve unexpected tag/source identity before running calculations. Existing checkouts are diagnosed, not reset or updated automatically. To reproduce the original pre-amendment v1.6.6 analysis, use its recorded `de84ed3337b900b36990f0230e5d497ff989061f` commit in a separate directory instead.

`GIT_LFS_SKIP_SMUDGE=1` leaves optional large Git LFS example data as committed pointer files during checkout, so engine help and the analytic demo do not require those downloads. Fetch needed example inputs later through the release's `examples/INPUTS.md` procedure; never pass a pointer file to a scientific calculation. See the [Git LFS setting](https://github.com/git-lfs/git-lfs/blob/main/docs/man/git-lfs-config.adoc).

If Git has globally configured LFS filters but `git-lfs` is absent from the
current PATH, `GIT_LFS_SKIP_SMUDGE=1` alone can still fail when Git starts the
missing filter process. For an engine checkout that does not yet need optional
LFS inputs, disable those filters for that command only: use
`git -c filter.lfs.process= -c filter.lfs.clean= -c filter.lfs.smudge= -c filter.lfs.required=false clone ...`
in place of `git clone ...`, or the same `-c` options before `checkout` (with
`-C <checkout-directory>` when needed). Do not change global Git configuration.
This leaves committed LFS pointers in place; fetch actual data with the
release's documented procedure before a scientific calculation.


For the Python workflow and analytic demo, select an installed Python 3.10+ interpreter and create a project environment. Check `python3 --version` first; if the system default is older, replace `python3` below with the absolute path to a supported interpreter:

```bash
python3 -m venv .venv-vaspberry-1.6.6-2e7067b
.venv-vaspberry-1.6.6-2e7067b/bin/python -m pip install -r VASPBERRY-1.6.6-2e7067b/requirements-transport.txt
.venv-vaspberry-1.6.6-2e7067b/bin/python /path/to/installed/vaspberry/scripts/vaspberry_agent.py doctor --source VASPBERRY-1.6.6-2e7067b
```

The released requirements are NumPy >=1.24 and Matplotlib >=3.7. Use the same interpreter for the helper and engine Python commands. Creating an environment requires `venv` support. The helper itself uses the Python standard library.

For native WAVECAR/topology calculations, read the checkout's `docs/BUILD.md`. Build from its root with the site's compatible Fortran compiler, Make, LP64 BLAS/LAPACK and (only when needed) MPI development/runtime stack. Examples are `make serial` for GNU serial or `make ifx-mpi` for Intel MPI. `make serial` produces `build/vaspberry`; validate with `make check-serial-help`. Compiler availability does not verify BLAS/LAPACK linkage. Platform-specific flags and libraries are documented in the build guide; Windows users need a supported Linux environment such as WSL2.

Do not run administrator package installation merely to try the demo: it needs no Fortran, MPI, VASP or WAVECAR. Public real-material examples may download large licensed/generated input files and need the release's `examples/INPUTS.md` procedure. Use the requested calculation to decide whether those inputs are relevant.
