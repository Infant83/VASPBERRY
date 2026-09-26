# Build and compiler portability

VASPBERRY's main program is fixed-form Fortran. Compile it once for the machine
where it will run, using a Fortran compiler and **LP64 BLAS/LAPACK**. MPI is
needed only for the parallel executable. Python, a Python environment manager,
VASP, and Wannier90 are not build dependencies. A calculation reads an existing
VASP `WAVECAR`; Python tools are optional for subsequent analysis and plotting.

The main route below uses **Intel oneAPI Fortran (`ifx`), Intel MPI, and
sequential oneMKL on Linux**. Choose the GNU route if those packages are not
available. All commands run from the source directory containing `Makefile`.
Extract or clone into a writable path without spaces; the Makefile rejects
whitespace in the absolute repository/build path.
Use the target for your chosen compiler; bare `make` builds both GNU serial
and MPI executables, and therefore also needs a GNU-compatible MPI SDK.

## Required build and execution environment

| Component | Build requirement | Requirement when running |
|---|---|---|
| Operating system and utilities | Linux or macOS; GNU Make, a POSIX shell and standard shell utilities (including `cmp` and `cksum`); C runtime and development/linker tools | Compatible operating system and CPU architecture |
| Fortran compiler | GNU `gfortran`, Intel `ifx`, or an existing Intel Classic `ifort` installation | The linked compiler runtime libraries |
| Numerical libraries | LP64 BLAS and LAPACK development libraries; sequential oneMKL is the Intel default | The corresponding shared libraries, unless linked statically |
| MPI, for parallel builds | An MPI SDK with Fortran support and its compiler wrapper (`mpiifx`, `mpiifort`, or `mpifort`) | The matching MPI runtime and launcher, available on every allocated node |
| Input and storage | No material input is needed to compile or display help | A compatible `WAVECAR`, a writable result directory, and memory/disk appropriate to its size |

The source uses POSIX C functions through `iso_c_binding`; a native Windows
compiler alone is not a supported build route. Use a Linux environment under
WSL2 as described below. Do not mix x86-64 and arm64 compiler/library packages.

No root privileges are needed to build VASPBERRY in a writable source directory.
Administrative commands below install system dependencies only. At a cluster,
use its installed compiler/library modules instead.

## Get the source

To follow the main development line, clone the repository normally; its default
branch is named `master`:

```bash
git clone https://github.com/Infant83/VASPBERRY.git
cd VASPBERRY
make help
```

To install the fixed release instead, download **Source code (tar.gz)** or
**Source code (zip)** from the [v1.4.2 release page](https://github.com/Infant83/VASPBERRY/releases/tag/v1.4.2),
extract it, and enter `VASPBERRY-1.4.2`. The same Makefile commands work without
Git metadata. For example:

```bash
curl -fL https://github.com/Infant83/VASPBERRY/archive/refs/tags/v1.4.2.tar.gz \
  -o VASPBERRY-1.4.2.tar.gz
tar -xzf VASPBERRY-1.4.2.tar.gz
cd VASPBERRY-1.4.2
make help
```

Git, Git LFS, and Python are not required to compile an extracted release.
For a Git checkout, the large Bi example may initially contain a Git LFS
pointer; obtain its actual input as described in the
[Bi example](../examples/Bi_Z2/README.md) before running that example. The
small MoS₂ input used below is included directly in the source tree.

<a id="support-matrix"></a>

## Platform and validation matrix

| Environment | Build route | Validation coverage |
|---|---|---|
| Ubuntu 22.04, Intel `ifx` 2025.0 and Intel MPI | `make ifx-mpi` | Serial and two-rank MPI runtime/help plus actual MoS₂ curvature and pair-export comparisons, using oneMKL |
| Ubuntu 22.04, retained Intel Classic `ifort` 2021.10 and Intel MPI | `make ifort-mpi` | Separate serial and MPI numerical comparison using oneMKL |
| Ubuntu 22.04/24.04, GNU 11/13 and Open MPI | `make gnu` | Compile/link, serial and MPI runtime/help; additional GNU numerical regression and Bi Z₂ validation |
| Ubuntu 22.04, GNU and MPICH | `make gnu MPIFC=mpifort.mpich` | Clean source-archive installation job checks serial and MPI native MoS₂ results |
| macOS 15 ARM64, Homebrew GNU/Open MPI/OpenBLAS | Explicit Homebrew paths below | Clean source-archive installation job checks the GNU build and serial/MPI native MoS₂ results |
| macOS Intel, existing GNU 15.1/Open MPI 5.0.7 and system LP64 BLAS/LAPACK | Existing compiler/MPI paths with `-llapack -lblas` | Local clean source-archive check using Apple system libraries; no claim of a newly provisioned Homebrew Intel environment |
| Debian, compatible GNU and Open MPI | Ubuntu/Debian package recipe below | Package recipe; Debian itself is not a separate CI runner |
| Windows with Ubuntu under WSL2 | Linux GNU recipe inside WSL2 | Documented route; WSL2 itself is not a tested release environment |
| Other Linux clusters and compiler/library versions | Site modules and Make overrides below | Run the local acceptance checks; compatibility is not implied by the existence of a recipe |

The [1.4.2 validation record](VALIDATION_1.4.2.md) identifies the release commit,
exact runner/compiler/library versions, and successful numerical checks.
A compiler/version omitted from that record is not claimed as tested.

<a id="intel-oneapi-build"></a>

## Intel oneAPI on Linux

Install or load **all three development components**: Intel Fortran (`ifx`),
the Intel MPI SDK, and oneMKL development libraries. An MPI runtime-only
installation cannot compile the executable because it lacks the wrapper and
headers. Intel supplies the compiler as a standalone package or through its
HPC Toolkit; see the [Linux compiler setup](https://www.intel.com/content/www/us/en/docs/fortran-compiler/get-started-guide/2025-2/get-started-on-linux.html)
and [MPI compiler-wrapper documentation](https://www.intel.com/content/www/us/en/docs/mpi-library/developer-reference-linux/2021-16/compiler-commands.html).

Load the environment in the same Bash terminal where you build and run. A
common component-layout installation uses:

```bash
source /opt/intel/oneapi/setvars.sh
command -v ifx mpiifx mpiexec.hydra
ifx --version
mpiifx -show
```

For a unified-layout installation, source the installed version's
`oneapi-vars.sh` instead. A personal installation may be under a different
prefix; cluster modules can provide the same settings. Follow Intel's
[environment setup layouts](https://www.intel.com/content/www/us/en/docs/oneapi/programming-guide/2025-1/use-the-setvars-and-oneapi-vars-scripts-with-linux.html).
Load this environment again in each new terminal and in your batch job script.

Build the parallel executable and check two-rank communication and native help:

```bash
make ifx-mpi
make check-ifx-mpi
mpiexec.hydra -n 4 build/vaspberry-ifx-mpi --help
```

This creates `build/vaspberry-ifx-mpi`. Four ranks means four MPI processes;
choose a count allowed by the scheduler allocation. The explicit Intel MPI
launcher avoids accidentally using an Open MPI `mpiexec` from another package.
A minimal real calculation is:

```bash
repo_dir="$PWD"
mkdir results-intel-kubo-01
(
  cd results-intel-kubo-01
  mpiexec.hydra -n 4 "$repo_dir/build/vaspberry-ifx-mpi" --task kubo \
    --wavecar "$repo_dir/examples/1H-MoS2/KPATH/2.band/WAVECAR" \
    --spinor 2 --bands 1:18 --bundle 1 --curvature-csv KUBO.csv \
    > vaspberry.log 2> vaspberry.err
)
```

It writes the supplied 48-point path's bundle curvature to `KUBO.csv`.
The [hands-on guide](HANDS_ON.md) explains the columns and plotting, followed
by pair export and occupation integration for transport. No Python command
is needed for this native calculation.

For serial operation, use `make ifx` and `make check-ifx`, then run
`build/vaspberry-ifx` without an MPI launcher. For a retained Classic setup:

```bash
make ifort-mpi
make check-ifort-mpi
mpiexec.hydra -n 4 build/vaspberry-ifort-mpi --help
```

Its serial targets are `make ifort` and `make check-ifort`, producing
`build/vaspberry-ifort`. Intel removed `ifort` from oneAPI 2025; the Classic
route is for existing installations, while new installations should use
`ifx`. See Intel's [porting guide](https://www.intel.com/content/www/us/en/developer/articles/guide/porting-guide-for-ifort-to-ifx.html).

These targets enable preprocessing, extended fixed-form source lines,
byte-based direct-access records, and sequential LP64 oneMKL. MPI builds
add `MPI_USE` and use `mpiifx` or `mpiifort`. `check-ifx-mpi` and
`check-ifort-mpi` also compile a small MPI runtime probe. Their launcher
variables are `INTEL_MPIEXEC` (default `mpiexec.hydra`) and
`INTEL_MPIEXEC_FLAGS`; use site-approved values when necessary.

<a id="gnu-build"></a>

## GNU on Ubuntu or Debian

Install the serial build dependencies first:

```bash
sudo apt-get update
sudo apt-get install build-essential gfortran libblas-dev liblapack-dev
make serial
make check-serial-help
build/vaspberry --help
```

`build-essential` supplies GNU Make and the system development/linking tools.
Use the ordinary LP64 BLAS/LAPACK packages, not packages with `64` integer
interfaces. The package distinction is also visible in the
[Debian LAPACK package information](https://packages.debian.org/trixie/liblapack-dev).

For Open MPI, install its runtime and development package, then build/check:

```bash
sudo apt-get install openmpi-bin libopenmpi-dev
make gnu
make check-gnu
mpiexec -n 4 build/vaspberry-mpi --help
```

`make gnu` creates both `build/vaspberry` (serial) and
`build/vaspberry-mpi`. `build/vaspberry-gfortran` remains a compatibility
symlink to the serial executable. If the machine provides fewer than two MPI
slots, a local Open MPI smoke check can use
`make check-gnu MPIEXEC_FLAGS=--oversubscribe`. Do not add that option to
scheduler jobs unless the site explicitly permits oversubscription.

The MPICH release installation check uses **Ubuntu 22.04**. For that route,
use the package-specific wrapper and launcher so the system's MPI alternatives
cannot select the wrong implementation:

```bash
sudo apt-get install mpich libmpich-dev
make gnu BUILD_DIR=build-mpich MPIFC=mpifort.mpich
make check-gnu BUILD_DIR=build-mpich MPIFC=mpifort.mpich \
  MPIEXEC=mpiexec.mpich
mpiexec.mpich -n 4 build-mpich/vaspberry-mpi --help
```

`make check-gnu` checks actual MPI communication, not just whether the
executable starts. Its runtime probe requires at least two ranks in the same
`MPI_COMM_WORLD` and then checks reduction and broadcast. If two processes
both report rank 0 in a one-process world, the check fails; a successful
`--help` invocation alone would not detect that broken parallel launch.

The release-preparation run on Ubuntu 24.04 with its packaged
**MPICH 4.2.0-5build3** encountered that failure with `mpiexec.mpich`: its
PMIx configuration was incompatible with the packaged Hydra launch path.
The problem is documented in the
[upstream MPICH issue](https://github.com/pmodels/mpich/issues/7064#issuecomment-2301026290).
For Ubuntu 24.04, use the Open MPI route above, or a corrected site-provided
MPICH compiler/runtime/launcher combination that passes `make check-gnu`.
Do not mix an MPICH executable with another MPI implementation's launcher as
an installation workaround.

Other distributions may name those commands differently. Use `mpifort` and
`mpiexec` from the same loaded MPI installation rather than copying an MPI
library list into the compiler command.

## GNU on macOS with Homebrew

The current hosted installation check uses macOS 15 on Apple Silicon.
Homebrew classifies Intel Mac installations as Tier 3; recent formulae may
require source builds there. Our Intel Mac check uses existing GNU 15.1/Open MPI 5.0.7 and Apple system
LP64 BLAS/LAPACK (Accelerate), rather than the OpenBLAS recipe below. It does
not establish support for freshly installing the latest Homebrew packages on
Intel hardware. See the linked Homebrew requirements.

Install the Xcode Command Line Tools if they are absent
(`xcode-select --install`) and set up [Homebrew](https://docs.brew.sh/Installation). Then
install its [GNU compiler](https://formulae.brew.sh/formula/gcc),
[Open MPI](https://formulae.brew.sh/formula/open-mpi), and
[OpenBLAS](https://formulae.brew.sh/formula/openblas):

```bash
brew install gcc open-mpi openblas
vaspberry_openblas="$(brew --prefix openblas)"
export FC="$(brew --prefix gcc)/bin/gfortran"
export MPIFC="$(brew --prefix open-mpi)/bin/mpifort"
export MPIEXEC="$(brew --prefix open-mpi)/bin/mpiexec"
export GNU_LIBS="-L${vaspberry_openblas}/lib -Wl,-rpath,${vaspberry_openblas}/lib -lopenblas"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
make gnu BUILD_DIR=build-macos
make check-gnu BUILD_DIR=build-macos
"$MPIEXEC" -n 4 build-macos/vaspberry-mpi --help
```

The formula provides the Fortran compiler as `gfortran`. Its Open MPI wrapper
must use the same GNU installation; inspect it with
`"$MPIFC" --showme:command`. The explicit OpenBLAS prefix is necessary because
Homebrew's OpenBLAS is keg-only. `-Wl,-rpath,...` retains that shared-library
path for execution. The recipe uses LP64 OpenBLAS for both BLAS and LAPACK;
it does not depend on whichever `-lblas` happens to resolve on the machine.
The thread settings prevent each MPI rank from spawning extra numerical
library threads during these checks.

`brew --prefix` handles the installation prefix without assuming `/usr/local`
or `/opt/homebrew`. Build and run in one architecture environment; avoid
mixing an arm64 compiler with x86-64 libraries from a Rosetta installation.
For serial operation, `make serial BUILD_DIR=build-macos` uses the same
exported compiler and library settings, and does not need MPI installed.
The final parallel `--help` command needs four available MPI slots; use
`-n 2` on a smaller machine.

## Windows through WSL2

The supplied Makefile targets Linux/Unix, so use Ubuntu inside WSL2 rather
than a Windows Command Prompt for compilation. Follow Microsoft's
[WSL installation instructions](https://learn.microsoft.com/en-us/windows/wsl/basic-commands)
(`wsl --install` from PowerShell), open the installed Ubuntu terminal, and
follow the GNU Ubuntu recipe above. Obtain and build the source inside that
Linux environment. Microsoft's [environment guide](https://learn.microsoft.com/en-us/windows/wsl/setup/environment)
recommends keeping Linux-tool workloads in the WSL filesystem.

This is a documented route, not a separately validated WSL2 platform. Native
Windows, MinGW, and MSYS2 builds are not claimed as supported by this release.

## Cluster modules and build overrides

Use the cluster's documented compiler, BLAS/LAPACK, and matching MPI module
set. Module names and batch launch commands differ by site; none is assumed
here. After loading the site's environment, check `command -v` for the chosen
compiler, wrapper, and launcher. For Open MPI, `mpifort --showme:command`
identifies the underlying compiler; for Intel MPI, use `mpiifx -show` or
`mpiifort -show`. Open MPI explains its wrappers and compiler selection in
its [wrapper documentation](https://docs.open-mpi.org/en/v5.0.x/man-openmpi/man1/ompi-wrapper-compiler.1.html).

The Makefile honors exported compiler/wrapper/launcher variables and explicit
command-line overrides. For example, with an Open MPI build compatible with
GNU 13:

```bash
OMPI_FC=gfortran-13 make gnu BUILD_DIR=build-gfortran13 \
  FC=gfortran-13 MPIFC=mpifort
OMPI_FC=gfortran-13 make check-gnu BUILD_DIR=build-gfortran13 \
  FC=gfortran-13 MPIFC=mpifort MPIEXEC=mpiexec
```

Changing the wrapper's underlying compiler does not make an incompatible MPI
installation compatible. Use an MPI SDK built for that compiler. Load the
same modules in the compute-node job environment and use the launch command
required by the site's scheduler; run multi-rank checks within an allocation.

| Make variable | Purpose |
|---|---|
| `FC`, `MPIFC` | GNU serial compiler and MPI Fortran wrapper |
| `IFX`, `MPIIFX`, `IFORT`, `MPIIFORT` | Intel compiler and MPI wrapper commands |
| `GNU_LIBS` | GNU LP64 link libraries and their search paths; default `-llapack -lblas` |
| `IFX_MKL_FLAGS`, `IFORT_MKL_FLAGS` | Intel library-link options; defaults select sequential oneMKL |
| `FFLAGS`, `LDFLAGS` | Additional compile/link flags, without replacing required source-format flags |
| `GNU_FLAGS`, `INTEL_FLAGS` | Full compiler flag presets; preserve preprocessing, fixed-form, and Intel byte-RECL settings if replacing them |
| `MPIEXEC`, `MPIEXEC_FLAGS` | GNU MPI check launcher and its additional arguments |
| `INTEL_MPIEXEC`, `INTEL_MPIEXEC_FLAGS` | Intel MPI check launcher and its additional arguments |
| `BUILD_DIR` | `build` or a direct-child directory named `build-*` |

For a local LP64 library outside the system search path, set `GNU_LIBS` to
its link flags or add its search/runtime path with `LDFLAGS`; there is no need
to edit the Makefile. Keep the executable in its build directory and use its
absolute path from result directories. There is no system-wide installation
step and no requirement to copy the executable into `/usr/local/bin`.

Use a distinct directory such as `build-ifx`, `build-openmpi`, or
`build-mpich` to retain executables for several toolchains. GNU build
configuration stamps trigger recompilation when compiler commands, flags,
source selection, MPI wrapper controls, or the Makefile change. If a compiler
or library is replaced
in place under the same command name, use a new directory or clean that build:

```bash
make clean BUILD_DIR=build-mpich
```

`BUILD_DIR` is restricted to a direct-child `build` or `build-*` directory;
source paths, the repository root, nested/outside paths, and symlinks are
rejected before creating or cleaning outputs.

Starting with 1.3.0, serial and MPI defaults share the current `vaspberry.f`,
including Kubo modes. The historical reduced serial source, which has no Kubo
mode, is available only by an explicit override:

```bash
make serial SERIAL_SOURCE=vaspberry_gfortran_serial.f BUILD_DIR=build-legacy
```

GNU builds retain `-cpp -ffixed-line-length-none` and
`-fallow-argument-mismatch` for the legacy interfaces; Intel uses
`-fpp -extend-source -assume byterecl`. Do not disable preprocessing or
remove the legacy-interface flag merely to silence compiler diagnostics.

## Installation troubleshooting

| Symptom | Check and correction |
|---|---|
| `make`, `gfortran`, `ifx`, or an MPI wrapper is not found | Install the development package, or load its environment/modules in this terminal. Intel MPI requires its SDK, not only its runtime. `make help` lists compiler-specific targets. |
| `mpif.h` is missing, or MPI routines fail to link | Build with the matching Fortran MPI wrapper and MPI development package; do not use the plain compiler for an MPI target. |
| `zgesvd_`/`zgetrf_` are unresolved, or `-llapack`/`-lblas` cannot be found | Install LP64 LAPACK/BLAS development libraries. For Homebrew use the explicit OpenBLAS recipe; for Intel load oneMKL. Keep libraries after the source/object in manual link commands. |
| The executable cannot load `libgfortran`, `libmpi`, or a numerical library | Restore the build environment on the running node. Inspect Linux dependencies with `ldd build/vaspberry-mpi`, or macOS dependencies with `otool -L build-macos/vaspberry-mpi`. Correct runtime paths rather than copying arbitrary library files. |
| MPI starts incorrectly, reports rank/PMI errors, or hangs | Check the linked MPI, wrapper, and launcher come from the same installation; use the scheduler's allocation and launcher. Do not launch an Intel MPI binary with Open MPI or vice versa. |
| MPI reports too few slots | Request enough scheduler ranks or reduce `-n`. Local Open MPI smoke tests can use the documented `--oversubscribe` override. |
| Compiler changes appear to leave an old executable | Use a fresh `BUILD_DIR`; command/flag changes are tracked, but replacement software under the same path may require rebuilding explicitly. |
| WAVECAR headers or record reads are invalid | Check the native byte-RECL input contract below and that the file is a complete WAVECAR, not a Git LFS pointer. This is an input-format issue, not a reason to change numerical thresholds. |

## What the MPI check covers

The Z2 work is divided into disjoint plaquette ranges. Every rank leaves
unowned field entries at zero; the fields and half-zone sums are reduced with
`MPI_DOUBLE_PRECISION`/`MPI_SUM` to rank 0. Rank 0 applies the common output
guards and broadcasts the integer pass/fail status before all ranks continue or
abort.

All remaining `REAL*8` collective buffers in the legacy Chern, spectrum,
optical-selectivity, and Kubo paths also use `MPI_DOUBLE_PRECISION`; the
optional implementation-specific `MPI_REAL8` name is not required. Output
array arguments use real arrays even in branches that do not inspect them, so
the build no longer relies on scalar-to-array placeholder aliasing.

The MPI program initializes MPI before parsing its command line. Consequently,
the `-h` path checks the initialization/finalization state and calls
`MPI_FINALIZE` before `STOP`. The normal completion path also finalizes MPI;
fatal paths use `MPI_ABORT` so that a rank-local error does not strand the other
ranks.

CI checks three distinct layers:

1. current serial and MPI builds compile and link against LP64 BLAS/LAPACK;
2. the MPI executable starts on two ranks and exposes the same command-line
   help as the serial executable;
3. a two-rank runtime test exercises the real(8) reduction and integer
   broadcast used by the Z2 path, while source-linked helper tests exercise the
   production reciprocal-G mapping, time reversal, overlap conditioning, and
   output guards.

These checks do not make serial and MPI floating-point output bitwise
identical: reduction order can change the final rounding. A scientific
comparison should instead require the same PASS/INVALID status, the same
integer n-field and Z2 parity, and agreement of continuous diagnostics within
their documented tolerances. A full calculation must still be repeated at a
denser k mesh to establish convergence.

The separate Bi validation workflow performs that comparison on the tracked
12 x 12 `WAVECAR`. It runs `build/vaspberry-gfortran` and a two-rank
`build/vaspberry-mpi` in different result directories, requires each result to
independently pass the v1.2 Z2 guards, compares all discrete field values
exactly, and compares finite continuous fields with relative tolerance
`1e-11` and absolute tolerance `1e-12`. Logs are not compared because MPI
rank ordering and timing are not numerical results.

The example runner selects MPI launch mode only when
`VASPBERRY_MPI_NPROCS` is a validated integer from 1 through 144. For example:

```bash
make mpi
VASPBERRY_BIN="$PWD/build/vaspberry-mpi" \
VASPBERRY_MPI_NPROCS=2 \
VASPBERRY_RESULT_DIR="$PWD/examples/Bi_Z2/results-z2-mpi" \
  ./examples/Bi_Z2/scripts/run_z2.sh
```

If `VASPBERRY_BIN` is omitted in this mode, the runner builds `vaspberry.f`
with `mpifort -DMPI_USE`; `VASPBERRY_MPIEXEC` may select a site-specific MPI
launcher. Serial and MPI runs must not share a result directory because the
VASPBERRY output names are fixed during each run.

## WAVECAR record length and numerical-library ABI

VASPBERRY reads `WAVECAR` with unformatted direct access. The writer and reader
must therefore agree on the meaning of `RECL`. GNU Fortran uses byte file
storage units for this build. The supplied Intel targets also use
`-assume byterecl` for this byte-RECL input contract; without it, the same
integer record length may be interpreted in four-byte units.
The source explicitly opens both header and data views with
`FORM='UNFORMATTED'` and `ACTION='READ'`.

A record-length mismatch can leave the first header partly plausible while
making `NKPOINT`, `NBANDS`, `ENCUT`, or lattice data invalid. The sound fix is
to match the reader to the writer's record convention or regenerate the
`WAVECAR` with the byte convention used by the supplied builds. Historical
four-byte-word-RECL files need a matching validated Intel build or a regenerated
byte-RECL input; the native GNU executable does not autodetect that convention.
The Python WAVECAR reader identifies both layouts, which does not change the
native executable's file contract. Do not tune scientific thresholds around
corrupt input.

The source calls `ZGESVD` and `ZGETRF` with default 32-bit Fortran integers.
Use the LP64 BLAS/LAPACK interface. Do not add GNU
`-fdefault-integer-8`, Intel `-i8`, or `-qmkl-ilp64` unless every BLAS, LAPACK,
MPI, and file-format interface is audited together.

The guarded output path also calls the POSIX C library functions `realpath`,
`free`, `strlen`, and `rename` through `iso_c_binding`. The current native
build and test contract is therefore Linux/Unix;
Windows users should follow the WSL2 route described above, whose platform
validation remains separate from the Linux CI checks.

## Manual compiler acceptance checklist

For a compiler not exercised in CI, record all of the following before calling
it supported:

- compiler and MPI-wrapper versions;
- the complete compile and link commands;
- successful serial and two-rank `-h` runs;
- successful source-linked Z2 helper tests;
- the same Z2 status and parity from serial and MPI runs on one reviewed
  `WAVECAR`;
- continuous diagnostic differences within the recorded tolerances; and
- a successful denser-mesh convergence calculation.
