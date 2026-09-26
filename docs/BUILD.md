# Build and compiler portability

VASPBERRY is a Linux/Unix fixed-form Fortran program. The main walkthrough
uses **Intel oneAPI Fortran (`ifx`), Intel MPI and sequential LP64 oneMKL**.
Existing Classic installations can use `ifort`; GNU serial and Open MPI
recipes remain available. The table separates available build recipes from
the compiler and numerical checks actually performed.

## Support matrix

| Build | Source | Command | Validation status |
|---|---|---|---|
| Intel `ifx` serial | `vaspberry.f` without `MPI_USE` | `make ifx` | Serial portability CI; matching serial reference in the Intel MPI numerical job |
| Intel `ifx` + Intel MPI | `vaspberry.f` | `make ifx-mpi` | `make check-ifx-mpi` and the required actual-input Intel MPI CI job |
| Intel Classic `ifort` | `vaspberry.f`, optionally with `MPI_USE` | `make ifort` or `make ifort-mpi` | Serial portability CI; `make check-ifort-mpi` and a separate required Intel MPI numerical job |
| GNU serial | `vaspberry.f` without `MPI_USE` | `make serial` | Local numerical regression and GNU 11/13 CI pass on Ubuntu 22.04/24.04 |
| GNU + Open MPI | `vaspberry.f` | `make mpi` | CI builds it and checks help, MPI collectives, native serial/MPI calculations and the separate Bi Z₂ case |

Starting with 1.3.0, default serial and MPI builds share the current source,
including the Kubo commands. To reproduce the historical reduced serial
implementation explicitly (it has no Kubo mode), use a separate build directory:

```bash
make serial SERIAL_SOURCE=vaspberry_gfortran_serial.f BUILD_DIR=build-legacy
```

The [1.4.1 validation record](VALIDATION_1.4.1.md) links publication checks
and separates compiler smoke tests from actual numerical comparisons.

Intel discontinued `ifort` in the oneAPI 2025 release and recommends `ifx` for
continued support. Retaining an `ifort` recipe helps older clusters, but new
installations should use `ifx`. See Intel's
[ifort-to-ifx porting guide](https://www.intel.com/content/www/us/en/developer/articles/guide/porting-guide-for-ifort-to-ifx.html)
and [Intel MPI compiler-wrapper list](https://www.intel.com/content/www/us/en/docs/mpi-library/developer-reference-linux/2021-16/compiler-commands.html).

## Intel oneAPI build

First load the oneAPI compiler, oneMKL, and Intel MPI environment provided by
the local installation. In Bash, a common system-wide setup is:

```bash
source /opt/intel/oneapi/setvars.sh
```

For a unified-layout installation, source its versioned `oneapi-vars.sh`
instead; cluster modules can supply the same environment. Use the actual
installation path. Intel documents both
[environment setup layouts](https://www.intel.com/content/www/us/en/docs/oneapi/programming-guide/2024-2/use-the-setvars-and-oneapi-vars-scripts-with-linux.html).

Build the MPI executable and inspect its arguments:

```bash
make ifx-mpi
make check-ifx-mpi
mpiexec -n 4 build/vaspberry-ifx-mpi --help
```

The executable is `build/vaspberry-ifx-mpi`. Four ranks means four MPI
processes; choose a count permitted by the scheduler allocation. Use the
Intel MPI launcher from the same environment as `mpiifx`. A minimal real
calculation from the repository root is:

```bash
repo_dir="$PWD"
mkdir results-intel-kubo-01
(
  cd results-intel-kubo-01
  mpiexec -n 4 "$repo_dir/build/vaspberry-ifx-mpi" --task kubo \
    --wavecar "$repo_dir/examples/1H-MoS2/KPATH/2.band/WAVECAR" \
    --spinor 2 --bands 1:18 --bundle 1 --curvature-csv KUBO.csv \
    > vaspberry.log 2> vaspberry.err
)
```

This writes the supplied 48-point path's bundle curvature to `KUBO.csv`.
The [hands-on guide](HANDS_ON.md) explains the columns and direct plotting,
then the separate pair-export and occupation-integration workflow.
Python is not required to compile or run this native calculation.

For a retained Classic installation, use its environment and substitute:

```bash
make ifort-mpi
make check-ifort-mpi
mpiexec -n 4 build/vaspberry-ifort-mpi --help
```

For serial operation, `make ifx` creates `build/vaspberry-ifx`; run it without
`mpiexec`. The corresponding Classic serial target is `make ifort`.

The Intel targets use preprocessing, extended fixed-form source lines,
byte-based direct-access records, and sequential LP64 oneMKL. Intel documents
`-qmkl=sequential` as the sequential oneMKL link option. The MPI targets use
`mpiifx` or the retained `mpiifort` wrapper and define `MPI_USE`.

`check-ifx-mpi` and `check-ifort-mpi` build the corresponding executable,
compile an MPI runtime probe, check its two-rank reduction/broadcast and
validate native help on two ranks. These check targets use
`INTEL_MPIEXEC=mpiexec.hydra`; override that Make variable and
`INTEL_MPIEXEC_FLAGS` only when the site requires another Intel MPI launcher.

The existing serial portability CI uses `ifx` 2025.0 and retained `ifort`
2021.10 against system LP64 BLAS/LAPACK. The separate
[Intel MPI validation workflow](../.github/workflows/intel-mpi-validation.yml)
provisions Intel MPI and oneMKL, exercises the default build/check targets,
then compares actual MoS₂ serial and two-rank MPI bundle curvature and pair
exports. Its ifx and ifort jobs are required publication checks for 1.4.1;
the release validation record identifies the successful source commit.

## GNU build

Install a Fortran compiler, LP64 BLAS/LAPACK, and, for the MPI executable, an
MPI development package. On Ubuntu:

```bash
sudo apt-get install gfortran libblas-dev liblapack-dev \
  openmpi-bin libopenmpi-dev
```

Build both executables and run the portable smoke checks:

```bash
make gnu
make check-gnu
```

The products are `build/vaspberry` and `build/vaspberry-mpi`.
`build/vaspberry-gfortran` is retained as a compatibility symlink to the serial
executable; existing scripts continue to work. A normal MPI calculation is launched, for example, as

```bash
mpiexec -n 4 build/vaspberry-mpi --task z2 --wavecar WAVECAR \
  --mesh 12,12 --spinor 2 --bands 1:10 --output NFIELD
```

Compiler commands can be overridden without editing the Makefile. This is
useful on module-based clusters or when checking another GNU release:

```bash
make clean
OMPI_FC=gfortran-13 make gnu FC=gfortran-13 MPIFC=mpifort
OMPI_FC=gfortran-13 make check-gnu FC=gfortran-13 MPIFC=mpifort \
  MPIEXEC=mpiexec MPIEXEC_FLAGS=--oversubscribe
```

`OMPI_FC` selects the underlying compiler for an Open MPI wrapper. Other MPI
implementations use different wrapper controls; use a wrapper built for, or
explicitly configured with, the selected Fortran compiler.

`BUILD_DIR` may be overridden only with the direct child `build` or a
`build-*` directory such as `build-gfortran13`. The Makefile rejects source
directories, the repository root, nested paths, and outside paths before
creating or recursively cleaning the directory. Restricting the target to a
build-named direct child also prevents `clean` from following an intermediate
symlink outside the repository.

`-fallow-argument-mismatch` is currently required by the legacy `mpif.h`
interface and heterogeneous external-procedure calls. Removing it requires an
interface modernization rather than merely dropping the flag.

Both production files use fixed form. GNU builds keep preprocessing enabled
and remove the fixed-form line limit with `-cpp -ffixed-line-length-none`;
Intel builds use the corresponding `-fpp -extend-source`. Do not compile the
MPI source without preprocessing, because the `MPI_USE` branches must be
selected consistently.

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
Fortran output names are fixed during each run.

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
`free`, `strlen`, and `rename` through `iso_c_binding`. The current build and
test contract is therefore Linux/Unix; Windows portability has not been
established.

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
