# Build and compiler portability

VASPBERRY is a Linux/Unix fixed-form Fortran program. The supported public
build paths are GNU Fortran in serial mode and GNU Fortran with Open MPI. CI
also compile- and help-tests the Intel serial source with provisioned `ifx` and
retained `ifort`; Intel MPI and numerical validation remain manual.

## Support matrix

| Build | Source | Command | Validation status |
|---|---|---|---|
| GNU serial | `vaspberry_gfortran_serial.f` | `make serial` | CI builds and help-tests the GNU versions supplied by Ubuntu 22.04 and 24.04 |
| GNU + Open MPI | `vaspberry.f` | `make mpi` | CI builds it and exercises two-rank help, `MPI_DOUBLE_PRECISION` reduction, and result-status broadcast |
| Intel `ifx` serial | `vaspberry_gfortran_serial.f` | `make ifx` | CI compile/help target with `ifx` 2025.0 and system LP64 BLAS/LAPACK; numerical test still required |
| Intel `ifx` + Intel MPI | `vaspberry.f` | `make ifx-mpi` | Build recipe reviewed; manual test required on a oneAPI host |
| Intel Classic `ifort` | serial source or `vaspberry.f` with Intel MPI | `make ifort` or `make ifort-mpi` | CI compile/help target for serial `ifort` 2021.10; Intel MPI and numerical tests remain manual |

Intel discontinued `ifort` in the oneAPI 2025 release and recommends `ifx` for
continued support. Retaining an `ifort` recipe helps older clusters, but new
installations should use `ifx`. See Intel's
[ifort-to-ifx porting guide](https://www.intel.com/content/www/us/en/developer/articles/guide/porting-guide-for-ifort-to-ifx.html)
and [Intel MPI compiler-wrapper list](https://www.intel.com/content/www/us/en/docs/mpi-library/developer-reference-linux/2021-16/compiler-commands.html).

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

The products are `build/vaspberry-gfortran` and
`build/vaspberry-mpi`. A normal MPI calculation is launched, for example, as

```bash
mpiexec -n 4 build/vaspberry-mpi -f WAVECAR -kx 12 -ky 12 \
  -s 2 -ii 1 -if 10 -z2 1
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

1. both production sources compile and link against LP64 BLAS/LAPACK;
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

## Intel oneAPI build

First load the oneAPI compiler, oneMKL, and Intel MPI environment provided by
the local installation. A common system-wide setup is:

```bash
source /opt/intel/oneapi/setvars.sh
```

Then build and check each executable available at the site:

```bash
make ifx
build/vaspberry-ifx -h

make ifx-mpi
mpiexec -n 2 build/vaspberry-ifx-mpi -h
```

For a retained Classic installation, substitute:

```bash
make ifort
make ifort-mpi
```

The Intel targets use preprocessing, extended fixed-form source lines,
byte-based direct-access records, and sequential LP64 oneMKL. Intel documents
`-qmkl=sequential` as the sequential oneMKL link option. The MPI targets use
`mpiifx` or the retained `mpiifort` wrapper and define `MPI_USE`.

Public CI provisions `ifx` 2025.0 and retained `ifort` 2021.10 through
`fortran-lang/setup-fortran`, links the serial executable against Ubuntu's
system LP64 BLAS/LAPACK, and checks `-h`. This verifies the fixed-form source
and serial link path without assuming oneMKL is installed on the runner. Intel
MPI wrappers are not provisioned there, so `ifx-mpi`/`ifort-mpi`, the oneMKL
recipes, and numerical comparisons remain explicitly manual checks.

## WAVECAR record length and numerical-library ABI

VASPBERRY reads `WAVECAR` with unformatted direct access. The writer and reader
must therefore agree on the meaning of `RECL`. GNU Fortran uses byte file
storage units for this build. Intel builds must keep `-assume byterecl`; without
it, the same integer record length may be interpreted in four-byte units.
The source explicitly opens both header and data views with
`FORM='UNFORMATTED'` and `ACTION='READ'`.

A record-length mismatch can leave the first header partly plausible while
making `NKPOINT`, `NBANDS`, `ENCUT`, or lattice data invalid. The sound fix is
to rebuild the reader with the writer's byte convention and regenerate the
`WAVECAR`, not to tune scientific thresholds around corrupt input.

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
