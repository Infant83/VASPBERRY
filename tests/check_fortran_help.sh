#!/bin/sh
set -eu

if [ "$#" -ne 1 ]; then
  echo "usage: $0 HELP_OUTPUT" >&2
  exit 2
fi

help_output=$1
version=$(tr -d '[:space:]' < VERSION)

# Concurrent MPI writers can splice help text mid-line. Require one writer.
banner_count=$(grep -F -o -- 'PROGRAM INSTRUCTION' "$help_output" | wc -l | tr -d '[:space:]')
if [ "$banner_count" -ne 1 ]; then
  echo "expected exactly one help banner, found $banner_count" >&2
  exit 1
fi

for required in \
  "Ver $version" \
  "build/vaspberry --task chern --wavecar WAVECAR" \
  "--mesh NX,NY" \
  "--bands FIRST:LAST or N" \
  "--task rejects conflicting legacy task options." \
  "No integer expectation on a path." \
  "All legacy flags remain accepted" \
  "-h               : Print this help and stop" \
  "-z2  1" \
  "-kubo_bundle 1" \
  "-kubo_pairs PATH" \
  "No gap division or occupations." \
  "Source occupations may vary with k." \
  "CSV options name exact files; parent dirs must exist." \
  "Legacy DAT files may be replaced: use a fresh cwd." \
  "Reference and output names: docs/NATIVE_COMMANDS.md." \
  "External gaps must exceed 1e-5 eV." \
  "Writes bundle CSV only; default 0." \
  "Fukui-Hatsugai n-field Z2 index" \
  "full, even Gamma-centered mesh" \
  "ICHARG=11 run with ISYM=-1" \
  "Nx,Ny >= 4 and kz=0 modulo G" \
  "occupied bands 1:NE, NBANDS>NE" \
  "each of four 2D TRIM once" \
  "PASS writes NFIELD.dat" \
  "Z2_FIELD.csv" \
  "Z2_FIELD.invalid.csv" \
  "needs the final PASS CSV" \
  "result_status=PASS" \
  "top/bottom half-BZ parities must" \
  "Check denser even meshes" \
  "Compilation: see Makefile and docs/BUILD.md." \
  "GNU: make serial; make mpi" \
  "Intel: make ifx; make ifx-mpi"
do
  if ! grep -F -- "$required" "$help_output" >/dev/null; then
    echo "missing help contract: $required" >&2
    exit 1
  fi
done

for forbidden in \
  "Legacy Fukui Z2 candidate" \
  "Run tools/wavecar_z2.py" \
  "ex-MPI)mpif90 -DMPI_USE -mkl"
do
  if grep -F -- "$forbidden" "$help_output" >/dev/null; then
    echo "obsolete help contract remains: $forbidden" >&2
    exit 1
  fi
done
