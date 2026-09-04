#!/bin/sh
set -eu

if [ "$#" -ne 1 ]; then
  echo "usage: $0 HELP_OUTPUT" >&2
  exit 2
fi

help_output=$1
version=$(tr -d '[:space:]' < VERSION)

for required in \
  "Ver $version" \
  "-h               : Print this help and stop" \
  "-z2  1" \
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
