#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
example_dir="$(cd -- "$script_dir/.." && pwd -P)"
repo_root="$(cd -- "$example_dir/../.." && pwd -P)"
wavecar="$example_dir/WAVECAR"
result_dir="${VASPBERRY_RESULT_DIR:-$example_dir/results-z2}"
mpi_nprocs="${VASPBERRY_MPI_NPROCS:-}"
expected_size=200421600
expected_sha256=a8d81854f2efc561938e478dde1be29a17ccc95d5d37325c122b9d9e82fa0838

if [[ -n "$mpi_nprocs" ]]; then
  if [[ ! "$mpi_nprocs" =~ ^[1-9][0-9]{0,2}$ ]] ||
     (( 10#$mpi_nprocs > 144 )); then
    echo "VASPBERRY_MPI_NPROCS must be an integer from 1 to 144." >&2
    exit 2
  fi
fi

if [[ ! -f "$wavecar" ]]; then
  echo "Missing examples/Bi_Z2/WAVECAR; fetch the Git LFS object first." >&2
  exit 1
fi
if grep -q '^version https://git-lfs.github.com/spec/v1' "$wavecar"; then
  echo "WAVECAR is still a Git LFS pointer; run git lfs pull first." >&2
  exit 1
fi

actual_size="$(wc -c < "$wavecar" | tr -d '[:space:]')"
if [[ "$actual_size" != "$expected_size" ]]; then
  echo "Unexpected WAVECAR size: $actual_size (expected $expected_size)." >&2
  exit 1
fi

if command -v sha256sum >/dev/null 2>&1; then
  actual_sha256="$(sha256sum "$wavecar" | awk '{print $1}')"
elif command -v shasum >/dev/null 2>&1; then
  actual_sha256="$(shasum -a 256 "$wavecar" | awk '{print $1}')"
else
  echo "A SHA-256 utility (sha256sum or shasum) is required." >&2
  exit 1
fi
if [[ "$actual_sha256" != "$expected_sha256" ]]; then
  echo "WAVECAR SHA-256 mismatch." >&2
  exit 1
fi

build_dir=""
cleanup() {
  if [[ -n "$build_dir" && -d "$build_dir" ]]; then
    rm -f -- "$build_dir/vaspberry"
    rmdir -- "$build_dir" 2>/dev/null || true
  fi
}
trap cleanup EXIT

if [[ -n "${VASPBERRY_BIN:-}" ]]; then
  requested_executable="$VASPBERRY_BIN"
  if [[ "$requested_executable" == */* ]]; then
    [[ -e "$requested_executable" ]] || {
      echo "VASPBERRY_BIN does not exist: $requested_executable" >&2
      exit 1
    }
    executable_dir="$({
      cd -- "$(dirname -- "$requested_executable")"
      pwd -P
    })"
    executable="$executable_dir/$(basename -- "$requested_executable")"
  else
    executable="$(command -v "$requested_executable" || true)"
  fi
  [[ -n "$executable" && -x "$executable" ]] || {
    echo "VASPBERRY_BIN is not executable: $requested_executable" >&2
    exit 1
  }
else
  build_dir="$(mktemp -d "${TMPDIR:-/tmp}/vaspberry-bi-z2.XXXXXX")"
  executable="$build_dir/vaspberry"
  compiler=gfortran
  source_file="$repo_root/vaspberry_gfortran_serial.f"
  compile_defines=()
  if [[ -n "$mpi_nprocs" ]]; then
    compiler=mpifort
    source_file="$repo_root/vaspberry.f"
    compile_defines=(-DMPI_USE)
  fi
  command -v "$compiler" >/dev/null 2>&1 || {
    echo "$compiler is required, or set VASPBERRY_BIN." >&2
    exit 1
  }
  "$compiler" -cpp -O2 -std=legacy -ffixed-line-length-none \
    -fallow-argument-mismatch \
    "${compile_defines[@]}" \
    -o "$executable" \
    "$source_file" \
    -llapack -lblas
fi

run_command=("$executable")
if [[ -n "$mpi_nprocs" ]]; then
  requested_launcher="${VASPBERRY_MPIEXEC:-mpiexec}"
  if [[ "$requested_launcher" == */* ]]; then
    [[ -e "$requested_launcher" ]] || {
      echo "VASPBERRY_MPIEXEC does not exist: $requested_launcher" >&2
      exit 1
    }
    launcher_dir="$({
      cd -- "$(dirname -- "$requested_launcher")"
      pwd -P
    })"
    launcher="$launcher_dir/$(basename -- "$requested_launcher")"
  else
    launcher="$(command -v "$requested_launcher" || true)"
  fi
  [[ -n "$launcher" && -x "$launcher" ]] || {
    echo "VASPBERRY_MPIEXEC is not executable: $requested_launcher" >&2
    exit 1
  }
  run_command=("$launcher" -n "$mpi_nprocs" "$executable")
fi

mkdir -p "$result_dir"
cd "$result_dir"

"${run_command[@]}" \
  -f "$wavecar" -o NFIELD -z2 1 \
  -kx 12 -ky 12 -s 2 -ii 1 -if 10 \
  >fortran.log 2>&1

test -s NFIELD.dat
test -s Z2_FIELD.csv
grep -Fq '# result_status=PASS' Z2_FIELD.csv
grep -Fq '# reportable_invariant=1' Z2_FIELD.csv
grep -Fq '# z2_invariant=1' Z2_FIELD.csv
grep -Fq '# half_top_z2_parity=1' Z2_FIELD.csv
grep -Fq '# half_bottom_z2_parity=1' Z2_FIELD.csv
grep -Fq '# half_bz_parity_consistent=1' Z2_FIELD.csv
grep -Fq '# half_top_nfield_sum=-3' Z2_FIELD.csv
grep -Fq '# half_bottom_nfield_sum=3' Z2_FIELD.csv
test ! -e INCOMPLETE
test ! -e Z2_FIELD.invalid.csv
test ! -e Z2_FIELD.tmp

if command -v python3 >/dev/null 2>&1 && \
   python3 -c 'import matplotlib, numpy' >/dev/null 2>&1; then
  python3 "$script_dir/plot_nfield.py" Z2_FIELD.csv \
    --output Z2_nfield_12x12.pdf
fi

if [[ -n "$mpi_nprocs" ]]; then
  echo "Bi n-field Z2 MPI calculation passed on $mpi_nprocs ranks: Z2=1"
else
  echo "Bi n-field Z2 serial calculation passed: Z2=1"
fi
echo "Outputs: $result_dir"
