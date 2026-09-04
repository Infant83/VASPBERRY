#!/usr/bin/env bash
set -euo pipefail

example_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd -- "$example_dir/.." && pwd -P)"
wavecar="$example_dir/WAVECAR"
result_dir="$example_dir/results-v1.1.1"
expected_size=200421600
expected_sha256=a8d81854f2efc561938e478dde1be29a17ccc95d5d37325c122b9d9e82fa0838

if [[ ! -f "$wavecar" ]]; then
  echo "Missing Bi_Z2/WAVECAR; fetch the Git LFS object first." >&2
  exit 1
fi
if grep -q '^version https://git-lfs.github.com/spec/v1' "$wavecar"; then
  echo "Bi_Z2/WAVECAR is still a Git LFS pointer; run git lfs pull first." >&2
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

command -v gfortran >/dev/null 2>&1 || {
  echo "gfortran is required." >&2
  exit 1
}
command -v python3 >/dev/null 2>&1 || {
  echo "python3 is required." >&2
  exit 1
}

build_dir="$(mktemp -d "${TMPDIR:-/tmp}/vaspberry-bi-z2.XXXXXX")"
cleanup() {
  rm -f -- "$build_dir/vaspberry"
  rmdir -- "$build_dir" 2>/dev/null || true
}
trap cleanup EXIT

gfortran -cpp -O2 -ffixed-line-length-none \
  -fallow-argument-mismatch \
  -o "$build_dir/vaspberry" \
  "$repo_root/vaspberry_gfortran_serial.f" \
  -llapack -lblas

mkdir -p "$result_dir"
cd "$result_dir"

set +e
# The legacy Fortran CLI parses -f with list-directed input, where a leading
# slash terminates the value. Use this relative path instead of $wavecar.
"$build_dir/vaspberry" \
  -f ../WAVECAR -o NFIELD -z2 1 \
  -kx 12 -ky 12 -s 2 -ii 1 -if 10 \
  >fortran.log 2>&1
fortran_status=$?

python3 "$repo_root/tools/wavecar_z2.py" "$wavecar" \
  --nx 12 --ny 12 --occupied-bands 10 \
  --axis both --output-dir . --plot \
  >python.log 2>&1
python_status=$?
set -e

if (( fortran_status != 0 || python_status != 0 )); then
  echo "Bi Z2 validation did not pass; inspect results-v1.1.1 logs." >&2
  echo "Fortran status: $fortran_status; Python status: $python_status" >&2
  exit 1
fi

grep -Fq '# result_status=PASS' Z2_FIELD.csv
grep -Fq '# numerical_self_consistency_checks_pass=1' Z2_FIELD.csv
grep -Fq '# reportable_invariant=0' Z2_FIELD.csv
test ! -e INCOMPLETE
test ! -e Z2_FIELD.invalid.csv
test ! -e Z2_FIELD.tmp

python3 - z2_diagnostics.json <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    result = json.load(handle)
if result.get("valid") is not True or result.get("z2") not in (0, 1):
    raise SystemExit("guarded Wilson-loop result is not reportable")
print(f"validated Z2={result['z2']}")
PY

echo "Bi Z2 validation passed; outputs are in $result_dir"
