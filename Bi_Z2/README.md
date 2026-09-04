# Bi Z2 validation case

Status: the VASP reference data and the spinor `WAVECAR` are available, but no
validated v1.1.1 Z2 result is committed yet. The tracked `WAVECAR` is managed
with Git LFS.

The files under [`legacy-pre-1.1.1/`](legacy-pre-1.1.1/) were produced by the
previous implementation. Its top and bottom half-zone candidates disagree
modulo 2. They are retained only as regression evidence and must not be cited
as the material's Z2 invariant.

`INCAR` and `POSCAR` received label-only corrections (`MoTe2`/`nanoribbon` to
the Bi buckled-honeycomb validation case); no numerical input value was
changed. VASP's license-restricted `POTCAR` is deliberately excluded. Its
non-proprietary provenance is recorded in
[`PSEUDOPOTENTIAL.md`](PSEUDOPOTENTIAL.md).

## Recalculate with v1.1.1

Install Git LFS, GNU Fortran, BLAS/LAPACK, Python 3, and the Python
requirements. From the repository root:

```bash
git lfs pull --include='Bi_Z2/WAVECAR'
python3 -m pip install -r requirements-transport.txt
./Bi_Z2/run_z2.sh
```

The script writes candidate outputs and logs to `results-v1.1.1/`. It runs both
the corrected Fortran reconstruction and the independent guarded Wilson-loop
workflow.

Do not publish a result as current unless all of the following hold:

1. `Z2_FIELD.csv` contains `result_status=PASS` and the top/bottom candidates
   agree modulo 2.
2. No `INCOMPLETE`, `Z2_FIELD.invalid.csv`, or staged temporary product
   remains.
3. `z2_diagnostics.json` contains `"valid": true` and `z2` is 0 or 1.
4. The outputs record the code revision and the full 12 x 12 spinor mesh.

The Fortran PASS establishes numerical self-consistency of its legacy
reconstruction only; its CSV intentionally records `reportable_invariant=0`.
The guarded Python result must also pass before a Z2 invariant is reported.

The raw pointwise n-field is gauge- and logarithm-branch-dependent, so
point-by-point mirror symmetry is not itself the acceptance criterion. The
relevant checks are half-zone parity agreement, time-reversal-odd physical
Berry flux, link conditioning, energy pairing, Kramers splitting, projector
residuals, the occupied-unoccupied gap, and mesh convergence.
