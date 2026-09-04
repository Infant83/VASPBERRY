# Bi Z2 validation case

Status: the VASP reference data and spinor `WAVECAR` are available. A v1.1.1
recalculation of the 12 x 12 mesh fixes the legacy top/bottom disagreement,
but the independent Wilson-loop guards do not yet validate a material Z2
invariant. The tracked `WAVECAR` is managed with Git LFS.

The files under [`legacy-pre-1.1.1/`](legacy-pre-1.1.1/) were produced by the
previous implementation. Its top and bottom half-zone candidates disagree
modulo 2. They are retained only as regression evidence and must not be cited
as the material's Z2 invariant.

The reviewed recalculation is stored under
[`diagnostics-v1.1.1-12x12/`](diagnostics-v1.1.1-12x12/). The corrected Fortran
reconstruction gives matching top/bottom legacy candidates, 1 and 1, and
passes its numerical self-consistency checks. The guarded Python workflow also
gives candidate 1 on both loop axes, but records `valid=false` and `z2=null`.
This is a diagnostic snapshot, not a validated Z2 result.

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

## Next validation calculation

The unresolved WCC checks are confined to the first pump interval next to
Gamma on this 12 x 12 mesh. Recalculate a full Gamma-centered 24 x 24 x 1
mesh with the same structure, SOC, `ISYM=-1`, occupied bands 1--10, and at
least band 11 as the unoccupied sentinel. Use `PREC=Accurate`,
`EDIFF<=1e-8`, and `LREAL=.FALSE.`. Keep the licensed Bi PAW dataset local;
do not add `POTCAR` to the repository.

Compare the 12 x 12 and 24 x 24 guarded results without relaxing the current
MOVE, GAP, or time-reversal thresholds. If the first interval remains
under-resolved, continue to a 48 x 48 mesh or an adaptive Wilson-loop workflow
that inserts pump lines near Gamma.
