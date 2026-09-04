# Corrected v1.1.1 field reference: 12 x 12

This snapshot was generated from the tracked Bi `WAVECAR` by commit
[`fa4d14b`](https://github.com/Infant83/VASPBERRY/commit/fa4d14b71cfccf816376488d33bc528f14126aa7)
in [GitHub Actions run 33850537368](https://github.com/Infant83/VASPBERRY/actions/runs/33850537368).
It is retained as the corrected n-field regression reference.

The v1.1.1 CSV still contains the historical policy strings
`LEGACY_FUKUI_Z2_CANDIDATE` and `reportable_invariant=0`. Those strings record
the software policy at the time the file was generated; they are not
retrofitted in an archived numerical product. The current `-z2` interface
publishes the Fukui-Hatsugai n-field result as a formal invariant when all of
its own guards pass.

## Numerical result

| Quantity | Value |
|---|---:|
| top half-zone n-field sum | -3 |
| bottom half-zone n-field sum | +3 |
| top parity | 1 |
| bottom parity | 1 |
| common 12 x 12 Z2 result | 1 |
| total Chern residual | 1.15938e-14 |
| maximum wrapped TR-odd flux residual | 1.81521e-14 rad |
| maximum paired n-field residual | 0 |
| minimum raw link singular value | 0.815455 |

The two complementary half-zone parities agree, and all numerical
self-consistency checks recorded in `Z2_FIELD.csv` pass. This is a result on
the archived 12 x 12 mesh. It does not certify the material's magnetic/TRS
state, the occupied-unoccupied gap, PAW augmentation, or convergence with
mesh density; those checks are outside this archived field calculation.

The pointwise n-field is gauge- and logarithm-branch-dependent. Pointwise
mirror equality is not required. The relevant internal criterion is that the
two half-zone sums have the same parity, together with the CSV guards.

## Files

- `NFIELD.dat`: repeated-zone n-field output from the corrected source.
- `Z2_FIELD.csv`: the 144-plaquette fundamental mesh and diagnostics.
- `fortran.log`: captured VASPBERRY output from the Actions run.

Use [`../scripts/run_z2.sh`](../scripts/run_z2.sh) to repeat the
post-processing calculation with the current source.
