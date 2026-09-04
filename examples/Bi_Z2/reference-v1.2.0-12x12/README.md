# Bi 12 x 12 Z2 reference field

`Z2_FIELD.csv` is the schema-2 PASS output produced by VASPBERRY 1.2.0
from the tracked Bi spinor `WAVECAR` with occupied bands 1-10:

```text
-z2 1 -kx 12 -ky 12 -s 2 -ii 1 -if 10
```

It was captured from the serial result in
[Bi n-field Z2 validation run 33881876978](https://github.com/Infant83/VASPBERRY/actions/runs/33881876978)
for commit
[`aef4d16`](https://github.com/Infant83/VASPBERRY/commit/aef4d163996d47a618b289f5251f22d26fa71981).

The 12 x 12 fundamental mesh contains 144 plaquettes. Its complementary
half-zone sums are -3 and +3, so both parities are 1 and the calculation gives
`Z2=1` on this mesh. The CSV also records a total Chern residual of
`1.1593839628498794e-14`, a maximum wrapped time-reversal-reconstruction flux
residual of `1.8152146452621309e-14 rad`, and a minimum link singular value of
`0.8154547168424493`.

The PASS status covers the numerical self-consistency checks named in the
CSV. It does not independently establish the raw input's physical
time-reversal symmetry or band gap, include PAW augmentation in the overlap,
or demonstrate k-mesh convergence.
