# Separate occupied-bundle Chern calculation

The unchanged fresh 12 × 12 Bi WAVECAR was evaluated separately in ordinary Fukui mode, with occupied bands 1–10 and both SOC components:

```sh
vaspberry-gfortran -f WAVECAR -o BERRYCURV \
  -kx 12 -ky 12 -s 2 -ii 1 -if 10
```

The native result is **C = 0.000000**. A separate complex128 Python overlap calculation gives **C = 2.83 × 10⁻¹⁶** without rounding, with minimum link singular value **0.79868**. The native run took 32.00 s; the Python check took 20.20 s, each on one CPU. Both use the original occupied states, not the time-reversal reconstructed states of the Z₂ routine.

`BERRYCURV.dat` and `fortran.log` are unchanged native outputs. `plaquettes.csv` and `python-fukui.npz` retain the higher-precision independent loop flux and curvature, including numerical residuals below the native printed precision. `summary.json` records the input and binary identity with a portable command; `independent-check.json` gives the full-precision checks.

The backend uses raw pseudo-wavefunction overlaps and does not add PAW overlap augmentation. This separate Chern result is consistent with time reversal; C = 0 alone does not distinguish a trivial insulator from the nontrivial Z₂ phase.
