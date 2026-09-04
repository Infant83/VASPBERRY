# 1H-MoS2 examples

This tree contains two distinct sampling workflows. Do not interchange their
inputs or numerical products.

## Full 12 x 12 Brillouin-zone map

The files in this directory are the original full-mesh Berry-curvature and
Chern example. `BERRYCURV.dat` contains 144 fundamental plaquettes repeated
over a 3 x 3 reciprocal display region. A full-mesh `WAVECAR` is not included.

To reproduce the calculation, generate a spinor `WAVECAR` on a full
Gamma-centered 12 x 12 x 1 mesh with symmetry reduction disabled, then run:

```bash
./build/vaspberry-gfortran \
  -f /path/to/full-mesh/WAVECAR \
  -kx 12 -ky 12 -ii 1 -if 18 -s 2
```

To redraw the stored curvature map:

```bash
python3 examples/1H-MoS2/contour.py
```

The plotting script resolves its input relative to its own location, so the
command can be issued from the repository root.

## K-Gamma-K' line data

`KPATH/1.scf/` contains a charge-density calculation, `KPATH/2.band/` contains
a 48-point K-Gamma-K' line calculation, and `KPATH/3.BC_kubo/` contains
band-resolved Kubo plot data on that line. The line-mode `WAVECAR` and
`EIGENVAL` are not a two-dimensional full-BZ mesh and cannot be used as Z2 or
Chern inputs.

The generic `KPATH/1.scf/sbatch_vasp.sh` template expects `VASP_BIN` to name
the local executable; scheduler resources remain site-specific.

## Data hygiene

All VASP `POTCAR` files have been removed from the current tree. Recreate them
locally according to [`PSEUDOPOTENTIAL.md`](PSEUDOPOTENTIAL.md). The old
`PdTe2` titles in user-facing `INCAR`, `POSCAR`, and `CONTCAR` files were
corrected to `1H-MoS2`; numerical input values were not changed. Raw VASP
outputs are retained unchanged as calculation records.
