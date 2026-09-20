# Bi: Z₂ invariant and integer field

Buckled honeycomb Bi has **Z₂ = 1** for the occupied SOC bands 1–10 in this
12 × 12 example. VASPBERRY evaluates the Fukui–Hatsugai integer field and
compares the parities of two complementary half Brillouin zones.

## Input and preparation

Run from the repository root:

```bash
make serial
python3 -m pip install -r requirements-transport.txt
python3 examples/fetch_inputs.py bi --output-dir results/inputs/bi
```

The public VASP `WAVECAR` contains 144 points on a full Γ-centered mesh and
18 SOC spinor bands. The [structure](../../Bi_Z2/inputs/POSCAR),
[band energies](../../Bi_Z2/archive-2016-run/EIGENVAL) and
[material calculation notes](../../Bi_Z2/README.md) accompany the input.

## Run VASPBERRY

```bash
repo_dir="$PWD"
mkdir results/bi-z2
(
  cd results/bi-z2
  "$repo_dir/build/vaspberry-gfortran" \
    -f "$repo_dir/results/inputs/bi/WAVECAR" -o NFIELD -z2 1 \
    -kx 12 -ky 12 -s 2 -ii 1 -if 10 > fortran.log
)
```

`-z2 1` selects the integer-field method, `-s 2` reads SOC spinors, and
`-ii 1 -if 10` selects the complete occupied subspace. The mesh dimensions
must match the WAVECAR. For MPI, use `make mpi` and prefix the MPI executable
with `mpiexec -n 4`.

The output `Z2_FIELD.csv` reports:

| Quantity | Reference |
|---|---:|
| Z₂ invariant | **1** |
| Half-zone integer sums | −3 and +3 |
| Half-zone parities | 1 and 1 |
| Minimum occupied–empty direct gap | 0.5924485 eV |
| Sampled global gap | 0.5100444 eV |
| Minimum link singular value | 0.815455 |

Use a result only when `result_status=PASS`, `reportable_invariant=1`, and the
two half-zone parities agree. The [method guide](../../../docs/Z2_FUKUI_HATSUGAI.md)
explains the numerical checks.

## Plot the integer field

```bash
python3 examples/features/z2/run.py \
  --plot-only results/bi-z2/Z2_FIELD.csv \
  --poscar examples/Bi_Z2/inputs/POSCAR \
  --figure results/bi-z2/nfield.png
```

![Bi Z2 integer field in its Cartesian first Brillouin zone](reference/figure.png)

The axes are Cartesian **kx and ky in Å⁻¹**, with equal scales and the
hexagonal Wigner–Seitz first Brillouin zone. Each tile shows its computed
integer n(k); the categorical colors are not interpolated. K is `(1/3, 2/3)`
and K′ is its time-reversed partner.

**The integer field is gauge and branch dependent.** Its local appearance is
not a measurable Berry-curvature distribution. Z₂ comes from the agreed
half-zone parity in the original calculation; folding the plot into the
hexagon does not redefine the half zones used for that calculation.

[Full field CSV](reference/Z2_FIELD.csv) ·
[Native n-field output](reference/NFIELD.dat) ·
[Numerical summary](reference/summary.csv) ·
[Figure PDF](reference/figure.pdf)

To calculate, validate and plot in one command:

```bash
python3 examples/features/z2/run.py \
  --wavecar results/inputs/bi/WAVECAR --output-dir results/bi-z2-checked
```

## Apply to your material

Use a nonmagnetic, time-reversal-symmetric insulator with a full, unshifted,
even SOC mesh and `ISYM=-1`. Select a fixed even-dimensional occupied bundle
separated from the unoccupied states throughout the BZ. Replace the input,
mesh and band range in the VASPBERRY command and use the matching POSCAR for
the figure. [VASP input templates](../../Bi_Z2/inputs/) provide a starting point.

Repeat the gap and invariant calculation on denser meshes. Passing numerical
checks does not independently establish physical time-reversal symmetry or
mesh convergence. WAVECAR overlaps omit PAW augmentation. The public Bi
wavefunctions reproduce the post-processing calculation; the original SCF
provenance is incomplete. Historical results remain in the
[material reference](../../Bi_Z2/reference-v1.2.0-12x12/) and
[supplementary reference collection](../../../docs/REFERENCE_MATERIALS.md).
