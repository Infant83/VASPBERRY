# Bi: Z₂ invariant and integer field

The [MoS₂–Bi comparison](comparison/) now contrasts a native **Z₂ = 0**
MoS₂ result with the fresh **Z₂ = 1** Bi result used in the technical report.
It includes the [MoS₂ native command and VASP preparation](mos2/) and a
plot-only helper for paired reduced-coordinate fields. The step-by-step
Bi tutorial below retains its distinct historical input and references.

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

## 1. Run native Fortran on WAVECAR

```bash
repo_dir="$PWD"
mkdir results/bi-z2
(
  cd results/bi-z2
  "$repo_dir/build/vaspberry" \
    --wavecar "$repo_dir/results/inputs/bi/WAVECAR" --output NFIELD --task z2 \
    --mesh 12,12 --spinor 2 --bands 1:10 > fortran.log
)
```

`--task z2` selects the integer-field method, `--spinor 2` reads SOC spinors, and
`--bands 1:10` selects the complete occupied subspace. The mesh dimensions
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

## 2. Plot the native integer field

```bash
python3 examples/features/z2/run.py \
  --plot-only results/bi-z2/Z2_FIELD.csv \
  --poscar examples/Bi_Z2/inputs/POSCAR \
  --figure results/bi-z2/nfield.png
```

The explicit `--plot-only` mode reads the completed CSV and does not run
VASPBERRY again.

![Bi Z2 integer field in reduced reciprocal coordinates](reference/figure.png)

The dimensionless reduced coordinates satisfy **k = q₁b₁ + q₂b₂**, where b₁
and b₂ are reciprocal lattice vectors. The square covers −1/2 ≤ q₁, q₂ ≤ 1/2,
with opposite edges periodically identified. Each tile shows the computed
integer n(k) of one native plaquette, without interpolation. The horizontal
line at q₂ = 0 separates the two half zones: the upper sum is **−3** and the
lower sum is **+3**, both giving **Z₂ = 1**.

**The integer field is gauge and branch dependent.** Its local appearance is
not a measurable Berry-curvature distribution. Z₂ comes from the agreed
half-zone parity in the original calculation. The reduced-coordinate display
shows those summation domains directly.

[Full field CSV](reference/Z2_FIELD.csv) ·
[Native n-field output](reference/NFIELD.dat) ·
[Numerical summary](reference/summary.csv) ·
[Figure PDF](reference/figure.pdf)

## Optional reproduction helper

To repeat the native calculation, validation and plotting in a fresh directory:

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
