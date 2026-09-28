# Graphene: projected-spin Kubo proxy and Chern numbers

Use the real SOC graphene inputs from the [parent example](../README.md).
This example compares the native **local differential curvature proxy** with
native **spin-sector Chern numbers** of the same selected energy subspace,
evaluated with the **Fukui–Hatsugai–Suzuki (FHS) method**, called the **Fukui method** below.
It includes all occupied spinor bands 1–8 and an optional local calculation for
the top occupied pair 7–8.

The `spin-kubo` task uses canonical momentum to approximate the projector's
k-derivative in the WAVECAR pseudo Gram metric. Its headers state the omitted
PAW, nonlocal/SOC-vertex, and basis-derivative terms. It is not the conventional
spin-current Kubo spin Hall conductivity. Python draws completed CSVs; the
native Fortran/MPI program computes every curvature and integral below.

## Run the native calculations

Follow the parent example to build with Intel ifx/Intel MPI and produce
`results/graphene-spin-chern/mesh12/` and `path/` using ordinary VASP.
From the repository root:

```sh
VB="$PWD/build/vaspberry-ifx-mpi"
WORK="$PWD/results/graphene-spin-chern"
EXAMPLE="$PWD/examples/materials/graphene-spin-chern"

# Any ordered k list or symmetry path: point curves, without a BZ integral.
(cd "$WORK/path" && mpiexec -n 4 "$VB" --task spin-kubo --bands 1:8)

# Full periodic mesh: the same point calculation plus the raw mesh sum.
(cd "$WORK/mesh12" && mpiexec -n 4 "$VB" --task spin-kubo --bands 1:8 --mesh 12,12)

# Independently evaluate the same occupied-sector Chern numbers with FHS.
(cd "$WORK/mesh12" && mpiexec -n 4 "$VB" --task spin-chern --bands 1:8 --mesh 12,12)
```

Matching `WAVECAR` and `OUTCAR` in each working directory and Cartesian `z` are
the defaults. Use `--wavecar FILE`, `--outcar FILE`, and `--spin-axis z` to make
them explicit. `--output NAME` changes the file prefix. The native program
refuses to overwrite an existing result, so choose a fresh prefix or folder
for repeats. GNU serial and MPI alternatives are described in the parent guide;
reference execution used GNU serial, with independent implementation checks.

| File | Native content |
|---|---|
| `SPIN_KUBO.csv` | Local sector and parent curvature vectors, parent-projected and spin-mixing terms, path distance, gaps and conditioning |
| `SPIN_KUBO_SPECTRUM.csv` | Projected Pauli spectrum at every stored point |
| `SPIN_KUBO_INTEGRAL.csv` | Raw noninteger C_est sums; created only with explicit `--mesh` |
| `SPIN_CHERN.csv`, `SPIN_BERRY.csv`, `SPIN_SPECTRUM.csv` | Chern numbers evaluated with the Fukui method, plaquette fluxes, and projected spectrum |

The native point file supplies xy, yz, zx curvature in Å². Fractional and Cartesian
k coordinates and cumulative path distance are included. These plain CSVs can
be used directly in other plotting or analysis tools.

## Draw your outputs or replay the saved figure

```sh
python3 "$EXAMPLE/kubo/plot.py" \
  --mesh-dir "$WORK/mesh12" --path-dir "$WORK/path" \
  --poscar "$EXAMPLE/inputs/POSCAR" --output-dir "$WORK/kubo-figure"
```

The mesh directory also supplies the Fukui method outputs. If both folders contain
`PAIR_KUBO.csv`, the pair curves and any completed `PAIR_CHERN.csv` are included
automatically. Detailed folder and prefix overrides remain available in `--help`.

For the complete saved four-panel figure, including the optional pair curve:

```sh
python3 "$EXAMPLE/kubo/plot.py" --input-dir "$EXAMPLE/kubo/reference" \
  --poscar "$EXAMPLE/inputs/POSCAR" \
  --title 'Graphene: occupied and selected-pair spin sectors' \
  --output-dir "$WORK/kubo-reference-figure"
```

The signed-log scale exposes the very large K-point values without clipping
or changing them. Lines connect the stored path samples; they do not resolve
the actual microscopic width of the peak. Full-BZ proxy values are drawn as
constant cells centered on native k points. Fluxes from the Fukui method are shown on their
own plaquettes in radians. The two maps have different meanings and units.
`comparison.csv` retains full numerical values, while the figure uses a finite
number of significant digits for display.

## Interpret this numerical example

For the 12×12 occupied subspace, the raw proxy gives
`C_est_spin=1.1359066841142153e12`; the Fukui method gives spin Chern number `C_spin=1`.
The sub-micro-eV gap makes direct point sampling particularly unsuitable for
integrating this sharply localized response on a coarse uniform mesh. Missing
derivative terms and the finite 16-band source are additional limitations.
The large cancellation residual in the summed proxy charge response is also
retained; it is not a physical prediction. Nothing is rounded to an integer,
rescaled, or adjusted to match the Chern number. A complete calculation and a converged
physical response are separate requirements.

The optional pair command is:

```sh
(cd "$WORK/path" && mpiexec -n 4 "$VB" --task spin-kubo --bands 7:8 --output PAIR)
(cd "$WORK/mesh12" && mpiexec -n 4 "$VB" --task spin-kubo --bands 7:8 --mesh 12,12 --output PAIR)
```

At the stored points the pair is energy/spin separated and the local proxy
succeeds. However, its actual full-BZ calculation with the Fukui method fails the link guard
(minimum singular value 2.35e-10 at k4). Its global isolation and sampling remain
unresolved, and **no Chern number is assigned to this graphene pair**.
The failure is retained in `reference/pair-fukui/`. The
[Bi example](../../bi-spin-hall/spin-chern-kubo/README.md) provides a positive
selected-pair Chern-number comparison. See the saved native headers and
`reference/provenance.json` for exact scope, checksums, and source identities.

## Local sampling and solver controls

Resolve the intrinsic-SOC valley peak, compare actual-wavefunction polygon flux with local Kubo quadrature, and check electronic convergence. See the [convergence example](convergence/).
