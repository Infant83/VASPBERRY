# Bi: occupied and selected-pair spin sectors

This example uses the same real buckled Bi geometry and self-consistent charge
density as the [parent Bi example](../README.md). It compares the occupied
space 1–10 with the selected top pair 9–10. Native `spin-kubo` produces local
projected-spin curvature from a canonical-momentum derivative proxy; native
`spin-chern` independently computes full-BZ Fukui invariants. This geometric
sector calculation is distinct from conventional spin-current Hall transport.

On the saved 6×6 mesh, occupied bands 1–10 give Fukui `C+=1,C−=-1,Cspin=1`.
The selected pair 9–10 gives `C+=2,C−=-2,Cspin=2`. The latter describes that pair;
it is not the occupied system's invariant or its Z2 index. One coarse mesh is
not evidence of mesh convergence. The pointwise proxy sums remain raw and
noninteger even where a Fukui integer is available.

## Prepare and run ordinary VASP

Build the current source with Intel ifx and Intel MPI as in
[BUILD.md](../../../../docs/BUILD.md). From the repository root, define:

```sh
VB="$PWD/build/vaspberry-ifx-mpi"
EXAMPLE="$PWD/examples/materials/bi-spin-hall"
WORK="$PWD/results/bi-spin-chern-kubo"
VASP_NCL=/absolute/path/to/your/licensed/vasp_ncl
BI_POTCAR=/absolute/path/to/your/licensed/PAW-PBE/Bi/POTCAR
mkdir -p "$WORK"
```

The default preparation restores the checksum-verified supplied SCF density;
its provenance is documented by the parent example. To use a newly converged
matching SCF calculation, pass `--charge /path/to/SCF/CHGCAR` to both commands.
Use the same density for path and mesh.

```sh
python3 "$EXAMPLE/prepare_vasp.py" --stage wavecar --mesh 6 6 --nbands 48 \
  --potcar "$BI_POTCAR" --output-dir "$WORK/mesh6"
python3 "$EXAMPLE/spin-chern-kubo/prepare_path.py" --nbands 48 \
  --potcar "$BI_POTCAR" --output-dir "$WORK/path"

(cd "$WORK/mesh6" && "$VASP_NCL" > vasp.stdout 2> vasp.stderr)
(cd "$WORK/path" && "$VASP_NCL" > vasp.stdout 2> vasp.stderr)
```

Use your VASP installation's normal MPI launcher and stack/thread settings as
needed. The local reference ordinary-VASP path ran serially with one thread
and a 64 MiB stack. The path is 49 points Γ–K–M–Γ with nodes 1, 17, 33, 49. For this
positive 60° primitive basis, K=(1/3,2/3),M=(0,1/2). The path helper preserves the
parent physical INCAR, POSCAR, POTCAR and CHGCAR; only the explicit KPOINTS list
changes from its ordinary full-mesh preparation.

## Compute curves and full-BZ comparisons with native MPI

```sh
(cd "$WORK/path" && mpiexec -n 4 "$VB" --task spin-kubo --bands 1:10)
(cd "$WORK/mesh6" && mpiexec -n 4 "$VB" --task spin-kubo --bands 1:10 --mesh 6,6)
(cd "$WORK/mesh6" && mpiexec -n 4 "$VB" --task spin-chern --bands 1:10 --mesh 6,6)

(cd "$WORK/path" && mpiexec -n 4 "$VB" --task spin-kubo --bands 9:10 --output PAIR)
(cd "$WORK/mesh6" && mpiexec -n 4 "$VB" --task spin-kubo --bands 9:10 --mesh 6,6 --output PAIR)
(cd "$WORK/mesh6" && mpiexec -n 4 "$VB" --task spin-chern --bands 9:10 --mesh 6,6 --output PAIR)
```

Each directory must contain its own matching actual WAVECAR and OUTCAR.
The Cartesian spin axis defaults to z. `SPIN_KUBO.csv` and `PAIR_KUBO.csv`
contain the path coordinate, all three Cartesian curvature components,
spin-mixing decomposition and local quality checks. `*_KUBO_INTEGRAL.csv`
contains raw full-mesh sums only when `--mesh` is supplied. The separate
Fukui `*_CHERN.csv` and `*_BERRY.csv` contain the invariants and plaquette flux.
The native program writes `result_status=PASS` only after its checks pass.
Choose new output prefixes before rerunning; existing files are not overwritten.

## Plot the completed CSVs

The reusable plotting program accepts arbitrary matching native output folders:

```sh
PLOT="$PWD/examples/materials/graphene-spin-chern/kubo/plot.py"
python3 "$PLOT" \
  --mesh-dir "$WORK/mesh6" --path-dir "$WORK/path" \
  --poscar "$EXAMPLE/inputs/POSCAR" --k-fractional 0.3333333333333333 0.6666666666666666 \
  --title 'Bi: occupied and selected-pair spin sectors' --output-dir "$WORK/figure"
```

The plot finds Fukui outputs in the mesh directory and automatically includes
the `PAIR_*` outputs when present. Folder/prefix overrides are listed by `--help`.

To replay the saved figure using only public CSVs:

```sh
python3 "$PLOT" --input-dir "$EXAMPLE/spin-chern-kubo/reference" \
  --poscar "$EXAMPLE/inputs/POSCAR" --k-fractional 0.3333333333333333 0.6666666666666666 \
  --title 'Bi: occupied and selected-pair spin sectors' --output-dir "$WORK/reference-figure"
```

Python reads and plots the native results without calculating topology.
Curvature uses a signed-log display; full-precision comparisons are saved in
`comparison.csv`. Straight lines join the actual path samples, and constant
native point-cell/plaquette values are clipped to the first Brillouin zone.

For this 6×6 source the occupied proxy has `C_est_spin≈1.9997748`, while the
pair proxy has `C_est_spin≈−118.76613`; the corresponding Fukui values are 1
and 2. These discrepancies are retained. Canonical momentum is an approximate
derivative, 48 stored bands are finite, and this mesh does not establish
quadrature convergence. Do not round the proxy, interpret it as a conventional
spin Hall conductivity, or claim full PAW validation. Native metadata states
these limitations explicitly. The modern same-density data and the public
2016 Bi fixture are different calculations and are not combined as convergence
points. Exact source and output hashes are in `reference/provenance.json`.

At shared path/mesh k points, bands 1–32 agree within 6e-11 eV, but the
highest stored empty band differs by up to 0.000473 eV. Reaching the VASP
total-energy stopping criterion does not establish convergence of every empty
eigenvalue. The raw 48-band proxy is retained with this additional limitation;
see `reference/empty-band-diagnostics.json`.

## Mesh and intermediate-state convergence

Separate mesh refinement, the retained intermediate-state cutoff, and VASP source-state quality on the same Hamiltonian. See the [convergence example](convergence/).
