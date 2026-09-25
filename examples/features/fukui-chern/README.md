# Bi: occupied-band Berry curvature and Chern number

This calculation uses the VASP SOC wavefunctions of buckled honeycomb Bi on a
12 × 12 mesh. Bands 1–10 form the occupied subspace. Their total Chern number
is **C = 0**, consistent with time-reversal symmetry. The same occupied
subspace has **Z₂ = 1** in the [Z₂ example](../z2/).

## Input and preparation

Run from the repository root:

```bash
make serial
python3 -m pip install -r requirements-transport.txt
python3 examples/fetch_inputs.py bi --output-dir results/inputs/bi
```

The input is the public Bi `WAVECAR`: 144 k points, 18 SOC spinor bands,
400 eV plane-wave cutoff. Its [structure](../../Bi_Z2/inputs/POSCAR),
[VASP eigenvalues](../../Bi_Z2/archive-2016-run/EIGENVAL) and
[calculation record](../../Bi_Z2/README.md) describe the material.
The download command checks the input automatically.

## 1. Run native Fortran on WAVECAR

```bash
repo_dir="$PWD"
mkdir results/bi-fukui
(
  cd results/bi-fukui
  "$repo_dir/build/vaspberry" \
    --wavecar "$repo_dir/results/inputs/bi/WAVECAR" --output BERRYCURV \
    --task chern --mesh 12,12 --spinor 2 --bands 1:10 > fortran.log
)
```

`--mesh 12,12` specifies the mesh already present in WAVECAR. `--spinor 2` reads
SOC spinors. `--bands 1:10` selects the complete occupied band bundle;
its internal Kramers degeneracies do not require separate band invariants.
The result is written to `BERRYCURV.dat`, and the log reports C = 0.

For MPI, build with `make mpi` and replace the executable by
`mpiexec -n 4 "$repo_dir/build/vaspberry-mpi"`.

## 2. Plot the native curvature file

```bash
python3 tools/plot_berry_curvature.py \
  --input results/bi-fukui/BERRYCURV.dat \
  --poscar examples/Bi_Z2/inputs/POSCAR \
  --output results/bi-fukui/curvature.png \
  --title "Bi: occupied bands 1–10, C = 0"
```

The plot uses Cartesian **kx and ky in Å⁻¹**, equal axis scales, and the
hexagonal Wigner–Seitz first Brillouin zone. Each colored cell retains the
native plaquette-averaged curvature; no smoothing or interpolation is applied.
The color scale is centered on zero. K is defined as `(1/3, 2/3)` in reciprocal
coordinates, and K′ is its time-reversed partner.

![Bi occupied-subspace curvature in its Cartesian first Brillouin zone](reference/figure.png)

For this occupied Bi bundle, all curvature values round to zero in
`BERRYCURV.dat`, which prints four decimal places. The blank white map is the
expected result at that precision. For a finite, sign-changing curvature map,
see the [MoS₂ example](../fukui-berry-curvature/). The ±10⁻⁴ Å² color scale shows the output
resolution, not a resolved signal. The direct occupied–empty gap is
**0.5924485 eV** and the sampled global gap is **0.5100444 eV**.

[Native curvature](reference/BERRYCURV.dat) ·
[Plaquette values](reference/plaquettes.csv) ·
[Band edges](reference/band_edges.csv) ·
[Numerical summary](reference/summary.csv) ·
[Figure PDF](reference/figure.pdf)

## Optional reproduction helper

To repeat the native calculation, validation and plotting in a fresh directory:

```bash
python3 examples/features/fukui-chern/run.py \
  --wavecar results/inputs/bi/WAVECAR --output-dir results/bi-fukui-checked
```

## Apply to your material

Use a full uniform mesh with an isolated band or a fixed band bundle separated
from all other bands. Set the mesh, spinor setting and band range from your
VASP calculation, then pass its matching POSCAR to the plotter. A band path
cannot replace the two-dimensional mesh needed for a Chern integral.

For an insulating occupied bundle, check the gap across the BZ and repeat on
denser meshes. An integer on one mesh is not a convergence test. WAVECAR
overlap calculations use pseudo-wavefunctions without PAW augmentation.
The supplied Bi result reproduces post-processing; its original preceding
SCF calculation is not fully archived. See the [material notes](../../Bi_Z2/README.md)
and [supplementary references](../../../docs/REFERENCE_MATERIALS.md).
