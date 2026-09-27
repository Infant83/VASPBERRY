# Bi spin-sector convergence: three separate controls

This follow-up varies the full Brillouin-zone mesh, the number of intermediate
states in `spin-kubo`, and the VASP source eigensystem separately. The structure,
PAW potential and fixed SCF density are identical to the [parent example](../../README.md).
The 2016 Bi dataset is not used in this convergence series.

`--bands 1:10` defines the occupied parent subspace. `--bands 9:10` defines only
the top occupied pair. `--sum-bands 48` retains source bands 1–48 in the response
sum, excluding the selected parent. It does **not** change the parent subspace,
source WAVECAR, VASP `NBANDS`, or the full-source energy-isolation check. Omitting
`--sum-bands` includes every stored source band. A cutoff that cuts a degenerate
boundary is rejected.

## Reproduce with native VASPBERRY

Build the current source with Intel ifx and Intel MPI following
[BUILD.md](../../../../../docs/BUILD.md). From the repository root:

```sh
VB="$PWD/build/vaspberry-ifx-mpi"
EXAMPLE="$PWD/examples/materials/bi-spin-hall"
WORK="$PWD/results/bi-spin-convergence"
VASP_NCL=/absolute/path/to/your/licensed/vasp_ncl
BI_POTCAR=/absolute/path/to/your/licensed/PAW-PBE/Bi/POTCAR
mkdir -p "$WORK"
```

Keep source `NBANDS=64` and the response cutoff 48 fixed while changing the mesh:

```sh
for n in 6 12 18; do
  python3 "$EXAMPLE/prepare_vasp.py" --stage wavecar --mesh "$n" "$n" \
    --nbands 64 --potcar "$BI_POTCAR" --output-dir "$WORK/mesh$n-source64"
  (cd "$WORK/mesh$n-source64" && "$VASP_NCL" > vasp.stdout 2> vasp.stderr)
  (cd "$WORK/mesh$n-source64" && mpiexec -n 4 "$VB" \
    --task spin-kubo --bands 1:10 --sum-bands 48 --mesh "$n,$n" --output OCC48)
  (cd "$WORK/mesh$n-source64" && mpiexec -n 4 "$VB" \
    --task spin-kubo --bands 9:10 --sum-bands 48 --mesh "$n,$n" --output PAIR48)
  (cd "$WORK/mesh$n-source64" && mpiexec -n 4 "$VB" \
    --task spin-chern --bands 1:10 --mesh "$n,$n" --output OCC)
  (cd "$WORK/mesh$n-source64" && mpiexec -n 4 "$VB" \
    --task spin-chern --bands 9:10 --mesh "$n,$n" --output PAIR)
done
```

Use the normal MPI launcher for your VASP installation. `OCC48_KUBO.csv` and
`PAIR48_KUBO.csv` contain pointwise Cartesian sector curvature; their
`*_KUBO_INTEGRAL.csv` files contain raw signed-area quadrature. `OCC_CHERN.csv`
and `PAIR_CHERN.csv` contain the separate Chern numbers evaluated with FHS. Each run requires its
own matching WAVECAR and OUTCAR. The source count and actual sum cutoff are
recorded separately in the native CSV metadata.

For intermediate-state convergence, make one 6×6 source with `NBANDS=80` and
reuse that unchanged WAVECAR for every response cutoff:

```sh
python3 "$EXAMPLE/prepare_vasp.py" --stage wavecar --mesh 6 6 --nbands 80 \
  --potcar "$BI_POTCAR" --output-dir "$WORK/mesh6-source80"
(cd "$WORK/mesh6-source80" && "$VASP_NCL" > vasp.stdout 2> vasp.stderr)
for n in 16 24 32 48 64 80; do
  (cd "$WORK/mesh6-source80" && mpiexec -n 4 "$VB" \
    --task spin-kubo --bands 1:10 --sum-bands "$n" --mesh 6,6 --output "OCC$n")
  (cd "$WORK/mesh6-source80" && mpiexec -n 4 "$VB" \
    --task spin-kubo --bands 9:10 --sum-bands "$n" --mesh 6,6 --output "PAIR$n")
done
```

To assess source eigenstate quality, compare `--sum-bands 48` between the
independent source64 and source80 calculations at the same mesh. Changing
`NBANDS` and the response cutoff together mixes two numerical controls.

## Read and plot the reference

At fixed source64 and response cutoff48, the actual results are:

| Mesh | Occupied proxy `C_est_spin` | Pair proxy `C_est_spin` | Occupied spin Chern number (FHS) | Pair spin Chern number (FHS) |
| --- | ---: | ---: | ---: | ---: |
| 6×6 | 1.9997800132 | −118.76612842 | 1 | 2 |
| 12×12 | 1.0382696725 | −29.02711625 | 1 | 2 |
| 18×18 | 1.0268699256 | −12.25375432 | 1 | 2 |

The occupied proxy changes by 0.01140 between 12×12 and 18×18 and remains
0.02687 above the spin Chern number. This is improved sampling behavior, not a
convergence certificate or a license to round. The pair proxy is clearly
unconverged. Its single sampled Γ point contributes −118.5620, −29.6405 and
−13.1736 on these meshes: the same large point value receives the shrinking
`1/N²` quadrature weight. These are **point weights**, not integrals over a
resolved Γ region. Common-point curvatures agree to less than `2.2e-8` of the
shared peak across the source64 meshes, supporting sampling sensitivity rather
than inconsistent electronic states as the source of this mesh change.
The saved [actual pair path](../reference/pair-path/SPIN_KUBO.csv) also changes
from `Ω⁺xy≈−11331.8 Å²` at Γ to `+138.5 Å²` at the next point, only
`0.05965 Å⁻¹` away. That path demonstrates rapid variation, but does not resolve
the peak width or supply a two-dimensional integral.

At a fixed 6×6 mesh and cutoff48, source64→source80 changes the occupied
integral by `1.6e-8` and the pair by `1.1e-6`. The first 48 eigenvalues in these
source64/source80 controls agree within `6.5e-11 eV`; however, among bands 49–64 the mismatch can
reach approximately `0.12 eV`. Therefore, sums retaining all stored upper empty
states are recorded as diagnostics, not as a converged improvement. On the
same source80 WAVECAR, cutoff32→48 changes the occupied value by about
`6.4e-4` and the pair by `−0.0040`; extending beyond48 also needs better empty
states. The fixed-source cutoff sweep and source-padding controls are both
preserved in the numerical table.

![Bi mesh and intermediate-band convergence](reference/figure.png)

The portable tables and full point data are in [reference/](reference/).
`convergence.csv` separates the mesh series, fixed-source cutoff series and
source-padding comparison. `source-quality.json` records common-k eigenvalue
differences and the final eigensolver residuals; total-energy convergence alone
does not establish the accuracy of every empty state.

```sh
python3 "$EXAMPLE/spin-chern-kubo/convergence/plot.py" \
  --input "$EXAMPLE/spin-chern-kubo/convergence/reference" \
  --output "$WORK/convergence-figure"
```

This optional plotting step uses Matplotlib, reads native-derived tables and writes
`figure.png`, `figure.pdf` and `figure.svg`. It does not rerun VASP or calculate
the Berry curvature. The figure displays raw estimates, including negative
and noninteger results, alongside the Chern numbers evaluated independently with FHS.

## Provenance and interpretation

The archived 6/12/18 calculations used partitions of the same complete mesh.
Each original partition was evaluated by native `spin-kubo` with its actual
matching OUTCAR. The reference aggregate verifies the complete periodic grid
and sums native point curvature with `(b1×b2)/(Nx Ny)/(2π)`. It is labelled
`POSTPROCESS_SUM_OF_NATIVE_POINT_CURVATURE`; it is not a synthetic native
integral file. Original source files and headers were not edited.

The independent NumPy FHS oracle uses the assembled complete WAVECAR and the
common spin frame verified across all original OUTCARs. Its outputs are marked
as independent validation, separately from the native full-mesh check. The
[original input pack](../../inputs/convergence/README.md) preserves the exact
partitions and eigensolver settings. The commands above give the simpler
ordinary full-mesh workflow with the same Hamiltonian.

An additional ordinary-VASP 12×12 source64 restart produced its own actual
complete OUTCAR and WAVECAR. Native full-mesh `spin-chern` confirmed the same
occupied/pair integers 1 and 2. The native cutoff48 integrals were
`1.03826968065` and `−29.02711596`, differing from the original partition
aggregates by `8.2e-9` and `3.0e-7`. The actual source, commands and CSVs are
saved in [native-mesh12-restart/](reference/native-mesh12-restart/), with
[comparison diagnostics](reference/native-restart-comparison.json). This check
uses real VASP output metadata; it does not combine fabricated OUTCAR headers.

No noninteger Kubo estimate is rounded into a Chern number. A stable spin Chern number
evaluated with FHS on these meshes establishes consistency within the sampled pseudo
wavefunction metric, while the curvature proxy still needs its own sampling,
source-state and operator-approximation checks. The pair invariant describes
bands 9–10 and is not the occupied system's Z2 index. See the full
[method and output contract](../../../../../docs/SPIN_KUBO.md).
