# Graphene: resolving the intrinsic-SOC curvature peak

This follow-up studies the same fixed planar graphene calculation as the
[parent example](../../README.md): 520 eV, 16 spinor bands, occupied bands 1–8,
and the same converged SCF density. It adds actual VASP points within
`1e-9` to `1e-5 Å^-1` of K and K′. The gap and curvature are neither enlarged
nor replaced by a model.

The outputs distinguish three calculations:

| Quantity | Calculation and scope |
|---|---|
| Local curvature | Native VASPBERRY `spin-kubo`, with its canonical-momentum derivative approximation |
| Polygon flux | Independent validation from overlaps of the actual WAVECAR projected-spin sectors; a small region around one valley |
| Dashed Dirac curves | A diagnostic local model fitted to the measured energy dispersion; not an additional VASP result |

Neither a local polygon flux nor a radial point sum is a full-BZ Chern number.
The separate occupied-space FHS calculation in the parent example gives
`C_plus=1`, `C_minus=-1` on the complete periodic BZ.

## Reproduce the ordinary VASP and native calculations

Complete the parent SCF preparation first. Use the same licensed carbon POTCAR
and SCF directory. From the repository root, with the Intel MPI executable
built as described in [the main example](../../README.md):

```sh
VB="$PWD/build/vaspberry-ifx-mpi"
WORK="$PWD/results/graphene-spin-chern"
EXAMPLE="$PWD/examples/materials/graphene-spin-chern/kubo/convergence"

python3 "$EXAMPLE/prepare.py" --scf-dir "$WORK/scf" \
  --potcar /path/to/C/POTCAR --output-dir "$WORK/local"
(cd "$WORK/local" && mpiexec -n 4 vasp_ncl)
(cd "$WORK/local" && mpiexec -n 4 "$VB" --task spin-kubo --bands 1:8)
```

`prepare.py` writes the ordinary VASP input files and `points.csv`, which names
both valleys, the Cartesian radii and angles, and each fractional k point. It
checks the SCF provenance and preserves the exact generating structure in the
CHGCAR header. The default has 146 points: both centers and nine rings of eight
angles per valley. These are local points, so do **not** add `--mesh` to this run.

VASPBERRY reads `WAVECAR` and matching `OUTCAR` from that directory and writes
`SPIN_KUBO.csv` and `SPIN_KUBO_SPECTRUM.csv`. The first file already contains the
curvature, selected-space gaps, and spin-sector decomposition. No Python step
is needed to obtain this native result.

The optional independent local overlap check and figure are:

```sh
python3 "$EXAMPLE/analyze.py" --input-dir "$WORK/local" \
  --output-dir "$WORK/local-analysis"
python3 "$EXAMPLE/plot.py" --input-dir "$WORK/local-analysis" \
  --output-dir "$WORK/local-figure"
```

`analyze.py` uses the actual coefficients to form orthonormal projected-spin
sectors and sums the phases of counterclockwise triangles from each valley
center to the boundary. It compares four- and eight-vertex polygons. It writes
`radial.csv`, `loops.csv`, and a source-hash `summary.json`. These diagnostic
polygon overlaps are computed by Python; they are not presented as a new
native task or as a native full-BZ `spin-chern` run.

`plot.py` draws angular means of the finished native samples, adds clearly
labeled dispersion fits, and writes PNG/PDF/SVG figures plus `model_diagnostics.csv` and
`integration_diagnostics.csv`. It requires Matplotlib. The local native point
integral uses the angular average and a finite trapezoid rule in log radius;
the radial sample count and incomplete spatial domain are retained explicitly.

## Check electronic convergence separately

The reference repeat samples both centers and three rings near the peak width,
with tighter EDIFF and more Davidson iterations:

```sh
python3 "$EXAMPLE/prepare.py" --scf-dir "$WORK/scf" \
  --potcar /path/to/C/POTCAR --output-dir "$WORK/local-strict" \
  --radii 3e-8,1e-7,3e-7 --ediff 1e-11 --nelmin 40
(cd "$WORK/local-strict" && mpiexec -n 4 vasp_ncl)
(cd "$WORK/local-strict" && mpiexec -n 4 "$VB" --task spin-kubo --bands 1:8)
python3 "$EXAMPLE/analyze.py" --input-dir "$WORK/local-strict" \
  --output-dir "$WORK/local-strict-analysis"
```

Changing EDIFF/NELMIN checks numerical eigenstate convergence at the same
Hamiltonian and same q points. It does not establish convergence of the SCF
density, cutoff, pseudopotential, retained empty bands, or omitted velocity
terms. Change one of those independently when studying physical accuracy.

## Saved measurements and integration strategy

The saved reference figure can be regenerated without licensed VASP files:

```sh
python3 "$EXAMPLE/plot.py" --input-dir "$EXAMPLE/reference" \
  --output-dir "$WORK/local-reference-figure"
```

The native CSVs remain the original unmodified files. Public references contain
plain data and source identities; licensed POTCAR, CHGCAR, WAVECAR and executables
are not included. See `reference/provenance.json` for their hashes and run
identities.

For a narrow-gap system, resolve the gap-derived width in Cartesian reciprocal
space first. Refine radial nodes near that width and compare the actual
WAVECAR polygon flux with the local point quadrature. A useful full-BZ adaptive
integration must also include the outer BZ with non-overlapping regions and
correct area weights. Agreement of a diagnostic Dirac fit is evidence about
the local peak, not a substitute for that integration or the omitted physical
velocity matrix elements.

## Measured resolution and approximation error

The completed 146-point local run gives the following values for K; K′ agrees
at the reported scale. These refer to the occupied projected-spin sectors.

| Quantity | K result | Interpretation |
|---|---:|---|
| Actual center gap | 0.912014 µeV | Same intrinsic-SOC Hamiltonian and SCF density |
| Local energy slope | 5.478965 eV Å | Fit to outer measured energy rings |
| Gap-derived width q₀ | 8.32287 × 10⁻⁸ Å⁻¹ | Diagnostic m/v, with m = gap/2 |
| Native center spin curvature | 6.82161 × 10¹³ Å² | Canonical-momentum proxy |
| Actual tiny-loop flux / polygon area, r = 10⁻⁹ Å⁻¹ | 7.21752 × 10¹³ Å² | Independent WAVECAR estimate of local pseudo curvature |
| Native / tiny-loop estimate | 0.94515 | About 5.5% smaller here; no correction factor is applied |
| Actual 8-vertex flux / 2π at r = 10⁻⁵ Å⁻¹ | 0.4956108 | One valley polygon, not the full BZ |
| Native local radial sum, five retained radii | 0.5208372 | Coarse subset of measured points |
| Native local radial sum, all nine radii | 0.4702077 | Refined finite point sum; convergence is not certified |

The native radial shape follows the diagnostic massive-Dirac shape closely:
the native/model angular-mean ratio ranges from 0.945055 to 0.945092 across both
valleys and all measured radii. The model's **positive chirality is a stated
choice matching the observed sector**; energy dispersion alone does not fix its
sign. The comparison supports the local peak shape for this calculation while
exposing the amplitude difference. It does not validate the exact PAW velocity
or establish a universal 5.5% correction.

Summing the two actual 8-vertex polygons gives 0.9912216. Their finite radius and
polygon boundaries exclude part of the BZ, so this is not rounded to one or
reported as a separately computed Chern number. Four vertices give a different
inscribed polygon; the 4→8 comparison measures boundary discretization as well
as angular sampling, and is not a proof of full angular convergence.

The coarse uniform 12×12 grid puts a huge weight on a peak whose measured width
is only about 8.3×10⁻⁸ Å⁻¹. This accounts for the enormous point-sum sensitivity;
the independently observed ~5.5% local amplitude difference is a separate
approximation issue. In particular, refinement does not imply that a canonical
Kubo integral must become the spin Chern number evaluated with FHS.

On the same 16-band WAVECAR, `--sum-bands 10`, `12`, and the default `16` give
maximum xy-sector relative changes of 1.13×10⁻¹⁰ (10→16) and 1.51×10⁻¹¹ (12→16).
This isolates the retained intermediate-state sum while keeping the VASP
Hamiltonian and stored wavefunctions fixed. It does not test additional source
bands above 16. The two-state conduction pair 9–10 is retained together.

```sh
(cd "$WORK/local" && mpiexec -n 4 "$VB" --task spin-kubo --bands 1:8 \
  --sum-bands 10 --output SUM10)
(cd "$WORK/local" && mpiexec -n 4 "$VB" --task spin-kubo --bands 1:8 \
  --sum-bands 12 --output SUM12)
```

The additional files are `SUM10_KUBO.csv`, `SUM10_KUBO_SPECTRUM.csv`, and the
corresponding `SUM12` files. The selected projector remains bands 1–8; only the
external intermediate-state sum changes. A cutoff that splits a numerically
degenerate boundary is rejected by VASPBERRY.

The 50-point electronic repeat converged in 47 Davidson iterations. Tightening
EDIFF from 10⁻⁹ to 10⁻¹¹ eV and increasing NELMIN from 20 to 40 changed the
matched gaps by at most 7.57×10⁻¹¹ eV (relative 7.81×10⁻⁵), the native xy spin
curvature by at most 0.01561%, and the actual local spin flux / 2π by at most
1.52×10⁻⁹. These observed solver changes are much smaller than the ~5.5%
center-curvature amplitude difference. The stricter run ended normally; the
finite total-energy numerical fluctuations and complete iteration history are
preserved in the research records. `reference/solver-comparison.json` contains
the exact metrics, and `reference/strict/` retains the corresponding native and
independent diagnostic CSVs.
