<a id="spin-sector-kubo-curvature-on-a-path-and-a-full-mesh"></a>

# Spin-sector Kubo-formula Berry curvature on a path and a full mesh

Use native `--task spin-kubo` to calculate the Kubo-formula Berry-curvature approximation
of the positive and negative projected-spin sectors at each supplied k point.
It supports internally degenerate selected bands and includes the change of
the spin-sector basis with k. VASPBERRY performs the calculation in Fortran;
Python is optional for drawing the resulting CSVs.

This native task is available in VASPBERRY 1.6.1.
Its canonical-momentum derivative approximation is intended for resolved
curvature and comparisons. A separate [spin Chern number calculation](SPIN_CHERN.md) uses the
Fukui–Hatsugai–Suzuki (FHS) link-variable method, also called the Fukui method,
to evaluate the geometric sector Chern numbers on a complete periodic mesh.

## 1. Build VASPBERRY

With [Intel Fortran and Intel MPI](BUILD.md#intel-oneapi-on-linux):

```bash
source /opt/intel/oneapi/setvars.sh
make ifx-mpi
VB="$PWD/build/vaspberry-ifx-mpi"
```

Use the same MPI environment when compiling and executing. A GNU serial
alternative is `make serial`, `VB="$PWD/build/vaspberry"`, with the MPI launcher
omitted. No Python installation is required to execute either native task.

## 2. Calculate a k-path

Start with an ordinary VASP SOC/noncollinear calculation containing the desired
k points. In that calculation's directory, with matching WAVECAR and OUTCAR:

```bash
mpiexec -n 4 "$VB" --task spin-kubo --bands 1:10
```

The example range `1:10` selects all ten occupied states of the Bi example;
adapt it to the material. The defaults are `WAVECAR`, `OUTCAR`, Cartesian
spin axis `z`, and output prefix `SPIN`. Paths are relative to the directory
where the command runs. Add `--wavecar PATH --outcar PATH` for files elsewhere.
OUTCAR supplies the source spin frame and matching dimensions/lattice;
WAVECAR supplies the wavefunctions and automatic two-component detection.

This command accepts an arbitrary k list, including repeated path endpoints.
It does not integrate the path or assign a Chern number. It writes:

| Native output | Data and use |
|---|---|
| `SPIN_KUBO.csv` | k coordinates, cumulative path distance, positive/negative/parent curvature in Å², its two contributions, ranks and gap/metric diagnostics |
| `SPIN_KUBO_SPECTRUM.csv` | Eigenvalues of the projected Pauli operator at every k; inspect the separation from zero |

For an xy-plane plot, use `k_distance_A-1` against `omega_plus_xy_A2` and
`omega_minus_xy_A2`. The `yz` and `zx` components are also saved. The separate
`parent_projected_*` and `spin_mixing_*` columns add to each sector's total;
the second term is required for a general k-dependent projected-spin split.
CSV data can be read by Origin, MATLAB, Julia, pandas or another plotting tool.
Skip lines starting with `#`, but retain them with the data as provenance.

The path distance follows the input order. It is useful for a band path; the
same column on an unordered full mesh is not a symmetry-path coordinate.

<a id="3-compare-a-full-mesh-with-fukui"></a>

<a id="3-compare-full-mesh-results-with-chern-numbers-from-the-fhs-method"></a>

## 3. Compare full-mesh results with Chern numbers from the Fukui method

Use a **separate full periodic VASP mesh**, with the same Hamiltonian, selected
band space and analysis axis. In its directory:

```bash
mpiexec -n 4 "$VB" --task spin-kubo --bands 1:10 --mesh 12,12
mpiexec -n 4 "$VB" --task spin-chern --bands 1:10 --mesh 12,12
```

Explicit `--mesh NX,NY` validates the full mesh and adds
`SPIN_KUBO_INTEGRAL.csv`. Its `C_est_plus`, `C_est_minus`, `C_est_charge` and
`C_est_spin` are **raw finite-grid approximation integrals**, never rounded
or certified as integers. They contract all three Cartesian curvature
components with the oriented reciprocal plaquette area, so tilted cells do
not silently become xy-plane integrals.

The second command independently writes `SPIN_CHERN.csv`, `SPIN_BERRY.csv`
and `SPIN_SPECTRUM.csv`. Compare its geometric `C_plus` and `C_minus` to the
approximate integrals. A k-path alone cannot supply a two-dimensional Chern
number. Refine sampling and the retained source bands before interpreting
integral differences; nonlocal/SOC velocity and PAW corrections can remain.
Do not round a poorly resolved integral to obtain the expected invariant.

Both commands preserve existing files. Use a new directory or `--output NAME`
for a new run. Without an explicit `--mesh`, spin-kubo writes no integral file,
even when the supplied points happen to form a mesh.

## Select a band pair or another spin axis

```bash
mpiexec -n 4 "$VB" --task spin-kubo --bands 9:10 --output PAIR
```

In the Bi example this range selects a degenerate pair and gives rank-one
positive and negative spin sectors. A sector Chern number can be associated
with one spin-resolved branch only when the sector has rank one and is globally
well defined.
For generic spin mixing, a projected-spin vector need not coincide with one
nondegenerate energy eigenstate. A rank-five sector is a five-state subspace,
not a single numbered band.

Use `--spin-axis x`, `y`, `z`, or a Cartesian unit vector such as `0,0,-1`.
Reversing the axis exchanges the sectors. The same axis and band range must
be used for path plots and the full-BZ comparison.

The default `--energy-gap-tol 1e-8` is in eV; `--spin-gap-tol 1e-6` is the
dimensionless distance of a projected Pauli eigenvalue from zero. These are
numerical rejection thresholds, not convergence targets. Do not lower them
to assign topology across a physically unresolved crossing.

## Real examples and interpretation

To test the intermediate-band sum without changing the VASP eigenstates,
add `--sum-bands N`. For example, on an 80-band WAVECAR:

```bash
mpiexec -n 4 "$VB" --task spin-kubo --bands 1:10 --mesh 6,6 --sum-bands 32 --output B32
mpiexec -n 4 "$VB" --task spin-kubo --bands 1:10 --mesh 6,6 --sum-bands 64 --output B64
```

This sums over stored bands **1 through N**, excluding the selected group;
the default uses every stored band. `--bands` still selects the subspace
whose curvature is calculated. The cap must contain that group and at least
one external state, and cannot exceed source `NBANDS`. A cap that cuts an
unresolved external degeneracy is rejected. Parent energy isolation is
always checked against **all stored bands**, including those omitted from
the sum. The CSV records `source_nbands` and `sum_band_max` separately.
Increasing this cap tests sum truncation; increasing VASP `NBANDS` and
converging the source eigenstates is a separate check.

- [Graphene SOC](../examples/materials/graphene-spin-chern/kubo/): a genuine
  Γ–K–M–Γ path, opposite occupied-sector curvature, and a tiny-gap example
  showing why a coarse pointwise integral need not resemble a quantized Chern number.
- [Bi bilayer](../examples/materials/bi-spin-hall/spin-chern-kubo/): compare
  the complete occupied space with an isolated topmost degenerate pair. The
  invariants belong to different projectors and need not be equal.

These examples retain native CSVs, optional plotting commands and the actual
VASP source identities. `result_status=PASS` certifies successful execution
of the declared numerical checks; it is not a full-PAW or continuum accuracy
certificate. A sampled energy/spin gap alone does not establish global band
connectivity: the graphene pair example records a failed link check in the Fukui method
and does not assign that pair a Chern number.

## Mathematical contract and approximation

Let C contain the selected raw pseudo-wavefunctions, O=C†C their Gram matrix,
and let X solve the generalized projected-spin eigenproblem
ΣX=OXΛ with X†OX=I. The orthonormal spin frame is U=CX. For each Cartesian
direction, the native canonical momentum matrix supplies

```math
D^a_{mn}=\frac{\hbar^2}{m_e}\langle\widetilde u_m|k_a+G_a|\widetilde u_n\rangle,
\qquad
\widetilde{\partial_a C_n}
=\sum_{m\notin P}\widetilde C_m\frac{D^a_{mn}}{E_n-E_m}.
```

The sum uses bands 1 through `sum_band_max` outside the selected group;
by default this is every retained WAVECAR band. No
internal energy denominator is evaluated. In the spin eigenframe define

```math
B_a=(1-UU^\dagger)\widetilde{\partial_a C}\,X,
\qquad
K_a=B_a^\dagger\sigma U+U^\dagger\sigma B_a.
```

The sector trace curvature is

```math
\Omega^s_{ab}=-2\operatorname{Im}\sum_{\alpha\in s}
\left[
(B_a^\dagger B_b)_{\alpha\alpha}
+\sum_{\beta\notin s}
\frac{(K_a)_{\beta\alpha}^*(K_b)_{\beta\alpha}}
{(\lambda_\alpha-\lambda_\beta)^2}
\right].
```

Only **opposite-sign spin eigenvalues** enter the second denominator.
Degeneracy inside either spin sector is allowed. The second contribution
cancels between the two sector traces; their sum equals the parent tangent
curvature. Exact-degenerate unitary basis rotations, source spin-frame
rotation and axis reversal are checked in the regression suite.

The geometric kernel is exact for a supplied correct tangent. The production
canonical tangent is an **approximation**: raw pseudo coefficients are not
a complete orthonormal physical PAW basis, and canonical momentum omits
nonlocal/SOC velocity corrections, PAW augmentation and basis-derivative
effects. The code records `exact_pseudo_projector_derivative=false` and
`geometric_integer_certified=false`. It does not reinterpret a global
whitening of unequal-energy states as an unchanged energy eigenbasis.

Independent analytic Hamiltonian tests compare all three native tensor
components with small geometric sector loops, include a nonzero spin-mixing
contribution, test nonorthogonal frame changes, and compare full-BZ integrals
with Chern numbers computed independently using the Fukui method. This validates the kernel; it does not remove the
material calculation's operator and sampling approximations.

The projected-spin construction follows
[Prodan](https://doi.org/10.1103/PhysRevB.80.125327).
This task is distinct from [conventional spin-current Hall conductivity](SPIN_HALL.md)
and from [PROCAR-weighted character attribution](POSTPROCESSING.md).

## Material convergence examples

- [Bi: mesh, retained intermediate states and source quality](../examples/materials/bi-spin-hall/spin-chern-kubo/convergence/) uses `--sum-bands` on unchanged WAVECAR files and compares the raw integral with sector Chern numbers independently evaluated using the Fukui method.
- [Graphene: resolving a tiny SOC peak](../examples/materials/graphene-spin-chern/kubo/convergence/) adds actual VASP valley points, native curvature, local overlap flux and solver controls. The local integrals do not cover the entire BZ.

The [technical report, Section 3.10](TECHNICAL_REPORT.md#310-follow-up-convergence-tests) interprets these numerical controls together.
