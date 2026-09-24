# Spin matrices, spin Hall response and quantum spin Hall topology

VASPBERRY evaluates the conventional intrinsic spin Hall conductivity of a
gapped two-dimensional occupied band group at zero temperature. It uses full
complex spin and velocity matrices from the same VASP eigenstates. The
separate Fukui–Hatsugai Z₂ calculation identifies the time-reversal topological
phase. These are different observables: a nontrivial Z₂ phase does not require
an exactly quantized conventional bulk spin Hall conductivity when SOC mixes
spin components.

The [actual Bi bilayer tutorial](../examples/materials/bi-spin-hall/) provides
VASP inputs, a fresh charge density, physical PAW matrix bundles, native
conductivity outputs, a Z₂ n-field, bands and ideal-edge figures. Its measured
sampling and interpolation errors are reported explicitly.

## What is implemented

| Command | Input and calculation |
|---|---|
| `spin-export` | Completed instrumented VASP run → validated PAW spin/full-velocity arrays |
| `spin-merge` | Same-density VASP k chunks → one checked full mesh and unchanged per-k wavefunctions |
| `spin-matrix` | SOC WAVECAR coefficients → all complex Pauli matrix elements and raw overlap; optional matching PAW augmentation |
| `spin-hall` | Audited full PAW spin/velocity export → occupied-bundle spin-curvature tensor and intrinsic 2D sheet response |
| `wannier-edge` | VASP-derived Wannier Hamiltonian → ideal finite-strip spectrum and boundary probabilities |
| Native Fukui–Hatsugai Z₂ | Occupied VASP spinor wavefunctions → integer n-field and half-zone parity |

The spin Hall command supports a full uniform 2D mesh, all occupied bands
starting at band 1, and a chemical potential strictly inside the sampled
global gap. It does not implement metallic occupations, finite temperature,
disorder, torque-dipole corrections, layer-current operators or 3D bulk spin
conductivity. Charge Hall and spin Hall have distinct current operators.
Weighting charge Berry curvature by a band's spin expectation is generally
insufficient in SOC systems.

## Spin from WAVECAR and PAW

For two-component coefficients $c_{n\mathbf Gs}$, the raw matrix is

$$
\widetilde\Sigma^a_{nm}
=\sum_{\mathbf Gss'}c^*_{n\mathbf Gs}(\sigma_a)_{ss'}c_{m\mathbf Gs'},
\qquad s_a=\frac{\hbar}{2}\sigma_a.
$$

Both spinor components and every off-diagonal band element are retained.
The coefficient overlap $\widetilde O_{nm}=\sum_{\mathbf Gs}
c^*_{n\mathbf Gs}c_{m\mathbf Gs}$ is also saved unchanged. WAVECAR contains
smooth pseudo-wavefunctions; it does not by itself contain the PAW on-site
augmentation. Independently normalizing its bands would change the operator
and does not supply that augmentation.

For the supported spin-independent PAW transformation, the physical matrix is

$$
\Sigma^a_{nm}=\widetilde\Sigma^a_{nm}
+\sum_{Aijss'}p^{A*}_{ni,s}\,\Delta Q^A_{ij}
(\sigma_a)_{ss'}p^A_{mj,s'},\qquad
\Delta Q^A_{ij}=\langle\phi^A_i|\phi^A_j\rangle
-\langle\widetilde\phi^A_i|\widetilde\phi^A_j\rangle.
$$

Here $p^A_{ni,s}=\langle\widetilde p^A_i|\widetilde u_{ns}\rangle$.
Replacing $\sigma_a$ by the spin identity gives the overlap correction.
This is the matrix representation of $T^\dagger\sigma_aT$, using the
same projectors, wavefunctions and partial waves as the VASP calculation.
The importer checks the corrected overlap against the identity and the
operator bounds $O\pm\Sigma^a\succeq0$. It does not repair failed matrices.

WAVECAR does not encode the spinor-axis convention. Obtain it from the
matching SAXIS and OUTCAR. `spin-matrix --spin-basis` records that declaration;
it does not rotate the coefficients. The current audited physical producer
uses Cartesian spin axes with `SAXIS = 0 0 1`.

```sh
python3 tools/vaspberry_kubo.py spin-matrix --wavecar WAVECAR \
  --spin-basis "Cartesian axes; SAXIS=(0,0,1)" \
  --augmentation augmentation.npz \
  --augmentation-metadata augmentation.json \
  --output-dir results/spin
```

Omit both augmentation arguments to inspect the raw pseudo matrices. Their
output explicitly records `operator_scope=raw_pseudo`. That output alone is
not accepted as a physical spin Hall input. PAW augmentation must match the
source wavefunction, k points, band order, energies, lattice and spin basis.

## Velocity and the finite-band spin current

Let $D_i=\hbar v_i$, measured in eV Å. The conventional spin current is
$J_i^a=\{s_a,v_i\}/2$. Its numerical representation uses
$K_i^a=\{\Sigma^a,D_i\}/2$ and accounts for $s_a=\hbar\sigma_a/2$
in the final unit conversion.

Forming this product needs the diagonal and exactly degenerate blocks of
the full velocity, in addition to occupied–empty optical matrix elements.
Reconstructing $D$ only by multiplying a stored interband connection by an
energy difference cannot recover blocks discarded at degeneracies. The
audited producer retains the velocity before that division, including the
PAW projector terms. Canonical plane-wave momentum alone is insufficient
for this physical-operator route.

In a finite source-band space $P$, the implemented product is

$$
K^{a,P}_i=\frac12\{P\Sigma^aP,PD_iP\}.
$$

It differs from $P\{\Sigma^a,D_i\}P/2$ by the omitted terms
$[P\Sigma^aQD_iP+PD_iQ\Sigma^aP]/2$, where $Q=1-P$.
Thus source-band convergence tests both the current-product closure and the
intermediate-state sum. The multiplication is performed over all retained
source bands **before** selecting occupied–empty transitions. The numerical
kernel also accepts independently supplied spin-current matrices; the present
VASP command uses the declared finite-band product. This approximation is
explicitly discussed by [Ryoo, Park and Souza](https://doi.org/10.1103/PhysRevB.99.235113).

## Kubo formula, signs and units

For the fixed occupied group $\mathcal V$,

$$
\Omega^{a}_{ij}(\mathbf k)=-2\,\mathrm{Im}
\sum_{n\in\mathcal V,m\notin\mathcal V}
\frac{K^a_{i,nm}D_{j,mn}}{(E_m-E_n)^2}.
$$

Equal-occupation internal terms cancel before division. Internal Kramers
degeneracies are allowed; a touching occupied–empty boundary is rejected.
The spin curvature-like tensor has units Å². With $e>0$ and normalized
full-zone weights, the physical sheet coefficient is

$$
\frac{\sigma^a_{ij}}{(\hbar/e)(e^2/h)}
=\frac{A_{\rm BZ}}{4\pi}\sum_{\mathbf k}w_{\mathbf k}\Omega^a_{ij}(\mathbf k).
$$

For comparison, charge Hall uses
$\sigma^{\rm charge}_{ij}/(e^2/h)=-A_{\rm BZ}\sum w\Omega^{\rm charge}_{ij}/(2\pi)$.
In a conserved $\sigma_z=+1$ sector, the normalized spin coefficient is
therefore minus one half of the normalized charge coefficient. The factor
comes from physical spin $\hbar/2$ and the electron charge $-e$; no additional
spin multiplicity or post hoc integer rounding is applied. The opposite
sectors can cancel in charge Hall while adding in spin Hall. This is the
conventional current response, not a spin Chern number.

## Calculate and reproduce

The [audited VASP producer recipe](../tools/vasp544_spin_bridge/) gives the
original extraction code, build steps and run wrapper. It requires a licensed
copy of the supported VASP source. The producer writes the NPZ arrays and accompanying metadata from VASP;
users do not author model matrices or numerical results in JSON. JSON records
the operator contract and provenance. Once the actual VASP-derived files are
available, a calculation is:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python3 tools/vaspberry_kubo.py spin-hall \
  --matrices physical-matrices.npz --metadata physical-matrices.json \
  --mesh 12 12 --occupied 10 --formats csv dat npz \
  --output-dir results/spin-hall
```

The mesh and occupied count shown are Bi-specific; set them for your own
system. The default chemical potential is the sampled gap midpoint. For a
closure diagnostic, use `--source-band-limit` on the same cache, preserving
the complete occupied group and without splitting a degenerate cutoff group.
Also rerun VASP with increased `NBANDS`, since changing the source calculation
and truncating a fixed set of matrices are different checks.
Check near-degenerate groups at the retained cutoff against the source
eigensolver accuracy, beyond the numerical zero-gap threshold. The output
records the minimum cutoff gap; numerical splitting of a Kramers pair does
not justify keeping only one partner.

`--k-chunk` bounds the temporary spin-current matrices. The default memory
estimate limit is 2048 MiB; `--time-limit` is checked between chunks. These
are preflight and cooperative checks, not operating-system resource caps.
The command is serial NumPy; it does not acquire MPI support from the native
Fortran executable. Independent VASP mesh chunks can be calculated separately. `spin-merge`
checks their full mesh, common fixed density, potentials, physical settings
and per-k gauge association before assembly:

```sh
python3 tools/vaspberry_kubo.py spin-merge \
  --run-dirs vasp-chunk00 vasp-chunk01 --mesh 12 12 --occupied 10 \
  --output-dir results/assembled
```

The command writes the physical matrix bundle and an assembled WAVECAR,
checking that every per-k record is byte-identical to its source. It preserves
the first chunk Fermi header, which is not a global chemical potential;
`spin-hall` selects its chemical potential from the assembled gap.

Outputs include every Cartesian spin/current/field component, local
curvature, weights, the charge-Hall comparison, gaps, units and operator
scope. CSV, whitespace DAT and compressed NPZ are independently selectable.
Metadata is always JSON. NPZ is read with `allow_pickle=False`. Existing
outputs are not overwritten; failed integrations retain `.partial/run.json`.
See [the format reference](OUTPUT_FORMAT.md).

## Edge spectrum and interpretation

```sh
python3 tools/vaspberry_kubo.py wannier-edge --operators results/operators \
  --width 40 --edge-cells 2 --periodic-axis 0 --open-axis 1 \
  --kpoints 201 --formats csv npz --output-dir results/edge
```

The strip retains hopping inside the specified number of Wannier unit cells
and removes hopping across its open edges. It is an ideal termination of the
supplied model, without surface relaxation or self-consistency. Check the
bulk interpolation against direct VASP bands and vary strip width before
interpreting an edge gap. Boundary weight summed over a degenerate group,
or a broadened edge spectral function, is invariant under rotations inside
that group; individual degenerate eigenvectors need not be edge-localized.

A gapped bulk with Z₂=1 and boundary branches connecting valence and
conduction states supports the QSH interpretation. It does not by itself
calculate a finite-device conductance. The Bi motivation follows
[Murakami](https://doi.org/10.1103/PhysRevLett.97.236805).

## References

- [Ryoo, Park and Souza, Phys. Rev. B 99, 235113 (2019)](https://doi.org/10.1103/PhysRevB.99.235113).
- [Qiao et al., Phys. Rev. B 98, 214402 (2018)](https://doi.org/10.1103/PhysRevB.98.214402).
- [Murakami, Phys. Rev. Lett. 97, 236805 (2006)](https://doi.org/10.1103/PhysRevLett.97.236805).
- [Kane and Mele, Phys. Rev. Lett. 95, 146802 (2005)](https://doi.org/10.1103/PhysRevLett.95.146802).
