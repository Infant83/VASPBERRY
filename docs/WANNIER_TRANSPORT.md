# Full-connection Wannier transport in VASPBERRY

`wannier-import`, `wannier-bands` and `wannier-hall` provide a second input
route for VASP-derived electronic structure. VASPBERRY performs the Fourier
sums, diagonalization, Berry-curvature evaluation and Hall integration itself.
Wannier90 prepares the localized electronic model; its `postw90` response is
used only as an independent comparison in the material example.

The [MnBi₂Te₄ tutorial](../examples/materials/mnbi2te4-qah/NATIVE_WANNIER.md)
provides actual VASP-derived operators and complete reproduction commands.

## Choose the input route

| Input | Operator and scope |
|---|---|
| `WAVECAR`, through `wavecar-hall` | Canonical-momentum approximation; chemical-potential and temperature scans |
| Standard optical VASP 5.4.4 run, through `waveder-hall` | PAW optical occupied-to-empty matrix elements; fixed insulating bundle at T=0 |
| Full Hamiltonian and position matrices, through `wannier-hall` | Full connection of the supplied finite Wannier model; fixed insulating bundle at T=0 |

The full connection needs **both** Hamiltonian and position matrices. A
Hamiltonian-only tight-binding calculation generally omits connection terms
of a VASP-derived Wannier basis. Adding k points to the WAVECAR route does
not supply its omitted PAW/nonlocal/SOC operator terms.

## Import actual operators

Use paired effective-model `seed_HH_R.dat` and `seed_AA_R.dat` files with
the matching VASP lattice. These are the documented Wannier90 effective-model
files, **not** the conventional `seed_hr.dat` Hamiltonian file. All real-space
and pair-translation weights must already be absorbed. The material tutorial
includes the exact full-operator export and restore procedure.

```sh
python3 tools/vaspberry_kubo.py wannier-import \
  --hh seed_HH_R.dat --aa seed_AA_R.dat --poscar POSCAR \
  --spinor-components 2 --spin-multiplicity 1 \
  --energy-reference "Unshifted source eigenvalues in eV" \
  --output-dir results/operators
```

The importer checks dimensions, fixed-width fields, lattice translations,
spin declarations and R/−R Hermiticity. Repeated source records are added;
they may represent double-precision residuals. It does not discard small
matrix elements, repair inconsistent matrices or infer a missing position
operator. The normalized NPZ/JSON cache can be reused for different meshes.
Record the source model's energy window and construction alongside it.

## Compute bands and Hall response

Define any continuous path by fractional-coordinate triples:

```sh
python3 tools/vaspberry_kubo.py wannier-bands --operators results/operators \
  --vertices 0 0 0 0.5 0 0 0.3333333333333333 0.3333333333333333 0 0 0 0 \
  --labels Gamma M K Gamma --points-per-segment 201 \
  --formats csv npz --output-dir results/bands
```

The Hall command requires the occupied model-band count and a chemical-potential
interval strictly inside the sampled global gap. Set these from your own model;
the following numerical values belong only to the supplied MnBi₂Te₄ example:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python3 tools/vaspberry_kubo.py wannier-hall --operators results/operators \
  --occupied 87 --mesh 80 80 \
  --refine 9 --refine-radius 0.18 --refine-center 0 0 \
  --mu-min 3.254290679571569 --mu-max 3.269681923569488 --mu-num 3 \
  --mu-reference 3.2619863015705284 \
  --workers 4 --batch-size 16 --formats csv dat npz \
  --output-dir results/hall
```

Without refinement options, the command uses the complete uniform mesh.
Refinement replaces selected coarse cells with odd centered submeshes;
cell areas supply the weights and the full periodic domain is preserved.
`--refine-center` can be repeated. Its two coordinates refer to the ordered
`--plane-axes` (default 0 1); the excluded fractional coordinate is zero.
For sheet units, the remaining real-cell vector must be perpendicular to
the selected plane. This is a two-dimensional model response, not a bulk
three-dimensional conductivity.

`--workers` distributes NumPy batches across threads. Set BLAS threads to one
when using several workers to avoid nested parallelism. The default
`--memory-limit-mib 4096` bounds a preflight memory estimate; it is not an OS
memory cap. `--time-limit 3600` is checked between batches. Failed calculations
retain a `.partial` directory with their status; completed output is renamed
only after source integrity and all requested gap conditions pass.

## Formula and outputs

With $A=i\langle u|\partial_k u\rangle$, let $o$ and $e$ denote occupied and
empty eigenstates, $D_a^{oe}=\langle o|\partial_aH|e\rangle/(E_e-E_o)$, and
$\bar A_a^{oe}=\langle o|A_a^{W}|e\rangle$. Let $U_o$ contain the occupied
Hamiltonian eigenvectors in the Wannier basis. The occupied trace is

$$
\begin{aligned}
\Omega_{ab}&=J0+J1+J2,\\
J0&=\mathrm{Re}\,\mathrm{Tr}[U_o^\dagger(\partial_a A_b^W-\partial_b A_a^W)U_o],\\
J1&=2\mathrm{Re}\sum_{oe}(\bar A_a^{oe}D_b^{oe*}-D_a^{oe}\bar A_b^{oe*}),\\
J2&=-2\mathrm{Im}\sum_{oe}D_a^{oe}D_b^{oe*}.
\end{aligned}
$$

Only occupied-to-empty denominators occur, so degeneracies internal to either
subspace are allowed.
The Fourier phase is $e^{2\pi i\mathbf q\cdot\mathbf R}$, while Cartesian
derivatives use $iR_a$ with real-space vectors in Å. The sheet response is
$\sigma/(e^2/h)=-A_{BZ}\sum_k w_k\Omega_{\rm normal}/(2\pi)$, including the
explicit spin multiplicity. No integer rounding or extra half factor is used.

- `hall/conductivity.csv`, `.dat`, `.npz`: common Hall table format.
- `hall/conductivity.json`: model, gap, sampling and operator conditions;
  integrated J0/J1/J2 for all three axial components `(yz,zx,xy)`.
- `curvature.npz`: actual k points, cell weights, parent cells, all local
  J0/J1/J2 components, total curvature and valence/conduction edge energies.
- `run.json`: completion state, timing and source integrity.
- Band output: `bands.csv` and/or `bands.npz`, plus units and path information
  in `bands.json`. Energies retain the source zero.

This route currently supports insulating T=0 bundles. A flat gap scan alone
does not test quadrature. Check mesh, refined-region size and model construction
separately. Full connection refers to the supplied finite subspace; it does
not certify that the source DFT or Wannier model is converged.

## References

- [Wang et al., Phys. Rev. B 74, 195118 (2006)](https://doi.org/10.1103/PhysRevB.74.195118).
- [Lopez et al., Phys. Rev. B 85, 014435 (2012), Eq. (51)](https://doi.org/10.1103/PhysRevB.85.014435).
