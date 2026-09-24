# VASPBERRY full-connection reference

These are the completed outputs of VASPBERRY's own NumPy solver, using the
supplied VASP-derived Hamiltonian and position matrices. VASPBERRY performs
the Fourier transforms, diagonalization, full occupied-bundle J0 + J1 + J2
curvature, and two-dimensional integration. No postw90 response is used as
an input. Follow [the native workflow](../../NATIVE_WANNIER.md) to reproduce
these calculations.

## Result

The final calculation gives **σxy = 1.0003849781827807 e²/h** at each of the
three calculated chemical potentials inside the gap. Its 27,600 weighted
points cover the full two-dimensional torus: an 80×80 base grid with 9×9
subdivision of cells selected within 0.18 Å⁻¹ of periodic Γ. The response
is not rounded or symmetrized. The independent direct-wavefunction Fukui
result is C = −1, which corresponds to positive sheet σxy in this convention.

| Base grid | Local subdivision | Radius (Å⁻¹) | Points | VASPBERRY σxy/(e²/h) |
|---|---|---:|---:|---:|
| 40×40 | 7×7 | 0.12 | 3,088 | 1.0001934913592478 |
| 60×60 | 7×7 | 0.12 | 6,528 | 1.0000548041196191 |
| 60×60 | 11×11 | 0.12 | 10,920 | 1.0000497725580064 |
| 60×60 | 11×11 | 0.18 | 21,720 | 1.0005919297781898 |
| 80×80 | 9×9 | 0.18 | 27,600 | 1.0003849781827807 |

All controls are retained, including the nonmonotonic wider-domain change.
The final same-domain change is 0.0002069515954091 e²/h, below the predefined
10⁻³ e²/h quadrature criterion. This tests integration of this fixed finite
model, rather than full convergence of the underlying material calculation.

## Files

| File | Contents |
|---|---|
| [bands.npz](bands.npz), [bands.json](bands.json) | All 138 native band energies at 3,601 Γ–M–K–Γ points; path, lattice, units and producer metadata |
| [convergence.csv](convergence.csv) | All five native controls and independently computed postw90 values |
| [metadata.json](metadata.json) | Material, source, comparison and integration scope |
| [Final Hall CSV](cases/n80-r9-R018/hall/conductivity.csv), [DAT](cases/n80-r9-R018/hall/conductivity.dat), [NPZ](cases/n80-r9-R018/hall/conductivity.npz) | Identical common Hall-response fields in e²/h and S |
| [Final Hall metadata](cases/n80-r9-R018/hall/conductivity.json) | Full operator, formula, gap, sampling and integrated J0/J1/J2 contributions |
| [Final curvature](cases/n80-r9-R018/curvature.npz) | Actual integration points, weights, energies and full curvature terms |
| [Final run record](cases/n80-r9-R018/run.json) | Completion, timing and source-integrity records |

Each directory under `cases/` contains the same six original output files:
`run.json`, `curvature.npz`, and `hall/conductivity.{csv,dat,npz,json}`.
They are copied without altering numerical values or producer metadata.
The table's electron count is **87 represented occupied model states**,
not the full physical electron count including the omitted deep bundle.

## Independent comparison and normalization

The unchanged [postw90 data](../wannier/README.md) are a separate reference.
All five calculations use identical points and weights. The native totals
and all J0/J1/J2 terms in each of the three Cartesian components agree with
that reference within its printed precision. The final published postw90
sheet value is 1.0003849765685957 e²/h; its raw difference from VASPBERRY is
1.61×10⁻⁹ e²/h.

VASPBERRY integrates the dimensionless sheet response directly using the
reciprocal area and converts S with exact modern SI constants. The original
postw90 binary uses a historical prefactor for its S/cm output, whereas the
published comparison converts that output with modern SI constants. The
ratio of the historical to modern prefactors is 0.999999996332018. This
small offset and the original printed output precision are included in the
comparison; both original datasets remain unchanged.

## Scope

This model has 138 spinor functions and 87 occupied model bands. The omitted
original VASP bands 1–36 have a separately checked C = 0; their local
curvature need not vanish. The structure is fixed and unrelaxed, the source
density and Wannier training meshes are coarse, and localization stopped
after 300 iterations before its requested spread tolerance. The parent
[material guide](../../README.md) gives the independent DFT band/curvature
checks and their limitations. A separate 160×160 energy scan finds a sampled
17.1014 meV model gap containing all three chemical potentials.

The response is a zero-temperature occupied-bundle result. Its constant
value inside the gap follows from unchanged occupations; quantization is
assessed from its value, the sampling controls and the independent Chern
invariant. No finite-temperature or metallic calculation is inferred from
these three points.
