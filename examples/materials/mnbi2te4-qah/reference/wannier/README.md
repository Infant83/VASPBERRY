# Full-connection Wannier reference data

These files describe the same fixed first-principles model as the supplied
[Hamiltonian and position matrices](../../inputs/wannier/operators/README.md).
The external producer is Wannier90/postw90 3.1.0. Its effective-model reader
uses the documented initialization fix; the physical full-connection formulas
are unchanged. The standard VASP/WAVEDER result is a separate
[coarse-grid diagnostic](../waveder-coarse/README.md).

| File | Meaning |
|---|---|
| `bands.csv` | Actual Wannier interpolation along Γ–M–K–Γ; all model bands entering the plotted energy window |
| `outputs/wannier90_band.dat.gz` | Complete original output for all 138 model bands |
| `local-validation.csv` | Band edges and full curvature at 12 coordinates, compared with actual VASP optical calculations |
| `conductivity.csv`, `.dat`, `.npz` | Identical final sheet-Hall data at three calculated chemical potentials, with raw Cartesian S/cm components retained |
| `convergence.csv` | All prescribed integration controls; `plot=true` selects the true-distance sequence shown in the figure |
| `quadrature.npz` | Final fractional coordinates, absolute integration weights, parent cells and lattice |
| `sampled-model-gap.json` | Energy-only check of the same exact model on a 160×160 mesh |
| `outputs/part*-ahc-fermiscan.dat`, `part*.wpout` | Original final partition responses and full J0/J1/J2 producer output |
| `metadata.json` | Model, numerical conventions, selected quadrature, limits and source records |

`mu_eV` retains the VASP energy zero. The display zero is the recorded DFT
midgap; no band alignment is fitted. postw90's x/y/z conductivity components
are σyz/σzx/σxy. The CSV reports the z component as `sigma_xy_S_per_cm` and
converts it to `sigma_xy_e2_over_h` using the complete simulation-cell height,
including the vacuum. No extra spin, surface or film-thickness factor is used.

Each parent cell of the two-dimensional torus retains its area. Cells whose
centres lie within the specified periodic Cartesian distance of Γ receive
an odd square submesh. All chemical potentials use the same quadrature, and
partition sums are added without separate normalization. `convergence.csv`
retains the narrower and wider refinement domains even when their deviation
from an integer is nonmonotonic. The final case is the prescribed wider,
finer calculation; it is not chosen for being numerically closest to one.

The three chemical potentials span the central 90% of the sampled DFT gap.
Their symbols in the figure are actual calculations. The connecting line is
a visual guide; it is not a densely sampled metallic Hall curve. At T = 0,
a response is constant in a gap because the occupied subspace is unchanged.
The unrounded value and quadrature comparisons establish the numerical
quantization check; flatness alone does not.

The basis contains 138 spinor Wannier functions with 87 occupied bands.
Original DFT bands 1–36 were excluded after a separate C = 0 check. Their
local curvature need not vanish, so the difference between the retained
model and the 123-band DFT occupied trace is not purely interpolation error.
Localization ended after 300 iterations without reaching its spread
threshold. The geometry, density, DFT training mesh and finite basis are
held fixed here. Convergence of this model's integral does not establish
convergence of the material calculation.
