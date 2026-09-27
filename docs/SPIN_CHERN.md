# Projected-spin Chern numbers from WAVECAR

The native `spin-chern` task calculates both projected-spin sectors of an
isolated spinor subspace. The main workflow selects the complete occupied
space. This native task is available in VASPBERRY 1.6.1.

## Execute VASPBERRY

Build the executable with [Intel Fortran and Intel MPI](BUILD.md#intel-oneapi-on-linux):

```bash
source /opt/intel/oneapi/setvars.sh
make ifx-mpi
vb_bin="$PWD/build/vaspberry-ifx-mpi"
```

In the directory containing a matching full-mesh `WAVECAR` and `OUTCAR`, run:

```bash
mpiexec -n 4 "$vb_bin" --task spin-chern \
  --mesh 12,12 --bands 1:8 --spin-axis z
```

Here `12,12` and `1:8` are example choices: use your source mesh and complete
occupied band range. The calculation runs inside VASPBERRY Fortran/MPI;
Python is optional for plotting its CSV results. Use a fresh working
directory for each calculation. `--spin-axis` accepts `x`, `y`, `z` or a
Cartesian unit vector such as `0,0,1`; the default is `z`.

`WAVECAR` is the default input in the working directory. It supplies the
wavefunctions, eigenvalues, k points and lattice, and VASPBERRY detects its
two-component layout automatically. `OUTCAR` supplies the matching spin-frame
transformation from SAXIS to Cartesian axes; it does not supply wavefunctions
or determine the spinor count. An absent or ambiguous spin frame is an error.
For files elsewhere, add `--wavecar /path/to/WAVECAR --outcar /path/to/OUTCAR`.
Both defaults are relative to the working directory.

`--bands` can also select another isolated subspace. Its reported charge sum
then belongs to that selected subspace, not the full occupied system. For an
insulator, use all occupied states and inspect the common global energy gap.
For a `1:N` selection, the CSV reports `has_global_energy_gap` separately from
`source_occupations_match_integer_projector`. Finite VASP smearing can make
the occupation diagnostic false without invalidating the selected subspace's
geometric calculation. This task requires `ISPIN=1` spinors and at least one
excluded source band.

## What is calculated

Within the complete occupied space, VASPBERRY forms the full matrix of
spin along the chosen Cartesian axis, including off-diagonal band elements.
Its positive and negative eigenspaces define the two sectors. The
Fukui–Hatsugai–Suzuki (FHS) method uses determinant link variables from
wavefunction overlaps to calculate their Chern numbers:

```math
C_{\mathrm{charge}}=C_++C_-,\qquad
\Delta C=C_+-C_-,\qquad C_{\mathrm{spin}}=\frac{\Delta C}{2}.
```

The files retain the raw numerical sums. For a general system, the declared
half-difference need not be an integer.

These sectors can contain internally degenerate states. The occupied space
must remain separated from excluded states, and the projected-spin spectrum
must stay separated from zero. Selecting individual bands by the sign of
their spin expectation does not construct these subspaces. Reversing the
analysis axis exchanges the sector labels.

The calculation uses the WAVECAR pseudo-wavefunction Gram matrix to define
an orthonormal occupied frame. This is geometric topology in the stated
pseudo representation; it does not add PAW augmentation. Its numerical
checks and mesh refinement must accompany interpretation. The projected-spin
construction follows [Prodan](https://doi.org/10.1103/PhysRevB.80.125327).

This result is distinct from [PROCAR-weighted charge attribution](POSTPROCESSING.md#add-layer-and-spin-character)
and from [conventional spin Hall conductivity](SPIN_HALL.md). The latter
uses a spin-current operator and has conductivity units. The two spinor
components of an SOC state are not separate VASP `ISPIN=2` channels.

To draw sector-resolved Kubo-formula Berry curvature along a path, use the separate native
[`spin-kubo` task](SPIN_KUBO.md). Its canonical-momentum derivative proxy
includes the changing spin-sector basis. Compare its raw full-mesh integral
with this task's geometric result; a path alone cannot determine a Chern number.

## Read and plot the results

| File | Numerical data | Useful plot |
|---|---|---|
| `SPIN_CHERN.csv` | `C_plus`, `C_minus`, `C_charge`, `delta_C`, `C_spin`, independently evaluated `C_parent` and sector ranks | Summary across meshes or material parameters |
| `SPIN_BERRY.csv` | Sector/parent fluxes in radians, plaquette centers, area and sector curvature in Å² | Positive/negative sector Berry-flux or curvature maps |
| `SPIN_SPECTRUM.csv` | Projected-spin eigenvalues at each sampled k point | Spin-spectrum and spin-gap maps |

Skip comment lines beginning `#` when importing the CSVs into Origin, a
spreadsheet or a Python plotting program. Keep the comments with the data:
they identify the representation, axis and numerical checks. See the
[output reference](OUTPUT_FORMAT.md#native-projected-spin-chern) for units.
Existing output files are rejected. `--output NAME` changes the prefix from
`SPIN` to `NAME`, producing `NAME_CHERN.csv`, `NAME_BERRY.csv` and
`NAME_SPECTRUM.csv`.

An energy gap is measured in eV. The projected-Pauli distance from zero is
dimensionless; it is not an energy broadening. A closed spin gap means this
axis-dependent spin-sector decomposition is unresolved, even when the
ordinary occupied-space Chern-number or Z₂ calculation remains meaningful.

The default tolerances are `--energy-gap-tol 1e-8` eV for selected-to-excluded
states and `--spin-gap-tol 1e-6` for the projected-Pauli distance from zero.
These numerical cutoffs do not establish VASP or mesh convergence.

| Error | What to check |
|---|---|
| Parent subspace is not energy isolated | Include the complete touching band group; inspect the VASP states and selected range. |
| Projected-spin zero gap is closed / sector rank changes | The chosen axis does not define resolved sectors on this mesh. Inspect the spin spectrum and sampling. |
| Singular link / plaquette branch / sector-sum mismatch | Refine and verify the periodic mesh; the selected-space invariant is unresolved. |
| Missing or inconsistent OUTCAR spin frame | Supply the matching OUTCAR with its printed SAXIS-to-Cartesian rotation, dimensions and lattice. |

The Gram metric and link singular values are also recorded. Ill-conditioned
input is rejected rather than assigned an invariant.
