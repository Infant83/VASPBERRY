# Apply the native WAVECAR workflow to your material

Start with a real-material tutorial and compare its native outputs with the
supplied reference. Then use the **same native Fortran command** with your
VASP files, band selection and k sampling. Fukui, Z₂, Kubo curvature, optical
transitions and wavefunction export read WAVECAR directly; Wannierization is
not required for these calculations.

The calculation produces text or CSV tables and, for wavefunctions, volumetric
amplitude grids. Use the supplied Python tools for plots and numerical
postprocessing. Hall transport adds occupation weighting and BZ integration
after native pair export. The optional tutorial `run.py` helpers combine
these steps with fixed reference checks; use the native CLI for general work.

## Select the input and observable

| Goal | Input sampling and decisions |
|---|---|
| Fukui curvature / Chern | A full periodic 2D mesh and an isolated band or isolated fixed-rank band bundle |
| Z₂ | A gapped nonmagnetic TR-symmetric spinor calculation; even, unshifted Gamma-centered full Nx×Ny×1 mesh with ISYM=-1 |
| Kubo curvature along a path | WAVECAR for the desired points; choose bands using your EIGENVAL and inspect near-degeneracies |
| Charge Hall | A full integration mesh, a valid occupied subspace or valid point curvature, an energy reference and a sufficient band window |
| Optical response | Initial/final bands and photon-energy range appropriate to your system; check transition strength before forming a ratio |
| Wavefunction | Matching WAVECAR/POSCAR/EIGENVAL, actual Gamma-point index, band and real-space grid |

A symmetry-reduced mesh or band path cannot replace the full periodic mesh
required by a Brillouin-zone integral. The example's numerical `--mesh`, `--bands`,
`--kpoint` and chemical potentials belong to its material.

## Replace the relevant settings

1. Put your VASP outputs in a stable input directory. Check WAVECAR format and
   record-length compatibility using the [build guide](../docs/BUILD.md).
2. Read your EIGENVAL/OUTCAR and choose occupied bands, selected states, spinor
   representation and the actual k sampling. SOC normally uses `--spinor 2`; this
   means two spinor components, not an extra factor of two in conductivity.
3. Copy the tutorial's explicit VASPBERRY command into a fresh result directory.
   Replace the WAVECAR path and all material-dependent indices/grid settings.
4. Inspect result status and diagnostics before plotting. Near-degenerate
   individual bands do not have a reliable separately resolved curvature;
   use a valid isolated subspace or a suitable formulation. A finite optical
   ratio from vanishing transition intensity is not reliable selectivity.
5. Recompute at denser k meshes and appropriate band/basis settings. Agreement
   with a reference tutorial verifies usage; it does not establish convergence
   for another material.

## Choose the correct transport route

The [Bi Hall tutorial](features/hall-valley/) uses occupied-subspace Fukui flux
in the insulating gap. It demonstrates a time-reversal-symmetric zero charge
Hall result. It is not a point-Kubo calculation, a valley-Hall effect, or a
nonzero Chern-insulator benchmark.

For point-curvature transport, the [Kubo/Hall guide](../docs/KUBO_TRANSPORT.md)
explains normalized curvature, occupations, intermediate-band windows and
user-defined reciprocal-space regions. The [MoS₂ Hall tutorial](features/kubo-hall/)
starts from a full VASP WAVECAR and shows native pair export, reusable pair
import, occupation-weighted integration and plotting as separate commands. It includes
chemical-potential and temperature scans, regional contributions and
independent mesh and band-window checks. Reuse those commands with your
material's band filling, energy reference and region definitions. Python
`pair-hall` performs the integration; `plot_hall.py` reads the finished table.

The native WAVECAR Kubo implementation
uses canonical momentum; full material velocity can require PAW, nonlocal,
SOC or other corrections. The general matrix interface accepts an explicitly
declared physical operator. It is a separate data-export route and does not
turn a JSON configuration into a VASP calculation.

The [output specification](../docs/OUTPUT_FORMAT.md) defines units and formats.
Use CSV, text DAT or NumPy NPZ tables in your own analysis; preserve normalization,
provenance, excluded points and region definitions. See [migration](../docs/MIGRATION.md)
before combining results with older doubled Kubo files.

## Optional extension: Wannier interpolation and edge spectra

For a validated VASP-derived Wannier representation, import **both** real-space
Hamiltonian and position operators with `wannier-import`, compute bands with
`wannier-bands`, and integrate an insulating occupied bundle with `wannier-hall`.
These calculations run inside VASPBERRY. The [general guide](../docs/WANNIER_TRANSPORT.md)
defines the input format, Cartesian conventions, full-zone refinement and
resource controls; the [MnBi₂Te₄ example](materials/mnbi2te4-qah/NATIVE_WANNIER.md)
provides actual inputs and reference results.

Replace the lattice, operator pair, spin convention, occupied count and energy
zero with your own system's values. Validate the Wannier bands and local
curvature against direct calculations, then converge the BZ integral. This
backend currently requires T=0 and all chemical potentials inside a sampled
global gap; use the documented μ/T workflows for their supported operators
when studying metallic occupations.

## Optional extension: physical spin Hall and quantum spin Hall topology

The [Bi spin Hall example](materials/bi-spin-hall/) starts from a fresh SCF
density, matching PAW inputs and the [serial VASP producer](../tools/vasp544_spin_bridge/).
For another material, regenerate the density and operators with its own
structure, potentials and converged electronic settings. Run `spin-export`
on each completed producer calculation, or `spin-merge` for same-density
chunks covering the full mesh, then use the general `spin-hall` command.
See [the spin Hall guide](../docs/SPIN_HALL.md) for the supported VASP revision,
matrix contract and units.

The current response is the conventional intrinsic spin current of a gapped
2D occupied group at T=0. Raw spin expectations multiplying charge Berry
curvature do not replace this operator calculation. Converge both the k mesh
and the finite source-band product; distinguish increased VASP `NBANDS` from
a retained-band cutoff on one set of matrices, and keep complete nearly
degenerate groups. Metallic occupations, finite-temperature spin response
and layer-current operators require a different supported implementation.

Calculate Z₂ separately for a time-reversal-symmetric occupied bundle.
A nontrivial Z₂ result does not require an integer conventional spin Hall
conductivity with SOC. `wannier-edge` gives the spectrum of an ideal strip
of the validated Hamiltonian; it does not include edge relaxation or a
finite-device conductance. The historical Bi_Z2 fixture remains a separate
topology tutorial, with different source provenance.
