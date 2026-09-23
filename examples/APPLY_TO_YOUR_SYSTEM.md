# Apply a VASPBERRY tutorial to your material

First reproduce the named tutorial with its supplied VASP outputs. Compare
native output files, numerical checks and the figure with `reference/`.
Then use the same **production command** with your own files and settings.
The wrappers provide fixed tutorial settings and compare the named reference
only for matching inputs. Use each production CLI for general calculations;
see individual READMEs for supported wrapper overrides.

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
required by a Brillouin-zone integral. The example's numerical `-kx`, `-ky`,
`-ii`, `-if`, `-k` and chemical potentials belong to its material.

## Replace the relevant settings

1. Put your VASP outputs in a stable input directory. Check WAVECAR format and
   record-length compatibility using the [build guide](../docs/BUILD.md).
2. Read your EIGENVAL/OUTCAR and choose occupied bands, selected states, spinor
   representation and the actual k sampling. SOC normally uses `-s 2`; this
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
starts from a full VASP WAVECAR and runs native pair export followed by the
bundled occupation and integration routines in one command. It includes
chemical-potential and temperature scans, regional contributions and
independent mesh and band-window checks. Reuse its general command with your
material's band filling, energy reference and region definitions.

The native WAVECAR Kubo implementation
uses canonical momentum; full material velocity can require PAW, nonlocal,
SOC or other corrections. The general matrix interface accepts an explicitly
declared physical operator. It is a separate data-export route and does not
turn a JSON configuration into a VASP calculation.

The [output specification](../docs/OUTPUT_FORMAT.md) defines units and formats.
Use CSV, text DAT or NumPy NPZ tables in your own analysis; preserve normalization,
provenance, excluded points and region definitions. See [migration](../docs/MIGRATION.md)
before combining results with older doubled Kubo files.
