# VASPBERRY
Berry curvature, Chern number, two-dimensional intrinsic charge Hall transport,
and two-dimensional Z2 calculations. VASPBERRY reads VASP `WAVECAR`
wavefunctions and exported interband matrices. It implements the discrete Brillouin-zone
method of [Fukui, Hatsugai, and Suzuki, J. Phys. Soc. Jpn. 74, 1674
(2005)](https://doi.org/10.1143/JPSJ.74.1674), together with circular
dichroism and real-space wavefunction output.

**Version 1.3.0 — development version.**

[Technical report](docs/TECHNICAL_REPORT.md) ([PDF](docs/TECHNICAL_REPORT.pdf)) · [Examples](examples/README.md) · [Installation](docs/BUILD.md)
See the [1.3.0 notes](docs/releases/v1.3.0.md) and
[migration guide](docs/MIGRATION.md) for the legacy Kubo factor-of-two fix.

# Download Git version
The examples below are on the unreleased development branch:

```bash
git clone --branch feature/general-kubo-1.3 https://github.com/Infant83/VASPBERRY.git
cd VASPBERRY
```

# Compile

GNU serial and OpenMPI builds are available through the `Makefile`:

```bash
make serial
make mpi
```

These targets write `build/vaspberry-gfortran` and `build/vaspberry-mpi`.
Commands below use the serial product; substitute the path if you build an
executable under another name.

Intel `ifx`/`mpiifx` and legacy `ifort` commands, the direct-access record
length requirement, and BLAS/LAPACK ABI constraints are documented in the
[`build guide`](docs/BUILD.md).

WAVECAR format and compiler compatibility are described in the
[build guide](docs/BUILD.md).

# Features

| Quantity / task | Main interface and guide | Runnable example and results |
|---|---|---|
| Fukui Berry curvature across the Brillouin zone | Fortran wavefunction overlaps | [MoS₂ curvature map](examples/features/fukui-berry-curvature/) |
| Chern number of an isolated band or band bundle | Fortran; [Python Fukui workflow](docs/VALLEY_TRANSPORT.md) | [Bi occupied Chern number](examples/features/fukui-chern/) |
| Two-dimensional Z₂ invariant | Fortran `-z2 1`; [Fukui–Hatsugai guide](docs/Z2_FUKUI_HATSUGAI.md) | [Z₂ / Bi](examples/features/z2/) |
| Circular dichroism / optical selectivity | Fortran `-cd`; [usage below](#usage) | [Circular dichroism](examples/features/circular-dichroism/) |
| Real-space wavefunction at Gamma | Fortran `-wf`; [usage below](#usage) | [Wavefunction](examples/features/wavefunction/) |
| WAVECAR Kubo Berry curvature | Fortran canonical-momentum implementation; [Kubo guide](docs/KUBO_TRANSPORT.md) | [Actual MoS₂ band-path calculation](examples/features/kubo-curvature/) |
| Exported interband-matrix curvature and transport | [Matrix interface and physical-operator contract](docs/KUBO_TRANSPORT.md) | [Developer numerical checks](validation/models/kubo-curvature/) |
| Two-dimensional intrinsic charge Hall response and reciprocal-space regions | Standardized curvature and occupations; [Hall guide](docs/KUBO_TRANSPORT.md) | [Actual Bi occupied-subspace Hall](examples/features/hall-valley/) |
| WAVECAR-direct Fukui export and guarded geometric transport | `tools/wavecar_fukui.py`; [valley-transport guide](docs/VALLEY_TRANSPORT.md) | [Actual Bi occupied-subspace calculation](examples/features/hall-valley/) |
| Historical Kubo normalization import | `tools/vaspberry_kubo.py import-legacy`; [migration guide](docs/MIGRATION.md) | [Import commands and identification checks](docs/MIGRATION.md#historical-doubled-output) |

These are general-purpose interfaces. Choose the observable and input operator
for the physical system; use the standardized output in your own analysis.
The [output specification](docs/OUTPUT_FORMAT.md) distinguishes point
curvature in Å² from geometric plaquette flux in radians. A finite-mesh Kubo
integral need not be an integer; a Fukui lattice integer still needs mesh and
band-isolation checks.

# Usage

## Start with actual VASP files

The [user tutorials](examples/README.md) reproduce calculations from public
**MoS₂ and Bi WAVECAR files**. Each explains the input files, the explicit
VASPBERRY command and its options, original output files, plotting, comparison
with a checked reference, and how to apply the command to your own material.

```bash
make serial
python3 -m pip install -r requirements-transport.txt
# A complete real MoS2 example; its WAVECAR is already in the repository:
python3 examples/features/kubo-curvature/run.py --output-dir results/mos2-kubo

# Fetch the real Bi WAVECAR if you do not have its Git LFS payload:
python3 examples/fetch_inputs.py bi --output-dir results/inputs/bi
python3 examples/features/z2/run.py --wavecar results/inputs/bi/WAVECAR \
  --output-dir results/bi-z2
```

The wrappers execute and record the production commands shown in each tutorial.
The [input catalog](examples/INPUTS.md) lists the available VASP files and their k sampling.
For another material, follow [the transfer guide](examples/APPLY_TO_YOUR_SYSTEM.md)
and set the mesh, bands, spinor setting and energies from your own VASP output.

## Apply the calculation to your input

Run `./build/vaspberry-gfortran -h`,
`python3 tools/wavecar_fukui.py --help`, or
`python3 tools/vaspberry_kubo.py --help` for all options.
The following Fortran commands assume a compatible `WAVECAR` in the working
directory. Use a separate output directory for each calculation.

```bash
# Berry flux / Chern: full 12 × 12 mesh, bands 1–18.
./build/vaspberry-gfortran -kx 12 -ky 12 -ii 1 -if 18

# Z₂: full even Gamma-centered SOC mesh, occupied bands 1–10.
./build/vaspberry-gfortran -f WAVECAR -o NFIELD -z2 1 \
  -kx 12 -ky 12 -s 2 -ii 1 -if 10

# Circular optical response for bands 11 → 12.
./build/vaspberry-gfortran -kx 12 -ky 12 -cd 1 -ii 11 -if 12

# Real-space band 18 at Gamma; matching POSCAR and EIGENVAL are also required.
./build/vaspberry-gfortran -wf 18 -k 1 -ng 40,40,40
```

Z₂ requires a nonmagnetic time-reversal-symmetric insulator and a full,
unshifted, even `Nx × Ny × 1` mesh with `Nx,Ny >= 4`, generated with `ISYM=-1`.
Use only a final `Z2_FIELD.csv` with `result_status=PASS`,
`reportable_invariant=1` and matching half-zone parities. Input gap, physical
time-reversal symmetry and mesh convergence require separate checks; see the
[Z₂ guide](docs/Z2_FUKUI_HATSUGAI.md).

For a semimetal, `-ne 18` can select a fixed 18-band geometric manifold when
the automatically counted occupations differ across k points. This does not
supply the changing occupations for metallic Hall transport. Use the
transport workflow with its occupation and band-window checks.

### Python workflows

The [MoS₂ Kubo tutorial](examples/features/kubo-curvature/) starts from a real
WAVECAR and exports canonical-momentum curvature along a band path. The
[Bi Hall tutorial](examples/features/hall-valley/) uses full-mesh occupied-subspace
Fukui flux in its insulating gap. Its unresolved Kramers partners are not
used as separately isolated Kubo bands.

For exported physical interband matrices, use the [matrix/Hall guide](docs/KUBO_TRANSPORT.md)
and [NPZ/JSON output specification](docs/OUTPUT_FORMAT.md). This interface is
separate from the native WAVECAR approximation.

The WAVECAR-direct workflow exports Fukui plaquette flux, four vertex energies,
adjacent-band gaps and link-quality diagnostics. A basic full-mesh export is:

```bash
python3 tools/wavecar_fukui.py WAVECAR --nx 12 --ny 12 \
  --energy-band 19 --output-dir results/direct-fukui
```

The [valley-transport guide](docs/VALLEY_TRANSPORT.md) gives
`--transport-t0`, `--transport-full-t0`, region definitions and the required
unoccupied sentinel band. These scans report cumulative `sigma_xy(mu)`;
a smooth energy derivative or material convergence needs further analysis.
Read the [migration guide](docs/MIGRATION.md) before comparing historical
Kubo data: version 1.3.0 corrects their factor-of-two normalization.

# Examples and material datasets

See the [technical report](docs/TECHNICAL_REPORT.md) for reciprocal-space maps,
band-resolved curvature, optical spectra and their physical interpretation.

- [Feature tutorials](examples/README.md): actual VASP files, production
  commands, numerical reference outputs, figures and application to your system.
- [Material catalog](examples/materials/): existing MoS₂/Bi data, sampling and
  reproduction requirements. Their original paths remain available.
- [1H-MoS₂](examples/1H-MoS2/): full-mesh stored maps and a separate band-path
  WAVECAR. A line-mode input cannot supply a full-BZ integral.
- [Bi buckled honeycomb](examples/Bi_Z2/): reviewed input templates,
  Z₂ result and full-mesh wavefunctions.
- [Standalone Hall CSV plotting](examples/kubo/): reusable plotting for
  standardized transport output.
- [Developer numerical checks](validation/models/): analytic/synthetic fixtures
  retained separately for testing formulas and regressions.

Method details, diagnostics and output contracts are kept in the linked
guides. An importable Python-library port remains a [roadmap](docs/ROADMAP.md)
item.

# Contributors
* Hyun-Jung Kim: Main developer
* Sun-Woo Kim: Circular dichroism and Kubo formula

# Citation of the code:
@software{Kim_VASPBERRY_2018,author = {Kim, Hyun-Jung},doi = {10.5281/zenodo.1402593},month = {8},title = {{VASPBERRY}},url = {https://github.com/Infant83/VASPBERRY},version = {1.0},year = {2018}}

@article{PhysRevLett.128.046401,
  title = {Circular Dichroism of Emergent Chiral Stacking Orders in Quasi-One-Dimensional Charge Density Waves},
  author = {Kim, Sun-Woo and Kim, Hyun-Jung and Cheon, Sangmo and Kim, Tae-Hwan},
  journal = {Phys. Rev. Lett.},
  volume = {128},
  issue = {4},
  pages = {046401},
  numpages = {6},
  year = {2022},
  month = {Jan},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevLett.128.046401},
  url = {https://link.aps.org/doi/10.1103/PhysRevLett.128.046401}
}
