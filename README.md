# VASPBERRY

**Topology and response functions directly from VASP wavefunctions.**
VASPBERRY reads `WAVECAR` and evaluates wavefunction overlaps and interband
matrix elements in Fortran. Use Fukui methods for Chern numbers and the 2D Z₂
invariant, and Kubo calculations for Berry-curvature maps, symmetry paths and
intrinsic charge Hall response. **No Wannier construction or interpolation is
required for this workflow.**

```text
VASP calculation → WAVECAR → VASPBERRY Fortran → numerical output → analysis / plots
```

VASP supplies the material's electronic structure. VASPBERRY postprocesses it;
band plots provide context for the calculated topology and response.

[Technical report](docs/TECHNICAL_REPORT.md) ([PDF](docs/TECHNICAL_REPORT.pdf)) ·
[Feature examples](examples/README.md) · [Build guide](docs/BUILD.md) ·
[Output formats](docs/OUTPUT_FORMAT.md)

**Latest release: [1.3.0](https://github.com/Infant83/VASPBERRY/releases/tag/v1.3.0).**
The commands below use current `master`, including the new readable Fortran
options. The original short options remain supported. See the
[unreleased changes](CHANGELOG.md#unreleased), [version policy](docs/RELEASING.md)
and [1.3.0 migration notes](docs/MIGRATION.md).

## Download and compile

```bash
git clone https://github.com/Infant83/VASPBERRY.git
cd VASPBERRY
make serial
./build/vaspberry --help
```

The serial executable is `build/vaspberry`; `build/vaspberry-gfortran` remains
a compatibility name. For MPI, run `make mpi` and use `build/vaspberry-mpi`.
See the [build guide](docs/BUILD.md) for GNU/Intel compilers and MPI checks.
Python is needed for the supplied postprocessing and plotting tools:

```bash
python3 -m pip install -r requirements-transport.txt
```

For the fixed 1.3.0 release, clone its tag and follow that version's commands:

```bash
git clone --branch v1.3.0 --single-branch https://github.com/Infant83/VASPBERRY.git VASPBERRY-1.3.0
```

Update an existing `master` checkout with `git pull --ff-only`. Release tags
stay fixed; source archives are available on the release page. Large Git LFS
inputs may need the [input-fetch procedure](examples/INPUTS.md).

## Features

| Task | Native command | Actual VASP example and results |
|---|---|---|
| Fukui Berry flux and Chern number of an isolated band or bundle | `--task chern` | [MoS₂ BZ map](examples/features/fukui-berry-curvature/), [Bi occupied bundle](examples/features/fukui-chern/) |
| 2D Fukui–Hatsugai Z₂ invariant and n-field | `--task z2` | [Bi topological insulator](examples/features/z2/) |
| Kubo Berry curvature on a BZ mesh or symmetry path | `--task kubo` | [MoS₂ occupied bundle and isolated-band maps, paths and bands](examples/features/kubo-curvature/) |
| Intrinsic charge Hall response versus chemical potential and temperature | `--task kubo-pairs`, then occupation-weighted postprocessing | [MoS₂ Hall and valley-region curves](examples/features/kubo-hall/) |
| Circular optical selectivity and transition spectra | `--task optical` / `--task spectrum` | [MoS₂ circular dichroism](examples/features/circular-dichroism/) |
| Real-space wavefunction at Γ | `--task wavefunction` | [MoS₂ wavefunction](examples/features/wavefunction/) |

Fukui plaquette flux and Kubo point curvature are different finite-grid
quantities. A Fukui integer needs band isolation and mesh checks. A Kubo
integral is not rounded to an integer; its k mesh and intermediate band
window must be converged. Native Kubo uses canonical momentum of the stored
pseudo-wavefunctions. Optional full-velocity comparisons assess the missing
PAW/nonlocal/SOC terms; see [operator choices](docs/OPERATOR_ROUTES.md).

## Usage

### First calculation: the supplied MoS₂ WAVECAR

Run from the repository root. This small example uses the actual SOC
band-path WAVECAR already in the repository, and writes the curvature of
the occupied bands 1–18 at its supplied k points:

```bash
mkdir -p results/mos2-path
./build/vaspberry --task kubo \
  --wavecar examples/1H-MoS2/KPATH/2.band/WAVECAR --spinor 2 --bands 1:18 --bundle 1 \
  --curvature-csv results/mos2-path/KUBO.csv
```

The CSV contains numerical point curvatures. The
[MoS₂ tutorial](examples/features/kubo-curvature/) separately provides a
matching full-BZ/path dataset and commands for the reference panels. A line
path alone cannot supply a BZ integral. Use a fresh output
path for each run.

### Apply the same commands to your system

Choose the mesh, band range and spinor count from your VASP calculation.
These illustrative commands assume the appropriate `WAVECAR` in the working
directory; the band indices are examples, not universal occupied counts.

```bash
# Fukui flux and Chern number: full 12 × 12 mesh, bands 1–18.
./build/vaspberry --task chern --wavecar WAVECAR \
  --mesh 12,12 --spinor 2 --bands 1:18 --output BERRYCURV

# Bi example: full even mesh, occupied spinor bands 1–10.
./build/vaspberry --task z2 --wavecar WAVECAR \
  --mesh 12,12 --spinor 2 --bands 1:10 --output NFIELD

# Occupied-bundle point curvature, excluding internal transitions.
./build/vaspberry --task kubo --wavecar WAVECAR --spinor 2 \
  --bands 1:18 --bundle 1 --curvature-csv KUBO_BUNDLE.csv

# Reusable all-band pair numerators for charge Hall postprocessing.
./build/vaspberry --task kubo-pairs --wavecar WAVECAR --spinor 2 \
  --pairs-csv PAIRS.csv
```

The native syntax groups a mesh as `NX,NY` and a band range as `FIRST:LAST`.
Specialized legacy options can still be used; for example, `--task kubo` is
`-kubo 2`, and `--bands 1:18` is `-ii 1 -if 18`. The optional
`--task kubo-integral` (`-kubo 1`) additionally evaluates the mesh integral.
Use `./build/vaspberry --help` for the full list. MPI uses the same arguments:

```bash
mpiexec -n 4 ./build/vaspberry-mpi --task chern --wavecar WAVECAR \
  --mesh 12,12 --spinor 2 --bands 1:18 --output BERRYCURV
```

For Z₂, use a nonmagnetic time-reversal-symmetric insulator and a full,
unshifted, even `Nx × Ny × 1` mesh with `Nx,Ny >= 4`, generated with
`ISYM=-1`. Report only a final `Z2_FIELD.csv` with `result_status=PASS`,
`reportable_invariant=1` and matching half-zone parities. See the
[Z₂ guide](docs/Z2_FUKUI_HATSUGAI.md) for input and convergence checks.

### Postprocess and plot saved results

The native program writes Berry-flux maps, point-curvature CSVs, pair data
and task-specific outputs. The supplied Python tools read those results:

- **Curvature maps and paths:** plot the computed samples, with optional
  interpolation for display; the [MoS₂ example](examples/features/kubo-curvature/)
  shows both maps and symmetry-path panels alongside the material's bands.
- **Charge Hall:** `import-pairs` checks and caches native pairs; `pair-hall`
  applies occupations, temperature and BZ/region integration. Reuse the cache
  for new scans. This is numerical postprocessing, in addition to plotting.
  See the [explicit three-stage workflow](docs/KUBO_TRANSPORT.md#native-pairs-to-charge-hall).
- **Figures and tables:** Hall output supports CSV, DAT and NPZ with JSON
  metadata; `tools/plot_hall.py` produces PNG, PDF and SVG. JSON records
  conditions and units, rather than replacing VASP input files.

Example `run.py` scripts and `wavecar-hall` remain optional convenience
wrappers that execute and record these stages. They are not needed to run
the Fortran program. For your own material, start with the
[transfer guide](examples/APPLY_TO_YOUR_SYSTEM.md).

## Optional extensions and supporting checks

The main tutorials above use ordinary VASP wavefunctions. Separate guides
cover extra operators and independent checks:

- [PAW optical and full-velocity comparisons](docs/OPERATOR_ROUTES.md),
  including [matched MoS₂ charge Hall data](examples/features/kubo-hall/operator-comparison/).
- [Conventional spin Hall response](docs/SPIN_HALL.md) with explicitly
  supplied spin and velocity matrices; [Bi results](examples/materials/bi-spin-hall/)
  complement its directly calculated Z₂ invariant.
- [Wannier-based supporting calculations](docs/WANNIER_TRANSPORT.md):
  externally prepared Wannier90 operators, with optional bundled NumPy
  response/strip calculations and independent postw90 comparisons. The Bi
  edge spectrum supports the Z₂ interpretation; the
  [MnBi₂Te₄ benchmark](examples/materials/mnbi2te4-qah/) illustrates a nonzero
  Chern invariant and integration convergence. These are optional validation
  routes, outside the WAVECAR/Fortran workflow.
- [Advanced geometric transport](docs/VALLEY_TRANSPORT.md),
  [matrix interfaces](docs/KUBO_TRANSPORT.md#applying-the-workflow-to-other-data)
  and [developer model checks](validation/models/).

The [technical report](docs/TECHNICAL_REPORT.md) connects these capabilities
to reference figures. The [material catalog](examples/materials/) records
input preparation and sampling. VASPBERRY implements the lattice method of
[Fukui, Hatsugai and Suzuki, JPSJ 74, 1674 (2005)](https://doi.org/10.1143/JPSJ.74.1674);
method-specific references are given in each guide.

## Contributors

* Hyun-Jung Kim: Main developer
* Sun-Woo Kim: Circular dichroism and Kubo formula

## Citation

```bibtex
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
```
