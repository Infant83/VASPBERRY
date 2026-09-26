# VASPBERRY

**Topology and response functions directly from VASP wavefunctions.**
VASPBERRY reads `WAVECAR` and evaluates wavefunction overlaps and interband
matrix elements in Fortran. Use Fukui methods for Chern numbers and the 2D Z₂
invariant, and Kubo calculations for Berry-curvature maps, symmetry paths and
intrinsic charge Hall response.

```text
VASP calculation → WAVECAR → VASPBERRY Fortran → numerical output → analysis / plots
```

VASP supplies the material's electronic structure. VASPBERRY postprocesses it;
band plots provide context for the calculated topology and response.

[Technical report](docs/TECHNICAL_REPORT.md) ([PDF](docs/TECHNICAL_REPORT.pdf)) ·
[Hands-on commands](docs/HANDS_ON.md) · [Feature examples](examples/README.md) · [Build guide](docs/BUILD.md) ·
[Postprocessing guide](docs/POSTPROCESSING.md) · [Output formats](docs/OUTPUT_FORMAT.md)

**Latest release: [1.5.0](https://github.com/Infant83/VASPBERRY/releases/tag/v1.5.0).**
The commands below are available in the fixed `v1.5.0` source. The original
short options remain supported. See the [release notes](docs/releases/v1.5.0.md),
[changelog](CHANGELOG.md), [version policy](docs/RELEASING.md)
and [migration notes](docs/MIGRATION.md).

## Install and run the Fortran program

The native program requires a **Fortran compiler, GNU Make, POSIX shell/tools
and LP64 BLAS/LAPACK libraries**. An MPI build additionally needs the matching
MPI **development package/compiler wrapper and runtime launcher**. These
native dependencies are sufficient to compile and run VASPBERRY on an existing
compatible `WAVECAR`. Python is used by the optional numerical postprocessing
and plotting tools.

| Environment | Required native development environment | Build and executable |
|---|---|---|
| Linux, Intel serial | Intel oneAPI `ifx` + oneMKL | `make ifx` → `build/vaspberry-ifx` |
| Linux, Intel MPI (main walkthrough) | `ifx` + Intel MPI SDK (`mpiifx`) + oneMKL | `make ifx-mpi` → `build/vaspberry-ifx-mpi` |
| Linux, GNU serial | GNU Fortran + LP64 BLAS/LAPACK | `make serial` → `build/vaspberry` |
| Linux, GNU MPI | GNU Fortran + Open MPI or MPICH development libraries + LP64 BLAS/LAPACK | `make mpi` → `build/vaspberry-mpi` |
| macOS, GNU | Compatible GNU Fortran/MPI libraries and LP64 BLAS/LAPACK | See the architecture-specific [build guide](docs/BUILD.md) for paths and available packages |

Retained Intel Classic installations use `make ifort` or `make ifort-mpi`.
On Windows use a Linux environment such as WSL2 and its Linux compiler stack;
there is no native Windows build. The [build guide](docs/BUILD.md) distinguishes
CI-validated environments from installation guidance and covers cluster modules,
macOS, MPICH, runtime libraries and common errors.

### Get the source

Clone the repository's default branch (`master`):

```bash
git clone https://github.com/Infant83/VASPBERRY.git
cd VASPBERRY
make help
```

For a fixed release, download and extract the source archive from
[v1.5.0](https://github.com/Infant83/VASPBERRY/releases/tag/v1.5.0), then run the
same build commands from the extracted directory containing `Makefile`.
**Archive builds do not require Git or Python.** No precompiled executable or
system-wide installation is needed; the build creates local files in `build/`.
Update an existing default-branch checkout with `git pull --ff-only`, then
rebuild. Release tags stay fixed. Large example inputs may separately need
the [input-fetch procedure](examples/INPUTS.md).

### Intel oneAPI and Intel MPI on Linux

Install/load the Fortran compiler, Intel MPI SDK and oneMKL development
components; having only runtime libraries is insufficient. In Bash, activate
the installed environment or use the site's equivalent modules:

```bash
source /opt/intel/oneapi/setvars.sh
make ifx-mpi
make check-ifx-mpi
mpiexec -n 4 ./build/vaspberry-ifx-mpi --help
```

Replace the setup path with the site's actual oneAPI installation. Use the
`mpiexec` from the same Intel MPI environment as `mpiifx`, and run within an
allocation that permits the requested rank count. `check-ifx-mpi` builds and
checks two-rank startup/communication/help. Intel Classic uses
`make check-ifort-mpi`. For a machine without MPI, use `make ifx`,
`make check-ifx`, and run `./build/vaspberry-ifx --help` directly.

### GNU on Ubuntu/Debian

The serial executable needs no MPI package:

```bash
sudo apt-get update
sudo apt-get install make gfortran libblas-dev liblapack-dev
make serial
make check-serial-help
./build/vaspberry --help
```

For Open MPI, add the development/runtime packages and build the MPI version:

```bash
sudo apt-get install openmpi-bin libopenmpi-dev
make mpi
make check-gnu
mpiexec -n 4 ./build/vaspberry-mpi --help
```

On a cluster, use its installed modules/libraries instead of these administrator
commands. Compiler, wrapper and library paths can be supplied to Make; see
[build overrides and platform recipes](docs/BUILD.md). The supplied builds
use byte-based WAVECAR records and **LP64**, not ILP64, numerical libraries.
Serial executables run directly; MPI executables run with their matching
launcher. Keep the compiler/MPI/library environment loaded when running.

Use your usual Python environment for the optional `tools/` commands. Their
full dependency versions are in [requirements-transport.txt](requirements-transport.txt);
Matplotlib is used by the supplied postprocessing and plotting tools.
Environment creation and package installation follow your site's usual practice.

## Features

| Task | Native command | Actual VASP example and results |
|---|---|---|
| Fukui Berry flux and Chern number of an isolated band or bundle | `--task chern` | [MoS₂ BZ map](examples/features/fukui-berry-curvature/), [Bi occupied bundle](examples/features/fukui-chern/) |
| 2D Fukui–Hatsugai Z₂ invariant and n-field | `--task z2` | [MoS₂ (Z₂ = 0) and Bi (Z₂ = 1)](examples/features/z2/comparison/) |
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

For **layer, atom, orbital and spin character**, combine a matching SOC
`PROCAR` with the WAVECAR and actual spin frame from `OUTCAR` using
the optional `[projection]` and `[group NAME]` sections in the
[postprocessing settings](docs/POSTPROCESSING.md). Named atom/orbital groups and a Cartesian spin
axis define the projections. The same saved native pair data support
chemical-potential/temperature scans of selected-band projected charge-Hall
contributions. Follow the [projection tutorial](examples/features/procar-character/)
for commands and an explicitly synthetic reproducibility fixture. These
projections explain state character; conventional spin-current Hall uses the
separate operator route described in the [spin guide](docs/SPIN_HALL.md).

## Usage

### First calculation: the supplied MoS₂ WAVECAR

Run from the repository root. This small example uses the actual SOC
band-path WAVECAR already in the repository, and writes the curvature of
the occupied bands 1–18 at its supplied k points:

```bash
mkdir -p results/mos2-path
mpiexec -n 4 ./build/vaspberry-ifx-mpi --task kubo \
  --wavecar examples/1H-MoS2/KPATH/2.band/WAVECAR --spinor 2 --bands 1:18 --bundle 1 \
  --curvature-csv results/mos2-path/KUBO.csv
```

`KUBO.csv` contains one row per k point/spin channel, with fractional
k coordinates, `omega_z_A2` (Ωxy in Å²) and `min_external_gap_eV`.
The occupied-bundle sum is already evaluated by Fortran: plot `omega_z_A2`
against the ordered `k_index` for this path. For a mesh map, convert the
fractional coordinates with the reciprocal lattice. The
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
mpiexec -n 4 ./build/vaspberry-ifx-mpi --task chern --wavecar WAVECAR \
  --mesh 12,12 --spinor 2 --bands 1:18 --output BERRYCURV

# Bi example: full even mesh, occupied spinor bands 1–10.
mpiexec -n 4 ./build/vaspberry-ifx-mpi --task z2 --wavecar WAVECAR \
  --mesh 12,12 --spinor 2 --bands 1:10 --output NFIELD

# Occupied-bundle point curvature, excluding internal transitions.
mpiexec -n 4 ./build/vaspberry-ifx-mpi --task kubo --wavecar WAVECAR --spinor 2 \
  --bands 1:18 --bundle 1 --curvature-csv KUBO_BUNDLE.csv

# Reusable all-band pair numerators for charge Hall postprocessing.
mpiexec -n 4 ./build/vaspberry-ifx-mpi --task kubo-pairs --wavecar WAVECAR --spinor 2 \
  --pairs-csv PAIRS.csv
```

The native syntax groups a mesh as `NX,NY` and a band range as `FIRST:LAST`.
Specialized legacy options can still be used; for example, `--task kubo` is
`-kubo 2`, and `--bands 1:18` is `-ii 1 -if 18`. The optional
`--task kubo-integral` (`-kubo 1`) additionally evaluates the mesh integral.
Use `mpiexec -n 4 ./build/vaspberry-ifx-mpi --help` or the [native command reference](docs/NATIVE_COMMANDS.md) for the full list. For example:

```bash
mpiexec -n 4 ./build/vaspberry-ifx-mpi --task chern --wavecar WAVECAR \
  --mesh 12,12 --spinor 2 --bands 1:18 --output BERRYCURV
```

For Z₂, use a nonmagnetic time-reversal-symmetric insulator and a full,
unshifted, even `Nx × Ny × 1` mesh with `Nx,Ny >= 4`, generated with
`ISYM=-1`. Report only a final `Z2_FIELD.csv` with `result_status=PASS`,
`reportable_invariant=1` and matching half-zone parities. See the
[Z₂ guide](docs/Z2_FUKUI_HATSUGAI.md) for input and convergence checks.

### Postprocess and plot saved results

Fortran produces the numerical outputs. Python has separate calculation
and plotting roles; neither stage requires a new VASP run when reusing the
same saved electronic structure.

| Stage | Input → output | What the output means / possible figure |
|---|---|---|
| Native `--task kubo` | `WAVECAR` → `KUBO.csv` | Gap-divided point curvature (Å²); plot Ω(k) maps or symmetry-path curves directly. |
| Native `--task kubo-pairs` | `WAVECAR` → `PAIRS.csv` | Energies, k coordinates and three interband pair numerators (eV² Å²). These are reusable intermediate matrix data, not yet a Hall conductivity. |
| Python `import-pairs` | `PAIRS.csv` + matching `WAVECAR` → `pairs.npz`, `pairs.json` | Checks normalization/mesh/spin metadata and caches arrays; no Hall integration. |
| Python `pair-hall` | Pair cache + μ/T/region choices → `conductivity.csv`, `.dat`, `.npz`, `.json` | Applies occupations, energy denominators and BZ/region integration; sheet σ and Δσ in e²/h, plus carriers per cell. |
| Python `plot_hall.py` | `conductivity.csv` → PNG/PDF/SVG | Draws σ(μ) or Δσ(μ) for selected temperatures and total/valley regions; no wavefunction calculation. |
| Python `procar_character.py` | Matching `PROCAR`, `WAVECAR`, `OUTCAR`, atom/orbital groups and spin axis → character tables; optional pair cache → projected Hall tables | Plots state-character maps and selected-band/group charge-Hall contributions. |

The [hands-on guide](docs/HANDS_ON.md) gives each stage's command, filenames,
and independent Matplotlib examples. The [output specification](docs/OUTPUT_FORMAT.md)
defines columns, units and metadata for reading the same files in Python,
Origin, gnuplot or other analysis tools. Follow the
[MoS₂ curvature](examples/features/kubo-curvature/),
[Hall](examples/features/kubo-hall/) and
[PROCAR](examples/features/procar-character/) examples for concrete figures.

For repeated Hall and PROCAR analysis, keep the input paths, μ/T scan and
optional atom groups or regions in **one commented settings file**:

```bash
python3 tools/vaspberry_post.py run analysis.ini
python3 tools/vaspberry_post.py plot results/run01
```

The first command launches the compiled Fortran executable and the requested
numerical postprocessing. The second draws the saved tables. `analysis.ini`
specifies `results/run01` relative to that file; it also sets the executable
and MPI rank count. No user-written JSON is needed. Start with the
[complete public Bi example](examples/features/simple-postprocess/) and
[settings guide](docs/POSTPROCESSING.md). Reuse a previous run's pair cache
with `run analysis-next.ini --reuse results/run01` when changing μ/T, without
repeating Fortran. The individual [Kubo commands](docs/KUBO_TRANSPORT.md#native-pairs-to-charge-hall)
and `procar_character.py` remain available for advanced controls.

## Optional extensions and supporting checks

The main tutorials above use ordinary VASP wavefunctions. Separate guides
cover extra operators and independent checks:

- [PAW optical and full-velocity comparisons](docs/OPERATOR_ROUTES.md),
  including [matched MoS₂ charge Hall data](examples/features/kubo-hall/operator-comparison/).
- [Conventional spin Hall response](docs/SPIN_HALL.md) with explicitly
  supplied spin and velocity matrices; [Bi results](examples/materials/bi-spin-hall/)
  complement its directly calculated Z₂ invariant.
- [Advanced geometric transport](docs/VALLEY_TRANSPORT.md),
  [matrix interfaces](docs/KUBO_TRANSPORT.md#applying-the-workflow-to-other-data)
  and [developer model checks](validation/models/).

The [technical report](docs/TECHNICAL_REPORT.md) connects these capabilities
to reference figures. The [material catalog](examples/materials/) records
input preparation and sampling. VASPBERRY implements the lattice method of
[Fukui, Hatsugai and Suzuki, JPSJ 74, 1674 (2005)](https://doi.org/10.1143/JPSJ.74.1674);
method-specific references are given in each guide.

## Contributors

* Hyun-Jung Kim: Main developer and maintainer; responsible for subsequent
  development and ongoing updates of the Kubo implementation.
* Sun-Woo Kim: Contributions to circular dichroism and the initial Kubo
  implementation.

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
