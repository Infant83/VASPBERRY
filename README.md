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
[Hands-on commands](docs/HANDS_ON.md) · [Feature examples](examples/README.md) · [Build guide](docs/BUILD.md) ·
[Output formats](docs/OUTPUT_FORMAT.md)

**Latest release: [1.4.1](https://github.com/Infant83/VASPBERRY/releases/tag/v1.4.1).**
The commands below are available in the fixed `v1.4.1` source. The original
short options remain supported. See the [release notes](docs/releases/v1.4.1.md),
[changelog](CHANGELOG.md), [version policy](docs/RELEASING.md)
and [migration notes](docs/MIGRATION.md).

## Download and compile

```bash
git clone --branch v1.4.1 --single-branch https://github.com/Infant83/VASPBERRY.git VASPBERRY-1.4.1
cd VASPBERRY-1.4.1
# Activate Intel oneAPI (or the equivalent compiler/MPI modules on your cluster).
source /opt/intel/oneapi/setvars.sh
make ifx-mpi
make check-ifx-mpi
mpiexec -n 4 ./build/vaspberry-ifx-mpi --help
```

The executable is `build/vaspberry-ifx-mpi`, built with Intel `mpiifx`,
Intel MPI and sequential LP64 oneMKL. Use the `mpiexec` from that same Intel
MPI installation. For Intel Classic, use `make ifort-mpi`,
`make check-ifort-mpi` and `build/vaspberry-ifort-mpi`. The
[build guide](docs/BUILD.md) also covers GNU/OpenMPI (`make gnu`) and
compiler-specific serial builds.
Python is needed for the supplied postprocessing and plotting tools:

```bash
python3 -m pip install -r requirements-transport.txt
```

For development work on the latest default branch:

```bash
git clone https://github.com/Infant83/VASPBERRY.git
```

Update an existing `master` checkout with `git pull --ff-only`. Release tags
stay fixed; source archives are available on the release page. Large Git LFS
inputs may need the [input-fetch procedure](examples/INPUTS.md).

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
`tools/procar_character.py`. Named atom/orbital groups and a Cartesian spin
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

For a shorter repeatable workflow, the existing `wavecar-hall` command runs
native Fortran pair export, cache import and Hall integration together. It
keeps `native/PAIRS.csv`, `pairs/pairs.npz` and `hall/conductivity.csv` under a
fresh output directory; plotting is a separate command. The
[Kubo tutorial](docs/KUBO_TRANSPORT.md#native-pairs-to-charge-hall) shows both
the shortcut and the explicit stages. Reuse the cache with `pair-hall` to
scan new μ/T values without repeating Fortran. For a new material, start with
the [transfer guide](examples/APPLY_TO_YOUR_SYSTEM.md).

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
