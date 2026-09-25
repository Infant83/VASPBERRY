# MoS₂ Berry curvature over the first Brillouin zone

The occupied bands of monolayer 1H-MoS₂ have opposite Berry curvature near the
K and K′ valleys. Time-reversal symmetry makes the total Chern number vanish
even though the local curvature is finite. This example calculates the **Fukui
plaquette curvature Ωz(kx, ky)** from an actual VASP spinor WAVECAR and plots
it in Cartesian reciprocal coordinates inside the hexagonal first Brillouin zone.

![MoS2 Fukui curvature map, matching bands and symmetry-path cut](reference/smooth/figure.png)

The reference above was calculated with the current native VASPBERRY routine
from a newly generated, complete **12 × 12 × 1 VASP mesh**. The public structure
and SCF charge density provide the starting inputs. The generated 149 MB
WAVECAR is not bundled; the [VASP preparation instructions](inputs/README.md)
give the executed setup and commands to generate it with your licensed VASP
and PAW datasets. The supplied K–Γ–K′ WAVECAR is a line calculation and cannot
replace this two-dimensional mesh.

## Calculation and inputs

| Quantity | Reference calculation |
|---|---|
| System | Nonmagnetic monolayer 1H-MoS₂, SOC |
| VASP calculation | Fixed public SCF density, `ICHARG=11`, complete Γ-centered 12 × 12 × 1 mesh |
| Basis | 400 eV cutoff, 26 stored spinor bands, 144 k points |
| Selected subspace | Occupied bands 1–18 |
| Sampled occupied-to-empty gap | 1.67355 eV; the sampled global gap is also positive |
| Fukui sampling | 144 independent plaquettes |
| Input files | [POSCAR](inputs/POSCAR), [INCAR](inputs/INCAR), [KPOINTS](inputs/KPOINTS), public [CHGCAR.gz](../../1H-MoS2/KPATH/1.scf/CHGCAR.gz) |
| VASP reference output | [OUTCAR](reference/vasp/OUTCAR), [EIGENVAL](reference/vasp/EIGENVAL), [OSZICAR](reference/vasp/OSZICAR) |

The NSCF step took about 144 seconds and 602 MB peak resident memory on the
reference machine. The native Fukui step took about one minute. These timings
indicate the scale of this small example; they are not a convergence study.

## 1. Generate the full-mesh WAVECAR

From the repository root, prepare the working directory using the matching
licensed Mo/S PAW-PBE datasets:

```bash
python3 examples/features/fukui-berry-curvature/prepare_vasp.py \
  --potcar /path/to/licensed/Mo-S/POTCAR \
  --output-dir results/mos2-fullmesh-vasp

cd results/mos2-fullmesh-vasp
/path/to/vasp_ncl > vasp.stdout.log 2> vasp.stderr.log
cd ../..
```

The [input guide](inputs/README.md) describes the charge density, potential
specification, MPI launch option and expected convergence checks. Complete
this VASP calculation before the next step.

## 2. Run native Fortran on WAVECAR

```bash
make serial
mkdir -p results/mos2-fukui-native
cd results/mos2-fukui-native
../../build/vaspberry \
  --wavecar ../mos2-fullmesh-vasp/WAVECAR \
  --task chern --spinor 2 --mesh 12,12 --bands 1:18 --output BERRYCURV \
  > vaspberry.log
cd ../..
```

`--spinor 2` reads both components of each SOC spinor. `--bands 1:18` uses the
determinant of the occupied-subspace overlap matrix, allowing degeneracy
within that subspace. `--mesh 12,12` describes the WAVECAR mesh; it does not
generate missing k points. This Fukui calculation does not use the Kubo
sum over empty intermediate states.

The sign convention is `phi = −Arg(product of link determinants)` around
`k → k+dk1 → k+dk1+dk2 → k+dk2 → k`. The reported curvature is
`Omega_z = phi / deltaS`, and the Chern number is `sum(phi)/(2*pi)`.
Since `deltaS` has units Å⁻², Ωz has units **Å²**.

## 3. Plot the native curvature

For a first look, plot the completed native file directly:

```bash
python3 -m pip install -r requirements-transport.txt
python3 tools/plot_berry_curvature.py \
  --input results/mos2-fukui-native/BERRYCURV.dat \
  --poscar examples/features/fukui-berry-curvature/inputs/POSCAR \
  --output results/mos2-fukui-native/curvature.png --title '1H-MoS2'
```

This reads the existing plaquette values and uses the matching POSCAR to draw
the Cartesian first Brillouin zone. It does not recalculate wavefunction overlaps.

### Add the band structure and symmetry-path cut

The figure pairs the BZ map with bands and curvature along **K–Γ–K′**.
The dashed line in the map is the path used in both right-hand panels.
K is (1/3, 2/3), Γ is (0, 0), and K′ is −K in reciprocal coordinates.

First prepare a matching path calculation from the completed full-mesh VASP
directory. It retains the same charge density, potentials, 400 eV cutoff and
26 bands; only the k-point list changes:

```bash
python3 examples/features/fukui-berry-curvature/prepare_path.py \
  --mesh-dir results/mos2-fullmesh-vasp \
  --output-dir results/mos2-path-vasp

cd results/mos2-path-vasp
/path/to/vasp_ncl > vasp.stdout.log 2> vasp.stderr.log
cd ../..

python3 tools/plot_berry_panels.py \
  --method fukui --input results/mos2-fukui-native/BERRYCURV.dat \
  --poscar examples/features/fukui-berry-curvature/inputs/POSCAR \
  --path-wavecar results/mos2-path-vasp/WAVECAR \
  --path-node-indices 1 25 49 --path-labels K Gamma Kprime \
  --map-style smooth --display-grid 401 \
  --title '1H-MoS2' --output results/mos2-fukui-native/panels.png
```

The [path input guide](inputs/path/README.md) provides the 49-point KPOINTS
and preparation details. The path NSCF step took about 55 seconds and 367 MB
peak resident memory on the reference machine.

The left panel uses **NumPy periodic bilinear interpolation** of the native
plaquette values onto a 401×401 Cartesian display grid, clipped to the first
BZ. The original samples and their integral are unchanged. Use
`--map-style cells` for the un-smoothed cell display. The lower-right
line is a **periodic bilinear cut of the 12×12 plaquette field** at the path
coordinates; the line does not add newly calculated pointwise curvature.
The upper-right panel uses the actual path VASP energies, relative to the
valence-band maximum. Occupied bands are blue and empty bands gray. Both
right-hand panels use cumulative Cartesian path distance.

To redraw the supplied results directly, replace `--path-wavecar` with
`--bands-csv examples/features/fukui-berry-curvature/reference/path/bands.csv`
and use the supplied `reference/BERRYCURV.dat` as `--input`. The plotter also
accepts PDF output. It writes the line data and a plotting record beside the
figure; exported band energies retain their original VASP zero.

For a standalone map, `tools/plot_berry_curvature.py` remains available.

## Optional reproduction helper

An optional runner repeats the native calculation and makes the standalone map.
It verifies the complete mesh, occupied
window and sampled gap, and checks the expected zero integral and
opposite-sign time-reversed curvature:

```bash
python3 examples/features/fukui-berry-curvature/run.py \
  --wavecar results/mos2-fullmesh-vasp/WAVECAR \
  --output-dir results/mos2-fukui
```

The runner writes both PNG and PDF figures, a numerical CSV and a machine-readable
result record. Its input is the VASP WAVECAR. JSON files record results and provenance.

## Reference results

| Output | Contents |
|---|---|
| [Panel PNG](reference/smooth/figure.png), [PDF](reference/smooth/figure.pdf) | BZ map, marked path, band structure and curvature cut shown above |
| [Path curve](reference/smooth/path_curvature.csv), [bands](reference/path/bands.csv) | Plotted line values and unchanged 26-band VASP energies |
| [Path EIGENVAL](reference/path/EIGENVAL), [OUTCAR](reference/path/OUTCAR) | Matching 49-point VASP calculation |
| [Standalone map](reference/figure.png), [PDF](reference/figure.pdf) | Original Cartesian first-BZ map |
| [BERRYCURV.dat](reference/BERRYCURV.dat) | Current native Fukui output; only the workstation path in its comment header is normalized |
| [summary.csv](reference/summary.csv) | 144 unique plaquette centers, Ωz and flux; fractional centers use a centered primitive cell |
| [result.json](reference/result.json) | Input identity, numerical checks, units and plotting convention |
| [crosscheck.json](reference/crosscheck.json) | Independent production Python FHS comparison using the same actual WAVECAR |
| [vasp/run.json](reference/vasp/run.json) | VASP input/output identities, convergence and measured resources |

The plaquette area is **0.031470265 Å⁻²**. Curvature ranges from
**−12.3125 to +12.3124 Å²**. Integrating the four-decimal native printed values
gives **C ≈ −5.01 × 10⁻⁷**; the native calculation reports **0.0000**, and the
independent production Python FHS calculation gives **5.65 × 10⁻¹⁶**.
Time-reversed native plaquettes cancel to within **2 × 10⁻⁴ Å²**.
These residuals reflect numerical precision and the finite VASP convergence
threshold. The opposite valleys remain individually visible although their
total Chern number is zero.

The native file repeats the 144 plaquettes over a 3 × 3 display region.
Integrate **one** set of unique plaquettes. The runner deduplicates the
display copies before computing the integral.

The reference is a verified calculation at one mesh and cutoff. For
quantitative material predictions, converge the SCF density, cutoff and
k mesh; a 12 × 12 color map does not establish a converged valley peak.

### Historical comparison

The original [contour script](../../1H-MoS2/contour.py) used SciPy
`griddata` with its default linear interpolation onto a 300×300 Cartesian
grid. The current smooth panels use NumPy bilinear interpolation with
explicit periodic wrapping; both are display operations. The earlier
[cell panels](reference/map-path/figure.png) are also retained.

The original public [BERRYCURV.dat](../../1H-MoS2/BERRYCURV.dat) is preserved.
Its separate [archived figure](reference/archive/figure.png) and
[reference record](reference/archive/result.json) have extrema ±12.233671 Å²
and a printed-value integral −4.01 × 10⁻⁸. The matching historical full-mesh
WAVECAR is unavailable. This older output and the new public-density NSCF
calculation show the same valley pattern; they are distinct calculations and
are not claimed to be bitwise reproductions of each other.

To replot that historical output without rerunning wavefunction overlaps:

```bash
python3 examples/features/fukui-berry-curvature/run.py \
  --archive-reference --output-dir results/mos2-fukui-archive
```

The historical file's `A^-2` curvature label is a typographical error; its
native `phi/deltaS` quantity has units Å². Its numerical file is unchanged.

## Apply the workflow to another material

Supply its complete uniform mesh, matching POSCAR and appropriate spinor
setting. Select an isolated band or isolated group of bands, then change
the band indices and mesh dimensions in the native command. The optional
runner above checks this specific 12 × 12 MoS₂ reference; use the general
native CLI and plotting tool for other systems. Plot with the matching
POSCAR so the Cartesian Brillouin zone follows the actual reciprocal lattice.

For a supplied-WAVECAR calculation of a topological invariant, see the
[Bi occupied-subspace example](../fukui-chern/). For a pointwise Kubo map and matching path, see the [Kubo example](../kubo-curvature/).
