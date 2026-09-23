# Examples: Berry curvature and response functions

These examples illustrate VASPBERRY calculations for monolayer MoS₂, a
Bi bilayer and a three-septuple-layer MnBi₂Te₄ film. Each guide gives the VASP files, calculation parameters, commands,
expected numerical results and figures. The [technical report](../docs/TECHNICAL_REPORT.md)
introduces the methods and discusses the physical interpretation
([PDF version](../docs/TECHNICAL_REPORT.pdf)).

## Fukui Berry curvature in the Brillouin zone

![MoS2 Fukui Berry curvature map, band structure and marked symmetry path](features/fukui-berry-curvature/reference/smooth/figure.png)

**Monolayer MoS₂.** Occupied-band Berry curvature from the Fukui method,
shown in Cartesian reciprocal coordinates with the first Brillouin-zone
boundary and the K–Γ–K′ path marked. The adjacent panels show the matching
band structure and a curvature cut along that path. The K and K′ valleys carry opposite
curvature, while the full-zone Chern number vanishes. The figure comes from
a new full-mesh VASP and VASPBERRY calculation. The [Fukui Berry-curvature tutorial](features/fukui-berry-curvature/)
explains how to prepare the VASP mesh, run VASPBERRY and draw this map.

## Calculation guides

| Quantity | System and sampling | Result |
|---|---|---|
| [Fukui Berry curvature](features/fukui-berry-curvature/) | MoS₂, full 12×12 mesh | BZ map, marked K–Γ–K′ cut and band structure |
| [Occupied Chern number](features/fukui-chern/) | Bi, full 12×12 mesh | C = 0 for occupied bands 1–10 |
| [Z₂ invariant](features/z2/) | Bi, full 12×12 mesh | Z₂ = 1 |
| [Kubo Berry curvature](features/kubo-curvature/) | MoS₂, full 12×12 mesh and matching path | Occupied-bundle BZ map, path curvature and band structure |
| [Single-band valley curvature](features/kubo-curvature/valleys/) | MoS₂, two 9×9 K/K′ patches | Isolated band-18 maps, line cuts and valence bands |
| [Kubo Hall conductivity](features/kubo-hall/) | MoS₂, full meshes; chemical potential and temperature scans | Total and regional Hall curves, band reference and convergence studies |
| [Charge Hall conductivity](features/hall-valley/) | Bi, full mesh, T = 0 | Zero charge-Hall response in the insulating gap |
| [Circular optical transitions](features/circular-dichroism/) | MoS₂, K–Γ–K′ path | Opposite polarization selectivity at K and K′ |
| [Real-space wavefunction](features/wavefunction/) | MoS₂, Γ point | Spinor-state density in the atomic unit cell |
| [Magnetic Chern insulator](materials/mnbi2te4-qah/) | MnBi₂Te₄, three septuple layers | Nonzero occupied Chern number, band structure and Hall-integration checks |

Full-zone integrals require a complete periodic mesh. The supplied MoS₂
band-path WAVECAR serves the Kubo, optical and wavefunction examples. The
MoS₂ full-mesh WAVECAR is generated once using the supplied VASP preparation
and public charge density; its size is about 149 MB. The map/path figures
use a matching 49-point, 26-band VASP path generated with the same setup.
The full Bi WAVECAR is
available for direct recalculation.

The MnBi₂Te₄ material guide supplies the actual fixed SCF charge density,
structure and NSCF inputs for a 21-atom magnetic film. It combines the Fukui
and PAW optical workflows with an independent Wannier reference. Its VASP
and Wannier preparation steps are separate from the eight-feature batch command.

## Getting started

```bash
make serial
python3 -m pip install -r requirements-transport.txt
```

The [input guide](INPUTS.md) lists the files. For Bi, download the wavefunctions:

```bash
python3 examples/fetch_inputs.py bi --output-dir results/inputs/bi
```

Each guide shows the native VASPBERRY command. A short Python helper also
combines the calculation and plotting, for example:

```bash
python3 examples/features/kubo-curvature/run.py --output-dir results/mos2-kubo
python3 examples/features/z2/run.py \
  --wavecar results/inputs/bi/WAVECAR --output-dir results/bi-z2
```

Results are written to the chosen directory. Compare the numerical tables
and figures with the guide, then follow [applying the calculation to your
material](APPLY_TO_YOUR_SYSTEM.md). Set the mesh, bands and energy window from
your own VASP output.

## Parallel calculation

Native Fortran calculations can use MPI. For the MoS₂ Kubo example:

```bash
make mpi
mkdir -p results/mos2-kubo-mpi
(
  cd results/mos2-kubo-mpi
  mpiexec -np 2 ../../build/vaspberry-mpi \
    -f ../../examples/1H-MoS2/KPATH/2.band/WAVECAR \
    -s 2 -kubo 2 -ii 17 -if 18 -kubo_csv KUBO.csv -o BERRYCURV
)
```

Python plotting and transport helpers run serially. See the
[build guide](../docs/BUILD.md) for compiler-specific instructions.

## Further material

- [Material preparation](materials/): structures, VASP settings and available files.
- [Technical report](../docs/TECHNICAL_REPORT.md): methods, figures and interpretation.
- [Supplementary references](../docs/REFERENCE_MATERIALS.md): analytic comparisons and teaching exercises.
- [Validation details](../docs/VALIDATION_1.3.0.md): numerical and implementation checks.

The six calculations with supplied WAVECAR files can also be run together:

```bash
python3 examples/run_examples.py \
  fukui-chern z2 hall-valley kubo-curvature circular-dichroism wavefunction \
  --bi-wavecar results/inputs/bi/WAVECAR --output-dir results/all-examples
```

After preparing the MoS₂ full mesh, include all eight calculations with
`--all --mos2-mesh-wavecar /path/to/full-mesh/WAVECAR`, together with
`--bi-wavecar` and a new `--output-dir`.

The Kubo Hall tutorial also accepts this full-mesh input on its own:

```bash
python3 examples/run_examples.py kubo-hall \
  --mos2-mesh-wavecar /path/to/full-mesh/WAVECAR \
  --output-dir results/mos2-transport
```

Its preparation guide explains how to enlarge the mesh and retain extra empty
states for convergence checks. The general command combines native matrix
export and Python integration; saved pair data can be reused for additional
chemical potentials or temperatures. CSV, text DAT and NumPy NPZ are selectable
table formats; PNG, PDF and SVG are supported for figures.
