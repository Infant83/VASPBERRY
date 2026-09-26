# Examples: calculate directly from VASP wavefunctions

VASPBERRY reads **VASP WAVECAR files directly**. Its native Fortran program
calculates Fukui curvature and Chern numbers, Z₂, Kubo curvature, circular
optical transitions and real-space wavefunctions. Use the resulting numerical files in your own analysis
or the supplied Python plotting tools.

Build with `make serial` and run `build/vaspberry --help` for the native CLI.
The older executable name `build/vaspberry-gfortran` and short flags remain
supported; the guides use descriptive task and option names.

The core guides follow the same sequence: **VASP input → native command →
output file → postprocessing and figure**. The [technical report](../docs/TECHNICAL_REPORT.md)
explains the methods and material results ([PDF](../docs/TECHNICAL_REPORT.pdf)).
The [report reproduction table](REPORT_REPRODUCTION.md) maps every report
figure to its calculation, numerical files and input-regeneration requirements.
Start with the [hands-on commands](../docs/HANDS_ON.md); use the compact
[native command reference](../docs/NATIVE_COMMANDS.md) when changing task options.
For a short Hall workflow, the [single-settings-file Bi example](features/simple-postprocess/)
uses one command to run Fortran plus integration, and one to draw the saved results.

## Choose a calculation

| Quantity | Real VASP input | Native result and tutorial |
|---|---|---|
| [Fukui Berry curvature](features/fukui-berry-curvature/) | MoS₂, full 12×12 mesh | `BERRYCURV.dat`; opposite K/K′ curvature in a Cartesian BZ map |
| [Occupied Chern number](features/fukui-chern/) | Bi, full 12×12 mesh | `BERRYCURV.dat`; C = 0 for bands 1–10 |
| [Z₂ invariant](features/z2/comparison/) | MoS₂ and Bi, full 12×12 SOC meshes | `Z2_FIELD.csv`; Z₂ = 0 and 1, with paired n-field maps and half-zone diagnostics |
| [Kubo Berry curvature](features/kubo-curvature/) | MoS₂, full mesh and matching K–Γ–K′ path | Native Kubo CSV; occupied-bundle curvature and band panels |
| [Single-band valley curvature](features/kubo-curvature/valleys/) | MoS₂, two 9×9 K/K′ patches | Native Kubo CSV; isolated band-18 maps and line cuts |
| [Charge and regional Hall response](features/kubo-hall/) | MoS₂, complete meshes and sufficient empty bands | Native `PAIRS.csv`, then occupation-weighted integration and Hall tables |
| [Circular optical transitions](features/circular-dichroism/) | Supplied MoS₂ K–Γ–K′ WAVECAR | Native left/right spectra; opposite valley selectivity |
| [Real-space wavefunction](features/wavefunction/) | Γ in the supplied MoS₂ path WAVECAR | Native real/imaginary amplitude grids; spinor-state density |

For an additional insulating charge-Hall check, the [Bi Hall guide](features/hall-valley/)
integrates occupied-subspace Fukui flux and obtains a zero charge response.
Kubo Hall occupation weighting and BZ integration currently use the bundled
Python tools after native matrix-element export; its guide shows both
stages explicitly. Python plotting reads the completed numerical outputs.

## First calculation: the supplied Bi wavefunctions

Run from the repository root. Download the actual Bi WAVECAR once, then
calculate the occupied-subspace Chern number with native Fortran:

```bash
make serial
python3 examples/fetch_inputs.py bi --output-dir results/inputs/bi
repo_dir="$PWD"
mkdir -p results/bi-fukui
(
  cd results/bi-fukui
  "$repo_dir/build/vaspberry" \
    --wavecar "$repo_dir/results/inputs/bi/WAVECAR" \
    --task chern --spinor 2 --mesh 12,12 --bands 1:10 --output BERRYCURV > vaspberry.log
)
```

Read `results/bi-fukui/vaspberry.log` and `BERRYCURV.dat`: the reference has
**C = 0**, consistent with time-reversal symmetry. Plot that output with:

```bash
python3 -m pip install -r requirements-transport.txt
python3 tools/plot_berry_curvature.py \
  --input results/bi-fukui/BERRYCURV.dat \
  --poscar examples/Bi_Z2/inputs/POSCAR \
  --output results/bi-fukui/curvature.png --title 'Bi: occupied bands'
```

The Bi map is zero at native output precision. For a finite valley-curvature
map, follow the MoS₂ guide below. The [input guide](INPUTS.md) distinguishes
supplied WAVECAR files from inputs that require a VASP preparation step.

## MoS₂: finite curvature with zero total Chern number

![MoS2 native Fukui map with matching bands and symmetry-path cut](features/fukui-berry-curvature/reference/smooth/figure.png)

The K and K′ valleys have opposite occupied-band curvature. The full-zone
Chern number vanishes. The left panel uses native Fukui plaquettes in
Cartesian reciprocal coordinates; the other panels show a matching VASP
band structure and an interpolated cut through the plaquette field.

The [Fukui guide](features/fukui-berry-curvature/) supplies the executed
VASP preparation, direct Fortran command, numerical references and plotting
commands. Generate its complete 12×12 WAVECAR once from the public SCF density;
the file is about 149 MB. The matching path has 49 points and 26 bands.
The separately supplied historical 48-point, 32-band path is used by the
optical and real-space examples. A path cannot replace a full integration mesh.

## Apply the commands to your material

Keep the native command structure and change the input path, spinor setting,
band selection and k mesh to match your VASP output. Start from
[applying the workflow to your system](APPLY_TO_YOUR_SYSTEM.md), then inspect
the [output conventions](../docs/OUTPUT_FORMAT.md). Verify the relevant gaps,
occupations and convergence before interpreting a material result.

Native calculations support MPI:

```bash
make mpi
mkdir -p results/mos2-kubo-mpi
(
  cd results/mos2-kubo-mpi
  mpiexec -np 2 ../../build/vaspberry-mpi \
    --wavecar ../../examples/1H-MoS2/KPATH/2.band/WAVECAR \
    --spinor 2 --task kubo --bundle 1 --bands 1:18 --curvature-csv KUBO.csv
)
```

The mesh and band choices still come from the supplied WAVECAR. See the
[build guide](../docs/BUILD.md) for compiler and MPI instructions. Python
integration and plotting run separately from the native MPI calculation.

## Material studies and optional extensions

These studies retain their VASP preparation, physical results and references.
Their additional operators or interpolation steps are described in each guide.

| Study | Scope |
|---|---|
| [MoS₂ stacking and valley selection](materials/mos2-stacking-valley/) | Monolayer and 1H/2H/3R bilayers; direct WAVECAR optics and optional standard-VASP PAW optical comparison |
| [MnBi₂Te₄ magnetic Chern insulator](materials/mnbi2te4-qah/) | Three-septuple-layer film; direct occupied Chern calculation, PAW optics and optional VASP-derived Wannier Hall/interpolation comparisons |
| [Bi quantum spin Hall and PAW spin response](materials/bi-spin-hall/) | Fresh SCF input, bulk Z₂, physical spin/full-velocity operators, and optional Wannier edge connectivity |
| [Matched Hall operator comparison](features/kubo-hall/operator-comparison/) | Canonical-momentum and full PAW velocity responses on the same electronic states |

Wannier interpolation is an additional route for validated dense-mesh and
edge-spectrum work; the direct WAVECAR calculations above remain complete
workflows. Physical PAW spin/full-velocity export uses a separate licensed
serial VASP producer. Its build requirements are independent of native
VASPBERRY MPI support. See the [material catalogue](materials/),
[spin Hall guide](../docs/SPIN_HALL.md) and
[Wannier guide](../docs/WANNIER_TRANSPORT.md).

## Optional reproduction helpers

Each tutorial may provide a `run.py` helper that combines calculation,
checks and plots using its fixed reference settings. The native-feature helpers launch
the Fortran executable and are optional conveniences for repeating examples.
The supplementary Bi gap-Hall helper instead evaluates Fukui overlaps in
Python as an independent transport check.
Their `--postprocess-only` or `--plot-only` modes, where documented, read
completed outputs without rerunning Fortran.

The six examples with supplied WAVECAR inputs can be repeated together:

```bash
python3 examples/run_examples.py \
  fukui-chern z2 hall-valley kubo-curvature circular-dichroism wavefunction \
  --bi-wavecar results/inputs/bi/WAVECAR --output-dir results/all-examples
```

After preparing the MoS₂ full mesh, use `--all` with
`--mos2-mesh-wavecar /path/to/full-mesh/WAVECAR`, `--bi-wavecar` and a new
`--output-dir` to include all eight catalogue calculations. Material studies
and their VASP/Wannier preparation are separate from this batch.

[Supplementary references](../docs/REFERENCE_MATERIALS.md) and
[validation details](../docs/VALIDATION_1.3.0.md) retain the model checks and
implementation tests supporting the real-material examples.
