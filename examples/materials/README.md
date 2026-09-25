# Material datasets and calculation routes

Choose the material data, then follow the linked native WAVECAR tutorial.
The core calculations read VASP wavefunctions directly and do not require
Wannierization. Each guide identifies its VASP inputs, native outputs and
reference figures.

## Direct WAVECAR examples

| Material and sampling | Available input | Calculation and result |
|---|---|---|
| [1H-MoS₂, full 12×12 mesh](../1H-MoS2/) | Public SCF density, structure and complete NSCF preparation; generate the 149 MB spinor WAVECAR with licensed VASP | [Fukui BZ curvature](../features/fukui-berry-curvature/), [bundle Kubo curvature](../features/kubo-curvature/), [native pair export and Hall integration](../features/kubo-hall/) |
| [1H-MoS₂, K–Γ–K′ path](../1H-MoS2/KPATH/) | Supplied 48-point, 32-band WAVECAR and EIGENVAL | [Path Kubo](../features/kubo-curvature/), [circular optics](../features/circular-dichroism/), [Γ wavefunction](../features/wavefunction/) |
| [Bi buckled honeycomb, full 12×12 mesh](../Bi_Z2/) | Downloadable SOC WAVECAR, structure and native references | [Occupied Chern number](../features/fukui-chern/), [Z₂](../features/z2/), [gap charge Hall check](../features/hall-valley/) |

Full-BZ Chern and Hall integrals require a complete periodic mesh. A path
WAVECAR supplies the states along that path only. The MoS₂ map/path guide
also prepares a matching 49-point, 26-band path from its full-mesh setup.

Download the historical Bi WAVECAR with:

```bash
python3 examples/fetch_inputs.py bi --output-dir results/inputs/bi
```

The file is about 200 MB. A Git LFS pointer is not the wavefunction payload;
`git lfs pull --include="examples/Bi_Z2/WAVECAR"` is an alternative when LFS is
available. The original Bi SCF provenance is incomplete, so this dataset
reproduces the WAVECAR postprocessing calculation rather than its preceding DFT run.

## Additional material studies

| Material | Direct calculation and additional capabilities |
|---|---|
| [MoS₂ monolayer and 1H/2H/3R bilayers](mos2-stacking-valley/) | Fixed structures and actual VASP band energies; native circular optics, with optional standard-VASP PAW optical strengths |
| [MnBi₂Te₄, three septuple layers](mnbi2te4-qah/) | Public fixed SCF density and 21-atom VASP preparation; occupied Chern and PAW optical diagnostics, plus supplied Wannier operators for interpolation and Hall comparisons |
| [Bi PAW spin Hall and ideal edges](bi-spin-hall/) | Fresh SCF density and physical spin/full-velocity preparation; gapped T=0 spin response and a separate VASP-derived Wannier edge model |

The operator-based extensions specify additional VASP outputs and their
producer requirements. The Bi PAW spin response uses matching licensed PAW
datasets and a serial instrumented VASP producer; its fresh SCF density is
separate from the historical Bi_Z2 fixture. Wannier models support optional
interpolation and ideal-edge calculations, alongside the direct WAVECAR examples.

VASP POTCAR files are not distributed. Follow each material's potential
specification and assemble them from a licensed library. See
[input acquisition](../INPUTS.md), [applying native commands to your system](../APPLY_TO_YOUR_SYSTEM.md),
and the [normalization migration guide](../../docs/MIGRATION.md) before
comparing historical Kubo magnitudes.
