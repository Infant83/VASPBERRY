# Material datasets

Material folders retain their existing paths. Feature tutorials in
[`../features/`](../features/) link to these datasets where appropriate.

| Material / sampling | Available data | Reproduction boundary | Related feature |
|---|---|---|---|
| [1H-MoS₂, full 12 × 12 mesh](../1H-MoS2/) | Public charge density, NSCF inputs and freshly calculated curvature | Generate the full-mesh SOC WAVECAR with VASP, then run VASPBERRY | [Fukui curvature](../features/fukui-berry-curvature/) |
| [1H-MoS₂, K–Γ–K′ line](../1H-MoS2/KPATH/) | Line WAVECAR, EIGENVAL and plot data | A band-path example; no full-BZ Chern or Hall integral | [Kubo](../features/kubo-curvature/), [optical](../features/circular-dichroism/), [wavefunction](../features/wavefunction/) |
| [MoS₂ monolayer and 1H/2H/3R bilayers](mos2-stacking-valley/) | Four fixed structures, actual VASP inputs/band energies, native and PAW optical references | Regenerate WAVECAR/WAVEDER with ordinary VASP; supplied tables reproduce all figures | Complete-group valley optical selection and stacking symmetry |
| [Bi buckled honeycomb, 12 × 12 mesh](../Bi_Z2/) | Input templates, wavefunctions, topological results and figures | New Fukui/Z₂/gap-Hall calculations from the Git LFS WAVECAR payload | [Fukui](../features/fukui-chern/), [Z₂](../features/z2/), [Hall](../features/hall-valley/) |
| [MnBi₂Te₄, three septuple layers](mnbi2te4-qah/) | Actual fixed SCF density, 21-atom structure, VASP preparation and full Wannier operators | VASP for wavefunctions/optical files; supplied operators for Wannier band interpolation and dense Hall integration | Occupied Chern, PAW optical diagnostic, full-connection QAH and independent postw90 comparison |
| [Bi PAW spin Hall and ideal edges](bi-spin-hall/) | Fresh SCF density, two-atom structure, physical spin/velocity preparation and a VASP-derived Wannier model | Matching licensed PAW data and the serial instrumented VASP producer for new operators; VASPBERRY for the gapped T=0 response and ideal strip | [Spin Hall, spin matrices and edge spectra](../../docs/SPIN_HALL.md) |

The historical Bi_Z2 `WAVECAR` is a Git LFS object of 200,421,600 bytes. A small text pointer
is not the wavefunction payload. If Git LFS is available, retrieve this one
file with `git lfs pull --include="examples/Bi_Z2/WAVECAR"`. The feature runner
checks for a pointer before attempting a WAVECAR calculation.

The archived Bi_Z2 SCF provenance is incomplete; its package demonstrates
post-processing rather than an end-to-end reproduction of the original DFT
calculation. The separate `bi-spin-hall` folder supplies a newly generated
density with its SCF inputs. It does not reuse the archived density as a
production spin Hall input. VASP `POTCAR` files are not distributed. Follow each material's
pseudopotential provenance notes and generate them using a licensed library.

For historical Kubo magnitudes, read the [normalization migration guide](../../docs/MIGRATION.md).
Reference examples never silently reinterpret old outputs as corrected data.

See [input acquisition](../INPUTS.md) and [applying the commands to your system](../APPLY_TO_YOUR_SYSTEM.md).
