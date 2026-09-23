# Material datasets

Material folders retain their existing paths. Feature tutorials in
[`../features/`](../features/) link to these datasets where appropriate.

| Material / sampling | Available data | Reproduction boundary | Related feature |
|---|---|---|---|
| [1H-MoS₂, full 12 × 12 mesh](../1H-MoS2/) | Public charge density, NSCF inputs and freshly calculated curvature | Generate the full-mesh SOC WAVECAR with VASP, then run VASPBERRY | [Fukui curvature](../features/fukui-berry-curvature/) |
| [1H-MoS₂, K–Γ–K′ line](../1H-MoS2/KPATH/) | Line WAVECAR, EIGENVAL and plot data | A band-path example; no full-BZ Chern or Hall integral | [Kubo](../features/kubo-curvature/), [optical](../features/circular-dichroism/), [wavefunction](../features/wavefunction/) |
| [Bi buckled honeycomb, 12 × 12 mesh](../Bi_Z2/) | Input templates, wavefunctions, topological results and figures | New Fukui/Z₂/gap-Hall calculations from the Git LFS WAVECAR payload | [Fukui](../features/fukui-chern/), [Z₂](../features/z2/), [Hall](../features/hall-valley/) |
| [MnBi₂Te₄, three septuple layers](mnbi2te4-qah/) | Actual fixed SCF density, 21-atom structure, VASP preparation and numerical references | Generate the optical full mesh with licensed VASP; validate Wannier interpolation separately | Occupied Chern, PAW optical Hall and magnetic-band comparison |

The Bi `WAVECAR` is a Git LFS object of 200,421,600 bytes. A small text pointer
is not the wavefunction payload. If Git LFS is available, retrieve this one
file with `git lfs pull --include="examples/Bi_Z2/WAVECAR"`. The feature runner
checks for a pointer before attempting a WAVECAR calculation.

The archived Bi SCF provenance is incomplete; its package demonstrates
post-processing rather than an end-to-end reproduction of the original DFT
calculation. VASP `POTCAR` files are not distributed. Follow each material's
pseudopotential provenance notes and generate them using a licensed library.

For historical Kubo magnitudes, read the [normalization migration guide](../../docs/MIGRATION.md).
Reference examples never silently reinterpret old outputs as corrected data.

See [input acquisition](../INPUTS.md) and [applying the commands to your system](../APPLY_TO_YOUR_SYSTEM.md).
