# VASP files for the examples

The examples use monolayer MoS₂, a Bi bilayer and a magnetic MnBi₂Te₄ film. Choose a full k mesh for
Brillouin-zone maps and topological integrals, or a band path for k-resolved
spectra.

| Dataset | Files | Sampling and use |
|---|---|---|
| MoS₂ band path | [WAVECAR, POSCAR and EIGENVAL](1H-MoS2/KPATH/2.band/) | SOC; 48 points along K–Γ–K′; 32 bands. Kubo curvature, optical spectra and Γ wavefunctions |
| MoS₂ full-zone curvature | [VASP inputs](features/fukui-berry-curvature/inputs/) and [calculated output](features/fukui-berry-curvature/reference/) | 12×12 mesh; occupied bands 1–18. Fukui curvature map; see the [Fukui tutorial](features/fukui-berry-curvature/) for available inputs and reproduction |
| MoS₂ Kubo transport | [VASP preparation and Hall commands](features/kubo-hall/) | Full periodic meshes, occupations at specified chemical potentials and temperatures; separate mesh and intermediate-band checks |
| Historical Bi_Z2 full mesh | [WAVECAR](Bi_Z2/WAVECAR), [band energies](Bi_Z2/archive-2016-run/EIGENVAL) | SOC; 12×12 mesh; 18 bands. Occupied Chern number, Z₂ and insulating-gap charge-Hall response |
| Historical Bi_Z2 preparation | [SCF and NSCF inputs](Bi_Z2/inputs/) | Templates for a new calculation; archived SCF provenance is incomplete |
| Bi physical spin Hall and ideal edges | [Material guide](materials/bi-spin-hall/), [fresh SCF density and inputs](materials/bi-spin-hall/inputs/scf/), [PAW dataset guide](materials/bi-spin-hall/PSEUDOPOTENTIAL.md) | Fresh fixed-density SOC calculation; full PAW spin/velocity matrices for gapped T=0 transport, and a separately validated Wannier Hamiltonian for the ideal strip |
| MnBi₂Te₄, three septuple layers | [Structure and VASP inputs](materials/mnbi2te4-qah/inputs/), [fixed SCF density](materials/mnbi2te4-qah/inputs/scf/), [full Hamiltonian/position operators](materials/mnbi2te4-qah/inputs/wannier/operators/) | 21 atoms; SOC+U; full 6×6 optical mesh; supplied actual VASP-derived operators for VASPBERRY bands and dense Hall |

## Download the historical Bi_Z2 wavefunctions

The MoS₂ band-path files are included in the repository. Retrieve the larger
Bi WAVECAR with:

```bash
python3 examples/fetch_inputs.py bi --output-dir results/inputs/bi
```

Use `results/inputs/bi/WAVECAR` in the commands. A Git LFS checkout can instead
use the file under `examples/Bi_Z2/` after:

```bash
git lfs pull --include='examples/Bi_Z2/WAVECAR'
```

## Prepare the MoS₂ full mesh

The [Fukui tutorial](features/fukui-berry-curvature/) regenerates the 12×12 SOC
WAVECAR from the public MoS₂ charge density and the supplied NSCF inputs.
It requires a VASP installation and a matching licensed Mo/S POTCAR. The
149 MB WAVECAR is generated locally; its calculated curvature and figures
are included for comparison.

The [Kubo Hall tutorial](features/kubo-hall/) reuses this density and structure.
Its preparation helper generates larger meshes and band windows. Extra empty
states help converge the lower states that enter the Kubo sum: the stored
VASP `NBANDS` and optional `--pair-band-max` integration cutoff are recorded
separately. Reaching the total-energy tolerance alone does not establish the
accuracy of the highest empty states.

## Prepare the magnetic MnBi₂Te₄ film

The [material guide](materials/mnbi2te4-qah/) starts from an actual compressed
SCF charge density and a fixed bulk-derived slab. Generate the matching licensed
Mn/Bi/Te POTCAR and run the supplied NSCF preparation command. The optical
workflow retains ordinary VASP WAVECAR, WAVEDER, INCAR and OUTCAR files from
the same calculation. The supplied full Hamiltonian/position operators support
VASPBERRY's own dense Wannier calculation without rerunning the external
response solver. The guide distinguishes this model input from the direct
WAVECAR/WAVEDER routes and retains postw90 as an independent comparison.

## Prepare the Bi PAW spin and velocity inputs

The [Bi spin Hall example](materials/bi-spin-hall/) starts from its own fresh
SCF density and fixed two-atom geometry. `prepare_vasp.py` restores the split
compressed CHGCAR, checks its identity and prepares a full explicit mesh with
your matching licensed Bi POTCAR. Its `--stage scf` option prepares regeneration
of that density. These inputs are separate from the historical Bi_Z2 archive.

Use the [audited producer recipe](../tools/vasp544_spin_bridge/) to build a
serial complex noncollinear executable from the supported locally licensed
VASP source. `spin-export` converts its completed same-run output into physical
spin and velocity matrices; `spin-merge` assembles compatible fixed-charge
chunks without rotating their eigenstates. WAVECAR alone provides raw pseudo
spin matrices and does not supply the required PAW augmentation or full velocity.

The [general spin guide](../docs/SPIN_HALL.md) describes the insulating T=0
response, source-band closure checks and independent k-mesh refinement. A
VASP-derived Wannier Hamiltonian supplies the separate ideal-strip spectrum;
validate its bulk bands and strip width before interpreting the edge branches.

## Preparing new VASP data

The material guides describe the structures, pseudopotentials and VASP input
settings. Generate POTCAR locally using your licensed VASP library. For a
full-zone calculation, ensure that WAVECAR contains every point of the stated
mesh rather than only its symmetry-irreducible part.

The archived outputs support the numerical examples. Where the original
preceding SCF calculation is incomplete, the supplied preparation templates
specify a new calculation rather than an exact recreation of its historical
wavefunctions. See the [MoS₂](1H-MoS2/README.md) and [Bi](Bi_Z2/README.md)
material descriptions.
