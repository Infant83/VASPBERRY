# VASP files for the examples

The examples use monolayer MoS₂ and a Bi bilayer. Choose a full k mesh for
Brillouin-zone maps and topological integrals, or a band path for k-resolved
spectra.

| Dataset | Files | Sampling and use |
|---|---|---|
| MoS₂ band path | [WAVECAR, POSCAR and EIGENVAL](1H-MoS2/KPATH/2.band/) | SOC; 48 points along K–Γ–K′; 32 bands. Kubo curvature, optical spectra and Γ wavefunctions |
| MoS₂ full-zone curvature | [VASP inputs](features/fukui-berry-curvature/inputs/) and [calculated output](features/fukui-berry-curvature/reference/) | 12×12 mesh; occupied bands 1–18. Fukui curvature map; see the [Fukui tutorial](features/fukui-berry-curvature/) for available inputs and reproduction |
| Bi full mesh | [WAVECAR](Bi_Z2/WAVECAR), [band energies](Bi_Z2/archive-2016-run/EIGENVAL) | SOC; 12×12 mesh; 18 bands. Occupied Chern number, Z₂ and insulating-gap Hall response |
| Bi VASP preparation | [SCF and NSCF inputs](Bi_Z2/inputs/) | Templates for a new calculation |

## Download the Bi wavefunctions

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
