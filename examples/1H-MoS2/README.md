# Monolayer 1H-MoS₂

These VASP data support complementary reciprocal-space and real-space
calculations for nonmagnetic monolayer MoS₂ with spin–orbit coupling.

| Dataset | Sampling and supplied files | Tutorial |
|---|---|---|
| Full Brillouin zone | Γ-centered 12 × 12 × 1 mesh; public SCF density and executable NSCF preparation, current native result and historical `BERRYCURV.dat`. Generate the 149 MB mesh WAVECAR using licensed VASP. | [Fukui Berry-curvature map](../features/fukui-berry-curvature/) |
| K–Γ–K′ path | 48 k points, 32 spinor bands; actual `KPATH/2.band/WAVECAR`, VASP energies and calculation records. | [Kubo curvature](../features/kubo-curvature/), [circular optical response](../features/circular-dichroism/) |
| Γ-point state | Records 24 and 25 of the path WAVECAR are Γ. | [Real-space wavefunction](../features/wavefunction/) |

## Fukui curvature over the first Brillouin zone

The occupied 1–18 subspace has opposite-sign Berry curvature near K and K′,
while the total Chern number vanishes. The [Fukui tutorial](../features/fukui-berry-curvature/)
provides the complete public-density NSCF preparation and current native
reference: 400 eV cutoff, 26 SOC bands, occupied bands 1–18, curvature about
±12.31 Å² and C = 0. Its [input guide](../features/fukui-berry-curvature/inputs/)
shows how to generate the required full-mesh WAVECAR.

The original full-mesh `BERRYCURV.dat` remains a separate historical
comparison. Its 144 independent plaquette values are repeated over a
3 × 3 reciprocal display region.

Plot the existing result in Cartesian reciprocal coordinates:

```bash
python3 tools/plot_berry_curvature.py \
  --input examples/1H-MoS2/BERRYCURV.dat \
  --poscar examples/1H-MoS2/POSCAR \
  --output results/mos2-berry-curvature.png \
  --title '1H-MoS2'
```

This command redraws the historical result. For a new calculation, generate
a full-mesh SOC WAVECAR and use the [step-by-step Fukui tutorial](../features/fukui-berry-curvature/).
The path WAVECAR cannot be substituted for a two-dimensional mesh.

## VASP preparation

The top-level `INCAR`, `KPOINTS` and `POSCAR` describe the full-mesh setup.
`KPATH/1.scf/` contains the charge-density calculation and compressed CHGCAR;
`KPATH/2.band/` contains the subsequent line calculation. `KPATH/3.BC_kubo/`
holds its historical band-resolved curvature output.

Construct `POTCAR` from a licensed PAW-PBE library in the order Mo, S, following
[`PSEUDOPOTENTIAL.md`](PSEUDOPOTENTIAL.md). To generate a complete mesh for
VASPBERRY, disable symmetry reduction with `ISYM=-1` and retain wavefunctions
with `LWAVE=.TRUE.`. The scheduler example in `KPATH/1.scf/sbatch_vasp.sh` takes
the executable from `VASP_BIN`; adjust resources for your system.
