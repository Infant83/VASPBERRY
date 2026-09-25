# MoS₂ stacking: band context and valley optical selection

This example extends the [monolayer circular-transition tutorial](../../features/circular-dichroism/)
to four small VASP systems. It demonstrates how the same VASPBERRY commands
apply to different layer registries and to degenerate band groups. VASP
supplies the band structure and wavefunctions; VASPBERRY evaluates the
polarization-resolved transitions. Ordinary WAVECAR is the starting point.
An optional comparison uses standard VASP WAVEDER and needs no source patch.

![VASP bands and valley optical selection for monolayer and bilayer MoS2](reference/figures/stacking-bands-selectivity.png)

[Calculated outputs, numerical checkpoints and channel spectra](reference/README.md)

The aligned bilayer is motivated by
[Yang et al., *Stacking-induced direct band gap in CVD-grown 1H MoS₂ bilayers*](https://doi.org/10.1038/s41586-026-11069-3).
The calculations here use explicitly documented, fixed idealized geometries
to illustrate software capabilities. The paper's actual slab input and
complete calculation settings were not available in the public source data.
This example does not reproduce its measured photoluminescence polarization
or infer an intervalley scattering rate from a Berry-curvature map.

## Systems and settings

| Case | Atoms | Registry | Mo-plane spacing | Occupied SOC states |
|---|---:|---|---:|---:|
| Monolayer | 3 | One trigonal-prismatic layer | — | 18 |
| 1H bilayer | 6 | Parallel layers, Mo on Mo and S on S | 6.74 Å | 36 |
| 2H bilayer | 6 | Antiparallel, inversion-symmetric control | 6.17 Å | 36 |
| 3R bilayer | 6 | Parallel layers translated by (1/3, 1/3) | 6.17 Å | 36 |

All cases use the same in-plane lattice, a = 3.1716343 Å, cell height
24 Å and intralayer Mo–S height 1.5968606 Å from the existing MoS₂ example.
Fractional translations refer to the supplied 60° real-space lattice basis.
The 1H and 3R separations use the experimental means reported in the paper's
supplement. The 2H control uses the same separation as 3R; it is not an
independently relaxed equilibrium structure. Atoms and cell vectors remain
fixed. There is no interlayer binding-energy or structural-stability claim.

The calculations use PBE+SOC, 400 eV, `PREC=Accurate`, `LASPH=.TRUE.` and
the documented [Mo/S PAW datasets](PSEUDOPOTENTIAL.md). Each case starts with
its own fresh 6×6 SCF density. The explicit path is Γ–M–K–Γ–K′, with 49
points. A separate full 6×6 mesh input recipe is supplied as an optional
extension; it is not part of the calculated reference results.
SCF and fixed-charge electronic tolerances are 10⁻⁷ and 10⁻⁸ eV.
Monolayer runs use 60 bands; bilayer SCFs use 96, and their path runs
use 64. The optional mesh templates use the same path-stage band counts.
This is a modest capability example, with no k-mesh, structural or
empty-band convergence claim.

The recorded calculations took about 59 minutes elapsed with up to four
single-thread VASP processes running concurrently. Individual VASP stages
took 13.5–30.3 minutes, and the largest measured per-process memory use was
1.85 GiB. These are timings from the reference machine, not performance
guarantees; the [source records](reference/) retain each completed stage.

## 1. Prepare and run ordinary VASP

Run these commands from the repository root. Install the postprocessing
dependencies and build native VASPBERRY:

```bash
python3 -m pip install -r requirements-transport.txt
make serial
```

Supply your licensed, concatenated Mo/S POTCAR. For one case:

```bash
python3 examples/materials/mos2-stacking-valley/prepare_vasp.py \
  --case monolayer --stage scf --potcar /path/to/Mo-S/POTCAR \
  --output-dir work/mos2-stacking/monolayer/scf

( cd work/mos2-stacking/monolayer/scf && /path/to/vasp_ncl > vasp.stdout )

python3 examples/materials/mos2-stacking-valley/prepare_vasp.py \
  --case monolayer --stage path --potcar /path/to/Mo-S/POTCAR \
  --scf-dir work/mos2-stacking/monolayer/scf \
  --output-dir work/mos2-stacking/monolayer/path

( cd work/mos2-stacking/monolayer/path && /path/to/vasp_ncl > vasp.stdout )
```

Repeat with `1h-bilayer`, `2h-bilayer` and `3r-bilayer`, using a separate
directory and fresh density for each structure. **For 3R, insert the
dipole-corrected SCF step below before preparing the path.** The helper checks the
completed SCF, matching structure and PAW datasets before preparing a
fixed-charge continuation. It never runs VASP itself or modifies its source.
Use your installation's supported MPI invocation if appropriate.

The 3R reference slab has a nonzero z dipole. It follows the standard VASP
recommendation to preconverge first, then enable the slab dipole correction
in a self-consistent warm restart. After the ordinary 3R SCF above:

```bash
python3 examples/materials/mos2-stacking-valley/prepare_vasp.py \
  --case 3r-bilayer --stage scf-dipole --potcar /path/to/Mo-S/POTCAR \
  --scf-dir work/mos2-stacking/3r-bilayer/scf \
  --output-dir work/mos2-stacking/3r-bilayer/scf-dipole

( cd work/mos2-stacking/3r-bilayer/scf-dipole && /path/to/vasp_ncl > vasp.stdout )

python3 examples/materials/mos2-stacking-valley/prepare_vasp.py \
  --case 3r-bilayer --stage path --potcar /path/to/Mo-S/POTCAR \
  --scf-dir work/mos2-stacking/3r-bilayer/scf-dipole \
  --output-dir work/mos2-stacking/3r-bilayer/path

( cd work/mos2-stacking/3r-bilayer/path && /path/to/vasp_ncl > vasp.stdout )
```

This restart uses `LDIPOL=.TRUE.`, `IDIPOL=3`, `DIPOL=0 0 0.5` and
`AMIN=0.01`, retaining the same nuclei and cell. Both subsequent 3R path
and mesh runs use the corrected density and settings. See VASP's
[LDIPOL](https://vasp.at/wiki/LDIPOL) and
[electrostatic-correction instructions](https://vasp.at/wiki/Electrostatic_corrections).
The other reference slabs have no detected z dipole in the density check.

The path/mesh templates include ordinary `LOPTICS=.TRUE.`,
`LPEAD=.FALSE.` and `LNABLA=.FALSE.` so the optional PAW comparison can be
performed from the same states. `--without-optics` prepares a WAVECAR-only
calculation. The present WAVEDER reader supports the documented VASP 5.4.4
longitudinal optical branch. See [PAW optical inputs](../../../docs/PAW_OPTICS.md)
for its exact contract.

For a full-zone check, repeat the path preparation with `--stage mesh
--mesh 6 6` and a new output directory. The mesh and path share the SCF
density, structure and electronic settings; their k-point lists differ.
A high-symmetry path alone cannot establish a global band gap or a BZ
integral. Increase `--mesh`, `--nbands` and `--path-points-per-segment` when
testing the corresponding sampling choices.

## 2. Evaluate transitions with VASPBERRY

The ordinary WAVECAR route is:

```bash
python3 examples/materials/mos2-stacking-valley/run.py \
  --case monolayer --run-dir work/mos2-stacking/monolayer/path \
  --output-dir results/mos2-stacking/monolayer
```

Add `--paw-optics` to evaluate the optional standard-WAVEDER comparison in
the same run. Use a fresh output directory when repeating a command.
The helper checks EIGENVAL against WAVECAR and writes the actual VASP path
bands, native circular spectra and a calculation record. It selects all
occupied initial states and a complete low-energy empty-state group: at
least two empty bands for monolayer and four for a bilayer. If those bands
touch further states, the upper boundary is expanded until it separates
the whole group by more than 2 meV at every point. `--final-band` may select
a larger complete group; at least one source band above it is retained as
a boundary check. Native and PAW calculations use the same band ranges.

This group sum is essential for the 2H inversion-symmetric control.
Polarization assigned to an arbitrarily chosen member of an unresolved
degenerate pair is not the observable of the complete subspace.

For another material, use the generic [native optical command](../../features/circular-dichroism/)
or [PAW optics command](../../../docs/PAW_OPTICS.md), with its actual filling
and band groups. This preparation helper intentionally retains the four
documented reference structures and PAW datasets.

## 3. Figures and interpretation

After calculating all four cases, draw their common comparison:

```bash
python3 examples/materials/mos2-stacking-valley/plot.py \
  --reference-dir results/mos2-stacking \
  --output-dir results/mos2-stacking-figures
```

`stacking-bands-selectivity` places each material's VASP bands above its K/K′
circular selectivity. Each band plot is referenced to that case's sampled
path VBM, not to a common absolute slab energy. Solid lines use native
WAVECAR momentum; dashed lines, when supplied, use PAW optical matrices.
Separate figures retain the full native and PAW channel spectra.
The figures focus on 1.3–2.2 eV; stored tables retain the complete 0.5–4 eV
photon grid.

The helicity convention is explicit: propagation along +z,
ε± = (x ± i y)/√2 and fields proportional to
Re[ε exp(i(q·r − ωt))]. Native `LEFT`/`RIGHT` use the same final-bra,
initial-ket contractions as `+`/`−`, up to a common normalization.
Both spectra use Gaussian σ = 0.05 eV. The plotted ratio is
η = (I₊ − I₋)/(I₊ + I₋); points below 0.1% of each k point's peak
total intensity are masked. Native output prints each intensity to four
decimal places. Its figure also requires the resulting conservative
uncertainty bound, 10⁻⁴/(I₊ + I₋), to be at most 0.005. This generic display
mask prevents a last-digit difference in weak channels from appearing as
resolved polarization. All original channels and ratios remain in the
reference tables; `native-eta-display.csv` records the bound and display
decision at every plotted K/K′ photon energy. No curve is forced to zero.

PAW transition strengths have units Å² and the broadened k-resolved
spectral density has units Å²/eV. No k weight, cell-volume prefactor or
absolute-absorption normalization is applied. Native intensities retain
their photon-energy/path normalization in arbitrary units. The η overlay
compares polarization selection; it is not an absolute intensity match
or a clean measure of the velocity correction under finite broadening.

These independent-particle transitions do not include excitons, scattering,
recombination or pumped populations. Their η is not the paper's measured
PL polarization. Time-reversal exchange of K/K′ helicity and cancellation
for the complete 2H degenerate subspace are the useful symmetry checks.
The earlier [Berry-curvature](../../features/fukui-berry-curvature/) and
[Hall examples](../../features/kubo-hall/) give the corresponding workflows;
their maps and convergence studies are not repeated here.

## Files and reuse

The `inputs/` tree contains POSCAR, INCAR and KPOINTS for each reference
stage, PAW identifiers and the geometry audit. POTCAR is supplied locally.
The [reference results](reference/) retain the actual VASP EIGENVAL files,
input snapshots, completion excerpts and normalized VASPBERRY outputs.
`bands.csv`/`.npz` preserve the VASP eigenvalues; `native-optical.csv`/`.npz`
preserve the native channels and η. Optional `paw-optical/` contains
transition and spectrum tables in CSV, DAT and NPZ with `optical.json`
units and operator metadata. Figures are PNG, PDF and SVG. Machine records
retain source associations and numerical settings separately from this guide.

To redraw the distributed numerical reference without running VASP:

```bash
python3 examples/materials/mos2-stacking-valley/plot.py \
  --output-dir results/mos2-stacking-reference-figures
```

Large WAVECAR and WAVEDER files are regenerated by the VASP preparation
steps above. The smaller archived outputs support figure reproduction;
they do not substitute for wavefunctions when recalculating transitions.

The [technical report](../../../docs/TECHNICAL_REPORT.md) uses this comparison
to illustrate VASPBERRY's optical capability across related input systems.
