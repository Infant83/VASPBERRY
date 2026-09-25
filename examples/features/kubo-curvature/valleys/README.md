# Single-band curvature near the MoS₂ valleys

The upper valence band of nonmagnetic SOC MoS₂ touches its partner at Γ,
so it is not an isolated band over the whole Brillouin zone. Near K and K′,
however, the spin–orbit splitting makes its **individual-band Berry curvature**
well defined. This example calculates band 18 on two small Cartesian patches
and checks its separation from every other stored band.

![Band-18 curvature and valence bands near K and Kprime](reference/figure.png)

The upper panels show Ω₁₈,z around K and K′. Each dashed line is the
δky = 0 cut shown below, alongside bands 17 and 18. The color maps and cuts
use the same NumPy bilinear interpolation of the actual **9 × 9 samples per
valley**. Symbols on the curvature lines mark those calculated points.
Interpolation changes only the display; it adds no new electronic-structure
information. Both patches extend ±0.12 Å⁻¹ from the valley center in kx and ky.

## Inputs and calculation

Start with the completed [full-mesh MoS₂ VASP example](../../fukui-berry-curvature/README.md).
The valley calculation retains its structure, fixed SCF charge density,
matching PAW datasets, **400 eV cutoff and 26 spinor bands**. Only KPOINTS
changes. Run these commands from the repository root:

```bash
python3 examples/features/kubo-curvature/valleys/prepare.py \
  --input-dir results/mos2-fullmesh-vasp \
  --output-dir results/mos2-valleys-vasp

cd results/mos2-valleys-vasp
/path/to/vasp_ncl > vasp.stdout.log 2> vasp.stderr.log
cd ../..

make serial
python3 examples/features/kubo-curvature/valleys/run.py \
  --wavecar results/mos2-valleys-vasp/WAVECAR \
  --output-dir results/mos2-valleys
```

The native command used by the runner is:

```bash
repo_dir="$PWD"
mkdir -p results/mos2-valleys-native
(
  cd results/mos2-valleys-native
  "$repo_dir/build/vaspberry" \
    --wavecar "$repo_dir/results/mos2-valleys-vasp/WAVECAR" \
    --spinor 2 --task kubo --bands 17:18 \
    --curvature-csv KUBO.csv --output BERRYCURV > vaspberry.log
)
```

It exports bands 17 and 18 with intermediate bands 1–26; the map shows
band 18. The runner compares CSV energies, coordinates and minimum gaps
against an independent WAVECAR read before plotting. It requires the target
band to remain separated by more than 10 meV at every sampled point.

The reference VASP step took **161 seconds**, with **653 MB** peak resident
memory on the reference machine. Timings are illustrative. The WAVECAR and
licensed PAW files are not bundled; the preparation command generates the
162-point KPOINTS file from your existing VASP inputs.

## Reference results

| Quantity | Result |
|---|---:|
| Isolated band-18 samples | 162 / 162 |
| Minimum separation from any other stored band | 0.130203 eV |
| Band 17–18 splitting at K | 0.147110 eV |
| Ω₁₈,z(K) | −6.414823 Å² |
| Ω₁₈,z(K′) | +6.414836 Å² |
| Maximum time-reversal pair residual | 0.000149 Å² |

The opposite signs describe valley-resolved anomalous motion. The patch
results are not a full-BZ Chern number or conductivity. They use the current
canonical-momentum Kubo approximation, which omits PAW augmentation and
nonlocal/SOC velocity corrections; the 26-band window is an example setting,
not a demonstration of empty-band convergence.

Reference files:

- [Native KUBO.csv](reference/KUBO.csv), [all VASP bands](reference/bands.csv)
  and [band-18 samples with isolation gaps](reference/summary.csv).
- [Plotted line values](reference/line.csv), [PDF figure](reference/figure.pdf)
  and [calculation checks](reference/result.json).
- VASP [OUTCAR](reference/vasp/OUTCAR), [EIGENVAL](reference/vasp/EIGENVAL)
  and [OSZICAR](reference/vasp/OSZICAR).

To redraw the supplied results without rerunning VASP:

```bash
python3 examples/features/kubo-curvature/valleys/plot.py \
  --output-dir results/mos2-valleys-redrawn
```

For another material, select a physically isolated band and verify its gap
throughout the region of interest. Use an occupied-band bundle when bands
touch inside the selected subspace, as in the [full-BZ Kubo example](../README.md).
