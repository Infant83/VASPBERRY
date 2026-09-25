# Calculated MoS₂ stacking references

These files are actual results from the four fixed-geometry PBE+SOC VASP
calculations described in the [reproduction guide](../README.md). The same
public `run.py` command produced the bands and native/optional PAW optical
tables. No model spectrum or symmetry averaging replaces the calculated data.

![VASP bands and circular optical selection for four MoS2 stackings](figures/stacking-bands-selectivity.png)

Top: direct VASP bands, referenced to each case's own path-sampled VBM.
Bottom: K/K′ circular selectivity with normal incidence. Solid lines use
ordinary WAVECAR momentum; dashed lines use standard WAVEDER optical
matrices from the same wavefunctions. The same band groups and Gaussian
σ = 0.05 eV are used. Weak or printing-precision-limited native ratios are
masked only in the figure, with the decision retained in
[the display table](figures/native-eta-display.csv).

## Numerical checkpoints

| Case | Initial → final bands | Path-sampled gap (eV) | PAW K peak (eV) | PAW η at K / K′ peak |
|---|---|---:|---:|---:|
| [Monolayer](monolayer/) | 1–18 → 19–20 | 1.67815 | 1.68 | −1 / +1 |
| [1H bilayer](1h-bilayer/) | 1–36 → 37–40 | 1.52841 | 1.68 | −1 / +1 |
| [2H bilayer](2h-bilayer/) | 1–36 → 37–40 | 1.26340 | 1.67 | +4.2×10⁻⁷ / +1.3×10⁻⁷ |
| [3R bilayer](3r-bilayer/) | 1–36 → 37–40 | 1.26563 | 1.68 | −0.99884 / +0.99884 |

The path-sampled gap is the lowest empty-state energy minus the highest
occupied-state energy on the 49-point path. It is not a full-zone gap
search. The idealized structures, including 1H, are software examples;
their band edges do not reproduce the paper's reported direct-gap result.
The reported η values refer to the strongest broadened K/K′ spectral peak
within the chosen band window, not the unbroadened sum over all transitions.
The [numerical summary](result-summary.json) includes both native and PAW
peak positions, ratios, direct K gaps and band-group boundary separations.

All four cases exchange helicity between K and K′. The complete 2H group
has negligible net selectivity, including both unresolved partners.
The [validation record](validation.json) gives the measured residuals,
independent WAVEDER contractions and CSV/DAT/NPZ comparisons. Native
four-decimal intensities can round a small contrast to zero or create an
unresolved last-digit difference; they are not a high-precision test of
the small PAW residuals.

Each `vasp/path-gmkgkp/path-tr-audit.json` also checks time-reversal-related
eigenvalues for the selected optical group and for every band visible in
the figure. The largest visible-band residual is below 1.3×10⁻⁸ eV.
This is a consistency check of the plotted states, not an empty-band
convergence claim for all stored high-energy states.

## Channel spectra and files

- [Native channel spectra](figures/stacking-native-spectra.png): retained
  photon/path normalization in arbitrary units.
- [PAW channel spectra](figures/stacking-paw-spectra.png): k-resolved
  length-matrix spectral density in Å²/eV.
- [Figure metadata](figures/plot.json): source associations, units and
  display masks. Each figure is supplied as PNG, PDF and SVG.

Each case directory contains:

| Files | Meaning |
|---|---|
| `bands.csv`, `bands.npz` | Original VASP eigenvalues, occupations and symmetry-path coordinates |
| `native-optical.csv`, `native-optical.npz` | Native circular channels and intensity-masked η |
| `paw-optical/transitions.*`, `paw-optical/spectra.*` | Optional PAW transition strengths and broadened spectra, in CSV/DAT/NPZ |
| `result.json`, `paw-optical/optical.json` | Calculation settings, units, source identities and commands |
| `vasp/` | Exact completed-stage INCAR, POSCAR, KPOINTS, EIGENVAL, OSZICAR and selected OUTCAR completion records |

The 3R records include both the initial SCF and the self-consistent dipole
restart used to generate its path density. Full-zone response meshes were
not calculated for this example. Licensed POTCAR and large WAVECAR/WAVEDER
files are regenerated with the guide's commands.

To redraw all three figures from these distributed data:

```bash
python3 examples/materials/mos2-stacking-valley/plot.py \
  --output-dir results/mos2-stacking-reference-figures
```

The figure covers the near-edge 1.3–2.2 eV interval; every spectrum table
retains the complete calculated 0.5–4 eV range. These are k-resolved
independent-particle spectra, without excitons or photoluminescence dynamics.
