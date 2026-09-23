# MoS₂ Hall reference results

These are outputs from actual VASP SOC calculations using the public MoS₂ structure
and SCF density. The [workflow](../README.md) gives the input preparation and
VASPBERRY commands. The `reference/calculation` and `reference/plot` folders are
the quick 12 × 12 / 26-band catalog baseline; `cases/` contains the separate
padded-source studies below.

## Calculation conditions

- Fixed density, `ICHARG=11`, `ISYM=-1`, 400 eV; 18 occupied spinor bands.
- Full Γ-centered 12 × 12, 24 × 24 and 36 × 36 meshes.
- Chemical potentials: 61 values, μ − Ev from −0.20 to +0.10 eV; 0 and 300 K.
- Periodic Cartesian disks of radius 0.35 Å⁻¹ around K and K′, with `rest` completing the BZ.
- Δσ(μ) = σ(μ) − σ(midgap); valley difference = K − K′, without dividing by two.
- The displayed path is an independent 26-band K–Γ–K′ calculation with the same structure, cutoff and fixed density. It does not supply the integration samples.

The main figure uses 36 × 36, 60 stored source bands and a retained pair cutoff
M = 40. The maximum absolute total charge response is 2.24e-09 e²/h;
the maximum absolute total Δσ is 1.8e-09 e²/h. Time-reversal cancellation
is measured from the raw calculation, without forcing the result to zero.

## Separate convergence checks

The metric compares the 300 K ΔσK − ΔσK′ curve on the common μ − Ev grid.
Relative L2 means the norm of the difference divided by the norm of the finer
curve. These checks address separate numerical choices; they do not establish
joint k-mesh and pair-cutoff convergence for every combination.

**k sampling:** 60 stored bands and M = 40 are fixed.

| Refinement | max absolute change (e²/h) | relative L2 | below 1% |
|---|---:|---:|:---:|
| 12 × 12 → 24 × 24 | 0.0211558 | 12.8897% | no |
| 24 × 24 → 36 × 36 | 0.0031028 | 1.4824% | no |

The last 24 × 24 → 36 × 36 mesh change is 1.482%, above the 1% criterion. The mesh study therefore remains unconverged by this criterion; the measured residual is retained explicitly.

**Pair cutoff:** the same 12 × 12, 96-band WAVECAR is used for every M.

| Refinement | max absolute change (e²/h) | relative L2 | below 1% |
|---|---:|---:|:---:|
| 26 → 40 | 0.0102059 | 4.2702% | no |
| 40 → 60 | 0.00854596 | 3.3388% | no |
| 60 → 80 | 0.00286355 | 1.2626% | no |
| 80 → 90 | 0.00051161 | 0.2024% | yes |

The last 80 → 90 cutoff refinement passes the stated 1% comparison criterion.
M = 40 remains a deliberately fixed cutoff for the mesh study; its amplitude
is not the result of the largest-cutoff calculation. The 0 K steps are retained
because occupations change discretely on a finite mesh. No artificial smoothing
is applied to the transport curves.

A further eigenstate check changes the stored source from 60 to 96 bands at
fixed 12 × 12 and M = 40: the curve changes by at most
3.66e-09 e²/h (relative L2 2.26e-06%).
This separates convergence of the retained eigenstates from convergence of the
intermediate-state sum.

## Numerical sanity and approximation

A loose occupied-energy stopping criterion left the highest empty eigenstates
poorly converged. Retaining all states in such files produced a spurious total
charge response, despite agreement between the native bundle and pair routines.
The published studies retain an accurately converged lower window within a larger
source. Every final VASP stage reached `EDIFF` and ended normally; checkpoint
stages are distinguished in the source conditions.

Changing the explicit numerical-degeneracy grouping tolerance from 10⁻⁷ to
10⁻⁶ eV changes absolute σ by at most 6.98e-12 e²/h and Δσ by at most
2.12e-13 e²/h in the checked 36 × 36 / M = 40 and 12 × 12 / M = 90 cases.
The recorded maximum occupation shift is 1.06e-08. CSV, DAT and NPZ contain
identical numerical values.

The native operator is canonical momentum without PAW augmentation or nonlocal/SOC
velocity corrections. Numerical refinement does not remove that approximation.
The μ scan is a rigid-band occupation scan; regional charge responses depend on
the stated geometric partition and are not a conserved valley-current calculation.

## Measured calculation cost

The calculations used one CPU thread per VASP worker on a 64 GiB host, with
eight concurrent workers as the shared scheduling target. Timings include
simultaneous calculations on the same machine. The table reports the latest
accepted warm NSCF stage; earlier seed calculations and checkpoint attempts are
retained in the source conditions. Worker-minutes sum elapsed time over all
workers, so the six 36 × 36 subsets can run concurrently within the resource cap.

| Mesh | stored bands | final NSCF worker-minutes | peak RSS per worker (GiB) |
|---|---:|---:|---:|
| 12 × 12 | 60 | 5.61 | 0.92 |
| 24 × 24 | 60 | 39.10 | 2.99 |
| 36 × 36 | 60 | 161.05 | 1.27 |
| 12 × 12 | 96 | 32.12 | 1.30 |

## Files and reproduction

Each `cases/` folder contains `conductivity.csv`, `.dat`, `.npz`, metadata,
`results.json` and `conditions.json`. Absolute σ, Δσ and carrier counts are included.
Shared `vasp/` folders contain the actual `INCAR`, `POSCAR`, `KPOINTS` and compressed
`OUTCAR`, `EIGENVAL`, `OSZICAR`, with the warm-restart iteration recipe. The 36 × 36
source was assembled from six independent 216-point fixed-density subsets;
all wavefunction records were checked byte-for-byte. Occupations and integration
weights are recomputed on the complete mesh, so chunk Fermi levels are not used.
Licensed `POTCAR`, large `WAVECAR` and raw interband-pair exports are regenerated
by the workflow rather than bundled here.

Regenerate both composite figures from the numerical tables:

```bash
python3 examples/features/kubo-hall/plot.py \
  --case examples/features/kubo-hall/reference/cases/12x12-source60-cap40 \
  --case examples/features/kubo-hall/reference/cases/12x12-source96-cap26 \
  --case examples/features/kubo-hall/reference/cases/12x12-source96-cap40 \
  --case examples/features/kubo-hall/reference/cases/12x12-source96-cap60 \
  --case examples/features/kubo-hall/reference/cases/12x12-source96-cap80 \
  --case examples/features/kubo-hall/reference/cases/12x12-source96-cap90 \
  --case examples/features/kubo-hall/reference/cases/24x24-source60-cap40 \
  --case examples/features/kubo-hall/reference/cases/36x36-source60-cap40 \
  --primary examples/features/kubo-hall/reference/cases/36x36-source60-cap40 \
  --output-dir results/mos2-hall-reference-figures
```

[Main figure (PDF)](figures/mos2-hall.pdf) ·
[Cutoff and temperature checks (PDF)](figures/mos2-hall-checks.pdf)
