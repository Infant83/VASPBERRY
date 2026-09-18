# Kubo Berry curvature: VASP MoS₂ WAVECAR → VASPBERRY → figure

This tutorial runs the current **Fortran VASPBERRY executable on an actual
VASP SOC WAVECAR**, then plots the two upper valence bands along K–Γ–K′.
The expected result is opposite-sign curvature near K and K′. Γ and nearby
unresolved individual-band states are explicitly omitted from the curvature
plot. The calculation is a line scan, so it does not produce a Chern number or
Hall conductivity.

## 1. Input files

The [public MoS₂ calculation](../../1H-MoS2/KPATH/2.band/) contains:

| File | Role in this tutorial |
|---|---|
| [WAVECAR](../../1H-MoS2/KPATH/2.band/WAVECAR) | Actual input to VASPBERRY: 48 k points, 32 spinor bands, 400 eV cutoff |
| [KPOINTS](../../1H-MoS2/KPATH/2.band/KPOINTS) | Two 24-point line segments; records 24 and 25 both represent Γ |
| [POSCAR](../../1H-MoS2/KPATH/2.band/POSCAR), [INCAR](../../1H-MoS2/KPATH/2.band/INCAR) | Original structure and VASP settings for interpreting the supplied data |
| [EIGENVAL](../../1H-MoS2/KPATH/2.band/EIGENVAL), [OUTCAR](../../1H-MoS2/KPATH/2.band/OUTCAR) | VASP energies and calculation record |

The WAVECAR is **60,521,760 bytes**, SHA256
`33f8546512856d6c04ad0a80454b18ac9b60e2af4b98f2f85b376ec49b1a8d9f`.
No VASP run or JSON model input is required to reproduce the VASPBERRY result.
POTCAR is not distributed; see the [material provenance](../../1H-MoS2/README.md)
before attempting to regenerate the VASP states themselves.

All commands below start in the repository root.

## 2. Run VASPBERRY directly

```sh
make serial
mkdir -p results/mos2-kubo-direct
cd results/mos2-kubo-direct
../../build/vaspberry-gfortran \
  -f ../../examples/1H-MoS2/KPATH/2.band/WAVECAR \
  -s 2 -kubo 2 -ii 17 -if 18 \
  -kubo_csv KUBO.csv -o BERRYCURV > vaspberry.log
cd ../..
```

`-s 2` reads the two components of each SOC spinor. `-ii 17 -if 18` selects the
bands whose curvature is printed; the Kubo sum uses **all 32 bands present in
WAVECAR** as intermediate states. `-kubo 2` selects line mode and avoids a BZ
integral. `-kubo_csv` writes physical, full-precision per-band point data.
Choose a fresh output directory and CSV filename for every run.

The main outputs are `KUBO.csv`, `BERRYCURV_KUBO.dat`, the band-specific DAT
files, and the log. CSV columns include spin, k index, band, fractional k,
energy in eV, `omega_z_A2`, and the closest other-band gap in eV.
The DAT file's generic `K-GRID`/`dk²` header retains the default grid parameters
in line mode; it is not evidence of a full mesh. The explicit
`integration=NONE_K_PATH` and `Chern_number=NOT_APPLICABLE_K_PATH` records define
this calculation's scope.

## 3. Reproduce the reference figure and checks

```sh
python3 examples/features/kubo-curvature/run.py \
  --wavecar examples/1H-MoS2/KPATH/2.band/WAVECAR \
  --binary build/vaspberry-gfortran \
  --output-dir results/mos2-kubo
```

This command **reruns VASPBERRY**, checks the CSV against the WAVECAR energies
and nearest gaps, then writes the figure. It never reads a model Hamiltonian.
Dependencies are NumPy and Matplotlib. Output directories are never overwritten.

| Reference | Expected result / meaning |
|---|---|
| [figure.png](reference/figure.png) | VASP valence energies and VASPBERRY Kubo curvature |
| [summary.csv](reference/summary.csv) | All 96 band/k rows; blank curvature means the individual band is unresolved |
| [KUBO.csv](reference/KUBO.csv) | Actual unmodified native per-band output, including diagnostic near-degenerate values |
| [result.json](reference/result.json) | Numerical checks, units, caveats, and output checksums |
| [provenance.json](reference/provenance.json) | WAVECAR/source/binary hashes and executed commands |

![MoS₂ native Kubo result](reference/figure.png)

For band 18, Ωz is approximately **−6.671251 Å² at K** and **+6.671253 Å² at K′**.
Each plotted band has **44 valid points out of 48** using a `10⁻⁵ eV` isolation
threshold. The source's near-Γ splitting reaches about `3.1 × 10⁻⁹ eV`;
dividing by that tiny gap gives unstable individual-band raw values. The
tutorial preserves them in the native CSV but masks them in the figure and
summary. It checks valid curvature against the stored reference with an
absolute tolerance of `10⁻⁵ Å²`, allowing compiler roundoff.

## 4. Apply the workflow to your VASP calculation

1. Supply your own SOC WAVECAR and its matching KPOINTS/structure records.
2. Choose the actual band numbers from your VASP energies. Keep `-s 2` for SOC
   spinors; scalar or collinear calculations require their appropriate setting.
3. Use `-kubo 2` for paths and inspect `min_gap_eV` before interpreting a band.
   Degenerate groups require a subspace treatment; do not reduce a safety
   threshold merely to retain every plotted point.
4. Plot the CSV using its actual k coordinates and units. Increase the number
   of intermediate bands through a new VASP calculation with larger NBANDS,
   and check convergence for your observable.

The convenience runner verifies the exact tutorial WAVECAR checksum. For
another system, adapt the explicit VASPBERRY command and the small plotting
script rather than expecting the MoS₂ reference checks to apply unchanged.

**Method scope:** this native implementation uses canonical momentum. PAW
augmentation and nonlocal/SOC velocity corrections are absent, so the numbers
are a reproducible bare-momentum approximation. They are already in the
standard `−2 Im` normalization; do not divide them by two. A path has no area
weights for a BZ integral. For an actual full-mesh occupied-subspace charge
response, continue with the [Bi Hall tutorial](../hall-valley/).

The QWZ implementation oracle now lives under
[developer model validation](../../../validation/models/kubo-curvature/).

For maintainers, package a successful replay into a **new** reference directory
with `python3 examples/features/kubo-curvature/export_reference.py --run-dir
results/mos2-kubo --output-dir results/mos2-kubo-reference`. Distributed hashes
cover the shipped files; the separately labeled full-run hashes retain the
original native-output and log record.
