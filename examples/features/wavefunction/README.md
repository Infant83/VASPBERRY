# A Γ-point state of monolayer MoS₂ in real space

VASPBERRY writes the real and imaginary amplitudes of a selected VASP state.
This example uses occupied band 18 at Γ and combines both spinor components
to obtain its pseudo-wavefunction density.

![Cartesian real-space density of a Gamma state in MoS2](reference/figure.png)

**Figure 1.** (a) State density integrated along z, shown in the actual oblique
unit cell with Cartesian x and y axes. Colors interpolate between the sampled
grid values. (b) The corresponding density averaged over the in-plane cell.
Both spinor components are included. [PDF](reference/figure.pdf) ·
[Plane-averaged numerical data](reference/summary.csv).

## 1. Calculate from VASP wavefunctions

Use `WAVECAR`, `POSCAR` and `EIGENVAL` from the supplied
[MoS₂ band-path calculation](../../1H-MoS2/KPATH/2.band/). Its 48-point
K–Γ–K′ path contains Γ at indices 24 and 25; this example selects index 24.
The wavefunctions include SOC and use a 400 eV plane-wave cutoff.

```bash
make serial
python3 -m pip install -r requirements-transport.txt
repo_dir="$PWD"
mkdir -p results/mos2-wavefunction
cp examples/1H-MoS2/KPATH/2.band/POSCAR results/mos2-wavefunction/
cp examples/1H-MoS2/KPATH/2.band/EIGENVAL results/mos2-wavefunction/
(
  cd results/mos2-wavefunction
  "$repo_dir/build/vaspberry" \
    --wavecar "$repo_dir/examples/1H-MoS2/KPATH/2.band/WAVECAR" \
    --task wavefunction --spinor 2 --mesh 48,1 --wavefunction-band 18 --kpoint 24 --real-grid 24,24,64 --imaginary 1 > stdout.log
)
```

## 2. Postprocess the native output

```bash
python3 examples/features/wavefunction/run.py \
  --output-dir results/mos2-wavefunction --postprocess-only
```

`--wavefunction-band 18` selects the band and `--kpoint 24` selects Γ. `--real-grid 24,24,64` sets the
real-space sampling along the cell vectors, and `--imaginary 1` includes imaginary
amplitudes. The postprocessor reads these amplitude grids, combines the two
spinor components and checks the result against the stored WAVECAR coefficients.
`--postprocess-only` does not rerun VASPBERRY.

## Optional reproduction helper

To calculate, check and plot in one step in a fresh directory:

```bash
python3 examples/features/wavefunction/run.py --output-dir results/mos2-wavefunction-auto
```

## Reading the result

`PARCHG-W-K024-E018-SPIN1` contains the real amplitudes and
`PARCHG-W-K024-E018-IM-SPIN1` the imaginary amplitudes. Each contains two
successive grids, one per spinor component. Despite the CHGCAR-like file
format, these are amplitudes rather than charge densities.

For the native output W(r) and cell volume V,

```math
\psi(\mathbf r)=W(\mathbf r)/V,\qquad
\rho(\mathbf r)=|\psi_\uparrow|^2+|\psi_\downarrow|^2.
```

The integrated pseudo-density is **0.825814**, consistent with the norm of
the stored plane-wave coefficients. It excludes the PAW augmentation and
therefore need not integrate to one. No renormalization is applied in the
figure. The reconstructed pseudo-wavefunction ψ agrees with a direct Fourier
sum of the stored coefficients to approximately 3.4 × 10⁻⁷ Å⁻³ᐟ².

The plotting step writes a z profile in `summary.csv` and figures in PNG
and PDF. The [original amplitude grids](reference/raw-output.tar.gz) can
also be viewed with a volumetric-data viewer.

## Other materials

Keep matching `WAVECAR`, `POSCAR` and `EIGENVAL` files together. Select the
band, actual Γ index and a real-space grid appropriate to the cutoff. The
helper accepts `--input-dir`, `--band` and `--k-index` for compatible SOC Γ
states in an xy-plane slab; other geometries can use the native amplitude
files with a suitable viewer.

A single state's phase, or its basis within a degenerate subspace, can change
between calculations. Sum densities over the entire subspace when comparing
a degenerate manifold. The [synthetic Fourier reference](../../../validation/models/wavefunction/)
illustrates the amplitude convention independently.
