# Band-resolved Berry curvature of monolayer MoS₂

This example calculates the two highest valence bands of 1H-MoS₂ along
K–Γ–K′ from the supplied VASP spinor wavefunctions. The band energies are even
under reversal of the path, while the Berry curvature near K and K′ has
opposite signs.

## Physical quantity

For a nondegenerate band, VASPBERRY evaluates

$$
\Omega_{n,z}(\mathbf{k})=-2\,\mathrm{Im}
\sum_{m\ne n}
\frac{D^x_{nm}(\mathbf{k})D^y_{mn}(\mathbf{k})}
{[E_n(\mathbf{k})-E_m(\mathbf{k})]^2},
\qquad D^a=\hbar v_a.
$$

The native WAVECAR implementation uses the canonical-momentum approximation
$v_a=p_a/m_e$. Its curvature is expressed in Å². PAW augmentation and
nonlocal/SOC velocity corrections are not included. The band-index sum uses
all intermediate states stored in WAVECAR.

## Input and calculation settings

The [MoS₂ dataset](../../1H-MoS2/KPATH/2.band/) provides the required VASP
output and the corresponding structure, sampling and calculation settings.

| Input or setting | Value |
|---|---|
| Calculation input | [WAVECAR](../../1H-MoS2/KPATH/2.band/WAVECAR) |
| Structure and settings | [POSCAR](../../1H-MoS2/KPATH/2.band/POSCAR), [INCAR](../../1H-MoS2/KPATH/2.band/INCAR) |
| Sampling | [KPOINTS](../../1H-MoS2/KPATH/2.band/KPOINTS): K–Γ and Γ–K′, 24 points per segment |
| VASP output | [EIGENVAL](../../1H-MoS2/KPATH/2.band/EIGENVAL), [OUTCAR](../../1H-MoS2/KPATH/2.band/OUTCAR) |
| Wavefunctions | SOC spinors; 32 bands; 400 eV plane-wave cutoff |
| Selected bands | 17 and 18, the two highest valence bands |
| Intermediate bands | 1–32 |

Γ occurs twice, at WAVECAR indices 24 and 25. The horizontal coordinate in the
figure is cumulative **Cartesian** distance along the path,
$s_j=\sum_{i<j}|\mathbf{k}_{i+1}-\mathbf{k}_i|$, in Å⁻¹. The reciprocal vectors
used in this conversion include $2\pi$.

## Run VASPBERRY

Start in the repository root and use a new output directory:

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

`-s 2` selects the two spinor components. `-kubo 2` evaluates the supplied
path; `-ii 17 -if 18` selects the bands to output. `-kubo_csv` produces the
full-precision point data used below. The other outputs are the combined
`BERRYCURV_KUBO.dat`, band-specific DAT files and the calculation log.

To run the calculation and generate both PNG and PDF figures:

```sh
python3 examples/features/kubo-curvature/run.py \
  --wavecar examples/1H-MoS2/KPATH/2.band/WAVECAR \
  --binary build/vaspberry-gfortran \
  --output-dir results/mos2-kubo
```

The plotting script requires NumPy and Matplotlib. For MPI, build with
`make mpi` and replace the executable in the native command with
`mpiexec -np 2 ../../build/vaspberry-mpi`. The two-process calculation
reproduces the serial point data for this input.

## Results

![Valence bands and Berry curvature of monolayer MoS2](reference/figure.png)

**Figure.** (a) The upper valence bands, with energies measured from the
maximum of band 18 on the sampled path, $E_{\mathrm{v}}=-1.274873881$ eV in
the original VASP energy reference. (b) Band-resolved Berry curvature in the
canonical-momentum approximation. Both panels use the same Cartesian path
distance, with K, Γ and K′ marked explicitly. Curvature is omitted where the
nearest other-band separation is at most $10^{-5}$ eV.

| Quantity | Band 17 | Band 18 |
|---|---:|---:|
| $\Omega_z(\mathrm{K})$ (Å²) | −5.489549 | −6.671251 |
| $\Omega_z(\mathrm{K}')$ (Å²) | +5.489548 | +6.671253 |
| Resolved points on the path | 44 of 48 | 44 of 48 |

[PNG figure](reference/figure.png) · [PDF figure](reference/figure.pdf) ·
[band and curvature table](reference/summary.csv) ·
[native Kubo output](reference/KUBO.csv)

The native CSV contains the band energy, fractional k coordinates, curvature
and nearest-band separation. `summary.csv` adds the Cartesian path distance
and marks unresolved curvature entries as blank. CSV energies retain the
original VASP reference; the energy shift is applied only in the figure.

## Interpretation and application to another system

The opposite signs at K and K′ are consistent with the time-reversal relation
between the two valleys. Near Γ, the selected states become nearly degenerate:
the smallest splitting is about $3.1\times10^{-9}$ eV. An individual-band Kubo
denominator is then unstable. Those points are omitted rather than interpreted
as a large physical curvature; degenerate states require a subspace treatment.

For another material, use its WAVECAR, choose the bands from the corresponding
VASP energies, and set the scalar/spinor option appropriately. Inspect band
separations before plotting individual-band curvature. Convergence of the
intermediate-state sum requires VASP calculations with increasing NBANDS.
The example script is configured for the supplied MoS₂ dataset; adapt the
native command and plotting selections for different systems.

A k path does not sample a Brillouin-zone area and cannot determine a Chern
number or Hall conductivity. The [Bi Hall example](../hall-valley/) illustrates
a full-mesh occupied-subspace calculation. The [material guide](../../1H-MoS2/README.md)
describes the supplied VASP dataset and preparation of new inputs.
