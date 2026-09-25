# Insulating charge-Hall response of a Bi bilayer

This example calculates the zero-temperature sheet Hall conductivity of a
Bi bilayer from the supplied full-mesh VASP spinor wavefunctions. The chemical
potential lies in the band gap, with ten occupied bands at every k point.
Time-reversal symmetry gives a vanishing total charge-Hall response, although
the same occupied subspace has a nontrivial [Z₂ invariant](../z2/).

## Physical quantity and method

For a two-dimensional insulator,

```math
\sigma_{xy}=-\frac{e^2}{h}C,\qquad
C=\frac{1}{2\pi}\int_{\mathrm{BZ}}\mathrm{Tr}\,\Omega_z(\mathbf{k})\,d^2k.
```

The Fukui method evaluates the occupied-subspace flux from determinants of
neighboring wavefunction overlaps. Internal Kramers degeneracies are retained
within the occupied subspace; a gap is required between occupied and empty
states. This is appropriate for the supplied Bi data, whose individual
Kramers partners have numerical splittings of order $`10^{-7}`$ eV.
The wavefunction overlaps use the stored PAW pseudo-wavefunctions.

## Input and calculation settings

Obtain the supplied Bi WAVECAR from the repository root:

```sh
python3 examples/fetch_inputs.py bi --output-dir results/inputs/bi
```

An existing `examples/Bi_Z2/WAVECAR` populated by Git LFS can also be used.
The calculation requires the full binary WAVECAR, approximately 200 MB.

| Input or setting | Value |
|---|---|
| Structure and VASP records | [Bi material guide](../../Bi_Z2/README.md), [archived calculation](../../Bi_Z2/archive-2016-run/README.md) |
| Sampling | Full Γ-centered 12 × 12 × 1 mesh, 144 k points |
| Wavefunctions | 18 SOC spinor bands |
| Occupied subspace | Bands 1–10 |
| First empty state | Band 11 |
| Temperature | 0 K |
| Chemical potential | −1.30 to −0.90 eV in the original VASP energy reference, 41 points |

The supplied [SCF and NSCF templates](../../Bi_Z2/inputs/README.md) provide a
starting point for a new VASP calculation. This example uses the archived
WAVECAR directly.

## Run the calculation

From the repository root:

```sh
python3 tools/wavecar_fukui.py results/inputs/bi/WAVECAR \
  --nx 12 --ny 12 --spinor-components 2 \
  --energy-band 11 --map occupied=1:10 \
  --transport-full-t0 10 \
  --mu-min -1.3 --mu-max -0.9 --mu-num 41 \
  --valley-k 0.6666666666666666,0.3333333333333333,0 \
  --valley-kp 0.3333333333333333,0.6666666666666666,0 \
  --output-dir results/bi-hall-direct
```

`--transport-full-t0 10` includes cumulative occupied subspaces through band
10 and checks that band 11 remains empty. The mesh dimensions describe the
points already stored in WAVECAR. No factor of two is added for the SOC
spinors. Link, subspace-gap and mesh-coverage criteria are applied throughout
the requested chemical-potential interval.

The current command also produces regional columns using the specified K/K′
centers. These define a geometric partition; the quantity interpreted here is
the **total charge conductivity**, not a physical valley contrast.

To run the calculation and generate the reference plot:

```sh
python3 examples/features/hall-valley/run.py \
  --wavecar results/inputs/bi/WAVECAR \
  --output-dir results/bi-hall
```

This Python WAVECAR calculation requires NumPy and Matplotlib and produces
both PNG and PDF figures. A Fortran executable is not needed for this route.

## Results

![Insulating sheet Hall response of the Bi bilayer](reference/figure.png)

**Figure.** Sheet charge-Hall conductivity as a function of chemical potential
relative to the sampled valence-band maximum. The shaded gap extends from
$`E_{\mathrm{v}}`$ to $`E_{\mathrm{c}}`$; gray regions indicate the neighboring
valence (VB) and conduction (CB) band ranges. The blue curve shows only the
41 calculated chemical potentials. The gap is obtained from the supplied
12 × 12 mesh.

| Quantity | Result |
|---|---:|
| Sampled valence maximum $`E_{\mathrm{v}}`$ | −1.348963526 eV |
| Sampled conduction minimum $`E_{\mathrm{c}}`$ | −0.838919163 eV |
| Sampled indirect gap $`E_{\mathrm{g}}`$ | 0.510044362 eV |
| Minimum direct occupied–empty gap | 0.592448500 eV |
| Chemical potentials relative to $`E_{\mathrm{v}}`$ | 0.048964–0.448964 eV |
| Occupied bands throughout the interval | 10 |
| Total Chern integral | Approximately 0 |
| Maximum $`\lvert\sigma_{xy}\rvert`$ in the calculation | $`1.17\times10^{-16}\,e^2/h`$ |

[PNG figure](reference/figure.png) · [PDF figure](reference/figure.pdf) ·
[conductivity table](reference/summary.csv) ·
[complete transport output](reference/transport_full_t0.csv) ·
[occupied-subspace flux](reference/fukui_occupied.csv)

The CSV chemical potentials retain the original VASP energy reference.
Only the figure uses $`\mu-E_{\mathrm{v}}`$. The small nonzero numerical value is
retained in the tables; it represents numerical residual, not a finite
anomalous Hall signal. The output is a two-dimensional sheet response in
$`e^2/h`$, not a three-dimensional conductivity per unit length.

## Interpretation and application to another system

The vanishing charge-Hall plateau is consistent with time-reversal symmetry.
It does not determine the spin-Hall response, and the geometric regional
columns do not establish a valley-Hall effect. The nontrivial topology of this
Bi example is characterized separately by Z₂.

For another material, use a full uniform two-dimensional WAVECAR mesh and
choose the occupied subspace and chemical-potential interval from its band
energies. Keep at least one additional unoccupied band to verify the upper
boundary. A path or a symmetry-reduced irreducible mesh is insufficient for
this integral.

In a metallic interval the occupied count changes across the mesh. Such a
calculation requires the active-subspace checks and a separate k-mesh
convergence study. For a point-Kubo treatment, use suitable velocity-matrix
data and the [Kubo transport formulation](../../../docs/KUBO_TRANSPORT.md).
The present 12 × 12 example demonstrates the insulating calculation; it does
not establish mesh or PAW convergence for a new material.
