# Berry curvature and response functions from VASP wavefunctions

[Download the PDF report](TECHNICAL_REPORT.pdf).

## Abstract

VASPBERRY evaluates geometric and response properties from VASP electronic
states. This report illustrates the discrete-wavefunction and Kubo approaches
with monolayer MoS₂, a Bi bilayer and a three-septuple-layer MnBi₂Te₄ film.
The examples cover reciprocal-space Berry
curvature, Chern and Z₂ indices, intrinsic charge Hall response, circular
optical transitions and real-space spinor densities. Berry-curvature maps use
Cartesian coordinates, band paths use distances along named high-symmetry
directions, and the Z₂ n-field is shown in reduced coordinates to display its
half-zone sums. The calculation guides provide the VASP files and executable commands.

## 1. Methods

### 1.1 Berry curvature from wavefunction overlaps

The Fukui–Hatsugai–Suzuki method evaluates the Berry phase around each cell
of a periodic k mesh using overlaps between neighboring wavefunctions.
For an isolated group of bands, the overlap is a matrix and its determinant
gives a gauge-invariant loop phase. With the connection
$\mathbf A=i\langle u|\nabla_\mathbf{k}u\rangle$ used here,

$$
\Phi_p=-\arg\!\left(U_1(\mathbf k)
U_2(\mathbf k+\Delta\mathbf k_1)
U_1(\mathbf k+\Delta\mathbf k_2)^{-1}
U_2(\mathbf k)^{-1}\right),\qquad
\overline\Omega_p=\frac{\Phi_p}{\Delta S_k}.
$$

The normalized link $U_j$ is the phase of the overlap determinant.
$\Delta S_k=|\mathbf b_1\times\mathbf b_2|/(N_1N_2)$ is the area of a mesh
cell, with reciprocal vectors containing $2\pi$. Thus $\Phi_p$ is a phase
and $\overline\Omega_p$ is an area-averaged curvature in Å². The Chern number
is $C=(2\pi)^{-1}\sum_p\Phi_p$. An isolated occupied group can include internal
degeneracies; its separation from excluded bands remains essential.
The discrete construction and its band-group extension are described by
[Fukui, Hatsugai and Suzuki](https://doi.org/10.1143/JPSJ.74.1674).

The map and integral convey different information. A material can have strong
local curvature at opposite valleys while its total Chern number is zero.
An integer lattice Chern result also needs a sufficiently resolved mesh and
a physically appropriate isolated band group.

### 1.2 Kubo curvature and Hall conductivity

For an isolated band and a Cartesian Hamiltonian derivative
$D_\alpha=\partial H/\partial k_\alpha$,

$$
\Omega_{n,xy}(\mathbf k)=-2\,\mathrm{Im}
\sum_{m\ne n}\frac{D_{x,nm}D_{y,mn}}{(E_n-E_m)^2}.
$$

This form resolves curvature at each sampled k point. Sharp curvature near
small gaps makes the k sampling and intermediate-band window particularly
important. Individual-band values are not used where the sampled bands
cannot be resolved separately. The native WAVECAR implementation uses
canonical momentum of the stored pseudo-wavefunctions; full material velocity
can require PAW, nonlocal and SOC terms. The exported-matrix interface permits
an explicitly supplied physical operator.

For an isolated occupied bundle $\mathcal V$, the trace curvature is

$$
\Omega^{\mathcal V}_{xy}=-2\,\mathrm{Im}
\sum_{n\in\mathcal V}\sum_{m\notin\mathcal V}
\frac{D_{x,nm}D_{y,mn}}{(E_n-E_m)^2}.
$$

Internal occupied-pair terms cancel in the sum of band curvatures. Excluding
them before division avoids singular individual-band terms at internal
degeneracies. The bundle must remain separated from the excluded bands.
This occupied-to-unoccupied formulation and its numerical advantage are
given by [Wang et al., Eq. (11) and Sec. III D](https://doi.org/10.1103/PhysRevB.74.195118).

For a two-dimensional system, the intrinsic charge sheet conductivity is

$$
\sigma_{xy}(\mu,T)=-\frac{e^2}{\hbar}
\sum_n\int_\mathrm{BZ}\frac{d^2k}{(2\pi)^2}\,
 f(E_{n\mathbf k}-\mu,T)\Omega_{n,xy}(\mathbf k).
$$

Within an insulating gap at zero temperature, this becomes
$\sigma_{xy}=-(e^2/h)C_\mathrm{occ}$ in the stated convention. General
background on Berry curvature and electronic transport is given in
[Xiao, Chang and Niu](https://doi.org/10.1103/RevModPhys.82.1959).
The [Bi Hall example](../examples/features/hall-valley/) evaluates the occupied-group flux directly;
it does not assign independent curvatures to unresolved Kramers partners.

For chemical potentials crossing band edges, the Kubo implementation uses
unordered interband pairs. With $N^{xy}_{nm}=-2\,\mathrm{Im}(D_{x,nm}D_{y,mn})$,
the occupation-weighted integrand is
$\sum_{n<m}(f_n-f_m)N^{xy}_{nm}/(E_n-E_m)^2$.
Equal-occupation pairs cancel before division. Changes from a gap reference
are evaluated with occupation differences before summation. This preserves
small doping responses without subtracting large filled-band baselines.
The temperature enters only the Fermi function; the electronic structure
remains fixed. Regional contributions integrate the same charge response over
specified parts of the BZ, with a K−K′ difference defined without a factor of
one-half. Such a partition is not a separately conserved valley-current operator.

### 1.3 PAW optical matrices and full-connection interpolation

For an insulating occupied bundle, the supported standard VASP longitudinal
optical calculation supplies matrix elements
$C_{cv,\alpha}=\langle u_c|\partial_{k_\alpha}u_v\rangle$ in Å. Here $v$ and
$c$ label occupied and empty states. Their contribution to the trace is

$$
\Omega^{\rm occ}_{xy}=-2\,\mathrm{Im}
\sum_{v,c} C_{cv,x}^{*} C_{cv,y}.
$$

The energy denominator is already contained in these derivatives. The
longitudinal PAW optical expression includes projector and augmentation terms,
as described by [Gajdoš et al.](https://doi.org/10.1103/PhysRevB.73.045112).
The `waveder-hall` route checks the actual VASP output, occupied filling and
global gap before integration. Its present scope is VASP 5.4.4, zero
temperature and a fixed insulating bundle; the occupied–empty separation
must exceed the producer's 2 meV degeneracy threshold. Mesh and empty-state
convergence remain separate requirements.

For narrow curvature peaks, VASPBERRY also interpolates the full Hamiltonian
and position matrices of a VASP-derived Wannier model. Its own Fourier,
diagonalization and occupied-trace routines evaluate $\Omega=J0+J1+J2$:
the basis-connection curl, the mixed connection/Hamiltonian-derivative term,
and the derivative-pair term. Only occupied-to-empty denominators occur;
internal degeneracies are allowed. The formulation follows
[Wang et al.](https://doi.org/10.1103/PhysRevB.74.195118) and
[Lopez et al., Eq. (51)](https://doi.org/10.1103/PhysRevB.85.014435).
Bands and local curvature are checked against direct VASP data. The separate
`postw90` calculation in [Wannier90](https://doi.org/10.1088/1361-648X/ab51ff)
provides an independent check of VASPBERRY's evaluation of the same finite model.

## 2. Materials and numerical settings

| Dataset | Wavefunctions and sampling | Quantities shown |
|---|---|---|
| Monolayer MoS₂, full zone | SOC; 12×12 mesh; 26 bands; 400 eV | Fukui and Kubo occupied-bundle curvature, bands 1–18 |
| Monolayer MoS₂, matching path | SOC; 49 K–Γ–K′ points; 26 bands; 400 eV | Band structure and Kubo curvature beside the maps |
| Monolayer MoS₂, local valleys | SOC; two 9×9 patches; 26 bands; 400 eV | Isolated band-18 curvature near K and K′ |
| Monolayer MoS₂, transport | SOC; 12×12 to 36×36 meshes; stored and retained band windows varied separately; 400 eV | Chemical-potential and temperature dependence of regional Hall response |
| Monolayer MoS₂, supplied path | SOC; 48 K–Γ–K′ points; 32 bands; 400 eV | Optical spectra and Γ-point state density |
| Buckled Bi bilayer | SOC; full 12×12 mesh; 18 bands; occupied bands 1–10 | Fukui–Hatsugai n-field and Z₂ invariant |
| MnBi₂Te₄, three septuple layers | SOC+U; 21 atoms; full 6×6 VASP mesh; 192 bands; occupied bands 1–123; 270 eV | Magnetic Chern invariant, bands and Hall-integration checks |

The [input guide](../examples/INPUTS.md) links the structures and VASP files.
The MoS₂ maps and their matching path use the same public SCF charge density,
structure, potentials, cutoff and 26-band window. VASP 5.4.4 generated the
12×12 mesh and 49-point path with `ICHARG=11` and `ISYM=-1`; only the k-point
list changes. The sampled occupied-to-empty direct gap is 1.674 eV. The
tutorial supplies preparation files and commands for both WAVECARs.
The optical and real-space examples retain the separate supplied 32-band path.

## 3. Results

### 3.1 MoS₂: Fukui Berry curvature over the full Brillouin zone

![Fukui Berry curvature map, band structure and symmetry-path cut](../examples/features/fukui-berry-curvature/reference/smooth/figure.png)

**Figure 1.** (a) Fukui Berry curvature of occupied bands 1–18 in the
Cartesian first Brillouin zone. The dashed line marks K–Γ–K′, with
K = (1/3, 2/3) and K′ = −K in reciprocal coordinates. Native 12×12 plaquette
values are displayed with periodic bilinear interpolation onto a 401×401
Cartesian grid. (b) Band structure along that path;
energies are relative to the valence-band maximum, with occupied bands in
blue and empty bands in gray. (c) A periodic bilinear cut of the plaquette
field along the marked path. This line visualizes panel (a), rather than
adding an independent pointwise curvature calculation. Both right-hand
panels share the same Cartesian path distance.

The curvature has opposite signs near the time-reversed K and K′ valleys.
The largest sampled magnitudes are approximately **12.31 Å²**, while the
full-zone Chern number is zero to numerical precision. The native text map
integrates to approximately −5 × 10⁻⁷ after decimal rounding. This illustrates the distinction between
valley-contrasting local geometry and vanishing total charge Chern number.
The connection between inversion breaking, spin–orbit coupling and valley
optical response in MoS₂ is discussed by
[Xiao et al.](https://doi.org/10.1103/PhysRevLett.108.196802).
The present map is a finite-mesh example rather than a convergence study.

[Calculation and plotting commands](../examples/features/fukui-berry-curvature/)

### 3.2 MoS₂: Kubo curvature of the occupied valence-band bundle

![MoS2 occupied-bundle Kubo curvature map, band structure and symmetry-path curvature](../examples/features/kubo-curvature/reference/bundle/figure.png)

**Figure 2.** (a) Canonical-momentum Kubo trace curvature of occupied bands
1–18 on the 12×12 mesh. Only transitions to the stored empty bands 19–26
enter the sum. The map uses the same periodic bilinear display interpolation
and dashed K–Γ–K′ path as Figure 1. (b) Band structure, with occupied bands
in blue and empty bands in gray. (c) Bundle curvature calculated directly
at all 49 path points. Internal valence-band degeneracies do not interrupt
the bundle curve.

The occupied-bundle curvature is approximately **−13.18 Å² at K** and
**+13.18 Å² at K′**. The occupied space stays separated from the empty states
by at least **1.674 eV** on both samplings, so every mesh and path point is
valid for the bundle calculation. The native option `-kubo_bundle 1` performs
the external-band sum directly, avoiding large cancelling internal terms.

Figures 1 and 2 now describe the same occupied band space. Their numerical
values need not coincide at this resolution: Fukui gives finite-plaquette
averages, while this Kubo result uses pointwise canonical momentum and a
finite empty-band window. The 12×12 mesh and 26 stored bands require
convergence for quantitative predictions. NumPy interpolation smooths only
the display; all integration and numerical checks use the original samples.

[Calculation and interpretation](../examples/features/kubo-curvature/)

### 3.2.1 MoS₂: an isolated single band near the valleys

![Isolated MoS2 valence-band curvature near K and Kprime](../examples/features/kubo-curvature/valleys/reference/figure.png)

**Figure 3.** Single-band Kubo curvature of the upper valence band, band 18,
around K and K′. Each local Cartesian patch extends ±0.12 Å⁻¹ and contains
9×9 VASP points. The upper maps use NumPy bilinear display interpolation;
the dashed lines mark the cuts below. The lower panels show the two upper
valence bands near both valleys and band-18 curvature at both valleys. Symbols mark
the calculated curvature samples.

Band 18 touches its Kramers partner at Γ, but is separated from every other
stored band by at least **0.130 eV** throughout these valley patches. The
17–18 splitting at K is **0.147 eV**, and the valley curvatures are
approximately **−6.415/+6.415 Å² at K/K′**. This is a suitable local single-band example of
the spin-split valley physics described by
[Xiao et al.](https://doi.org/10.1103/PhysRevLett.108.196802).
The patch data describe local curvature; a full-zone integral requires
the appropriate isolated band bundle and complete BZ sampling.

[VASP preparation, calculation and reference data](../examples/features/kubo-curvature/valleys/)

### 3.2.2 MoS₂: intrinsic charge and regional valley Hall response

![MoS2 regional Hall response, bands and mesh refinement](../examples/features/kubo-hall/reference/figures/mos2-hall.png)

**Figure 4.** (a) Cartesian first BZ with periodic K/K′ disks of radius
0.35 Å⁻¹ and the marked K–Γ–K′ line. (b) The matching 26-band dispersion,
with the chemical-potential interval shaded. (c) Regional Hall changes at
300 K on the 36×36 mesh, using 60 stored bands and pairs within bands 1–40.
(d) Mesh refinement at that fixed band window. Energies are relative to each
source's valence-band maximum; the band path uses the same density and cutoff.

The 61-point scan covers $E_v-0.20$ to $E_v+0.10$ eV at 0 and 300 K.
We define $\Delta\sigma(\mu)=\sigma(\mu)-\sigma(\mu_{\rm ref})$ with a midgap
reference, and the valley difference as $\Delta\sigma_K-\Delta\sigma_{K'}$.
The total charge response remains below **$2.24\times10^{-9}\,e^2/h$**,
consistent with time-reversal symmetry. Numerical tables also give absolute
conductivities and carrier counts.

Two numerical choices were varied separately at 300 K. The relative L2
difference below is the norm of the curve change divided by the finer curve's
norm, on a common $\mu-E_v$ grid.

| Test | Refinement | Relative L2 change |
|---|---|---:|
| k mesh; 60 stored bands, pair cutoff 40 | 12×12 → 24×24 | 12.89% |
| k mesh; same band window | 24×24 → 36×36 | 1.482% |
| Pair cutoff; one 12×12, 96-band source | 60 → 80 | 1.263% |
| Pair cutoff; same source | 80 → 90 | 0.202% |

The final mesh change exceeds the stated 1% criterion. The last cutoff
increment is below 1%; these separate tests do not establish joint
convergence. The [reference guide](../examples/features/kubo-hall/) includes
the commands and additional checks on empty-state accuracy and numerical
degeneracy grouping.

This rigid-band response uses canonical momentum, omitting PAW augmentation
and nonlocal/SOC velocity terms. Numerical refinement does not remove those
approximations. The [supplementary plot](../examples/features/kubo-hall/reference/figures/mos2-hall-checks.pdf)
retains the finite-mesh steps at 0 K.

### 3.3 MoS₂: valley-selective circular transitions

![Circular optical spectra in MoS2](../examples/features/circular-dichroism/reference/figure.png)

**Figure 5.** Left- and right-circular spectra at K and K′ and the selectivity
$\eta=(I_L-I_R)/(I_L+I_R)$ along the path. Normal incidence is used, with
0.05 eV Gaussian broadening and transitions from occupied bands 1–18 to
bands 19–20. White regions in the selectivity map exclude negligible intensity.

The first peak occurs near **1.67 eV**. Its selectivity is −1 at K and +1
at K′ in VASPBERRY's channel convention. The opposite helicities provide a
clear illustration of valley-selective transitions. Intensities are in
arbitrary units and retain the native sampling normalization; an absolute
absorption spectrum requires a full-zone calculation and the appropriate
optical matrix elements.

[Optical calculation guide](../examples/features/circular-dichroism/)

### 3.4 Bi: Fukui–Hatsugai Z₂ invariant and n-field map

For the Bi bilayer, VASPBERRY evaluates the
[Fukui–Hatsugai lattice n-field](https://doi.org/10.1143/JPSJ.76.053702)
from the occupied SOC bands 1–10 on a full Γ-centered 12×12 mesh.
With the time-reversal-compatible gauge used in this construction, the
Z₂ invariant is the parity of the integer-field sum over a half Brillouin zone:

$$
\nu=\left[\sum_{p\in B_{1/2}} n(p)\right]\bmod 2.
$$

![Bi Fukui-Hatsugai integer n-field and Z2 invariant](../examples/features/z2/reference/figure.png)

**Figure 6.** Fukui–Hatsugai integer field $n(\mathbf k)$ calculated from the
Bi VASP spinor wavefunctions. Red, white and blue denote +1, 0 and −1,
respectively. The native plaquettes are displayed without interpolation in
dimensionless reduced coordinates, $\mathbf k=q_1\mathbf b_1+q_2\mathbf b_2$,
with opposite edges periodically identified. The line at $q_2=0$ separates
the upper and lower half zones. The occupied-band calculation gives **Z₂ = 1**.

The upper and lower half-zone sums are **−3 and +3**, respectively. Both are odd,
so their parities agree at **$\nu=1$**, identifying the nontrivial
time-reversal-symmetric insulating phase for this calculation.
The sampled minimum direct and global gaps are **0.592 eV** and **0.510 eV**.

The local n-field depends on the gauge and logarithm branch; its half-zone
parity is the invariant. The reduced-coordinate map shows the original mesh
half-zones used for these sums. The separate occupied **C = 0** and zero charge-Hall
response are consistent with this nontrivial Z₂ result.

[Z₂ calculation and plotting guide](../examples/features/z2/) ·
[Numerical n-field data](../examples/features/z2/reference/Z2_FIELD.csv)

### 3.5 MoS₂: a real-space Γ state

![MoS2 Gamma state in Cartesian real space](../examples/features/wavefunction/reference/figure.png)

**Figure 7.** Projected pseudo-wavefunction density of band 18 at Γ and its
plane average along z. The in-plane plot uses the Cartesian geometry of the
oblique primitive cell. Both spinor components are retained.

The integrated pseudo-density is approximately **0.8258**, consistent with
the plane-wave coefficient norm. The PAW augmentation is not part of this
plot, so a unit density integral is not imposed. The example illustrates
conversion of complex amplitudes into a physical-coordinate density map.

[Wavefunction calculation guide](../examples/features/wavefunction/)

### 3.6 MnBi₂Te₄: Chern number and VASPBERRY Hall response

This example connects an occupied-band topological invariant to a Hall integral
computed by VASPBERRY. The three-septuple-layer MnBi₂Te₄ film has 21 atoms
with +z/−z/+z Mn magnetization. Its fixed, unrelaxed geometry retains bulk
vertical spacings from Yan et al., with in-plane $a=4.336$ Å, PBE+SOC and
Mn $U_{\rm eff}=5.34$ eV following Otrokov et al. It is distinct from the
relaxed films in that study.

![MnBi2Te4 bands and VASPBERRY Hall convergence](../examples/materials/mnbi2te4-qah/reference/figures/mnbi2te4-qah.png)

**Figure 8.** (a) VASPBERRY-interpolated bands along Γ–M–K–Γ, with direct
VASP checks; the inset resolves the K–Γ–M gap near Γ. (b) VASPBERRY
full-connection sheet Hall response at three chemical potentials inside the
gap, at zero temperature; open markers show the independent `postw90`
comparison. The dashed line is the topological expectation. (c) Deviation
$|\sigma_{xy}/(e^2/h)-1|$ under integration refinement. Labels give the base
mesh, local subdivision and Γ-region radius in Å⁻¹. Both solvers use the
same model and weights; their agreement is distinct from mesh convergence.

The VASP gap is **17.10 meV**. Fukui links of the **complete occupied
valence-band bundle, bands 1–123, give C = −1**, corresponding to
$\sigma_{xy}=+e^2/h$ in our convention. A direct 6×6 PAW optical integral
instead gives **163.109 $e^2/h$**: the Γ-point curvature reaches approximately
$-1.52\times10^4$ Å², and this coarse mesh overweights the narrow peak.
An integer Fukui result does not establish convergence of the pointwise integral.

For dense integration, VASPBERRY's `wannier-hall` evaluates all J0/J1/J2
terms from the supplied VASP-derived Hamiltonian and position matrices.
`wannier-bands` independently generates the plotted dispersion. The model
contains 138 orbitals and **87 occupied states**; the omitted 36 deep VASP
bands form a separately checked C = 0 bundle. Near Γ, direct VASP band edges
agree within 1.26 meV and local curvature within 1.44%; the latter comparison
retains the full-123 versus model-87 subspace difference. The maximum
band-edge difference across all twelve validation locations is 18.14 meV.

VASPBERRY's final 27,600-point result is
**$\sigma_{xy}=1.000384978\,e^2/h$** at all three gap chemical potentials.
It uses an 80×80 base mesh and 9×9 subdivisions in the 0.18 Å⁻¹ Γ region.
The last refinement changes the result by **$2.07\times10^{-4}\,e^2/h$**,
below the stated $10^{-3}\,e^2/h$ comparison criterion. The nonmonotonic
region-expansion control is retained; no integer rounding is applied.
Independent `postw90` totals and all three component decompositions agree
within their original output precision. A separate 160×160 energy scan
places all three chemical potentials inside the sampled model gap.

This validates VASPBERRY's full-connection evaluation of a fixed finite model.
The source 6×6 mesh, structure and basis still require convergence for
quantitative material predictions; localization stopped before its spread
tolerance. The [material tutorial](../examples/materials/mnbi2te4-qah/NATIVE_WANNIER.md)
provides the actual inputs, VASPBERRY commands and outputs. The independent
external reference is retained separately. This input route requires both
Wannier operators; it does not infer missing PAW/nonlocal/SOC velocity terms
from a WAVECAR alone.

## 4. Practical use

Begin with the [example guide](../examples/README.md), reproduce the chosen
quantity, and then select the sampling and band window appropriate to your
own material. Use Cartesian reciprocal-space maps to display lattice symmetry,
named k paths for band-resolved quantities, and energy or chemical potential
for spectra and transport. State units, represented bands, temperature and
broadening in the caption.

For quantitative conclusions, converge the VASP electronic structure,
k mesh and relevant band window. The examples show how the methods are used;
the [supplementary analytic references](REFERENCE_MATERIALS.md) check formulas
and conventions, while [validation details](VALIDATION_1.3.0.md) document the
numerical implementation. Detailed execution records remain alongside the
numerical files for reproducibility.

## References

1. T. Fukui, Y. Hatsugai and H. Suzuki, *Chern Numbers in Discretized Brillouin
   Zone: Efficient Method of Computing (Spin) Hall Conductances*,
   [J. Phys. Soc. Jpn. **74**, 1674–1677 (2005)](https://doi.org/10.1143/JPSJ.74.1674).
2. D. Xiao, M.-C. Chang and Q. Niu, *Berry Phase Effects on Electronic
   Properties*, [Rev. Mod. Phys. **82**, 1959–2007 (2010)](https://doi.org/10.1103/RevModPhys.82.1959).
3. D. Xiao, G.-B. Liu, W. Feng, X. Xu and W. Yao, *Coupled Spin and Valley
   Physics in Monolayers of MoS₂ and Other Group-VI Dichalcogenides*,
   [Phys. Rev. Lett. **108**, 196802 (2012)](https://doi.org/10.1103/PhysRevLett.108.196802).
4. T. Fukui and Y. Hatsugai, *Quantum Spin Hall Effect in Three Dimensional
   Materials: Lattice Computation of Z₂ Topological Invariants and Its
   Application to Bi and Sb*,
   [J. Phys. Soc. Jpn. **76**, 053702 (2007)](https://doi.org/10.1143/JPSJ.76.053702).
5. X. Wang, J. R. Yates, I. Souza and D. Vanderbilt, *Ab initio calculation
   of the anomalous Hall conductivity by Wannier interpolation*,
   [Phys. Rev. B **74**, 195118 (2006)](https://doi.org/10.1103/PhysRevB.74.195118).
6. M. Gajdoš et al., *Linear optical properties in the projector-augmented
   wave methodology*,
   [Phys. Rev. B **73**, 045112 (2006)](https://doi.org/10.1103/PhysRevB.73.045112).
7. M. M. Otrokov et al., *Unique Thickness-Dependent Properties of the van der
   Waals Interlayer Antiferromagnet MnBi₂Te₄ Films*,
   [Phys. Rev. Lett. **122**, 107202 (2019)](https://doi.org/10.1103/PhysRevLett.122.107202).
8. J.-Q. Yan et al., *Crystal growth and magnetic structure of MnBi₂Te₄*,
   [Phys. Rev. Materials **3**, 064202 (2019)](https://doi.org/10.1103/PhysRevMaterials.3.064202).
9. G. Pizzi et al., *Wannier90 as a community code: new features and applications*,
   [J. Phys.: Condens. Matter **32**, 165902 (2020)](https://doi.org/10.1088/1361-648X/ab51ff).

10. M. Lopez, D. Vanderbilt, T. Thonhauser and I. Souza, *Wannier-based
    calculation of the orbital magnetization in crystals*,
    [Phys. Rev. B **85**, 014435 (2012)](https://doi.org/10.1103/PhysRevB.85.014435).
