# MnBi₂Te₄: a magnetic Chern-insulator example

This example starts from an actual VASP calculation of a three-septuple-layer
MnBi₂Te₄ slab. The main workflow computes the occupied-bundle Fukui invariant from VASP
wavefunctions, imports the actual VASP-derived Hamiltonian and position
operators, and evaluates their full Berry-curvature Hall response with
VASPBERRY. A separate postw90 calculation checks the interpolated result.

## Material and reference calculation

The 21-atom slab has alternating out-of-plane Mn moments (+z, −z, +z),
PBE+SOC, Mn U_eff = 5.34 eV, and a 270 eV plane-wave cutoff. There are
123 occupied spinor bands; their multiplicity is one. The supplied
[structure](inputs/POSCAR) is **fixed and unrelaxed**, constructed from the
published bulk fractional coordinates and vertical spacings of Yan et al.
The in-plane lattice constant a = 4.336 Å follows Otrokov et al.; the
surface-to-surface vacuum is 15 Å. It is not the original authors' relaxed
three-layer structure.

The structural coordinates are from [Yan et al., Physical Review Materials
3, 064202 (2019)](https://doi.org/10.1103/PhysRevMaterials.3.064202).
The magnetic-film setting follows [Otrokov et al., Physical Review Letters
122, 107202 (2019)](https://doi.org/10.1103/PhysRevLett.122.107202).
The gap below is this example's calculated value, rather than the relaxed
film gap reported in that paper.

| Quantity | Direct VASP reference |
|---|---:|
| Full NSCF mesh | 6×6×1 |
| Stored spinor bands | 192 |
| Sampled global gap | 17.1014 meV |
| Fukui C, occupied bands 1–123 | −1 |
| Fukui C, deep bands 1–36 | approximately 0 |
| Standard-WAVEDER Hall integral, σ_xy/(e²/h) | +163.1092, unconverged |

The Γ-point curvature is very narrow. The 6×6 optical integral greatly
overweights this peak and is retained as a convergence diagnostic. A
chemical-potential scan inside a gap is flat at zero temperature because
the occupied subspace is unchanged; flatness alone does not demonstrate
quantization. With the convention used here, a converged C = −1 insulator
has σ_xy = +e²/h.

## 1. Prepare the VASP inputs

Use the PAW datasets listed in [PSEUDOPOTENTIAL.md](PSEUDOPOTENTIAL.md).
VASP and licensed POTCAR files must be supplied locally. The bundled
[self-consistent density](inputs/scf/README.md) fixes the electronic source
state so the NSCF and post-processing steps can be repeated directly.
**This standard-WAVEDER workflow currently requires VASP 5.4.4**: the
audited `waveder-hall` adapter rejects other producer versions. Confirm the
version before running the optical calculation.

Run from the repository root:

```sh
python examples/materials/mnbi2te4-qah/prepare_vasp.py \
  --potcar /path/to/licensed/Mn-Bi-Te/POTCAR \
  --mesh 6 --nbands 192 --output-dir work/mbt3-optics
```

The output directory contains ordinary `INCAR`, `POSCAR`, `KPOINTS`,
`CHGCAR`, and `POTCAR` files. Run your noncollinear VASP executable there,
for example:

```sh
cd work/mbt3-optics
mpirun -np 8 /path/to/vasp_ncl > vasp.stdout.log
```

The input uses `ICHARG=11`, `ISYM=-1`, `LOPTICS=.TRUE.`,
`LPEAD=.FALSE.`, and `LNABLA=.FALSE.`. Retain the completed `WAVECAR`,
`WAVEDER`, `OUTCAR`, and input files. These standard optical files were
produced with unmodified VASP 5.4.4 for the reference. The supplied density
uses a coarse 3×3 SCF mesh; converge the density, structure, basis, magnetic
state, and NSCF sampling separately before using this setup for a material
prediction.

## 2. Calculate the occupied Fukui invariant and optical diagnostic

Return to the repository root and run:

```sh
python examples/materials/mnbi2te4-qah/run.py \
  --run-dir work/mbt3-optics --mesh 6 --output-dir work/mbt3-results
```

This invokes the general [wavefunction Fukui tool](../../../tools/wavecar_fukui.py)
for the occupied bundle and the common
[`waveder-hall` command](../../../tools/vaspberry_kubo.py) for the
zero-temperature occupied/empty optical sum. The output contains Fukui
maps and link-quality diagnostics, conductivity CSV/DAT/NPZ tables, and a
summary with the actual sampled band edges. No rounding or factor of two
is applied to the optical integral.

If a full mesh is divided into fixed-charge VASP runs, `--run-dir` accepts
all of those directories. Supply the complete corresponding
`--wavecar` for the Fukui links. The example checks the mesh coverage,
energies, plane-wave bases, and all complex coefficients against those
source runs. Independently recomputed or interpolated wavefunctions do not
satisfy this assembly contract.

## 3. Compute the full Hall response with VASPBERRY

Follow [the native Wannier workflow](NATIVE_WANNIER.md): restore the actual
Hamiltonian and position matrices, import them with the general
`wannier-import` command, then run `wannier-hall`. VASPBERRY performs the
Fourier transforms, diagonalization, occupied-bundle J0 + J1 + J2 calculation,
and numerical integration in its own NumPy backend. Postw90 is not called
by this solver.

The [supplied operators](inputs/wannier/operators/README.md) come from this
actual VASP-derived Wannier model. Both H and all three position-connection
components are needed; a WAVECAR alone does not contain the latter. No fitted
analytic Hamiltonian or smoothed Fukui map replaces the response.

Three distinctions matter for this example:

- The WAVECAR momentum route is a canonical-momentum approximation and omits
  PAW/nonlocal velocity terms.
- The standard-WAVEDER 6×6 result contains the audited optical matrix elements,
  but its grid badly undersamples the narrow Γ peak.
- A dense result calculated by postw90 establishes an independent reference.
  The VASPBERRY result must come from the native solver's own output.

All five native integration runs below were completed and independently
checked against postw90 on identical points and weights, including all three
curvature components and their J0/J1/J2 contributions.

The model has 138 spinor Wannier functions and 87 occupied model bands.
The omitted original bands 1–36 form a separately checked C = 0 bundle.
Its local curvature need not vanish. Near Γ, the full Wannier and
192-band PAW optical curvature differ by at most 1.44% at the seven
sampled points; corresponding band-edge differences are at most 1.26 meV.
Across all 12 comparison points, the largest band-edge difference is
18.14 meV. These differences include the different finite subspaces and
the omitted deep-band curvature, as well as interpolation error.

This is a **fixed, validated numerical model**: localization stopped after
300 optimization iterations before the requested spread tolerance was
reached. The quadrature study tests integration of that model. It does
not establish convergence with respect to the source 6×6 DFT mesh,
SCF density, basis, structural relaxation, or Wannier optimization.

### VASPBERRY results and convergence

The completed VASPBERRY full-connection integral is **σ_xy = 1.000384978 e²/h**
at each of three calculated chemical potentials spanning the central 90% of
the sampled gap. A separate 160×160 energy scan of the same model finds a
17.1014 meV global sampled gap, with all three points inside it.

| Base grid | Local subdivision | Refined radius (Å⁻¹) | Points | σ_xy/(e²/h) |
|---|---|---:|---:|---:|
| 40×40 | 7×7 | 0.12 | 3,088 | 1.000193491 |
| 60×60 | 7×7 | 0.12 | 6,528 | 1.000054804 |
| 60×60 | 11×11 | 0.12 | 10,920 | 1.000049773 |
| 60×60 | 11×11 | 0.18 | 21,720 | 1.000591930 |
| 80×80 | 9×9 | 0.18 | 27,600 | 1.000384978 |

Each refined region replaces whole coarse cells by the stated two-dimensional
subdivision. The full periodic integration domain and total weight remain
unchanged; all chemical potentials use the same points and weights. The
last refinement changes σ_xy by 0.000207 e²/h, below the predefined
10⁻³ e²/h criterion. The wider-domain check is deliberately retained even
though it moves the value farther from an integer. These checks establish
numerical quantization to the reported quadrature accuracy for this fixed
model; they do not establish full material convergence. No result was rounded
to an integer or selected for being closest to one.

The [native reference data](reference/vaspberry-wannier/README.md) contain
matching common CSV/DAT/NPZ Hall tables, all sampled curvatures, run records
and all 138 bands. The [convergence table](reference/vaspberry-wannier/convergence.csv)
retains every native control beside its independent postw90 comparison.
The final native run took 9.6 minutes with four NumPy workers and about
2.2 GiB peak resident memory on the reference computer; timings vary by host.
The independently calculated postw90 value is 1.000384977 e²/h, agreeing
within its recorded output precision.

## 4. Plot the reference

```sh
python examples/materials/mnbi2te4-qah/plot.py \
  --output-dir work/mbt3-reference-figures
```

The figure combines VASP-derived Wannier bands with direct VASP sample markers and a
gap enlargement, the three native Hall samples inside the gap, and unrounded
integration convergence. Hollow markers show the independent postw90
calculation. VASPBERRY integrates the sheet response directly from reciprocal
area; the external S/cm reference uses the full simulation-cell height for
comparison, without a guessed slab thickness.

![MnBi₂Te₄ bands, sheet Hall response, and quadrature convergence](reference/figures/mnbi2te4-qah.png)

Blue band curves and filled Hall/convergence markers are calculated by
VASPBERRY from the actual VASP-derived H and position operators. Open band
circles are direct VASP samples; the gap inset follows K–Γ–M. Hollow Hall
and convergence markers are the independent postw90 values on the same
grids. All five prescribed quadrature controls are shown. The independent
occupied-bundle Fukui invariant is C = −1; the Hall values are not rounded.

## 5. Run the independent postw90 crosscheck

The independent reference uses **Wannier90/postw90** to check the
VASPBERRY full-connection calculation on the same finite operators. It retains the Hamiltonian and all
three position-connection matrices from the actual VASP-derived model;
postw90 evaluates the full J0 + J1 + J2 expression. This is not a fitted
analytic Hamiltonian or an integral of a smoothed Fukui map.

Follow [the external-reference workflow](WANNIER_REFERENCE.md) to restore
the supplied numerical operators, build the documented reader, and run
`wannier_reference.py prepare`, `run`, and `collect`. The output contains
ordinary conductivity CSV/DAT/NPZ tables with their own external-producer
metadata. The [operator archive](inputs/wannier/operators/README.md) is
about 65 MiB compressed and restores about 201 MiB; the original large
VASP/Wannier intermediate files are not needed for this route.

The physical reference was checked with unmodified official Wannier90
3.1.0. Its portable effective-model reader needs the documented
[one-line dimension initialization fix](inputs/wannier/toolchain/README.md).
The patch does not change the Berry formulas. All band energies and
curvature components at 12 test points, a 50-point Hall test, and the
complete 3,088-point integral agree with the original full-connection
producer at the printed precision.

The original VASP-to-Wannier export used a separately checked compatibility
correction for spin-labelled projections in the older VASP 5.4.4/Wannier90
2.1 interface. The standard `WAVEDER` calculation above used unmodified
VASP. No licensed VASP source, patch, executable, or POTCAR is distributed.
The supplied operators make the dense external benchmark reproducible
without that private interface. Generating a new model from another VASP
installation requires checking its own spinor-projection interface and
exported band/projection counts.


## Reference data

- [VASP-derived Wannier bands, VASPBERRY Hall response and all five convergence runs](reference/vaspberry-wannier/README.md).
- [Occupied-bundle Fukui map](reference/fukui/fukui_occupied.csv) and
  [link diagnostics](reference/fukui/diagnostics.json).
- [Actual VASP band samples](reference/direct-dft/sample-bands.csv),
  including nine independent points away from the training mesh.
- [PAW optical curvature samples](reference/direct-dft/optical-curvature.csv)
  and [sample metadata](reference/direct-dft/metadata.json).
- [Coarse optical Hall tables and interpretation](reference/waveder-coarse/README.md).

These tables distinguish direct calculations from any subsequent plotting
interpolation. The reference band energies retain the VASP energy zero;
the plotted zero may be shifted to the recorded midgap energy.

The external reference additionally provides [band data](reference/wannier/bands.csv),
[local VASP comparisons](reference/wannier/local-validation.csv), and
[model settings and limitations](reference/wannier/metadata.json).
