# Compute the full Wannier Hall response with VASPBERRY

This route uses **VASPBERRY's own NumPy solver** to Fourier-transform the
Hamiltonian and position matrices, diagonalize the Hamiltonian, evaluate the
occupied-bundle Berry curvature, and integrate the sheet Hall response.
It retains all three terms J0 + J1 + J2. It does not invoke postw90.

The inputs are the actual VASP-derived operators supplied with this material
example. Wannier90 was used to construct the finite representation; the Hall
evaluation below is performed by VASPBERRY. The independently computed
[postw90 reference](WANNIER_REFERENCE.md) provides a separate check of the
same operators and integration grid.

## Why these additional inputs are needed

The usual WAVECAR momentum matrix omits PAW/nonlocal velocity corrections.
The standard WAVEDER route includes the audited optical matrix elements, but
the 6×6 MnBi₂Te₄ mesh misses the narrow Γ-region integral. Its large response
is a sampling diagnostic.

Dense interpolation therefore uses both the Hamiltonian and the full position
connection from the VASP-derived Wannier representation. A WAVECAR alone does
not contain these position matrices. Keeping only Hamiltonian derivatives
would give J2 and would omit J0 and J1.

This implementation uses a zero-temperature occupied projector. It avoids
energy denominators within the occupied bundle, allowing internal band
degeneracies, and rejects a gap closure between the selected occupied and
excluded states. The example does not claim a finite-temperature or metallic
transport calculation.

## 1. Restore the actual operator files

Run all commands from the repository root. Install the transport dependencies
(including NumPy and Matplotlib), and limit each numerical-library worker to
one thread:

```sh
python -m pip install -r requirements-transport.txt
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
```

The output directory below must be new:

```sh
mkdir -p work
python examples/materials/mnbi2te4-qah/inputs/wannier/operators/restore.py \
  examples/materials/mnbi2te4-qah/inputs/wannier/operators \
  work/mbt3-full-operators
```

The restored files are `wannier90_HH_R.dat` and `wannier90_AA_R.dat`, together
about 201 MiB. The archive is about 65 MiB. Both files are needed; the three
Cartesian position components are in the AA file.

## 2. Import the operators into VASPBERRY

```sh
python tools/vaspberry_kubo.py wannier-import \
  --hh work/mbt3-full-operators/wannier90_HH_R.dat \
  --aa work/mbt3-full-operators/wannier90_AA_R.dat \
  --poscar examples/materials/mnbi2te4-qah/inputs/POSCAR \
  --spinor-components 2 --spin-multiplicity 1 \
  --energy-reference "Unshifted VASP eigenvalues in eV; DFT midgap 3.2619863015705284 eV" \
  --output-dir work/mbt3-operator-cache
```

The importer checks the real-space matrix records and their consistency and
writes the VASPBERRY operator cache. These files use the documented effective
HH/AA format, including additive records. They already contain the original
Wigner–Seitz weights and translated vectors. Do not divide by those weights
again, apply another distance correction, or rename an arbitrary `_hr.dat`
file into this format.

This model has 138 spinor functions and **87 occupied model bands**. The
direct VASP Fukui calculation uses 123 occupied bands. Their difference is
the separately checked deep bundle, original bands 1–36, with C = 0. The
model's local curvature can still differ from the full 123-band trace.

## 3. Interpolate the VASP-derived band structure

The following command evaluates all 138 bands of the VASP-derived Wannier
Hamiltonian along Γ–M–K–Γ. VASP supplies the electronic structure;
VASPBERRY evaluates its interpolation for the response calculation. There are
1,201 points per segment; shared endpoints are counted once, giving 3,601
points in the complete path:

```sh
python tools/vaspberry_kubo.py wannier-bands \
  --operators work/mbt3-operator-cache \
  --vertices 0 0 0 0.5 0 0 0.3333333333333333 0.3333333333333333 0 0 0 0 \
  --labels Gamma M K Gamma --points-per-segment 1201 \
  --formats npz --output-dir work/mbt3-native-bands
```

`bands.npz` contains the actual k points, path distances, energies and lattice;
`bands.json` records the producer, units and source. Energies retain the input
zero. The example figure subtracts the recorded midgap only for presentation.

## 4. Integrate with the VASPBERRY solver

```sh
python tools/vaspberry_kubo.py wannier-hall \
  --operators work/mbt3-operator-cache --occupied 87 \
  --mesh 80 80 --refine 9 --refine-radius 0.18 --refine-center 0 0 \
  --mu-min 3.254290679571569 --mu-max 3.269681923569488 --mu-num 3 \
  --mu-reference 3.2619863015705284 \
  --workers 4 --batch-size 16 --time-limit 3600 \
  --formats csv dat npz --output-dir work/mbt3-native-hall
```

The base grid covers the full two-dimensional torus. Cells around periodic Γ
within the stated 0.18 Å⁻¹ selection radius are subdivided into 9×9 cells,
preserving their area. The final grid has 27,600 points. Each chemical
potential uses the same points and weights. The three chemical potentials
span the central 90% of the sampled DFT gap, in the unchanged VASP energy zero.

Use `--workers 1` on a smaller computer, and reduce `--batch-size` when memory
is limited. More workers and batches use more memory; choose them together.
The time limit is an execution bound, not a convergence criterion. Keep the
recorded run status and outputs when checking a result.

## 5. Inspect the outputs and convergence

The output directory contains:

- `hall/conductivity.csv`, `.dat`, and `.npz`: matching Hall-response tables.
- `hall/conductivity.json`: model, occupied gap, mesh, full-connection formula
  and integrated J0/J1/J2 components.
- `curvature.npz`: the computed curvature data for the integration points.
- `run.json`: completion status, timing and source-integrity records.

The full-connection output identifies VASPBERRY as the solver. The independent
postw90 tables remain separately identified by their producer. With this
example's conventions, C = −1 corresponds to sheet σxy = +e²/h, with no extra
spin or surface factor. A flat response inside a gap follows from unchanged
occupations; its value and integration convergence establish whether it is
quantized numerically.

Repeat the command into new output directories for the following controls:

| Base mesh | Local subdivision | Refinement radius (Å⁻¹) |
|---|---|---:|
| 40×40 | 7×7 | 0.12 |
| 60×60 | 7×7 | 0.12 |
| 60×60 | 11×11 | 0.12 |
| 60×60 | 11×11 | 0.18 |
| 80×80 | 9×9 | 0.18 |

These separately check base sampling, local resolution and the refined domain.
Retain every control, including nonmonotonic changes. All five native runs
have been completed; their [outputs and convergence table](reference/vaspberry-wannier/README.md)
are supplied. The final VASPBERRY value is **1.0003849781827807 e²/h** at each
of the three chemical potentials. The final same-domain refinement changes
it by 0.0002069516 e²/h, below the predefined 10⁻³ e²/h integration criterion.
The final run took 9.6 minutes with four NumPy workers, with about 2.2 GiB
peak resident memory on the reference host.

The independent postw90 result is **1.0003849765685957 e²/h**. Every native
control agrees with its independent counterpart within the external output
precision, accounting for the documented historical constant convention.
This is a comparison of independently evaluated responses; no postw90
response table or executable is an input to the native solver.

## Applying this to another system

Prepare matching Hamiltonian and position matrices from that system's
validated first-principles Wannier representation, import them with its own
lattice and energy reference, choose its occupied bundle, and repeat the
gap and integration checks. Bands, local curvature and gauge conventions
need independent validation. A band fit alone does not validate the position
connection.

This supplied model is fixed and unrelaxed, uses a coarse 3×3 source density
and 6×6 Wannier training mesh, and stopped localization after 300 iterations
before the requested spread tolerance. Agreement between two solvers tests
the calculation of this model; it does not establish full material convergence.

For the mathematical background see
[Lopez et al., PRB85,014435(2012), Eq.(51)](https://doi.org/10.1103/PhysRevB.85.014435)
and [Wang et al., PRB74,195118(2006)](https://doi.org/10.1103/PhysRevB.74.195118).
The [general Wannier-transport guide](../../../docs/WANNIER_TRANSPORT.md)
describes the input contract, formula and common output schema.
