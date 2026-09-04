# Two-dimensional Z2 invariant by the Fukui-Hatsugai n-field method

VASPBERRY's `-z2 1` option evaluates the two-dimensional class-AII invariant
directly from a VASP spinor `WAVECAR`. It implements the lattice n-field
construction of Fukui and Hatsugai on a full, uniform, Gamma-centered mesh.

## Method implemented by VASPBERRY

Let `Psi(k)` contain the occupied cell-periodic spinors at a mesh point. For
each positive mesh direction, VASPBERRY forms the non-Abelian overlap matrix
and its normalized determinant link,

```math
U_mu(k) = det[Psi(k)^dagger Psi(k+mu)] /
          |det[Psi(k)^dagger Psi(k+mu)]|.
```

The four link arguments around a plaquette give a principal plaquette phase
`phi_p` and an integer field

```math
phi_p = wrap(sum_(links in boundary p) arg U_link),
```

```math
n_V(p) = [sum_(links in boundary p) arg U_link - phi_p] / (2 pi).
```

VASPBERRY constructs the time-reversal gauge required on the half Brillouin
zone. In its plane-wave convention,

```math
psi_k(r) = sum_G C_(G,k) exp[i(k+G).r],
```

the reciprocal-vector mapping used to represent the time-reversed state is

```math
G_target = -G_source + round(-k_source-k_target).
```

This last expression is a code-level consequence of the VASPBERRY plane-wave
basis convention; it is not quoted as an equation from the method papers.
The spinor action is
`Theta(C_up,C_down)=(-conjg(C_down),conjg(C_up))`, which differs from the
equivalent `i sigma_y K` convention only by an overall phase and satisfies
`Theta^2=-1`.

For either complementary half-zone,

```math
nu = sum_(p in B_half) n_V(p) mod 2.
```

Fukui and Hatsugai write `D_L=-sum n_12` (Ref. 1). With the plaquette
orientation used by VASPBERRY, the stored field is `n_V=-n_12`; therefore the
program's `sum n_V mod 2` is the same invariant. The top and bottom values are
not two independent invariants. They are complementary evaluations of the
same two-dimensional invariant and must agree modulo two.

## Valid input contract

The physical calculation must satisfy all of the following conditions:

- a two-dimensional, nonmagnetic, time-reversal-symmetric insulator;
- a full, unshifted, Gamma-centered `Nx x Ny x 1` mesh with even `Nx` and
  `Ny`, both at least 4, `kz=0`, and `ISYM=-1`;
- VASP spinors produced with `ISPIN=1` and `LSORBIT=.TRUE.`;
- an occupied subspace `1:Nocc` of even rank with a positive global insulating
  gap and separation from band `Nocc+1` at every represented k point; and
- `NBANDS>Nocc`, with at least the first unoccupied band retained in the
  `WAVECAR`.

VASPBERRY checks the mesh reconstruction, reciprocal-G bijections, link
conditioning, n-field integrality, time-reversal-odd wrapped flux, zero total
Chern residual, and equality of the two half-zone parities. A PASS establishes
the numerical self-consistency of the time-reversal reconstruction used by the
Z2 routine. It does not independently establish that the raw `WAVECAR` is a
nonmagnetic, time-reversal-symmetric insulating input, certify its
occupied-unoccupied gap, or perform a mesh-convergence study. Those physical
input and convergence checks remain the user's responsibility.

Setting `-ne Nocc` only overrides occupation inference when the isolated
occupied subspace has already been established. It cannot turn a metal or
semimetal into a valid insulating Z2 calculation.

## Preparing the WAVECAR with VASP

VASPBERRY itself reads only `WAVECAR`. The VASP inputs below are needed to
create it; `OUTCAR` and `EIGENVAL` should be retained to verify provenance and
the band gap.

### 1. Converge the charge density

Run an SOC calculation that writes a converged `CHGCAR`. Use the same
structure, PAW datasets, exchange-correlation functional, cutoff, and SOC
settings that will be used in the next step. The Bi example provides a
template under `examples/Bi_Z2/inputs/01_scf/`.

For example, run the following from the SCF directory after supplying the
local `POTCAR`:

```bash
VASP_BIN=/path/to/vasp_ncl
mpiexec -n 4 "$VASP_BIN" > vasp.out
```

The SOC-capable binary path, MPI launcher, and rank count are site-specific.
The example SCF template sets `LORBIT=11` so that site-projected x/y/z moments
can be inspected in `OUTCAR`; these PAW-sphere projections are qualitative
evidence rather than an independent proof of time-reversal symmetry.

### 2. Write a full-mesh spinor WAVECAR

Read that `CHGCAR` in a fixed-charge calculation with `ICHARG=11`. The
important settings are:

```text
ICHARG  = 11
ISPIN   = 1
LSORBIT = .TRUE.
MAGMOM  = 6*0.0   # two-atom nonmagnetic Bi example
ISYM    = -1
LREAL   = .FALSE.
GGA_COMPAT = .FALSE.
LASPH   = .TRUE.
LMAXMIX = 2       # Bi s/p PAW datasets
LWAVE   = .TRUE.
LCHARG  = .FALSE.
```

This CHGCAR-only `ICHARG=11` recipe is for the PBE/GGA workflow represented by
the example templates. Do not carry it over unchanged to hybrid or
kinetic-energy-density meta-GGA calculations, which require their
version-appropriate restart workflow and additional information. For hybrid
functionals, follow VASP's
[`hybrid-functional workflow`](https://vasp.at/wiki/Band-structure_calculation_using_hybrid_functionals)
and retain the required `WAVECAR` information instead.

For the 12 x 12 example, use:

```text
Z2 full Gamma-centered mesh
0
Gamma
12 12 1
0 0 0
```

Then run the fixed-charge stage from its own directory:

```bash
VASP_BIN=/path/to/vasp_ncl
mpiexec -n 4 "$VASP_BIN" > vasp.out
```

`ISYM=0` is not sufficient: VASP can still relate `k` and `-k` and reduce the
sampling. `ISYM=-1` disables that reduction. After the run, confirm that
`OUTCAR` reports `NKPTS=Nx*Ny` and that the selected occupied bands remain
separated from the next band at every represented k point, with a positive
sampled global gap. Check the gap again with denser sampling or interpolation.
Also verify that the total and, where evaluated, site-projected magnetization
are consistent with the assumed nonmagnetic state. Zero initial, total, or
projected moments do not by themselves prove physical time-reversal symmetry.

## Command line

From the repository root, compile VASPBERRY as described in
[`BUILD.md`](BUILD.md), then run the Bi 12 x 12 case as follows:

```bash
./build/vaspberry-gfortran \
  -f examples/Bi_Z2/WAVECAR \
  -o NFIELD \
  -z2 1 \
  -kx 12 -ky 12 \
  -s 2 \
  -ii 1 -if 10
```

The same options are accepted by the serial and MPI executables. With MPI,
for example:

```bash
mpiexec -n 4 ./build/vaspberry-mpi \
  -f examples/Bi_Z2/WAVECAR -o NFIELD -z2 1 \
  -kx 12 -ky 12 -s 2 -ii 1 -if 10
```

Run `./build/vaspberry-gfortran -h` to print the option contract and an
abbreviated example. If an executable is built directly under another name,
substitute that path in the same commands.

## Outputs and decision rule

On PASS, the calculation writes:

- `NFIELD.dat`: the traditional repeated-zone n-field view; and
- `Z2_FIELD.csv`: the `Nx x Ny` fundamental mesh, result metadata, numerical
  checks, flux, curvature, link quality, and integer n-field.

Rejected or interrupted Z2 runs retain `Z2_FIELD.invalid.csv` when the
guarded output preflight has started. Such a file is non-reportable, and an
`NFIELD.dat` is not a valid result unless the final `Z2_FIELD.csv` exists with
`result_status=PASS`.

A result is reportable by this method only when the CSV says
`result_status=PASS`, `reportable_invariant=1`, and both half-zone parities
agree. The common value is the Z2 invariant:

| top parity | bottom parity | interpretation |
|---:|---:|---|
| 0 | 0 | `Z2=0`, topologically trivial |
| 1 | 1 | `Z2=1`, topologically nontrivial |
| 0 | 1 | invalid; do not select either half |
| 1 | 0 | invalid; do not select either half |

Repeat the calculation on denser even Gamma-centered meshes. A PASS on one
mesh is a method result for that mesh, not proof of mesh convergence.

The pointwise n-field depends on the gauge and principal-logarithm branch.
Individual colored plaquettes are therefore not Berry-curvature hot spots,
and point-by-point mirror symmetry is not a validity condition. The invariant
comes from the half-zone integer sum modulo two together with the stated
guards.

## Plotting

The Bi helper reads a current schema-v2 PASS CSV directly:

```bash
python3 -m pip install -r requirements-transport.txt
```

```bash
python3 examples/Bi_Z2/scripts/plot_nfield.py \
  examples/Bi_Z2/reference-v1.2.0-12x12/Z2_FIELD.csv \
  --output examples/Bi_Z2/Z2_nfield_12x12.pdf
```

The figure shows the current 12 x 12 integer n-field and annotates the two
half-zone sums and parities. The helper also writes a PNG copy when the
requested output has a `.pdf` suffix.

The checked schema-v2 Bi reference contains 144 fundamental plaquettes. Its
top and bottom sums are -3 and +3, respectively, so both parities are 1 and
the result is `Z2=1` on this 12 x 12 mesh. This fixed-mesh PASS is not a
substitute for the physical input and mesh-convergence checks above.

## Implementation and numerical notes

- Serial and MPI builds call the same Z2 state, overlap, link, and field
  routines. The MPI path distributes plaquettes and reduces the resulting
  arrays before rank 0 publishes the atomic outputs.
- The time-reversal construction uses an explicit one-to-one reciprocal-G
  mapping for both spinor components and the same `Theta` operation for
  Kramers partners. Coefficient norms are accumulated after promoting their
  real and imaginary components.
- Link quality is measured by the minimum singular value of each overlap
  matrix. The determinant phase is accumulated from an LU factorization,
  avoiding rank-dependent determinant-magnitude underflow.
- Guarded output publication uses an invalid sentinel, same-directory staging,
  and an atomic final rename. A PASS `Z2_FIELD.csv` is the final commit marker.
- The direct reader uses pseudo-wavefunction coefficients from `WAVECAR` and
  does not add PAW augmentation terms. The backend and this limitation are
  recorded in the CSV; Ref. 5 discusses augmentation in PAW overlap matrix
  elements.
- Compiler and direct-access record-length requirements are documented in
  [`BUILD.md`](BUILD.md).

## References

1. T. Fukui and Y. Hatsugai, "Quantum Spin Hall Effect in Three Dimensional
   Materials: Lattice Computation of Z2 Topological Invariants and Its
   Application to Bi and Sb," *J. Phys. Soc. Jpn.* **76**, 053702 (2007),
   [doi:10.1143/JPSJ.76.053702](https://doi.org/10.1143/JPSJ.76.053702).
2. T. Fukui, Y. Hatsugai, and H. Suzuki, "Chern Numbers in Discretized
   Brillouin Zone: Efficient Method of Computing (Spin) Hall Conductances,"
   *J. Phys. Soc. Jpn.* **74**, 1674-1677 (2005),
   [doi:10.1143/JPSJ.74.1674](https://doi.org/10.1143/JPSJ.74.1674).
3. L. Fu and C. L. Kane, "Time reversal polarization and a Z2 adiabatic spin
   pump," *Phys. Rev. B* **74**, 195312 (2006),
   [doi:10.1103/PhysRevB.74.195312](https://doi.org/10.1103/PhysRevB.74.195312).
4. C. L. Kane and E. J. Mele, "Z2 Topological Order and the Quantum Spin Hall
   Effect," *Phys. Rev. Lett.* **95**, 146802 (2005),
   [doi:10.1103/PhysRevLett.95.146802](https://doi.org/10.1103/PhysRevLett.95.146802).
5. A. Ferretti, A. Calzolari, B. Bonferroni, and R. Di Felice, "Maximally
   localized Wannier functions constructed from projector-augmented waves or
   ultrasoft pseudopotentials," *J. Phys.: Condens. Matter* **19**, 036215
   (2007),
   [doi:10.1088/0953-8984/19/3/036215](https://doi.org/10.1088/0953-8984/19/3/036215).

VASP input semantics used above follow the official documentation for
[`ICHARG`](https://vasp.at/wiki/ICHARG),
[`ISYM`](https://vasp.at/wiki/ISYM),
[`LSORBIT`](https://vasp.at/wiki/LSORBIT),
[`MAGMOM`](https://vasp.at/wiki/MAGMOM),
[`LORBIT`](https://vasp.at/wiki/LORBIT),
[`LASPH`](https://vasp.at/wiki/LASPH), and the element-dependent
[`LMAXMIX`](https://vasp.at/wiki/LMAXMIX) restart setting. Backend scope
follows the VASP
[`WAVECAR`](https://vasp.at/wiki/WAVECAR) and
[`PAW formalism`](https://vasp.at/wiki/Projector-augmented-wave_formalism)
documentation. Functional-specific restart constraints are summarized in the
official VASP
[`meta-GGA workflow`](https://vasp.at/wiki/Band-structure_calculation_using_meta-GGA_functionals).

An independent Wilson-loop calculation may be used as an optional cross-check;
it is not part of `-z2`. See [Yu et al., *Phys. Rev. B* **84**, 075119 (2011)](https://doi.org/10.1103/PhysRevB.84.075119).
