# Bi 12 x 12 Fukui-Hatsugai Z2 example

This example demonstrates the official VASPBERRY `-z2 1` path on a buckled
honeycomb Bi test structure. It evaluates the two-dimensional Z2 invariant by
the lattice n-field method.

## What is reproducible here

The tracked 12 x 12 x 1 SOC spinor [`WAVECAR`](WAVECAR) is a post-processing
fixture. It has 144 k points, 18 spinor bands, and SHA-256
`a8d81854f2efc561938e478dde1be29a17ccc95d5d37325c122b9d9e82fa0838`.
Bands 1-10 are used as the occupied subspace.

The archived VASP files are in [`archive-2016-run/`](archive-2016-run/).
`OUTCAR` and `vasprun.xml` show a VASP 5.4.1 fixed-charge run from 2016 with
`ICHARG=11`, `PREC=Normal`, `EDIFF=1e-6`, `ISPIN=1`, `LSORBIT=.TRUE.`,
`ISYM=-1`, and `LREAL=.FALSE.`. The complete inputs and output of the earlier
SCF run that produced its input `CHGCAR` are unavailable. The archive is
therefore sufficient to reproduce VASPBERRY post-processing, but not an
end-to-end DFT calculation. `archive-2016-run/INCAR.unverified` is retained as
an old tracked file and is explicitly not asserted to be the input that made
the archived `WAVECAR`.

The recommended, internally consistent input templates for a new calculation
are under [`inputs/`](inputs/). They are not claimed to reproduce every
numerical detail of the 2016 run.

## Files used by each stage

| File | Needed by | Role |
|---|---|---|
| `POSCAR`, local `POTCAR`, `INCAR`, `KPOINTS` | VASP SCF | produce a converged `CHGCAR` |
| converged `CHGCAR` plus the same structure and potentials | VASP fixed-charge run | produce the full-mesh SOC `WAVECAR` |
| `WAVECAR` | VASPBERRY | calculate links, n-field, and Z2 |
| `OUTCAR`, `EIGENVAL` | user verification | confirm settings, band count, and occupied-unoccupied separation |

The licensed Bi PAW data are excluded. Non-proprietary provenance is recorded
in [`PSEUDOPOTENTIAL.md`](PSEUDOPOTENTIAL.md).

## Reproduce the archived 12 x 12 post-processing result

Install Git LFS, GNU Fortran, and LP64 BLAS/LAPACK. From the repository root:

```bash
git lfs pull --include='examples/Bi_Z2/WAVECAR'
./examples/Bi_Z2/scripts/run_z2.sh
```

The script compiles the serial GNU source in a temporary directory, checks the
`WAVECAR` size and checksum, and uses the same options as this command with the
standard Makefile product:

```bash
./build/vaspberry-gfortran \
  -f examples/Bi_Z2/WAVECAR \
  -o NFIELD -z2 1 \
  -kx 12 -ky 12 -s 2 -ii 1 -if 10
```

It writes `NFIELD.dat`, `Z2_FIELD.csv`, and `fortran.log` under
`examples/Bi_Z2/results-z2/`. The same calculation can be run with the MPI
executable; see [`../../docs/BUILD.md`](../../docs/BUILD.md) and the main
[`Z2 guide`](../../docs/Z2_FUKUI_HATSUGAI.md).

## Read the result

The corrected reference is stored under
[`reference-v1.1.1-12x12/`](reference-v1.1.1-12x12/). Its half-zone sums are

```text
top    = -3 -> parity 1
bottom = +3 -> parity 1
```

That directory is a v1.1.1 numerical snapshot and retains its historical
candidate-policy strings. The current runner recalculates the same field with
the v1.2.0 result schema under `results-z2/`.

Thus the two complementary half-zone evaluations agree and the n-field method
gives `Z2=1` on this 12 x 12 mesh. The numerical self-consistency diagnostics
also pass: the total Chern residual is `1.159e-14`, the maximum wrapped
time-reversal-odd flux residual is `1.815e-14 rad`, and the minimum link
singular value is `0.815455`.

On the archived 144-point `EIGENVAL`, the minimum direct band-10/band-11
separation is `0.592449 eV` and the sampled global gap is `0.510045 eV`.
These values describe this fixed mesh and do not establish gap or k-mesh
convergence.

This does not replace the physical input checks or a mesh-convergence study.
Verify that the material is nonmagnetic and time-reversal symmetric, that
bands 1-10 are separated from band 11 throughout the Brillouin zone, and that
the common parity is stable on denser even Gamma-centered meshes.

The files in [`legacy-pre-1.1.1/`](legacy-pre-1.1.1/) were produced by the
incomplete implementation. Their top and bottom parities disagree and they
must not be used as a Z2 result.

## Draw the 12 x 12 field

Install NumPy and Matplotlib from the repository root if they are not already
available:

```bash
python3 -m pip install -r requirements-transport.txt
```

```bash
python3 examples/Bi_Z2/scripts/plot_nfield.py \
  examples/Bi_Z2/reference-v1.1.1-12x12/Z2_FIELD.csv \
  --legacy examples/Bi_Z2/legacy-pre-1.1.1/NFIELD.dat \
  --output examples/Bi_Z2/Z2_nfield_12x12.pdf
```

Rendered reference: [`Z2_nfield_12x12.pdf`](Z2_nfield_12x12.pdf) and
[`Z2_nfield_12x12.png`](Z2_nfield_12x12.png).

The pointwise n-field is gauge- and logarithm-branch-dependent. The plot is a
discrete regression view; Z2 is determined by the half-zone integer sum
modulo two, not by requiring every point to look mirror-symmetric.

## New VASP calculation

1. Copy `inputs/POSCAR`, `inputs/01_scf/INCAR`, and
   `inputs/01_scf/KPOINTS` to a clean VASP directory. Add a locally licensed
   `POTCAR` in the `Bi` species order and run the site's SOC-capable VASP
   executable until `CHGCAR` is converged.
2. Copy the converged `CHGCAR` with the same `POSCAR` and `POTCAR` into a new
   directory containing `inputs/02_z2_nscf/INCAR` and `KPOINTS`. This is an
   `ICHARG=11` fixed-charge run, not a second SCF run; it writes the full-mesh
   spinor `WAVECAR` consumed by VASPBERRY.
3. Confirm `NKPTS=144`, an even occupied rank, `NBANDS>Nocc`, a direct
   separation between the selected occupied bundle and the next band at every
   k point, and a positive sampled global gap. Check both on denser sampling.
4. Run VASPBERRY with the command above and repeat on denser even meshes.

Method references and the exact result contract are given in
[`docs/Z2_FUKUI_HATSUGAI.md`](../../docs/Z2_FUKUI_HATSUGAI.md).
