# Real graphene SOC: native projected-spin Chern example

This example starts from ordinary VASP calculations of planar graphene. It uses
unscaled intrinsic spin–orbit coupling, an eight-band occupied spinor projector,
and the native VASPBERRY `spin-chern` task. No artificial SOC multiplier is used. The saved reference has `C_plus=1`, `C_minus=-1`, and
`C_spin=(C_plus-C_minus)/2=1`; its charge Chern number is zero to roundoff.

The gap is approximately **0.912 micro-eV for the specified calculation**.
This is a numerical example for the recorded PBE carbon potential and geometry,
not a claim that graphene's physical gap has converged with respect to its
pseudopotential, cutoff, charge-density mesh, or exchange-correlation treatment.
In particular, the classic [Gmitra et al. calculation](https://arxiv.org/abs/0904.3315)
reported about 24 micro-eV with higher-orbital contributions. Do not replace the
actual eigenvalues with that literature value.

## 1. Build and define local paths

Run the commands below from the repository root. Build the current source with
Intel oneAPI and Intel MPI, following [BUILD.md](../../../docs/BUILD.md). Adapt
the oneAPI initialization path to your installation:

```bash
source /opt/intel/oneapi/setvars.sh
make ifx-mpi
EXAMPLE="$PWD/examples/materials/graphene-spin-chern"
WORK="$PWD/results/graphene-spin-chern"
VB="$PWD/build/vaspberry-ifx-mpi"
VASP_NCL=/absolute/path/to/your/licensed/vasp_ncl
C_POTCAR=/absolute/path/to/your/licensed/PAW-PBE/C/POTCAR
mkdir -p "$WORK"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
```

The preparation helper uses only the Python standard library. Plotting uses
Matplotlib in your usual Python environment; see `requirements-transport.txt`
if needed. A GNU serial alternative is `make serial`
with `VB="$PWD/build/vaspberry"` and the MPI launcher omitted. The saved local
reference was cross-checked with GNU serial/MPI and Intel ifort serial builds;
Intel ifx/MPI execution was not available on the reference macOS host.

Use the exact potential identity in [pseudopotential.json](pseudopotential.json).
The preparation helper checks its SHA-256. The public example contains neither
POTCAR nor VASP executable/source. The reference used ordinary serial VASP 5.4.4
on macOS with a 64 MiB process stack (`ulimit -s 65536` in the VASP subshells);
choose the stack setting required by your own installation.

## 2. Obtain the self-consistent SOC charge density

```sh
python3 "$EXAMPLE/prepare_vasp.py" --stage scf --mesh 9 \
  --potcar "$C_POTCAR" --output-dir "$WORK/scf"
(cd "$WORK/scf" && ulimit -s 65536 && "$VASP_NCL" > vasp.stdout 2> vasp.stderr)
```

The fixed geometry has `a=2.46 Å`, `c=20 Å`, and two equivalent carbon sites.
The SCF settings include `ENCUT=520`, `EDIFF=1e-11`, `LSORBIT=.TRUE.`,
`ISYM=2`, `NBANDS=16`, `MAGMOM=6*0`, and `SIGMA=0.001 eV`.
Check normal termination and the `EDIFF is reached` line in `OUTCAR`.
The next preparation step checks those conditions and verifies the recorded
SCF input hashes before accepting its density.

## 3. Generate a full periodic SOC mesh

```sh
python3 "$EXAMPLE/prepare_vasp.py" --stage mesh --mesh 12 \
  --potcar "$C_POTCAR" --scf-dir "$WORK/scf" --output-dir "$WORK/mesh12"
(cd "$WORK/mesh12" && ulimit -s 65536 && "$VASP_NCL" > vasp.stdout 2> vasp.stderr)
```

This is an explicit 12×12×1 list of all 144 points, `ICHARG=11`, `ISYM=-1`,
`EDIFF=1e-9`, and at least 20 electronic iterations. A mesh dimension divisible
by three includes both K=(1/3,1/3) and K′=(2/3,2/3). Use the same prepared SCF
density for the 6×6, 9×9, and 12×12 comparisons.

The helper restores the exact *generating* POSCAR header in the copied CHGCAR.
VASP 5.4.4 writes this header to six decimal places; at a sub-micro-eV scale the
rounded coordinates and cell gave a spurious approximately 5.25 micro-eV
no-SOC splitting in the initial diagnostic. The helper accepts only exact or
six-decimal-rounded header values matching the known source structure and
leaves every density/augmentation byte unchanged. Source/output/payload hashes
are saved in `charge-header-repair.json`. This operation does not modify the
SCF density and must not be applied to a different structure.

## 4. Run VASPBERRY through MPI

```sh
(cd "$WORK/mesh12" && mpiexec -n 4 "$VB" --task spin-chern --bands 1:8 --mesh 12,12)
```

`WAVECAR`, `OUTCAR`, Cartesian spin axis `z`, and output prefix `SPIN` are the
defaults. An equivalent explicit command is:

```sh
mpiexec -n 4 "$VB" --task spin-chern --wavecar "$WORK/mesh12/WAVECAR" \
  --outcar "$WORK/mesh12/OUTCAR" --bands 1:8 --mesh 12,12 \
  --spin-axis z --energy-gap-tol 1e-8 --spin-gap-tol 1e-6 \
  --output "$WORK/mesh12/EXPLICIT"
```

| Native output | Meaning |
|---|---|
| `SPIN_CHERN.csv` | Sector Chern numbers, their sum, difference, and half-difference |
| `SPIN_BERRY.csv` | Sector/parent plaquette flux in radians, cell area, and curvature in Å² |
| `SPIN_SPECTRUM.csv` | Projected Pauli eigenvalues, sector ranks, and energy/metric diagnostics |

All three files have `result_status=PASS` only after native validation. They can
be read directly by pandas, NumPy, Julia, MATLAB, or a plotting program. The
Fortran code calculates the projector, generalized projected-spin eigenproblem,
periodic overlap links, and Chern sums. Python below only reads completed CSVs
and draws them.

The source Fermi smearing partly occupies both sides of the tiny gap. `1:8`
selects the isolated zero-temperature energy projector, and the metadata
explicitly records that source occupations do not match an integer projector.
The implementation uses the WAVECAR pseudo-wavefunction Gram metric;
`physical_paw_validated=false` remains explicit. The result is a projected-spin
invariant, not a Kubo spin Hall conductivity, a finite-temperature occupation
integral, or a full PAW spin-operator validation. See [SPIN_CHERN.md](../../../docs/SPIN_CHERN.md).

## 5. Plot raw native outputs and optional real VASP bands

A compact plot needs only the three native CSVs:

```sh
python3 "$EXAMPLE/plot.py" --input-dir "$WORK/mesh12" \
  --output-dir "$WORK/sector-plots"
```

To add an actual Γ–K–M–Γ band path and show the small K-point gap:

```sh
python3 "$EXAMPLE/prepare_vasp.py" --stage path \
  --potcar "$C_POTCAR" --scf-dir "$WORK/scf" --output-dir "$WORK/path"
(cd "$WORK/path" && ulimit -s 65536 && "$VASP_NCL" > vasp.stdout 2> vasp.stderr)
python3 "$EXAMPLE/plot.py" --input-dir "$WORK/mesh12" \
  --path-wavecar "$WORK/path/WAVECAR" \
  --path-node-indices 1 17 33 49 --path-labels Gamma K M Gamma \
  --gap-k-index 17 --output-dir "$WORK/figure"
```

The plot exports unchanged full-precision WAVECAR eigenvalues to `bands.csv`
and shifts the display origin to the supplied path's midgap. OUTCAR's usual
four-decimal band printing cannot resolve this gap. Sector maps show constant
native plaquette fluxes clipped to the first Brillouin zone, without smoothing
or interpolation. These coarse meshes establish the integer result here;
they do not resolve the microscopic width or peak height of Dirac curvature.

To reproduce the published plot without a VASP license or fresh calculation:

```sh
python3 "$EXAMPLE/plot.py" --input-dir "$EXAMPLE/reference/native-12x12" \
  --bands-csv "$EXAMPLE/reference/bands.csv" \
  --path-node-indices 1 17 33 49 --path-labels Gamma K M Gamma \
  --gap-k-index 17 --output-dir "$WORK/reference-figure"
```

## 6. Repeat the convergence and gap-closed control

| Fixed-charge mesh | EDIFF (eV) | Minimum direct gap (µeV) | C+ | C− | Cspin |
|---|---:|---:|---:|---:|---:|
| 6×6 | 1e-9 | 0.9120148 | 1 | −1 | 1 |
| 9×9 | 1e-9 | 0.9119943 | 1 | −1 | 1 |
| 12×12 | 1e-9 | 0.9120028 | 1 | −1 | 1 |
| 6×6, tighter EDIFF | 1e-10 | 0.9120202 | 1 | −1 | 1 |
| 6×6, SOC off | 1e-9 | 0.0001940 | gap-closed rejection | — | — |

For topology, repeat steps 3–4 with `--mesh 6`/`--mesh 6,6`, then 9/9,9 and
12/12,12. Tighten the 6×6 fixed-charge calculation to `--ediff 1e-10` in a new
folder. Saved numerical results and source hashes are in
[reference/summary.json](reference/summary.json) and
[reference/provenance.json](reference/provenance.json).

A matched spinor calculation without SOC should fail the energy-gap guard:

```sh
python3 "$EXAMPLE/prepare_vasp.py" --stage mesh --mesh 6 --no-soc \
  --potcar "$C_POTCAR" --scf-dir "$WORK/scf" --output-dir "$WORK/no-soc6"
(cd "$WORK/no-soc6" && ulimit -s 65536 && "$VASP_NCL" > vasp.stdout 2> vasp.stderr)
(cd "$WORK/no-soc6" && mpiexec -n 4 "$VB" --task spin-chern --bands 1:8 --mesh 6,6)
```

A nonzero exit and no completed invariant are the expected outcome for that
control. Do not lower the energy-gap tolerance to assign topology to a
numerically gapless state. The saved control, tighter electronic tolerance,
and mesh changes are separate checks; none substitutes for converging the
material Hamiltonian when using this workflow for a research prediction.

## Local projected-spin curvature

The [native spin-sector Kubo example](kubo/README.md) adds real path curves and
raw full-mesh differential-proxy sums beside these Fukui invariants. Its CSVs
retain the approximation and sampling limits explicitly.
