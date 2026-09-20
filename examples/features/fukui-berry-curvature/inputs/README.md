# VASP input for the MoS₂ 12 × 12 curvature example

These files generate a complete spinor WAVECAR for VASPBERRY from the public
MoS₂ SCF charge density. The setup was executed with VASP 5.4.4 and then passed
to the current native Fukui routine. The generated WAVECAR is about 149 MB and
is not bundled in the repository.

## Prepare the calculation

The preparation command copies `INCAR`, `KPOINTS` and `POSCAR`, decompresses
[`CHGCAR.gz`](../../../1H-MoS2/KPATH/1.scf/CHGCAR.gz), and copies your licensed
Mo/S PAW-PBE file. It checks that the potential datasets match the
[material specification](../../../1H-MoS2/PSEUDOPOTENTIAL.md).

From the repository root:

```bash
python3 examples/features/fukui-berry-curvature/prepare_vasp.py \
  --potcar /path/to/licensed/Mo-S/POTCAR \
  --output-dir results/mos2-fullmesh-vasp
```

The generated directory contains:

| File | Purpose |
|---|---|
| `POSCAR` | Public 1H-MoS₂ structure; in-plane lattice length 3.171634303 Å, cell height 15 Å. |
| `POTCAR` | Mo and S PAW-PBE datasets supplied by the user. |
| `CHGCAR` | Stored public SOC SCF density and on-site information. |
| `INCAR` | Fixed-charge SOC calculation: `ICHARG=11`, `ENCUT=400`, `NBANDS=26`, `ISYM=-1`, `LWAVE=.TRUE.`. |
| `KPOINTS` | Complete Γ-centered 12 × 12 × 1 mesh. |
| `input_manifest.json` | Machine-readable preparation record. |

`PREC=Normal`, `LASPH=.FALSE.`, `GGA_COMPAT=.FALSE.` and `LMAXMIX=2` follow
the supplied charge-density calculation. This is an NSCF continuation of that
public calculation, not a new self-consistent ground-state optimization.
The actual VASP outputs are retained in the [reference directory](../reference/vasp/).

## Run your VASP executable

Use the licensed noncollinear/SOC executable configured for your machine:

```bash
cd results/mos2-fullmesh-vasp
/path/to/vasp_ncl > vasp.stdout.log 2> vasp.stderr.log
cd ../..
```

On an MPI installation, use the corresponding launcher and processor count,
for example `mpirun -np 4 /path/to/vasp_ncl`. The supplied inputs do not fix
site-specific parallelization settings.

Check that the electronic iterations converge and the VASP run finishes
normally. The expected WAVECAR has **144 k points, 26 spinor bands, 18 occupied
bands and a sampled gap of about 1.67355 eV**. The demonstrated one-process run
took about **144 seconds** and used a peak **602 MB** of resident memory on an
Intel Core i9 development machine. These are measured example costs, not
performance guarantees for another machine.

The reference binary used a process-local 64 MiB stack limit and one thread
per numerical library. Apply the runtime settings appropriate to your local
VASP build; the executable and licensed potentials are not distributed here.

## Run VASPBERRY on the resulting file

```bash
make serial
python3 examples/features/fukui-berry-curvature/run.py \
  --wavecar results/mos2-fullmesh-vasp/WAVECAR \
  --output-dir results/mos2-fukui
```

The [main tutorial](../README.md) gives the equivalent direct Fortran command,
output columns, physical interpretation and plotting instructions.

## Reference calculation and convergence

The recorded WAVECAR is 149,002,560 bytes. Its calculation record is
available in [result.json](../reference/result.json).

Different compilers, VASP versions and diagonalization paths can produce
different wavefunction phases and numerical representations. Compare the occupied-subspace
curvature and gap, rather than requiring coefficient-by-coefficient equality.
The 12 × 12 mesh is a reproducible tutorial setting. Converge the underlying
SCF calculation and k mesh separately when preparing quantitative material
results.
