# Reproduce the Bi spin Hall convergence study

This supplement preserves the actual fixed-charge VASP inputs used for the
mesh and source-band study. The ordinary 48-band quickstart remains the
simplest first calculation. This study additionally distinguishes the number
of eigenstates calculated by VASP from the number retained by VASPBERRY.

## Calculation sequence

| Case | Mesh | Stored bands | Bands retained in the comparison | Initialization | Minimum electronic iterations |
|---|---:|---:|---:|---|---:|
| `mesh6-source48` | 6×6 | 48 | 48 | Cold | 20 |
| `mesh12-source48-continuation` | 12×12 | 48 | 48 | Same-k 48-band checkpoint | 8, after 9 cold iterations |
| `mesh6-source64` | 6×6 | 64 | 48 | Same-k converged 48-band WAVECAR | 8 |
| `mesh12-source64` | 12×12 | 64 | 48 | Same-k converged 48-band WAVECAR | 8 |
| `mesh6-source80` | 6×6 | 80 | 48, 56 and 64 | Same-k converged 64-band WAVECAR | 8 |
| `mesh18-source64` | 18×18 | 64 | 48 | Cold, no input WAVECAR | 12 |

Every case uses the same new SCF charge density, structure, potential, 400 eV
cutoff and physical settings. The occupied bundle comprises bands 1–10.
`EDIFF = 1E-8` must be reached in addition to the listed `NELMIN`. A minimum
iteration count alone does not establish eigenstate or transport convergence.

The top states of a finite VASP calculation can be less accurate than lower
states. The 64-band runs therefore supply the retained 48-band mesh study;
the 80-band runs supply the retained 64-band closure check. Increasing
`--source-band-limit` changes both the spin-current product and its
intermediate-state sum. It does not improve the stored VASP eigenstates.
Increasing VASP `NBANDS` and increasing the retained limit are separate tests.

The `inputs/` directories contain byte-for-byte copies of the actual INCAR
and KPOINTS files. Their strided partitions use qx as the inner index and
signed fractional coordinates. The quickstart preparation helper instead
uses contiguous partitions and qy as the inner index. Both cover the same
full mesh, but their individual WAVECAR files must not be interchanged.

## 1. Prepare the common licensed inputs

Run commands from the VASPBERRY repository root. Set `PACK` to this supplement
and `WORK` to a new working directory. Build the opt-in producer by following
`tools/vasp544_spin_bridge/README.md` first.

```bash
PACK=examples/materials/bi-spin-hall/inputs/convergence
WORK=results/bi-study
VASP_BIN=/path/to/isolated-vasp544/bin/vasp_ncl
PRODUCER_MANIFEST=/path/to/isolated-vasp544/vaspberry-spin-producer.json

# This preparation restores the supplied density and validates the potential.
# The charge-source directory is only a source of these common inputs.
python3 examples/materials/bi-spin-hall/prepare_vasp.py \
  --stage spin --potcar /path/to/matching/Bi/POTCAR \
  --mesh 6 6 --nbands 48 --output-dir "$WORK/charge-source"
```

The exact charge, structure and potential identities are recorded in
`manifest.json`. No POTCAR, executable, VASP implementation, or WAVECAR is
included here. To regenerate the SCF density, use the material example's
SCF preparation procedure. All subsequent chunks must share that completed
SCF calculation. A newly generated density need not be byte-identical to the
archived density, and its results form a new calculation series.

For each case, create a new directory and copy its two supplied text inputs
and the common files. For example:

```bash
mkdir -p "$WORK/mesh6-source48"
for part in 00 01; do
  dest="$WORK/mesh6-source48/chunk$part"
  mkdir "$dest"
  cp "$PACK/inputs/mesh6-source48/chunk$part/INCAR" "$dest/INCAR"
  cp "$PACK/inputs/mesh6-source48/chunk$part/KPOINTS" "$dest/KPOINTS"
  cp "$PACK/POSCAR" "$dest/POSCAR"
  cp "$WORK/charge-source/CHGCAR" "$dest/CHGCAR"
  cp "$WORK/charge-source/POTCAR" "$dest/POTCAR"
  python3 tools/run_vasp_spin_producer.py \
    --run-dir "$dest" --binary "$VASP_BIN" \
    --producer-manifest "$PRODUCER_MANIFEST" \
    --timeout 3600 --memory-limit-mib 4096
done
```

This loop is sequential. Independent one-thread chunks can be run concurrently
within available CPU and memory limits. Keep a record of every attempt and
require successful electronic convergence, normal completion and both matrix
exports before using a result. A zero process exit code alone is insufficient.

## 2. Enlarge the same-k source space

To prepare `mesh6-source64/chunk00`, copy its supplied INCAR and KPOINTS and
the common files into a new directory as above. Also copy **the completed
`mesh6-source48/chunk00/WAVECAR`** to that directory. Repeat for chunk01.
The supplied INCAR sets `ISTART=1`, `NBANDS=64` and `NELMIN=8`. VASP reads the
48 stored states and initializes additional states itself. Check that the
output confirms a successful WAVECAR read and the expected old band count.

For `mesh6-source80`, use the corresponding completed `mesh6-source64`
WAVECAR. The supplied settings request 80 bands and eight minimum warm
iterations. The physical input files and ordered KPOINTS must be identical
between a restart and its source. Do not copy coefficients to a different
k-point list, reorder records informally, or treat a cropped 48-band matrix
as a new 64-band VASP calculation.

### Preparing the 12×12 source

For a fresh reproduction, run all four supplied
`mesh12-source48-checkpoint/chunk00` through `chunk03` inputs **to normal
electronic convergence** with their default `NELMIN=20`. These are cold
48-band inputs. Then prepare each `mesh12-source64` directory with its
matching completed 48-band WAVECAR and the supplied warm 64-band inputs.

The archived calculation used a resource checkpoint instead: the original
cold 48-band attempts were gracefully stopped after nine iterations; each
was then continued with eight minimum warm iterations using
`mesh12-source48-continuation`. Only the converged continuation was accepted
as the source of the 64-band run. `execution.json` preserves both stages and
marks the initial attempts as checkpoint-only. A fresh uninterrupted 20-step
source follows the same physical setup but need not reproduce this historical
initialization path bit for bit. Do not deliberately interrupt a new run at
a wall-clock time to imitate the archived checkpoint.

## 3. Calculate the new 18×18 mesh

Prepare the nine `mesh18-source64/chunk00` through `chunk08` directories using
their supplied inputs and the common charge, potential and structure. Run
each with the same producer command. These are **cold** runs:
`ISTART=0`, `NBANDS=64`, `NELMIN=12`, and no input WAVECAR. The smaller
`NELMIN` differs intentionally from the quickstart's 20; the actual
electronic residuals and independent retained-state comparisons are recorded
alongside the calculation. A WAVECAR from the 12×12 mesh is not used.

## 4. Merge and evaluate the response

Once all nine 18×18 chunks have passed the producer checks:

```bash
python3 tools/vaspberry_kubo.py spin-merge \
  --run-dirs "$WORK"/mesh18-source64/chunk0? \
  --mesh 18 18 --occupied 10 --output-dir "$WORK/assembled18"

python3 tools/vaspberry_kubo.py spin-hall \
  --matrices "$WORK/assembled18/physical-matrices.npz" \
  --metadata "$WORK/assembled18/physical-matrices.json" \
  --mesh 18 18 --occupied 10 --source-band-limit 48 \
  --formats csv dat npz --output-dir "$WORK/hall18-cap48"
```

Use the corresponding mesh and source directories for the 6×6 and 12×12
series. In the 80-band 6×6 series, compare limits 48, 56 and 64 while holding
the underlying matrix files fixed. The merge command checks the source
Hamiltonian identities and complete mesh union, and preserves the original
per-k eigenvector gauge. The current spin Hall integration is an insulating,
zero-temperature occupied-bundle calculation.

The reported conventional spin Hall unit is `(ħ/e)(e²/h)` for the sheet.
Spin–orbit coupling need not conserve the conventional spin current, so this
response need not be an integer. The source-band product approximation,
eigenstate accuracy and k-mesh integration error must each be assessed.

## Execution record

`execution.json` contains the observed per-chunk wall time, sampled process
RSS, electronic iteration count, final residuals, exit status and completion
checks. It identifies running or checkpoint-only cases explicitly. These
records do not label a material conductivity converged merely because the
underlying VASP runs completed.

The measured workstation had 10 physical cores and 64 GiB RAM. Each VASP
process used one scientific thread. The target was six concurrent VASP
processes; the actual peak was seven. Timings under this mixed load are
observations, not performance guarantees. The execution record preserves the
complete resource history.
