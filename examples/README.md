# Reproduce VASPBERRY calculations from VASP files

These tutorials start with **actual public VASP WAVECAR files**, run the current
VASPBERRY programs, and compare the new output with committed numerical results
and figures. Each tutorial shows the underlying VASPBERRY command, what its
options mean, how to plot the result, and what to change for your own material.

## 1. Choose a calculation

| Tutorial | Actual VASP input | Result to reproduce |
|---|---|---|
| [Berry curvature / Chern](features/fukui-chern/) | Bi, full 12×12 SOC mesh, occupied bands 1–10 | Occupied-subspace Chern = 0 and native curvature map |
| [Z₂](features/z2/) | The same full-mesh Bi WAVECAR | Z₂ = 1, matching half-zone parities and n-field figure |
| [Kubo Berry curvature](features/kubo-curvature/) | MoS₂, SOC K–Γ–K′ path | Canonical-momentum curvature versus k; unresolved band degeneracies are marked |
| [Charge Hall in an insulating gap](features/hall-valley/) | Bi, full mesh and occupied-subspace Fukui flux | Charge Hall = 0 in the verified gap, plus diagnostics |
| [Circular optical response](features/circular-dichroism/) | MoS₂, the same actual SOC band-path WAVECAR | Left/right optical response and meaningful selectivity |
| [Gamma wavefunction](features/wavefunction/) | MoS₂ WAVECAR with its matching POSCAR/EIGENVAL | A selected actual spinor state and real-space figure |

Bi is time-reversal symmetric: a zero charge Chern/Hall result is the expected
reference, while its Z₂ invariant is nontrivial. The Hall tutorial uses the
occupied-subspace Fukui route because individual Bi bands have unresolved
Kramers degeneracies. It does not present those bands as valid point-Kubo input.
The MoS₂ line is useful for k-resolved plots; it cannot provide a full-BZ Hall
integral or Chern number.

## 2. Prepare the real input files

The MoS₂ band-path WAVECAR, POSCAR and EIGENVAL are already in
[`1H-MoS2/KPATH/2.band/`](1H-MoS2/KPATH/2.band/). The Bi WAVECAR is stored with
Git LFS. Obtain its actual payload using either:

```bash
git lfs pull --include='examples/Bi_Z2/WAVECAR'
```

or the public, checksum-pinned download helper (Git LFS is not required):

```bash
python3 examples/fetch_inputs.py bi --output-dir results/inputs/bi
```

The latter creates `results/inputs/bi/WAVECAR` (200,421,600 bytes) and a download
record. A tiny Git LFS pointer is not the input payload. Missing inputs stop
the calculation; the runners never substitute an analytic model or an old
output. Reference comparisons identify the input by checksum. See
[input files and provenance](INPUTS.md).

## 3. Run VASPBERRY and compare the reference

Build once and install the plotting/analysis dependencies:

```bash
make serial
python3 -m pip install -r requirements-transport.txt
```

Open the selected tutorial above and follow its explicit production commands.
For example, the Bi Z₂ calculation itself is:

```bash
mkdir -p results/bi-z2-manual
(
cd results/bi-z2-manual
../../build/vaspberry-gfortran -f ../inputs/bi/WAVECAR -o NFIELD \
  -z2 1 -kx 12 -ky 12 -s 2 -ii 1 -if 10
)
```

It writes `NFIELD.dat` and `Z2_FIELD.csv`. The
[Z₂ tutorial](features/z2/) explains their required PASS status, expected
values, plotting and comparison. Its `run.py` performs the same steps and
records the commands, input hashes and checks:

```bash
# From the repository root:
python3 examples/features/z2/run.py \
  --wavecar results/inputs/bi/WAVECAR --output-dir results/bi-z2
```

Each feature's `reference/` contains original VASPBERRY numerical output (or a
lossless compressed copy), a comparison summary, provenance, and a figure
produced from that output. These are **fresh calculations using the stated
real VASP input**, not illustrative model curves. See its README for exact
comparison tolerances and for formatted-output precision limits.

### Run the native calculation with MPI

For native Fortran calculations, build the MPI executable and use the same
physical options. For example, the actual MoS₂ Kubo calculation with two ranks:

```bash
make mpi
mkdir -p results/mos2-kubo-mpi
(
cd results/mos2-kubo-mpi
mpiexec -np 2 ../../build/vaspberry-mpi \
  -f ../../examples/1H-MoS2/KPATH/2.band/WAVECAR \
  -s 2 -kubo 2 -ii 17 -if 18 -kubo_csv KUBO.csv -o BERRYCURV
)
```

Follow the [Kubo tutorial](features/kubo-curvature/) for output interpretation
and the unresolved-state mask. Python tutorial wrappers, plotting and the
Python occupied-subspace Hall path are serial programs; do not launch a
wrapper once per MPI rank. See the [build guide](../docs/BUILD.md).

## 4. Apply it to your own system

Follow [applying a tutorial to your material](APPLY_TO_YOUR_SYSTEM.md).
Replace the input path in the underlying production command and set the mesh,
spinor setting, bands and energy window from your own VASP calculation.
The Bi and Kubo tutorial wrappers require their named reference input. The
optical and wavefunction wrappers also accept compatible user data, with the
checks and limits stated in their README. For general use, adapt the documented
production command; do not change a reference checksum to bypass a check.

## Optional: run the complete reference set

After fetching Bi through the helper:

```bash
python3 examples/run_examples.py --list
python3 examples/run_examples.py --all \
  --bi-wavecar results/inputs/bi/WAVECAR --output-dir results/all-examples
```

With a full Git LFS checkout, omit `--bi-wavecar`. Each output directory must
be new. The batch runner records commands, elapsed times, exit codes and logs;
it also requires a matching successful result manifest. No real calculation
is silently replaced with a stored-field check. CI runs this same real-input
reference set and retains the outputs as an artifact.

## Materials, plots and developer checks

- [Material catalog](materials/): retained MoS₂/Bi files and preparation notes.
- [Standalone Hall CSV plotting](kubo/): plots standardized transport output.
- [Developer numerical checks](../validation/models/): analytic models,
  synthetic states and historical field checks, kept outside user tutorials.
- [Supplementary reference map](../docs/REFERENCE_MATERIALS.md): connect these
  examples and method checks for hands-on documents and release technical reports.

VASP `POTCAR` files are not distributed. The supplied outputs reproduce
VASPBERRY postprocessing; consult each material's provenance before claiming
an end-to-end reproduction of its original VASP SCF/NSCF calculation.
