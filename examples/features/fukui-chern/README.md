# Fukui Berry curvature and Chern number from a real Bi WAVECAR

**Start with VASP's WAVECAR, run VASPBERRY, and compare the output with the
stored reference.** This example uses the complete occupied band bundle
1–10 of buckled honeycomb Bi. Its total Chern number is **0**; the same
material has a nontrivial Z₂ invariant in the [next tutorial](../z2/).

## 1. Obtain the actual VASP output

Run these commands from the repository root. The calculation input is the
**real Bi `WAVECAR`**, not a JSON configuration or model Hamiltonian.

```bash
git lfs pull --include='examples/Bi_Z2/WAVECAR'
shasum -a 256 examples/Bi_Z2/WAVECAR
```

Expected SHA-256:

```text
a8d81854f2efc561938e478dde1be29a17ccc95d5d37325c122b9d9e82fa0838
```

The payload is **200,421,600 bytes**. A 134-byte Git LFS pointer is not a usable
WAVECAR. If Git LFS is unavailable, download the same public file separately:

```bash
mkdir -p inputs
curl -fL https://media.githubusercontent.com/media/Infant83/VASPBERRY/a692e21482d24767859d02b7885dd70b234c6a2c/examples/Bi_Z2/WAVECAR -o inputs/Bi-WAVECAR
shasum -a 256 inputs/Bi-WAVECAR
```

Then pass `--wavecar inputs/Bi-WAVECAR` to the runner below. It verifies the
exact size and checksum before executing VASPBERRY. It never substitutes a
model or a stored result for a missing input.

### Input files and provenance

| File | Role |
|---|---|
| [`Bi_Z2/WAVECAR`](../../Bi_Z2/WAVECAR) | Actual VASPBERRY input: 144 points on a 12 × 12 × 1 mesh, 18 SOC spinor bands |
| [`archive-2016-run/EIGENVAL`](../../Bi_Z2/archive-2016-run/EIGENVAL) | VASP band energies, cross-checked against WAVECAR by the runner |
| [`archive-2016-run/OUTCAR`](../../Bi_Z2/archive-2016-run/OUTCAR) | Archived VASP 5.4.1 settings and output |
| [`inputs/POSCAR`](../../Bi_Z2/inputs/POSCAR) | Bi structure for a new calculation |
| [`inputs/01_scf/`](../../Bi_Z2/inputs/01_scf/) and [`inputs/02_z2_nscf/`](../../Bi_Z2/inputs/02_z2_nscf/) | Recommended SCF and full-mesh SOC input templates |
| [`PSEUDOPOTENTIAL.md`](../../Bi_Z2/PSEUDOPOTENTIAL.md) | Potential provenance; licensed POTCAR is not redistributed |

The archived WAVECAR came from a 2016 fixed-charge calculation. The full
preceding SCF provenance is unavailable. This tutorial reproduces **VASPBERRY
post-processing of that public VASP result**. The recommended VASP input
templates are not claimed to recreate every digit of the 2016 WAVECAR.

## 2. Build VASPBERRY and plotting dependencies

```bash
make serial
python3 -m pip install -r requirements-transport.txt
```

## 3. Run the actual calculation and make the figure

```bash
python3 examples/features/fukui-chern/run.py --output-dir results/bi-fukui
```

The wrapper checks the VASP input, runs the executable, checks the output,
and draws the figure. It runs precisely this native calculation inside a
new output directory (where `WAVECAR` points to the downloaded input):

```bash
/path/to/vaspberry/build/vaspberry-gfortran \
  -f WAVECAR -o BERRYCURV \
  -kx 12 -ky 12 -s 2 -ii 1 -if 10 > fortran.log 2>&1
```

`-kx/-ky` describe the **mesh already present in WAVECAR**; changing them
does not generate additional VASP k points. `-s 2` specifies two spinor
components. `-ii 1 -if 10` selects the occupied bundle, including its internal
Kramers degeneracies; it does not compute ten independent band invariants.
`-o BERRYCURV` chooses the native output basename.

For an MPI run, build with `make mpi` and replace the executable in that
native command by `mpiexec -n 4 /path/to/vaspberry/build/vaspberry-mpi`.
The convenience Python runner itself uses serial execution.

## 4. Compare your outputs

| Produced file | Stored reference | What to check |
|---|---|---|
| `BERRYCURV.dat` | [Native output](reference/BERRYCURV.dat) | Chern number 0; selected bands 1–10; 12 × 12 mesh |
| `fortran.log` | [Native log](reference/fortran.log) | Successful completion and `Chern Number = -0.000000` |
| `plaquettes.csv` | [First-BZ table](reference/plaquettes.csv) | 144 plaquettes selected once from the native 3 × 3 display tiling |
| `band_edges.csv` | [Band edges](reference/band_edges.csv) | WAVECAR bands 10 and 11, matching archived EIGENVAL |
| `summary.csv` | [Summary](reference/summary.csv) | Zero Chern number and positive sampled gaps |
| `figure.png` | [Reference figure](reference/figure.png) | Native curvature map and occupied/empty band edges |
| `result.json` | [Run provenance](reference/result.json) | PASS, source/input/output hashes, actual command and limitations |

![Bi Fukui reference from an actual WAVECAR calculation](reference/figure.png)

The WAVECAR gives **0.592448500 eV minimum direct gap** and
**0.510044362 eV sampled global gap**. Subtracting the individually rounded
energies in archived EIGENVAL gives 0.592449 and 0.510045 eV respectively.
The curvature of this complete Bi occupied bundle is near zero.
`BERRYCURV.dat` prints four decimal places, so **all curvature values in this
file round to zero**. The map preserves that native output precision; it
does not imply an exact zero at unlimited precision. This example teaches a
real C = 0 calculation and does not serve as a nonzero-Chern material demo.

The reference files were generated by the current VASPBERRY executable from
the checked public WAVECAR. Figure or file bytes can differ between plotting
libraries/compilers; compare physical values and declared tolerances, not PNG
checksums. The runner accepts `|C| ≤ 10⁻⁵` and maximum absolute native printed
curvature ≤ 5.1 × 10⁻⁵ Å² for this exact fixture. It also verifies unique
coverage of all 144 plaquettes against the actual WAVECAR mesh, checks Cartesian
coordinates against its lattice, and compares each curvature with the stored
map at 5.1 × 10⁻⁵ Å² tolerance. Duplicate or missing coordinates fail.
The malformed-output regression checks can be run with
`python3 -m unittest discover -s tests -p test_real_fukui_example.py`.

## 5. Apply the commands to your system

1. Produce a converged VASP calculation and an **unreduced, uniform 2D mesh**
   WAVECAR. Keep the lattice, potential and SOC settings consistent; use
   `ISYM=-1` so every required mesh point is written.
2. Replace the WAVECAR path, set `-kx/-ky` to its actual mesh, and set `-s`
   for its scalar or spinor representation.
3. Choose a band or fixed band bundle isolated from the remaining bands at
   every k point. For an insulator, select the complete occupied bundle.
   Check the band edges and gap; a metallic occupied-band count that changes
   across k cannot be represented by one insulating occupied bundle.
4. Run the native command above and compare Chern number, overlaps and gap
   stability on denser meshes. A K–Γ–K′ line cannot give a full-BZ Chern number.

The checked-input runner intentionally reproduces this exact Bi fixture;
for a different material use the native command with your own parameters.
WAVECAR overlaps use pseudo-wavefunctions without PAW augmentation. This
reference does not establish mesh, basis or PAW convergence for another
system. See the [main usage guide](../../../README.md#usage) and
[Bi material provenance](../../Bi_Z2/README.md).
