# Selected WAVEDER bands: executable software example

This example generates **synthetic format and contraction fixtures**. Its
WAVECAR, WAVEDER, INCAR and OUTCAR are not outputs of a VASP calculation and
do not model a material. The purpose is to check the public commands, required
matrix coverage, μ/T weighting and figure labels. It establishes no mesh or
NBANDS convergence. For actual source preparation, use the
[standard WAVEDER protocol](../../../docs/WAVEDER_KUBO_PROTOCOL.md).

Run from the repository root with Python, NumPy and Matplotlib available.
Use a new output directory each time:

```sh
python3 examples/features/waveder-selected/run.py \
  --output-dir results/selected-demo --binary build/vaspberry
```

The native binary is optional; omit `--binary` to run the Python INI, plot,
direct-sum oracle and weighted-pair rejection checks alone. Build the current
binary first to include the native checks. Each command, exit status and
stdout/stderr log is retained. `result.json` records the outcome; expected
rejections remain visibly separate from successful calculations.

The generator builds a 2×2 grid with four source bands and a declared
Hermitian optical connection. Bands 1–2 form an exactly degenerate group;
their erased internal connection is zero. The `full/` source has all four
ket columns. The `rectangular/` source keeps just columns 1–2 without altering
the four source energies or NBANDS. `expected.json` retains the direct
all-intermediate-band curvature sums, input hashes and fixture definition.

| Check | Expected behavior |
|---|---|
| Full source, selected 3:4, Python μ/T scan | 27 μ/T rows agree with the independent direct-band sum; selected-contribution metadata and figure label. |
| Full or rectangular source, native selected 3:4 trace | Same geometric trace. Reverse orientations supply every external pair; the missing 3↔4 internal pair cancels at equal unit weights. |
| Full source, native selected 3:4 with `--per-band 1` | Separate bands agree with direct sums over all four intermediate states. |
| Rectangular source, native selected 3 | Rejected: the required 3↔4 pair is missing. |
| Full source, native selected 1 | Rejected: it splits the producer-degenerate 1–2 group. |
| Rectangular source, selected 3:4 at 300 K | Rejected: unequal Fermi weights require the missing 3↔4 pair. |

For the INI stages individually, generate a new fixture and run the copied
[configuration](selected.ini):

```sh
python3 examples/features/waveder-selected/make_fixture.py \
  --output-dir results/selected-input
python3 tools/vaspberry_post.py check results/selected-input/analysis.ini
python3 tools/vaspberry_post.py run results/selected-input/analysis.ini
python3 tools/vaspberry_post.py plot results/selected-input/hall-selected
```

The INI uses `[hall] bands = 3:4`, not `occupied`. It computes T=0, 100 and
300 K on the full synthetic source. The generated `hall-selected/hall/`
contains conductivity tables and metadata; `hall-selected/figures/charge-hall/`
contains PNG, PDF, SVG and `plot.json`. The figure title and saved metadata
identify **Selected-band contribution: 3,4**. Its `total` region means the
full sampled reciprocal region of this contribution, not total material AHC.

Native selected geometry can be inspected without an occupation scan:

```sh
build/vaspberry --task kubo --input-dir results/selected-input/rectangular \
  --bands 3:4 --mesh 2,2 --curvature-csv results/selected-input/TRACE.csv
build/vaspberry --task kubo --input-dir results/selected-input/full \
  --bands 3:4 --per-band 1 --curvature-csv results/selected-input/BANDS.csv
```

The trace uses unit geometric weights and does not claim a physical total
Hall response. Native `--mesh` uses one comma-separated value; Python uses
two values. An explicit native selection, even `--bands 1:2`, remains the
selected-geometric mode. Omitting native `--bands` instead infers the complete
insulating occupied space and enables that separate physical Hall contract.

For actual completed optical inputs, replace the directory and choose
source-appropriate bands, grid and energy values. Selectors accept a single
band, an inclusive range or a list such as `31,33:34`; they never truncate
the intermediate sum over all source NBANDS. Python `--bands` and
`--occupied` are mutually exclusive. A selected contribution does not become
the total AHC at finite temperature, and missing required matrix elements
are never padded. All four inputs must belong to the same completed standard
VASP 5.4.4 optical run. `--input-dir` defaults to the invocation directory;
per-file overrides change only their file and retain cwd-relative paths.
