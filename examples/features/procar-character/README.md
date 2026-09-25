# Atom, orbital and spin character of charge-Hall contributions

This v1.4.0 workflow combines noncollinear `LORBIT=11` PROCAR projections with
the **same final WAVECAR** used for native Kubo pair export. It supports arbitrary
named atom groups, optional orbital selectors, a Cartesian spin-analysis axis,
chemical-potential/temperature scans and the existing region JSON format.

The result is a **selected-band, character-weighted charge-Hall attribution**.
It is useful for comparing which layer, atomic plane, orbital or spin projection
contributes to a canonical-momentum Hall trend. It does not compute an independent
layer-current, orbital-current or spin-current operator. Wannier is unnecessary.

## Run the small analytic example

Run from the repository root, using a Python environment with
[`requirements-transport.txt`](../../../requirements-transport.txt). Choose fresh
output directories: commands refuse to overwrite previous results.

```sh
python examples/features/procar-character/make_fixture.py \
  --output-dir results/procar-demo/input

python tools/procar_character.py project \
  --procar results/procar-demo/input/PROCAR \
  --wavecar results/procar-demo/input/WAVECAR \
  --outcar results/procar-demo/input/OUTCAR \
  --groups results/procar-demo/input/groups.json \
  --axis 1 0 0 --output-dir results/procar-demo/character

python tools/vaspberry_kubo.py import-pairs \
  --csv results/procar-demo/input/PAIRS.csv \
  --wavecar results/procar-demo/input/WAVECAR \
  --spinor-components 2 --spin-multiplicity 1 --mesh 2 2 \
  --energy-reference 'synthetic analytic fixture zero' \
  --output-dir results/procar-demo/pairs

python tools/procar_character.py hall \
  --character-dir results/procar-demo/character \
  --pairs-dir results/procar-demo/pairs --bands 1 \
  --mu-min -1.2 --mu-max 0.2 --mu-num 71 --mu-reference 0 \
  --temperatures 0 300 --output-dir results/procar-demo/hall

python tools/procar_character.py plot \
  --character-dir results/procar-demo/character \
  --hall-dir results/procar-demo/hall \
  --group lower --band 1 --temperature 300 \
  --output-dir results/procar-demo/figures
```

The fixture has two atoms, two orbitals, two bands and a 2×2 mesh. Its OUTCAR
rotation maps the spin-frame z axis to Cartesian +x. The lower atom has raw
charge 0.2 and +x character; the upper atom has charge 0.6 and −x character.
`PAIRS.csv` assigns `N_xy=4 eV² Å²` and the band gap is 2 eV, so the selected
band has `Omega_z=1 Å²`. These are **analytic format fixtures**, not a VASP
material result; the assigned pair numerator is not calculated from the toy
WAVECAR coefficients. No private material inputs or VASP license are needed.

At `T=0`, `mu=0`, band 1 gives:

| Quantity | Expected value in e²/h |
| --- | ---: |
| Unweighted selected-band charge Hall | −2π |
| Lower charge / lower plus | −0.4π |
| Lower minus / upper plus | 0 |
| Upper charge / upper minus | −1.2π |
| All-ion projected charge | −1.6π |
| Charge projection residual | −0.4π |

[`run_example.py`](run_example.py) runs the same CLI stages, verifies those
numbers, checks joint-spin and residual closure, then saves a verification JSON.
For CI or a one-command check:

```sh
python examples/features/procar-character/run_example.py \
  --output-dir results/procar-demo-checked
```

## Use your own VASP calculation

Keep `PROCAR`, `WAVECAR` and `OUTCAR` from one completed static, noncollinear
`LORBIT=11` calculation. The current reader supports the four-block
charge/m₁/m₂/m₃ format used by VASP 5.4.4; unsupported phase/collinear formats
fail explicitly. The `project` command checks every k point, band energy,
occupation, OUTCAR lattice and final-state table. It reads the actual printed
SAXIS-to-Cartesian matrix; `--axis 0 0 1` always means **Cartesian** +z.

For a layer/orbital definition, create `groups.json`, for example:

```json
{
  "groups": [
    {"name": "lower", "ions": [1, 3]},
    {"name": "upper", "ions": [2, 4]},
    {"name": "lower_d", "ions": [1], "orbitals": ["dxy", "dyz", "dz2", "dxz", "x2-y2"]}
  ]
}
```

Atom IDs are 1-based POSCAR/PROCAR indices. Determine these lists from the
actual structure; a pair of magnetic atomic planes inside one monolayer is
not automatically a bilayer. Orbital names must match the PROCAR header exactly.
Groups may overlap; do not sum overlapping groups as a partition.

First build the native executable with `make serial`. For your **full uniform
2D mesh**, export and import the occupation-independent pair cache:

```sh
mkdir -p results/my-sample-native
build/vaspberry --task kubo-pairs --wavecar /absolute/path/to/WAVECAR \
  --spinor 2 --pairs-csv results/my-sample-native/PAIRS.csv

python tools/vaspberry_kubo.py import-pairs \
  --csv results/my-sample-native/PAIRS.csv \
  --wavecar /absolute/path/to/WAVECAR \
  --spinor-components 2 --spin-multiplicity 1 --mesh 24 24 \
  --energy-reference 'unchanged VASP eigenvalue zero' \
  --output-dir results/my-sample-pairs

python tools/procar_character.py project \
  --procar /absolute/path/to/PROCAR --wavecar /absolute/path/to/WAVECAR \
  --outcar /absolute/path/to/OUTCAR --groups groups.json \
  --axis 0 0 1 --output-dir results/my-sample-character

python tools/procar_character.py hall \
  --character-dir results/my-sample-character --pairs-dir results/my-sample-pairs \
  --bands 19 20 --mu-min -0.5 --mu-max 0.5 --mu-num 101 \
  --mu-reference 0 --temperatures 0 100 300 \
  --output-dir results/my-sample-character-hall
```

Replace the illustrative 24×24 mesh, bands, energy range and reference with
your system's values. Add `--regions regions.json` for periodic valley circles
or selected k IDs, using the same schema as
[`kubo-hall`](../kubo-hall/README.md). `--plane-axes` on pair import selects a
different reciprocal plane when needed. The native GNU build needs a compatible
byte-RECL WAVECAR; see [native commands](../../../docs/NATIVE_COMMANDS.md).

Within `procar_character.py`, only `project` reads the VASP files.
Once both caches exist, rerun `hall`
with a new output directory to change selected bands, mu, temperature or
regions; rerun `project` to change atom/orbital groups or the spin-analysis axis.
Change plots independently using the saved CSV/NPZ tables.

## Outputs and interpretation

- `character.csv`, `character.npz`, `character.json`: state/group charge `q`,
  projected Cartesian Pauli weights, `m_axis`, and the **joint** projectors
  `plus=(q+m_axis)/2`, `minus=(q-m_axis)/2`. They are not group charge multiplied
  by the whole-state polarization. `S_axis/hbar=m_axis/2` if that quantity is desired.
- `projection_diagnostics.csv`: all-ion charge, raw `1−q` residual, all-ion
  Cartesian Pauli weights and the independently rounded printed state charge.
  No residual spin is inferred. PAW weights are neither clipped nor normalized
  to unity; small negative values from projection/rounding remain visible.
- `character_hall.csv`, `character_hall.npz`, `character_hall.json`: raw and
  reference-subtracted attribution for each group, spin component, region, mu
  and temperature. Reserved diagnostic groups `$unweighted`, `$all_projected`
  and `$unprojected_residual` expose the unweighted selected-band response,
  all-ion projected response and their difference.
- `selected_curvature.npz`: physical Cartesian curvature in Å², selected band
  IDs, full-source virtual-band gap diagnostics, coordinates and energies.
- `character.*` and `character_hall.*` in the plot directory: PNG/PDF/SVG.

The integration uses `−A_BZ/(2π) Σ_kn w_k f_n Omega_n·normal c_ng` in e²/h.
The current native pair normalization already contains `−2 Im`; no legacy
factor of one-half is applied. All stored source bands enter the virtual-state
sum. `--bands` specifies the states whose character-weighted contributions are
reported, not a declaration that omitted occupied bands have zero response.
Use ordinary `pair-hall` alongside this analysis for the total charge response.

Every selected state must be separated from **every other stored band** by
more than `--degeneracy-threshold-eV` (default `1e-5`). Equal occupations do not
remove gauge dependence from character-weighted degenerate states, so this
workflow rejects them and offers no coalescence shortcut. Inspect/refine the
selected window or use the ordinary pair/bundle total instead. The saved
WAVECAR hash must match between both caches; neither header agreement nor hashes
alone can prove that user-supplied PROCAR and OUTCAR came from that run.

Canonical momentum is a useful working approximation for comparative trends;
mesh, virtual-band cutoff and missing PAW/nonlocal/SOC/U velocity corrections
still require convergence/accuracy checks for quantitative claims. Mu/T scans
change occupations of fixed bands, not self-consistent doping or magnetism.
See [Kubo transport](../../../docs/KUBO_TRANSPORT.md),
[spin-response scope](../../../docs/SPIN_HALL.md), and the
[technical report](../../../docs/TECHNICAL_REPORT.md).
