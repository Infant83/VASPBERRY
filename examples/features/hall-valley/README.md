# Charge Hall response: VASP Bi WAVECAR → occupied-subspace transport

This tutorial recalculates a **Bi bilayer's insulating charge-Hall response
from its actual VASP SOC WAVECAR**. It demonstrates a complete 2D k mesh,
occupied-state selection, validity checks, output CSV, and a reference figure.
The expected total charge response is zero: Bi preserves time reversal.
Its nontrivial Z₂ invariant does not imply a nonzero charge-Hall conductivity.

Bi's Kramers pairs have numerical splittings of order `10⁻⁷ eV`. Individual-band
point-Kubo denominators are therefore inappropriate for this input. This
tutorial uses the production **WAVECAR Fukui occupied-subspace transport
command**, which retains the occupied pairs together. The separate
[MoS₂ tutorial](../kubo-curvature/) shows the native Kubo line calculation.

## 1. Obtain the actual VASP input

From the repository root:

```sh
python3 examples/fetch_inputs.py bi --output-dir results/inputs/bi
```

This retrieves the public Git LFS WAVECAR payload and checks it. Alternatively,
use an already hydrated `examples/Bi_Z2/WAVECAR` from a Git LFS checkout.
The 134-byte text pointer alone cannot be used for a calculation.

| Input | Reference |
|---|---|
| WAVECAR | 200,421,600 bytes; SHA256 `a8d81854f2efc561938e478dde1be29a17ccc95d5d37325c122b9d9e82fa0838` |
| Sampling | Full Γ-centered 12 × 12 × 1 mesh: 144 k points |
| States | 18 SOC spinor bands; occupied subspace 1:10; band 11 is the unoccupied sentinel |
| Structure and VASP records | [Bi material guide](../../Bi_Z2/README.md), [archived calculation](../../Bi_Z2/archive-2016-run/README.md) |
| VASP regeneration templates | [SCF and NSCF inputs](../../Bi_Z2/inputs/README.md); these are documented templates, not proof of the archived run's exact input provenance |

The WAVECAR is sufficient for this VASPBERRY calculation. The μ values below
use its unchanged energy zero; they are **not inferred Fermi-energy offsets**.

## 2. Run the production command

```sh
python3 tools/wavecar_fukui.py results/inputs/bi/WAVECAR \
  --nx 12 --ny 12 --spinor-components 2 \
  --energy-band 11 --map occupied=1:10 \
  --transport-full-t0 10 \
  --mu-min -1.3 --mu-max -0.9 --mu-num 41 \
  --valley-k 0.6666666666666666,0.3333333333333333,0 \
  --valley-kp 0.3333333333333333,0.6666666666666666,0 \
  --output-dir results/bi-hall-direct
```

The sampled occupied maximum is **−1.348963526 eV**, and band 11 starts at
**−0.838919163 eV**, giving a sampled global gap of **0.510044362 eV**.
All 41 selected μ values lie inside this gap, so ten spinor bands are occupied
everywhere. No extra factor of two is applied for spin.

`--transport-full-t0 10` computes cumulative subspaces through band 10 and
checks band 11 remains empty. Internal degeneracies of the complete occupied
subspace are allowed; its boundary must remain gapped. Default link, gap,
plane-wave coverage, and phase checks stay enabled. The scan also checks
occupation changes between the selected μ points.

The current transport CLI requires two centers for its regional output.
Here they are the hexagonal reciprocal-cell K/K′ centers and define only a
**geometric partition**. The physical observable shown is the **total charge
Hall response**. Do not interpret the tiny regional residuals as a valley Hall
effect in this Bi example.

## 3. Reproduce the figure and compare the results

```sh
python3 examples/features/hall-valley/run.py \
  --wavecar results/inputs/bi/WAVECAR \
  --output-dir results/bi-hall
```

The runner executes the production command above and writes a new figure.
It needs NumPy and Matplotlib; this Python WAVECAR route does not require a
Fortran binary. The reference takes roughly 10–20 seconds on the development
machine; actual time depends on the machine. Output directories are never
overwritten.

| Reference file | Contents |
|---|---|
| [figure.png](reference/figure.png) | Zero charge-Hall plateau and occupied-subspace flux residual map |
| [summary.csv](reference/summary.csv) | μ, quality flags, occupied count, total Chern integral and sheet conductivity |
| [transport_full_t0.csv](reference/transport_full_t0.csv) | Unmodified production transport output, including geometric partitions |
| [transport_full_t0_diagnostics.json](reference/transport_full_t0_diagnostics.json) | Active-subspace guards, continuous-range checks, boundary gaps, and mesh diagnostics |
| [fukui_occupied.csv](reference/fukui_occupied.csv) | Every occupied-subspace plaquette, its coordinates, flux, gap and link checks |
| [result.json](reference/result.json), [provenance.json](reference/provenance.json) | Reference checks, input/source checksums, and actual commands |

![Bi insulating charge-Hall reference](reference/figure.png)

Expected checks are:

- **41/41 μ points: PASS**, with ten occupied bands throughout the scan.
- Total `|σxy| < 10⁻¹⁰ e²/h` (typically of order `10⁻¹⁶ e²/h`). No rounding is
  applied to the CSV value.
- Minimum occupied-subspace boundary gap: approximately **0.5924485 eV**.
- Minimum common plane-wave coverage: approximately **0.9837398**.
- Minimum occupied-subspace link singular value: approximately **0.8156088**.

The output is the **2D sheet response in e²/h**, not a 3D conductivity per unit
length. The flux map shows small numerical residuals, not a finite predicted
anomalous Hall signal. Coarse-mesh reproduction is not a convergence study.

## 4. Apply the workflow to your material

Use a full, uniform two-dimensional WAVECAR mesh and declare its actual
dimensions. A line-mode or symmetry-reduced irreducible mesh is insufficient.
Choose `MAX_BAND` and the μ interval from the matching band energies; retain
at least one higher unoccupied sentinel band. An insulating interval is the
simplest first check. Preserve SOC pairs in a complete occupied subspace.

For a metallic interval the occupancy varies across the mesh, and this
plaquette method must pass its active-subspace guards and a separate mesh
convergence study. If it rejects a range, do not use `--allow-invalid-transport`
to turn a diagnostic into a physical reference. For point-Kubo Hall integration,
use appropriate velocity matrix data and the [general Kubo transport guide](../../../docs/KUBO_TRANSPORT.md),
including its completeness and degeneracy requirements.

The historical folder name `hall-valley` is retained for stable links. This
specific material reference demonstrates total charge Hall and a geometric
partition; it does not demonstrate a nonzero physical valley or spin Hall
effect. The QWZ metal/region oracle is kept under
[developer model validation](../../../validation/models/hall-valley/).

For maintainers, `python3 examples/features/kubo-curvature/export_reference.py
--run-dir results/bi-hall --output-dir results/bi-hall-reference` packages a
successful replay in a new directory. Distributed-file hashes and the original
full-run output hashes have separate, explicitly documented path scopes.
