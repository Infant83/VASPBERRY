# MoS₂: compare two current operators on the same states

The [ordinary WAVECAR Hall example](../README.md) is the starting workflow.
This optional comparison keeps the MoS₂ structure, density, wavefunctions,
k points, occupations and retained-band window fixed, and changes the current
operator from pseudo-wavefunction canonical momentum to the full PAW velocity
of the supported producer. VASP provides the electronic states; VASPBERRY
calculates both responses from them.

The full-velocity matrices include PAW, nonlocal and SOC terms within the
producer's declared Hamiltonian scope. They offer a more complete operator
description. This comparison does not establish a converged material value,
nor imply that the more complete operator converges faster with k sampling.
The [operator guide](../../../../docs/OPERATOR_ROUTES.md) distinguishes ordinary
WAVECAR, unmodified standard WAVEDER and optional instrumented VASP inputs.

## Recalculate from VASP

Use the [parent MoS₂ preparation](../README.md#1-prepare-and-run-vasp) for the
licensed Mo/S PAW inputs and actual public charge density. Choose a full
12×12 mesh and 60 stored bands for this comparison. Converge the electronic
states, retaining extra empty states above the pair window.

For the optional full matrices, follow the
[producer build and run instructions](../../../../tools/vasp544_spin_bridge/README.md)
in an isolated supported VASP source tree. They specify the extra export
settings, matching density and physical inputs, completed-run wrapper and
validation. This setup is additional to the ordinary WAVECAR calculation.
The comparison driver then performs both postprocessing routes on the **same
final WAVECAR**:

```bash
python3 examples/features/kubo-hall/operator-comparison/run.py \
  --run-dir results/mos2-full-velocity-vasp \
  --binary build/vaspberry-gfortran \
  --pair-band-max 40 50 \
  --output-dir results/mos2-operator-comparison
```

The driver is specific to this neutral MoS₂ example (18 occupied spinor
states); all computational kernels are general. It validates the producer,
calls `velocity-pairs`, runs the native canonical-pair exporter on the matching
WAVECAR, and calls `pair-hall` for each operator and cutoff. The same periodic
K/K′ disks of radius 0.35 Å⁻¹ and the same 0/300 K occupations are used.

Each output directory contains:

- `physical/`: validated full matrices and their operator metadata;
- `canonical-pairs/` and `paw-pairs/`: reusable, explicitly distinguished pair caches;
- `canonical-capM/` and `paw-capM/`: equivalent CSV, DAT and NPZ conductivity tables;
- `result.json` and `logs/`: input association, actual commands and execution status.

Failed stages remain recorded and do not become reference results. Neither
matrix Hermiticity nor charge cancellation is imposed by averaging.

## Reuse the numerical outputs

The supplied [reference results](reference/) contain both actual pair caches,
all four conductivity scans and the figure. Recalculate their numerical
tables without a VASP installation:

```bash
python3 examples/features/kubo-hall/operator-comparison/reproduce.py \
  --output-dir results/mos2-operator-reproduced
```

This checks the regenerated CSV, DAT and NPZ against the reference values.
The final VASP input files and compressed OUTCAR, EIGENVAL and OSZICAR are
in `reference/vasp/`. The physical matrix validation records are in
`reference/physical-validation/`; regenerate the full matrices with the
producer recipe when repeating that earlier stage. Licensed POTCAR and the
large WAVECAR and raw matrix streams are not bundled. The ordinary parent
SCF preparation and the final input files document their regeneration.

Once the two pair caches are available, changing the chemical potential,
temperature, retained-band cutoff or valley partition requires only the
ordinary `pair-hall` command. No VASP executable or source modification is
needed for that integration. Use the same options for both caches and keep
the operator identity in every resulting table's metadata.

Plot a completed comparison with:

```bash
python3 examples/features/kubo-hall/operator-comparison/plot.py \
  --input-dir results/mos2-operator-comparison \
  --output-dir results/mos2-operator-comparison-figures
```

The figure shows regional Δσ and the unhalved K−K′ difference. For this
nonmagnetic system the total charge Hall response should vanish by time
reversal; nonzero regional contributions make the operator difference
visible. A region partition is not a separately conserved valley-current
operator. Separate k-mesh and source/retained-band studies remain necessary
for each operator.

## Reference observations

![Matched MoS2 Hall operators](reference/figures/mos2-operator-comparison.png)

At 300 K and μ−Ev = −0.20 eV, the K−K′ change is −0.227790 e²/h for
canonical momentum and −0.333633 e²/h for full PAW velocity, with `M=40`.
The source has 60 bands; only the lower verified windows are retained.
Increasing `M` from 40 to 50 changes the whole 300 K difference curve by
1.445% and 0.465%, respectively. Relative L2 uses the `M=50` curve as its
denominator. At `M=50`, the operator difference is 29.44% relative to PAW.
These are fixed-mesh comparisons, not a joint convergence result.

Maximum total charge Hall residuals over both cutoffs and temperatures are
1.071×10⁻⁸ and 5.118×10⁻⁹ e²/h. Independent ordered-matrix sums agree with
the production integrator within 1.54×10⁻¹⁵ e²/h. The first 40 and 50 states
satisfy time-reversal energy pairing within 2.18×10⁻⁷ and 5.10×10⁻⁷ eV;
the highest source states are less accurate and are excluded. No response
is repaired or forced to vanish. The [numerical validation record](reference/validation.json)
retains the format, source, partition and operator checks.

The accepted warm VASP export took 11.4 minutes on one CPU thread, with
2.88 GiB sampled peak RSS. This follows the ordinary density/eigenstate
preparation and includes optical exports; it is not the cost of a fresh
DFT calculation. Both VASPBERRY workflows and all four Hall scans took
44 seconds in that run. Timings depend on hardware and concurrent work.
