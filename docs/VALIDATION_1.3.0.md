# Local validation of the 1.3.0 development candidate

Initial checks on 2026-09-19; scientific examples updated on 2026-09-20. This is a local source validation, not a hosted release
or a completed material-convergence study.

## Reproducible public checks

```bash
python3 -m pip install -r requirements-transport.txt
python3 -m unittest discover -s tests
make check-gnu
```

The initial Kubo/Hall validation run passed **224 tests**, including compiled
Fortran physics/parser checks. Compiler-dependent tests explicitly skip if
their compiler is unavailable; a skip is not compiler validation. GNU serial,
GNU/OpenMPI build/help/runtime and Intel Classic serial build checks passed
locally on that source. Hosted CI and Intel MPI remain separate checks.

Coverage includes existing Fukui/Z2/transport tests; independent two-band
curvature and QWZ sign/phase oracles; Hermiticity and missing vertices; full-mesh
geometry and oriented planes; occupations, partial bands and spin multiplicity;
schema types and hashes; one-time normalization migration; imported energy-gap
checks; regional/band sum rules; and small doping responses on large baselines.

The [method guide](KUBO_TRANSPORT.md) includes runnable public matrix → curvature
→ Hall commands. Those commands and the standalone
[CSV plotting example](../examples/kubo/README.md) passed local smoke checks.

## Actual VASP tutorials

The user-facing examples now start from the supplied **real MoS₂ and Bi VASP
WAVECAR files**. All six were recalculated with the current production
Fortran or Python WAVECAR path on 2026-09-19, with original outputs, figures,
checksums and numerical comparisons retained in each `reference/` directory.
The combined local run passed in approximately 76 seconds; timings depend on
hardware. The public Bi download was independently retrieved and verified
against its pinned 200,421,600-byte SHA-256 payload.

```bash
make serial
python3 examples/fetch_inputs.py bi --output-dir results/inputs/bi
python3 examples/run_examples.py \
  fukui-chern z2 hall-valley kubo-curvature circular-dichroism wavefunction \
  --bi-wavecar results/inputs/bi/WAVECAR --output-dir results/feature-examples
python3 -m unittest discover -s tests -v
```

| Actual input and calculation | Checked result |
|---|---|
| Bi 12×12 occupied bands 1:10, native Fukui | C = 0; sampled direct/global gaps 0.592448500 / 0.510044362 eV |
| Bi 12×12, native Z₂ | Z₂ = 1; two half-zone integer sums −3 and +3 |
| Bi occupied-subspace Fukui transport, T=0 | All 41 chemical potentials within the gap pass; max charge Hall residual 1.17e−16 e²/h |
| MoS₂ actual 48-point path, native Kubo | 96 band/point rows; 88 meet the stated 1e−5 eV isolation threshold; K/K′ opposite signs |
| MoS₂ actual path, native `-cd 2` | 9,648 k/energy rows; independent momentum/spectral sum agrees within 4.991e−5 a.u., within text-output rounding |
| MoS₂ Γ band 18, native `-wf` | Both spinor components; independent Fourier amplitude error 3.42e−7 Å^(−3/2), integrated pseudo-density 0.8258141586 versus coefficient norm 0.8258141556 |

The full local Python suite passed **238 tests** with no skips, including an
actual MoS₂ native Kubo batch calculation, immutable reference hashes and
archive contents, input-pointer rejection, malformed Fukui mesh/coordinate rejection, retained failures, and compiled
Fortran regressions. See the individual tutorials for comparison tolerances.
The CI `feature-examples` job downloads the actual Bi input, recalculates all
six examples and retains outputs/logs. Local success alone does not establish
hosted CI success.

An actual two-rank OpenMPI MoS₂ Kubo run also passed. Its 96-row high-precision
CSV is byte-identical to the serial reference, including the 88 valid / 8
unresolved-state mask. This check covers the native Kubo path; Python tutorial
wrappers and Python transport do not gain MPI support from it.

The [feature index](../examples/README.md) links the exact inputs, production
commands, reference outputs and [own-material guide](../examples/APPLY_TO_YOUR_SYSTEM.md).
Existing `examples/1H-MoS2/` and `examples/Bi_Z2/` inputs and historical outputs
remain unchanged. Bi is a zero charge-Hall sanity check, not a nonzero-Chern
or valley-Hall material demonstration. The MoS₂ path cannot supply a BZ
integral. Native canonical momentum does not include all PAW/nonlocal/SOC
velocity terms.

## Cartesian figures and full-zone MoS₂ example

The 2026-09-20 update adds a seventh real-material tutorial: occupied-subspace
Fukui curvature on a full 12×12 MoS₂ mesh. VASP 5.4.4 regenerated the 26-band
SOC WAVECAR from the public charge density and matching licensed potentials.
The fixed-charge run converged in eight electronic iterations (about 144 s,
peak resident memory 602 MB on the measured machine). Its minimum direct and
global gaps are approximately 1.674 eV.

Native VASPBERRY gives curvature extrema −12.3125 and +12.3124 Å²; integration
of the rounded text map gives C ≈ −5.0e−7. An independent Python overlap
calculation gives C ≈ 5.65e−16, with minimum link singular value 0.714 and
minimum plane-wave coverage 0.978. The historical MoS₂ map is retained as a
separate reference and is not claimed to be identical to the regenerated NSCF
calculation.

The final local suite passed **252 tests** with no skips. The first-BZ plotter
preserves native plaquette values and Cartesian lattice geometry. Geometry tests cover hexagonal and square cells, skew and rotated
bases, reversed orientation, area preservation, malformed meshes and lattice
mismatches. Existing Bi native outputs and the six earlier numerical reference
tables remain unchanged; the figure updates affect presentation only.

The seventh tutorial requires a locally generated VASP mesh. Public CI runs
the six tutorials with downloadable wavefunctions and tests the seventh
workflow's discovery/input checks. The [scientific report](TECHNICAL_REPORT.md)
and [MoS₂ tutorial](../examples/features/fukui-berry-curvature/) show the
resulting figures and the complete VASP preparation commands.

## Separate developer checks and historical wavefunction fix

The earlier QWZ, analytic optical, synthetic wavefunction and stored-field
checks are preserved under [validation/models](../validation/models/), outside
the user tutorials. All six still pass after the move (about 8.4 seconds
locally). They check signs, formulas, sum rules and regressions; they are not
actual VASP material demonstrations.

The synthetic wavefunction fixture previously exposed a phase-array extent
mismatch (`npmax` versus active `ncnt`) and an uninitialized POSCAR-header loop
flag. Both complete Fortran sources pass a regression compiled with
`-fcheck=all -finit-integer=1`, checking real/imaginary grids and headers.
The actual MoS₂ tutorial independently checks the retained amplitude convention.
No production numerical kernel changed during the real-input tutorial update.

## Supplemental historical comparisons

Retained local 14×14 SOC data were replayed without changing original inputs
or results. The private material files are not distributed with this source.

| Comparison | Result |
|---|---|
| Historical canonical **raw/2**, 40 bands, 1,002 chemical potentials, two temperatures | Total/per-band Hall reproduced; maximum absolute difference below `5.3e-14 e^2/h` |
| Historical PAW optical curvature, 196×40 states | Curvature arrays exactly equal |
| Historical PAW Hall spectrum, 2,006 rows | Maximum absolute difference below `1.6e-14 e^2/h` |
| Material Fukui bands 31/32 | Same +1/−1 and identical sampled output contents |

Small Hall differences are consistent with summation order. Reproducibility
targets the previously corrected canonical result; the old doubled curvature
is retained only as historical input with an explicit migration record.

## Scope

These checks validate formulas, data contracts, migration and the checked
execution paths. They do not establish k-mesh/intermediate-band convergence
for an arbitrary material or the physical completeness of an exporter.
Experimental matrix inputs retain their labels, and the bare-momentum path
retains its missing-velocity-term limitations.

The new Hall path computes intrinsic charge sheet response on uniform full
2D meshes. Invalid individual-band curvature is rejected. Spin/layer currents,
adaptive integration and general 3D bulk response need separate implementations.
Existing Fukui subspace workflows remain available for their documented scope.
