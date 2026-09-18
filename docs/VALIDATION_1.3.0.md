# Local validation of the 1.3.0 development candidate

Checked on 2026-09-19. This is a local source validation, not a hosted release
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

## Feature examples and wavefunction follow-up

The subsequent feature-example update, also checked on 2026-09-19, passed
**230 tests** with no skips in the local Python discovery run. A fresh serial
build ran all six defaults in approximately 8.8 seconds on the validation
machine (timing is illustrative, not a performance guarantee):

```bash
make serial
python3 examples/run_examples.py --all --output-dir results/feature-examples
python3 -m unittest discover -s tests -v
```

| Example | Checked result |
|---|---|
| Fukui/Chern | QWZ phases +1, −1, 0; empty/filled lower-band Hall signs |
| Matrix Kubo | Pointwise oracle error below 7.3e-16 Å²; Chern integral error below 9.4e-10 |
| Hall/regions | Direct occupation/curvature oracle and region/band sums agree below 8e-15 |
| Z₂ | Stored Bi schema-2 field validates with matching parity 1; no new Bi WAVECAR run |
| Optical selectivity | Current Fortran angular dependence agrees within 4.7e-7 |
| Gamma wavefunction | Current Fortran complex amplitudes agree within 5e-7 |

The new wavefunction fixture exposed a pre-existing phase-array extent mismatch
(`npmax` versus active `ncnt`) and an uninitialized POSCAR-header loop flag.
Both complete Fortran sources now pass a regression compiled with
`-fcheck=all -finit-integer=1`, checking the real/imaginary grids and header.
This retains the existing amplitude output convention, documented by the
[wavefunction example](../examples/features/wavefunction/).

The [feature index](../examples/README.md) links the exact inputs, compact
reference results, figures and provenance. Existing MoS₂ and Bi material
files remain unchanged. The CI `feature-examples` job repeats the six defaults
and retains results/logs; local success does not imply hosted CI success.

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
