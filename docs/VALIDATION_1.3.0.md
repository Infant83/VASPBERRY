# Local validation of the 1.3.0 development candidate

Checked on 2026-09-19. This is a local source validation, not a hosted release
or a completed material-convergence study.

## Reproducible public checks

```bash
python3 -m pip install -r requirements-transport.txt
python3 -m unittest discover -s tests
make check-gnu
```

The final local Python discovery run passed **224 tests**, including compiled
Fortran physics/parser checks. Compiler-dependent tests explicitly skip if
their compiler is unavailable; a skip is not compiler validation. GNU serial,
GNU/OpenMPI build/help/runtime and Intel Classic serial build checks passed
locally on the final source. Hosted CI and Intel MPI remain separate checks.

Coverage includes existing Fukui/Z2/transport tests; independent two-band
curvature and QWZ sign/phase oracles; Hermiticity and missing vertices; full-mesh
geometry and oriented planes; occupations, partial bands and spin multiplicity;
schema types and hashes; one-time normalization migration; imported energy-gap
checks; regional/band sum rules; and small doping responses on large baselines.

The [method guide](KUBO_TRANSPORT.md) includes runnable public matrix → curvature
→ Hall commands. Those commands and the standalone
[CSV plotting example](../examples/kubo/README.md) passed local smoke checks.

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
