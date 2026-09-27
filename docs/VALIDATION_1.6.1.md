# VASPBERRY 1.6.1 validation scope

This compatible feature release adds native projected-spin `spin-chern` and
`spin-kubo`, including `--sum-bands` for a response cutoff on unchanged source
states. Existing ordinary Kubo, Fukui, Z₂ and postprocessing interfaces remain
available. [Release notes](releases/v1.6.1.md) explain the commands and outputs.

## Publication and installation

Publication requires thirteen successful jobs on the exact release commit:

- Two Python versions run the full regression suite, including compiled spin
  geometry, error conditions, MPI and source metadata checks.
- Feature examples, the existing Fortran Z₂ regression, GNU 11/13 portability
  and Intel ifx/ifort serial portability run in the main CI workflow.
- The separate actual Bi workflow checks serial/two-rank MPI Z₂ against the
  retained material reference.
- Both Intel MPI compilers execute native numerical checks.
- Clean source-archive installations use GNU/OpenMPI on macOS ARM64 and
  GNU/MPICH on Linux, including the new native include files.

The new spin tasks also receive numerical validation with prebuilt compiler/MPI
executables. The publication guard verifies workflow names, all required jobs,
source SHA, run attempt, previous tag identity and the current default branch.
A failed or skipped required job prevents publication. The GitHub release body
records the successful runs and exact source; it does not certify material
convergence. Historical tags and reference data retain their original identity.

## Scientific implementation

The selected pseudo-wavefunction Gram matrix, full projected-spin matrix and
spin-frame rotation define the two sectors. Energy isolation, spin gap, sector
ranks, link conditioning and finite results are checked explicitly. Splitting
an unresolved energy multiplet or a retained-sum boundary is rejected.

Independent exact-Hamiltonian tests compare all Cartesian sector curvatures
with small geometric loops, include a nonzero spin-mixing derivative, and check
axis reversal and selected-frame changes. An exact-model integral approaches
the independent Fukui integer as its mesh is refined. Real Bi cutoff checks
also agree with an independent NumPy calculation. Default/full-band results
are unchanged by introducing the optional retained sum cutoff.

These checks validate the implemented geometric kernel and declared tangent.
The ordinary WAVECAR Kubo route uses canonical momentum and a finite pseudo
basis. It is not the complete physical PAW/nonlocal/SOC velocity, nor a
conventional spin-current Hall conductivity.

## Material examples and remaining research

The [technical report](TECHNICAL_REPORT.md#310-follow-up-convergence-tests)
([PDF](TECHNICAL_REPORT.pdf)) separates these follow-ups:

| Follow-up | Verified result | Remaining physical calculation |
|---|---|---|
| [Bi](../examples/materials/bi-spin-hall/spin-chern-kubo/convergence/) | Occupied/pair spin Chern 1/2 persists on 6×6, 12×12, 18×18 meshes; a genuine full-mesh restart agrees with partition aggregation. | Pair Kubo integration near Γ is unresolved; the occupied estimate improves without a convergence certificate. High empty-state quality needs separate checks. |
| [Graphene](../examples/materials/graphene-spin-chern/kubo/convergence/) | Actual local VASP samples resolve the tiny-SOC peak. A stricter solver repeat changes native curvature by at most 0.0156%; the local geometric comparison differs by about 5.5%. | Radial point quadrature is still sample-sensitive; local polygons and fitted diagnostics are not full-BZ invariants. |
| Graphene selected pair 7–8 | The failed sector-link check is preserved, with no assigned global pair Chern number. | A sampled energy/spin gap alone does not certify resolved global connectivity. |

These limitations are part of the documented numerical scope, not hidden
successful-integer claims. Users can calculate and export local curvature,
inspect diagnostics, choose complete isolated groups, and refine sampling
without changing the output interface. No empirical rescaling or rounding is
used to force a Kubo estimate to match a Chern integer.

Saved example files retain the producer labels, source hashes and failed-case
records from the development calculations that generated them. They are not
relabeled as freshly executed release outputs. The release validation checks
the final source independently; reproductions should record their own version,
inputs, compiler, sampling and retained band range.
