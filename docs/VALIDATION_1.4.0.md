# Validation of VASPBERRY 1.4.0

This record separates release software checks from convergence of a material
prediction. The immutable release's source SHA and hosted validation URLs
are attached to its [publication record](https://github.com/Infant83/VASPBERRY/releases/tag/v1.4.0).

## Release checks

The release workflow requires the full CI matrix and actual Bi serial/MPI
validation on its own exact source commit. CI covers Python 3.10/3.12,
GNU 11/13 serial/MPI, actual-source Z₂ guards, real Bi/MoS₂ feature commands,
and Intel ifx/ifort compilation and help. Intel numerical and MPI behavior
are not claimed by those compile-only checks.

Final local preparation passed all 510 unit/regression tests with no skips,
GNU serial/MPI compilation and two-rank smoke checks. The projection subset
contains 25 tests. The public analytic workflow ran its five CLI stages,
produced 5,396 Hall rows and six plot files, and agreed with its independent
analytic values to 2.22×10⁻¹⁶ e²/h. An independent oblique-cell, three-band
calculation also checked arbitrary-axis projections and virtual-state sums.

Documentation checks covered 85 Markdown files, 717 local links, 26 anchors
and 135 explicit Python-command references. The technical report was rebuilt
as a 20-page PDF and visually inspected, preserving its eleven report figures.
The immutable source receives separate hosted checks before publication.

## Independently checked native changes

The pair-export occupation exemption applies only to explicit pair export:
its numerators do not depend on occupation. Otherwise identical input files
with changed occupation records give byte-identical pair outputs, including
serial/MPI comparison. Fixed-subspace task checks remain active.

Canonical velocity was tested against analytic plane-wave states for scalar,
SOC and explicit collinear spin channels, including extrema coordinates and
array-bounds checks. Actual MoS₂ comparison against an independent SI
calculation confirms the scale. Older velocity files require regeneration.
The nine Kubo calculation/export subroutines were independently compared with
the pre-fix baseline and retain their numerical kernels.

## Reproduction evidence retained from preparation

- Report-to-example commands regenerated all eleven existing report figures.
  Optional model panels were redrawn from their saved model results; they
  are not described as newly produced Wannier calculations.
- MoS₂ Fukui, bundle Kubo mesh/path, Z₂=0 and fresh Bi Z₂=1 were compared
  with their stored numerical references. MnBi₂Te₄ native four-rank occupied
  Chern calculation gave C=−1 with all 324 printed numerical rows preserved.
- A new 24×24/source60 native pair export followed by pair-window40 integration
  reproduced the MoS₂ reference NPZ exactly when using its full-precision μ.
- Eight saved Hall cases, 4,880 rows, gave equal CSV/DAT/NPZ numerical values.
  Reusing fixed pair data with different μ/T/region choices and independent
  readers produced additional curves and maps. Twelve matched PAW/canonical
  comparison tables were reproduced byte for byte.

These runs retain their individual source/binary hashes and timestamps.
Several actual-material checks precede the final version-label and projection
changes; the exact release commit's hosted checks provide separate evidence.
No new VASP production calculation was needed for release preparation.

## Projection interpretation

The [public analytic example](../examples/features/procar-character/) is a
small generated fixture with known projections and pair numerators. It is
explicitly separate from the report's VASP material benchmarks. The workflow
preserves raw projection coverage and uses the selected bands' normalized
curvature. It rejects ambiguous selected-state degeneracies and mismatched
source or spin-frame inputs. The supported results are state character and
projected charge-Hall attribution, not conventional spin-current conductivity.

Material convergence remains open where noted in the report: the Bi SHC
12×12→18×18 difference is 6.7544%, and the coarse MnBi₂Te₄ optical integral
is not a converged Hall prediction. Historical numerical references remain
unchanged and retain their original producer versions.
