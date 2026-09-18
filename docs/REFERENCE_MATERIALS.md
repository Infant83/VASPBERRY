# Supplementary references for hands-on guides and technical reports

This index connects the actual VASP tutorials with the analytic, synthetic and
historical reference checks. Both sets are retained with their numerical
outputs and figures for future hands-on documents and release technical
reports.

Start a user exercise with the [actual VASP tutorials](../examples/README.md):
obtain the VASP files, run VASPBERRY, compare the saved output and figure, then
[adapt the settings to another material](../examples/APPLY_TO_YOUR_SYSTEM.md).
Use the supplementary checks below to explain the numerical method, sign
conventions and regression evidence behind that exercise.

## Reference map

| Topic | Actual VASP exercise and reference | Supplementary reference | Use in a guide or report |
|---|---|---|---|
| Fukui / Chern | [Bi occupied bundle](../examples/features/fukui-chern/) | [QWZ phases](../validation/models/fukui-chern/) | Explain lattice flux, known Chern phases and Hall sign; Bi supplies the actual C = 0 workflow |
| Kubo curvature | [MoS₂ band path](../examples/features/kubo-curvature/) | [Analytic matrix curvature](../validation/models/kubo-curvature/) | Explain the curvature formula, units and comparison with an independent exact model |
| Charge Hall / regions | [Bi insulating-gap response](../examples/features/hall-valley/) | [Model occupations and regions](../validation/models/hall-valley/) | Explain occupation weighting and regional/band sum rules; Bi demonstrates a zero charge-Hall plateau |
| Z₂ | [Bi WAVECAR recalculation](../examples/features/z2/) | [Stored Bi field validation](../validation/models/z2/) | Explain field diagnostics, schema compatibility and plotting separately from recalculation |
| Optical response | [MoS₂ spectra](../examples/features/circular-dichroism/) | [Synthetic angular selectivity](../validation/models/circular-dichroism/) | Explain polarization conventions and an analytic angular-response check |
| Wavefunction | [MoS₂ Γ spinor](../examples/features/wavefunction/) | [Synthetic scalar Γ state](../validation/models/wavefunction/) | Explain Fourier reconstruction, amplitude conventions and the array-bounds regression |

Each linked directory documents its inputs, commands and `reference/` files.
The model and synthetic fixtures use JSON settings to generate test data;
the actual-material exercises use VASP outputs. The stored Z₂ field comes
from a material calculation, but validating that field alone is not a new
WAVECAR calculation.

## Suggested hands-on sequence

1. Follow one real-material tutorial through native output and its plotted
   reference. Keep the input checksums and the generated run record.
2. Add the matching supplementary reference as a method note or optional
   exercise: identify the expected answer and inspect the numerical error.
3. Explain which settings the participant must choose for their own material,
   including mesh, band window, spinor representation and energy reference.
4. Keep model and material figures separately captioned. Compare the relevant
   method or convention; the different datasets need not share numerical
   curves or parameter values.

For example, the QWZ references help explain nonzero Chern phases and numerical
integration. The supplied Bi example has occupied C = 0 and Z₂ = 1, while the
MoS₂ path has insufficient sampling for a Brillouin-zone Hall integral. Model
figures therefore complement these exercises without supplying missing
material results or a physical valley identification.

An exact model-vertex check does not establish the completeness of the native
material velocity operator, and a synthetic scalar state does not establish
SOC/PAW accuracy. Use these checks for their stated numerical purpose.

## Evidence to retain in a release technical report

- **Method checks:** independent expected values, signs, units, tolerances and
  residuals from [developer checks](../validation/models/).
- **Actual-input reproduction:** input provenance/checksums, exact production
  commands, original output, comparison results and figures from the
  [VASP tutorials](../examples/README.md).
- **Implementation and compatibility:** the exact source commit, build/runtime
  environment, successful and failed checks, serial/MPI comparison where
  measured, and [historical normalization migration](MIGRATION.md).
- **Interpretation:** label the input kind, operator approximation, sampling,
  unresolved states and what has or has not been converged.

Use the [1.3.0 validation record](VALIDATION_1.3.0.md) and
[release notes](releases/v1.3.0.md) as the current evidence index. Historical
comparisons that use undistributed material inputs remain separately labeled;
they are not exercises a public checkout can reproduce. Historical
reference hashes and source locations describe the recorded run; retain them
when moving or reusing a figure. A new run should have its own output folder
and provenance. In figure captions, name the dataset, method, plotted units,
important parameters and source commit, and link the numerical data used.
