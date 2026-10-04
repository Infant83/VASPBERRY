# VASPBERRY feature poster

[View the one-page poster (PDF)](VASPBERRY_FEATURE_POSTER.pdf).

Revised 5 October 2026; example and capability overview for VASPBERRY 1.6.6.
The six panels cover Berry-curvature maps and paths, charge and projected-spin
Chern numbers, the 2D Z₂ n-field, regional Hall response, circular dichroism,
and real-space wavefunctions. Original method references are linked in the PDF.

The poster illustrates documented example calculations. The Bi maps show
discrete FHS plaquette flux; the selected-pair spin Chern number is distinct
from the full occupied-space invariant. The projected-spin results use the
WAVECAR pseudo-Gram metric without PAW augmentation.

The MnBi₂Te₄ bands provide context from a fixed VASP-derived Wannier model
whose localization did not meet its convergence threshold. The occupied
Chern number comes from a separate native VASP calculation.

The Hall panel shows a valley-resolved change using the WAVECAR approximation.
The Kubo/optical figures labelled bare momentum are archived examples of that
approximation; same-run WAVEDER is the standard charge-Kubo input in version
1.6.6. Display smoothing, where labelled, does not add calculated k points.

Citation and GitHub-star counts are the 4 October 2026 snapshot. Citations are
summed across the code-related 2016 PRB and 2022 PRL papers, not a count of code
users.

See the [technical report](TECHNICAL_REPORT.md),
[feature examples](../examples/README.md), and
[standard Kubo protocol](WAVEDER_KUBO_PROTOCOL.md) for inputs, validation and
limitations. The poster is a visual overview, not an additional convergence
study.

## Sharing images

The [full-resolution PNG](VASPBERRY_FEATURE_POSTER.png) is displayed in the
repository README and links to the PDF.

The [social-preview PNG](VASPBERRY_SOCIAL_PREVIEW.png) is 1280×640 and under
1 MB. It preserves the entire poster with side margins. Repository owners
can upload it under Settings → Social preview. This setting is separate
from the README image; committing this file alone does not change link cards.
