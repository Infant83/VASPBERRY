# Actual Bi reference outputs

These files come from the fresh two-atom Bi PBE+SOC calculation documented in
the [material guide](../README.md). The physical inputs and operator exports
are in [`../inputs/`](../inputs/).

| Directory or file | Calculation |
|---|---|
| `z2/` | Native occupied-band Fukui–Hatsugai n-field, Z₂ = 1 |
| `fukui/` | Separate full occupied-bundle Chern calculation, C = 0 |
| `hall-6x6-b48/`, `hall-12x12-b48/` | Full 48-band PAW matrix integration on each full mesh |
| `hall-6x6-b48-cap*/` | Retained-band study on the same 48-band source |
| `hall-6x6-b64-cap48/` | Enlarged VASP eigensystem, same retained 48-band space |
| `hall-12x12-b64-cap48/` | Matching accurate 48-band space for the mesh study |
| `hall-18x18-b64-cap48/` | Direct 324-point extension with the same source/retained counts |
| `hall-6x6-b80-cap*/` | Accurate retained 48, 56 and 64 states from an 80-band eigensystem |
| `hall-6x6-b64/` | Diagnostic using all 64 states; its upper empty states are **not converged** |
| `bands/`, `direct-dft/` | Native Wannier dispersion and separate actual VASP markers |
| `edge-w20/`, `edge-w40/` | Native ideal-strip energies and boundary probabilities |
| `diagnostics/` | Independent physical/numerical audits and measured limitations |
| `convergence.csv` | Curves used by the response figure; units `(ħ/e)(e²/h)` |
| `figures/` | PNG, PDF and SVG generated from these numerical outputs |

`complete: true` in a command's metadata means that calculation finished and
passed its input/numerical checks. It does not certify material convergence.
The 12×12→18×18 change remains 6.75% of the final value, so these responses
are not converged spin Hall values. The conventional spin-current
definition and finite-band product approximation are recorded in each output.

All local response and tensor files are native VASPBERRY outputs. The large
edge CSV/DAT tables are stored as lossless `.gz` files; their uncompressed
hashes remain those listed by `edge.json`, with archive details in
`text-archives.json`. NPZ files are directly readable by NumPy and the plotter.
For example, `gzip -dc edge-w40/edge.csv.gz > edge.csv` restores the original
table exactly.

The figures use periodic bilinear interpolation for the BZ display and a
10 meV Gaussian broadening for the edge spectrum. Neither operation changes
the original samples, eigenvalues or conductivity integrals. Replot with
`plot.py` from the material directory, following the guide's repository-root
command.
