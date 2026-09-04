# Legacy pre-v1.1.1 outputs

These files were generated with the pre-v1.1.1 implementation and are kept
only as regression evidence.

`NFIELD.dat` reports top = -1 and bottom = 0, which disagree modulo 2. That
disagreement is a known symptom of the incomplete legacy calculation,
including the missing reciprocal-G shift at nonzero time-reversal invariant
momenta. The legacy Z2 plots and data are not validated results and must not be
cited as the material's invariant.

The pointwise n-field is gauge- and logarithm-branch-dependent. The failure is
the inconsistent half-zone parity together with the time-reversal diagnostics,
not an expectation that every raw n-field value be pointwise symmetric.

`BERRYCURV.dat` and `BERRYCURV.eps` came from the ordinary Berry-curvature
path and are not direct evidence of the Z2/TRIM defect. They are co-located
here so that every result and plotting helper from the older run remains a
clearly separated historical bundle. `auto-Z2.sh`, `contour.py`, and
`gnuZ2.gp` are archived for provenance and are not supported run scripts.

Use [`../run_z2.sh`](../run_z2.sh) to recalculate with the current code.
