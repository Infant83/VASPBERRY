# Full Hamiltonian and position operators

These numerical inputs reproduce the **same full Berry connection** as the
Wannier90 3.1.0 reference model. They contain the Hamiltonian in eV and all three
position-matrix components in Å. Both files are required for full J0+J1+J2 AHC.

Restore the nine compressed parts into a new directory:

```sh
python restore.py . /path/to/new-operator-directory
```

The restoration checks the ordered parts and complete matrices and refuses to
overwrite an existing directory. It produces `wannier90_HH_R.dat` and
`wannier90_AA_R.dat`. Together they require about 201 MiB after decompression;
the compressed inputs total about 65 MiB.

Use these files with the [reader described in the toolchain guide](../toolchain/README.md)
and the example's prepared `wannier90.win`. The required settings are
`effective_model = true` and `use_ws_distance = false`: the original model's
pair-dependent distance correction is already incorporated into the matrices.
Do not apply it a second time. The original `transl_inv = false` definition is
preserved, including Hermitization of the full position matrix.

The files keep every nonzero matrix element and include all translations.
No curvature rounding, spin doubling or fitted Hall prefactor is applied.
Two additive rows may describe one matrix element, as supported by the reader;
this preserves precision within its fixed-width input format.

The packaged model has 138 functions and 87 occupied states after removal of a
separately checked deep bundle with Chern number zero. Its 300-step localization
reached the iteration limit. Conversion checks in `validation.json` establish
equivalence to the source model; material accuracy and integration convergence
must be assessed from the example's scientific results.
