# Wannier90 toolchain for the numerical reference

The reference uses [official Wannier90 3.1.0](https://github.com/wannier-developers/wannier90/releases/tag/v3.1.0).
The full numerical operators can be used without VASP, a checkpoint or overlap
files. A small input-initialization correction is required in this release's
`effective_model` reader.

## Build the reader

Extract a fresh official 3.1.0 source tree. From that tree, apply:

```sh
patch -p1 < /path/to/toolchain/effective-reader-modern31.patch
```

Prepare `make.inc` for your compiler and BLAS/LAPACK installation, following
Wannier90's installation guide, with `COMMS = serial`. Then build:

```sh
make -j1 post
```

The patch initializes `num_bands = num_wann` only for an effective model.
Without it, this release leaves that dimension uninitialized and fails an
unrelated disentanglement-input check. Operator readers, interpolation,
derivatives and Berry-response formulas are unchanged. This is a **patched
3.1.0 reader**, and comparisons use unmodified 3.1.0 with the original checkpoint
and overlaps as the numerical reference.

The patches derive from Wannier90 and are supplied under its GNU GPL version 2;
the license text is included. No executable is distributed.

## Regenerate the compact matrices

This step is optional when using the supplied numerical inputs. In a second,
fresh official 3.1.0 tree, apply `export-only-modern31.patch` and build `post`.
Run that export build with the original model checkpoint, overlaps and
eigenvalues, `use_ws_distance = true`, `transl_inv = false`, and a small full-AHC
or bands-plus-curvature task. It writes `wannier90_exact_operators.bin` after
forming the full Hermitian Berry connection; it does not change the response.

```sh
python export_compact_operators.py \
  --input /path/to/wannier90_exact_operators.bin \
  --output /path/to/new-operator-package
```

The converter requires NumPy. The export stream from the documented build uses
little-endian 32-bit integers, 64-bit reals and 128-bit complex values. The
converter checks dimensions, coverage and finite values. It expands the exact
pair-dependent translations, divides by both original degeneracy factors, and
retains all nonzero Hamiltonian and position elements.

The ordinary `_tb.dat` file is insufficient to preserve this reference's
connection definition: its diagonal position terms use a logarithmic overlap,
whereas the selected postw90 definition uses the overlap itself. Extracting
the actual `AA_R` avoids changing that convention. Both the export hook and
reader fix were checked against the original full response.
