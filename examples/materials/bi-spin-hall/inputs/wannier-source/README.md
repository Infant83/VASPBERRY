# Bi bilayer: reproducible 16-orbital band and edge model

These numerical inputs were generated from the same fresh Bi VASP calculation used for the direct spin Hall example. They contain no licensed pseudopotential, wavefunction, VASP source or executable.

The model includes s, px, py and pz orbitals on each of two Bi atoms, with both spinor components: 16 Wannier orbitals and 10 occupied bands. No deep bands are excluded. The source has 48 VASP bands on a full 6 × 6 mesh. The energy zero is the unchanged VASP eigenvalue zero.

## Use the supplied operators in VASPBERRY

Run the following from the VASPBERRY repository root. All output directories must be new.

```sh
case_dir=examples/materials/bi-spin-hall/inputs/wannier-source
python3 "$case_dir/restore.py" --kind operators --output-dir bi-wannier-operators
python3 tools/vaspberry_kubo.py wannier-import \
  --hh bi-wannier-operators/wannier90_HH_R.dat \
  --aa bi-wannier-operators/wannier90_AA_R.dat \
  --poscar "$case_dir/POSCAR" \
  --spinor-components 2 --spin-multiplicity 1 \
  --energy-reference "Unchanged Bi VASP eigenvalue zero in eV" \
  --output-dir bi-wannier-cache
```

The resulting cache is the input to `wannier-bands` and `wannier-edge` in the [material tutorial](../../README.md). The full position connection is included for a complete operator cache; the edge calculation uses its Hamiltonian. No small matrix entries were discarded.

The operators incorporate the original real-space degeneracy weights and pair-dependent Wigner–Seitz translations. The importer therefore reads them as the effective `HH_R`/`AA_R` pair, not as a conventional `_hr.dat` file.

## Reproduce the localization from actual VASP overlaps

The `raw/` archives contain the actual complex 48-band `.mmn`, `.amn` and `.eig` exports. They are first-principles numerical inputs, not analytic-model parameters. From the repository root, use official Wannier90 3.1.0:

```sh
python3 "$case_dir/reproduce_model.py" \
  --wannier90 /path/to/wannier90.x --output-dir bi-localization
```

The helper checks archives, records inputs and actual exits, and uses one serial CPU. It reproduces the recorded sequence of 1000 iterations followed by three 5000-iteration continuations from the preceding checkpoint. `bi-localization/stage03` contains the final checkpoint and conventional Wannier90 Hamiltonian/position outputs.

### Regenerate the exact effective operators

Build a separate official Wannier90 3.1.0 `postw90.x` with the open-source export-only patch following the [toolchain instructions](../../../mnbi2te4-qah/inputs/wannier/toolchain/README.md). That hook writes the actual full Hermitian position connection without changing it. The supplied `operator-export.win` requests a short bands-and-curvature task with the same model settings and `use_ws_distance = true`, `transl_inv = false`.

From the repository root:

```sh
mkdir bi-operator-export
cp bi-localization/stage03/wannier90.chk bi-operator-export/
cp bi-localization/stage03/wannier90.mmn bi-operator-export/
cp bi-localization/stage03/wannier90.eig bi-operator-export/
cp "$case_dir/operator-export.win" bi-operator-export/wannier90.win
(cd bi-operator-export && /path/to/export-only/postw90.x wannier90)
python3 examples/materials/mnbi2te4-qah/inputs/wannier/toolchain/export_compact_operators.py \
  --input bi-operator-export/wannier90_exact_operators.bin \
  --output bi-regenerated-operators
```

Use the generated `wannier90_HH_R.dat` and `wannier90_AA_R.dat` with `wannier-import` as above. This expansion preserves the pair-dependent translations and degeneracies. Ordinary `_r.dat` or `_tb.dat` must not be substituted for this particular full position connection without validation.

For a licensed VASP → Wannier regeneration, use the [export guide](VASP_EXPORT.md), guarded source instrumenter and explicit full-mesh VASP inputs.

The VASP producer was 5.4.4 with the audited correction for spin-labelled projection columns. This compatibility correction and its producer are recorded explicitly; the data are not represented as outputs of an unmodified historical interface. A different VASP/Wannier version combination must validate its 16 projection columns, complete neighbor coverage and source-state association.

## Accuracy and interpretation

The fixed 16,000-iteration model reaches a last spread change of about 1.57 × 10⁻⁸ Å². It does not meet the requested strict 10⁻⁸ Å² stopping condition. Source frozen-band energies agree within 4 × 10⁻¹² eV. Against 144 directly computed VASP points on the independent 12 × 12 grid, the frontier band error is at most 30.7 meV; the sampled band gap differs by 1.725 meV. The 1000-to-16,000-iteration change of frontier energies is at most 1.404 meV (1.641 meV across all 16 bands). Time-reversal eigenvalue residuals are below 7.2 × 10⁻⁸ eV; no symmetry averaging was applied. Full evidence is in `validation.json`.

This is a finite model for band and ideal-edge illustration, with its interpolation error stated. The conventional spin Hall response is calculated separately from direct PAW spin and full velocity matrices. It is not inferred by attaching mean spins to this model's charge Berry curvature.

For the supplied +60° cell, a conventional adjacent path is Γ(0,0)–M(1/2,0)–K(2/3,1/3)–Γ. Edge weights describe selected Wannier cells and depend on the ideal termination; they are not an all-electron spatial density or a spin texture. Check strip width and combine the edge branches with the independently calculated bulk Z₂ result.
