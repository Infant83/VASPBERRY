# Re-export the Bi Wannier input from licensed VASP

The supplied raw overlaps already reproduce the recorded localization. This guide is for regenerating those quantities from your own VASP calculation. It requires your licensed VASP source and pseudopotential; neither is distributed here.

The reference used VASP 5.4.4 with the Wannier90 **2.1 serial library** for export, followed by standalone official **Wannier90 3.1.0** for localization and postprocessing. These are separate programs. The spin/velocity producer described in the material tutorial does not automatically enable or correct the Wannier interface.

## Correct the historical spin-labelled projection interface

In the supported historical interface, Wannier90 returns an already spin-labelled list of projection columns. The original VASP loop multiplied that list by two a second time. The correction consumes each returned column once, selects component 1 for spin +1 and component 2 for spin −1, and checks dimensions before writing. Both the VASP `SAXIS` and projection spin axes must be +z. Other axes are explicitly rejected; MPI and the older v1 interface are outside this tested recipe.

The supplied `apply_vasp_wannier_compat.py` stores 18 guarded source edits. It contains the newly written insertion text and short API replacements, without an unchanged licensed source body. It accepts the tested original `mlwf.F` revision only and verifies that its output equals the interface source used for these numerical data. A modified or different revision is rejected.

From the repository root, with the paths below replaced by your licensed source tree:

```sh
case_dir=examples/materials/bi-spin-hall/inputs/wannier-source
python3 "$case_dir/apply_vasp_wannier_compat.py" \
  --source /licensed/vasp-5.4.4/src/mlwf.F \
  --output /licensed/vasp-5.4.4/src/mlwf.F.spin-labelled
```

Keep your original file as a backup, then use the generated file as `src/mlwf.F` in a separate export build. Enable both `-DVASP2WANNIER90` and `-DVASP2WANNIER90v2`, and link the serial `libwannier.a` built from Wannier90 2.1 using a compatible compiler. Perform a clean noncollinear VASP rebuild so all units see the same preprocessor options. This correction changes the Wannier projection interface, not the Hamiltonian, SCF equations or optical/spin operators.

The optional exact projection cache reuses raw overlaps only for identical site, radial function, decay parameter and spin. Individual angular rotations and hybrid combinations are still evaluated separately. Cached and uncached exports were checked to give identical projection matrices. `LREUSE_MMN` must remain false in the recipe below, so the full overlap file is freshly calculated.

## Prepare the full 6 × 6 export

First run the fresh SOC SCF described in the [material tutorial](../../README.md). Use its `CHGCAR`, identical `POSCAR` and identical licensed `POTCAR` for the fixed-density NSCF calculation. The supplied `INCAR.nscf-wannier` keeps the reference physical settings, requests 48 bands and 20 initial diagonalization iterations, and enables the corrected Wannier export. The full automatic Gamma mesh has the order expected by the supplied `.win` file.

From the repository root:

```sh
mkdir bi-vasp-wannier
cp "$case_dir/POSCAR" bi-vasp-wannier/POSCAR
cp /path/to/completed-bi-scf/CHGCAR bi-vasp-wannier/CHGCAR
cp /path/to/completed-bi-scf/POTCAR bi-vasp-wannier/POTCAR
cp "$case_dir/INCAR.nscf-wannier" bi-vasp-wannier/INCAR
cp "$case_dir/KPOINTS" bi-vasp-wannier/KPOINTS
cp "$case_dir/wannier90.win" bi-vasp-wannier/wannier90.win
(cd bi-vasp-wannier && /path/to/corrected/vasp_ncl > vasp.stdout.log 2> vasp.stderr.log)
```

Check normal VASP completion, electronic convergence and the actual band gap before accepting the model input. A process exit code of zero alone does not demonstrate a successful export. Require 16 projection columns, 48 source bands, all 36 k points, complete overlap-neighbor coverage, finite values, and agreement between `.eig` energies and the source VASP energies. The test for projection count and spin axis runs before the expensive projection calculation. Inspect lower conduction states independently when changing the band count or iteration controls.

The reference itself first diagonalized the full source, then used `ALGO=None` solely to export the already converged states. That is a state-preserving export step, not a convergence method. Do not use it with an unconverged input WAVECAR. An assembled explicit-mesh WAVECAR must be reordered into the automatic mesh order before that route; merely replacing its KPOINTS header is incorrect. The fresh automatic-mesh calculation above avoids that ordering issue.

A fresh calculation can select a different gauge inside degenerate bands; compare band spectra, occupied subspaces and physical outputs rather than requiring raw complex matrix entries to match the supplied archive. To reproduce the recorded model exactly, use the supplied raw archives and [localization helper](README.md). To construct a model from your new export, use your new `.amn`, `.mmn` and `.eig` with the same recorded windows and localization sequence, then regenerate the exact effective operators as described there.
