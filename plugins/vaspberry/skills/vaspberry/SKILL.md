---
name: vaspberry
description: Use VASPBERRY to inspect VASP outputs and calculate Berry curvature, Chern numbers, supported two-dimensional Z2 invariants and intrinsic charge Hall response with numerical validation and citations. Also install or diagnose VASPBERRY and run its public analytic demonstration.
---

# VASPBERRY for Codex

Use the released VASPBERRY code and its own validation for scientific results. This plugin provides a workflow and local helpers; the calculation engine is a separate checkout. Plugin 0.1.3 was checked with the 2026-10-09 help amendment of source v1.6.6, commit `2e7067b8b8b6d3ab444eef5bc025f4ce5720ca4d`; that is a validation baseline, not a requirement to install an old version. It neither runs VASP nor supplies licensed VASP inputs.

## Find the engine and choose a calculation

Use the user's checkout or requested version when supplied. Otherwise check `VASPBERRY_ROOT`, then the current directory. For a new installation, resolve the **latest stable GitHub release at that time** using [setup.md](references/setup.md), announce the release selected, then install it and record its tag and exact commit. Use published releases rather than assuming the default branch is a tested release. If the user requests an update, check again and install the selected release in a separate versioned directory; preserve previous calculations and modified checkouts. This is an agent workflow when invoked, not a background update service.

Select a Python 3.10+ interpreter and use it for the helper and engine; `python3` below means that interpreter, not necessarily the system default. Resolve helper paths from the directory containing this installed SKILL.md, whether installed through the plugin or as a standalone skill. The examples use `/path/to/installed/vaspberry` for that directory; do not resolve helper paths from the user's working directory:

```bash
python3 /path/to/installed/vaspberry/scripts/vaspberry_agent.py doctor --source /path/to/VASPBERRY
```

`doctor` returns JSON containing the actual source version, commit, dirty state, Python dependencies, available compiler/launcher commands, and existing binaries. `is_target_release_commit` only compares the checkout to the last validated v1.6.6 baseline; false does not itself reject a newer release. Compiler presence is not a successful build. It does not certify input compatibility. If installation is needed, follow [setup.md](references/setup.md). Read the selected release's README, build guide, changelog and actual `--help` before choosing commands. For newer versions, use that release's supported procedures and checks; do not claim this plugin's v1.6.6 verification has tested them.

Identify the requested observable, dimensionality, mesh/path, band or occupied subspace, spin mode and energy reference. Inventory the supplied files; do not infer missing input data. Read only the matching guide in the checkout:

For material calculations, invariant classification or manuscript text, read the bundled [scientific validation guide](references/scientific-validation.md). It defines the terminology, route-specific input conditions and evidence required for a claim. Use the [reporting guide](references/reporting.md) when writing Methods, Results, captions or a scientific delivery.

| Request | Guide and route |
|---|---|
| Fukui Chern number / plaquette flux | `docs/NATIVE_COMMANDS.md`, `docs/VALLEY_TRANSPORT.md`; native `--task chern` |
| Standard Kubo charge curvature / Hall | `docs/WAVEDER_KUBO_PROTOCOL.md`, `docs/POSTPROCESSING.md`; same-run WAVEDER, WAVECAR, INCAR and OUTCAR |
| 2D time-reversal Z2 invariant | `docs/Z2_FUKUI_HATSUGAI.md`, native `--help z2`, `examples/features/z2/README.md`; Fukui–Hatsugai (FH) lattice n-field route |
| Projected-spin-sector topology or Kubo contribution | `docs/SPIN_CHERN.md` or `docs/SPIN_KUBO.md`, native task help; matching OUTCAR spin frame |
| Conventional spin Hall / Wannier | `docs/SPIN_HALL.md`, `docs/OPERATOR_ROUTES.md`; verify the required operator route first |
| First use without VASP data | `demo` command below; explicitly an analytic model |

## Calculate and check

Preserve the user's inputs. Use a new result directory and save exact commands, stdout/stderr, exit status, source identity and input hashes. Use existing native task help and the release's settings examples instead of inventing flags. For the supported **charge Hall INI frontend**, use:

```bash
python3 /path/to/VASPBERRY/tools/vaspberry_post.py check /absolute/path/analysis.ini
python3 /path/to/VASPBERRY/tools/vaspberry_post.py run /absolute/path/analysis.ini
python3 /path/to/VASPBERRY/tools/vaspberry_post.py plot /absolute/path/result
```

This frontend is not a universal Chern/Z2/spin-topology driver: those native tasks use their own commands, output schemas and validity guards from the matching guide. Explicit relative input paths in the INI resolve from the INI directory; an omitted `[run] input_dir` defaults to the invocation working directory. Command-line interface (CLI) paths resolve from the launch directory. For native `--task wavefunction`, POSCAR/EIGENVAL discovery also uses the invocation directory, not automatically the `--input-dir` directory. The INI `plot` command accepts a successful INI run with `run.json`; native tables, pair caches and `workflow.json` use their documented plotting routes. `check` for WAVEDER evaluates the requested spectrum and can cost about as much as `run`. It does not prove convergence. `run.json` retains stage status; preserve it with `logs/` and any failure output. Run Python once when the INI itself configures Message Passing Interface (MPI) parallel execution.

Scientific constraints:

- In v1.6.6 standard charge Kubo requires supported same-run WAVEDER optical input; missing pairs or unsupported producer data must remain errors. The documented supported producer is the standard VASP 5.4.4 longitudinal optical branch, not arbitrary WAVEDER files from every VASP version.
- Use `--kubo-source wavecar` or the equivalent INI setting only when the user knowingly chooses the canonical-momentum approximation. Explain omitted PAW augmentation/nonlocal velocity terms. Never silently fall back to this route.
- Hall integration needs a complete appropriate 2D BZ mesh. A band path or symmetry-reduced mesh is not a total Hall integral. A selected-band/region contribution is not automatically total AHC. Preserve units and energy zero.
- FHS finite-cell flux differs from point Kubo curvature. Do not round a Kubo integral to an integer. A Fukui integer alone does not establish band isolation or mesh convergence.
- Before reporting a 2D Z2 invariant, establish physical time-reversal symmetry, the complete even-rank occupied spinor subspace and an insulating global gap. For v1.6.6, use the documented full even Gamma-centered mesh and retained first-unoccupied sentinel band; require `result_status=PASS` and `reportable_invariant=1`. A time-reversal reconstruction's internal PASS does not establish the raw material's symmetry or gap.
- Direct WAVECAR overlap routes also omit projector augmented-wave (PAW) augmentation. State that limitation for Chern/Z2/projected-spin results as well as the separate canonical-momentum Kubo approximation; do not label these outputs all-electron PAW validation.
- Report execution status separately from numerical validity and material convergence. Show gap/link/coverage diagnostics applicable to the route, and say which mesh, band-window and source-state checks remain undone. Never invent a material result when inputs are absent.
- PROCAR spin/character projections are not conventional spin-current Hall operators. Spin-sector Kubo integrals retain the approximation and convergence qualifications in the source guides.

## Demonstration without proprietary inputs

```bash
python3 /path/to/installed/vaspberry/scripts/vaspberry_agent.py demo --source /path/to/VASPBERRY --output /new/result/directory
```

This calls the release's analytic Qi-Wu-Zhang Fukui runner and saves logs, checksums, figures and citations. It checks Chern values +1, -1 and 0 for the supplied three gapped model masses and empty/full lower-band Hall limits. Clearly label it an analytic model calculation. It does not exercise the WAVECAR parser, PAW treatment, or any material calculation.

## Credit and delivery

Generate software credit from the checkout actually used:

```bash
python3 /path/to/installed/vaspberry/scripts/vaspberry_agent.py cite --source /path/to/VASPBERRY --output /new/citation/directory --method fhs
```

Omit `--method` when no calculation-method reference is selected; supported selectors are `fhs`, `fh-z2`, and `berry-review`. Only select methods used or background actually discussed. Follow the author's README citation guidance and [citation notes](references/citations.md) for the PRB and PRL papers. Keep these scholarly project references distinct from numerical-method references. Consult the checkout's technical report for other method references and add verified references where needed. Use the generated software credit, BibTeX, CFF and provenance files in the delivery. GitHub version/commit credit is sufficient for identifying the software; Zenodo deposition is optional. Never attach the historical v1.0 DOI to a later release.

In the final scientific answer name **VASPBERRY by Hyun-Jung Kim**, link the exact source version/commit and applicable method references, and link generated figures/tables and citation files. Preserve source modifications and approximation/convergence qualifications. Credit is a scientific reporting practice; this plugin does not claim to impose additional license terms.

The citation command prepares references; it is not evidence that a calculation ran. Separate measured outputs, user-supplied claims, analytic expectations and unchecked manuscript placeholders. Report execution, numerical checks, convergence evidence and physical interpretation separately, without turning a program PASS into a validated material claim.

Preserve the engine's existing WAVETRANS credit, including R. M. Feenstra and M. Widom, its source URL, and the description of wavefunction-reading and G-matching routines. Preserve other contributor and third-party notices when installing, copying or adapting code.
