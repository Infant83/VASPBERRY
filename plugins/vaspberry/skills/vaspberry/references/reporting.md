# Scientific delivery and manuscript text

Use [scientific-validation.md](scientific-validation.md) for definitions and claim conditions. Write in the user's requested language, with conventional scientific names rather than new labels for code-specific stages. Expand an abbreviation at first use in each standalone text. Use standard symbols C, Z2, Ω, σxy, μ and T with their definitions, units and conventions. “Supported,” “tested” and “converged” must each name their actual scope.

## Build the report from evidence

Use the run's numerical tables, settings, logs and provenance as the source. Separate what was executed, what was only supplied by the user, what is an analytic expectation and what remains unchecked. The `cite` command creates references without executing a material analysis; do not describe bibliography preparation as a completed calculation.

A scientific delivery should let the reader identify:

- The requested observable, input/model and calculation route, with links to the actual data, commands, figures and logs.
- Execution outcome, numerical diagnostics, convergence evidence and physical interpretation as distinct statements. List consequential missing evidence beside the result it limits.
- Software version/commit, modifications, input hashes, relevant build/environment and source-state settings. “Latest release” alone is not a reproducible identity; selected-file hashes alone are not a complete modified source snapshot.
- What the plotted/table value represents: selected bands or full occupied response; valley/region or full BZ; sheet or bulk units; curvature or flux; raw σ or Δσ; operator/PAW approximation.
- Software credit, applicable method citations and the author's recommended PRB/PRL research references, with the roles explained in [citations.md](citations.md).

Do not introduce missing parameter values to make a paragraph look complete. Ask for missing essential data or mark explicit placeholders in a draft. A Methods draft with placeholders and a qualified result may be useful; it is not a ready-to-submit manuscript.

## Methods content to retain

For real VASP-based work, retain the system/geometry and dimensionality, electronic-structure producer/version, exchange-correlation functional, PAW dataset identifiers, cutoff, source-state and mesh settings, SOC/magnetic state, occupation and band/subspace selection. Then describe the actual VASPBERRY algorithm/operator route, mesh and reciprocal orientation, energy reference, temperature, required partner bands, regularization/tolerances and approximation. Record the performed convergence/reference checks and software identity. Fields that do not apply to the chosen route need not be added.

For an analytic demonstration, state the Hamiltonian/model parameters and their units, selected band/subspace, grid, connection convention, production routines exercised, measured checks and limiting scope. Do not call analytic states VASP wavefunctions or describe the demo as a material validation. Read the actual input/result files; do not assume that future releases retain the same model parameters.

Example structure, to fill only from actual records:

> We evaluated [observable] for [system/model] using VASPBERRY [version, commit] with [algorithm/operator route]. The calculation used [occupied/selected subspace], [mesh and geometry] and [energy/occupation conventions]. [Relevant approximation and validity evidence]. Refinement from [actual settings] changed [measured quantity] by [measured amount]; [checks still missing] remain outside the present validation.

This is an optional drafting aid, not text to paste with invented entries. Link supporting tables/logs in the delivery and retain fuller reproducibility information in the supplementary material when preparing a paper.

## Results, captions and claim strength

Report values with units and precision supported by numerical evidence. Preserve unrounded integrals in data; do not round a Kubo result into a topological integer. For a valid integer invariant, distinguish its class from finite-grid residuals and show the diagnostics that justify interpreting it.

Use “the executed model check passed,” “the sampled invariant is numerically self-consistent,” or “the response was stable over the reported mesh sequence” when that is the actual evidence. Reserve “the material is a [phase]” for a supported physical classification with the relevant symmetry, gap, occupied-subspace and convergence evidence. Avoid “proves,” “exact,” “fully validated” or “publication-ready” when the underlying assumptions or checks are missing. A computation does not establish experimental observation.

Each standalone caption should identify the observable/component, system, axes/units, energy reference or k coordinates, mesh/band/region selection, occupation parameters and relevant approximation. Define special symbols in the caption or its immediately accessible legend. State whether maps show computed cell flux, normalized cell-average curvature or point curvature; distinguish calculated samples from interpolation.

Keep numerical convergence and model/approximation limitations separate. If an independent comparison was not done, say so without presenting agreement between quantities sharing the same underlying implementation as independent confirmation.

## User guidance when a task cannot be completed

Identify the missing artifact or failed physical/numerical condition, explain what claim it prevents, and provide the smallest supported next step. Do not silently change an observable, substitute a demonstration for a requested material calculation, or use a different operator route. Present the analytic demo as an explicit option when no material inputs are available. Keep failed commands and diagnostics accessible so the user can reproduce or repair them.
