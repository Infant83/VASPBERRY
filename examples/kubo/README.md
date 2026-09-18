# Public Kubo/Hall example and CSV plotting

For a complete input → calculation → checked result workflow, start with the
[Kubo curvature](../features/kubo-curvature/) and [Hall/region](../features/hall-valley/)
feature examples. This directory retains the reusable CSV plotter.

The [Kubo transport guide](../../docs/KUBO_TRANSPORT.md#public-demonstration)
provides a complete QWZ analytic-model workflow: `demo` → `matrix` → `hall`.
It needs no VASP input, licensed potential or material-specific setting. Run
those commands from the repository root to create `results/qwz-hall/`.

This standalone plotting example reads the resulting standardized CSV using
only Python's standard library and the existing Matplotlib dependency:

```bash
python3 examples/kubo/plot_hall.py \
  --input results/qwz-hall/conductivity.csv \
  --output results/qwz-hall/hall.png
```

Select channels and an observable for a different calculation or figure:

```bash
python3 examples/kubo/plot_hall.py \
  --input results/qwz-hall/conductivity.csv \
  --regions total --bands 0 1 --temperatures 0 \
  --quantity delta_sigma_e2_over_h --relative-mu \
  --title "Change from the reference chemical potential" \
  --output results/qwz-hall/hall-change.pdf
```

`--regions` accepts names already present in the CSV, including a named
regional contrast produced by `hall --difference`. The plotter does not
invent valley boundaries or insert a factor of one-half. `--bands 0` means
the sum of represented bands; individual band rows require the producer's
`--band-resolved` option. Temperatures default to all available values.

Available quantities are `sigma_e2_over_h`, `delta_sigma_e2_over_h`, `sigma_S`,
`electrons_per_cell`, and `delta_electrons_per_cell`. `sigma_S` is sheet
conductance in siemens, not bulk conductivity in S/m. Occupation columns in
regional-difference channels are signed differences. PNG, PDF and SVG output
can be selected by filename suffix; existing files are not overwritten.

If a same-stem JSON sidecar is present, the figure annotates partial-band and
experimental/approximate-operator status. Keep the CSV/NPZ/JSON and source
provenance with the figure. Plotting does not revalidate the physics, recompute
curvature, smooth the data, enforce an integer or establish convergence. The
QWZ example's sign and units are model checks, not a material benchmark.
