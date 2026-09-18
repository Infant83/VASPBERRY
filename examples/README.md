# Feature examples

Start with a physical quantity, run a small example, then replace its model
or input with your own system. The [material catalog](materials/) points to
the existing MoS₂ and Bi datasets. Their paths and stored results are retained.

Every feature directory provides an `input.json`, a `run.py`, a README and a
small `reference/` directory containing numerical results and a figure.
The machine-readable index is [`catalog.json`](catalog.json).

| Feature | Input / workflow | Checked result | Reference figure |
|---|---|---|---|
| [Fukui curvature and Chern](features/fukui-chern/) | Public two-band eigenstates → production plaquette-flux and geometric transport helpers | Three gapped phases; flux sums and filled-band Hall sign | [Figure](features/fukui-chern/reference/figure.png) |
| [Kubo curvature](features/kubo-curvature/) | Public analytic matrices → production `demo` and `matrix` commands | Pointwise analytic curvature and finite-mesh Chern integral | [Figure](features/kubo-curvature/reference/figure.png) |
| [Hall and region decomposition](features/hall-valley/) | Point curvature → production `hall` command | Chemical-potential scan, occupation and region sum rules | [Figure](features/hall-valley/reference/figure.png) |
| [Z₂](features/z2/) | Stored Bi field validation; optional WAVECAR recalculation | Schema/status checks and equal half-zone parities | [Figure](features/z2/reference/figure.png) |
| [Circular dichroism](features/circular-dichroism/) | Generated Gamma WAVECAR → current Fortran optical routine | Optical selectivity versus an analytic angular dependence | [Figure](features/circular-dichroism/reference/figure.png) |
| [Wavefunction](features/wavefunction/) | Generated Gamma WAVECAR → current Fortran real-space routine | Real/imaginary parts versus a known plane-wave sum | [Figure](features/wavefunction/reference/figure.png) |

## Run examples

From the repository root, using Python 3.10+:

```bash
python3 -m pip install -r requirements-transport.txt
python3 examples/run_examples.py --list
python3 examples/run_examples.py fukui-chern kubo-curvature hall-valley z2 \
  --output-dir results/python-examples
```

For the two Fortran examples, build the current source first. To run all six:

```bash
make serial
python3 examples/run_examples.py --all --output-dir results/all-examples
# Or select a binary you have already built:
python3 examples/run_examples.py --all --binary build/vaspberry-gfortran \
  --output-dir results/all-examples-second-run
```

Every output directory must be new. The batch runner retains each command,
elapsed time, exit status, stdout and stderr in `run_manifest.json` and log
files. A successful process must also report a matching `feature_id`, schema
version 1 and `status=PASS` in `result.json`. Failures are retained, reported
and produce a nonzero exit code.
Each feature records its numerical checks and input/source provenance in
`result.json`; `summary.csv` and `figure.png` are compact presentation outputs.
Larger arrays and generated wavefunctions remain in the run directory.

Use an individual `run.py --help` for its input and optional recalculation
arguments. Follow the feature README to customize `input.json`; generic
calculation tools remain independent of these small tutorial runners.

The CI `feature-examples` job runs these same six defaults and retains generated
results and logs as an artifact, including when a workflow fails.

## What the results demonstrate

Five default workflows calculate new results from public analytic or synthetic
inputs. **Z₂ defaults to checking and plotting a stored Bi result.** It is
listed separately as `stored_reference_validation`; the feature's `--wavecar`
mode performs a new Fortran calculation when the full input is available.

The examples test numerical routines, conventions and data interchange.
They do not establish k-mesh, band-window or DFT convergence for a material.
The Hall example's regions are illustrative partitions of a model Brillouin
zone; a physical valley definition must come from the system being studied.
The optical example reports dimensionless selectivity, not an absolute
absorption rate. Wavefunction plotting is restricted to Gamma.

Fukui plaquette flux and Kubo point curvature are different output quantities;
see the [output format](../docs/OUTPUT_FORMAT.md). For old Kubo data, use the
explicit [normalization migration workflow](../docs/MIGRATION.md).

## Material files and reusable plotting

- [Material catalog](materials/): available MoS₂/Bi inputs and reproduction limits.
- [1H-MoS₂](1H-MoS2/): existing full-mesh reference maps.
- [MoS₂ band path](1H-MoS2/KPATH/): line sampling, unsuitable for full-BZ integrals.
- [Bi](Bi_Z2/): reviewed templates and Git LFS wavefunction.
- [Standalone Hall CSV plotter](kubo/): retained at its original path.

No VASP `POTCAR` is distributed. No new material-specific or unpublished
collaboration data is needed by the default workflows.
