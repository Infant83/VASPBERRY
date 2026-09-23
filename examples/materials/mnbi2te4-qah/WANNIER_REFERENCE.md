# Reproduce the full Wannier Hall reference

This workflow uses the supplied VASP-derived Hamiltonian and position matrices
to evaluate the full **J0+J1+J2** anomalous Hall response. VASP, a WAVECAR,
a Wannier checkpoint and overlap files are not required for this step.

Build the serial [Wannier90 3.1.0 reader with the documented initialization
fix](inputs/wannier/toolchain/README.md). Python requires NumPy. Run the following
commands from the repository root.

## Prepare the integration

```sh
python examples/materials/mnbi2te4-qah/wannier_reference.py prepare \
  --output-dir work/mbt3-wannier-hall \
  --base 80 --refine 9 --radius 0.18 --chunks 8
```

Preparation restores and verifies both full operators, then writes eight
ordinary `wannier90.win` and `kpoint.dat` inputs. To reuse matrices restored
earlier, add `--operators-dir /path/to/restored-operators`.

The base mesh covers the full two-dimensional Brillouin zone at kz = 0. Cells
within 0.18 Å⁻¹ of Γ are replaced by 9×9 submeshes while retaining their area.
All partitions together have weight one; their separate integrals must be
summed without renormalization. The calculation uses zero-temperature
occupations. Three actual chemical potentials span the
central 90% of the sampled DFT gap, using the unchanged VASP energy zero.

## Run

```sh
python examples/materials/mnbi2te4-qah/wannier_reference.py run \
  --output-dir work/mbt3-wannier-hall \
  --postw90 /path/to/patched-3.1.0/postw90.x --workers 8
```

`--workers` defaults to one and accepts up to eight concurrent serial processes.
Each process uses one OpenMP/BLAS thread. Select a worker count appropriate for
your computer. The default time limit is 3,600 seconds per partition; adjust
`--timeout-seconds` when needed. Actual commands, executable and input records,
exit status, elapsed time and logs remain with each partition.

The runner rejects stale outputs and checks both the process exit code and
postw90's normal completion message. A failed or interrupted calculation keeps
its logs. Prepare a new output directory for another attempt.

## Collect

```sh
python examples/materials/mnbi2te4-qah/wannier_reference.py collect \
  --output-dir work/mbt3-wannier-hall
```

After validating every partition, the `results` directory contains matching
`conductivity.csv`, `conductivity.dat` and `conductivity.npz`, plus
`integration.json` describing the external producer, operators and quadrature.
The columns report μ, μ relative to the DFT midgap, the three Cartesian Hall
components in S/cm, and sheet σxy in e²/h. Conversion uses the full simulation
cell height, including vacuum. No additional spin or surface factor is applied.

These are explicitly external full-connection Wannier results. The original
pair-dependent distance correction is already incorporated into the matrices;
the prepared effective-model input therefore sets `use_ws_distance = false`.
Both H and all three position components are retained.

Repeat preparation with other `--base`, `--refine` and `--radius` values to
check integration convergence. Compare the actual gap, band and local-curvature
checks in the [material example](README.md) to assess the finite model itself.
A flat three-point response inside a gap alone does not establish quantization.
