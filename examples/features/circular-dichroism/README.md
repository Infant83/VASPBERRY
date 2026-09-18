# Circular optical selectivity

Run the production Fortran `-cd 1` calculation on a **640-byte synthetic WAVECAR**, then compare its angle dependence to an analytic answer. No VASP installation, PAW dataset, or downloaded material file is needed.

![Actual Fortran selectivity and analytic curve](reference/figure.png)

## Run

From the repository root:

```sh
python3 -m pip install -r requirements-transport.txt
make serial
python3 examples/features/circular-dichroism/run.py --output-dir /tmp/vaspberry-circular-example
```

The output directory must not exist. Use `--binary /path/to/vaspberry` to select another compiled production executable. The script reads [input.json](input.json); `--input /path/to/input.json` selects a modified copy. You may change the finite, distinct `theta_degrees` samples between 0° and 180°. The fixed scalar Gamma fixture and `phi_degrees = 0` are validated; unsupported changes are rejected. The script does not build or modify the executable and leaves generated inputs and raw Fortran outputs in the output directory.

For every angle `0, 15, …, 180`, the script creates a numbered `theta-NNN/` case directory and executes this command there (replace `$BINARY` with the absolute executable path):

```sh
"$BINARY" -f ../WAVECAR.synthetic -s 1 -kx 1 -ky 1 -ii 1 -if 2 -cd 1 -theta 0 -phi 0 -kp 1 -o selectivity
```

Only the `-theta` value changes between those 13 commands. Each command and exit code is recorded in `result.json`.

## Why the answer is known

The scalar Gamma-point basis has a cubic cell of side `2π Å` and seven plane waves, with cutoff 5 eV. Only `G = 0, +x, +y` have nonzero coefficients:

```text
valence:    (1, 1, 1) / √3
conduction: (1, ω, ω*) / √3,  ω = exp(2πi/3)
```

These two states are normalized and mutually orthogonal before complex64 storage. With the code's `p_cv = <c|p|v>` convention and azimuth `φ = 0`, their optical selectivity is

```text
η(θ) = (|P+|² − |P−|²) / (|P+|² + |P−|²)
     = √3 cos(θ) / [1 + cos²(θ)].
```

The checked result is `+0.866025` at 0°, zero at 90°, and `−0.866025` at 180°. All 13 standard-output values agree with the analytic formula within `4.61e-7`. The legacy `.dat` column has four decimal places; its separate tolerance is `5.1e-5`. The script checks both representations rather than inventing extra precision.

## Outputs and scope

- `summary.csv`: measured selectivity, analytic value, and both rounding errors.
- `figure.png`: actual measured points with the analytic curve.
- `result.json`: pass/fail status, commands, input/source/binary hashes, and errors.
- `theta-NNN/selectivity.dat`, `stdout.log`, `stderr.log`: original production outputs.

The committed [reference result](reference/result.json) and [summary](reference/summary.csv) come from a real run of VASPBERRY 1.3.0. Compressed original `.dat` and stdout files are retained in `reference/raw/`; executable paths in the reference JSON are expressed relative to the repository root for portability.

This is a single synthetic scalar transition, **not a material prediction or an absolute absorption rate**. Its coefficients and energies are authored model data rather than converged eigenstates of a material Hamiltonian. The ratio uses bare momentum and cancels common amplitude factors; it does not test PAW/nonlocal optical corrections, occupations, photon-intensity prefactors, spinor/SOC behavior, or Brillouin-zone integration. The signs follow the program's `P+`/`P−` convention and should not be relabeled as an experimental helicity convention without specifying the viewing direction.
