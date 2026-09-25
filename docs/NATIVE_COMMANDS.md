# Native Fortran command reference

`build/vaspberry` reads VASP WAVECAR and writes numerical results for subsequent
analysis and plotting. It does not require Wannierization. Start with
`make serial`, then `build/vaspberry --help`; see [build details](BUILD.md).
The [hands-on guide](HANDS_ON.md) provides complete material workflows.

## A runnable first command

From the repository root, use the included MoS₂ path WAVECAR:

```bash
make serial
VB_ROOT="$PWD"
mkdir results-native-path-01
cd results-native-path-01
"$VB_ROOT/build/vaspberry" --task kubo \
  --wavecar "$VB_ROOT/examples/1H-MoS2/KPATH/2.band/WAVECAR" \
  --spinor 2 --bands 1:18 --bundle 1 --curvature-csv KUBO.csv
cd "$VB_ROOT"
```

This writes 48 occupied-bundle curvature rows to `KUBO.csv`. It uses every
stored empty band as an intermediate state; it does not integrate a path.
The [curvature example](../examples/features/kubo-curvature/) supplies the
matched full-mesh calculation and figure commands. To repeat this command,
choose another fresh directory.

## Choose the result

All commands accept `--wavecar PATH --spinor 2` for SOC spinors. The parameters
below use the MoS₂ occupied rank 18 as an example, not a universal band count.

| Task | Additional arguments | Numerical result and required sampling |
|---|---|---|
| `chern` | `--mesh 12,12 --bands 1:18` | Fukui occupied-subspace flux in `BERRYCURV.dat`; Chern number in its header. Full periodic 2D mesh. |
| `z2` | `--mesh 12,12 --bands 1:18` | `Z2_FIELD.csv` and `NFIELD.dat` after PASS. Full even Γ-centered mesh, SOC, time-reversal symmetry, occupied rank even and gap open; [full requirements](Z2_FUKUI_HATSUGAI.md). |
| `kubo` | `--bands 1:18 --bundle 1 --curvature-csv KUBO.csv` | Bundle trace Ωxy in Å², one row per source k/spin. Mesh or path. External gaps must exceed 1e−5 eV. |
| `kubo` | `--bands 18 --curvature-csv BAND18.csv` | Individual-band Ωxy, energy and minimum gap, plus legacy DAT companions. Inspect isolation before interpretation. |
| `kubo-pairs` | `--pairs-csv PAIRS.csv` | Every source-band pair's undivided momentum numerator, energies and k coordinates; all yz/zx/xy components. No band selector. |
| `kubo-integral` | `--mesh 12,12 --bands 1:18 --bundle 1 --curvature-csv KUBO.csv` | Same bundle CSV plus a finite-grid integral in stdout. Requires the full mesh; this diagnostic is not an occupation-weighted transport scan. |
| `optical` | `--mesh 12,12 --bands 18:19 --theta 0 --phi 0` | Selected transition 18→19 circular selectivity in `BERRYCURV.dat` (legacy filename). |
| `spectrum` | `--mesh 12,12 --bands 1:20 -ien 1 -fen 3 -nediv 201 -sigma 0.05` | Broadened left/right transition spectra in `OPT_TRANS_RATE_LEFT.dat` and `OPT_TRANS_RATE_RIGHT.dat`. Energies/broadening in eV. [Optical workflow](../examples/features/circular-dichroism/). |
| `wavefunction` | `--wavefunction-band 18 --kpoint 24 --real-grid 24,24,64 --imaginary 1` | CHGCAR-style `PARCHG-W-K024-E018-SPIN1` and its `-IM-SPIN1` companion. Index 24 is Γ in the included 48-point path only. Matching POSCAR/EIGENVAL must be in the working directory. [Wavefunction workflow](../examples/features/wavefunction/). |
| `velocity` | `--bands 18` | Canonical-momentum velocity expectation in `VEL_EXPT.dat`, x/y components in m/s; a single-band diagnostic. |

`kubo-line` is a synonym of `kubo`; both preserve the supplied k-point order.
For `kubo`, `--bands 1:18` without `--bundle 1` means eighteen individual-band
curvatures, not the degenerate occupied-subspace trace. Optical band endpoints
instead select the initial/final transition; `spectrum` uses the selected band
window and source occupied boundary. `velocity` uses the first selected band,
so specify a single band.

## A small set of common arguments

| Argument | Meaning |
|---|---|
| `--task NAME` | One calculation. Explicit task names reject conflicting legacy task flags. Default is Fukui `chern`. |
| `--wavecar PATH` | Input WAVECAR, default `WAVECAR` in the working directory. Quote paths containing spaces. |
| `--spinor 2` / `--spinor 1` | Two-component SOC/noncollinear state / scalar state. This is not a spin-degeneracy factor. Always state it explicitly for your data. |
| `--mesh NX,NY` | Mesh dimensions for mesh-based algorithms. It does not generate k points, interpolate, or convert a path into a mesh. |
| `--bands FIRST:LAST` / `--bands N` | Inclusive one-based band range / one band. Explicit selection avoids guessing the occupied boundary. |
| `--bundle 1` | Trace curvature of the selected group, allowing internal degeneracies; requires `--curvature-csv`. |
| `--curvature-csv PATH` / `--pairs-csv PATH` | Exact CSV output filename. Parent directory must already exist. |
| `--output LABEL` | Legacy output naming label; **not an output directory** and not a replacement for a CSV path. |

Every option takes a separate value: use `--mesh 12,12`, not
`--mesh=12,12` or `--mesh 12 12`. Lists contain no spaces. Modern names and
legacy flags can be mixed where needed; optical energy-grid controls above
retain their short names. There is no need to learn every historical flag.

GNU native builds expect ordinary byte-record, single-precision complex
WAVECAR coefficients (`RTAG=45200`). Older files whose RECL is in compiler
words are not automatically converted by the native reader. Retain the VASP
producer/compiler provenance and verify compatibility before a large run.

## Where outputs go

Use a fresh working directory per native calculation and an absolute WAVECAR
path. Legacy DAT and wavefunction outputs can replace existing files. Kubo
CSV exporters instead refuse an existing filename; they never append another
calculation to it. A complete pair export ends with `result_status=PASS`;
bundle and Z₂ validation status must also be checked before use.

Legacy labels are retained for compatibility: `--output sample` produces
`BERRYCURV.sample.dat` for Fukui, `sample.dat` for optical selectivity,
`CIRC_DICHROISM_W.sample_LEFT/RIGHT.dat` for spectra, and
`VEL_EXPT.sample.dat` for velocity. Individual-band Kubo adds `.EIG-N`, and
its sum has the base label. Scalar `ISPIN=2` DAT outputs add `.UP`/`.DN`.
Wavefunction filenames are selected by k and band; `--output` does not rename
them. Bundle/pair exports use only their explicit CSV path.

The native program currently stores most paths in 256-character fields.
Prefer short working paths. The CSV interfaces are the preferred precision
and interchange route for Kubo results; legacy DAT tables have fewer printed
digits. [Output schemas](OUTPUT_FORMAT.md) define units, columns, index bases,
operator assumptions, and the NPZ/JSON caches used by postprocessing.

## Charge transport and reusable figures

The main charge-Hall workflow is:

```text
VASP full-mesh WAVECAR
  -> Fortran --task kubo-pairs --pairs-csv PAIRS.csv
  -> Python import-pairs -> pairs.npz + pairs.json
  -> Python pair-hall -> conductivity.csv / .dat / .npz + .json
  -> plot_hall.py or your preferred plotting program
```

Pair export does not use source occupations or divide by energy gaps;
k-dependent metallic/smeared occupations therefore need no manual `-ne`
override. Postprocessing supplies μ, temperature, spin multiplicity and
regions. Its complete, uniform 2D mesh and degeneracy checks still apply.
Once the cache exists, varying μ, temperature, regions or the retained
intermediate-band cutoff does not require another native calculation.
See the [three-stage commands](KUBO_TRANSPORT.md#native-pairs-to-charge-hall)
and [actual MoS₂ Hall example](../examples/features/kubo-hall/).

For atom/layer/orbital and chosen-axis spin character, retain the matching
SOC `PROCAR` and `OUTCAR` as well. The separate
[PROCAR workflow](../examples/features/procar-character/) accepts user-defined
groups, combines their projections with the saved canonical pair cache for
selected isolated bands, and writes charge-attribution tables and plots.
No additional native export or Wannier model is needed.

Native WAVECAR Kubo uses canonical momentum of pseudo-wavefunctions. It
does not supply missing PAW/nonlocal/SOC velocity terms. These outputs support
the documented intrinsic charge-response approximation, not a complete
longitudinal conductivity, relaxation-time transport or spin-Hall calculation.
For additional physical operators and spin Hall, use the separately documented
[operator routes](OPERATOR_ROUTES.md). Wannier is an optional model route.

Regenerate historical native `-vel 1` outputs: the previous velocity conversion
omitted the factor c² needed with electron mass expressed in eV/c², and its
extrema used an invalid k index. Rounded legacy zero values cannot be repaired
afterward. The corrected diagnostic prints scientific-notation m/s values and
still represents only canonical momentum, not the full material velocity.

For MPI, run `make mpi` and replace the executable with
`mpiexec -n 4 "$VB_ROOT/build/vaspberry-mpi"`; calculation arguments stay the same.
