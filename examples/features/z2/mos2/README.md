# MoS₂: a trivial Z₂ insulator with valley Berry curvature

The native Fukui–Hatsugai calculation gives **Z₂ = 0** for occupied SOC
bands 1–18 of monolayer 1H-MoS₂ on the actual full 12×12 mesh. Both half-zone
integer sums are zero. This contrasts with the fresh Bi bilayer's **Z₂ = 1**;
the [paired example](../comparison/) displays both native fields.

## 1. Reuse the full-mesh VASP preparation

Use the same full-mesh WAVECAR as the
[MoS₂ Fukui-curvature example](../../fukui-berry-curvature/).
It contains 26 stored spinor bands at a 400 eV cutoff, generated from the
supplied fixed SCF density with `ICHARG=11` and `ISYM=-1`. A band-path
WAVECAR cannot replace this full periodic mesh.

From the repository root:

```sh
python3 examples/features/fukui-berry-curvature/prepare_vasp.py \
  --potcar /path/to/licensed/Mo-S/POTCAR \
  --output-dir results/mos2-fullmesh-vasp
cd results/mos2-fullmesh-vasp
/path/to/vasp_ncl > vasp.stdout.log 2> vasp.stderr.log
cd ../..
```

If this matching calculation is already complete, reuse its WAVECAR.
The [VASP input guide](../../fukui-berry-curvature/inputs/README.md) records
the structure, charge density and PAW datasets. VASP and licensed POTCAR
files must be supplied locally; the Z₂ step needs no Wannier model.

## 2. Run native VASPBERRY

```sh
make serial
repo_dir="$PWD"
mkdir results/mos2-z2
(
  cd results/mos2-z2
  "$repo_dir/build/vaspberry" --task z2 \
    --wavecar "$repo_dir/results/mos2-fullmesh-vasp/WAVECAR" \
    --mesh 12,12 --spinor 2 --bands 1:18 --output NFIELD \
    > fortran.log
)
```

These task aliases refer to the current repository. For an MPI build, use
`make mpi` and `mpiexec -n 4 "$repo_dir/build/vaspberry-mpi"` with the same
arguments. All numerical links and the integer field are calculated by the
Fortran executable; Python is only needed for the figure and row checks.

| Native reference quantity | Result |
|---|---:|
| Occupied bundle | 1–18 |
| Half-zone sums | 0 and 0 |
| Half-zone parities | 0 and 0 |
| Z₂ | **0** |
| Sampled direct / global gap | 1.67355 / 1.67355 eV |
| Minimum link singular value | 0.713968 |

Require `result_status=PASS`, `reportable_invariant=1` and agreeing
half-zone parities. The [native CSV](reference/Z2_FIELD.csv),
[text output](reference/NFIELD.dat) and [result record](reference/result.json)
retain the source association and diagnostics. Follow the
[comparison plotting command](../comparison/#plot-the-completed-fields)
to place this result beside the fresh Bi field without rerunning VASPBERRY.

## Additional numerical checks

The full occupied bundle also gives **Z₂ = 0**, with half-zone sums 0 / 0,
in the [24×24 native field](reference/mesh-24x24/Z2_FIELD.csv).
That source stores 60 bands, while the
12×12 source stores 26; both select occupied bands 1–18 and retain the same
fixed SCF density. This checks invariant stability under those input changes.
It does not isolate a pure k-mesh change or establish transport convergence.
The minimum link singular value changes from 0.713968 to 0.741853, and
the largest plaquette phase decreases from 0.387476 to 0.121679 rad.

The 12×12 and 24×24 serial/MPI field CSVs are byte-identical at each mesh. A separate raw
occupied-subspace check tests time-reversal correspondence before native
reconstruction; the archived result metadata retain its definition and values.

## Interpretation

MoS₂ has C = 0 and Z₂ = 0 in this example, while Bi has C = 0 and Z₂ = 1.
The finite valley-contrasting Berry curvature and regional response shown
in the MoS₂ examples do not require nontrivial Z₂ topology. Individual
positive or negative n-field tiles also do not imply a nontrivial invariant:
their half-zone parity matters, and the local tile pattern is gauge and
branch dependent. No smoothing or integer adjustment is used.

The native consistency checks validate the time-reversal reconstruction;
they do not independently establish raw-input time-reversal symmetry or
physical convergence. For another material, verify the insulating gap,
occupied subspace and time-reversal symmetry, and refine the mesh. These
WAVECAR overlaps omit PAW augmentation; see the
[method guide](../../../../docs/Z2_FUKUI_HATSUGAI.md).
