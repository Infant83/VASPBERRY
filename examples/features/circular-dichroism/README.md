# MoS₂: VASP WAVECAR → circularly polarized optical spectra

This tutorial runs the production `vaspberry-gfortran` executable on the **actual
VASP WAVECAR supplied with the repository**. Start from the VASP output, run the
command below, then compare your spectrum and figure with the saved reference.

## 1. Input files

The calculation reads [`WAVECAR`](../../1H-MoS2/KPATH/2.band/WAVECAR), from the
existing [1H-MoS₂ band-path dataset](../../1H-MoS2/README.md). Its
[`INCAR`](../../1H-MoS2/KPATH/2.band/INCAR),
[`KPOINTS`](../../1H-MoS2/KPATH/2.band/KPOINTS),
[`POSCAR`](../../1H-MoS2/KPATH/2.band/POSCAR), and
[`EIGENVAL`](../../1H-MoS2/KPATH/2.band/EIGENVAL) provide the VASP setup and context;
`-cd 2` itself reads `WAVECAR`.

| Property of the supplied WAVECAR | Value |
|---|---|
| System | Monolayer 1H-MoS₂, SOC |
| Stored states | 48 k points, 32 bands, two spinor components |
| Cutoff | 400 eV |
| Path | K → Γ → K′; Γ appears at indices 24 and 25 |
| Occupation | Bands 1–18 occupied; 19–32 empty |
| Input size | 60,521,760 bytes |
| SHA-256 | `33f8546512856d6c04ad0a80454b18ac9b60e2af4b98f2f85b376ec49b1a8d9f` |

No JSON input is needed. This is a replay of the supplied VASP output; it does
not run VASP or generate a model WAVECAR.

## 2. Run VASPBERRY

Run these commands from the repository root, using a new output directory:

```bash
make serial
repo_dir="$PWD"
mkdir -p results
mkdir results/mos2-optical
(
  cd results/mos2-optical
  "$repo_dir/build/vaspberry-gfortran" \
    -f "$repo_dir/examples/1H-MoS2/KPATH/2.band/WAVECAR" \
    -s 2 -kx 48 -ky 1 -cd 2 -if 20 \
    -ien 1 -fen 3 -nediv 201 -sigma 0.05 \
    -theta 0 -phi 0 -o optical > stdout.log 2> stderr.log
)
python3 examples/features/circular-dichroism/run.py \
  --output-dir results/mos2-optical --postprocess-only
```

The last command needs NumPy and Matplotlib. It reads the Fortran output,
produces `summary.csv` and `figure.png`, and checks the saved reference.
Its report records `execution_mode=postprocess_existing_outputs` and marks the
original producer revision `unknown`: a current plotting script cannot certify
which executable created pre-existing files. Its script/reader hashes describe
the local analysis only. The automated replay below records
`execution_mode=calculate_and_postprocess`, the actual command and binary hash,
and the local Fortran source revision; that is how the shipped reference was made.

| Flag | Meaning here |
|---|---|
| `-s 2` | Read both SOC spinor components; this is not a spin-degeneracy multiplier. |
| `-cd 2` | Compute left/right optical spectra from all occupied states. |
| `-if 20` | Include final bands 19 and 20; `-cd 2` automatically uses occupied bands 1–18. |
| `-ien 1 -fen 3 -nediv 201` | Sample photon energies from 1 to 3 eV. |
| `-sigma 0.05` | Gaussian broadening parameter, in eV. |
| `-theta 0 -phi 0` | Light incident along the out-of-plane z direction. |
| `-kx 48 -ky 1` | Bookkeeping for these 48 path points; this does not create a 2D mesh. |

For an automated replay of exactly those steps:

```bash
python3 examples/features/circular-dichroism/run.py \
  --output-dir results/mos2-optical-replay
```

## 3. Outputs and expected result

| Output | What to inspect |
|---|---|
| `CIRC_DICHROISM_W.optical_LEFT_KP1.dat`, `...RIGHT_KP1.dat` | Photon-energy spectra at the first k point, K. |
| `...LEFT_KP48.dat`, `...RIGHT_KP48.dat` | Spectra at K′. There are 48 files per polarization. |
| `...LEFT.dat`, `...RIGHT.dat` | Sums over this path. They are **not full-BZ absorption integrals**. |
| `summary.csv` | 9,648 rows: k index, fractional k, path distance, photon energy, both intensities and selectivity. |
| `figure.png` | Left/right k-resolved spectra and their normalized difference. |
| `result.json` | Input/source/binary checksums, command, validation and reference comparison. |

The lowest K peak occurs at **1.67 eV** on this photon-energy grid. At that point
the printed left/right intensities are **0.0000 / 0.0402** (arbitrary units), so
the selectivity is **−1**. The corresponding K′ peak has selectivity **+1** in
VASPBERRY's channel convention. These are reference values for the supplied
files and selected band window, not a convergence claim for optical absorption.

The current `-cd 2` output convention includes a factor `1/NKPOINTS` in each
k-resolved spectrum, as well as the Gaussian broadening and `1/(photon energy)²`.
Keep that convention in mind when comparing intensity scales from different
path samplings. It cancels in the left/right selectivity ratio.

The postprocessor leaves selectivity blank where the summed intensity is too
small to resolve reliably. In particular, a single fixed band-pair ratio can be
misleading when that transition is spin-forbidden or when the energy ordering of
the conduction states changes along the path. Summing the two selected empty
states before forming the ratio avoids displaying those near-zero denominators.

Reference files:
[numeric CSV](reference/summary.csv) · [provenance and checks](reference/result.json) ·
[raw Fortran outputs](reference/raw-output.tar.gz) · [figure](reference/figure.png)

![MoS2 optical reference](reference/figure.png)

## 4. Apply this workflow to your material

Use the `WAVECAR` from your own VASP calculation and change the flags to match
its spin representation, k points, empty-band range and photon-energy window.
For SOC use `-s 2`; for scalar wavefunctions use `-s 1`. Include enough empty
states to converge the chosen photon-energy range. Use a path for k-resolved
plots, or a suitable full-BZ mesh for an integrated observable. Do not interpret
the sum of path points as a Brillouin-zone integral.

The supplied helper deliberately checks the 18-occupied-band setup: a different
material generally needs the direct CLI command and an adapted postprocessor.
The production routine evaluates bare momentum matrix elements of the stored
pseudo-wavefunctions and prints intensities in arbitrary units; the example is
not an absolute absorption coefficient or a PAW-corrected optical response.

The earlier synthetic numerical check is retained under
[`validation/models/circular-dichroism`](../../../validation/models/circular-dichroism/).
