# Actual VASP input files

The user tutorials operate on the public VASP outputs below. They do not need
a JSON model input. JSON files written by the runners are provenance and
validation records.

| Input | Location | Sampling / content |
|---|---|---|
| MoS₂ WAVECAR | [KPATH/2.band/WAVECAR](1H-MoS2/KPATH/2.band/WAVECAR) | 60,521,760 bytes; SOC, 48 path points, 32 bands |
| Matching MoS₂ geometry | [POSCAR](1H-MoS2/KPATH/2.band/POSCAR) | Needed for the real-space output header |
| Matching MoS₂ eigenvalues | [EIGENVAL](1H-MoS2/KPATH/2.band/EIGENVAL) | Band/occupation and Gamma-point selection |
| Bi WAVECAR | [Bi_Z2/WAVECAR](Bi_Z2/WAVECAR) | Git LFS: 200,421,600-byte SOC payload, 144 full-mesh points, 18 bands |
| Bi archived eigenvalues | [EIGENVAL](Bi_Z2/archive-2016-run/EIGENVAL) | Independent view of the original band-edge gap |
| Bi new-run templates | [inputs](Bi_Z2/inputs/) | Reviewed SCF then full-mesh NSCF templates for a new calculation |

The Bi payload SHA256 is
`a8d81854f2efc561938e478dde1be29a17ccc95d5d37325c122b9d9e82fa0838`.
Use `git lfs pull --include='examples/Bi_Z2/WAVECAR'`, or:

```bash
python3 examples/fetch_inputs.py bi --output-dir results/inputs/bi
```

The helper fetches the public file pinned to repository commit
`a692e21482d24767859d02b7885dd70b234c6a2c`, verifies both size and SHA256, and
retains a failed download record if verification fails. A completed output
folder contains an actual binary WAVECAR and `download.json`. A downloaded
file can be checked/copied without network using `--source /path/to/WAVECAR`.

Each tutorial's `reference/result.json` records hashes of the exact inputs,
commands and software used for its reference. Material files are shared across
feature tutorials; large wavefunction binaries are not duplicated in each
feature directory.

## Reproduction boundary

The full-mesh MoS₂ WAVECAR for its historical curvature map is not supplied.
That map remains an archived result, not a new full-mesh tutorial. Do not use
the supplied line-mode WAVECAR in its place.

The Bi original SCF provenance is incomplete. Supplied VASP outputs allow the
postprocessing calculation to be repeated, but the archived CHGCAR is not
claimed to come from the reviewed new-run template. No POTCAR is distributed;
follow the material pseudopotential provenance using a licensed VASP library.
