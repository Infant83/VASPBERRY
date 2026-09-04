# Pseudopotential provenance

The VASP `POTCAR` is intentionally not distributed with this repository.
Reconstruct it from a licensed VASP PAW-PBE potential library, preserving the
species order in [`inputs/POSCAR`](inputs/POSCAR).

- Species order: Bi
- Dataset: `PAW_PBE Bi 08Apr2002`
- Exchange tag: `PE`
- POMASS: 208.980
- ZVAL: 5.000
- ENMAX: 105.037 eV
- ENMIN: 78.777 eV
- Reference SHA-256:
  `d6b6753ed5db3f0e277fb15e6dbb6699c4bc829850a481068e9e7236faeca489`

The metadata above is recorded from the reference calculation's `OUTCAR` and
from a locally licensed copy. See the official
[VASP POTCAR documentation](https://vasp.at/wiki/POTCAR) for access and usage
requirements.
