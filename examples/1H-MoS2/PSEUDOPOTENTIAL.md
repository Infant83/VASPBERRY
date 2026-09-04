# Pseudopotential provenance

VASP `POTCAR` files are intentionally not distributed with this repository.
Reconstruct the concatenated file from a licensed VASP PAW-PBE potential
library, preserving the species order in `POSCAR`.

- Species order: `Mo S`
- Mo dataset recorded in `OUTCAR`: `PAW_PBE Mo 08Apr2002`, `ZVAL=6`,
  `ENMAX=224.584 eV`
- S dataset recorded in `OUTCAR`: `PAW_PBE S 06Sep2000`, `ZVAL=6`,
  `ENMAX=258.689 eV`
- SHA-256 of the formerly tracked concatenated file:
  `b8f5bccb188e7a8b1082445abe307e62f27d3cd9f9f7b128325dacee3e9192a0`

The hash is provenance only and does not provide the licensed data. See the
official [VASP POTCAR documentation](https://vasp.at/wiki/POTCAR).
