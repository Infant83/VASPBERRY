# Pseudopotential setup

VASP `POTCAR` files are intentionally not distributed with this repository.
Reconstruct the concatenated file from a licensed VASP PAW-PBE potential
library, preserving the species order in `POSCAR`.

- Species order: `Mo S`
- Mo dataset recorded in `OUTCAR`: `PAW_PBE Mo 08Apr2002`, `ZVAL=6`,
  `ENMAX=224.584 eV`
- S dataset recorded in `OUTCAR`: `PAW_PBE S 06Sep2000`, `ZVAL=6`,
  `ENMAX=258.689 eV`

See the official [VASP POTCAR documentation](https://vasp.at/wiki/POTCAR).
