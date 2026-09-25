# Reference PAW datasets

Use the following licensed VASP PBE datasets, concatenated in the POSCAR
species order **Mo, then S**:

| Species | Dataset title | Valence charge |
|---|---|---:|
| Mo | `PAW_PBE Mo 08Apr2002` | 6 |
| S | `PAW_PBE S 06Sep2000` | 6 |

The preparation helper checks the supplied combined POTCAR against the
reference identity recorded in [inputs/provenance.json](inputs/provenance.json).
The potential files are not distributed. All four structures use these same
datasets; a neutral SOC monolayer has 18 occupied spinor states and a bilayer
has 36. This count belongs to these potentials and structures. Determine
NELECT and the occupied group again when applying VASPBERRY to other inputs.
