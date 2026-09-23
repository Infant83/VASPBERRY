# PAW dataset specification

Construct the licensed PBE `POTCAR` in **Mn, Bi, Te** order:

| Element | VASP dataset title | Valence electrons per atom |
|---|---|---:|
| Mn | PAW_PBE Mn 06Sep2000 | 7 |
| Bi | PAW_PBE Bi 08Apr2002 | 5 |
| Te | PAW_PBE Te 08Apr2002 | 6 |

The cell contains Mn₃Bi₆Te₁₂ and 123 valence electrons. With spin–orbit
coupling, the insulating occupied bundle is bands 1–123; no additional factor
of two is applied.

`prepare_vasp.py` checks that the concatenated file matches the reference
datasets before using the supplied charge density. The datasets and VASP
executables are licensed software and are not included in this example.
