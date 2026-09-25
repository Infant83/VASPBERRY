# Bi potential used in this example

Use the licensed **PAW_PBE Bi 08Apr2002** dataset, with five valence electrons
per atom and recommended ENMAX 105.037 eV. The calculation uses ENCUT 400 eV.
Two atoms give ten occupied SOC spinor states. The spinor spectrum already
enumerates spin; no factor of two is added.

Obtain POTCAR through your VASP licence. It is not redistributed here. The
preparation helper checks the potential against the source record to avoid
silently mixing the supplied density with a different PAW dataset. For a
different dataset, regenerate the SCF density and reconverge the calculation.
Machine-readable input identities are in `inputs/provenance.json`.
