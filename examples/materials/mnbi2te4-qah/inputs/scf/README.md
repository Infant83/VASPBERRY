# Self-consistent density used by the example

`CHGCAR.gz` is the actual converged SOC density used for the reference
calculations. It is supplied so that the NSCF and post-processing steps start
from the same electronic state.

The density uses a 3×3 Γ-centered grid, PBE+SOC, Mn U_eff = 5.34 eV, a 270 eV cutoff,
and alternating +z/−z/+z Mn moments. The fixed 21-atom slab is unrelaxed.
The final electronic tolerance is 10⁻⁵ eV and the net moment is 5.0581 μB.

The original calculation started with `INCAR.initial` and was checkpointed
after 32 iterations. `INCAR.continuation` resumed that state with
`ALGO=Normal`, reaching the tolerance in four further iterations. `OUTCAR`
and `OSZICAR` retain this completed continuation. The checkpoint was a local
resource decision; it is not a special physical parameter.

For a new self-consistent calculation, use the same structure and licensed
PAW datasets and converge the magnetic state and k mesh for the intended
material study. The example's current numerical reference deliberately keeps
the supplied coarse-grid density fixed while testing post-processing.
