# Reviewed input templates

These templates define a clean two-stage workflow for generating a Z2-ready
spinor `WAVECAR`:

1. `01_scf/`: converge and write `CHGCAR`.
2. `02_z2_nscf/`: read that charge density with `ICHARG=11` on the full,
   unshifted 12 x 12 x 1 Gamma-centered mesh and write `WAVECAR`.

Use the common [`POSCAR`](POSCAR) and reconstruct `POTCAR` locally from a
licensed VASP PAW-PBE library in the Bi species order. Keep the PAW datasets,
exchange-correlation functional, `ENCUT`, structure, and SOC settings
identical between stages. Test all numerical settings for the system being
studied; these files are starting templates, not universal convergence
parameters.

For the two-atom nonmagnetic test structure, `MAGMOM=6*0.0` explicitly sets
the three initial moment components per atom to zero. Confirm the final
magnetization in `OUTCAR`; the initialization alone does not establish
time-reversal symmetry.

The CHGCAR-only `ICHARG=11` stage is a semilocal-functional recipe. A hybrid
Hamiltonian also depends on regular-mesh orbitals, so a hybrid calculation
requires VASP's separate hybrid-functional workflow and its `WAVECAR`
information.
