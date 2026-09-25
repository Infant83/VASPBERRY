# Fresh Bi SCF density

This density was generated for the spin Hall example using the supplied fixed
two-atom geometry, PBE+SOC, 400 eV, Accurate precision, a Γ-centered 12×12×1
mesh and 32 spinor bands. Crystallographic symmetry reduces the SCF mesh to
19 points. The nonmagnetic calculation reached EDIFF 1e−8 eV in 22 electronic
iterations. Subsequent operator calculations disable symmetry and retain
the full mesh.

`prepare_vasp.py` restores the split compressed CHGCAR and checks its
integrity. To regenerate it, choose `--stage scf` with your licensed matching
Bi POTCAR. The historical Bi_Z2 fixture remains a separate calculation;
its old density was used only for exporter development tests.

The geometry is fixed, without relaxation or substrate. Numerical response
convergence at this geometry does not establish structure or functional
convergence for a material prediction.
