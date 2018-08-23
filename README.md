# VASPBERRY
Berry curvature and Chern number calculations with the output (WAVECAR) of VASP code.
VASPBERRY is written for the post-processing pur pose of the VASP outputs, i.e., WAVECAR the Bloch wavefunction infor mation. VASPBERRY can compute Berr y cur vature and Cher n number via Fukui's method [S ee J. Phys. S oc. Jap. 74, 1674 (2005)]. In addition Circular dichroism also can be evaluated. Since it directly reads the wavefunction coefficients, one can also obtain real space wavefunction character psi_nk(r) by simple command.

# Compile
* Serial version : 
ifort -fpp -assume byterecl -mkl -o vaspberry vaspberry.f
* Multicore version : 
mpif90 -DMPI_USE -mkl -fpp -assume byterecl -o vaspberry vaspberry.f

# Features
* Berry curvature calculation
* Compute Chern number for certain band(s)
* Circular dichroism (optical selectivity response to the circulary polarized light)
* Wavefunction plot (Gamma point only in the current version)

# Example
* 1H-MoS2 : Berry curvature and Chern number
* Quantum Anomalous Hall effect (Trypheny-lead lattice) : See H.-J. Kim, C. Li, J. Feng, J.-H. Cho, and Z. Zhang, PRB 93, 041404(R) (2016) (the example files will be provided upon request)
