# VASPBERRY
Berry curvature and Chern number calculations with the output (WAVECAR) of VASP code.

# Compile
* Serial version : 
ifort -fpp -assume byterecl -mkl -o vaspberry vaspberry.f
* Multicore version : 
mpif90 -DMPI_USE -mkl -fpp -assume byterecl -o vaspberry vaspberry.f

# Features
* Berry curvature calculation
* Compute Chern number for certain band(s)
* Circular dichroism (optical selectivity due to the circulary polarized right)
* Wavefunction plot (Gamma point only in the current version)

# Example
* 1H-MoS2 : Berry curvature and Chern number
* Quantum Anomalous Hall effect (Trypheny-lead lattice) : See H.-J. Kim, C. Li, J. Feng, J.-H. Cho, and Z. Zhang, PRB 93, 041404(R) (2016) (the example files will be provided upon request)
