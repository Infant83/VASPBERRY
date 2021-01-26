# VASPBERRY
Berry curvature and Chern number calculations with the output (WAVECAR) of VASP code.
VASPBERRY is written for the post-processing purpose of the VASP outputs, i.e., WAVECAR the Bloch wavefunction information. VASPBERRY can compute Berry curvature and Chern number via Fukui's method [See J. Phys. Soc. Jap. 74, 1674 (2005)]. In addition Circular dichroism also can be evaluated. Since it directly reads the wavefunction coefficients, one can also obtain real space wavefunction character psi_nk(r) by simple command.

# Compile
* Serial version : 
    > ifort -fpp -assume byterecl -mkl -o vaspberry vaspberry.f
* Multicore version : 
    > mpif90 -DMPI_USE -mkl -fpp -assume byterecl -o vaspberry vaspberry.f

* Note for gfortran:
    For gfortran, please use vaspberry_gfortran_serial.f for the compilation. This only support non-parallel calculations.
    For the compilation, for example
    > gfortran -L/usr/local/lib/lapack/ -l lapack -o vaspberry vaspberry_gfortran_serial.f

# Features
* Berry curvature calculation
* Compute Chern number for certain band(s)
* Circular dichroism (optical selectivity response to the circulary polarized light)
* Wavefunction plot (Gamma point only in the current version)

# Usage
* Instruction and possible options
> ./vaspberry -h
* Berry curvature calculation and Chern number (ex, k-grid: 12x12, multiband berry curvature from 1-th to 18-th band)
> ./vaspberry -kx 12 -ky 12 -ii 1 -if 18
* Circular dichroism [ex, transition rate from 11-th to 12-th state by right(+) polarized light]
> ./vaspberry -kx 12 -ky 12 -cd -ii 11 -if 12
* Real space wavefunction plot [ex, to plot 18-th state with 1-st k-point (if it is gamma point), with 40x40x40 grid for density file]
> ./vaspberry -wf 18 -k 1 -ng 40,40,40
* If your system is semimetallic, there can be following error messages: "error. !!! ne(k) /= ne(k') !!!". This is due to that the number of occupied states for certain k-point (ne(k)) counted based on the calculated Fermi level is differ over the Brillouin zone. In this case, one can explicitly specify the number of electrons (NE) of your system, so that VASPBERRY calculate berry curvature with "NE" bands. 
> ./vaspberry -kx 12 -ky 12 -ii 1 -if 18 -ne 18

# Example
* 1H-MoS2 : Berry curvature and Chern number
* Quantum Anomalous Hall effect (Trypheny-lead lattice) : See H.-J. Kim, C. Li, J. Feng, J.-H. Cho, and Z. Zhang, PRB 93, 041404(R) (2016) (the example files will be provided upon request)
