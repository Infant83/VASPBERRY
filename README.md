# VASPBERRY
Berry curvature and Chern number calculations with the output (WAVECAR) of VASP code.
VASPBERRY is written for the post-processing purpose of the VASP outputs, i.e., WAVECAR the Bloch wavefunction information. VASPBERRY can compute Berry curvature and Chern number via Fukui's method [See J. Phys. Soc. Jap. 74, 1674 (2005)]. In addition Circular dichroism also can be evaluated. Since it directly reads the wavefunction coefficients, one can also obtain real space wavefunction character psi_nk(r) by simple command.

# Download Git version
* git clone --branch master  https://github.com/Infant83/VASPBERRY.git

# Compile
* Serial version : 
    > ifort -fpp -assume byterecl -mkl -o vaspberry vaspberry.f
* Multicore version : 
    > mpif90 -DMPI_USE -mkl -fpp -assume byterecl -o vaspberry vaspberry.f

* Note for gfortran:
    For gfortran, please use vaspberry_gfortran_serial.f for the compilation. This only support non-parallel calculations.
    For the compilation, for example
    > gfortran -L/usr/local/lib/lapack/ -l lapack -o vaspberry vaspberry_gfortran_serial.f

> **WAVECAR compatibility:** The VASP writer and VASPBERRY reader must
> use the same direct-access `RECL` convention. Byte-based `RECL` is
> recommended; with Intel Fortran, compile both using
> `-assume byterecl`. If the first WAVECAR header is readable but
> `NKPOINT`, `NBANDS`, or `ENCUT` is zero/invalid, or the lattice
> vectors are `NaN`, suspect a 4-byte `RECL` mismatch. Rebuilding VASP
> with byte `RECL` and regenerating `WAVECAR` is the preferred fix.

# Features
* Berry curvature calculation
* Compute Chern number for certain band(s) which are well-isolated over the BZ
* Circular dichroism (optical selectivity response to the circulary polarized light)
* Wavefunction plot (Gamma point only in the current version)
* Guarded WAVECAR-direct Fukui export and valley-transport analysis (Python)


## Legacy Fortran Z2 field diagnostics (version 1.1.1)

The `-z2 1` path now reconstructs spinor time reversal in the actual
plane-wave G basis. At a represented TRIM `Lambda`, the Kramers partner uses
`G_target = -G_source - 2*Lambda`; the additional reciprocal-lattice shift
is essential at the three M points. A validated run writes `Z2_FIELD.csv`
for the fundamental mesh and keeps the older `NFIELD.dat` compatibility
output.

`Z2_FIELD.csv` separates the physical Berry flux, which must obey
`F(-k) = -F(k)`, from the Fukui integer n-field, whose pointwise pattern
depends on gauge and logarithm branch. The file is rejected if its TR partner
map, flux oddness, integer residual, phase margin, total Chern number, or
half-zone parity checks fail.

This direct reader uses pseudo-wavefunction coefficients from `WAVECAR`;
the CSV therefore records
`WAVECAR_PSEUDO_NO_PAW_AUGMENTATION`. For quantitative PAW-corrected
neighbor overlaps, use a VASP-generated Bloch-overlap route such as
`wannier90.mmn` and verify the VASP version and symmetry settings. The
Fortran parity remains a legacy candidate; the guarded
`tools/wavecar_z2.py` Wilson-loop result is the reportable route only when
all diagnostics and mesh-convergence checks pass.

# Usage
* Instruction and possible options
> ./vaspberry -h
* Berry curvature calculation and Chern number (ex, k-grid: 12x12, multiband berry curvature from 1-th to 18-th band)
> ./vaspberry -kx 12 -ky 12 -ii 1 -if 18
* Circular dichroism [ex, transition rate from 11-th to 12-th state by right(+) polarized light]
> ./vaspberry -kx 12 -ky 12 -cd 1 -ii 11 -if 12
* Real space wavefunction plot [ex, to plot 18-th state with 1-st k-point (if it is gamma point), with 40x40x40 grid for density file]
> ./vaspberry -wf 18 -k 1 -ng 40,40,40
> NOTE: current version only support gamma point for wavefunction plot. (there is some problem in boundary region in off-gamma k-point)
* If your system is semimetallic, there can be following error messages: "error. !!! ne(k) /= ne(k') !!!". This is due to that the number of occupied states for certain k-point (ne(k)) counted based on the calculated Fermi level is differ over the Brillouin zone. In this case, one can explicitly specify the number of electrons (NE) of your system, so that VASPBERRY calculate berry curvature with "NE" bands. 
> ./vaspberry -kx 12 -ky 12 -ii 1 -if 18 -ne 18

* Energy-resolved and valley-resolved Hall transport is available through the
  opt-in Python 3.10+ tools. The direct reader exports high-precision plaquette flux,
  four vertex energies, adjacent-band gaps, and link-quality diagnostics from
  the same full-mesh WAVECAR:
> python -m pip install -r requirements-transport.txt
> python tools/wavecar_fukui.py WAVECAR --nx 12 --ny 12 --energy-band 19 --output-dir direct_raw --valley-k 0.6666667,0.3333333,0 --valley-kp 0.3333333,0.6666667,0 --plot

  To request a guarded zero-temperature chemical-potential scan, add for example:
> python tools/wavecar_fukui.py WAVECAR --nx 12 --ny 12 --energy-band 19 --output-dir direct_transport --transport-t0 --mu-min 0.40 --mu-max 0.55 --mu-num 151 --valley-k 0.6666667,0.3333333,0 --valley-kp 0.3333333,0.6666667,0 --plot

  To scan every represented occupied subspace through a VBM at band 18,
  including hole occupations, use the separate full-window mode. Band 19 is
  retained as the mandatory unoccupied sentinel:
> python tools/wavecar_fukui.py WAVECAR --nx 12 --ny 12 --energy-band 18 --output-dir full_valence_scan --transport-full-t0 18 --mu-min -8.0 --mu-max 0.40 --mu-num 1001 --valley-k 0.6666667,0.3333333,0 --valley-kp 0.3333333,0.6666667,0 --plot --plot-domain first-bz

  This output is cumulative `sigma_xy(mu)`. It does not claim a smooth
  `d sigma_xy/d mu` spectral density or a selected-mu k-resolved transport map;
  on a 12x12 mesh those require separate convergence and representation work.

  See the [valley-transport guide](docs/VALLEY_TRANSPORT.md) before interpreting
  a single-band or transport result. The
  guarded transport path refuses active plaquettes with nearly singular links,
  insufficient band isolation, or branch-risk phases. The three-manifold mode
  validates its fully occupied valence baseline globally and keeps its window
  above the VBM. The full-window mode instead validates the `MAX_BAND` reference
  bundle globally and keeps `mu_max` below the `MAX_BAND+1` sentinel.
  Existing Fortran Fukui output and default calculations are unchanged; the
  opt-in spin-polarized Kubo path includes an independent accumulator fix.

# Example
* 1H-MoS2 : Berry curvature and Chern number
* Quantum Anomalous Hall effect (Trypheny-lead lattice) : See H.-J. Kim, C. Li, J. Feng, J.-H. Cho, and Z. Zhang, PRB 93, 041404(R) (2016) (the example files will be provided upon request)
* Circular dichroism : See S.-W. Kim, H.-J. Kim, S. Cheon, and T.-H. Kim, Phys. Rev. Lett. accepted (2021) (the example will be provided upon reasonable request).

# Contributors
* Hyun-Jung Kim: Main developer (Infant@kias.re.kr, currently at LG Display)
* Sun-Woo Kim: Circular dichroism, Kubo formular (kimsunwoo821@gmail.com, Department of Physics, Sungkyunkwan University)

# Citation of the code:
@software{Kim_VASPBERRY_2018,author = {Kim, Hyun-Jung},doi = {10.5281/zenodo.1402593},month = {8},title = {{VASPBERRY}},url = {https://github.com/Infant83/VASPBERRY},version = {1.0},year = {2018}}

@article{PhysRevLett.128.046401,
  title = {Circular Dichroism of Emergent Chiral Stacking Orders in Quasi-One-Dimensional Charge Density Waves},
  author = {Kim, Sun-Woo and Kim, Hyun-Jung and Cheon, Sangmo and Kim, Tae-Hwan},
  journal = {Phys. Rev. Lett.},
  volume = {128},
  issue = {4},
  pages = {046401},
  numpages = {6},
  year = {2022},
  month = {Jan},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevLett.128.046401},
  url = {https://link.aps.org/doi/10.1103/PhysRevLett.128.046401}
}
