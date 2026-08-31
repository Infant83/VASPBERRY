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

The `-z2 1` path reconstructs spinor time reversal in the actual
plane-wave reciprocal-G basis. For a time-reversed image,

\[
\mathbf G_t=-\mathbf G_s+
\operatorname{round}(-\mathbf k_s-\mathbf k_t).
\]

At a represented TRIM \(\Lambda\), this becomes
\(\mathbf G_t=-\mathbf G_s-2\Lambda\). The reciprocal-lattice shift is
therefore required at nonzero two-dimensional TRIMs (the M points of
hexagonal MoS₂); for TRIM partners folded onto the same stored representative, the former
shift-free \(-G\) rule was valid only at Gamma.
This is an independently confirmed code-level defect. Quantitatively
attributing a previously observed MoS₂ field asymmetry to it still requires
rerunning this corrected path on the corresponding full-mesh WAVECAR.

Each link is now checked by its minimum singular value rather than by
`abs(det S)`, which depends exponentially on band rank. Determinant phases
are accumulated from a separate LU factorization without multiplying the
determinant, and the single-precision WAVECAR coefficient norm is accumulated
after promoting its real and imaginary components to double precision.

A run that passes the legacy reconstruction checks atomically publishes
`Z2_FIELD.csv`; a completed rejection writes
`Z2_FIELD.invalid.csv`. After path/ownership preflight succeeds, stale
products are removed; failures before finalization leave an `INCOMPLETE`
sentinel. The legacy NFIELD is first completed as a temporary file, then the
sentinel is removed, and `Z2_FIELD.csv` is published last as the PASS commit
marker. If that final rename itself fails, neither a regular PASS CSV nor the
sentinel is present: the staged `Z2_FIELD.tmp` remains and the program exits
nonzero.
Before cleanup, reserved basenames and POSIX `realpath()` identities reject
direct, relative, absolute, and symbolic-link aliases of the input WAVECAR.
Existing files are deleted only after their Z2 schema or legacy NFIELD markers
are recognized; otherwise preflight stops without modifying any file. In that
case, any older output is not a result of the rejected run. This also keeps a
hard-linked binary input from being mistaken for a deletable Z2 product. The
Z2 `-o` base is limited to 71 characters so the checked `.dat` and `.tmp`
paths cannot truncate. Concurrent Z2 runs in one working directory are not
supported. The CSV
contains the wrapped Berry flux, curvature, minimum link singular value,
Fukui integer n-field, thresholds, and per-plaquette diagnostics on the
fundamental mesh. The physical flux check is

\[
\operatorname{wrap}[\Phi(-\mathbf k)+\Phi(\mathbf k)]=0
\quad(\mathrm{mod}\ 2\pi).
\]

The pointwise n-field remains gauge- and logarithm-branch-dependent; only its
integer consistency and half-zone parity are used.

These checks test the numerical self-consistency of the time-reversal gauge
constructed by the legacy routine. Because the B-minus states and the even
TRIM partners are generated with the time-reversal operator, they do not
independently establish time-reversal symmetry of the input WAVECAR. Raw
\(E(\mathbf k)-E(-\mathbf k)\), TRIM Kramers splitting, occupied-projector
residuals, the occupied-unoccupied gap, and mesh convergence must be checked
separately. The guarded `tools/wavecar_z2.py` workflow performs those tests
and reports an invariant only when every guard passes.

The direct reader uses pseudo-wavefunction coefficients from `WAVECAR`;
the CSV records `WAVECAR_PSEUDO_NO_PAW_AUGMENTATION`. Missing PAW
augmentation can change the complex finite-neighbor overlap matrices,
including their phases and conditioning, but it does not generate the omitted
reciprocal-G shift. Assuming correct spinor time reversal, reciprocal-G
folding, link orientation, and a consistently TR-paired mesh, omitting a
TR-covariant augmentation term is not by itself an exact TR-covariance
breaker; it can still amplify coarse-mesh or branch-cut failures. For a PAW-aware cross-check, use
VASP-generated Bloch overlaps such as `wannier90.mmn`, with a VASP version
and symmetry setup appropriate for noncollinear spinors. See
[Ferretti et al.](https://doi.org/10.1088/0953-8984/19/3/036215) and the
[VASP PAW formalism](https://vasp.at/wiki/Projector-augmented-wave_formalism).

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
