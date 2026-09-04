# VASPBERRY
Berry curvature, Chern number, and two-dimensional Z2 calculations from VASP
`WAVECAR` wavefunctions. VASPBERRY implements the discrete Brillouin-zone
method of [Fukui, Hatsugai, and Suzuki, J. Phys. Soc. Jpn. 74, 1674
(2005)](https://doi.org/10.1143/JPSJ.74.1674), together with circular
dichroism and real-space wavefunction output.

# Download Git version
* git clone --branch master  https://github.com/Infant83/VASPBERRY.git

# Compile

GNU serial and OpenMPI builds are available through the `Makefile`:

```bash
make serial
make mpi
```

These targets write `build/vaspberry-gfortran` and `build/vaspberry-mpi`.
Commands below use the serial product; substitute the path if you build an
executable under another name.

Intel `ifx`/`mpiifx` and legacy `ifort` commands, the direct-access record
length requirement, and BLAS/LAPACK ABI constraints are documented in the
[`build guide`](docs/BUILD.md).

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
* Two-dimensional Z2 invariant by the Fukui-Hatsugai lattice n-field method
* Circular dichroism (optical selectivity response to the circulary polarized light)
* Wavefunction plot (Gamma point only in the current version)
* Guarded WAVECAR-direct Fukui export and valley-transport analysis (Python)


## Fukui-Hatsugai Z2 invariant

`-z2 1` is a first-class two-dimensional Z2 option. It requires a
nonmagnetic, time-reversal-symmetric insulating spinor `WAVECAR` on a full,
unshifted, even Gamma-centered `Nx x Ny x 1` mesh generated with `ISYM=-1`.
For the bundled Bi 12 x 12 example:

```bash
./build/vaspberry-gfortran \
  -f examples/Bi_Z2/WAVECAR -o NFIELD -z2 1 \
  -kx 12 -ky 12 -s 2 -ii 1 -if 10
```

A PASS writes `NFIELD.dat` and `Z2_FIELD.csv`. A field result is reportable
only when `result_status=PASS`, `reportable_invariant=1`, and the top/bottom
half-zone parities agree; their common value is Z2. Rejected or interrupted Z2
runs retain `Z2_FIELD.invalid.csv` when the guarded output preflight has
started, and any `NFIELD.dat` without the final PASS CSV is non-reportable.
The input gap, physical time-reversal symmetry, and convergence with mesh
density must still be checked for the system under study.

See the [`Fukui-Hatsugai Z2 guide`](docs/Z2_FUKUI_HATSUGAI.md) for the method,
input preparation, output schema, plotting, scope, and references. Run
`./build/vaspberry-gfortran -h` for the built-in option summary.

An importable Python-library port is deferred to a later release; its proposed
validation boundary is recorded in the [`roadmap`](docs/ROADMAP.md).

# Usage
* Instruction and possible options
> ./build/vaspberry-gfortran -h
* Berry curvature calculation and Chern number (ex, k-grid: 12x12, multiband berry curvature from 1-th to 18-th band)
> ./build/vaspberry-gfortran -kx 12 -ky 12 -ii 1 -if 18
* Z2 invariant (12x12 full Gamma-centered SOC mesh, occupied bands 1--10)
> ./build/vaspberry-gfortran -f WAVECAR -o NFIELD -z2 1 -kx 12 -ky 12 -s 2 -ii 1 -if 10
* Circular dichroism [ex, transition rate from 11-th to 12-th state by right(+) polarized light]
> ./build/vaspberry-gfortran -kx 12 -ky 12 -cd 1 -ii 11 -if 12
* Real space wavefunction plot [ex, to plot 18-th state with 1-st k-point (if it is gamma point), with 40x40x40 grid for density file]
> ./build/vaspberry-gfortran -wf 18 -k 1 -ng 40,40,40
> NOTE: current version only support gamma point for wavefunction plot. (there is some problem in boundary region in off-gamma k-point)
* If your system is semimetallic, there can be following error messages: "error. !!! ne(k) /= ne(k') !!!". This is due to that the number of occupied states for certain k-point (ne(k)) counted based on the calculated Fermi level is differ over the Brillouin zone. In this case, one can explicitly specify the number of electrons (NE) of your system, so that VASPBERRY calculate berry curvature with "NE" bands. 
> ./build/vaspberry-gfortran -kx 12 -ky 12 -ii 1 -if 18 -ne 18

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

# Examples
* [1H-MoS2](examples/1H-MoS2/): Berry curvature, Chern, and Kubo plot data
* [Bi buckled honeycomb layer](examples/Bi_Z2/): 12 x 12
  Fukui-Hatsugai n-field Z2 example, input templates, reference result, and
  legacy regression evidence
* See the [examples index](examples/README.md) for the sampling requirements
  and the distinction between full-BZ and line-mode data.
* Quantum Anomalous Hall effect (Trypheny-lead lattice) : See H.-J. Kim, C. Li, J. Feng, J.-H. Cho, and Z. Zhang, PRB 93, 041404(R) (2016) (the example files will be provided upon request)
* Circular dichroism : See S.-W. Kim, H.-J. Kim, S. Cheon, and T.-H. Kim, Phys. Rev. Lett. accepted (2021) (the example will be provided upon reasonable request).

# Contributors
* Hyun-Jung Kim: Main developer
* Sun-Woo Kim: Circular dichroism and Kubo formula

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
