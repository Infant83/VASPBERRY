# Valley-resolved intrinsic Hall transport

This document specifies the physically validated post-processing path for
chemical-potential-dependent Hall transport from VASPBERRY Fukui data. It also
separates quantities that the current code can calculate from quantities that
require additional energy information.

## 1. Observable and sign convention

For a two-dimensional system, the intrinsic Berry-curvature contribution is

\[
C_{\mathrm{occ}}(\mu,T)=
\sum_n\int_{\mathrm{BZ}}\frac{d^2k}{2\pi}
f(E_{n\mathbf{k}}-\mu,T)\,\Omega_{n,z}(\mathbf{k}),
\]

\[
\frac{\sigma_{xy}(\mu,T)}{e^2/h}=-C_{\mathrm{occ}}(\mu,T).
\]

The postprocessor writes both quantities. This makes the electron Hall-sign
convention explicit instead of hiding it in a plot.

The reported unit is the 2D sheet value \(e^2/h\). Conversion to S/m or S/cm
requires a separately declared effective layer thickness.

This is an oriented sheet response on the plane spanned by the ordered pair
\((\mathbf b_1,\mathbf b_2)\). The symbol \(\sigma_{xy}\) is literal when those
vectors span the Cartesian \(xy\) slab plane and
\(\mathbf b_1\times\mathbf b_2\) points along \(+z\). For a nonstandard plane
or a left-handed reciprocal basis, the reported component and sign must instead
be interpreted relative to that declared orientation.

## 2. What VASPBERRY's Fukui output represents

For each uniform-mesh cell, VASPBERRY follows the loop

\[
\mathbf{k},\quad
\mathbf{k}+\Delta\mathbf{k}_1,\quad
\mathbf{k}+\Delta\mathbf{k}_1+\Delta\mathbf{k}_2,\quad
\mathbf{k}+\Delta\mathbf{k}_2
\]

and calculates the gauge-invariant plaquette flux

\[
\phi_{n,p}=-\operatorname{Im}\log
\prod_{i=1}^{4}\langle u_{n\mathbf{k}_i}|u_{n\mathbf{k}_{i+1}}\rangle.
\]

The legacy file prints

\[
\Omega_{n,p}=\phi_{n,p}/\Delta S_k,
\qquad
\Delta S_k=\frac{|\mathbf b_1\times\mathbf b_2|}{N_xN_y}.
\]

Despite an old header that labels the curvature as `A^-2`, the physical unit
of \(\Omega\) is \(\text{Angstrom}^2\); \(\Delta S_k\) is in
\(\text{Angstrom}^{-2}\).

**The fractional coordinates printed in `BERRYCURV.dat` are the plaquette
centers**, not the lower-left vertices. The four energy vertices are therefore

\[
\mathbf q_p \pm \frac{\Delta\mathbf k_1}{2}
             \pm \frac{\Delta\mathbf k_2}{2}.
\]

## 3. Correct finite-mesh occupation weighting

For a Fukui plaquette, the implemented quadrature is

\[
\bar f_{n,p}(\mu,T)=\frac14\sum_{i=1}^{4}
 f(E_{n\mathbf{k}_{p,i}}-\mu,T),
\]

\[
C_{\mathrm{occ}}(\mu,T)
\simeq\sum_{n,p}\frac{\phi_{n,p}}{2\pi}\bar f_{n,p}(\mu,T).
\]

The code intentionally does **not** use \(f(\bar E_p)\), because in general

\[
\frac14\sum_i f(E_i)\ne f\!\left(\frac14\sum_iE_i\right).
\]

At zero temperature, a cell crossed by the Fermi surface can consequently have
an occupation fraction of 0, 1/4, 1/2, 3/4, or 1. This is a convergent mesh
quadrature, not an exact reconstruction of the Fermi contour. A dense mesh and
mesh-convergence study remain mandatory near sharp valley pockets.

## 4. CBM-only contribution versus total conductivity

A single isolated CBM band gives the doping-induced active-band contribution,

\[
\Delta C_{\mathrm{CBM}}(\mu,T)=
\sum_p\frac{\phi_{\mathrm{CBM},p}}{2\pi}\bar f_{\mathrm{CBM},p}.
\]

It is not automatically the total Hall conductivity. For a chemical-potential
window in which all valence bands remain fully occupied, their combined Fukui
Chern number can be supplied as a constant baseline:

\[
C_{\mathrm{occ}}(\mu,T)=C_{\mathrm{core}}+
\sum_{n\in\mathrm{active}}\Delta C_n(\mu,T).
\]

Use `--core-chern` for this baseline. This avoids recomputing every deep valence
band separately. A scalar baseline is included only in the **total** response;
it cannot be assigned to K or K' without a corresponding spatially resolved
valence-manifold curvature map.

## 5. Conditions for a single-band calculation

Fermi-weighting a single-band Fukui flux is physically meaningful only if that
band is isolated and nondegenerate over the relevant Brillouin zone and energy
window. The postprocessor reports the minimum direct gap to all other EIGENVAL
bands and rejects gaps below `--isolation-tolerance`.

The default threshold is only a numerical guard. It is not a publication-level
definition of isolation. The threshold, k mesh, band tracking, temperature, and
the physical separation from neighboring bands must all be tested.

If the CBM is degenerate, crosses another band, or becomes entangled in the
target window, an occupation-weighted gauge-covariant multiband/density-matrix
or independently validated Kubo/Wannier treatment is required. A non-Abelian
manifold curvature cannot be multiplied by one band's Fermi occupation.

## 6. Valley decomposition

The K and K' contributions are defined by user-supplied valley centers and a
Cartesian radius in inverse Angstrom. Distances are minimized over neighboring
reciprocal images, which is necessary for oblique and hexagonal cells. The
nearest-image search uses a singular-value bound and is not limited to the
surrounding 3x3 copies; a singular or pathologically ill-conditioned in-plane
basis is rejected with a request to use a reduced reciprocal basis.

The tool rejects overlapping K/K' disks. It reports

- total intrinsic Hall conductivity;
- K contribution;
- K' contribution;
- contribution outside both disks;
- the diagnostic contrast K minus K'.

The contrast is a declared numerical diagnostic, not a universally normalized
definition of valley conductivity. The valley radius must be converged and tied
to the actual band/curvature distribution.

## 7. Required inputs and validation gates

The legacy path needs:

1. one **single-band Fukui** `BERRYCURV...dat` per active band;
2. the matching full-uniform-BZ `EIGENVAL` on the same k mesh;
3. reciprocal vectors and mesh metadata in the curvature header.

Before integration, `tools/vaspberry_transport.py` checks:

- periodic 3x3 legacy copies are grouped modulo reciprocal lattice vectors;
- all copies have the same Berry curvature before being collapsed;
- the unique grid agrees with header `NKPOINT` and `K-GRID`;
- the mesh is a complete uniform \(N_x\times N_y\) product on one \(k_z\) plane;
- both mesh dimensions contain at least two points, and the curvature and
  EIGENVAL \(k_z\) planes agree modulo a reciprocal vector;
- EIGENVAL has the same full mesh and uniform weights summing to one;
- all four EIGENVAL vertices exist for every printed Fukui center;
- requested band numbers, total band counts, and labeled UP/DN channels agree
  with the curvature headers, filenames, and selected EIGENVAL spin column;
- an `ISPIN=2` curvature file used with EIGENVAL retains an unambiguous `.UP` or
  `.DN` filename (or explicit `SPIN` header); a renamed unlabeled channel is
  rejected rather than inferred from global `--spin`;
- header \(\Delta S_k\) agrees with
  \(|\mathbf b_1\times\mathbf b_2|/(N_xN_y)\);
- the unique-grid full-occupation Chern sum agrees with the file header;
- the full-occupation single-band Fukui sum is within 0.005 of its nearest
  integer (the legacy rounded \(C\simeq-4.0\times10^{-8}\) example passes);
- the maximum plaquette flux stays below a safety fraction of the principal-log
  branch boundary; otherwise a denser mesh is required;
- every requested active band passes the numerical isolation check.

Nonfinite coordinates, curvature, energies, occupations, k weights, chemical
potentials, temperature, and numerical thresholds are rejected before plotting
or integration. The `map` and `cut` commands require a two-dimensional full-BZ
mesh, while `line` requires an ordered path; these sampling topologies are not
silently reinterpreted.

At the start of each `sigma` attempt, the tool removes only the three explicitly
requested `--csv`, `--plot`, and optional `--summary` targets. Thus a rejected
rerun cannot leave an older artifact that looks like the failed run's result;
no parent directory or wildcard target is removed.
Before that cleanup, all target paths are resolved without requiring them to
exist. The targets must be mutually distinct and must not alias any `--input`
or `--eigenval` path; a collision is rejected before any file is unlinked.

The supplied external archive currently contains `WAVECAR` but no matching
`EIGENVAL`. Therefore a real \(\sigma_{xy}(\mu,T)\) must use either energies
exported directly from that WAVECAR or an EIGENVAL generated on the identical
full mesh. A WAVECAR-only Berry-curvature map must not be presented as an
energy-resolved conductivity result.

## 8. Direct WAVECAR export and cumulative-subspace transport

The postprocessors require Python 3.10 or later and the packages listed in
`requirements-transport.txt`.

`tools/wavecar_fukui.py` is the preferred path when the matching full-mesh
`EIGENVAL` is unavailable. It reads energies and spinor coefficients directly
from the same WAVECAR, supports both byte RECL and legacy four-byte-word RECL,
and reproduces VASPBERRY's plane-wave ordering and boundary gauge shifts.
The current direct reader supports standard single-precision complex
coefficients (`RTAG=45200`) and rejects other WAVECAR coefficient encodings.

For a target CBM band 19, the default export calculates five maps:

- `V`: cumulative valence manifold 1:18;
- `VC`: cumulative manifold 1:19;
- `VCC`: cumulative manifold 1:20;
- `C`: single band 19 (diagnostic only unless isolated);
- `Cpair`: conduction pair 19:20 (symmetry and conditioning diagnostic).

```bash
python tools/wavecar_fukui.py WAVECAR --nx 12 --ny 12 \
  --energy-band 19 --output-dir direct_raw \
  --valley-k 0.6666667,0.3333333,0 \
  --valley-kp 0.3333333,0.6666667,0 --plot
```

Each `fukui_*.csv` row is one original-BZ plaquette. The stable interchange
columns are:

- `cell_id`, `ix`, `iy`, `band`, `band_first`, `band_last`;
- `k0` through `k3` (one-based WAVECAR vertex indices);
- `e0_eV` through `e3_eV` for band 19 and `energy_above1_ev` through
  `energy_above4_ev` for band 20;
- the four target-band gaps above and below, plus the cumulative subspace's
  `subspace_min_gap_above_ev`;
- `q1`, `q2`, `q3`, `phi_rad`, `omega_A2`, and `dS_A-2`;
- four link determinant magnitudes, four minimum singular values,
  `link_quality_min_sv`, and plane-wave intersection coverage.

`diagnostics.json` records lattice/reciprocal vectors, mesh, RECL convention,
energy gaps, full-occupation Chern sums, phase range, and worst conditioning.
The direct reader rejects nonpositive plane-wave counts and nonfinite k points,
energies, occupations, or wave-function coefficients before they can enter an
overlap. JSON diagnostics are emitted in strict form; a diagnostic nonfinite
value is represented as `null`, never the nonstandard `NaN` or `Infinity`.
The mesh detector accepts a uniform shifted Monkhorst mesh and records its
fractional offset. It rejects line-mode, incomplete, duplicate, or nonuniform
point sets, nonfinite coordinates, and either in-plane dimension below two.
Map plots preserve that offset in both their image extent and periodic K/K'
markers. `--min-pw-coverage` (default 0.90) is a hard guard on the common
plane-wave fraction of every neighboring overlap.

The default map ranges follow `--energy-band N` dynamically: 1:`N-1`, 1:`N`,
1:`N+1`, `N`, and `N`:`N+1`. Raw export and plotting require band `N+1` for
the reported upper energy/gap and pair maps. `--transport-t0` additionally
requires band `N+2` for the upper-manifold isolation and energy-window check.
For a collinear
`ISPIN=2` WAVECAR, use `--spinor-components 1 --spin 1` or `--spin 2`; an SOC
`ISPIN=1` WAVECAR uses `--spinor-components 2 --spin 1`.

At zero temperature, the guarded cumulative-subspace scan assigns each of the
four vertices to 1:18, 1:19, or 1:20 according to its band-19 and band-20
energies, then applies the Sawahata-style average

\[
\bar\phi_p(\mu)=\frac{1}{4}\sum_{i=1}^4
\phi_p^{(N_i(\mu))}.
\]

This follows the local Berry-phase strategy tested against Wannier/Kubo results
by Sawahata *et al.*, Phys. Rev. B **107**, 024404 (2023),
https://doi.org/10.1103/PhysRevB.107.024404.

```bash
python tools/wavecar_fukui.py WAVECAR --nx 12 --ny 12 \
  --energy-band 19 --output-dir direct_transport \
  --transport-t0 --mu-min 0.40 --mu-max 0.55 --mu-num 151 \
  --valley-k 0.6666667,0.3333333,0 \
  --valley-kp 0.3333333,0.6666667,0 --plot
```

The K coordinates above are only examples. With no radius, every cell is
assigned to the nearest periodic K or K' center. `--valley-radius R` instead
creates disjoint K/K'/outside masks in inverse Angstrom.
K and K' must be finite and periodically distinct. Exact Voronoi ties are put
in the outside mask rather than assigned by floating-point accident; their
count and the exclusive/exhaustive mask totals are written to
`transport_t0_diagnostics.json`. A radius partition is rejected if its two
masks overlap or if either valley contains no plaquette center.
Periodic distances and the K-to-K' line image use an exact two-dimensional
closest-lattice-vector search after Gauss reduction, so highly skew or
non-reduced reciprocal bases are not limited to an unsafe fixed image stencil.

The transport command first treats the fully occupied 1:`N-1` valence manifold
as a hard global baseline. Its minimum link singular value, minimum direct gap
to band `N`, and maximum absolute plaquette phase must pass the same declared
quality thresholds at every cell. This baseline gate cannot be bypassed by
`--allow-invalid-transport`.

Because that baseline is held fully occupied, `--mu-min` must be strictly above
the global maximum of band `N-1`; otherwise the requested window contains holes
that the three-manifold formula does not represent. Similarly, `--mu-max` must
stay below the global minimum of band `N+2`.

After the baseline passes, the command checks active-cell link singular values,
adjacent-band gaps, and the principal-branch phase margin. If a plaquette uses
the 1:`N` (`V+1`) flux at even one vertex, the `E(N+1)-E(N)` isolation gap is
checked at **all four plaquette vertices**. Likewise, use of the 1:`N+1`
(`V+2`) flux triggers an all-four-vertex `E(N+2)-E(N+1)` gap check. This avoids
accepting a mixed-occupation plaquette whose unused corner is degenerate and
therefore makes the cumulative-subspace flux ill-conditioned. If any active
cell fails, the command writes `transport_t0_diagnostics.json` but refuses
`transport_t0.csv`. The
`--allow-invalid-transport` switch exists only to inspect an explicitly labeled
`INVALID` diagnostic scan; such a file is not a validated physical result.
Before a transport attempt, the tool removes only its three prior transport
artifacts (`transport_t0.csv`, `transport_t0_diagnostics.json`, and
`wavecar_fukui_sigma_mu.png`) so a rejected rerun cannot leave an older scan or
plot masquerading as the new result. All chemical-potential, valley, radius,
and quality-threshold inputs are checked for finite values before use.
Before deleting or writing anything, it also compares the resolved WAVECAR path
with every output planned for that invocation (raw CSVs, diagnostics, plots,
and transport files) and refuses any collision, preventing the input from being
unlinked or overwritten.

With `--plot`, the direct tool writes:

- `wavecar_fukui_kresolved.png`: energy, cumulative curvature, link quality,
  and neighbor-gap maps with K/K' marked. Every diverging curvature panel uses
  symmetric limits about zero so white always means zero;
- `wavecar_fukui_line_K_Kp.png`: the periodic shortest K-to-K' energy and
  curvature line cut;
- `wavecar_fukui_sigma_mu.png`: the T=0 net/K/K'/contrast scan. Rejected
  chemical potentials are drawn with dotted low-alpha curves and red markers.

## 9. Command examples

Plot a legacy curvature map:

```bash
python tools/vaspberry_transport.py map \
  --input 1H-MoS2/BERRYCURV.dat \
  --output berry_map.png
```

Plot a Kubo line result together with its band energy:

```bash
python tools/vaspberry_transport.py line \
  --input 1H-MoS2/KPATH/3.BC_kubo/BERRYCURV_KUBO.EIG-10.dat \
  --eigenval 1H-MoS2/KPATH/2.band/EIGENVAL \
  --band 10 --output kubo_path.png
```

Interpolate a periodic shortest-image K-to-K' cut from a full-BZ map:

```bash
python tools/vaspberry_transport.py cut \
  --input BERRYCURV_BAND-19.dat --eigenval EIGENVAL --band 19 \
  --k-center 0.333333 0.666667 \
  --kp-center 0.666667 0.333333 \
  --output valley_cut.png
```

Calculate the active CBM contribution and add a fully occupied valence baseline:

```bash
python tools/vaspberry_transport.py sigma \
  --input BERRYCURV_BAND-19.dat --bands 19 \
  --eigenval EIGENVAL \
  --mu-min 0.0 --mu-max 0.5 --mu-points 501 --temperature 10 \
  --core-chern 0.0 \
  --k-center 0.333333 0.666667 \
  --kp-center 0.666667 0.333333 \
  --valley-radius 0.20 \
  --csv sigma_xy.csv --plot sigma_xy.png --summary sigma_xy.json
```

The K coordinates above are examples only. Use the convention and reciprocal
basis of the actual calculation.

## 10. Kubo status

Legacy `-kubo 2` line-mode output is accepted for visualization only. It is not
accepted by the validated transport integrator because a one-dimensional path
does not define a Brillouin-zone integral. In addition, the checked-in band-10
example has a zero minimum direct gap and an approximately
\(6.27\times10^4\ \text{Angstrom}^2\) spike at Gamma. This exposes the present
machine-epsilon degeneracy handling and makes that example unsuitable as a
quantitative transport reference.

Before any full-grid legacy Kubo result is used quantitatively, its degeneracy
regularization, normalization, PAW/nonlocal/SOC velocity terms, unoccupied-band
convergence, and agreement with Fukui/Wannier benchmarks must be established.

## 11. Test coverage

Run:

```bash
PYTHONPATH=tools python -m unittest discover -s tests -v
```

The tests cover replica invariance and disagreement detection, exact Chern/sign
normalization, four-vertex occupation averaging at zero and finite temperature,
full and shifted/shuffled EIGENVAL meshes, nonuniform and duplicate mesh
rejection, K/K'/outside partition closure, core-Chern baselines, manifold/Kubo
rejection, and the checked-in MoS2 9-copy-to-unique-grid Chern reconstruction.

An independent Qi-Wu-Zhang two-band oracle additionally fixes the orientation
and sign convention: for right-handed \(\mathbf b_1\times\mathbf b_2=+\hat z\)
and VASPBERRY's \(-\operatorname{Im}\log\) plaquette phase, the occupied
lower band at \(m=+1\) has \(C=-1\), and therefore
\(\sigma_{xy}=+e^2/h\). An independent analytic point-curvature oracle verifies
monotonic Riemann-sum convergence from \(8\times8\) through \(32\times32\)
without weakening the nearest-integer gate applied to actual Fukui plaquette
fluxes. Separate synthetic
K/K' packets verify exact cancellation for degenerate, opposite-curvature
valleys and the expected finite total response when a controlled valley-energy
splitting places the chemical potential between the two minima.

## 12. Remaining physical limitations

This is an intrinsic, rigid-band Berry-curvature calculation. It does not include
disorder, side jump, skew scattering, phonons, carrier self-consistency, or
intervalley relaxation. Electron/hole density should be calculated alongside
\(\sigma_{xy}(\mu,T)\) before connecting a chemical-potential window to an
experiment. For scalar non-SOC calculations, whether a factor-two spin
degeneracy is required must be determined from how the WAVECAR bands were stored;
SOC and explicitly collinear spin channels must not be doubled automatically.
