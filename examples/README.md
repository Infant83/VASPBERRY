# Examples

The examples are grouped by physical task. Run commands from the repository
root unless an example says otherwise.

| Example | Sampling | Purpose | Main input |
|---|---|---|---|
| [`Bi_Z2/`](Bi_Z2/) | full even Gamma-centered 2D mesh | Fukui-Hatsugai n-field Z2 invariant | SOC spinor `WAVECAR` |
| [`1H-MoS2/`](1H-MoS2/) | full 12 x 12 mesh | Berry curvature and Chern regression | stored `BERRYCURV.dat` |
| [`1H-MoS2/KPATH/`](1H-MoS2/KPATH/) | K-Gamma-K' line | band and Kubo plot data | line-mode `WAVECAR`/`EIGENVAL` |

The line-mode MoS2 data are not a full Brillouin-zone mesh and must not be
used for a Chern or Z2 calculation.

VASP `POTCAR` files are not distributed. Reconstruct them locally from a
licensed VASP potential library in the element order given by `POSCAR`; see
the pseudopotential provenance file in each example.
