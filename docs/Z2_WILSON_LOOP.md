# Wilson-loop Z2 validation

`tools/wavecar_z2.py` evaluates the occupied spinor subspace of a complete,
Gamma-centered two-dimensional WAVECAR. It uses SVD polar-unitary links and
hybrid Wannier charge centers (WCCs). The tool reports an integer `z2` only
when all gap, time-reversal, link-quality, topology, and fixed-grid checks pass.

The existing Fortran `-z2` routine and its `NFIELD.dat` output are retained for
compatibility. `NFIELD.dat` is gauge-dependent and its parity is a **legacy Z2
candidate**, not a validated invariant. A successful run of this Python tool is
the reportable path in version 1.1.0.

## Run

Install the same Python dependencies used by the transport tools, then run:

```bash
python tools/wavecar_z2.py WAVECAR \
  --nx 12 --ny 12 --occupied-bands 18 --axis both \
  --output-dir z2_result --plot
```

Requirements:

- `ISPIN=1` with a two-component SOC/noncollinear spinor WAVECAR;
- a full, uniform, Gamma-centered 2D mesh on `kz=0` modulo a reciprocal lattice
  vector, with even `nx` and `ny`, each at least four;
- an even occupied rank and the next band as an unoccupied gap sentinel; and
- both Wilson-loop directions for a reportable result. `--axis x` or
  `--axis y` can be used for diagnosis, but it deliberately returns INVALID.

The default outputs are:

| File | Contents |
|---|---|
| `z2_wilson_wcc.csv` | Both loop directions, pump coordinate, WCC index, and WCC position in cycles |
| `z2_diagnostics.json` | Thresholds, measured diagnostics, individual checks, `candidate_z2`, and reportable `z2` |
| `z2_wilson_wcc.png` | Optional (`--plot`) full-pump TR comparison and half-BZ largest-gap diagnostic |

The program protects the input WAVECAR from output-name collisions and writes
CSV/JSON/PNG through temporary files. A hard failure removes all planned
outputs rather than leaving a partial result.

| Exit status | Meaning |
|---:|---|
| `0` | Every strict check passed; the integer `z2` is reportable for this mesh |
| `2` | Calculation completed but a check failed; `z2` is `null` and `candidate_z2` is diagnostic only |
| `1` | Input, I/O, or computation failed; planned outputs are removed |

Do not substitute `candidate_z2` for a failed `z2`. In particular, a stable
parity on one coarse mesh does not establish convergence.

## Calculation and checks

For neighboring occupied frames, the overlap matrix is replaced by its polar
factor,

\[
M=U\Sigma V^\dagger,\qquad \widetilde M=UV^\dagger.
\]

Multiplication around a closed reciprocal-space string gives the Wilson loop.
Its eigenphases divided by \(2\pi\) are the WCC set. The code evaluates the
largest-gap crossing parity over the closed half Brillouin zone following the
set-based Wilson-loop construction. It never connects sorted WCC indices into
putative physical branches: Kramers partners may exchange their sorted labels.

The following checks are independent hard gates:

- positive direct and indirect gaps between bands `N` and `N+1`;
- band-energy and occupied-projector agreement under \(\mathbf{k}\to-\mathbf{k}\);
- complete time-reversed G-basis matching, retained raw wavefunction norm, and
  Kramers energy pairing at the four TR-invariant momenta;
- zero occupied-subspace Chern number and odd Berry flux under time reversal;
- nonsingular overlap matrices, plane-wave coverage, and polar unitarity;
- equality of the **unordered** WCC sets at pump coordinates \(q\) and \(-q\),
  plus Kramers-paired WCCs at \(q=0\) and \(q=1/2\);
- agreement of candidate parities from the two loop directions; and
- fixed-grid position (POS), neighboring-gap (GAP), and WCC-motion (MOVE)
  controls.

The GAP and MOVE formulas follow the neighboring-line controls used by Z2Pack,
including its union-largest-gap WCC matcher. POS compares the full loop with a
loop sampled at every second k point. This program does not adaptively insert
new pump lines, so these are **fixed-grid convergence guards**, not the same as
an adaptive Z2Pack convergence run and not statistical error estimates. A
denser full-BZ WAVECAR remains necessary when any of them fails.

Small numerical TR residuals show that the supplied wavefunctions are
consistent with the implemented time-reversal transformation. They do not, by
themselves, establish that the physical Hamiltonian has time-reversal symmetry;
that also requires the magnetic structure, VASP inputs, and symmetry setting.

## Reading the figure

The horizontal coordinate in a WCC plot is the Wilson-loop pump parameter, not
a Cartesian reciprocal-space map. Therefore this plot should not be folded
into a hexagonal first Brillouin zone. Hexagonal plotting is appropriate for a
k-resolved Berry-curvature map, while the Z2 figure must display the WCC sets at
`q` and `-q` and their Kramers pairing at the two TR pump lines.

Blue points are unordered WCC eigenvalues. Red crosses mark the largest-gap
centers used by the parity calculation. An invalid plot has a neutral
`INVALID (candidate Z2 = ...)` title and lists the failed axis-specific guards;
it does not visually promote the candidate to an invariant.

## References

- T. Fukui and Y. Hatsugai, “Quantum Spin Hall Effect in Three Dimensional
  Materials: Lattice Computation of Z2 Topological Invariants and Its
  Application to Bi and Sb,” *J. Phys. Soc. Jpn.* **76**, 053702 (2007).
  [doi:10.1143/JPSJ.76.053702](https://doi.org/10.1143/JPSJ.76.053702)
- R. Yu, X. L. Qi, A. Bernevig, Z. Fang, and X. Dai, “Equivalent expression of
  Z2 topological invariant for band insulators using the non-Abelian Berry
  connection,” *Phys. Rev. B* **84**, 075119 (2011).
  [doi:10.1103/PhysRevB.84.075119](https://doi.org/10.1103/PhysRevB.84.075119)
- D. Gresch *et al.*, “Z2Pack: Numerical implementation of hybrid Wannier
  centers for identifying topological materials,” *Phys. Rev. B* **95**,
  075146 (2017).
  [doi:10.1103/PhysRevB.95.075146](https://doi.org/10.1103/PhysRevB.95.075146)
