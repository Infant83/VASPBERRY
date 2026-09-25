# Circular transition strengths from standard VASP optical output

VASPBERRY offers two complementary routes for optical valley selection:

| Input and command | Calculated quantity | VASP source changes |
|---|---|---|
| Existing `WAVECAR`, native `-cd 2` | Bare-momentum spectra and helicity selectivity; native photon and sampling normalization | None |
| Standard optical run, `waveder-optics` | PAW occupied-to-empty circular transition strengths in Å² and Gaussian spectral densities in Å²/eV | None |

Start with the [ordinary WAVECAR example](../examples/features/circular-dichroism/).
The optional PAW route below includes the longitudinal optical matrix elements
already calculated by VASP. It does not require the instrumented producer used
for full velocity and spin matrices. Increasing a bare-momentum calculation's
band count cannot supply missing PAW operator terms; compare the two routes
with matched states, transition windows and explicitly stated normalizations.

These calculations describe independent-particle transitions. Neither gives
photoluminescence polarization, exciton dynamics or an absolute absorption
coefficient. VASP calculates the electronic structure; VASPBERRY analyzes its
matrix elements and polarization channels.

## Prepare one standard optical run

The currently validated binary producer is **VASP 5.4.4**, with:

```text
LOPTICS = .TRUE.
LPEAD = .FALSE.
LNABLA = .FALSE.
LREAL = .FALSE.
ISYM = -1
NSW = 0
```

Retain `WAVECAR`, `WAVEDER`, `INCAR` and the completed `OUTCAR` together.
The adapter requires electronically converged states and a completed optical
calculation. Hybrid, meta-GGA, PEAD, gamma-real and other unvalidated branches
are rejected. Choose the source `NBANDS` to contain accurately converged empty
states **and at least one band beyond the selected final window**. The extra
band checks whether the final window cuts a degenerate group.

Use the [VASP WAVEDER description](https://vasp.at/wiki/WAVEDER) and
[LOPTICS instructions](https://vasp.at/wiki/LOPTICS) to prepare the optical run.
The standard output contains occupied-to-empty elements in Å; these already
include their transition-energy denominator. The adapter does not divide by
that energy again. The supported longitudinal PAW matrix follows
[Gajdoš et al.](https://doi.org/10.1103/PhysRevB.73.045112).

Paths and selected points are accepted. A full BZ mesh is needed for a BZ
integral, but this command produces point spectra and performs **no BZ
integration**. Printed OUTCAR k weights are retained for reference, never
applied or normalized by this command. The occupied/empty partition must be
consistent with source occupations and `NELECT`; metallic fractional filling
is outside this adapter's scope.

## Run the analysis

Set `N_OCC`, `V_FIRST`, `V_LAST`, `C_FIRST` and `C_LAST` from the actual VASP
run. All band indices are one-based and inclusive. For example:

```bash
python3 tools/vaspberry_kubo.py waveder-optics \
  --run-dir path/to/optics-run --occupied "$N_OCC" \
  --spinor-components 2 --energy-reference 'unchanged VASP eigenvalue zero' \
  --initial "$V_FIRST" "$V_LAST" --final "$C_FIRST" "$C_LAST" \
  --beam-vector 0 0 1 --polarization-axis 1 0 0 \
  --photon-min 1 --photon-max 3 --photon-num 401 --sigma-eV 0.05 \
  --formats csv dat npz --output-dir results/paw-optics
```

The standalone command `python3 tools/waveder_optics.py` takes the same options.
Use `--spinor-components 1` for a scalar calculation and `--spin` to select a
collinear spin channel. Intensities refer to that one channel; the command
applies no additional spin factor. A scalar spin-degenerate total has a factor
of two, while SOC spinors already represent both spin components. This common
factor does not affect helicity selectivity.

## Polarization, signs and degeneracies

The Cartesian beam direction is normalized to **b**. The projection of
`--polarization-axis` onto its transverse plane defines **e₁**, and
**e₂** = **b** × **e₁**. Thus:

```math
\boldsymbol\epsilon_\pm=(\mathbf e_1\pm i\mathbf e_2)/\sqrt2,
\qquad
\mathbf E(\mathbf r,t)=\mathrm{Re}
  [E_0\boldsymbol\epsilon_\pm e^{i(\mathbf q\cdot\mathbf r-\omega t)}].
```

With final-state bra $`f`$, initial-state ket $`i`$, and the stored length
matrix $`C^a_{fi}`$, a complete initial group $`I`$ and final group $`F`$ give

```math
I_\pm(\mathbf k;I,F)=\sum_{i\in I,f\in F}
  \left|\sum_a\epsilon_{\pm,a}C^a_{fi}(\mathbf k)\right|^2,
\qquad
\eta=\frac{I_+-I_-}{I_++I_-}.
```

The polarization vector is **not conjugated** in this contraction. The common
phase between the stored connection and dipole matrix cancels in the squared
modulus. Plus/minus are algebraic channel labels; observer-dependent
left/right naming is avoided. Reversing the beam swaps the channels for a
fixed transverse axis. Time-reversal partners exchange the channels when
complete corresponding groups are compared.

Adjacent states separated by at most `--degeneracy-threshold-eV` form
transitive groups on the full stored band axis. The default and minimum is
the standard producer's **0.002 eV** threshold. Both selected windows must
include complete groups at every point. For a centrosymmetric SOC bilayer,
select whole Kramers pairs rather than an arbitrary eigenvector in a pair.
If a window fails this check, expand it or request more VASP bands.

Strengths are summed before forming η. A group's transition energy is the
difference of the two group mean energies. This treats splittings inside the
declared threshold as unresolved and preserves group-rotation invariance;
the minimum/maximum transition energies and largest group spread are saved.
Converge this approximation if such small splittings matter for the spectrum.

## Outputs and interpretation

| Output | Contents |
|---|---|
| `transitions.csv`, `.dat`, `.npz` | One row per k and pair of complete groups; band bounds, centroid and range of transition energies, I₊/I₋ in Å², η and validity mask |
| `spectra.csv`, `.dat`, `.npz` | One row per k and photon energy; summed normalized-Gaussian spectral densities in Å²/eV, η and validity mask |
| `optical.json` | Polarization convention, source and operator scope, units, selected windows, broadening, grouping and validation records |

Each Gaussian has standard deviation `--sigma-eV` and unit area over the
whole energy axis. The finite displayed energy range is not renormalized.
The spectra contain no extra photon-energy factor, path-length normalization,
cell-volume factor or k weight. Their absolute intensity is therefore not the
native `-cd 2` intensity scale. Compare selection rules or explicitly convert
the two normalizations before comparing numerical intensities.

The default mask requires the total intensity to exceed $`10^{-10}`$ of the
largest total at the same k, using separate maxima for transition and spectral
tables. Zero intensity is always invalid. Such rows have `eta_valid=false`
and `eta=NaN`; do not plot them as η = 0. Adjust
`--relative-intensity-floor` when interpreting very weak transitions.

The adapter verifies WAVECAR/OUTCAR dimensions, geometry, energies,
occupations, spin layout and WAVEDER coverage. WAVEDER does not store energies
or coordinates: keeping the same run's files together remains the user's
source-association assertion. Passing these checks establishes consistency,
not k-mesh, empty-band or broadening convergence. The tests also check group
unitary invariance, time-reversal exchange, circular/linear sum rules and the
shared Hall-curvature sign convention.
