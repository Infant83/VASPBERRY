# Point-curvature and Hall output formats

Version 1.3.0 introduces `vaspberry.band-curvature` schema version **1** and
`vaspberry.hall-spectrum` schema version **1**. Schema versions describe file
contracts and are independent of the software version. These formats do not
replace the existing plaquette-flux format of the Fukui–Hatsugai–Suzuki (FHS)
method, also called the Fukui method, for Chern numbers. The
`VASPBERRY_Z2_FIELD` format belongs to the separate Fukui–Hatsugai (FH)
Z₂ n-field method and is also unchanged.

## Standard native WAVEDER curvature (from 1.6.5)

`--task kubo` now defaults to optical connections. With `--bands` omitted its schema is
`VASPBERRY_WAVEDER_KUBO_OCCUPIED_V1`, default file `KUBO_WAVEDER.csv`, with
columns `spin,k_index,kx_frac,ky_frac,kz_frac,omega_z_A2,min_external_gap_eV`.
Coordinates are fractional; `omega_z_A2` is Cartesian Ωxy in Å².
Comments record `source_operator=VASP_5.4.4_LONGITUDINAL_PAW_OPTICAL`,
`complete_occupied_bundle=true`, `occupied_bands_spin_N`, source counts and
completion status. The same numerical trace columns do not make this the
same physical operator as the canonical trace schema. When a
full mesh is explicitly supplied and validated, `total_chern` and
`sheet_hall_e2_over_h` report unrounded integrals. The latter is the negative
physical-spin-weighted Chern integral, in e²/h. A path has no integral.
An explicit `--bands` singleton, range or list always writes selected
geometry, even when the selection equals the occupied space:

| Native schema | Columns | Interpretation |
|---|---|---|
| `VASPBERRY_WAVEDER_KUBO_BUNDLE_V1` | The seven occupied-trace columns above | Selected geometric trace, one row per source k/spin. |
| `VASPBERRY_WAVEDER_KUBO_BAND_V1` | `spin,k_index,band,kx_frac,ky_frac,kz_frac,energy_eV,omega_z_A2,min_gap_eV` | Separate selected bands with explicit `--per-band 1`, when individually resolvable. |

Selected metadata records `band_ids`, rank, all source NBANDS, the stored
NDBANDS rectangle, required-pair/producer-cluster checks, and
`occupation_weighting=NONE`. Its explicit full-mesh integrals are
`selected_chern_spin_N` and `selected_chern_sum`, with unit geometric weights;
physical spin multiplicity is recorded but is not applied. There is no
`sheet_hall_e2_over_h` total-Hall claim. Missing or producer-erased pairs
with unequal selection weights stop the calculation. If all source bands
are selected, the zero trace is labeled `zero_trace_scope=TRUNCATED_WAVEDER_BASIS`;
it is not an infinite-basis material invariant. Output readers preserve these
scopes and plot labels use the actual selected IDs.

No legacy canonical DAT companion is produced. The historical band/bundle
and pair schemas below apply to the explicit `--kubo-source wavecar` route.
Preserve the metadata to distinguish operators; do not relabel old CSVs.

## Read the file that matches your question

The Intel and GNU executables write the same numerical formats. The
[Intel MPI hands-on guide](HANDS_ON.md) gives the producing commands.

| File | Contains | Calculation still needed? | Typical plot |
|---|---|---|---|
| Native trace CSV (`KUBO_WAVEDER.csv` for standard WAVEDER; `KUBO.csv` for the explicit WAVECAR approximation) | k coordinates and computed Ωxy in `omega_z_A2` (Å²), selected-bundle gap (eV); metadata identifies the operator | None for a curvature plot | Ω versus path samples; a full-mesh Ω map with the matching lattice |
| Native individual-band CSV | Band ID, energy, coordinates, Ωxy and minimum gap for each selected isolated band | None for a curvature plot | Separate band curves, filtered by `band` and `spin` |
| Native `SPIN_CHERN.csv`, `SPIN_BERRY.csv`, `SPIN_SPECTRUM.csv` | Projected-spin sector invariants, plaquette flux and spin spectrum | None for plotting | Sector flux maps, spin-gap maps and mesh comparisons; [spin Chern number guide](SPIN_CHERN.md) |
| Native `SPIN_KUBO.csv`, `SPIN_KUBO_SPECTRUM.csv`, optional `SPIN_KUBO_INTEGRAL.csv` | Positive/negative/parent point curvature, its decomposition and raw explicit-mesh integrals | None for plotting; use separate `spin-chern` for geometric integers | Spin-sector curves on a k path and full-mesh maps; [spin-sector Kubo-formula guide](SPIN_KUBO.md) |
| Native `PAIRS.csv` | Undivided interband products `numerator_*_eV2_A2`, two band energies and gap | Yes: occupations, squared-gap denominators and BZ integration | Intermediate pair analysis; it is not a conductivity table |
| `pairs.npz` + `pairs.json` | Validated typed cache of those same pairs, mesh and provenance | Yes: `pair-hall` calculates transport | Reusable input for μ/T/region scans |
| `conductivity.csv`, `.dat` or `.npz` + `.json` | σ, Δσ, μ, T, region and represented carrier count | None for plotting existing rows | Hall versus μ, temperature comparison, regional contrast or Hall versus carrier count |
| `character.csv/npz/json` | State/group raw charge and spin projections | `procar_character.py hall` for Hall attribution | Atom/layer/orbital/spin character of states |
| `character_hall.csv/npz/json` | Completed selected-band charge-Hall attributions | None for plotting existing rows | Group/spin-projection contribution versus μ |

**Import, integrate and plot are different operations.** `import-pairs`
validates/repackages the native data; `pair-hall` performs a physical
occupation-weighted Kubo-formula response calculation; `plot_hall.py` only plots the completed
conductivity table. `wavecar-hall` orchestrates the first two together with
the native executable and retains their intermediate files. Similarly,
PROCAR `project` reads state projections, `hall` calculates attribution,
and `plot` draws figures.

Native CSVs have comment lines beginning `#`, followed by named comma-separated
columns. Skip those comments when importing into Origin, a spreadsheet, R or
another plotter. For the bundle file, select `spin=1` and plot `k_index`
against `omega_z_A2`. Python is optional for that native-result plot. A
Cartesian BZ map additionally needs the matching reciprocal lattice; a line
path does not become a full-zone map by interpolation. Keep the comments or
JSON sidecars with any exported table.

<a id="native-projected-spin-chern"></a>

## Native projected-spin Chern numbers

The native [`spin-chern` task](SPIN_CHERN.md), available in 1.6.1, writes three CSVs
with schema `VASPBERRY_SPIN_CHERN_V1`. `--output NAME` changes the default
`SPIN` filename prefix. Existing outputs are rejected. Each file has `#`
comment metadata, named columns and a final `# result_status=PASS` on success.
They are separate from the spin-current conductivity tensors below.

| File | Data and interpretation |
|---|---|
| `SPIN_CHERN.csv` | One row: `C_plus,C_minus,C_charge,delta_C,C_spin,C_parent,rank_plus,rank_minus`; `C_charge=C_plus+C_minus`, `delta_C=C_plus-C_minus`, `C_spin=delta_C/2`. `C_parent` is evaluated independently from the parent links. |
| `SPIN_BERRY.csv` | `k_index,k1,k2,k3,flux_plus_rad,flux_minus_rad,flux_parent_rad,plaquette_area_A-2,omega_plus_A2,omega_minus_A2`. Fractional coordinates are plaquette centers; `omega` is flux divided by area. |
| `SPIN_SPECTRUM.csv` | `k_index,k1,k2,k3,spin_eigen_index,projected_pauli_eigenvalue,rank_plus,rank_minus,external_energy_gap_eV,min_gram_eigenvalue,gram_condition`. One row per selected-space spin eigenvalue at each source k point. |

The projected-Pauli eigenvalues are dimensionless. Their distance from zero
is distinct from the parent energy gap in eV. Metadata records
`representation=WAVECAR_PSEUDO_GRAM`, `physical_paw_validated=false`,
`scope=selected_subspace`, the axis/frame, thresholds and minimum gaps,
Gram/link conditioning, plane-wave overlap coverage and maximum flux.
For a selection starting at band 1, `global_indirect_gap_eV` and
`has_global_energy_gap` are also given. The separate
`source_occupations_match_integer_projector` checks whether the source
occupations match the selected filled space. Finite VASP smearing can make
that diagnostic false; it does not reject the geometric calculation of an
isolated selected subspace.
The two sectors are subspaces of the complete selected spinor space, not
independent `ISPIN=2` channels or per-band PROCAR spin weights.

<a id="native-projected-spin-kubo-curvature"></a>

## Native projected-spin Kubo-formula Berry curvature

The native [`spin-kubo` task](SPIN_KUBO.md), available in 1.6.1, writes schema
`VASPBERRY_SPIN_SECTOR_KUBO_V1`. It uses the same spin axis, band selection
and pseudo Gram metric as `spin-chern`, with a canonical-momentum derivative
approximation. It accepts an arbitrary k list by default. Only an explicit
`--mesh NX,NY` requests a validated complete-mesh integral. The default
prefix is `SPIN`; `--output NAME` changes it. Existing outputs are rejected.

| File / column family | Meaning |
|---|---|
| `SPIN_KUBO.csv`: `k_index,k1,k2,k3` | One-based source index and fractional source k coordinates; these are points, not plaquette centers |
| `kx_A-1,ky_A-1,kz_A-1,k_distance_A-1` | Cartesian coordinates and cumulative distance along the supplied order, in Å⁻¹ |
| `omega_plus_*_A2,omega_minus_*_A2,omega_parent_*_A2` | Total sector and parent trace curvature; `*` is `yz`, `zx` or `xy`, in Å² |
| `parent_projected_plus_*_A2,parent_projected_minus_*_A2` | The parent tangent contribution resolved in the spin eigenframe |
| `spin_mixing_plus_*_A2,spin_mixing_minus_*_A2` | Contribution from the k dependence of the projected-spin sectors; adds to the preceding contribution and cancels in their sum |
| `rank_plus,rank_minus` | Dimensions of the two spin sectors, not counts of independent collinear channels |
| `external_energy_gap_eV,min_abs_spin_eigenvalue` | Sampled isolation from excluded bands (eV) and distance of the projected Pauli spectrum from zero (dimensionless) |
| `min_gram_eigenvalue,gram_condition` | Selected pseudo-wavefunction metric conditioning |
| `selected_energy_min_eV,selected_energy_max_eV` | Energy bounds of the selected source bands at that point |
| `SPIN_KUBO_SPECTRUM.csv` | One row per selected spin eigenvalue: `k_index,k1,k2,k3,spin_eigen_index,projected_pauli_eigenvalue,rank_plus,rank_minus` |
| `SPIN_KUBO_INTEGRAL.csv` | Explicit-mesh only: `C_est_plus,C_est_minus,C_est_charge,delta_C_est,C_est_spin,C_est_parent,rank_plus,rank_minus` |

The dimensionless `C_est` columns integrate the point curvature with the
**oriented** reciprocal plaquette area. `C_est_charge=C_est_plus+C_est_minus`,
`delta_C_est=C_est_plus-C_est_minus` and `C_est_spin=delta_C_est/2`.
They are raw estimates and are never rounded to integers. The integral file
does not exist after an ordinary path calculation. A mesh's cumulative
input-order distance is not a high-symmetry band-path coordinate.

Metadata records the source filenames, spin frame, selected bands, numerical
thresholds and the derivative/representation contract:
`derivative=CANONICAL_MOMENTUM_DERIVATIVE_PROXY`,
`representation=PSEUDO_GRAM`, `exact_pseudo_projector_derivative=false`,
`physical_paw_validated=false`, `geometric_integer_certified=false` and
`conventional_spin_current_hall=false`. `source_nbands` identifies the original
WAVECAR band count; `sum_band_max` is the upper band included in the sum,
and `summed_external_bands` counts its states outside the selected group.
With `--sum-bands N`, a reduced sum also records `sum_boundary_gap_eV`.
Parent isolation continues to use all source states. All retained bands
outside the selected group enter the sum by default; fixed-rank projectors do not apply
source occupations. A final `result_status=PASS` reports completed numerical
checks, not operator completeness or k-mesh convergence.

Plot `k_distance_A-1` against `omega_plus_xy_A2` and `omega_minus_xy_A2`
for an xy-plane path, or the appropriate tensor component for another plane.
The [graphene](../examples/materials/graphene-spin-chern/kubo/) and
[Bi](../examples/materials/bi-spin-hall/spin-chern-kubo/) examples show these
plots beside full-zone Chern numbers and Berry fluxes from the Fukui method, retaining failures and unconverged
integrals rather than assigning unsupported Chern numbers.

## PROCAR character and charge-Hall attribution

The [PROCAR workflow](../examples/features/procar-character/) uses
`tools/procar_character.py` and writes version-1 schemas separate from the
ordinary Hall-spectrum schema. Matching SOC `PROCAR`, `WAVECAR` and `OUTCAR`
are required. The current parser accepts the four-block noncollinear
`LORBIT=11` format. Atom IDs and band IDs are one-based; atom groups and
optional orbital labels come from the user's group JSON.

`project` writes `character.csv`, `character.npz` and `character.json`
(`vaspberry.procar-character`). For K k points, N bands and G groups:

| Array or column | Meaning |
|---|---|
| `kpoints_fractional`, `energies_eV`, `occupations` | K×3 coordinates and K×N source energies/occupations |
| `weights`, `lattice_A`, `reciprocal_inv_A` | Source K weights and direct/reciprocal row vectors; Å and Å⁻¹ with 2π |
| `raw_cartesian` | K×N×G×4 raw projected charge and Cartesian Pauli weights, ordered `charge,mx,my,mz` |
| `characters` | K×N×G×4 values ordered `charge,pauli_axis,plus,minus` |
| `all_ions_cartesian`, `printed_state_cartesian` | K×N×4 sum over the per-ion totals and the independently rounded printed state total, respectively |
| CSV `k_id`, `band_id`, `group` | One row per source state and user group, with coordinates, energy and occupation |
| CSV `charge`, `pauli_axis`, `plus`, `minus`, `mx`, `my`, `mz` | The same dimensionless raw projector weights in readable form |

The chosen Cartesian unit axis is `axis_cartesian`. The actual
SAXIS-to-Cartesian transformation is read from OUTCAR and retained as
`spin_to_cartesian`. For group charge q and axis-projected Pauli weight m,
`plus=(q+m)/2`, `minus=(q-m)/2`, and projected spin in units of ħ is m/2.
No weight is clipped or normalized to unity. PAW atom/orbital projections
can be incomplete, and user groups may overlap; they are not automatically
a partition of the full wavefunction. The metadata records group definitions,
input hashes, matching checks and an explicit same-run association.
`projection_diagnostics.csv` records `k_id`, `band_id`, `all_ions_charge`,
`charge_residual=1-all_ions_charge`, Cartesian `mx,my,mz` and
`printed_state_charge`. The residual is bookkeeping relative to unit weight,
not a separately calculated interstitial projector.

`hall` consumes that cache and matching Kubo pairs, then writes:

| File | Contract |
|---|---|
| `character_hall.csv`, `.npz` | Identical long-form columns `mu_eV`, `mu_minus_reference_eV`, `temperature_K`, `region`, `group`, `component`, `attribution_e2_over_h`, `delta_attribution_e2_over_h` |
| `character_hall.json` | `vaspberry.character-hall`: selected bands, full source virtual-band window, minimum isolation gap, operator, axis, groups, regions, reference μ, formula and input/output checksums |
| `selected_curvature.npz` | K×B×3 `omega_A2` for the B selected bands, K×B `min_gap_eV` and `energies_eV`, k/band IDs and fractional coordinates |

The response is the charge-Hall integral weighted by each selected band's
raw character. Δ attribution applies `f(mu,T)-f(mu_reference,T)` before
summation. The `pauli_axis` channel has charge-attribution units e²/h; it
is not a conventional spin Hall coefficient. Selected-band outputs omit
all other bands, even if those bands are fully occupied. Every selected
band must be isolated from all other source states at each sampled k;
unresolved degeneracies are rejected, not coalesced for this decomposition.
Use ordinary `pair-hall` for the full occupation-weighted charge response.
Three reserved diagnostic groups have only a `charge` component:
`$unweighted`, `$all_projected`, and `$unprojected_residual`. The latter two
sum to the unweighted selected-band baseline. This closure does not make
arbitrary overlapping user groups a partition, and no unprojected spin is
inferred from charge residuals.

The `plot` command reads these saved products and writes `character.*` and
`character_hall.*` in PNG, PDF and SVG, plus `plots.json`. It does not run
VASP or the VASPBERRY executable. The CSVs and NPZs are also ordinary numerical inputs
for an independent plotting program; keep their JSON metadata beside them.

## Spin operators and insulating spin Hall tensors

All schemas below are version 1. Arrays use float64 or complex128; integer
band indices are int64 and one-based. Spinor multiplicity is one.

`spin-export` validates the opt-in VASP producer and writes
`physical-matrices.npz/json` (`vaspberry.spin-velocity`). Its arrays are:

| Array | Shape / units |
|---|---|
| `energies_eV` | K×N eigenvalues, unchanged source zero |
| `kpoints_fractional`, `weights` | K×3 coordinates and K normalized weights |
| `lattice_A`, `band_indices` | 3×3 direct row vectors in Å; all N source band IDs |
| `spin_pauli` | K×3×N×N complex dimensionless Pauli matrices, Cartesian x,y,z |
| `velocity_eVA` | K×3×N×N complex physical ħv matrices, eV Å |
| `overlap` | K×N×N complete physical PAW overlap |

The [n,m] convention is bra n, ket m. Source identity, Cartesian axes,
PAW operator scope, no normalization, complete source-band coverage and
availability of diagonal/degenerate velocity blocks are mandatory metadata.
The archive checksum binds that metadata to its arrays. Source consistency
and numerical validation do not establish band or mesh convergence.

`spin-matrix` writes `spin.npz/json` (`vaspberry.spin-operators`). It preserves
`spin_pauli`, `overlap`, `pseudo_spin_pauli`, `pseudo_overlap`, k points,
eigenvalues, lattice and band indices. Scope is either `raw_pseudo` or
`paw_augmented`; the latter requires a matching validated augmentation input.
`augmentation.npz/json` (`vaspberry.spin-augmentation`) holds
`delta_spin_pauli`, `delta_overlap` and the matching coordinates, energies,
lattice and band indices. The raw pseudo overlap is never forced to identity.

`spin-hall` writes `spin_hall.json` (`vaspberry.spin-hall`) and selected
CSV/DAT/NPZ products:

- `spin_hall.csv/dat`: 27 rows with `spin`, `current`, `electric_field`,
  `sigma_hbar_over_e_e2_over_h`, `sigma_hbar_over_e_S`.
- `spin_curvature.csv/dat`: each k point and all 27 Cartesian components,
  including fractional coordinates, normalized weight and `omega_spin_A2`.
- `spin_hall.npz`: K×3×3×3 `omega_spin_A2`, K×3×3 `omega_charge_A2`,
  3×3×3 integrated spin tensors, 3×3 `charge_sigma_e2_over_h`, coordinates,
  weights, lattice and scalar `mu_eV`.
- `run.json`: completion, elapsed time and processed points. A failed
  integration retains a `.partial` directory instead of a complete output.

Tensor order is **spin, current, electric field**. The spin sheet units are
(ħ/e)(e²/h) or (ħ/e)S; these are two representations of the same coefficient.
They are not e²/h charge conductivity. Metadata records both source and
retained band counts, the finite-band current-product approximation, occupied
count, sampled gaps, chemical potential, mesh and plane normal. In-plane
current/field contractions give sheet transport; no vacuum thickness is used.
See [formulas and sign conventions](SPIN_HALL.md).

`wannier-edge` writes `edge.json` (`vaspberry.wannier-edge`) and selected
CSV/DAT/NPZ. NPZ contains K parallel fractional coordinates, K×(width·N)
`energies_eV`, `left_edge_weight`, `right_edge_weight`, and the lattice.
The weights are probabilities in the chosen boundary Wannier cells, not
real-space spin densities. CSV/DAT has one row per k point and strip state.
The metadata states the cut axes, width, edge-cell count and ideal-termination
scope. Sum weights over a degenerate group or use an edge spectral function.

## Native Kubo pair export and reusable cache

`-kubo 2 -kubo_pairs PAIRS.csv` writes
`VASPBERRY_BARE_MOMENTUM_KUBO_PAIRS_V1`. Each row describes one unordered
pair `n_band < m_band` at a source k point and spin channel:

| Columns | Meaning |
|---|---|
| `spin`, `k_index`, `n_band`, `m_band` | One-based source indices |
| `kx_frac`, `ky_frac`, `kz_frac` | Source reduced coordinates |
| `energy_n_eV`, `energy_m_eV`, `gap_eV` | Raw energies and absolute energy difference |
| `numerator_yz_eV2_A2`, `numerator_zx_eV2_A2`, `numerator_xy_eV2_A2` | Undivided `−2 Im(D_a,nm D_b,mn)` in eV² Å² |

Comments declare the real and reciprocal lattices, source dimensions,
component order, units, normalization and operator. No occupation, spin
multiplicity or denominator is included in the numerators. Only a final
`result_status=PASS` footer certifies that the native exporter completed.
The pair schema remains V1. Native pair export writes `PATH.partial` and
publishes `PATH` only after successful completion, flush and close; existing
final or partial paths are refused.
The importer additionally verifies every selected-spin row against WAVECAR;
a partial export is not an integration input.

The three real pair numerators are not the full complex velocity matrices.
For an isolated band, curvature includes `Nab/(E_n-E_m)^2`; occupation-weighted
pair transport additionally uses `f_n-f_m` and the k/BZ weights. The separate
native `--task kubo --kubo-source wavecar --curvature-csv KUBO.csv` mode has already performed the
curvature sum, which is why its Å² column is directly plottable.

`import-pairs` writes `pairs.npz` and `pairs.json`, schema
`vaspberry.kubo-pairs` version 1. For K k points, B stored bands and
P=B(B−1)/2 pairs, the arrays are:

| Array | Shape / meaning |
|---|---|
| `k_ids`, `band_ids` | Integer one-based consecutive source IDs |
| `pair_n`, `pair_m` | Integer zero-based band-array positions, lexicographic n<m |
| `kpoints_fractional`, `weights` | K×3 coordinates and K normalized weights |
| `energies_eV` | K×B raw eigenvalues |
| `numerator_eV2_A2` | K×P×3 real values, component order yz,zx,xy |
| `lattice_A`, `reciprocal_inv_A` | 3×3 row vectors, reciprocal includes 2π |

Numerical floating arrays are float64. JSON records full-mesh geometry,
energy reference, multiplicity, source operator and the matching NPZ checksum.
Load NPZ with `allow_pickle=False`. A full uniform 2D mesh is required for
Hall integration; line paths and irreducible meshes are not accepted.
This import/cache step does not evaluate Fermi occupations or output a Hall
coefficient. The subsequent `pair-hall` stage writes the completed transport
tables described under [Hall directory](#hall-directory).

`--pair-band-max M` selects a finite pair space from the complete cache without
changing its source arrays. Hall metadata keeps `source_nbands`, records
`pair_band_window=[1,M]` and the cutoff gap. A reduced virtual-state window is
`truncated_pair_space`, and the full stored window is `all_source_bands`.
If an occupied upper cutoff is explicitly allowed with `--allow-partial-bands`,
the scope is instead `partial_band_contribution`. M never replaces the source
VASP band count.

`pair-hall` and `bundle-hall` use the existing Hall-spectrum schema below.
`wavecar-hall` places native logs in `native/`, reusable arrays in `pairs/`,
and final numerical results in `hall/`. Its `workflow.json` records overall
status and an error on failure; native subprocess logs give execution details.
A failed workflow does not create a valid Hall result.

`kubo-hall` and `waveder-hall` also use the Hall-spectrum schema.
With `--occupied N` the method is `standard_waveder_occupied_bundle_T0`;
metadata identifies the longitudinal
PAW optical operator, producer settings, occupied count, stored empty-band
coverage, source files and checked gap. These tables describe a constant
zero-temperature response within the same global insulating gap.

With `--bands` the method is `standard_waveder_selected_weighted_pairs` and
`scope=selected_band_contribution`. Metadata includes `selected_band_ids`,
all `intermediate_band_ids`, `rectangular_coverage`, `producer_cluster_checks`,
and `total_AHC_certified=false`, `total_delta_AHC_certified=false`. Requested
μ/T values are accepted only if every noncanceling pair is available and
not erased by the producer. Finite-temperature tails cannot waive coverage.
The `total` region means the full sampled reciprocal region of the selected
contribution; it does not mean total material AHC. `plot_hall.py` checks the
table hash against the selected metadata and labels exported figures
`Selected-band contribution: ...`; `plot.json` retains that scope and IDs.

## Full-connection Wannier operators and results

`wannier-import` writes `operators.npz` and `operators.json`, schema
`vaspberry.wannier-operators` version 1. The input is the paired effective-model
HH_R/AA_R format with translation weights already absorbed. It is not the
usual Hamiltonian-only hr.dat format. The arrays are:

| Array | Shape / units |
|---|---|
| `irvec` | R×3 int64 integer translations |
| `hamiltonian_eV` | R×N×N complex128 Hamiltonian in eV |
| `position_A` | R×3×N×N complex128 position connection, Cartesian x,y,z in Å |
| `lattice_A` | 3×3 float64 direct row vectors in Å |

Metadata declares the positive Fourier phase, unit real-space degeneracy,
model dimensions, energy zero, spin convention and source integrity. An
imported cache establishes format consistency, not model convergence.

`wannier-hall` writes the common Hall schema in `hall/`, with method
`vaspberry_wannier_full_connection_T0`. The zero-temperature scan has one
fixed occupied model bundle; electron counts describe this represented model,
not omitted deep DFT bands. `sigma_terms_e2_over_h` is 3×3, with rows J0/J1/J2
and columns yz/zx/xy. The scalar Hall table projects their sum onto the
oriented plane normal. Spin multiplicity is explicit.

The accompanying `curvature.npz` contains K×3 `kpoints_fractional`, K
normalized area `weights`, K integer `parent_cell` IDs, K×3×3
`omega_terms_A2`, K×3 `omega_A2`, K×2 valence/conduction `band_edges_eV`,
and the direct/reciprocal lattices. A complete parent-cell partition replaces
coarse cells with refined children; a partial point cloud is not renormalized.
Internal occupied/empty degeneracies require no individual-band curvature.

`wannier-bands` writes CSV and/or NPZ with all represented model bands.
`bands.npz` stores K×3 `kpoints_fractional`, K `distance_inv_A`, K×N
`energies_eV`, direct lattice and path vertices/ticks. `bands.json` uses
schema `vaspberry.wannier-bands` version 1 and records path labels, source,
units and output formats. See [the workflow guide](WANNIER_TRANSPORT.md).

## Native Kubo individual-band CSV

For explicit `--kubo-source wavecar` with `--task kubo`, `kubo-line` or
`kubo-integral`, select `--bands N` for one
band, or `--bands FIRST:LAST --per-band 1` for separate bands. Add
`--curvature-csv PATH` to save the full-precision CSV as well as the legacy
DAT outputs; this individual-band CSV has no automatic filename.
Its schema is `VASPBERRY_BARE_MOMENTUM_KUBO_V3`, with columns:

```text
spin,k_index,band,kx_frac,ky_frac,kz_frac,energy_eV,omega_z_A2,min_gap_eV
```

IDs are one-based. Coordinates are fractional in the source reciprocal
basis, energy and gap are in eV, and `omega_z_A2` is Cartesian Ωxy in Å².
Every selected band must be isolated from every other stored band by more
than 1e-5 eV at all source k points and spin channels. The program checks
this before writing numerical results. A multi-band trace has the distinct
schema below and does not contain separate band energies or curves.

Version 3 records `source_nkpoints`, `source_nspin`, `source_nbands`,
`band_min`, `band_max`, `band_rank` and `expected_rows`. The expected count is the number of source k points
× spin channels × selected bands. The header starts with
`result_status=INCOMPLETE`; only a terminal `result_status=PASS` after all
successful writes marks completion. Readers check complete, unique coverage
and finite numerical values. Older V2 files retain their original metadata
and are accepted as `legacy_unverified`, without this completion guarantee.

The exporter writes `PATH.partial` while calculating and publishes the final
`PATH` only after successful completion, flush and close. An existing final
or partial filename is refused. A leftover partial file is a failed-run
artifact, not a numerical result; inspect its log before rerunning elsewhere.

## Native Kubo bundle CSV

Native `--task kubo --kubo-source wavecar --bands FIRST:LAST --curvature-csv PATH`, with two or
more selected bands and no `--per-band 1`, writes
`VASPBERRY_BARE_MOMENTUM_KUBO_BUNDLE_V2`. It contains one row per spin channel
and source k point, with comment metadata followed by a CSV header. The
file stores the trace curvature of the selected band range, with unit
occupation. It has no per-band energy column or legacy `.dat` companion.
The default trace filename is `KUBO.csv` when `--curvature-csv` is omitted.
The numerical columns and physical bundle meaning are unchanged. Version 2
adds explicit completion and coverage metadata; the removed native
`--bundle` flag is not needed to produce this format.

| Column | Meaning |
|---|---|
| `spin` | One-based WAVECAR spin-channel index; SOC spinors are already represented |
| `k_index` | One-based source k-point index |
| `kx_frac`, `ky_frac`, `kz_frac` | Fractional coordinates in the source reciprocal basis |
| `omega_z_A2` | Selected-bundle Cartesian Ωxy in Å² |
| `min_external_gap_eV` | Smallest selected-to-excluded source-band gap at this point |

| Metadata | Meaning |
|---|---|
| `result_kind` | `ISOLATED_BUNDLE_TRACE` |
| `result_status` | Initial `INCOMPLETE`; terminal `PASS` only after all spin/k results are written successfully |
| `source_nkpoints`, `source_nspin`, `expected_rows` | Source dimensions and required row count; trace rows = k points × spin channels |
| `normalization` | `STANDARD_MINUS_TWO_IM`, with `A_i=i<u\|d/dk_i u>` |
| `operator` | `WAVECAR_BARE_MOMENTUM_NO_PAW_NONLOCAL_VELOCITY` |
| `band_min`, `band_max`, `band_rank` | Selected one-based contiguous range and its size |
| `source_nbands` | All retained source bands, including excluded intermediate states |
| `intermediate_bands` | `EXTERNAL_TO_SELECTED_BUNDLE_WITHIN_SOURCE_NBANDS` |
| `internal_transitions` | `EXCLUDED_ANALYTICALLY` |
| `gap_threshold_eV` | `1e-5`; every external gap must exceed this value |
| `occupation_weighting` | `NONE`; no Fermi or multiplicity weighting is added |
| `no_external_states` | `true` only when every retained source band is selected |

All spin channels and k points are checked before opening the new output
file. Internal degeneracies are allowed; unresolved external gaps reject the
run. Selecting all retained source bands produces zero in that truncated
basis, with `min_external_gap_eV=NA` and
`zero_trace_scope=TRUNCATED_WAVECAR_BASIS`. This is not evidence that omitted
higher states make no physical contribution.

Readers require the complete expected spin/k coverage, unique row IDs,
finite coordinates/curvature and a terminal PASS. The `NA` external-gap
sentinel is reserved for the declared all-source-band trace above. Older
trace V1 files remain readable as `legacy_unverified`; an old header PASS
cannot establish that all records were successfully produced. Preserve the
original files and rerun when the V2 guarantee is required.
Trace export uses the same `.partial` staging and final-publication rule as
the individual-band exporter above.

Keep the matching POSCAR or WAVECAR to recover the reciprocal lattice for
Cartesian plots. This CSV does not carry energies, integration weights or
the lattice, and it is not the `vaspberry.band-curvature` NPZ/JSON input of
the `hall` command. The [native bundle guide](KUBO_TRANSPORT.md#native-curvature-of-a-band-bundle)
explains the calculation and the occupied-state interpretation.

## Curvature directory

The producer writes a new directory containing:

| File | Purpose |
|---|---|
| `curvature.npz` | Typed numerical arrays; load with `allow_pickle=False` |
| `curvature.json` | Meaning, units, normalization, sampling, operator provenance and NPZ SHA256 |
| `curvature.csv` | One readable row per k point and selected band |

Keep the NPZ and JSON together. `data_npz_sha256` binds the JSON to the exact
NPZ bytes. The reader checks that checksum before accepting the arrays.
Output commands require a fresh directory; do not overwrite an earlier run
to conceal a failed or differently normalized attempt.

### NPZ arrays

K is the number of represented k points, N the selected output bands, and M
the number of declared intermediate-band IDs. IDs are global, positive,
one-based integers; array positions are zero-based NumPy indices.

| Key | dtype and shape | Meaning |
|---|---|---|
| `k_ids` | `int64[K]` | Unique source k IDs |
| `band_ids` | `int64[N]` | Unique output band IDs |
| `intermediate_band_ids` | `int64[M]` | Declared spectral-sum band window |
| `kpoints_fractional` | `float64[K,3]` | Reciprocal fractional coordinates |
| `weights` | `float64[K]` | Nonnegative k weights summing to one |
| `energies_eV` | `float64[K,N]` | Energies in the declared reference |
| `omega_A2` | `float64[K,N,3]` | Cartesian curvature vector, Å² |
| `valid_nondegenerate` | `bool[K,N]` | Numerical validity for the represented band calculation |
| `min_gap_eV` | `float64[K,N]` | Minimum checked gap; matrix mode checks all exported band energies independently of the m sum |
| `lattice_A` | `float64[3,3]` | Direct lattice vectors as rows, Å |
| `reciprocal_inv_A` | `float64[3,3]` | Reciprocal vectors as rows, including 2π, Å⁻¹ |

The components of `omega_A2` are `yz,zx,xy`, equivalently `x,y,z` of the
curvature vector. A value is finite only where both the band's validity mask
and the component-availability mask are true; all other curvature entries are
NaN, not zero. JSON itself contains no NaN values. The CSV uses blank cells
for unavailable/invalid curvature entries.

Minimum gaps are finite and nonnegative, and a valid individual-band entry
requires a strictly positive known gap. For matrix output,
`kernel_diagnostics.isolation_check_scope` records the all-exported-energy
check, while `all_source_band_energies_available` indicates whether the input
contained every source band energy. A truncated m sum does not hide a known
near-degenerate band. Imported legacy gaps retain their producer-declared
scope; they are not newly reconstructed by the importer.

### JSON contract

| Key | Required meaning |
|---|---|
| `schema`, `version` | `vaspberry.band-curvature`, `1` |
| `complete` | `true` for a finished data product; not a physical convergence certificate |
| `normalization` | `physical_Omega=-2Im;A=i<u|grad_k u>` |
| `axes` | `["k", "band", "cartesian"]` |
| `components` | `["yz", "zx", "xy"]` |
| `available_components` | Three booleans identifying supplied Cartesian components |
| `units` | Energy `eV`, curvature `Angstrom^2`, lattice `Angstrom`, reciprocal `1/Angstrom` |
| `weights_convention` | `sum_one` |
| `reciprocal_convention` | `2pi` |
| `source_nbands` | Original source band count, not merely the selected n window |
| `spin_multiplicity` | `1` for spinors or a single spin channel; `2` for declared scalar spin degeneracy |
| `method`, `energy_reference` | Nonempty descriptions of origin and energy zero |
| `sampling` | `points` or validated `uniform_full_2d`, as described below |
| `source_operator` | Operator identity and available accuracy/provenance information |
| `provenance` | Source files, hashes, generating command and other producer records |
| `convergence_status` | `not_established_by_format_validation` |
| `data_npz_sha256` | Exact NPZ checksum |

The exact `units` object is:

```json
{"energy":"eV","curvature":"Angstrom^2","lattice":"Angstrom","reciprocal":"1/Angstrom"}
```

For example, `sampling={"kind":"uniform_full_2d","mesh":[32,32],
"plane_axes":[0,1]}` denotes a full tensor-product mesh in the first two
reciprocal directions. `plane_axes` is zero-based and ordered, so reversing
it reverses the surface normal. The validator checks the grid, duplicate
points and uniform weights. A shifted full grid can be represented. A path,
partial patch or irreducible mesh must not be declared a full 2D grid.

`sampling={"kind":"points"}` permits storage and inspection of point data,
but the Hall integrator rejects it. Normalized weights alone do not establish
that a dataset covers the Brillouin zone.

## Generic interband-matrix input

The `matrix` command reads `vaspberry.interband-matrix`, schema version **1**,
as implemented in `tools/exported_matrix_kubo.py`. This is an input contract,
not the point-curvature output schema.

Required NPZ arrays are `k_ids`, `band_ids`, `kpoints_fractional`, `weights`,
`energies_eV`, `D_eVA`, `coverage`, `lattice_A`, and `reciprocal_inv_A`.
The scalar/vector arrays follow the corresponding types above, now with B
matrix bands. `D_eVA` is `complex128[K,3,B,B]` and stores
`<n|D_a|m>` with axes k, Cartesian, bra band, ket band. `coverage` is a Boolean
array of the same shape. Unavailable elements have both real and imaginary
parts NaN and must be uncovered; covered zeros remain valid measured zeros.

The input JSON binds the NPZ with `matrix_npz_sha256`. It declares the exact
units (`matrix: eV*Angstrom` in place of `curvature`), axes, bra/ket convention,
Cartesian order, diagonal coverage status, source k/band counts, normalized
weights, 2π reciprocal convention, and an explicitly Hermitian operator.
The operator object includes `kind`, `definition`, `included_terms`,
`missing_terms`, `validation_evidence`, `accuracy_status`, and `hermitian`.
Provenance includes `run_id`, `exporter_revision`, and named `source_hashes`.
Both matrix directions must be covered consistently; the reader does not
manufacture the reverse direction or symmetrize a failed matrix.

`accuracy_status=experimental` requires explicit diagnostic opt-in.
`accuracy_status=validated` requires written evidence but is still a
producer declaration. Numerical format/Hermiticity checks cannot independently
validate the physical exporter.

## Hall directory

Numerical outputs are independently selectable with `--formats csv dat npz`.
`hall` preserves its CSV+NPZ default; the new pair/bundle workflows default to
all three. `conductivity.json` is always written. CSV, tab-separated DAT
(commented column header) and NPZ use identical long-form rows, one per chemical
potential, temperature, region and band channel. JSON includes the selected
formats, source metadata, reciprocal area/normal, region definitions,
occupation-window checks and algorithm settings.

`plot_hall.py` reads any of these numerical formats and writes PNG, PDF and/or
SVG. Pair and bundle integrations use `band_id=0` for the total represented
subspace; they do not assign a gauge-dependent curvature to each degenerate band.

| Column / NPZ key | Meaning |
|---|---|
| `mu_eV` | Chemical potential in the input energy reference |
| `mu_minus_reference_eV` | μ minus the declared comparison reference |
| `temperature_K` | Fermi-function temperature |
| `region` | `total`, a user region, automatic `rest`, or a requested difference |
| `band_id` | Global band ID; `0` denotes the sum of represented bands |
| `sigma_e2_over_h` | Oriented intrinsic charge sheet response in e²/h |
| `sigma_S` | The same sheet conductance in siemens, not S/m |
| `delta_sigma_e2_over_h` | Response minus its value at the reference μ, at the same temperature |
| `electrons_per_cell` | Represented-state occupation sum, including declared spin multiplicity |
| `delta_electrons_per_cell` | Occupation change from the reference μ |

For a named difference channel, response and occupation columns are signed
differences, not standalone particle counts. At zero temperature the convention
is full occupation when `E==mu`. Temperature broadens occupations; it is not
a scattering-rate or spectral-denominator broadening.

The metadata's `delta_formula` records direct integration with
`f(mu,T)-f(mu_reference,T)`. Delta outputs avoid cancellation from subtracting
large filled-band baselines: zero-temperature sums start at the reference
occupation boundary, while finite-temperature chunks use occupation differences.
The reported mathematical quantity remains the response change at fixed
electronic states.

The point-curvature `hall` command requires a complete represented source-band
window by default. Pair integration also permits a reduced virtual-state
window 1:M, with its cutoff recorded as described above. Both require negligible
occupation of the highest included band over the requested μ/T range.
`--allow-partial-bands` labels an incomplete occupation window `partial_band_contribution`;
it does not add omitted occupied bands. A successful window check still does
not prove that the source calculation includes enough unoccupied states for
curvature convergence.

## Read and plot with other tools

The output tables are the numerical result. The supplied plotters are
optional readers: no VASPBERRY executable, VASP installation or Wannier
model is needed to plot a completed table. Keep its JSON sidecar with the
table so that the operator approximation, units, reference energy, band
window and region definitions remain available for interpretation.

CSV opens directly in spreadsheet applications, Origin, R or Python. DAT
has a tab-separated column header prefixed by `# `; strip that prefix when
using a reader that expects an uncommented header. NPZ contains the same
named columns and uses no Python object arrays. For example, run this from
the repository root after installing `requirements-transport.txt`:

```bash
python3 - <<'PY'
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

source = Path('examples/features/kubo-hall/reference/cases/24x24-source60-cap40')
with np.load(source / 'conductivity.npz', allow_pickle=False) as archive:
    d = {name: archive[name] for name in archive.files}

fig, ax = plt.subplots(layout='constrained')
for region in ['K', 'Kprime', 'valley']:
    selected = ((d['region'] == region) & (d['temperature_K'] == 300)
                & (d['band_id'] == 0))
    order = np.argsort(d['mu_minus_reference_eV'][selected])
    ax.plot(d['mu_minus_reference_eV'][selected][order],
            d['delta_sigma_e2_over_h'][selected][order], label=region)
ax.set(xlabel='μ − μref (eV)', ylabel='Δσ (e²/h)')
ax.legend()
out = Path('results/custom-hall-plot')
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / 'regional-hall.svg')
plt.close(fig)
PY
```

This selects already calculated 300 K rows. The same arrays support absolute
σ versus μ, several temperatures, a μ/T map, or response versus carrier
count. For the last plot, join `delta_electrons_per_cell` from **`region=total`
and `band_id=0`** to the desired response rows at the same μ and temperature.
A difference region's carrier column is a signed regional contrast, not the
total number of doped carriers. Counts are electrons per cell, not cm⁻²;
conversion to sheet density requires the physical in-plane cell area.

The equivalent plain CSV can be read without VASPBERRY or NumPy:

```python
import csv
with open('results/mos2-hall/hall/conductivity.csv') as handle:
    rows = list(csv.DictReader(handle))
curve = sorted((float(r['mu_eV']), float(r['sigma_e2_over_h']))
               for r in rows if r['region'] == 'total'
               and int(r['band_id']) == 0 and float(r['temperature_K']) == 300)
```

For native `PAIRS.csv` or bundle CSV, skip `#` comment lines before passing
rows to a CSV reader. Retain and inspect those comments separately: they
identify the operator. Current schemas require terminal `result_status=PASS`;
older curvature schemas remain explicitly legacy-unverified. Pair
numerators are not curvature or conductivity until occupations and energy
denominators have been applied. `import-pairs` and `pair-hall` perform those
validated numerical stages before plotting.

## PAW circular transition strengths

`waveder-optics` writes `optical.json` with schema
`vaspberry.circular-transition-strength`, version 1. Requested CSV, DAT and
NPZ files contain the same long-form columns; CSV/DAT use `NaN` for invalid
ratios, and NPZ has an explicit boolean `eta_valid` array.

- `transitions.*`: k index, fractional coordinates, unapplied source k weight,
  initial/final inclusive band bounds, `transition_eV` (group centroid
  difference), `transition_min_eV`, `transition_max_eV`, `I_plus_A2`,
  `I_minus_A2`, `eta`, `eta_valid`.
- `spectra.*`: k index, fractional coordinates, unapplied source k weight,
  `photon_eV`, `I_plus_A2_per_eV`, `I_minus_A2_per_eV`, `eta`, `eta_valid`.

Both initial and final windows contain complete near-degenerate groups;
their strengths are summed before forming the ratio. Spectra use normalized
Gaussians of the declared standard deviation and group-centroid transition
energies. No k weight, spin multiplicity, photon-energy or cell-volume factor
is applied. These are point transition strengths, not absolute absorption or
photoluminescence. Metadata specifies the algebraic circular-polarization
convention, input operator, band grouping, finite-intensity mask and source
validation. See the [PAW optical guide](PAW_OPTICS.md) for formulas and usage.

## User-defined reciprocal-space regions

The JSON has one top-level `regions` list. Each region has a unique `name`
and either a periodic circle (`center_fractional`, `radius_inv_A`) or an
explicit list of global `k_ids`. `total` and `rest` are reserved.

```json
{
  "regions": [
    {"name":"region_A", "center_fractional":[0.25,0.25,0.0], "radius_inv_A":0.15},
    {"name":"region_B", "center_fractional":[0.75,0.75,0.0], "radius_inv_A":0.15}
  ]
}
```

These are illustrative locations, not material-specific valley coordinates.
Circle radii are Cartesian reciprocal distances in Å⁻¹ with periodic images
restricted to the chosen plane. Centers must lie in that plane. Regions must
contain sampled points and be disjoint; `rest` is their complement. An explicit
ID region instead looks like `{"name":"patch","k_ids":[1,2,3]}`.

A regional contrast is a partition-dependent charge-response difference.
It is not automatically a conserved valley-current observable, a spin/layer
current or an integer Chern number.
