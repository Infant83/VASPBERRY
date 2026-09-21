# Point-curvature and Hall output formats

Version 1.3.0 introduces `vaspberry.band-curvature` schema version **1** and
`vaspberry.hall-spectrum` schema version **1**. Schema versions describe file
contracts and are independent of the software version. These formats do not
replace the existing Fukui plaquette or `VASPBERRY_Z2_FIELD` formats.

## Native Kubo bundle CSV

Native `-kubo_bundle 1 -kubo_csv PATH` writes
`VASPBERRY_BARE_MOMENTUM_KUBO_BUNDLE_V1`. It contains one row per spin channel
and source k point, with comment metadata followed by a CSV header. The
file stores the trace curvature of bands `-ii` through `-if`, with unit
occupation. It has no per-band energy column or legacy `.dat` companion.

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
| `result_status` | `PASS` after the spin/k isolation check |
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

`hall` writes `conductivity.csv`, `conductivity.npz`, and `conductivity.json`.
The CSV and NPZ use the same long-form rows, one per chemical potential,
temperature, region and band channel. JSON includes hashes for both files,
the source curvature metadata, reciprocal area/normal, region definitions,
occupation-window checks and algorithm settings.

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

The default integration checks require a complete represented source-band
window and negligible occupation of its highest band over the requested μ/T
range. `--allow-partial-bands` labels the result `partial_band_contribution`;
it does not add omitted occupied bands. A successful window check still does
not prove that the source calculation includes enough unoccupied states for
curvature convergence.

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
