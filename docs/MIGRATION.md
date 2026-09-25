# Migrating to the 1.3.0 source

Version 1.3.0 is available as a versioned source release. Preserve the exact producer
version/commit and original files when migrating a calculation. The existing
Fukui and Z2 command-line interfaces remain available.

## Legacy Kubo factor of two

The older Fortran Kubo path formed circular matrix elements `Px+iPy` and
`Px-iPy` without the `1/sqrt(2)` normalization. Their squared difference is

```math
|P_x+iP_y|^2-|P_x-iP_y|^2=4\mathrm{Im}(P_xP_y^*).
```

Consequently, its negated numerator contained `-4 Im` instead of the standard
`-2 Im` for `A=i<u|grad_k u>`. The corrected 1.3.0 Kubo path uses the standard
factor. At otherwise identical inputs and numerical settings, the affected
legacy curvature magnitude is halved. This is a normalization correction,
not a guarantee of physical velocity completeness or mesh convergence.

| Data source | Normalization action |
|---|---|
| Confirmed affected older Fortran Kubo output | Import with the explicit historical doubled convention |
| Corrected 1.3.0 Fortran Kubo output | Import as standard normalization |
| Generic interband-matrix Kubo result using `-2 Im` | Already normalized; do not divide by two again |
| Fukui plaquette flux / Z2 n-field | Unaffected by this Kubo correction |
| Unknown or modified producer | Establish its formula from source/provenance before importing |

The `import-legacy` command requires an explicit source convention. It does
not infer a factor from a plausible Chern number, a filename, or a rounded
output. Keep the original data and the generated normalization metadata so a
later reader can tell which transformation was applied. Historical example
files remain historical and must not be silently overwritten or relabeled as
new calculations.

### Corrected Fortran output

The optional `-kubo_csv PATH` writes full-precision per-band rows with spin,
k index, band, fractional coordinates, energy, minimum gap and `omega_z_A2`.
Its comments identify `STANDARD_MINUS_TWO_IM`. Request all bands present in
the WAVECAR when preparing the current importer's input; a selected n-band
subset is not a complete import file. Use `-kubo 1` for a full 2D mesh and
`-kubo 2` only for point/path output. For example, after a full-grid calculation
has produced `kubo.csv`, a **scalar, spin-degenerate** 12×12 input can be imported
as follows:

```bash
python3 tools/vaspberry_kubo.py import-legacy \
  --csv kubo.csv --wavecar WAVECAR --normalization physical \
  --spin 1 --spinor-components 1 --spin-multiplicity 2 \
  --mesh 12 12 --plane-axes 0 1 \
  --energy-reference "Unchanged WAVECAR energy zero" \
  --degeneracy-threshold-eV 1e-6 \
  --output-dir imported-kubo
```

Replace the mesh and gap threshold with the actual calculation's choices.
SOC spinors instead use `--spinor-components 2 --spin-multiplicity 1`.
For explicit collinear spin channels, select `--spin` and use multiplicity 1;
each imported channel remains a charge-response contribution. The importer
checks coordinates and energies against the supplied WAVECAR.

### Historical doubled output

To prepare an older result, preserve its original files and make a CSV with
these columns:

```text
k_index,band,kx_frac,ky_frac,kz_frac,energy_eV,min_gap_eV,omega_legacy_A2
```

An optional `spin` column selects channels. There must be exactly one row per
source k point and band for the selected channel; the CSV must contain actual
matching energies and gap information. Old plot tables without these fields
are not self-contained transport input, and the importer does not guess the
missing values. Use `--normalization legacy-double` to apply 0.5 exactly once.
For an already normalized file the column is `omega_z_A2` and the flag is
`--normalization physical`. A standard-normalization producer marker prevents
accidentally importing a new Fortran CSV as historical doubled data.

The normalized output retains `normalization_migration`, including
`input_convention`, `factor_applied` and producer comments. Small-gap entries
are marked invalid; the Hall command rejects invalid band curvature. This is
not a replacement for the producer's degeneracy or physical-operator treatment.

The correction does not supply nonlocal, PAW, SOC or Hubbard-projector terms
missing from a canonical-momentum approximation. It also does not change a
line-mode calculation into a two-dimensional integration mesh. See the
[Kubo method guide](KUBO_TRANSPORT.md).

## Choose the correct output kind

- Keep Fukui flux in radians with its cell and vertex-energy interpretation.
  Its old output is not an input to the point-curvature importer.
- Use the standardized point-curvature NPZ/JSON pair for the new `hall`
  command. Arrays and metadata travel together; the JSON binds the NPZ bytes
  by checksum.
- Carry the actual energy reference and occupations. A selected-band Chern
  sum is not a replacement for a chemical-potential-dependent metal response.
- Preserve invalid or experimental status. Converting storage formats does
  not turn a failed matrix check into a validated physical calculation.

## Version and release records

`VERSION`, CLI version output and current citation metadata identify 1.3.0.
Historical changelog entries, the 1.2.0 release notes and archived reference
data retain their original versions. The CFF schema version is independent
of the software version.

The 2018 DOI `10.5281/zenodo.1402593` identifies VASPBERRY V1.0. It is not a
DOI for version 1.3.0. Use the immutable
[`v1.3.0` release](https://github.com/Infant83/VASPBERRY/releases/tag/v1.3.0)
for version-pinned source, or `master` for current source. Build binaries locally.
Cite the software version, exact commit and method references appropriate to
the calculation; see the [version policy](RELEASING.md).
