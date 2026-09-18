# Developer checks and supplementary references

These analytic models, synthetic WAVECAR fixtures and stored-field checks are
for testing numerical routines. The [user tutorials](../../examples/README.md)
start from actual public VASP calculations and reproduce VASPBERRY outputs.

These references are also retained as supplementary teaching and technical
report material. The [reference map](../../docs/REFERENCE_MATERIALS.md) pairs
each check with an actual VASP exercise and explains how to reuse its results
in a hands-on guide or release report.

The six former feature checks are retained here, including their historical
reference files and source hashes. Their JSON files configure generated
fixtures or validation; they are not VASP calculation inputs. Historical
provenance retains the original file locations at its recorded source commit.

```bash
python3 validation/models/run_checks.py --list
make serial
python3 validation/models/run_checks.py --all --output-dir results/developer-checks
```

| Check | Purpose |
|---|---|
| [Fukui/Chern](fukui-chern/) | Known QWZ phases and flux/Hall sign |
| [Matrix Kubo](kubo-curvature/) | Independent analytic curvature oracle |
| [Hall/regions](hall-valley/) | Occupation and arbitrary-region sum rules |
| [Stored Z2 field](z2/) | Historical schema-2 consistency and plotting |
| [Optical](circular-dichroism/) | Synthetic scalar angular-selectivity oracle |
| [Wavefunction](wavefunction/) | Synthetic scalar Gamma amplitude/bounds regression |

Identify analytic models and synthetic inputs explicitly when reusing results.
The stored Z2 check validates an existing material result but does not itself
rerun the wavefunctions. Public user examples do not substitute these checks
when an input WAVECAR is missing.
