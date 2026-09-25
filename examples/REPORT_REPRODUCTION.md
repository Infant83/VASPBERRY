# Reproduce the technical report

The [technical report](../docs/TECHNICAL_REPORT.md) uses the numerical
examples below. **Redrawing a figure** reads distributed numerical output;
**recalculating it** starts from WAVECAR or the stated operator matrices.
The table keeps those two operations separate. VASP and licensed POTCAR
files are needed only when regenerating electronic-structure inputs that
are not distributed. No main native workflow requires Wannierization.

For the first complete calculation, follow the [hands-on guide](../docs/HANDS_ON.md).
The [native command reference](../docs/NATIVE_COMMANDS.md) groups the task
selectors and their arguments for reuse with your own material.

Run commands from the repository root, after `make serial` and installing
`requirements-transport.txt`. Each linked guide gives the complete command
sequence and the expected result. Use a fresh output directory for a new
calculation. Keep mesh and path WAVECARs in separate directories.

| Report result | Calculation and exact commands | Numerical files used by the figure | What is available without a new VASP calculation? |
|---|---|---|---|
| §3.1 / Fig. 1: MoS₂ Fukui map | [Fukui guide](features/fukui-berry-curvature/) — `--task chern --mesh 12,12 --bands 1:18` | `features/fukui-berry-curvature/reference/BERRYCURV.dat` and `reference/path/bands.csv` | Redraw map, bands and cut with `tools/plot_berry_panels.py`. Regenerate the full mesh and matching 49-point path to recalculate. |
| §3.2 / Fig. 2: MoS₂ occupied Kubo bundle | [Kubo guide](features/kubo-curvature/) — `--task kubo --bundle 1 --bands 1:18` on both WAVECARs | `features/kubo-curvature/reference/bundle/{KUBO_mesh.csv,KUBO_path.csv,bands.csv}` | Redraw all panels. Recalculation uses the same regenerated mesh/path as Fig. 1. |
| §3.2.1 / Fig. 3: isolated valley band | [Valley guide](features/kubo-curvature/valleys/) — 162-point VASP patches, native bands 17–18 | `features/kubo-curvature/valleys/reference/{KUBO.csv,summary.csv,bands.csv}` | `valleys/plot.py` redraws the reference. Recalculate after generating the two 9×9 patches. |
| §3.2.2 / Fig. 4: μ/T-dependent Hall response | [Hall guide](features/kubo-hall/) — native `kubo-pairs`, `import-pairs`, `pair-hall`, `plot_hall.py` | [Eight mesh/cutoff case tables](features/kubo-hall/reference/) | Redraw both composite figures. Regenerate each full WAVECAR for fresh pairs; a single 24×24 run does not reproduce the complete convergence study. |
| §3.2.3 / Fig. 5: optional matched operators | [Operator comparison](features/kubo-hall/operator-comparison/) | `reference/canonical-pairs`, `reference/paw-pairs`, and four conductivity tables | Both supplied pair caches support new integrations and plots. Generating new full-velocity operators requires the separate supported producer. |
| §3.3 / Fig. 6: stacking optical selection | [Stacking guide](materials/mos2-stacking-valley/) — fresh SCF and path for each of four geometries | Each case's `bands.csv`, `native-optical.csv`, optional `paw-optical/` | `plot.py` redraws all panels and channel spectra. Recalculate transitions after four ordinary VASP path runs; 3R additionally needs the documented dipole SCF restart. |
| §3.4 / Fig. 7: MoS₂/Bi Z₂ comparison | [Paired guide](features/z2/comparison/), [MoS₂](features/z2/mos2/), [fresh Bi](materials/bi-spin-hall/#3-fukui-z₂-and-the-n-field) | `features/z2/mos2/reference/Z2_FIELD.csv` and `materials/bi-spin-hall/reference/z2/Z2_FIELD.csv` | `features/z2/compare.py` redraws both fields. Fresh Bi has an ordinary `--stage wavecar` preparation; PAW spin matrices are optional. The public historical Bi WAVECAR is a different source. |
| §3.5 / Fig. 8: Γ-state density | [Wavefunction guide](features/wavefunction/) — native `wavefunction`, then `--postprocess-only` | `reference/raw-output.tar.gz`, `summary.csv`; bundled 48-point MoS₂ WAVECAR | Complete native calculation, Fourier check and figure reproduction without VASP. |
| §3.6: MnBi₂Te₄ Chern number | [Magnetic-film guide](materials/mnbi2te4-qah/) — native `chern --mesh 6,6 --bands 1:123` | `reference/native-fukui/BERRYCURV.dat`, stdout and plaquette comparison | Inspect/replot occupied flux. Recalculate after ordinary SOC+U VASP from the supplied SCF density. The separate coarse WAVEDER integral is explicitly unconverged. |
| Appendix A.3 / Fig. 9: optional Bi spin Hall | [Bi guide](materials/bi-spin-hall/) — `restore_matrices.py`, `spin-hall` | Public 6×6 and 12×12 physical matrices; conductivity and convergence tables | Recalculate the supplied matrix cases and redraw the whole convergence figure. New mesh/source calculations need licensed VASP and the supported PAW producer. |
| Appendix B.2 / Fig. 10: optional Bi ideal edges | [Bi guide](materials/bi-spin-hall/#4-confirm-the-bulk-z₂-result-with-an-ideal-edge) — `wannier-bands`, `wannier-edge` | Public Hamiltonian/position operators; bulk/edge reference arrays | Recalculate the finite model and redraw figures. This supporting model is separate from native Z₂ and PAW spin Hall. |
| Appendix B.3 / Fig. 11: optional MnBi₂Te₄ finite model | [Native Wannier guide](materials/mnbi2te4-qah/NATIVE_WANNIER.md) | Public archived operators; five quadrature cases and direct-DFT sample bands | Recalculate the fixed model and redraw figures. It does not establish convergence of the coarse direct-VASP optical integral. |

Sections 1.1–1.3 define the Fukui, Kubo and output contracts used by the
main examples. Section 2 identifies their material settings. Appendix A.1
is exercised by the stacking WAVEDER spectra and the MnBi₂Te₄ coarse optical
integral; A.2 by the Bi physical spin matrices; B.1 by the two supplied
Wannier-model integrations. Section 4 explains reuse and convergence.

The original [circular-optics tutorial](features/circular-dichroism/) also
recalculates directly from the bundled 48-point MoS₂ WAVECAR. The
[archived Bi Chern](features/fukui-chern/), [Z₂](features/z2/) and
[gap-Hall](features/hall-valley/) tutorials use the separately downloadable
historical Bi WAVECAR. They are introductory checks, not the fresh Bi input
used in report Fig. 7. `examples/run_examples.py --all` runs the eight
catalogue tutorials; it does not run every material study or reproduce all
report figures automatically.

The additional [PROCAR character example](features/procar-character/) gives
a small analytic fixture for the layer/orbital/spin-attribution commands in
report §1.4. It supplements the report's workflow discussion; it is not an additional
VASP material result or one of Figures 1–11. Transfer its `project`, `hall`
and `plot` stages to matching material files after defining the atom groups
and choosing isolated bands and a spin axis.

For another material, retain numerical outputs and their units, band
selection, operator identity and energy zero. The generic native CLI,
`tools/vaspberry_kubo.py` and plotters accept new inputs; material `run.py`
helpers deliberately check their particular reference settings. CSV/DAT
can be plotted in Python, gnuplot, Origin or a spreadsheet; NPZ preserves
multidimensional arrays. Keep original samples for physical integration:
display interpolation changes the figure, not k-mesh convergence.
