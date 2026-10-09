# Scientific validation and scope

Read this before a material calculation or a topological/transport claim. These are scientific decision criteria; the selected engine release's input contract and diagnostics determine the commands and supported routes. Version-specific conditions below describe the validated v1.6.6 release, not a promise about every future version.

## Terms and quantities

| Term | Meaning in this workflow |
|---|---|
| VASP and its files | Vienna Ab initio Simulation Package, the separate electronic-structure producer. WAVECAR stores wavefunction coefficients/eigenstate data; WAVEDER carries optical derivative matrix data for a supported producer. INCAR specifies calculation controls; OUTCAR records producer settings/results. Retain the actual associated files and provenance. |
| Brillouin zone (BZ) | Reciprocal-space primitive cell with periodic boundaries; k is the Bloch wavevector (crystal momentum ℏk, with ℏ the reduced Planck constant), and Gamma (Γ) is k=0. A high-symmetry band path samples lines, not the area needed for a 2D Hall integral. |
| Berry curvature, Ω | Local geometric curvature of Bloch states. A Cartesian k-space component has units of length squared; reduced-coordinate components require their own Jacobian. |
| Plaquette flux, φ | Berry phase integrated around one mesh cell, dimensionless and reported in radians. Dividing by the physical cell area can give a finite-cell average, not an independently evaluated point curvature. |
| Chern number, C | Integral of curvature divided by 2π over a closed 2D BZ for an isolated band/subspace. Occupied-subspace C and a selected-band C answer different questions. |
| FHS | Fukui–Hatsugai–Suzuki discrete link-variable method for Chern numbers. Integer-valued output on a grid does not alone establish that the continuum material phase was resolved. |
| FH / lattice n-field | Fukui–Hatsugai Z2 construction: parity of a lattice integer field over a half BZ, with the required time-reversal gauge. The local n-field is gauge dependent; the qualified invariant is not. This is distinct from FHS Chern evaluation. |
| Time-reversal symmetry (TRS) | Symmetry under the physical time-reversal operator. Spinful time reversal squares to −1 and implies Kramers pairing. Zero reported magnetic moment alone does not prove TRS. |
| Band bundle / occupied subspace | A band bundle is a fixed-rank selected span of Bloch states; it need not be occupied. The physical occupied subspace contains all occupied states. At internal degeneracies, use a supported multi-band subspace rather than assigning an isolated-band invariant. |
| Direct / global insulating gap | Direct: separation at a given k point. Sampled global: minimum conduction energy minus maximum valence energy over the sampled BZ. Positive direct gaps can coexist with a negative global gap. |
| Link conditioning / branch margin | Quality of neighboring-state overlap links (or subspace overlaps) and distance of cell phases from the principal-log branch boundary. Poor links or near-π cells can indicate an unresolved grid or invalid subspace. |
| Intrinsic charge Hall conductivity / AHC | Charge-current response calculated from the specified velocity operators and band occupations; AHC means anomalous Hall conductivity. The intrinsic Berry/Kubo contribution does not include every extrinsic transport mechanism. |
| Kubo route | Linear-response evaluation using the supported velocity/current matrix elements, interband energy denominators and occupation weights. Its input/operator approximation and degeneracy handling must be specified. |
| Chemical potential μ / temperature T | Parameters of the electronic occupations. A fixed-band μ scan is a rigid-band analysis, not a self-consistent doping calculation; changing occupation T is not a new finite-temperature structure calculation. |
| SOC / PAW | Spin–orbit coupling / projector augmented-wave method. WAVECAR pseudo-wavefunction coefficients alone do not reconstruct all PAW overlap or velocity augmentation. |
| Spin Chern / conventional SHC | Projected-spin-sector topology / spin Hall conductivity defined with a specified spin-current operator, conventionally jᵢˢ={vᵢ,s}/2. Here vᵢ is velocity along axis i, s is the chosen spin angular-momentum component and braces denote the anticommutator. They are distinct observables, especially when spin is not conserved. |
| Valley contribution | Integral over an explicitly defined reciprocal-space region. Its value depends on the region and weighting; it is not automatically an invariant or the total Hall response. |

## Establish the calculation before running

Identify the observable, material/model, dimensionality or 2D slice, magnetic/SOC state, occupied/selected bands, energy zero, BZ geometry and intended claim. Check that the supplied artifacts support that route. A band plot or PROCAR (orbital/spin projection data) is not a replacement for wavefunctions or required velocity/spin-current matrix elements.

For VASP data, retain the input/control files and output metadata that establish a common producer calculation: structure/lattice, functional, PAW dataset identities, cutoff, charge-density workflow, VASP version/build, spin/SOC settings, mesh, eigenvalues, band count and occupation. Reconcile dimensions, reciprocal coordinates and source hashes before combining files. Follow the release's exact same-run requirements; matching filenames are insufficient. Do not distribute licensed PAW data to create a provenance record; record identifiers and permitted hashes.

Treat a dirty engine checkout explicitly as modified: keep the source patch or a complete reproducible source snapshot, including relevant untracked files. Helper hashes cover selected files, not the whole dependency graph. For native work also record executable hash, compiler/build flags, relevant BLAS/LAPACK (linear-algebra libraries), MPI launcher and rank count. If a source commit cannot be recovered, identify the archive and its hash rather than inventing a commit.

## Route-specific conditions

| Route | Evidence required before interpretation |
|---|---|
| FHS Chern | Full closed periodic 2D grid; a fixed isolated band or isolated subspace; valid links and branch diagnostics; correct reciprocal orientation. Check separation from excluded bands across the sampled BZ. Degeneracies inside a supported subspace do not imply separation from its complement. Direct WAVECAR overlaps use the pseudo-wavefunction approximation without PAW augmentation. |
| FH Z2 | Physically TR-symmetric, nonmagnetic 2D spinor insulator with the complete occupied even-rank subspace, not an arbitrary subset. For v1.6.6: unshifted Gamma-centered Nx×Ny×1 grid, even Nx/Ny ≥4, kz=0, ISYM=−1, ISPIN=1, LSORBIT enabled, occupied bands 1:Nocc and NBANDS>Nocc. Retain the first unoccupied sentinel band and verify both direct isolation and a positive sampled global insulating gap. |
| Standard charge Kubo/WAVEDER | Use the supported optical producer and complete required matrix-pair coverage from the same run; v1.6.6 documents the standard VASP 5.4.4 longitudinal optical branch. Check the selected/partner band windows, occupation weights and exclusions. Producer provenance determines which PAW/velocity terms are present; do not add them twice or declare arbitrary WAVEDER files equivalent. |
| WAVECAR Kubo approximation | Requires explicit informed choice. Canonical momentum from pseudo-wavefunctions omits the documented PAW augmentation/nonlocal velocity terms; it cannot be silently substituted for a missing standard optical input. |
| Projected-spin sectors | Specify the spin axis/frame, occupied subspace and projected-spin spectral gap. Keep the pseudo-overlap/operator approximation labels and inspect the sector/band gaps. A spin-sector integral is not conventional spin-current SHC. |
| Conventional spin Hall / Wannier | Verify the supported exported operator route, units and spin-current definition in the release's operator guide. In v1.6.6 the conventional spin-Hall route is an insulating 2D, T=0 response using validated same-state spin/full-velocity matrices; do not infer unsupported metallic, finite-T, layer-current or 3D results. PROCAR projections alone are insufficient. The v1.6.6 `wannier-hall` route also assumes a fixed insulating occupied bundle in 2D at T=0. A full VASP-derived Wannier connection requires the documented Hamiltonian and position-matrix inputs; a band fit alone does not validate curvature/operators. |

In the FH grid condition, Nx and Ny are in-plane grid dimensions, Nocc is the number of occupied spinor states at each k, and NBANDS is the producer's total band count. The VASP controls ISYM, ISPIN and LSORBIT specify symmetry handling, spin mode and SOC, respectively; ISYM=−1 disables symmetry handling. These input tags do not themselves prove the physical TRS of the states.

For FH Z2, inspect the final native `Z2_FIELD.csv` validity fields, not only the printed parity or NFIELD diagnostic. In v1.6.6 require `result_status=PASS` and `reportable_invariant=1`, agreement of the complementary half-zone parities and passing associated diagnostics; an `.invalid.csv` output is not reportable. An internal TR reconstruction can be consistent even when raw input has not been shown to satisfy physical TRS. A PASS cannot supply that missing evidence. An occupation override cannot turn a metal into an insulator. This plugin does not infer a 3D strong/weak Z2 classification from its supported 2D FH route, or a Z2 index from zero Chern number. A local n-field map is not a physical Berry-curvature map.

For a valley analysis, record periodic region centers, shape/mask and radius, any normalization, the complementary/rest contribution and sensitivity to the region definition. Do not introduce an automatic factor of one half in a K−K′ diagnostic. Orbital/layer character attribution is not a corresponding current operator. Keep projected-spin-sector Chern numbers, any defined spin-Chern combination, spin-sector Kubo contributions and conventional spin-current conductivity separately labeled.

Distinguish geometric band weights from physical occupation/spin multiplicity. In particular, do not multiply explicitly represented SOC spinors by an extra factor of two. Use the release's selected/occupied weighting contract and inspect recorded metadata.

## Numerical and physical validation

Keep four separate conclusions: command execution; method-specific numerical checks; observed convergence under the studied approximations; physical interpretation. Record PASS, FAIL or not evaluated for each applicable check. Describe the scope of PASS rather than treating it as a universal material certificate.

- Refine the relevant BZ mesh and compare the target response, gap and diagnostics over a systematic sequence. A stable integer on a coarse mesh is insufficient if narrow avoided crossings or Berry hotspots remain unresolved. Report mesh dimensions/offsets, results, changes, declared tolerance and which refinements were actually run. Additional values cannot be manufactured from interpolation of a plot.
- For Kubo response, converge the required partner/intermediate-band space and the occupied or selected contribution scope. Intermediate states may lie below or above a selected band; do not restrict them to unoccupied states unless the route explicitly requires that restriction. Examine sensitivity to any broadening/regularization and to actual degeneracies. A degeneracy-rejection tolerance or optical-producer clustering threshold is not a physical lifetime broadening. Where the route does not implement lifetime broadening, do not invent it. Do not drop missing matrix pairs or relax input/validity guards merely to obtain a requested number.
- Check the parent electronic-structure convergence relevant to the claim (functional, cutoff, source mesh, structure/magnetism, SOC and adequate empty states). A denser postprocessing mesh does not replace source-state convergence. If new VASP inputs are needed, say what is missing; the plugin does not run VASP.
- Use appropriate analytic/reference cases and limiting occupations where available. Empty-band zero response and insulating occupied-band limits are useful numerical checks. A known-result test of shared production routines does not independently validate the WAVECAR parser, PAW reconstruction or real-material applicability.
- Compare an independent method/backend when a disputed phase or approximation-sensitive publication claim needs it, after matching subspaces, units and conventions. State when that comparison was not performed or is outside the supported inputs.

There is no universal mesh size or tolerance that certifies every material. Adopt the release's validity guards and a declared study-specific accuracy target. When only user-reported scalar values are available, label them as unverified reported values and explain the additional evidence needed. C=0.998 or a near-quantized coarse Hall value alone does not establish a quantum anomalous Hall (QAH) phase.

## Units, signs and figure interpretation

Record the Berry connection convention, reciprocal-plane orientation and Hall tensor/current convention. In the v1.6.6 public Qi–Wu–Zhang (QWZ) model, A=i⟨u|∇k u⟩ with positive b1×b2 orientation gives the filled lower-band relation σxy=−C e²/h at T=0. Here u is the periodic part of a Bloch state, b1 and b2 are the ordered in-plane reciprocal basis vectors, and σxy is the specified xy Hall tensor component. Use the actual output convention rather than changing a sign to fit a remembered phase diagram. Define e as the positive elementary-charge magnitude and h as Planck's constant when writing this relation.

Distinguish 2D sheet conductivity in e²/h or siemens from 3D bulk conductivity in S/m. Do not convert a slab response using the simulation's vacuum length without a stated physical thickness convention. Distinguish σ from reference-subtracted Δσ. A selected-band or valley table must remain labeled as a contribution unless the full physical occupied response and mesh coverage have been established.

For each plot identify the quantity, tensor/spin component, units, k coordinates or μ energy reference, temperature, mesh/band/region scope and approximation. Do not label plaquette flux as point curvature; disclose area normalization and interpolation if applied. A path curvature plot is a path diagnostic, not a full-BZ integral.

## Primary references and release contracts

- FHS Chern method: [Fukui, Hatsugai and Suzuki, JPSJ 74, 1674–1677 (2005)](https://doi.org/10.1143/JPSJ.74.1674).
- FH Z2 method: [Fukui and Hatsugai, JPSJ 76, 053702 (2007)](https://doi.org/10.1143/JPSJ.76.053702). The paper's scope does not extend this plugin's 2D implementation to an unimplemented 3D workflow.
- Berry/transport background: [Xiao, Chang and Niu, Rev. Mod. Phys. 82, 1959 (2010)](https://doi.org/10.1103/RevModPhys.82.1959).
- The v1.6.6 implementation contracts are `docs/Z2_FUKUI_HATSUGAI.md`, `docs/VALLEY_TRANSPORT.md`, `docs/WAVEDER_KUBO_PROTOCOL.md`, `docs/SPIN_CHERN.md`, `docs/SPIN_KUBO.md`, `docs/SPIN_HALL.md`, `docs/WANNIER_TRANSPORT.md`, `docs/OPERATOR_ROUTES.md` and `docs/OUTPUT_FORMAT.md` in the exact selected engine checkout. Read the matching route; terminology here does not replace its input schema or guards.
