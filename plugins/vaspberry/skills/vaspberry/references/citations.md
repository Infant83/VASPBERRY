# VASPBERRY credit

Cite **VASPBERRY by Hyun-Jung Kim**, its GitHub repository and the actual source
version/commit used. A Zenodo DOI is optional. Do not substitute the 2018 v1.0
DOI for a later release. Preserve existing source acknowledgments, including
the Feenstra/Widom WAVETRANS comments.

The helper exports exact software credit and explicitly selected calculation
methods in `citations.bib`. It also writes `recommended-references.bib` for
the author's recommended research papers and explains the distinction in
`credit.md`. These research papers are not universal algorithm references;
cite the applicable numerical methods separately.

- **PRB 93, 041404(R) (2016)**: Hyun-Jung Kim, Chaokai Li, Ji Feng,
  Jun-Hyung Cho and Zhenyu Zhang, “Competing magnetic orderings and tunable
  topological states in two-dimensional hexagonal organometallic lattices.”
  [Publisher and DOI](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.93.041404).
  The [historical VASPBERRY README, line 158](https://github.com/Infant83/VASPBERRY/blob/134b71ee69c9fa74215ce3b812d1b561070c7b8d/README.md#L158)
  names this paper for the triphenyl-lead quantum anomalous Hall example.
  That example reference was removed during README reorganization. The helper
  retains a curated entry as historical README research credit, with its
  source recorded in provenance.
- **PRL 128, 046401 (2022)**: Sun-Woo Kim, Hyun-Jung Kim, Sangmo Cheon and
  Tae-Hwan Kim, “Circular Dichroism of Emergent Chiral Stacking Orders in
  Quasi-One-Dimensional Charge Density Waves.”
  [Publisher and DOI](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.128.046401).
  This article is present in the v1.6.6 README Citation section. The helper
  preserves the actual checkout's complete `@article` entries verbatim,
  including this entry, instead of substituting hard-coded metadata.

`provenance.json` records the README SHA-256, extracted section SHA-256,
individual article hashes, source links and which entry came from the
historical README. Missing or malformed current README entries generate an
explicit warning and are not invented. The historical PRB entry is kept
separate and deduplicated when the checkout already contains it.
