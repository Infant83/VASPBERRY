# MoS₂ bands and Kubo curvature along K–Γ–K′

This 49-point path supplements the 12 × 12 Brillouin-zone maps. The path
calculation uses the same structure, charge density, PAW datasets, spin–orbit
coupling, 400 eV cutoff and 26 bands as the mesh calculation. Only the k-point
list changes. The band energies and Kubo path curvature are calculated directly
by VASP and VASPBERRY; they are not interpolated from the mesh.

First complete the [full-mesh VASP calculation](../README.md). From the
repository root, prepare a separate working directory:

```bash
python3 examples/features/fukui-berry-curvature/prepare_path.py \
  --mesh-dir results/mos2-fullmesh-vasp \
  --output-dir results/mos2-path-vasp
```

The helper copies `INCAR`, `POSCAR`, `CHGCAR` and your licensed `POTCAR`
unchanged and supplies the path [`KPOINTS`](KPOINTS). Run the same licensed
noncollinear VASP executable used for the mesh:

```bash
cd results/mos2-path-vasp
/path/to/vasp_ncl > vasp.stdout.log 2> vasp.stderr.log
cd ../..
```

For an MPI installation, use the launcher and processor count configured for
that executable. Check electronic convergence and normal completion before
using the resulting `WAVECAR`.

## Path and expected output

| Point | Fractional reciprocal coordinates | Record |
|---|---|---|
| K | (1/3, 2/3, 0) | 1 |
| Γ | (0, 0, 0) | 25 |
| K′ | (−1/3, −2/3, 0) | 49 |

Each segment contains 25 points; Γ is included once. For the supplied lattice,
each segment is 1.320704 Å⁻¹ long. Plot energies relative to the highest occupied
band energy along this path; the CSV retains the original VASP energy zero.
The expected occupied–unoccupied gap is about 1.67355 eV. The demonstrated
single-process calculation took about 55 seconds with a peak resident memory
of 367 MB; costs depend on the VASP build and hardware.

The band dispersion is shared by the Fukui and Kubo figures. The Fukui curve
is extracted from the occupied-subspace plaquette field; its sampling is
described in the [Fukui tutorial](../../README.md). Individual-band Kubo
curvature is undefined at degeneracies, so near-degenerate points must be
masked as described in the [Kubo tutorial](../../../kubo-curvature/README.md).

The compact [reference outputs](../../reference/path/) contain the band data
and VASP convergence records. The licensed potential files, VASP executable
and generated `WAVECAR` are not distributed with this example.
