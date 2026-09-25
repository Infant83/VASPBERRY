# Native occupied-bundle Fukui reference

The native Fortran executable reads the actual 6×6 MnBi₂Te₄ VASP WAVECAR
and evaluates the complete occupied bands 1–123 (192 stored spinor bands).
It reports **Chern Number = −1.000000**. No Wannier model or Kubo integral
enters this result.

- [Original native curvature map](BERRYCURV.dat)
- [Native stdout](vaspberry.stdout.log)
- [All 36 plaquette comparisons](flux-comparison.csv)
- [Calculation metadata and checks](result.json)
- [Resource measurements](resources.log)

Every unique native plaquette agrees with the existing independent
[Python overlap reference](../fukui/fukui_occupied.csv) within the native
four-decimal map printing precision. The largest map flux difference is
3.32814e-06 rad; the six-decimal stdout comparison differs by at most
3.23692e-08 rad. Periodic display copies are removed before integration and
comparison; no curvature value is altered. The integral reconstructed from
the printed map is -1.000003324824.

The comparison retains the native Fortran reciprocal-area normalization
when reconstructing flux from curvature. The map's extrema are independently
checked against its numerical rows and recorded in `result.json`.

The recorded command used four MPI ranks and the current task-based CLI.
The serial and MPI builds select the same Fukui numerical routine:

```sh
mpirun -np 4 build/vaspberry-mpi --task chern --wavecar WAVECAR \
  --mesh 6,6 --spinor 2 --bands 1:123 --output BERRYCURV
```

The portable command, source identity and source-input association
are recorded in `result.json`. The MPI run completed in 371.88 seconds.
The launcher's resource measurement reported 15.41 MiB maximum
resident memory; this is not an independently measured aggregate MPI RSS.
These timings are host dependent. Task aliases require the current repository
version; they are not part of the immutable v1.3.0 source release.

The existing Python reference remains unchanged and supplies the additional
link singular-value and plane-wave-coverage diagnostics. The native agreement
checks implementation consistency for this lattice invariant; it does not
establish convergence of the coarse pointwise Hall integral or the physical
material model. Follow the [material workflow](../../) to regenerate the
licensed VASP source files locally.
