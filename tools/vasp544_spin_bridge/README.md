# PAW spin and velocity producer

This opt-in instrumentation writes the **full band matrices** needed for a
conventional intrinsic spin Hall calculation. The Python instrumenter and the
two Fortran insertion files are original VASPBERRY code. They do not contain a
VASP implementation module, executable, or PAW potential.

You need a licensed VASP source installation to build the producer. The
instrumenter supports the audited VASP 5.4.4 source revision and refuses other
or previously modified versions. Use a separate copy of your source tree.
Your site's VASP licence and build requirements still apply.

## Build and run

From the VASPBERRY repository root:

```bash
python tools/vasp544_spin_instrument.py /path/to/isolated-vasp544
```

Configure that copy for a **serial, complex `ncl` build**, using your site's
compiler and BLAS/LAPACK setup, then build it:

```bash
make -j1 -C /path/to/isolated-vasp544 ncl
```

Use a serial build for the first instrumentation rebuild. The response module
must be regenerated before its optical consumer; older VASP make dependencies
can leave stale module files or compile them out of order under parallel make.

The instrumenter records the supported revision and all modifications in
`vaspberry-spin-producer.json`. Both output flags default to false, preserving
the ordinary calculation path. This producer does not support MPI or the
Gamma-only executable.

Prepare a new, fixed-charge SOC calculation with a converged `CHGCAR` and the
same structure, potentials, functional, cutoff and FFT-grid settings as its
SCF source. Include the following settings:

```text
ICHARG = 11
ISPIN = 1
LSORBIT = .TRUE.
SAXIS = 0 0 1
ISYM = -1
LREAL = .FALSE.
NSW = 0
LWAVE = .TRUE.
LCHARG = .FALSE.
LOPTICS = .TRUE.
LPEAD = .FALSE.
LNABLA = .FALSE.
LBERRY_EXPORT = .TRUE.
LSPIN_EXPORT = .TRUE.
```

Choose the k mesh, band count and electronic convergence settings for your
system. Empty-state accuracy matters even for an insulating occupied bundle.
Hybrid and meta-GGA functionals, Hubbard U, spin-dependent PAW overlap, non-Cartesian spin
axes, finite-q response and spin spirals are outside the supported contract.

Run with the supplied wrapper so input/output association and completion are
recorded. It uses one scientific thread, a one-hour default time limit and a
4 GiB default memory limit; limits are explicit options.

```bash
python tools/run_vasp_spin_producer.py \
  --run-dir /path/to/new-nscf \
  --binary /path/to/isolated-vasp544/bin/vasp_ncl \
  --producer-manifest /path/to/isolated-vasp544/vaspberry-spin-producer.json

python tools/vaspberry_kubo.py spin-export \
  --run-dir /path/to/new-nscf --output-dir spin-matrices
```

The converter requires a completed, converged run and validates its recorded
files before writing a new output directory. A zero process exit code alone
does not count as success. Outputs are `physical-matrices.npz/.json`,
`augmentation.npz/.json`, and `producer-audit.json`. The physical matrices feed
the [`spin-hall` workflow](../../docs/SPIN_HALL.md). Independent
fixed-charge k-point chunks can be combined with `spin-merge`; their union must
cover the declared full uniform mesh with no duplicates.

## Operators and conventions

Arrays retain every stored band and the original eigenvector gauge.
`spin_pauli[k,a,n,m]` is the dimensionless matrix of Pauli component `a`;
physical spin is `hbar/2` times this matrix. `velocity_eVA[k,a,n,m]` is `hbar*v`
in eV Å, with Cartesian component `a` and bra/ket indices `n,m`.

The spin matrix includes both the smooth plane-wave contraction and the PAW
on-site overlap correction. Their sum is checked against the physical overlap
metric. Neither the pseudo coefficients nor the resulting matrices are
renormalized or symmetrized.

For velocity, the instrumentation preserves the commutator **before** the
optical routine divides by an energy difference or clears a degenerate block.
It restores the removed energy derivative times the full PAW overlap and adds
the projector-derivative and one-center dipole contribution. In the producer's
derivative-with-respect-to-`ik` convention,

```text
D_mn = i [Xi_mn + eprime_n O_mn + (E_n - E_m) T_mn].
```

The stream stores the first two terms together as `xi_restored`. This formula
provides diagonal and degenerate blocks directly, with no replacement by zero.
The converter compares separated pairs with the optical connection and checks
the independently exported diagonal energy derivative. These are consistency
checks within the producer. Finite-k energy slopes and independent response
benchmarks provide additional physical validation.

### Optical comparison near degeneracies

With `LSPIN_EXPORT=.TRUE.`, the optical comparison uses the standard **0.002 eV**
cluster threshold. The full velocity is captured before that optical division
and retains its diagonal and all degenerate blocks. The cluster threshold
therefore controls which optical elements can be used as a consistency check;
it does not truncate the exported full velocity or change its eigenstates.

This matters when numerical splittings are very small: reconstructing a
derivative state with a tiny energy denominator can amplify contraction and
orthogonality roundoff into other optical entries. The comparison still
requires an absolute residual below `1e-7` eV Å on every covered separated
element, along with the unchanged diagonal, commutator and Hermiticity tests.
There is no tolerance relaxation or matrix repair.

Legacy `LBERRY_EXPORT`-only optical exports retain their `1e-10` eV threshold.
Both policies remain readable, and the stream header, build manifest and
producer audit record the actual policy. Original Bi reference files keep
their recorded legacy producer settings. When rebuilding the producer, check
the full velocity against the previous build on the same fixed eigenstates.

The spin Hall backend currently constructs the anticommutator from these
**finite band matrices**. This is a projected-product approximation: excluded
states can contribute to the spin-current matrix. Increase the source band
count and check convergence separately from the k mesh. Passing the converter
does not establish material or response convergence, nor does a nontrivial Z2
invariant require the conventional spin Hall conductivity to be quantized.

## Reproduction boundary

The insertion points are in the optical and linear-response modules: before
optical evaluation, immediately after the raw commutator is assembled, after
the PAW projector/dipole correction is available, and after the final optical
matrix is computed. The instrumenter validates the complete source files and
all insertion points before writing either file. The two insertion files are
readable, original contraction/export routines. They require the locally
licensed VASP data types and numerical routines at compilation time.

For users without a VASP source licence, published numerical matrix bundles
allow the complete VASPBERRY post-processing stage to be reproduced. They do
not substitute for generating a new material's PAW operators.
