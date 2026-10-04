# Choose the input operator for a VASPBERRY response calculation

For charge Kubo calculations, begin with the same-run standard `WAVEDER`,
`WAVECAR`, `INCAR` and `OUTCAR` from a supported standard VASP optical calculation.
Version 1.6.5 selects this route by default. An explicit
`--kubo-source wavecar` selects the canonical-momentum approximation and
prints a warning. Missing or unsupported optical input never selects it
automatically. Other operator files enable separately validated observables.
VASP computes the electronic states; VASPBERRY evaluates their topology and
response. A separately prepared Wannier representation is retained for optional
supporting checks; it is not required by the native workflow.

## Available routes

| Route | Input and approximation | Present scope | VASP source modification |
|---|---|---|---|
| Native `--task kubo` (default) / `kubo-hall` | Same-run WAVEDER, WAVECAR, INCAR and OUTCAR; PAW longitudinal optical connections | Selected geometric band/bundle curvature and validated selected μ/T charge contributions; occupied T=0 compatibility | None; supported VASP 5.4.4 branch |
| Explicit `--kubo-source wavecar` | WAVECAR canonical momentum of pseudo-wavefunctions | Historical band/bundle curvature and reusable pairs; approximate μ/T charge scans and spin-sector proxies | None |
| `waveder-hall` | Same-run WAVEDER, WAVECAR, INCAR and OUTCAR; stored PAW longitudinal optical connections | Selected μ/T charge contributions with required-pair checks; occupied T=0 response and regional integrals | None; supported standard VASP 5.4.4 optical branch |
| `waveder-optics` | Same standard optical files; complete initial and final degenerate groups | k-resolved circular transition strengths and spectra; paths are allowed, with no inferred BZ integral | None; the same supported standard optical branch |
| Full velocity export / `velocity-pairs` → `pair-hall` | Audited full complex PAW velocity matrices, converted to reusable charge-response pairs | Charge and regional Hall response with chemical-potential and temperature scans | Optional instrumentation of the supported VASP 5.4.4 source copy |
| Full velocity and spin export / `spin-hall` | Audited full complex PAW velocity and spin matrices from the same eigenstates | Conventional spin response and a companion charge response for a gapped 2D occupied bundle at T=0 | Same optional producer |
| `wannier-hall` | VASP-derived Hamiltonian and position matrices | Dense full-connection integration of a finite Wannier model; gapped 2D occupied bundle at T=0 | Follow the source's VASP/Wannier interface recipe; the Hall postprocessor itself uses the supplied operators |

The [Kubo-formula response guide](KUBO_TRANSPORT.md) gives the full input contracts and commands.
The [output guide](OUTPUT_FORMAT.md) specifies the CSV, DAT, NPZ and metadata
formats. JSON stores configuration and provenance; it does not replace the
VASP calculation with an analytic model.

## 1. WAVECAR approximation: explicit opt-in

Follow the [MoS₂ charge-Hall example](../examples/features/kubo-hall/) to
prepare a full VASP mesh and run
`--task kubo-pairs --kubo-source wavecar` with the Fortran
executable, then `import-pairs` and `pair-hall` on the saved output. The
[explicit commands](KUBO_TRANSPORT.md#native-pairs-to-charge-hall) separate
wavefunction evaluation from postprocessing. The bundled Python stage applies occupations,
integrates the BZ and writes reusable numerical outputs. No custom VASP
producer or user-written integration script is needed. The same pair cache
can be reused for other chemical potentials, temperatures or region choices.

The canonical-momentum approximation does not include the full PAW,
nonlocal and SOC velocity terms. Increasing k sampling or the intermediate
band window improves the numerical evaluation of this approximation, but
does not supply the missing operator terms. Report that distinction with
the results. An integer Chern number, or cancellation of the total Hall
response by time reversal, is not by itself an accuracy test of local or
regional curvature magnitudes.

## 2. Standard charge-Kubo protocol: WAVEDER, no source patch

For a gapped T=0 occupied bundle, standard VASP `WAVEDER` already supplies
PAW longitudinal optical transition elements. The supported route uses
VASP 5.4.4 with `LOPTICS=.TRUE.`, `LPEAD=.FALSE.`, `LNABLA=.FALSE.`,
`NSW=0`, `ISYM=-1` and `LREAL=.FALSE.`. Preserve WAVEDER, WAVECAR, INCAR and
OUTCAR from that same completed run. Required non-canceling pairs must survive the producer's 2 meV treatment.
A complete producer-degenerate group may be selected if its internal response
weights are equal; an unresolved split or unequal weighting is rejected.

Use `kubo-hall --bands SELECTOR` for an occupation-weighted selected-band
contribution, with all source bands retained in the intermediate sum.
`--occupied N` remains the separate insulating T=0 compatibility mode; its
chemical potentials must lie inside the common gap. The two selectors are
mutually exclusive. The [complete command and supported settings](KUBO_TRANSPORT.md#standard-waveder-insulating-paw-hall-response)
include multi-run mesh assembly. The supplied [MnBi₂Te₄ optical example](../examples/materials/mnbi2te4-qah/)
illustrates this unmodified-VASP route.

For a PAW charge-Hall calculation within this scope, this is the default
route in 1.6.5. Native `--task kubo` writes point curvature; `kubo-hall`
provides Hall tables. Use `--input-dir DIR` to select the source directory;
omission means the invocation working directory. Each per-file option changes
only that file, so a WAVECAR override never redirects the optical inputs.
The [protocol](WAVEDER_KUBO_PROTOCOL.md) gives both commands and path rules.
The stored matrix may be rectangular. Every pair with unequal selected/
occupation weights must be stored in at least one orientation; equal-weight
internal pairs cancel. The adapter neither pads missing pairs nor truncates
the intermediate sum to the target set. The PAW terms
already enter VASP's matrix construction; the adapter does not add them again.
If the requested filling or operator lies outside the supported scope, stop
or choose a separately validated route. Do not silently substitute the
WAVECAR canonical-momentum approximation. File acceptance alone also does
not certify complete physical SOC/+U response terms.

This file stores an interband optical connection with an energy denominator
already applied. VASPBERRY does not divide by that denominator again. It
does not reconstruct the diagonal or degenerate velocity blocks that this
producer discards. Selected metallic or finite-temperature scans are
accepted only if their required matrix pairs and producer-cluster weights
are resolvable. Their result is a selected contribution, not automatically
the total AHC; the output records this distinction.

For circular optical selection, the ordinary native `-cd 2` WAVECAR workflow
also remains available. The optional [`waveder-optics` route](PAW_OPTICS.md)
uses the standard PAW optical elements and complete initial/final groups to output strengths in Å²
and point spectra in Å²/eV. It accepts paths as well as meshes, with no BZ
integration or absolute-absorption normalization. Its helicity selectivity
is an independent-particle transition property, not a PL polarization.

## 3. Optional full velocity and spin matrices

The [producer instructions](../tools/vasp544_spin_bridge/README.md) explain
how to instrument a separate, licensed VASP 5.4.4 source copy, build its
serial complex `ncl` executable and run the supplied wrapper. The original
insertion routines export the full velocity before the optical division,
including diagonal and degenerate blocks, together with physical PAW spin
and overlap matrices. This enables the spin-current calculation that
cannot be formed from occupied–empty optical transitions alone.

The procedure is:

1. Prepare and validate the ordinary SCF calculation for your material.
2. Follow the producer recipe in an isolated source tree; retain its manifest.
3. Prepare a fixed-charge SOC run using the same density, structure, PAW
   datasets, functional, cutoff and FFT settings. Choose accurate source
   bands and the required k sampling.
4. Run `run_vasp_spin_producer.py`, then `spin-export`, or assemble matching
   fixed-charge chunks with `spin-merge`.
5. For charge response, convert the full velocity matrix to a reusable pair
   cache as below, then use the ordinary `pair-hall` integration commands.
   For spin response, follow the [spin Hall commands](SPIN_HALL.md#calculate-and-reproduce)
   or [Bi material tutorial](../examples/materials/bi-spin-hall/).

With the actual `NX`, `NY` and energy reference chosen for your full mesh:

```bash
python3 tools/vaspberry_kubo.py velocity-pairs \
  --matrices spin-matrices/physical-matrices.npz \
  --metadata spin-matrices/physical-matrices.json \
  --mesh "$NX" "$NY" \
  --energy-reference 'unchanged VASP eigenvalue zero' \
  --output-dir results/full-velocity-pairs
```

The cache records the full-velocity operator provenance and source band
count. Pass it to `pair-hall --pairs-dir results/full-velocity-pairs` with
the same chemical-potential, temperature, retained-band and region arguments
as the WAVECAR comparison. This charge calculation uses velocity only;
it does not multiply charge curvature by spin expectation values. Supported
degeneracy rules and required convergence tests are the same as in the
[pair-Hall workflow](KUBO_TRANSPORT.md#wavecar-to-charge-hall-in-one-command).

The current audited producer has a specific supported Hamiltonian and build
contract; for example, it excludes hybrid/meta-GGA, Hubbard U and MPI.
Consult its guide before applying it to a new material. The distributed
numerical Bi matrices allow the VASPBERRY integration to be repeated without
building that producer. These extra steps are needed for generating this
full physical-operator input, not for the ordinary WAVECAR calculations.

The [generic matrix interface](KUBO_TRANSPORT.md#applying-the-workflow-to-other-data)
accepts other declared operators when they satisfy its schema and band-gap
rules. Its individual-band curvature command must not be applied to an
unresolved Kramers partner as if that partner were isolated.

## 4. Compare operators separately from convergence

Use the same geometry, density, k points, original eigenstates, source bands,
retained band window, occupations and region definitions whenever possible.
A side-by-side comparison then measures the effect of changing the operator.
The [MoS₂ matched-operator example](../examples/features/kubo-hall/operator-comparison/)
applies both charge workflows to the same final VASP wavefunctions.
For a nonmagnetic material, inspect nonzero regional or local responses as
well as the near-zero total. Keep complete degenerate groups.

Repeat k-mesh and source/retained-band tests for each operator. PAW and
nonlocal terms can improve the physical description, but do not guarantee
faster or monotonic k convergence. Sharp curvature peaks still require
adequate sampling. The source eigenstates and chosen electronic-structure
approximation must also be suitable for the material.

For a narrow-gap insulator, full-connection [Wannier interpolation](WANNIER_TRANSPORT.md)
can make dense integration affordable after independent checks against VASP
bands and operator-sensitive observables. Both Hamiltonian and position
matrices are required. A band fit alone does not validate the response.
