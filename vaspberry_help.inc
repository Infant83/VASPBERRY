! Native help is processed before the calculation parser or input I/O.
! This include is shared by current serial/MPI and historical serial builds.
      subroutine vaspberry_help_request(ver_tag,legacy_build)
      implicit none
      character*(*) ver_tag
      logical legacy_build
      integer narg,j,topic_length,ios
      character*256 arg,topic
      narg=command_argument_count()
      if(narg.eq.0)return
! Only option slots are inspected. A filename/value equal to --help is data.
      do j=1,narg,2
       call get_command_argument(j,arg)
       if(trim(arg).eq.'--help'.or.trim(arg).eq.'-h')then
        if(j.ne.1.or.narg.gt.2)then
         call vaspberry_help_error(
     &    'use --help [TOPIC] alone; no calculation options or extra arguments')
        endif
        topic='overview'
        if(narg.eq.2)then
         call get_command_argument(2,topic,length=topic_length,
     &                             status=ios)
         if(topic_length.gt.len(topic).or.ios.ne.0)then
          call vaspberry_help_error('help topic is too long or unreadable')
         endif
         if(len_trim(topic).eq.0)then
          call vaspberry_help_error('help topic must not be empty')
         endif
        endif
        call vaspberry_help_show(ver_tag,topic,legacy_build)
       endif
      enddo
! Reject removed options before any input I/O, including historical builds.
      do j=1,narg,2
       call get_command_argument(j,arg)
       if(trim(arg).eq.'--bundle'.or.
     &    trim(arg).eq.'-kubo_bundle')then
        call vaspberry_help_error('--bundle/-kubo_bundle was removed; '//
     &   'use --task kubo --bands FIRST:LAST for the selected-space trace, '//
     &   'or add --per-band 1 for individual bands')
       endif
      enddo
      end subroutine vaspberry_help_request

      subroutine vaspberry_help_rank(print_help)
      implicit none
      logical print_help
#ifdef MPI_USE
      include 'mpif.h'
      integer ierr,help_rank
      logical initialized,finalized
#endif
      print_help=.true.
#ifdef MPI_USE
      call MPI_INITIALIZED(initialized,ierr)
      if(initialized)then
       call MPI_FINALIZED(finalized,ierr)
       if(.not.finalized)then
        call MPI_COMM_RANK(MPI_COMM_WORLD,help_rank,ierr)
        print_help=help_rank.eq.0
       endif
      endif
#endif
      end subroutine vaspberry_help_rank

      subroutine vaspberry_help_finish(exit_status)
      implicit none
      integer exit_status
#ifdef MPI_USE
      include 'mpif.h'
      integer ierr
      logical initialized,finalized
      call MPI_INITIALIZED(initialized,ierr)
      if(initialized)then
       call MPI_FINALIZED(finalized,ierr)
       if(.not.finalized)call MPI_FINALIZE(ierr)
      endif
#endif
      if(exit_status.ne.0)stop 1
      stop
      end subroutine vaspberry_help_finish

      subroutine vaspberry_help_error(message)
      implicit none
      character*(*) message
      logical print_help
      call vaspberry_help_rank(print_help)
      if(print_help)then
       write(0,'(A)')'*** error - '//trim(message)
       write(0,'(A)')'Use --help task, --help options or --help all.'
      endif
      call vaspberry_help_finish(1)
      end subroutine vaspberry_help_error

      subroutine vaspberry_help_show(ver_tag,requested,legacy_build)
      implicit none
      character*(*) ver_tag,requested
      character*256 topic
      logical legacy_build,print_help,known
      topic=trim(requested)
      if(topic(1:2).eq.'--')topic=topic(3:)
      select case(trim(topic))
      case('tasks')
       topic='task'
      case('legacy')
       topic='all'
      case('kubo-line')
       topic='kubo'
      case('-f')
       topic='wavecar'
      case('-s')
       topic='spinor'
      case('-kx')
       topic='mesh'
      case('-ky')
       topic='mesh'
      case('-ii')
       topic='bands'
      case('-if')
       topic='bands'
      case('-is')
       topic='bands'
      case('-o')
       topic='output'
      case('-kubo_csv')
       topic='curvature-csv'
      case('-kubo_pairs')
       topic='pairs-csv'
      case('-kubo_bundle')
       topic='bundle'
      case('-wf')
       topic='wavefunction-band'
      case('-k')
       topic='kpoint'
      case('-ng')
       topic='real-grid'
      case('-im')
       topic='imaginary'
      case('-theta')
       topic='theta'
      case('-phi')
       topic='phi'
      case('-ien')
       topic='ien'
      case('-fen')
       topic='fen'
      case('-nediv')
       topic='nediv'
      case('-sigma')
       topic='sigma'
      end select
      known=.false.
      select case(trim(topic))
      case('overview')
       known=.true.
      case('task')
       known=.true.
      case('options')
       known=.true.
      case('chern')
       known=.true.
      case('spin-chern')
       known=.true.
      case('spin-kubo')
       known=.true.
      case('z2')
       known=.true.
      case('kubo')
       known=.true.
      case('kubo-integral')
       known=.true.
      case('kubo-pairs')
       known=.true.
      case('optical')
       known=.true.
      case('spectrum')
       known=.true.
      case('wavefunction')
       known=.true.
      case('velocity')
       known=.true.
      case('bands')
       known=.true.
      case('mesh')
       known=.true.
      case('wavecar','input-dir')
       known=.true.
      case('spinor')
       known=.true.
      case('output')
       known=.true.
      case('curvature-csv')
       known=.true.
      case('pairs-csv')
       known=.true.
      case('bundle','per-band')
       known=.true.
      case('spin-axis')
       known=.true.
      case('outcar','incar','waveder','kubo-source')
       known=.true.
      case('energy-gap-tol')
       known=.true.
      case('spin-gap-tol')
       known=.true.
      case('sum-bands')
       known=.true.
      case('wavefunction-band')
       known=.true.
      case('kpoint')
       known=.true.
      case('real-grid')
       known=.true.
      case('imaginary')
       known=.true.
      case('theta')
       known=.true.
      case('phi')
       known=.true.
      case('ien')
       known=.true.
      case('fen')
       known=.true.
      case('nediv')
       known=.true.
      case('sigma')
       known=.true.
      case('build')
       known=.true.
      case('outputs')
       known=.true.
      case('examples')
       known=.true.
      case('all')
       known=.true.
      end select
      if(.not.known)then
       call vaspberry_help_error('unknown help topic: '//trim(requested))
      endif
      if(legacy_build)then
       known=.false.
       select case(trim(topic))
       case('overview')
        known=.true.
       case('task')
        known=.true.
       case('options')
        known=.true.
       case('chern')
        known=.true.
       case('z2')
        known=.true.
       case('optical')
        known=.true.
       case('wavefunction')
        known=.true.
       case('velocity')
        known=.true.
       case('bands')
        known=.true.
       case('mesh')
        known=.true.
       case('wavecar','input-dir')
        known=.true.
       case('spinor')
        known=.true.
       case('output')
        known=.true.
       case('wavefunction-band')
        known=.true.
       case('kpoint')
        known=.true.
       case('real-grid')
        known=.true.
       case('imaginary')
        known=.true.
       case('build')
        known=.true.
       case('outputs')
        known=.true.
       case('examples')
        known=.true.
       case('all')
        known=.true.
       end select
       if(.not.known)then
        call vaspberry_help_error(
     &   'help topic unavailable in this legacy serial build: '//trim(topic))
       endif
      endif
      if(topic.eq.'all')call help_all(ver_tag)
      call vaspberry_help_rank(print_help)
      if(print_help)then
       write(6,'(A)')'          **** PROGRAM INSTRUCTION ***'
       write(6,'(A)')trim(ver_tag)
       write(6,'(A)')'Topic: '//trim(topic)
       write(6,'(A)')'In examples, replace vaspberry with your executable path.'
       write(6,'(A)')''
       if(legacy_build)then
        call vaspberry_help_legacy_topic(topic)
       else
        call vaspberry_help_topic(topic)
       endif
      endif
      call vaspberry_help_finish(0)
      end subroutine vaspberry_help_show

      subroutine vaspberry_help_topic(topic)
      implicit none
      character*(*) topic
      select case(trim(topic))
      case('overview')
       write(6,'(A)')'Usage: vaspberry --task TASK [OPTIONS]'
       write(6,'(A)')'       vaspberry --help [TOPIC]   (or -h [TOPIC])'
       write(6,'(A)')''
       write(6,'(A)')'Choose a calculation'
       write(6,'(A)')'  Topology:   chern, spin-chern, z2'
       write(6,'(A)')'  Curvature:  kubo, kubo-integral, kubo-pairs, spin-kubo'
       write(6,'(A)')'  Other:      optical, spectrum, wavefunction, velocity'
       write(6,'(A)')''
       write(6,'(A)')'Start here'
       write(6,'(A)')'  --help task       Compare calculations and input requirements'
       write(6,'(A)')'  --help kubo       Read one calculation''s example, input and output'
       write(6,'(A)')'  --help bands      Check an option''s values and meaning'
       write(6,'(A)')'  --help options    Browse options by purpose'
       write(6,'(A)')'  --help outputs    Understand numerical data and plotting routes'
       write(6,'(A)')'  --help examples   Copy complete native commands'
       write(6,'(A)')'  --help build      Build Intel/GNU serial or MPI executables'
       write(6,'(A)')'  --help all        Full reference, including legacy syntax'
       write(6,'(A)')'  --help legacy     Alias of --help all'
       write(6,'(A)')''
       write(6,'(A)')'Inputs default to cwd; --input-dir DIR changes omitted source filenames.'
       write(6,'(A)')'Standard Kubo requires same-run WAVECAR, WAVEDER, INCAR and OUTCAR.'
       write(6,'(A)')'Use a fresh result directory. Help never reads inputs or writes results.'
       write(6,'(A)')'More: docs/NATIVE_COMMANDS.md and docs/HANDS_ON.md'
      case('task')
       write(6,'(A)')'Select one calculation with --task NAME; default: chern.'
       write(6,'(A)')''
       write(6,'(A)')'Topology'
       write(6,'(A)')'  chern          Chern number and Berry flux; full periodic 2D mesh'
       write(6,'(A)')'  spin-chern     Projected-spin sector Chern numbers; full mesh'
       write(6,'(A)')'  z2             FH n-field Z2 invariant; full even SOC mesh'
       write(6,'(A)')''
       write(6,'(A)')'Curvature and integration'
       write(6,'(A)')'  kubo           WAVEDER selected/occupied curvature; default optical source'
       write(6,'(A)')'  kubo-integral  Full-mesh Kubo-formula integration diagnostic'
       write(6,'(A)')'  kubo-pairs     Canonical pair export; explicit --kubo-source wavecar'
       write(6,'(A)')'  spin-kubo      Approximate projected-spin sector local curvature'
       write(6,'(A)')''
       write(6,'(A)')'Other calculations'
       write(6,'(A)')'  optical        Selected-transition circular selectivity'
       write(6,'(A)')'  spectrum       Broadened left/right transition spectra'
       write(6,'(A)')'  wavefunction   Real-space wavefunction reconstruction at Gamma'
       write(6,'(A)')'  velocity       Canonical-momentum velocity expectation'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  --help NAME gives the command, inputs, outputs and validity conditions.'
       write(6,'(A)')'  kubo-line is a synonym of kubo.'
       write(6,'(A)')'  chern/spin-chern use the Fukui method (FHS); Z2 uses the distinct FH method.'
       write(6,'(A)')'  Explicit --task rejects conflicting legacy calculation flags.'
      case('options')
       write(6,'(A)')'Each calculation option takes one separate value; no --name=value.'
       write(6,'(A)')''
       write(6,'(A)')'Inputs and selection'
       write(6,'(A)')'  task  input-dir  wavecar  waveder  incar  outcar  kubo-source'
       write(6,'(A)')'  bands  mesh  spinor'
       write(6,'(A)')''
       write(6,'(A)')'Outputs'
       write(6,'(A)')'  output  curvature-csv  pairs-csv  per-band'
       write(6,'(A)')''
       write(6,'(A)')'Spin analysis'
       write(6,'(A)')'  spin-axis  energy-gap-tol  spin-gap-tol  sum-bands'
       write(6,'(A)')''
       write(6,'(A)')'Wavefunction and optics'
       write(6,'(A)')'  wavefunction-band  kpoint  real-grid  imaginary  theta  phi'
       write(6,'(A)')'  Optical controls retain short flags: -ien, -fen, -nediv, -sigma.'
       write(6,'(A)')''
       write(6,'(A)')'Usage'
       write(6,'(A)')'  --help bands or --help --bands opens an option topic.'
       write(6,'(A)')'  List values have no spaces: --mesh 12,12 --bands 1:18'
       write(6,'(A)')'  Match example band IDs and meshes to your own VASP data.'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  build  outputs  examples  all'
      case('chern')
       write(6,'(A)')'Chern number and plaquette Berry flux via the Fukui method (FHS).'
       write(6,'(A)')''
       write(6,'(A)')'Run'
       write(6,'(A)')'  vaspberry --task chern --wavecar WAVECAR '//achar(92)
       write(6,'(A)')'    --mesh 12,12 --bands 1:18'
       write(6,'(A)')'  Replace example bands and mesh with those of your VASP calculation.'
       write(6,'(A)')''
       write(6,'(A)')'Inputs'
       write(6,'(A)')'  WAVECAR containing a full periodic uniform 2D mesh.'
       write(6,'(A)')'  --bands N selects one band; FIRST:LAST selects a band subspace.'
       write(6,'(A)')'  Omitted bands infer an occupied boundary from source occupations.'
       write(6,'(A)')''
       write(6,'(A)')'Outputs'
       write(6,'(A)')'  BERRYCURV.dat: plaquette values and a Chern-number header.'
       write(6,'(A)')'  --output LABEL writes BERRYCURV.LABEL.dat, in cwd.'
       write(6,'(A)')'  DAT coordinates are periodically repeated for plotting by default.'
       write(6,'(A)')'  Existing DAT files may be replaced; use a fresh result directory.'
       write(6,'(A)')''
       write(6,'(A)')'Conditions'
       write(6,'(A)')'  Use an isolated subspace and check mesh/link convergence.'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  --help bands | --help mesh | --help outputs'
       write(6,'(A)')'  docs/NATIVE_COMMANDS.md'
      case('spin-chern')
       write(6,'(A)')'Projected-spin sector Chern numbers via the Fukui method (FHS).'
       write(6,'(A)')''
       write(6,'(A)')'Run'
       write(6,'(A)')'  vaspberry --task spin-chern --input-dir /data/soc '//achar(92)
       write(6,'(A)')'    --mesh 12,12 --bands 1:8 --spin-axis z'
       write(6,'(A)')'  Replace the example mesh and bands with those of your input.'
       write(6,'(A)')''
       write(6,'(A)')'Inputs'
       write(6,'(A)')'  ISPIN=1 WAVECAR with two-component spinors and matching OUTCAR.'
       write(6,'(A)')'  A full periodic 2D mesh and an explicit isolated band subspace.'
       write(6,'(A)')'  OUTCAR supplies the source spin frame; inputs default to cwd.'
       write(6,'(A)')''
       write(6,'(A)')'Outputs'
       write(6,'(A)')'  SPIN_CHERN.csv: C_plus, C_minus, C_charge, C_spin.'
       write(6,'(A)')'  SPIN_BERRY.csv: sector plaquette data.'
       write(6,'(A)')'  SPIN_SPECTRUM.csv: projected-spin eigenvalues and gap diagnostics.'
       write(6,'(A)')'  C_spin=(C_plus-C_minus)/2. --output NAME replaces the SPIN prefix.'
       write(6,'(A)')'  Output files must be new.'
       write(6,'(A)')''
       write(6,'(A)')'Common options'
       write(6,'(A)')'  --spin-axis x|y|z|A,B,C   Cartesian unit axis; default z'
       write(6,'(A)')'  --energy-gap-tol VALUE   Selected energy gap, eV; default 1e-8'
       write(6,'(A)')'  --spin-gap-tol VALUE     Dimensionless spin gap; default 1e-6'
       write(6,'(A)')''
       write(6,'(A)')'Conditions'
       write(6,'(A)')'  The selected energy group must be isolated from excluded states.'
       write(6,'(A)')'  The projected-Pauli spectrum must stay separated from zero.'
       write(6,'(A)')'  Uses the pseudo-wavefunction metric, not a PAW-complete metric.'
       write(6,'(A)')'  Check mesh, energy-gap and spin-gap convergence before interpretation.'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  --help spin-axis | --help outcar | --help bands'
       write(6,'(A)')'  docs/SPIN_CHERN.md'
      case('spin-kubo')
       write(6,'(A)')'Approximate projected-spin sector Berry curvature from a Kubo formula.'
       write(6,'(A)')''
       write(6,'(A)')'Run'
       write(6,'(A)')'  vaspberry --task spin-kubo --kubo-source wavecar '//achar(92)
       write(6,'(A)')'    --input-dir /data/soc --bands 1:10 --spin-axis z'
       write(6,'(A)')''
       write(6,'(A)')'Inputs'
       write(6,'(A)')'  ISPIN=1 WAVECAR with two-component spinors and matching OUTCAR.'
       write(6,'(A)')'  An explicit isolated energy group with an open projected-spin gap.'
       write(6,'(A)')'  WAVEDER lacks the spin matrices this task requires; use explicit wavecar.'
       write(6,'(A)')''
       write(6,'(A)')'Outputs'
       write(6,'(A)')'  SPIN_KUBO.csv: local sector Berry curvature at supplied k points.'
       write(6,'(A)')'  SPIN_KUBO_SPECTRUM.csv: projected-spin spectrum and diagnostics.'
       write(6,'(A)')'  Explicit --mesh NX,NY also writes SPIN_KUBO_INTEGRAL.csv.'
       write(6,'(A)')'  --output NAME changes the SPIN prefix. Output files must be new.'
       write(6,'(A)')''
       write(6,'(A)')'Common options'
       write(6,'(A)')'  --spin-axis x|y|z|A,B,C   Analysis axis; default z'
       write(6,'(A)')'  --sum-bands N            Intermediate states 1:N; default all'
       write(6,'(A)')'  --mesh NX,NY             Validate a full mesh and estimate integrals'
       write(6,'(A)')''
       write(6,'(A)')'Conditions'
       write(6,'(A)')'  A path gives local curvature only; --mesh creates no missing k points.'
       write(6,'(A)')'  Uses a canonical-momentum response tangent with source-energy'
       write(6,'(A)')'  denominators and projected-spin derivatives in the pseudo metric.'
       write(6,'(A)')'  No integer invariant or conventional spin-current Hall is reported.'
       write(6,'(A)')'  Check energy/spin gaps, mesh and intermediate-band convergence.'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  --help sum-bands | --help outcar | --help energy-gap-tol'
       write(6,'(A)')'  --help spin-gap-tol | docs/SPIN_KUBO.md'
      case('z2')
       write(6,'(A)')'Fukui-Hatsugai (FH) n-field Z2 invariant; distinct from FHS Chern numbers.'
       write(6,'(A)')''
       write(6,'(A)')'Run'
       write(6,'(A)')'  vaspberry --task z2 --wavecar WAVECAR '//achar(92)
       write(6,'(A)')'    --mesh 12,12 --bands 1:10'
       write(6,'(A)')''
       write(6,'(A)')'Inputs'
       write(6,'(A)')'  Gapped, nonmagnetic, time-reversal-symmetric SOC WAVECAR.'
       write(6,'(A)')'  ISPIN=1 with two-component spinors; select all occupied bands 1:NE.'
       write(6,'(A)')'  NE is the occupied-band count: even, with NBANDS>NE.'
       write(6,'(A)')'  Use ICHARG=11, ISYM=-1 and a full even Gamma-centered 2D mesh.'
       write(6,'(A)')'  Nx,Ny>=4; kz=0 modulo G. Include each of the four 2D'
       write(6,'(A)')'  time-reversal-invariant momenta (TRIM) exactly once.'
       write(6,'(A)')''
       write(6,'(A)')'Outputs'
       write(6,'(A)')'  PASS writes Z2_FIELD.csv and NFIELD.dat.'
       write(6,'(A)')'  --output LABEL changes only the DAT label; CSV remains Z2_FIELD.csv.'
       write(6,'(A)')'  Rejected runs may leave Z2_FIELD.invalid.csv.'
       write(6,'(A)')''
       write(6,'(A)')'Conditions'
       write(6,'(A)')'  Report z2_invariant only for result_status=PASS and agreeing half-BZ'
       write(6,'(A)')'  parities. Check denser even meshes before physical interpretation.'
       write(6,'(A)')'  Legacy selector: -z2 1. Half-mesh -hf 1 is rejected.'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  --help mesh | --help bands | docs/Z2_FUKUI_HATSUGAI.md'
      case('kubo')
       write(6,'(A)')'Berry curvature from the Kubo formula; default source: WAVEDER.'
       write(6,'(A)')''
       write(6,'(A)')'Run'
       write(6,'(A)')'  vaspberry --task kubo --input-dir /data/optics'
       write(6,'(A)')'  Omit --bands for verified fixed occupied T0 filling.'
       write(6,'(A)')''
       write(6,'(A)')'Inputs'
       write(6,'(A)')'  Same-run WAVECAR, WAVEDER, INCAR and OUTCAR; default cwd.'
       write(6,'(A)')'  --input-dir DIR changes omitted filenames; explicit relative paths stay'
       write(6,'(A)')'  relative to invocation cwd. --wavecar overrides only WAVECAR.'
       write(6,'(A)')'  Supported: VASP 5.4.4 longitudinal static non-PEAD LOPTICS.'
       write(6,'(A)')''
       write(6,'(A)')'Outputs'
       write(6,'(A)')'  KUBO_WAVEDER.csv: k coordinates, omega_z_A2 (Angstrom^2), exact'
       write(6,'(A)')'  band IDs, source metadata and diagnostics. CSV paths must be new.'
       write(6,'(A)')'  Point/path data alone give no BZ integral. Explicit --mesh validates'
       write(6,'(A)')'  a full grid and adds unrounded integration metadata.'
       write(6,'(A)')''
       write(6,'(A)')'Common options'
       write(6,'(A)')'  --bands N or FIRST:LAST   Selected geometric band/subspace contribution'
       write(6,'(A)')'  --bands 31,33:34          Noncontiguous WAVEDER selection'
       write(6,'(A)')'  --per-band 1             Resolved individual bands; requires --bands'
       write(6,'(A)')'  --mesh NX,NY             Full stored 2D mesh for integration'
       write(6,'(A)')'  --curvature-csv PATH     Exact output filename'
       write(6,'(A)')''
       write(6,'(A)')'Conditions'
       write(6,'(A)')'  Selected boundaries must not split full-source transitive 0.002 eV'
       write(6,'(A)')'  degenerate clusters. Internal trace degeneracies are allowed.'
       write(6,'(A)')'  Required optical pairs must be stored; no invented missing elements'
       write(6,'(A)')'  or second energy denominator. Check k/NBANDS convergence.'
       write(6,'(A)')'  Explicit-band integrals are geometric contributions; only verified'
       write(6,'(A)')'  omitted-band occupied integration labels total physical Hall.'
       write(6,'(A)')'  Native integration sums both ISPIN=2 channels; Python kubo-hall'
       write(6,'(A)')'  selects --spin (default 1). Sum separate spin 1/2 runs for charge.'
       write(6,'(A)')''
       write(6,'(A)')'Approximation and compatibility'
       write(6,'(A)')'  --kubo-source wavecar explicitly selects canonical momentum, omitting'
       write(6,'(A)')'  full PAW/nonlocal velocity terms. No automatic WAVECAR fallback.'
       write(6,'(A)')'  Canonical trace default: KUBO.csv; isolated/per-band gaps >1e-5 eV.'
       write(6,'(A)')'  Standard Kubo rejects wavefunction/optical/folding controls:'
       write(6,'(A)')'  --kpoint/--real-grid/--imaginary, --theta/--phi, legacy aliases,'
       write(6,'(A)')'  -ishift, -ien/-fen/-nediv/-sigma, -kp/-skp. No silent broadening.'
       write(6,'(A)')'  kubo-line is a synonym; legacy -kubo uses the same source barrier.'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  --help bands | --help waveder | --help kubo-source'
       write(6,'(A)')'  --help curvature-csv | docs/WAVEDER_KUBO_PROTOCOL.md'
      case('kubo-integral')
       write(6,'(A)')'Full-mesh Kubo integral; default source: standard WAVEDER.'
       write(6,'(A)')''
       write(6,'(A)')'Run'
       write(6,'(A)')'  vaspberry --task kubo-integral --input-dir /data/optics '//achar(92)
       write(6,'(A)')'    --mesh 12,12'
       write(6,'(A)')''
       write(6,'(A)')'Inputs'
       write(6,'(A)')'  Same-run WAVECAR, WAVEDER, INCAR and OUTCAR.'
       write(6,'(A)')'  An explicit --mesh matching the full stored periodic 2D grid.'
       write(6,'(A)')'  Source, path and gap requirements: --help kubo.'
       write(6,'(A)')''
       write(6,'(A)')'Outputs'
       write(6,'(A)')'  KUBO_WAVEDER.csv: omega_z_A2 in Angstrom^2 and integration metadata.'
       write(6,'(A)')'  --curvature-csv PATH changes the filename; CSV paths must be new.'
       write(6,'(A)')'  The finite-grid geometric integral is unrounded and dimensionless.'
       write(6,'(A)')''
       write(6,'(A)')'Conditions'
       write(6,'(A)')'  Explicit --bands gives a unit-weight selected Chern contribution,'
       write(6,'(A)')'  not total Hall. A path cannot be integrated as a full mesh.'
       write(6,'(A)')'  Canonical diagnostic requires explicit --kubo-source wavecar.'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  Occupation-weighted Hall: python3 tools/vaspberry_kubo.py kubo-hall'
       write(6,'(A)')'  --help mesh | --help curvature-csv | docs/WAVEDER_KUBO_PROTOCOL.md'
      case('kubo-pairs')
       write(6,'(A)')'Canonical-momentum pair export; explicit --kubo-source wavecar.'
       write(6,'(A)')''
       write(6,'(A)')'Run'
       write(6,'(A)')'  vaspberry --task kubo-pairs --kubo-source wavecar '//achar(92)
       write(6,'(A)')'    --wavecar WAVECAR --pairs-csv PAIRS.csv'
       write(6,'(A)')''
       write(6,'(A)')'Inputs'
       write(6,'(A)')'  WAVECAR; omit --bands to export all unordered source-band pairs n<m.'
       write(6,'(A)')'  Standard WAVEDER cannot promise all unordered source-band pairs.'
       write(6,'(A)')'  Use --task kubo for standard selected/occupied curvature.'
       write(6,'(A)')''
       write(6,'(A)')'Outputs'
       write(6,'(A)')'  PAIRS.csv: energies, fractional k coordinates and Cartesian pair'
       write(6,'(A)')'  numerators (yz,zx,xy), in eV^2 Angstrom^2.'
       write(6,'(A)')'  No gap division or occupation weighting; pairs are not conductivity.'
       write(6,'(A)')'  PATH and PATH.partial must be new; require result_status=PASS.'
       write(6,'(A)')''
       write(6,'(A)')'Conditions'
       write(6,'(A)')'  Canonical momentum omits full PAW/nonlocal velocity terms.'
       write(6,'(A)')'  Source occupations may vary with k; no physical filling is assumed.'
       write(6,'(A)')'  Incompatible with --bands, --curvature-csv and --per-band.'
       write(6,'(A)')'  Full-BZ sampling is needed for later transport integration.'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  python3 tools/vaspberry_kubo.py import-pairs --help'
       write(6,'(A)')'  python3 tools/vaspberry_kubo.py pair-hall --help'
       write(6,'(A)')'  docs/OUTPUT_FORMAT.md | docs/WAVEDER_KUBO_PROTOCOL.md'
      case('optical')
       write(6,'(A)')'Circular optical selectivity for one initial-to-final transition.'
       write(6,'(A)')''
       write(6,'(A)')'Run'
       write(6,'(A)')'  vaspberry --task optical --wavecar WAVECAR '//achar(92)
       write(6,'(A)')'    --mesh 12,12 --bands 18:19 --theta 0 --phi 0'
       write(6,'(A)')''
       write(6,'(A)')'Inputs'
       write(6,'(A)')'  WAVECAR on the stored 2D mesh; --bands endpoints are initial/final'
       write(6,'(A)')'  bands, not a subspace trace. Defaults infer VBM/CBM.'
       write(6,'(A)')'  Angles in degrees; (theta,phi)=(0,0) means light along z.'
       write(6,'(A)')''
       write(6,'(A)')'Outputs'
       write(6,'(A)')'  BERRYCURV.dat: normalized circular selectivity; legacy filename.'
       write(6,'(A)')'  --output LABEL writes LABEL.dat. Existing DAT files may be replaced.'
       write(6,'(A)')''
       write(6,'(A)')'Conditions'
       write(6,'(A)')'  Uses pseudo-wavefunction momentum matrix elements.'
       write(6,'(A)')'  Choose endpoints explicitly for reproducibility. Legacy flag: -cd 1.'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  --help spectrum | --help theta | --help phi'
      case('spectrum')
       write(6,'(A)')'Broadened left/right optical transition spectra (legacy -cd 2).'
       write(6,'(A)')''
       write(6,'(A)')'Run'
       write(6,'(A)')'  vaspberry --task spectrum --wavecar WAVECAR '//achar(92)
       write(6,'(A)')'    --mesh 12,12 --bands 1:20 '//achar(92)
       write(6,'(A)')'    -ien 1 -fen 3 -nediv 201 -sigma 0.05'
       write(6,'(A)')''
       write(6,'(A)')'Inputs'
       write(6,'(A)')'  WAVECAR; occupied bands 1:VBM and empty bands VBM+1:LAST.'
       write(6,'(A)')'  LAST is the upper --bands endpoint; FIRST does not cap occupied bands.'
       write(6,'(A)')''
       write(6,'(A)')'Outputs'
       write(6,'(A)')'  OPT_TRANS_RATE_LEFT.dat and OPT_TRANS_RATE_RIGHT.dat plus per-k _KP#'
       write(6,'(A)')'  files: LEFT/RIGHT intensity diagnostics in arbitrary units.'
       write(6,'(A)')'  --output LABEL uses CIRC_DICHROISM_W.LABEL_LEFT/RIGHT.dat.'
       write(6,'(A)')'  Energies <=1e-11 eV are omitted. DAT files may be replaced.'
       write(6,'(A)')''
       write(6,'(A)')'Common options'
       write(6,'(A)')'  -ien 0 -fen 10 -nediv 1000 -sigma 0.01   Default spectrum controls'
       write(6,'(A)')'  Photon energies and Gaussian standard deviation are in eV.'
       write(6,'(A)')''
       write(6,'(A)')'Conditions'
       write(6,'(A)')'  Use -fen > -ien, -nediv >= 2 and finite positive -sigma.'
       write(6,'(A)')'  This is a transition diagnostic, not a calibrated absolute rate.'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  --help bands | --help theta | --help ien | --help sigma'
       write(6,'(A)')'  examples/features/circular-dichroism/'
      case('wavefunction')
       write(6,'(A)')'Reconstruct real-space wavefunction amplitudes at Gamma.'
       write(6,'(A)')''
       write(6,'(A)')'Run'
       write(6,'(A)')'  vaspberry --task wavefunction --wavecar WAVECAR '//achar(92)
       write(6,'(A)')'    --wavefunction-band 18 --kpoint 24 '//achar(92)
       write(6,'(A)')'    --real-grid 24,24,64 --imaginary 1'
       write(6,'(A)')'  Index 24 is Gamma only in the included MoS2 path example.'
       write(6,'(A)')''
       write(6,'(A)')'Inputs'
       write(6,'(A)')'  WAVECAR and matching POSCAR/EIGENVAL in the working directory.'
       write(6,'(A)')'  POSCAR/EIGENVAL are required; --input-dir does not relocate them.'
       write(6,'(A)')'  Choose your file''s Gamma point. Indices are one-based.'
       write(6,'(A)')''
       write(6,'(A)')'Outputs'
       write(6,'(A)')'  PARCHG-W-K###-E###-SPIN#; optional -IM-SPIN# companion.'
       write(6,'(A)')'  CHGCAR-style container with real/imaginary amplitudes, one block'
       write(6,'(A)')'  per spinor component; this is not a charge-density calculation.'
       write(6,'(A)')'  --output does not rename these files. Existing files may be replaced.'
       write(6,'(A)')''
       write(6,'(A)')'Common options'
       write(6,'(A)')'  --kpoint 1        Default first source k point; must be Gamma'
       write(6,'(A)')'  --real-grid      Default basis-derived grid'
       write(6,'(A)')'  --imaginary 0    Default no imaginary companion'
       write(6,'(A)')''
       write(6,'(A)')'Conditions'
       write(6,'(A)')'  The retained reconstruction supports Gamma; use a fresh directory.'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  --help wavefunction-band | --help real-grid | --help imaginary'
       write(6,'(A)')'  examples/features/wavefunction/'
      case('velocity')
       write(6,'(A)')'Canonical-momentum velocity expectation for one selected band.'
       write(6,'(A)')''
       write(6,'(A)')'Run'
       write(6,'(A)')'  vaspberry --task velocity --wavecar WAVECAR --bands 18'
       write(6,'(A)')''
       write(6,'(A)')'Inputs'
       write(6,'(A)')'  WAVECAR; choose one band. This routine uses the first selected band.'
       write(6,'(A)')''
       write(6,'(A)')'Outputs'
       write(6,'(A)')'  VEL_EXPT.dat: x/y velocity components in m/s.'
       write(6,'(A)')'  --output LABEL writes VEL_EXPT.LABEL.dat; scalar ISPIN=2 adds UP/DN.'
       write(6,'(A)')'  DAT files may be replaced; use a fresh working directory.'
       write(6,'(A)')''
       write(6,'(A)')'Conditions'
       write(6,'(A)')'  Bare momentum in the pseudo-wavefunction metric; excludes full'
       write(6,'(A)')'  PAW/nonlocal velocity terms. Legacy selector: -vel 1.'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  --help bands | --help output | --help outputs'
      case('bands')
       write(6,'(A)')'--bands N | FIRST:LAST: one-based inclusive selection, no duplicates.'
       write(6,'(A)')''
       write(6,'(A)')'Selection and defaults'
       write(6,'(A)')'  N:N is a singleton, like N; multiple selected IDs form a subspace.'
       write(6,'(A)')'  WAVEDER Kubo also accepts N,M:P, for example --bands 31,33:34'
       write(6,'(A)')'  Inferred defaults follow source occupations, not a fixed band 18.'
       write(6,'(A)')'  Spin tasks require explicit selection; pair export requires omission.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning by task'
       write(6,'(A)')'  chern          One band or the selected link-variable subspace'
       write(6,'(A)')'  WAVEDER kubo   Selected geometry; omitted bands verify occupied T0'
       write(6,'(A)')'  wavecar kubo   One band individually; multiple selected bands as a trace'
       write(6,'(A)')'  spin tasks     Explicit isolated energy subspace'
       write(6,'(A)')'  z2             All occupied bands 1:NE; even rank and NBANDS>NE'
       write(6,'(A)')'  optical        Initial/final endpoints'
       write(6,'(A)')'  spectrum       LAST caps empty bands; occupied bands always start at 1'
       write(6,'(A)')'  velocity       First selected band only'
       write(6,'(A)')'  kubo-pairs     Omit --bands; all source-band pairs are exported'
       write(6,'(A)')''
       write(6,'(A)')'Compatibility'
       write(6,'(A)')'  Prefer one --bands selection. Repeated/mixed selectors print a warning.'
       write(6,'(A)')'  A later --bands replaces earlier selection; -ii/-if/-is retain'
       write(6,'(A)')'  historical endpoint/singleton precedence after scalar/range selection.'
       write(6,'(A)')'  Legacy flags cannot modify an earlier comma-list; --bands can replace'
       write(6,'(A)')'  legacy selection. Inspect effective band IDs before interpretation.'
       write(6,'(A)')'  With explicit --kubo-source wavecar, --task kubo uses the range rule;'
       write(6,'(A)')'  pure legacy -kubo -ii/-if keeps individual bands.'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  --help per-band | --help kubo-source'
       write(6,'(A)')'  --sum-bands separately controls spin-kubo intermediate states.'
      case('mesh')
       write(6,'(A)')'--mesh NX,NY : dimensions of the stored uniform periodic 2D mesh.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  12,12; match the actual stored mesh.'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  Stored 2D meshes; explicit integration rules below.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Example: --mesh 12,12 ; -kx 12 -ky 12 sets the same mesh dimensions.'
       write(6,'(A)')'  For WAVEDER or spin Kubo, only explicit --mesh requests integration;'
       write(6,'(A)')'  legacy -kx/-ky alone set sizes and leave point curvature unintegrated.'
       write(6,'(A)')'  WAVEDER kubo-integral requires --mesh, not only -kx/-ky.'
       write(6,'(A)')'  Default native dimensions are 12,12; specify the actual VASP mesh.'
       write(6,'(A)')'  This option does not generate k points, interpolate, or make a path a mesh.'
       write(6,'(A)')'  chern/spin-chern require a matching full mesh; omitted dimensions use 12,12.'
       write(6,'(A)')'  z2 additionally needs an even Gamma-centered mesh, Nx,Ny>=4.'
       write(6,'(A)')'  kubo and kubo-pairs accept arbitrary source k points without --mesh.'
       write(6,'(A)')'  spin-kubo evaluates points by default; an explicit mesh additionally'
       write(6,'(A)')'  requests SPIN_KUBO_INTEGRAL.csv and checks complete mesh geometry.'
       write(6,'(A)')'  For transport use a full BZ mesh and verify sampling convergence.'
      case('kubo-source')
       write(6,'(A)')'--kubo-source waveder|wavecar : Kubo operator source; default waveder.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  waveder'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  kubo, kubo-integral, kubo-pairs, spin-kubo'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  waveder: same-run optical connections; selected geometry or verified'
       write(6,'(A)')'    occupied trace.'
       write(6,'(A)')'  wavecar: explicit canonical-momentum approximation; prints a warning.'
       write(6,'(A)')'  Missing WAVEDER or unsupported requests stop; fallback is never automatic.'
       write(6,'(A)')'  kubo-pairs and spin-kubo presently require explicit wavecar.'
       write(6,'(A)')'  This choice changes the operator, not just input speed or file preference.'
       write(6,'(A)')'  Reference: docs/WAVEDER_KUBO_PROTOCOL.md.'
      case('input-dir')
       write(6,'(A)')'--input-dir DIR : base for omitted WAVECAR/WAVEDER/INCAR/OUTCAR filenames.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  Invocation cwd; explicit paths remain relative to cwd.'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  Omitted WAVECAR/WAVEDER/INCAR/OUTCAR paths; outputs stay in cwd.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Default is the invocation working directory. DIR must exist and be a'
       write(6,'(A)')'    directory.'
       write(6,'(A)')'  Example: --task kubo --input-dir /data/optics --curvature-csv result.csv'
       write(6,'(A)')'  Each explicit --wavecar/--waveder/--incar/--outcar overrides only that file.'
       write(6,'(A)')'  Explicit relative file paths are relative to cwd, never relative to DIR.'
       write(6,'(A)')'  Option order does not change precedence. --wavecar never redirects other'
       write(6,'(A)')'    inputs.'
       write(6,'(A)')'  Does not change cwd, outputs, or unrelated input lists. Empty/overlong paths'
       write(6,'(A)')'    fail.'
       write(6,'(A)')'  Standard Kubo prints all selected input paths before opening the source'
       write(6,'(A)')'    files.'
       write(6,'(A)')'  Wavefunction POSCAR/EIGENVAL remain required in cwd, not under DIR.'
      case('waveder')
       write(6,'(A)')'--waveder PATH : optical connection input for --kubo-source waveder.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  WAVEDER under --input-dir (default cwd).'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  Standard WAVEDER kubo/kubo-integral.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Default is WAVEDER under --input-dir (default cwd), independent of'
       write(6,'(A)')'    --wavecar.'
       write(6,'(A)')'  Explicit relative PATH is relative to cwd; keep same-run INCAR and OUTCAR.'
       write(6,'(A)')'  Supported rectangular coverage is sufficient; do not pad to NBANDS squared.'
       write(6,'(A)')'  Standard VASP 5.4.4 longitudinal LOPTICS branch only; PAW terms are already'
       write(6,'(A)')'  included upstream. The stored connection already contains the gap'
       write(6,'(A)')'    denominator.'
       write(6,'(A)')'  Each required pair needs a stored orientation; both must agree if present.'
       write(6,'(A)')'  No invented missing elements or second energy denominator.'
      case('incar')
       write(6,'(A)')'--incar PATH : same-run INCAR for the standard WAVEDER Kubo route.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  INCAR under --input-dir (default cwd).'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  Standard WAVEDER kubo/kubo-integral.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Default is INCAR under --input-dir (default cwd), independent of --wavecar.'
       write(6,'(A)')'  Explicit relative PATH is relative to cwd. It validates the producer branch.'
       write(6,'(A)')'  Use LOPTICS=T, LPEAD=F, LNABLA=F, NSW=0, ISYM=-1, LREAL=F.'
       write(6,'(A)')'  Keep the actual input; do not synthesize settings to make validation pass.'
      case('wavecar')
       write(6,'(A)')'--wavecar PATH : overrides only WAVECAR; default is --input-dir/WAVECAR.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  WAVECAR under --input-dir (default cwd).'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  Tasks reading VASP wavefunctions.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  --input-dir defaults to cwd. Explicit relative PATH is relative to cwd.'
       write(6,'(A)')'  Changing --wavecar never changes WAVEDER/INCAR/OUTCAR defaults.'
       write(6,'(A)')'  Example: --wavecar /data/project/WAVECAR ; legacy equivalent: -f PATH.'
       write(6,'(A)')'  Quote paths with spaces. Prefer short paths (native fields: 256 chars).'
       write(6,'(A)')'  Reads full-complex single-precision coefficients with RTAG=45200.'
       write(6,'(A)')'  GNU and supplied Intel targets use byte RECL. The native reader detects'
       write(6,'(A)')'  byte or four-byte-word producer layouts and reports physical byte stride.'
       write(6,'(A)')'  Ambiguous, truncated or invalid layouts fail before allocation; no'
       write(6,'(A)')'    conversion.'
       write(6,'(A)')'  --spinor auto checks scalar/two-component layout throughout the file.'
       write(6,'(A)')'  Help never attempts to open or inspect this input.'
      case('spinor')
       write(6,'(A)')'--spinor auto|1|2 : infer or assert the WAVECAR component layout.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  auto'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  WAVECAR component layout.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Default auto checks full-complex coefficient counts for every spin/k point.'
       write(6,'(A)')'  1 or 2 asserts that the file has that detected number of components.'
       write(6,'(A)')'  This is not a spin degeneracy multiplier and does not enable SOC.'
       write(6,'(A)')'  ISPIN counts channels; it does not by itself identify spinor components.'
       write(6,'(A)')'  Spin-sector and Z2 tasks require two-component spinors.'
       write(6,'(A)')'  Use matching OUTCAR for spin-frame conventions. Legacy equivalent: -s.'
      case('output')
       write(6,'(A)')'--output LABEL : legacy output label or spin-task prefix, not a directory.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  Task-specific filename; spin prefix SPIN. This is not a directory.'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  Legacy DAT labels and spin CSV prefixes.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Examples for LABEL=sample:'
       write(6,'(A)')'    chern       BERRYCURV.sample.dat'
       write(6,'(A)')'    optical     sample.dat'
       write(6,'(A)')'    spectrum    CIRC_DICHROISM_W.sample_LEFT/RIGHT.dat'
       write(6,'(A)')'    velocity    VEL_EXPT.sample.dat'
       write(6,'(A)')'    z2          sample.dat; CSV remains Z2_FIELD.csv'
       write(6,'(A)')'    spin-chern  sample_CHERN/BERRY/SPECTRUM.csv'
       write(6,'(A)')'    spin-kubo   sample_KUBO.csv and sample_KUBO_SPECTRUM.csv'
       write(6,'(A)')'    spin-kubo   sample_KUBO_INTEGRAL.csv only with explicit full --mesh'
       write(6,'(A)')'  Spin-task default prefix is SPIN. Scalar ISPIN=2 DAT files add UP/DN.'
       write(6,'(A)')'  --curvature-csv and --pairs-csv specify their own exact CSV paths.'
       write(6,'(A)')'  Wavefunction names are set by k/band; --output does not rename them.'
       write(6,'(A)')'  Use a fresh working directory; DAT outputs can replace existing files.'
       write(6,'(A)')'  Legacy equivalent: -o LABEL. See --help outputs for data usage.'
      case('curvature-csv')
       write(6,'(A)')'--curvature-csv PATH: exact Kubo curvature CSV filename.'
       write(6,'(A)')''
       write(6,'(A)')'Default and scope'
       write(6,'(A)')'  WAVEDER: KUBO_WAVEDER.csv. Canonical wavecar trace: KUBO.csv.'
       write(6,'(A)')'  Applies to kubo/kubo-integral (legacy -kubo 2/1), not pair export.'
       write(6,'(A)')''
       write(6,'(A)')'WAVEDER schemas'
       write(6,'(A)')'  Omitted bands: VASPBERRY_WAVEDER_KUBO_OCCUPIED_V1.'
       write(6,'(A)')'  Selected trace: VASPBERRY_WAVEDER_KUBO_BUNDLE_V1.'
       write(6,'(A)')'  Individual bands: VASPBERRY_WAVEDER_KUBO_BAND_V1.'
       write(6,'(A)')'  Records omega_z_A2, exact band IDs, optical metadata and gap diagnostics.'
       write(6,'(A)')''
       write(6,'(A)')'Canonical WAVECAR schemas'
       write(6,'(A)')'  Trace V2: spin,k,min_external_gap and selected range/source NBANDS.'
       write(6,'(A)')'  Band V3: spin,k,band,energy,min_gap,omega_z_A2.'
       write(6,'(A)')'  Individual-band CSV is opt-in; a selected trace always writes CSV.'
       write(6,'(A)')''
       write(6,'(A)')'File policy and units'
       write(6,'(A)')'  Parent directory must exist; PATH and PATH.partial must be new.'
       write(6,'(A)')'  Declared counts and terminal result_status=PASS mark complete output.'
       write(6,'(A)')'  omega_z_A2 is normalized curvature in Angstrom^2 for both sources.'
       write(6,'(A)')'  Preserve comment metadata with the numeric table.'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  --help kubo | --help per-band | docs/OUTPUT_FORMAT.md'
       write(6,'(A)')'  Legacy equivalent: -kubo_csv PATH'
      case('pairs-csv')
       write(6,'(A)')'--pairs-csv PATH : exact reusable native pair-export CSV filename.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  No default path; required for kubo-pairs.'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  kubo-pairs only.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Use with --task kubo-pairs --kubo-source wavecar; legacy -kubo_pairs PATH.'
       write(6,'(A)')'  The parent directory must exist; PATH and PATH.partial must be new.'
       write(6,'(A)')'  PATH.partial is renamed only after all records are written and closed.'
       write(6,'(A)')'  Omit --bands, --curvature-csv and --per-band.'
       write(6,'(A)')'  All unordered source-band pairs are exported, with no gap division or'
       write(6,'(A)')'  occupation weighting. Require the final result_status=PASS marker.'
       write(6,'(A)')'  See --help kubo-pairs for quantities, units and postprocessing.'
      case('per-band')
       write(6,'(A)')'--per-band 0|1: individual Kubo output when 1; default 0.'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  kubo, kubo-line, kubo-integral; incompatible with pair/spin/FHS tasks.'
       write(6,'(A)')''
       write(6,'(A)')'WAVEDER (default source)'
       write(6,'(A)')'  Requires explicit --bands and resolved full-source 0.002 eV clusters.'
       write(6,'(A)')'  Every required partner among all NBANDS must be stored.'
       write(6,'(A)')'  Always writes CSV; default KUBO_WAVEDER.csv. No occupation weights.'
       write(6,'(A)')''
       write(6,'(A)')'Explicit canonical WAVECAR'
       write(6,'(A)')'  Individual bands require all gaps >1e-5 eV at every source k/spin,'
       write(6,'(A)')'  including gaps to other selected bands; degeneracy fails before output.'
       write(6,'(A)')'  Selected traces require only external gaps >1e-5 eV and exclude'
       write(6,'(A)')'  internal transitions before division. No range is enlarged automatically.'
       write(6,'(A)')'  Trace: one CSV row per source k/spin, default KUBO.csv.'
       write(6,'(A)')'  Individual bands: CSV is opt-in; existing DAT output remains.'
       write(6,'(A)')'  Omitted bands retain occupied inference and individual output.'
       write(6,'(A)')'  Full-source trace is zero in the truncated basis, not a convergence test.'
       write(6,'(A)')'  These canonical modes apply no occupations.'
       write(6,'(A)')''
       write(6,'(A)')'Run'
       write(6,'(A)')'  vaspberry --task kubo --input-dir /data/optics '//achar(92)
       write(6,'(A)')'    --bands 1:18 --per-band 1 --curvature-csv BANDS.csv'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  --help bands | --help kubo-source | --help curvature-csv'
      case('bundle')
       write(6,'(A)')'--bundle and -kubo_bundle were removed; old commands fail explicitly.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  Removed option; migration help only.'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  No accepted calculation mode.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Use --task kubo --bands FIRST:LAST for the selected-space trace.'
       write(6,'(A)')'  Defaults: WAVEDER KUBO_WAVEDER.csv; explicit wavecar trace KUBO.csv.'
       write(6,'(A)')'  For either source use explicit --bands with --per-band 1 for resolved bands.'
       write(6,'(A)')'  WAVEDER also accepts disjoint selections such as --bands 31,33:34.'
      case('spin-axis')
       write(6,'(A)')'--spin-axis x|y|z|A,B,C : Cartesian analysis axis for spin-sector tasks.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  z'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  spin-chern, spin-kubo'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Default z; for example --spin-axis 0,0,1.'
       write(6,'(A)')'  For A,B,C supply a finite unit vector (norm squared within 1e-8 of 1).'
       write(6,'(A)')'  Used only by spin-chern and spin-kubo; it does not rotate VASP input.'
       write(6,'(A)')'  Matching OUTCAR gives the source SAXIS-to-Cartesian spin frame.'
       write(6,'(A)')'  Changing the analysis axis may close the projected-spin gap.'
      case('outcar')
       write(6,'(A)')'--outcar PATH: same-run source metadata; default --input-dir/OUTCAR.'
       write(6,'(A)')''
       write(6,'(A)')'Paths'
       write(6,'(A)')'  --input-dir defaults to cwd; explicit relative PATH stays relative to cwd.'
       write(6,'(A)')'  --wavecar does not redirect OUTCAR or any other auxiliary input.'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  Standard WAVEDER Kubo: validates the optical producer and run settings.'
       write(6,'(A)')'  spin-chern/spin-kubo: spin frame, SAXIS-to-Cartesian rotation and'
       write(6,'(A)')'  matching dimensions/lattice. Preserve the WAVECAR/OUTCAR pair.'
       write(6,'(A)')'  Charge Kubo with explicit wavecar does not use an OUTCAR override.'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  --help input-dir | --help incar | --help spin-axis'
       write(6,'(A)')'  OUTCAR does not replace or infer wavefunction coefficients.'
      case('energy-gap-tol')
       write(6,'(A)')'--energy-gap-tol VALUE : positive selected-to-excluded gap tolerance.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  1e-8 eV'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  spin-chern, spin-kubo'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Used by spin-chern and spin-kubo; units eV, default 1e-8.'
       write(6,'(A)')'  The selected energy subspace must stay isolated from excluded states.'
       write(6,'(A)')'  This numerical cutoff is not broadening or evidence of convergence.'
       write(6,'(A)')'  Charge thresholds: standard WAVEDER >0.002 eV; wavecar >1e-5 eV.'
      case('spin-gap-tol')
       write(6,'(A)')'--spin-gap-tol VALUE : positive projected-Pauli distance from zero.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  1e-6 (dimensionless)'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  spin-chern, spin-kubo'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Used by spin-chern and spin-kubo; dimensionless, default 1e-6.'
       write(6,'(A)')'  It is not an energy broadening or an eV gap.'
       write(6,'(A)')'  A closed spin gap makes the axis-dependent sector decomposition unresolved.'
       write(6,'(A)')'  Inspect the spin spectrum and convergence instead of relaxing it blindly.'
      case('sum-bands')
       write(6,'(A)')'--sum-bands N : spin-kubo-only intermediate-state convergence control.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  All stored source bands.'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  spin-kubo only.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Include source bands 1:N outside the selected group in the Kubo sum.'
       write(6,'(A)')'  Default is all stored source bands. N must include the selected group'
       write(6,'(A)')'  and cannot exceed source NBANDS. VASP eigenstates remain fixed.'
       write(6,'(A)')'  At least one intermediate state outside the selected group must remain.'
       write(6,'(A)')'  Example: --task spin-kubo --kubo-source wavecar --bands 1:10 --sum-bands 20'
       write(6,'(A)')'  This does not truncate source WAVECAR, change --bands, or apply to'
       write(6,'(A)')'  charge kubo, kubo-pairs, or spin-chern. Compare increasing cutoffs.'
      case('wavefunction-band')
       write(6,'(A)')'--wavefunction-band N : positive one-based band to reconstruct.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  Required; no default band.'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  wavefunction only.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Required with --task wavefunction. Legacy equivalent: -wf N.'
       write(6,'(A)')'  This is the reconstruction selector, separate from --bands.'
       write(6,'(A)')'  Related: --help wavefunction, --help kpoint, --help real-grid.'
      case('kpoint')
       write(6,'(A)')'--kpoint N : positive one-based source k index for wavefunction output.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  1 (choose the Gamma index).'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  wavefunction only.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Default 1. Choose the Gamma point for the existing reconstruction routine.'
       write(6,'(A)')'  This option does not select k points for curvature/transport tasks.'
       write(6,'(A)')'  Legacy equivalent: -k N. Related: --help wavefunction.'
      case('real-grid')
       write(6,'(A)')'--real-grid NX,NY,NZ : positive real-space wavefunction grid dimensions.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  Basis-derived real-space grid.'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  wavefunction only.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Default uses the basis-derived grid. Example: --real-grid 24,24,64'
       write(6,'(A)')'  Used for wavefunction reconstruction, not reciprocal-space sampling.'
       write(6,'(A)')'  Legacy equivalent: -ng NX,NY,NZ. Related: --help wavefunction.'
      case('imaginary')
       write(6,'(A)')'--imaginary 0|1 : write an imaginary wavefunction companion if 1.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  0'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  wavefunction only.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Default 0. Used with wavefunction reconstruction.'
       write(6,'(A)')'  Legacy equivalent: -im 1. The -IM-SPIN# file accompanies the real output.'
       write(6,'(A)')'  Related: --help wavefunction.'
      case('theta')
       write(6,'(A)')'--theta DEG : light incidence polar angle in degrees, default 0.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  0 degrees'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  optical, spectrum'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Used by optical/spectrum; legacy equivalent: -theta DEG.'
       write(6,'(A)')'  (theta,phi)=(0,0) is along z; (90,0) along x; (90,90) along y.'
       write(6,'(A)')'  This is a light-direction control, not the spin analysis axis.'
      case('phi')
       write(6,'(A)')'--phi DEG : light incidence azimuthal angle in degrees, default 0.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  0 degrees'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  optical, spectrum'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Used by optical/spectrum; legacy equivalent: -phi DEG.'
       write(6,'(A)')'  (theta,phi)=(0,0) is along z; (90,0) along x; (90,90) along y.'
       write(6,'(A)')'  Related: --help theta and --help optical.'
      case('ien')
       write(6,'(A)')'-ien EV : initial photon energy for the optical spectrum, default 0 eV.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  0 eV'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  spectrum'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Together with -fen EV (default 10 eV), sets the energy window.'
       write(6,'(A)')'  Example: --task spectrum --bands 1:20 -ien 1 -fen 3 -nediv 201'
       write(6,'(A)')'  These controls retain their legacy short-option names.'
       write(6,'(A)')'  Related: --help spectrum, --help nediv, --help sigma.'
       write(6,'(A)')'  Use a finite increasing photon-energy window with -nediv >= 2.'
      case('fen')
       write(6,'(A)')'-fen EV : final photon energy for the optical spectrum, default 10 eV.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  10 eV'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  spectrum'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Use together with -ien EV (default 0 eV) and -nediv N (default 1000).'
       write(6,'(A)')'  Units eV; this option retains its legacy short name.'
       write(6,'(A)')'  Related: --help spectrum and --help ien.'
       write(6,'(A)')'  The final energy must exceed the initial energy.'
      case('nediv')
       write(6,'(A)')'-nediv N : number of photon-energy points for the optical spectrum.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  1000 points'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  spectrum'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Default 1000; combine with -ien and -fen, for example -nediv 201.'
       write(6,'(A)')'  This option retains its legacy short name.'
       write(6,'(A)')'  Related: --help spectrum and --help sigma.'
       write(6,'(A)')'  Use N >= 2 for the retained spectrum grid.'
      case('sigma')
       write(6,'(A)')'-sigma EV : Gaussian broadening for optical transition spectra.'
       write(6,'(A)')''
       write(6,'(A)')'Default'
       write(6,'(A)')'  0.01 eV'
       write(6,'(A)')''
       write(6,'(A)')'Applies to'
       write(6,'(A)')'  spectrum'
       write(6,'(A)')''
       write(6,'(A)')'Meaning'
       write(6,'(A)')'  Default 0.01 eV; for example -sigma 0.05'
       write(6,'(A)')'  It does not change spin-gap/energy-gap tolerances or Kubo denominators.'
       write(6,'(A)')'  This option retains its legacy short name. Related: --help spectrum.'
       write(6,'(A)')'  Supply a finite positive Gaussian standard deviation.'
      case('build')
       write(6,'(A)')'Build and launch the native VASPBERRY executable.'
       write(6,'(A)')''
       write(6,'(A)')'Intel oneAPI / Intel MPI'
       write(6,'(A)')'  Load the site compiler, Intel MPI and oneMKL environment first.'
       write(6,'(A)')'  make ifx          -> build/vaspberry-ifx'
       write(6,'(A)')'  make ifx-mpi      -> build/vaspberry-ifx-mpi'
       write(6,'(A)')'  mpiexec -n 4 build/vaspberry-ifx-mpi --help task'
       write(6,'(A)')'  Classic: make ifort / make ifort-mpi when installed at the site.'
       write(6,'(A)')''
       write(6,'(A)')'GNU alternatives'
       write(6,'(A)')'  make serial       -> build/vaspberry (vaspberry-gfortran alias)'
       write(6,'(A)')'  make mpi          -> build/vaspberry-mpi'
       write(6,'(A)')''
       write(6,'(A)')'Runtime'
       write(6,'(A)')'  Use matching compiler/MPI libraries and launcher, BLAS/LAPACK and byte RECL.'
       write(6,'(A)')'  Serial and MPI use the same vaspberry.f; help needs no VASP inputs.'
       write(6,'(A)')''
       write(6,'(A)')'Historical serial source'
       write(6,'(A)')'  make serial SERIAL_SOURCE=vaspberry_gfortran_serial.f '//achar(92)
       write(6,'(A)')'    BUILD_DIR=build-legacy'
       write(6,'(A)')'  That source has no modern --task parser or Kubo/spin tasks.'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  docs/BUILD.md'
      case('outputs')
       write(6,'(A)')'Numerical VASPBERRY files can be analysed and plotted independently.'
       write(6,'(A)')''
       write(6,'(A)')'Files -> data -> plots'
       write(6,'(A)')'  Kubo --curvature-csv: omega_z_A2 (Angstrom^2), k/spin IDs and gaps;'
       write(6,'(A)')'    plot maps or path curves. Band CSV also carries band energies.'
       write(6,'(A)')'  Fukui BERRYCURV.dat: plaquette flux/curvature and Chern header;'
       write(6,'(A)')'    plot a BZ map. Legacy coordinates repeat periodically by default.'
       write(6,'(A)')'  Pair CSV: energies and undivided numerators; integrate later with'
       write(6,'(A)')'    occupations/gaps. Pair records alone are not conductivity.'
       write(6,'(A)')'  Spin CSVs: sector summaries, k maps and projected-spin spectrum.'
       write(6,'(A)')'  Z2_FIELD.csv + NFIELD.dat: status/parity diagnostics and n-field map.'
       write(6,'(A)')'  Wavefunction files: real/imaginary amplitudes in a CHGCAR container.'
       write(6,'(A)')''
       write(6,'(A)')'Postprocessing routes'
       write(6,'(A)')'  Standard optical Hall: python3 tools/vaspberry_kubo.py kubo-hall --help'
       write(6,'(A)')'  Saved pair cache: import-pairs, then pair-hall in that same tool.'
       write(6,'(A)')'  Hall tables -> tools/plot_hall.py -> mu-dependent Hall curves.'
       write(6,'(A)')'  Fukui DAT -> tools/plot_berry_curvature.py -> Berry-curvature map.'
       write(6,'(A)')'  Python is optional for reading native tables in Origin or another tool.'
       write(6,'(A)')''
       write(6,'(A)')'Conditions and file policy'
       write(6,'(A)')'  Kubo/spin CSVs require new files; parent directories must exist.'
       write(6,'(A)')'  Legacy DAT and wavefunction files may be replaced; use a fresh cwd.'
       write(6,'(A)')'  Require final PASS where provided; incomplete files are not results.'
       write(6,'(A)')'  Transport needs a full BZ mesh/weights; a path is for plotting only.'
       write(6,'(A)')'  CSV readers may skip # comments, but preserve metadata with analysis.'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  docs/OUTPUT_FORMAT.md | docs/HANDS_ON.md'
      case('examples')
       write(6,'(A)')'Copy native commands, then adapt paths, bands and meshes to your data.'
       write(6,'(A)')''
       write(6,'(A)')'Setup'
       write(6,'(A)')'  Run in a fresh result directory. Set an absolute executable path:'
       write(6,'(A)')'    VB=/absolute/path/to/build/vaspberry-ifx'
       write(6,'(A)')'  For MPI, choose the MPI executable and prefix with mpiexec -n N.'
       write(6,'(A)')''
       write(6,'(A)')'Chern number'
       write(6,'(A)')'  "$VB" --task chern --wavecar /data/WAVECAR '//achar(92)
       write(6,'(A)')'    --mesh 12,12 --bands 1:18'
       write(6,'(A)')''
       write(6,'(A)')'Standard Kubo'
       write(6,'(A)')'  "$VB" --task kubo --input-dir /data/optics '//achar(92)
       write(6,'(A)')'    --bands 1:18 --curvature-csv KUBO_WAVEDER.csv'
       write(6,'(A)')'  Needs same-run WAVECAR, WAVEDER, INCAR and OUTCAR.'
       write(6,'(A)')''
       write(6,'(A)')'Explicit canonical pair approximation'
       write(6,'(A)')'  "$VB" --task kubo-pairs --kubo-source wavecar '//achar(92)
       write(6,'(A)')'    --wavecar /data/WAVECAR --pairs-csv PAIRS.csv'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  Backslashes continue commands; each option takes one separate value.'
       write(6,'(A)')'  --help TASK gives specific requirements, units and output filenames.'
       write(6,'(A)')'  docs/NATIVE_COMMANDS.md | docs/HANDS_ON.md | examples/features/'
      end select
      end subroutine vaspberry_help_topic

      subroutine vaspberry_help_legacy_topic(topic)
      implicit none
      character*(*) topic
      select case(trim(topic))
      case('overview')
       write(6,'(A)')'Historical serial source: vaspberry_gfortran_serial.f'
       write(6,'(A)')''
       write(6,'(A)')'Usage and navigation'
       write(6,'(A)')'  This build accepts legacy calculation flags, not modern --task options.'
       write(6,'(A)')'  Help: --help [TOPIC] or -h [TOPIC], as a separate command.'
       write(6,'(A)')'    --help task       supported legacy calculations'
       write(6,'(A)')'    --help bands      legacy band selection'
       write(6,'(A)')'    --help options    supported option topics'
       write(6,'(A)')'    --help all        original complete legacy reference'
       write(6,'(A)')'  Input defaults to WAVECAR. Help does not read input or write results.'
       write(6,'(A)')'  Example: vaspberry -f WAVECAR -kx 12 -ky 12 -ii 1 -if 18'
       write(6,'(A)')'  No Kubo or projected-spin tasks are implemented in this historical source.'
       write(6,'(A)')'  For current functionality build vaspberry.f; see --help build.'
       write(6,'(A)')'  Reference: docs/NATIVE_COMMANDS.md (current source) and --help all.'
      case('task')
       write(6,'(A)')'This historical build uses legacy flags, not --task NAME.'
       write(6,'(A)')''
       write(6,'(A)')'Usage and navigation'
       write(6,'(A)')'    chern         default Berry-flux/Chern-number calculation'
       write(6,'(A)')'    z2            -z2 1, full even Gamma-centered SOC mesh'
       write(6,'(A)')'    optical       -cd 1, selected-band circular selectivity'
       write(6,'(A)')'    wavefunction  -wf N, existing real-space reconstruction'
       write(6,'(A)')'    velocity      -vel 1, historical velocity diagnostic'
       write(6,'(A)')'  No Kubo, projected-spin or optical-spectrum task is implemented here.'
       write(6,'(A)')'  Use --help TOPIC for supported tasks; --help all gives the old reference.'
       write(6,'(A)')'  Build the current vaspberry.f for the modern interface.'
      case('options')
       write(6,'(A)')'Supported help topics for this historical serial source:'
       write(6,'(A)')''
       write(6,'(A)')'Usage and navigation'
       write(6,'(A)')'    bands mesh input-dir wavecar spinor output'
       write(6,'(A)')'    wavefunction-band kpoint real-grid imaginary'
       write(6,'(A)')'  Calculation syntax uses legacy flags: -f, -kx/-ky, -ii/-if/-is, -s, -o,'
       write(6,'(A)')'  -wf, -k, -ng, -im. Except --input-dir, long names here are help topics only.'
       write(6,'(A)')'  More help: task build outputs examples all.'
       write(6,'(A)')'  Current-source Kubo/spin/optical-spectrum options are unavailable here.'
      case('chern')
       write(6,'(A)')'Legacy Chern-number/Berry-flux calculation on a full periodic 2D mesh.'
       write(6,'(A)')''
       write(6,'(A)')'Run'
       write(6,'(A)')'    vaspberry -f WAVECAR -kx 12 -ky 12 -ii 1 -if 18'
       write(6,'(A)')''
       write(6,'(A)')'Inputs, outputs and conditions'
       write(6,'(A)')'  Fukui-Hatsugai-Suzuki (FHS), also called the Fukui method.'
       write(6,'(A)')'  Writes BERRYCURV.dat; select an isolated band/subspace and check'
       write(6,'(A)')'    convergence.'
       write(6,'(A)')'  Modern --task/--mesh/--bands syntax is unavailable in this source.'
      case('z2')
       write(6,'(A)')'Legacy interface to the shared Fukui-Hatsugai (FH) n-field Z2 method.'
       write(6,'(A)')''
       write(6,'(A)')'Run'
       write(6,'(A)')'    vaspberry -f WAVECAR -z2 1 -kx 12 -ky 12 -ii 1 -if 10'
       write(6,'(A)')''
       write(6,'(A)')'Inputs, outputs and conditions'
       write(6,'(A)')'  Requires nonmagnetic gapped time-reversal-symmetric SOC, two components,'
       write(6,'(A)')'  full even Gamma-centered mesh Nx,Ny>=4, kz=0 modulo G, four TRIM once.'
       write(6,'(A)')'  Select all occupied bands with even rank, and keep empty source bands.'
       write(6,'(A)')'  Use ICHARG=11, ISYM=-1. Report only a final result_status=PASS CSV.'
       write(6,'(A)')'  See --help all for full requirements. Modern --task syntax is unavailable.'
      case('optical')
       write(6,'(A)')'Historical selected-transition circular selectivity: -cd 1.'
       write(6,'(A)')''
       write(6,'(A)')'Run'
       write(6,'(A)')'    vaspberry -f WAVECAR -cd 1 -kx 12 -ky 12 -ii 18 -if 19'
       write(6,'(A)')''
       write(6,'(A)')'Inputs, outputs and conditions'
       write(6,'(A)')'  The endpoints select initial/final bands. Default output:'
       write(6,'(A)')'    CIRC_DICHROISM.dat.'
       write(6,'(A)')'  This source has no spectrum or angle-option implementation.'
       write(6,'(A)')'  Use the current vaspberry.f for modern optical functionality.'
      case('wavefunction')
       write(6,'(A)')'Historical reconstruction: -wf BAND -k KPOINT -ng NX,NY,NZ -im 0|1.'
       write(6,'(A)')''
       write(6,'(A)')'Usage'
       write(6,'(A)')'  Historical reconstruction: -wf BAND -k KPOINT -ng NX,NY,NZ -im 0|1.'
       write(6,'(A)')''
       write(6,'(A)')'Inputs, outputs and conditions'
       write(6,'(A)')'  Matching POSCAR and EIGENVAL are required in cwd; --input-dir does not'
       write(6,'(A)')'    relocate them.'
       write(6,'(A)')'  Output holds wavefunction amplitudes in a CHGCAR container, not charge'
       write(6,'(A)')'    density.'
       write(6,'(A)')'  Select the appropriate Gamma index; defaults are first k, inferred grid,'
       write(6,'(A)')'  and no imaginary companion. Writes PARCHG-W-K###-E###-SPIN# files.'
       write(6,'(A)')'  Use a fresh directory. Modern --task syntax is unavailable here.'
      case('velocity')
       write(6,'(A)')'Historical velocity diagnostic: -vel 1 -is BAND.'
       write(6,'(A)')''
       write(6,'(A)')'Usage'
       write(6,'(A)')'  Historical velocity diagnostic: -vel 1 -is BAND.'
       write(6,'(A)')''
       write(6,'(A)')'Inputs, outputs and conditions'
       write(6,'(A)')'  Writes VEL_EXPT.dat; choose a single band and fresh working directory.'
       write(6,'(A)')'  This retained historical implementation is not the corrected current'
       write(6,'(A)')'  velocity task. Prefer the current vaspberry.f for quantitative use.'
      case('bands')
       write(6,'(A)')'Legacy selection: -is N or -ii FIRST -if LAST (one-based inclusive).'
       write(6,'(A)')''
       write(6,'(A)')'Meaning and compatibility'
       write(6,'(A)')'  Example: -ii 1 -if 18. Optical -cd 1 uses initial/final endpoints.'
       write(6,'(A)')'  The default boundary is inferred from source occupations.'
       write(6,'(A)')'  Modern --bands is not accepted by this historical calculation parser.'
      case('mesh')
       write(6,'(A)')'Legacy mesh dimensions: -kx NX -ky NY; defaults 12 and 12.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning and compatibility'
       write(6,'(A)')'  They must match the stored full periodic 2D mesh; no k points are generated.'
       write(6,'(A)')'  Z2 needs an even Gamma-centered mesh with Nx,Ny>=4.'
       write(6,'(A)')'  Modern --mesh is not accepted by this historical calculation parser.'
      case('input-dir')
       write(6,'(A)')'--input-dir DIR : directory for an omitted WAVECAR; default cwd.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning and compatibility'
       write(6,'(A)')'  DIR must exist. An explicit -f PATH overrides only WAVECAR, in any order.'
       write(6,'(A)')'  Relative -f PATH remains relative to invocation cwd; outputs stay in cwd.'
       write(6,'(A)')'  Example: vaspberry --input-dir /data/run -kx 12 -ky 12 -ii 1 -if 18'
       write(6,'(A)')'  Historical path fields support 75 characters; overlong/empty paths fail.'
      case('wavecar')
       write(6,'(A)')'-f PATH selects the input WAVECAR; default WAVECAR.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning and compatibility'
       write(6,'(A)')'  Requires full-complex RTAG=45200 data and matching reader/writer RECL units.'
       write(6,'(A)')'  GNU builds use byte RECL. This source stores most paths in 75 characters.'
       write(6,'(A)')'  Modern --wavecar is not accepted by this historical calculation parser.'
      case('spinor')
       write(6,'(A)')'-s auto|1|2 : infer or assert the component count, default auto.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning and compatibility'
       write(6,'(A)')'  Checks every spin/k layout; not a spin degeneracy factor or SOC switch.'
       write(6,'(A)')'  Modern --spinor is not accepted by this historical calculation parser.'
      case('output')
       write(6,'(A)')'-o LABEL changes legacy output naming, not the working directory.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning and compatibility'
       write(6,'(A)')'  Chern -o sample writes BERRYCURV.sample.dat; names depend on the task.'
       write(6,'(A)')'  Default BERRYCURV.dat, CIRC_DICHROISM.dat, VEL_EXPT.dat or NFIELD.dat'
       write(6,'(A)')'  depends on the selected mode; wavefunction names use k and band indices.'
       write(6,'(A)')'  Use a fresh directory because DAT files may be replaced.'
       write(6,'(A)')'  Modern --output is not accepted by this historical calculation parser.'
      case('wavefunction-band')
       write(6,'(A)')'-wf N selects the one-based wavefunction band for reconstruction.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning and compatibility'
       write(6,'(A)')'  This historical source does not accept --wavefunction-band.'
      case('kpoint')
       write(6,'(A)')'-k N selects the one-based source k index for wavefunction reconstruction.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning and compatibility'
       write(6,'(A)')'  Default 1; use the Gamma point. This source does not accept --kpoint.'
      case('real-grid')
       write(6,'(A)')'-ng NX,NY,NZ sets the real-space reconstruction grid; no spaces.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning and compatibility'
       write(6,'(A)')'  Default is basis-derived. This source does not accept --real-grid.'
      case('imaginary')
       write(6,'(A)')'-im 1 writes the imaginary wavefunction companion; default 0.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning and compatibility'
       write(6,'(A)')'  This historical source does not accept --imaginary.'
      case('build')
       write(6,'(A)')'Build and launch the native VASPBERRY executable.'
       write(6,'(A)')''
       write(6,'(A)')'Intel oneAPI / Intel MPI'
       write(6,'(A)')'  Load the site compiler, Intel MPI and oneMKL environment first.'
       write(6,'(A)')'  make ifx          -> build/vaspberry-ifx'
       write(6,'(A)')'  make ifx-mpi      -> build/vaspberry-ifx-mpi'
       write(6,'(A)')'  mpiexec -n 4 build/vaspberry-ifx-mpi --help task'
       write(6,'(A)')'  Classic: make ifort / make ifort-mpi when installed at the site.'
       write(6,'(A)')''
       write(6,'(A)')'GNU alternatives'
       write(6,'(A)')'  make serial       -> build/vaspberry (vaspberry-gfortran alias)'
       write(6,'(A)')'  make mpi          -> build/vaspberry-mpi'
       write(6,'(A)')''
       write(6,'(A)')'Runtime'
       write(6,'(A)')'  Use matching compiler/MPI libraries and launcher, BLAS/LAPACK and byte RECL.'
       write(6,'(A)')'  Serial and MPI use the same vaspberry.f; help needs no VASP inputs.'
       write(6,'(A)')''
       write(6,'(A)')'Historical serial source'
       write(6,'(A)')'  make serial SERIAL_SOURCE=vaspberry_gfortran_serial.f '//achar(92)
       write(6,'(A)')'    BUILD_DIR=build-legacy'
       write(6,'(A)')'  That source has no modern --task parser or Kubo/spin tasks.'
       write(6,'(A)')''
       write(6,'(A)')'More'
       write(6,'(A)')'  docs/BUILD.md'
      case('outputs')
       write(6,'(A)')'Historical output is legacy DAT/wavefunction data, not Kubo/spin CSV.'
       write(6,'(A)')''
       write(6,'(A)')'Meaning and compatibility'
       write(6,'(A)')'  Use a fresh working directory; files may be overwritten.'
       write(6,'(A)')'  For Z2, require the final PASS status in Z2_FIELD.csv.'
       write(6,'(A)')'  Build the current vaspberry.f for reusable Kubo and projected-spin output.'
       write(6,'(A)')'  See --help build and docs/OUTPUT_FORMAT.md.'
      case('examples')
       write(6,'(A)')'Historical serial source examples (use your actual bands and mesh):'
       write(6,'(A)')''
       write(6,'(A)')'Meaning and compatibility'
       write(6,'(A)')'    vaspberry -f WAVECAR -kx 12 -ky 12 -ii 1 -if 18'
       write(6,'(A)')'    vaspberry -f WAVECAR -z2 1 -kx 12 -ky 12 -ii 1 -if 10'
       write(6,'(A)')'    vaspberry -f WAVECAR -cd 1 -kx 12 -ky 12 -ii 18 -if 19'
       write(6,'(A)')'  Use a fresh working directory. See --help TOPIC for task prerequisites.'
       write(6,'(A)')'  Modern --task commands require building the current vaspberry.f.'
      end select
      end subroutine vaspberry_help_legacy_topic
