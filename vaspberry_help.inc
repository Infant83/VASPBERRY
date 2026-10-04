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
       write(6,'(A)')'Usage: vaspberry --help [TOPIC]   (or -h [TOPIC])'
       write(6,'(A)')'Help is a separate command; it never opens WAVECAR or writes results.'
       write(6,'(A)')''
       write(6,'(A)')'Start here:'
       write(6,'(A)')'  --help task       choose a calculation'
       write(6,'(A)')'  --help kubo       Berry curvature from the Kubo formula'
       write(6,'(A)')'  --help bands      select one band or a band range'
       write(6,'(A)')'  --help options    find an option topic'
       write(6,'(A)')'  --help examples   copy a complete native command'
       write(6,'(A)')'  --help build      Intel/GNU serial and MPI builds'
       write(6,'(A)')'  --help outputs    use the numerical files in other tools'
       write(6,'(A)')'  --help all        complete legacy reference (also: --help legacy)'
       write(6,'(A)')''
       write(6,'(A)')'Intel MPI: make ifx-mpi'
       write(6,'(A)')'  mpiexec -n 4 build/vaspberry-ifx-mpi --help task'
       write(6,'(A)')'Example (join lines using the shown backslash):'
       write(6,'(A)')'  mpiexec -n 4 build/vaspberry-ifx-mpi --task kubo --bands 1:18 '//achar(92)
       write(6,'(A)')'    --curvature-csv KUBO_WAVEDER.csv'
       write(6,'(A)')'Standard Kubo needs same-run WAVEDER/WAVECAR/INCAR/OUTCAR.'
       write(6,'(A)')'Inputs default to cwd; --input-dir DIR changes omitted source filenames.'
       write(6,'(A)')'Use --kubo-source wavecar only to opt into the warned approximation.'
       write(6,'(A)')'Use a fresh working directory for each calculation.'
       write(6,'(A)')'Topic names accept -- prefixes: --help --bands'
       write(6,'(A)')'Further guidance: docs/NATIVE_COMMANDS.md and docs/HANDS_ON.md'
      case('task')
       write(6,'(A)')'Select one calculation with --task NAME; default: chern.'
       write(6,'(A)')'  chern          Chern number and Berry flux on a periodic 2D mesh'
       write(6,'(A)')'  spin-chern     projected-spin sector Chern numbers on a full mesh'
       write(6,'(A)')'  spin-kubo      approximate spin-sector local Berry curvature'
       write(6,'(A)')'  z2             FH n-field Z2 invariant on an even full mesh'
       write(6,'(A)')'  kubo           WAVEDER selected/occupied curvature; default PAW optical source'
       write(6,'(A)')'  kubo-integral  full-mesh Kubo-formula integration diagnostic'
       write(6,'(A)')'  kubo-pairs     approximate pair export; requires --kubo-source wavecar'
       write(6,'(A)')'  optical        selected-transition circular selectivity'
       write(6,'(A)')'  spectrum       broadened left/right optical transition spectra'
       write(6,'(A)')'  wavefunction   real-space wavefunction reconstruction'
       write(6,'(A)')'  velocity       canonical-momentum velocity expectation'
       write(6,'(A)')'kubo-line is a synonym of kubo. See --help NAME for each task.'
       write(6,'(A)')'Fukui-Hatsugai-Suzuki (FHS) is called the Fukui method here.'
       write(6,'(A)')'chern and spin-chern use it; Z2 uses the distinct FH n-field method.'
       write(6,'(A)')'Explicit --task rejects conflicting legacy calculation flags.'
      case('options')
       write(6,'(A)')'Each calculation option takes one separate value; no --name=value.'
       write(6,'(A)')'Read a topic with --help NAME or --help --NAME:'
       write(6,'(A)')'  task  input-dir  wavecar  waveder  incar  kubo-source  bands  mesh  spinor  output'
       write(6,'(A)')'  curvature-csv  pairs-csv  per-band'
       write(6,'(A)')'  spin-axis  outcar  energy-gap-tol  spin-gap-tol  sum-bands'
       write(6,'(A)')'  wavefunction-band  kpoint  real-grid  imaginary  theta  phi'
       write(6,'(A)')'Optical controls retain legacy names: ien fen nediv sigma.'
       write(6,'(A)')'  --help ien       photon-energy window'
       write(6,'(A)')'  --help sigma     Gaussian broadening'
       write(6,'(A)')'More help: build outputs examples all'
       write(6,'(A)')'List values have no spaces: --mesh 12,12 --bands 1:18.'
       write(6,'(A)')'Task examples use illustrative bands/meshes; match your own VASP data.'
      case('chern')
       write(6,'(A)')'Chern number and plaquette Berry flux via the Fukui method (FHS).'
       write(6,'(A)')'Requires the full periodic uniform 2D mesh stored in WAVECAR.'
       write(6,'(A)')'  vaspberry --task chern --wavecar WAVECAR '//achar(92)
       write(6,'(A)')'    --mesh 12,12 --bands 1:18'
       write(6,'(A)')'Replace the example mesh and bands with your VASP calculation.'
       write(6,'(A)')'--bands N selects one band; FIRST:LAST selects a subspace.'
       write(6,'(A)')'Use an isolated subspace and check mesh/link convergence.'
       write(6,'(A)')'Default occupied boundary is inferred from source occupations;'
       write(6,'(A)')'explicit band selection is clearer for reproducible calculations.'
       write(6,'(A)')'Writes BERRYCURV.dat: plaquette values and a Chern-number header.'
       write(6,'(A)')'--output LABEL changes the legacy naming label, not the directory.'
       write(6,'(A)')'Related: --help mesh, --help bands, --help outputs.'
      case('spin-chern')
       write(6,'(A)')'Projected-spin sector Chern numbers via the Fukui method (FHS).'
       write(6,'(A)')'  vaspberry --task spin-chern --wavecar WAVECAR --outcar OUTCAR '//achar(92)
       write(6,'(A)')'    --mesh 12,12 --bands 1:8 --spin-axis z'
       write(6,'(A)')'Needs a full periodic 2D mesh, two-component spinors, and matching'
       write(6,'(A)')'OUTCAR spin-frame metadata. The selected energy group must be isolated;'
       write(6,'(A)')'the projected-Pauli spectrum must stay separated from zero.'
       write(6,'(A)')'Writes SPIN_CHERN.csv, SPIN_BERRY.csv, SPIN_SPECTRUM.csv.'
       write(6,'(A)')'--output NAME replaces the default SPIN prefix; files must be new.'
       write(6,'(A)')'Reports C_plus, C_minus, C_charge and C_spin=(C_plus-C_minus)/2.'
       write(6,'(A)')'These use the pseudo-wavefunction metric, not a PAW-complete metric.'
       write(6,'(A)')'Check mesh, energy-gap and spin-gap convergence before interpretation.'
       write(6,'(A)')'Related: spin-axis outcar energy-gap-tol spin-gap-tol.'
       write(6,'(A)')'Reference: docs/SPIN_CHERN.md.'
      case('spin-kubo')
       write(6,'(A)')'Approximate projected-spin sector Berry curvature from a Kubo formula.'
       write(6,'(A)')'Requires --kubo-source wavecar; WAVEDER lacks the required spin matrices.'
       write(6,'(A)')'  vaspberry --task spin-kubo --kubo-source wavecar --wavecar WAVECAR --outcar OUTCAR '//achar(92)
       write(6,'(A)')'    --bands 1:10 --spin-axis z'
       write(6,'(A)')'The default evaluates supplied k points, including a path.'
       write(6,'(A)')'An explicit --mesh NX,NY additionally requests full-mesh integrals;'
       write(6,'(A)')'it does not interpolate or create missing points.'
       write(6,'(A)')'--sum-bands N restricts intermediate source bands to 1:N (default all).'
       write(6,'(A)')'Writes SPIN_KUBO.csv and SPIN_KUBO_SPECTRUM.csv; with an explicit'
       write(6,'(A)')'full mesh also SPIN_KUBO_INTEGRAL.csv. --output sets the prefix.'
       write(6,'(A)')'Uses a canonical-momentum response tangent with source-energy denominators'
       write(6,'(A)')'and projected-spin derivatives in the pseudo metric. No integer invariant.'
       write(6,'(A)')'This is not the conventional spin-current Kubo spin Hall conductivity.'
       write(6,'(A)')'Requires an isolated selected group, open spin gap and matching OUTCAR.'
       write(6,'(A)')'Related: sum-bands spin-axis outcar energy-gap-tol spin-gap-tol.'
       write(6,'(A)')'Reference: docs/SPIN_KUBO.md.'
      case('z2')
       write(6,'(A)')'Fukui-Hatsugai (FH) n-field Z2 invariant; distinct from FHS Chern numbers.'
       write(6,'(A)')'  vaspberry --task z2 --wavecar WAVECAR --mesh 12,12 --bands 1:10'
       write(6,'(A)')'Requires a gapped, nonmagnetic, time-reversal-symmetric SOC calculation,'
       write(6,'(A)')'ISPIN=1 with two-component spinors, and all occupied bands 1:NE.'
       write(6,'(A)')'NE must be even and NBANDS must exceed NE.'
       write(6,'(A)')'Use ICHARG=11 with ISYM=-1 and a full even Gamma-centered 2D mesh:'
       write(6,'(A)')'Nx,Ny >= 4, kz=0 modulo G, each of the four 2D TRIM exactly once.'
       write(6,'(A)')'PASS writes Z2_FIELD.csv and NFIELD.dat. Report z2_invariant only'
       write(6,'(A)')'when result_status=PASS and the two half-BZ parities agree.'
       write(6,'(A)')'Rejected runs may leave Z2_FIELD.invalid.csv; check denser even meshes.'
       write(6,'(A)')'Legacy selector: -z2 1. Half-mesh -hf 1 is rejected.'
       write(6,'(A)')'Reference: docs/Z2_FUKUI_HATSUGAI.md.'
      case('kubo')
       write(6,'(A)')'Kubo Berry curvature from PAW optical connections; default source: waveder.'
       write(6,'(A)')'  vaspberry --task kubo --input-dir /data/optics --bands 1:18 '//achar(92)
       write(6,'(A)')'    --curvature-csv KUBO.csv'
       write(6,'(A)')'WAVECAR/WAVEDER/INCAR/OUTCAR default to cwd, or --input-dir DIR.'
       write(6,'(A)')'--wavecar only overrides WAVECAR; it never redirects auxiliary inputs.'
       write(6,'(A)')'Overrides: --wavecar PATH --waveder PATH --incar PATH --outcar PATH.'
       write(6,'(A)')'Explicit relative file paths remain relative to invocation cwd, in any option order.'
       write(6,'(A)')'Supported: standard VASP 5.4.4 LOPTICS, non-PEAD, longitudinal static branch.'
       write(6,'(A)')'Explicit --bands 31,33:34 selects a bundle independently of all NBANDS partners.'
       write(6,'(A)')'Omit --bands for verified fixed occupied T0 filling. Internal trace degeneracies are allowed.'
       write(6,'(A)')'The selected boundary must not split full-source transitive 0.002 eV clusters.'
       write(6,'(A)')'--per-band 1 with explicit --bands writes only resolved individual bands.'
       write(6,'(A)')'Each required pair needs a stored orientation and strict consistency if both exist;'
       write(6,'(A)')'no square padding, invented missing elements, or second energy denominator.'
       write(6,'(A)')'Default trace CSV: KUBO_WAVEDER.csv; --curvature-csv changes the filename.'
       write(6,'(A)')'CSV records omega_z_A2 in Angstrom^2, exact band IDs and optical source metadata.'
       write(6,'(A)')'Source paths give point curvature only. Explicit --mesh validates a full grid'
       write(6,'(A)')'and adds unrounded integration metadata. Selected integrals are geometric contributions.'
       write(6,'(A)')'Only the omitted-band verified occupied integral labels total physical Hall.'
       write(6,'(A)')'Native integral sums both ISPIN=2 channels. Python kubo-hall selects --spin'
       write(6,'(A)')'(default 1); sum separately calculated spin 1 and 2 for collinear charge.'
       write(6,'(A)')'Missing/unsupported input stops. No automatic WAVECAR fallback.'
       write(6,'(A)')'Modern standard Kubo rejects wavefunction/optical/folding controls it cannot use:'
       write(6,'(A)')'--kpoint/--real-grid/--imaginary, --theta/--phi and their legacy aliases,'
       write(6,'(A)')'-ishift, -ien/-fen/-nediv/-sigma, -kp/-skp. No silent k selection or broadening.'
       write(6,'(A)')'Explicit approximation: --task kubo --kubo-source wavecar --bands 18'
       write(6,'(A)')'WAVECAR uses canonical momentum without full PAW/nonlocal velocity terms.'
       write(6,'(A)')'It permits isolated single/per-band or selected traces with gaps >1e-5 eV.'
       write(6,'(A)')'Canonical trace default is KUBO.csv. --per-band 1 selects separate bands.'
       write(6,'(A)')'Both routes require new CSV paths and complete finite output with PASS.'
       write(6,'(A)')'Check k/NBANDS convergence. kmesh cannot restore missing operator terms.'
       write(6,'(A)')'kubo-line is a synonym; legacy -kubo selectors obey the same source barrier.'
       write(6,'(A)')'Reference: docs/WAVEDER_KUBO_PROTOCOL.md; related: kubo-source curvature-csv.'
      case('kubo-integral')
       write(6,'(A)')'Full-mesh Kubo integral; default source is standard WAVEDER.'
       write(6,'(A)')'  vaspberry --task kubo-integral --mesh 12,12 --bands 1:18 '//achar(92)
       write(6,'(A)')'    --curvature-csv KUBO.csv'
       write(6,'(A)')'Requires explicit full --mesh and same-run WAVEDER/WAVECAR/INCAR/OUTCAR.'
       write(6,'(A)')'Explicit --bands gives unit-weight selected Chern contribution, not total Hall.'
       write(6,'(A)')'A path is not an integration mesh. The finite-grid integral is never rounded.'
       write(6,'(A)')'For occupation-weighted Hall tables use python3 tools/vaspberry_kubo.py kubo-hall.'
       write(6,'(A)')'The historical canonical diagnostic requires --kubo-source wavecar.'
       write(6,'(A)')'Operator, source validation and gap assumptions: --help kubo.'
      case('kubo-pairs')
       write(6,'(A)')'Canonical-momentum pair export requires explicit --kubo-source wavecar.'
       write(6,'(A)')'  vaspberry --task kubo-pairs --kubo-source wavecar --pairs-csv PAIRS.csv'
       write(6,'(A)')'The default WAVEDER route does not promise all unordered source-band pairs.'
       write(6,'(A)')'Use standard --task kubo for selected or automatic occupied curvature.'
       write(6,'(A)')'WAVECAR approximation: omit --bands; all n<m pairs are exported.'
       write(6,'(A)')'Writes energies, fractional k coordinates and three Cartesian numerators'
       write(6,'(A)')'(yz,zx,xy), in eV^2 Angstrom^2. No gap division or occupation weighting.'
       write(6,'(A)')'Source occupations may vary with k. No physical filling is assumed.'
       write(6,'(A)')'Not compatible with --curvature-csv or --per-band. Use a new file.'
       write(6,'(A)')'Numerators omit PAW/nonlocal velocity terms; a warning identifies this choice.'
       write(6,'(A)')'Final result_status=PASS marks complete export; pairs are not conductivity.'
       write(6,'(A)')'Postprocess full-BZ data: python3 tools/vaspberry_kubo.py import-pairs, then pair-hall.'
       write(6,'(A)')'Reference: docs/OUTPUT_FORMAT.md and docs/WAVEDER_KUBO_PROTOCOL.md.'
      case('optical')
       write(6,'(A)')'Circular optical selectivity for one initial-to-final transition.'
       write(6,'(A)')'  vaspberry --task optical --wavecar WAVECAR '//achar(92)
       write(6,'(A)')'    --mesh 12,12 --bands 18:19 '//achar(92)
       write(6,'(A)')'    --theta 0 --phi 0'
       write(6,'(A)')'The two endpoints select the initial and final bands, not a trace bundle.'
       write(6,'(A)')'Default endpoints are inferred VBM/CBM; explicit choices are reproducible.'
       write(6,'(A)')'Angles are in degrees; (theta,phi)=(0,0) means light along z.'
       write(6,'(A)')'Writes BERRYCURV.dat by default (the retained legacy filename).'
       write(6,'(A)')'--output LABEL produces LABEL.dat. Use a fresh working directory.'
       write(6,'(A)')'Uses pseudo-wavefunction momentum matrix elements.'
       write(6,'(A)')'For a broadened spectrum use --help spectrum. Legacy selector: -cd 1.'
      case('spectrum')
       write(6,'(A)')'Broadened left/right optical transition spectra (legacy -cd 2).'
       write(6,'(A)')'  vaspberry --task spectrum --wavecar WAVECAR '//achar(92)
       write(6,'(A)')'    --mesh 12,12 --bands 1:20 '//achar(92)
       write(6,'(A)')'    -ien 1 -fen 3 -nediv 201 -sigma 0.05'
       write(6,'(A)')'Uses occupied bands 1:VBM and empty bands VBM+1:LAST.'
       write(6,'(A)')'LAST is the upper --bands endpoint; FIRST does not limit occupied bands.'
       write(6,'(A)')'Photon energies and Gaussian broadening are in eV.'
       write(6,'(A)')'Defaults: -ien 0, -fen 10, -nediv 1000, -sigma 0.01.'
       write(6,'(A)')'Writes OPT_TRANS_RATE_LEFT.dat and OPT_TRANS_RATE_RIGHT.dat'
       write(6,'(A)')'plus per-k _KP# files, in arbitrary units; zero-energy rows are omitted.'
       write(6,'(A)')'--output LABEL instead uses CIRC_DICHROISM_W.LABEL_LEFT/RIGHT.dat.'
       write(6,'(A)')'This is a broadened transition diagnostic, not a calibrated absolute rate.'
       write(6,'(A)')'Related: bands theta phi ien fen nediv sigma.'
       write(6,'(A)')'Reference: examples/features/circular-dichroism/.'
      case('wavefunction')
       write(6,'(A)')'Reconstruct the existing Gamma-point real-space wavefunction output.'
       write(6,'(A)')'  vaspberry --task wavefunction --wavecar WAVECAR '//achar(92)
       write(6,'(A)')'    --wavefunction-band 18 --kpoint 24 '//achar(92)
       write(6,'(A)')'    --real-grid 24,24,64 --imaginary 1'
       write(6,'(A)')'The example k index 24 is Gamma only in the included MoS2 path WAVECAR;'
       write(6,'(A)')'choose the Gamma index in your own file. Indices are one-based.'
       write(6,'(A)')'Matching POSCAR/EIGENVAL should be present in the working directory.'
       write(6,'(A)')'Writes CHGCAR-style PARCHG-W-K###-E###-SPIN# files and, if requested,'
       write(6,'(A)')'the -IM-SPIN# companion. --output does not rename these files.'
       write(6,'(A)')'Existing files may be overwritten; use a fresh working directory.'
       write(6,'(A)')'Defaults: first k point, basis-derived real grid, no imaginary companion.'
       write(6,'(A)')'Legacy selector: -wf N; related: wavefunction-band kpoint real-grid.'
       write(6,'(A)')'Reference: examples/features/wavefunction/.'
      case('velocity')
       write(6,'(A)')'Canonical-momentum velocity expectation for a selected band.'
       write(6,'(A)')'  vaspberry --task velocity --wavecar WAVECAR --bands 18'
       write(6,'(A)')'Use a single band: this routine uses the first selected band.'
       write(6,'(A)')'Writes x/y components in m/s to VEL_EXPT.dat.'
       write(6,'(A)')'--output LABEL uses VEL_EXPT.LABEL.dat; scalar ISPIN=2 adds UP/DN.'
       write(6,'(A)')'This is the bare-momentum expectation in the pseudo-wavefunction metric,'
       write(6,'(A)')'not a full PAW/nonlocal velocity operator. Legacy selector: -vel 1.'
       write(6,'(A)')'Use a fresh working directory; existing DAT files may be replaced.'
      case('bands')
       write(6,'(A)')'--bands N | FIRST:LAST | N,M:P : one-based inclusive selection (no duplicates).'
       write(6,'(A)')'Examples: --bands 18 ; --bands 1:18 ; WAVEDER kubo also accepts 31,33:34.'
       write(6,'(A)')'chern: one band or the full selected link-variable subspace.'
       write(6,'(A)')'WAVEDER kubo: explicit selection is geometric; omitted bands verifies occupied T0.'
       write(6,'(A)')'Explicit wavecar: N is individual; FIRST:LAST is the group trace.'
       write(6,'(A)')'spin-chern/spin-kubo: the selected isolated energy subspace.'
       write(6,'(A)')'z2: all occupied bands 1:NE, with even rank and NBANDS>NE.'
       write(6,'(A)')'optical: initial/final endpoints. spectrum: LAST caps empty bands;'
       write(6,'(A)')'occupied bands still start at 1 regardless of FIRST.'
       write(6,'(A)')'velocity: only the first selected band; choose one explicitly.'
       write(6,'(A)')'kubo-pairs: omit this option; every source-band pair is exported.'
       write(6,'(A)')'Default boundary is inferred from source occupations, not fixed at 18.'
       write(6,'(A)')'Legacy selectors: -is N or -ii FIRST -if LAST.'
       write(6,'(A)')'Prefer one --bands selection. Repeated/mixed selectors print a warning.'
       write(6,'(A)')'A later --bands replaces the earlier selection; later legacy -ii/-if/-is'
       write(6,'(A)')'keeps historical endpoint/singleton precedence after a scalar/range selector.'
       write(6,'(A)')'Legacy -ii/-if/-is cannot modify an earlier comma-list; a later --bands'
       write(6,'(A)')'can replace legacy selectors. Inspect effective band IDs before interpretation.'
       write(6,'(A)')'With --task kubo these follow the new range rule; pure legacy -kubo'
       write(6,'(A)')'with -ii/-if keeps individual bands. Explicit --bands uses the new rule.'
       write(6,'(A)')'This selector does not set the number of intermediate states;'
       write(6,'(A)')'--sum-bands is a separate spin-kubo-only convergence control.'
      case('mesh')
       write(6,'(A)')'--mesh NX,NY : dimensions of the stored uniform periodic 2D mesh.'
       write(6,'(A)')'Example: --mesh 12,12 ; -kx 12 -ky 12 sets the same mesh dimensions.'
       write(6,'(A)')'For WAVEDER or spin Kubo, only explicit --mesh requests integration;'
       write(6,'(A)')'legacy -kx/-ky alone set sizes and leave point curvature unintegrated.'
       write(6,'(A)')'WAVEDER kubo-integral requires --mesh, not only -kx/-ky.'
       write(6,'(A)')'Default native dimensions are 12,12; specify the actual VASP mesh.'
       write(6,'(A)')'This option does not generate k points, interpolate, or make a path a mesh.'
       write(6,'(A)')'Required for chern, spin-chern and full-mesh integration tasks.'
       write(6,'(A)')'z2 additionally needs an even Gamma-centered mesh, Nx,Ny>=4.'
       write(6,'(A)')'kubo and kubo-pairs accept arbitrary source k points without --mesh.'
       write(6,'(A)')'spin-kubo evaluates points by default; an explicit mesh additionally'
       write(6,'(A)')'requests SPIN_KUBO_INTEGRAL.csv and checks complete mesh geometry.'
       write(6,'(A)')'For transport use a full BZ mesh and verify sampling convergence.'
      case('kubo-source')
       write(6,'(A)')'--kubo-source waveder|wavecar : Kubo operator source; default waveder.'
       write(6,'(A)')'waveder: same-run optical connections; selected geometry or verified occupied trace.'
       write(6,'(A)')'wavecar: explicit canonical-momentum approximation; prints a warning.'
       write(6,'(A)')'Missing WAVEDER or unsupported requests stop; fallback is never automatic.'
       write(6,'(A)')'kubo-pairs and spin-kubo presently require explicit wavecar.'
       write(6,'(A)')'This choice changes the operator, not just input speed or file preference.'
       write(6,'(A)')'Reference: docs/WAVEDER_KUBO_PROTOCOL.md.'
      case('input-dir')
       write(6,'(A)')'--input-dir DIR : base for omitted WAVECAR/WAVEDER/INCAR/OUTCAR filenames.'
       write(6,'(A)')'Default is the invocation working directory. DIR must exist and be a directory.'
       write(6,'(A)')'Example: --task kubo --input-dir /data/optics --curvature-csv result.csv.'
       write(6,'(A)')'Each explicit --wavecar/--waveder/--incar/--outcar overrides only that file.'
       write(6,'(A)')'Explicit relative file paths are relative to cwd, never relative to DIR.'
       write(6,'(A)')'Option order does not change precedence. --wavecar never redirects other inputs.'
       write(6,'(A)')'Does not change cwd, outputs, or unrelated input lists. Empty/overlong paths fail.'
       write(6,'(A)')'Standard Kubo prints all selected input paths before opening the source files.'
      case('waveder')
       write(6,'(A)')'--waveder PATH : optical connection input for --kubo-source waveder.'
       write(6,'(A)')'Default is WAVEDER under --input-dir (default cwd), independent of --wavecar.'
       write(6,'(A)')'Explicit relative PATH is relative to cwd; keep same-run INCAR and OUTCAR.'
       write(6,'(A)')'Supported rectangular coverage is sufficient; do not pad to NBANDS squared.'
       write(6,'(A)')'Standard VASP 5.4.4 longitudinal LOPTICS branch only; PAW terms are already'
       write(6,'(A)')'included upstream. The stored connection already contains the gap denominator.'
      case('incar')
       write(6,'(A)')'--incar PATH : same-run INCAR for the standard WAVEDER Kubo route.'
       write(6,'(A)')'Default is INCAR under --input-dir (default cwd), independent of --wavecar.'
       write(6,'(A)')'Explicit relative PATH is relative to cwd. It validates the producer branch.'
       write(6,'(A)')'Use LOPTICS=T, LPEAD=F, LNABLA=F, NSW=0, ISYM=-1, LREAL=F.'
       write(6,'(A)')'Keep the actual input; do not synthesize settings to make validation pass.'
      case('wavecar')
       write(6,'(A)')'--wavecar PATH : overrides only WAVECAR; default is --input-dir/WAVECAR.'
       write(6,'(A)')'--input-dir defaults to cwd. Explicit relative PATH is relative to cwd.'
       write(6,'(A)')'Changing --wavecar never changes WAVEDER/INCAR/OUTCAR defaults.'
       write(6,'(A)')'Example: --wavecar /data/project/WAVECAR ; legacy equivalent: -f PATH.'
       write(6,'(A)')'Quote paths with spaces. Prefer short paths (native fields: 256 chars).'
       write(6,'(A)')'Reads full-complex single-precision coefficients with RTAG=45200.'
       write(6,'(A)')'GNU and supplied Intel targets use byte RECL. The native reader detects'
       write(6,'(A)')'byte or four-byte-word producer layouts and reports physical byte stride.'
       write(6,'(A)')'Ambiguous, truncated or invalid layouts fail before allocation; no conversion.'
       write(6,'(A)')'--spinor auto checks scalar/two-component layout throughout the file.'
       write(6,'(A)')'Help never attempts to open or inspect this input.'
      case('spinor')
       write(6,'(A)')'--spinor auto|1|2 : infer or assert the WAVECAR component layout.'
       write(6,'(A)')'Default auto checks full-complex coefficient counts for every spin/k point.'
       write(6,'(A)')'1 or 2 asserts that the file has that detected number of components.'
       write(6,'(A)')'This is not a spin degeneracy multiplier and does not enable SOC.'
       write(6,'(A)')'ISPIN counts channels; it does not by itself identify spinor components.'
       write(6,'(A)')'Spin-sector and Z2 tasks require two-component spinors.'
       write(6,'(A)')'Use matching OUTCAR for spin-frame conventions. Legacy equivalent: -s.'
      case('output')
       write(6,'(A)')'--output LABEL : legacy output label or spin-task prefix, not a directory.'
       write(6,'(A)')'Examples for LABEL=sample:'
       write(6,'(A)')'  chern       BERRYCURV.sample.dat'
       write(6,'(A)')'  optical     sample.dat'
       write(6,'(A)')'  spectrum    CIRC_DICHROISM_W.sample_LEFT/RIGHT.dat'
       write(6,'(A)')'  velocity    VEL_EXPT.sample.dat'
       write(6,'(A)')'  spin-chern  sample_CHERN/BERRY/SPECTRUM.csv'
       write(6,'(A)')'  spin-kubo   sample_KUBO.csv and sample_KUBO_SPECTRUM.csv'
       write(6,'(A)')'Spin-task default prefix is SPIN. Scalar ISPIN=2 DAT files add UP/DN.'
       write(6,'(A)')'--curvature-csv and --pairs-csv specify their own exact CSV paths.'
       write(6,'(A)')'Wavefunction names are set by k/band; --output does not rename them.'
       write(6,'(A)')'Use a fresh working directory; DAT outputs can replace existing files.'
       write(6,'(A)')'Legacy equivalent: -o LABEL. See --help outputs for data usage.'
      case('curvature-csv')
       write(6,'(A)')'--curvature-csv PATH : exact Kubo curvature CSV filename.'
       write(6,'(A)')'Standard WAVEDER default is KUBO_WAVEDER.csv. Omitted --bands uses schema'
       write(6,'(A)')'VASPBERRY_WAVEDER_KUBO_OCCUPIED_V1; omega_z_A2 and external gap.'
       write(6,'(A)')'Explicit selection: WAVEDER_KUBO_BUNDLE_V1 or WAVEDER_KUBO_BAND_V1 (VASPBERRY_ prefix).'
       write(6,'(A)')'Example: --task kubo --kubo-source wavecar --bands 18 --curvature-csv BAND18.csv.'
       write(6,'(A)')'Used only with kubo or kubo-integral (legacy -kubo 2 or 1).'
       write(6,'(A)')'The parent directory must exist; PATH and PATH.partial must be new.'
       write(6,'(A)')'Trace V2/band V3 includes declared counts and terminal result_status=PASS.'
       write(6,'(A)')'Band rows include spin,k,band,energy,min_gap,omega_z_A2.'
       write(6,'(A)')'A multi-band range exports selected-group trace and external-gap data.'
       write(6,'(A)')'omega_z_A2 is normalized Berry curvature in Angstrom^2.'
       write(6,'(A)')'Preserve comment metadata with the table; see docs/OUTPUT_FORMAT.md.'
       write(6,'(A)')'Not compatible with --pairs-csv. Legacy equivalent: -kubo_csv PATH.'
      case('pairs-csv')
       write(6,'(A)')'--pairs-csv PATH : exact reusable native pair-export CSV filename.'
       write(6,'(A)')'Use with --task kubo-pairs --kubo-source wavecar; legacy -kubo_pairs PATH.'
       write(6,'(A)')'The parent directory must exist; PATH and PATH.partial must be new.'
       write(6,'(A)')'PATH.partial is renamed only after all records are written and closed.'
       write(6,'(A)')'Omit --bands, --curvature-csv and --per-band.'
       write(6,'(A)')'All unordered source-band pairs are exported, with no gap division or'
       write(6,'(A)')'occupation weighting. Require the final result_status=PASS marker.'
       write(6,'(A)')'See --help kubo-pairs for quantities, units and postprocessing.'
      case('per-band')
       write(6,'(A)')'--per-band 0|1 : individual Kubo output when 1; default 0.'
       write(6,'(A)')'WAVEDER requires explicit --bands and resolved full-source 0.002 eV clusters.'
       write(6,'(A)')'Every required partner among all NBANDS must be stored in either orientation.'
       write(6,'(A)')'WAVEDER always writes CSV (default KUBO_WAVEDER.csv); no occupation weights.'
       write(6,'(A)')'Requires --task kubo, kubo-line or kubo-integral.'
       write(6,'(A)')'Example: --task kubo --kubo-source wavecar --bands 1:18 --per-band 1 --curvature-csv BANDS.csv.'
       write(6,'(A)')'Without it, a selected multi-band range is the selected-space trace.'
       write(6,'(A)')'For explicit wavecar the following historical rules apply:'
       write(6,'(A)')'Individual bands require all gaps >1e-5 eV at every source k/spin,'
       write(6,'(A)')'including gaps to other selected bands. Degeneracy fails before output.'
       write(6,'(A)')'No band range is enlarged automatically. A trace requires only external'
       write(6,'(A)')'gaps >1e-5 eV; internal transitions are excluded before division.'
       write(6,'(A)')'The trace writes CSV only, one row per source k/spin; default KUBO.csv.'
       write(6,'(A)')'For the full source space its trace is zero in the truncated basis;'
       write(6,'(A)')'this is not a convergence result. These canonical modes apply no occupations.'
       write(6,'(A)')'For individual output CSV is opt-in; existing DAT output remains.'
       write(6,'(A)')'Omitted bands retain occupied inference and individual output; use --bands.'
       write(6,'(A)')'Rejected with pair export, spin-sector, FHS and all other tasks.'
      case('bundle')
       write(6,'(A)')'--bundle and -kubo_bundle were removed; old commands fail explicitly.'
       write(6,'(A)')'Use --task kubo --bands FIRST:LAST for the selected-space trace.'
       write(6,'(A)')'Defaults: WAVEDER KUBO_WAVEDER.csv; explicit wavecar trace KUBO.csv.'
       write(6,'(A)')'For either source use explicit --bands with --per-band 1 for resolved bands.'
       write(6,'(A)')'WAVEDER also accepts disjoint selections such as --bands 31,33:34.'
      case('spin-axis')
       write(6,'(A)')'--spin-axis x|y|z|A,B,C : Cartesian analysis axis for spin-sector tasks.'
       write(6,'(A)')'Default z; for example --spin-axis 0,0,1.'
       write(6,'(A)')'For A,B,C supply a finite unit vector (norm squared within 1e-8 of 1).'
       write(6,'(A)')'Used only by spin-chern and spin-kubo; it does not rotate VASP input.'
       write(6,'(A)')'Matching OUTCAR gives the source SAXIS-to-Cartesian spin frame.'
       write(6,'(A)')'Changing the analysis axis may close the projected-spin gap.'
      case('outcar')
       write(6,'(A)')'--outcar PATH : same-run optical producer record, or spin-frame metadata.'
       write(6,'(A)')'Defaults to OUTCAR under --input-dir (default cwd) for Kubo and spin tasks.'
       write(6,'(A)')'Explicit relative PATH is relative to cwd; --wavecar never redirects it.'
       write(6,'(A)')'Needed for the printed SAXIS-to-Cartesian spin-frame rotation and'
       write(6,'(A)')'matching source dimensions/lattice. Preserve the WAVECAR/OUTCAR pair.'
       write(6,'(A)')'Not used to replace or infer wavefunction coefficients.'
      case('energy-gap-tol')
       write(6,'(A)')'--energy-gap-tol VALUE : positive selected-to-excluded gap tolerance.'
       write(6,'(A)')'Used by spin-chern and spin-kubo; units eV, default 1e-8.'
       write(6,'(A)')'The selected energy subspace must stay isolated from excluded states.'
       write(6,'(A)')'This numerical cutoff is not broadening or evidence of convergence.'
       write(6,'(A)')'Charge thresholds: standard WAVEDER >0.002 eV; wavecar >1e-5 eV.'
      case('spin-gap-tol')
       write(6,'(A)')'--spin-gap-tol VALUE : positive projected-Pauli distance from zero.'
       write(6,'(A)')'Used by spin-chern and spin-kubo; dimensionless, default 1e-6.'
       write(6,'(A)')'It is not an energy broadening or an eV gap.'
       write(6,'(A)')'A closed spin gap makes the axis-dependent sector decomposition unresolved.'
       write(6,'(A)')'Inspect the spin spectrum and convergence instead of relaxing it blindly.'
      case('sum-bands')
       write(6,'(A)')'--sum-bands N : spin-kubo-only intermediate-state convergence control.'
       write(6,'(A)')'Include source bands 1:N outside the selected group in the Kubo sum.'
       write(6,'(A)')'Default is all stored source bands. N must include the selected group'
       write(6,'(A)')'and cannot exceed source NBANDS. VASP eigenstates remain fixed.'
       write(6,'(A)')'Example: --task spin-kubo --kubo-source wavecar --bands 1:10 --sum-bands 20.'
       write(6,'(A)')'This does not truncate source WAVECAR, change --bands, or apply to'
       write(6,'(A)')'charge kubo, kubo-pairs, or spin-chern. Compare increasing cutoffs.'
      case('wavefunction-band')
       write(6,'(A)')'--wavefunction-band N : positive one-based band to reconstruct.'
       write(6,'(A)')'Required with --task wavefunction. Legacy equivalent: -wf N.'
       write(6,'(A)')'This is the reconstruction selector, separate from --bands.'
       write(6,'(A)')'Related: --help wavefunction, --help kpoint, --help real-grid.'
      case('kpoint')
       write(6,'(A)')'--kpoint N : positive one-based source k index for wavefunction output.'
       write(6,'(A)')'Default 1. Choose the Gamma point for the existing reconstruction routine.'
       write(6,'(A)')'This option does not select k points for curvature/transport tasks.'
       write(6,'(A)')'Legacy equivalent: -k N. Related: --help wavefunction.'
      case('real-grid')
       write(6,'(A)')'--real-grid NX,NY,NZ : positive real-space wavefunction grid dimensions.'
       write(6,'(A)')'Default uses the basis-derived grid. Example: --real-grid 24,24,64.'
       write(6,'(A)')'Used for wavefunction reconstruction, not reciprocal-space sampling.'
       write(6,'(A)')'Legacy equivalent: -ng NX,NY,NZ. Related: --help wavefunction.'
      case('imaginary')
       write(6,'(A)')'--imaginary 0|1 : write an imaginary wavefunction companion if 1.'
       write(6,'(A)')'Default 0. Used with wavefunction reconstruction.'
       write(6,'(A)')'Legacy equivalent: -im 1. The -IM-SPIN# file accompanies the real output.'
       write(6,'(A)')'Related: --help wavefunction.'
      case('theta')
       write(6,'(A)')'--theta DEG : light incidence polar angle in degrees, default 0.'
       write(6,'(A)')'Used by optical/spectrum; legacy equivalent: -theta DEG.'
       write(6,'(A)')'(theta,phi)=(0,0) is along z; (90,0) along x; (90,90) along y.'
       write(6,'(A)')'This is a light-direction control, not the spin analysis axis.'
      case('phi')
       write(6,'(A)')'--phi DEG : light incidence azimuthal angle in degrees, default 0.'
       write(6,'(A)')'Used by optical/spectrum; legacy equivalent: -phi DEG.'
       write(6,'(A)')'(theta,phi)=(0,0) is along z; (90,0) along x; (90,90) along y.'
       write(6,'(A)')'Related: --help theta and --help optical.'
      case('ien')
       write(6,'(A)')'-ien EV : initial photon energy for the optical spectrum, default 0 eV.'
       write(6,'(A)')'Together with -fen EV (default 10 eV), sets the energy window.'
       write(6,'(A)')'Example: --task spectrum --bands 1:20 -ien 1 -fen 3 -nediv 201.'
       write(6,'(A)')'These controls retain their legacy short-option names.'
       write(6,'(A)')'Related: --help spectrum, --help nediv, --help sigma.'
      case('fen')
       write(6,'(A)')'-fen EV : final photon energy for the optical spectrum, default 10 eV.'
       write(6,'(A)')'Use together with -ien EV (default 0 eV) and -nediv N (default 1000).'
       write(6,'(A)')'Units eV; this option retains its legacy short name.'
       write(6,'(A)')'Related: --help spectrum and --help ien.'
      case('nediv')
       write(6,'(A)')'-nediv N : number of photon-energy points for the optical spectrum.'
       write(6,'(A)')'Default 1000; combine with -ien and -fen, for example -nediv 201.'
       write(6,'(A)')'This option retains its legacy short name.'
       write(6,'(A)')'Related: --help spectrum and --help sigma.'
      case('sigma')
       write(6,'(A)')'-sigma EV : Gaussian broadening for optical transition spectra.'
       write(6,'(A)')'Default 0.01 eV; for example -sigma 0.05.'
       write(6,'(A)')'It does not change spin-gap/energy-gap tolerances or Kubo denominators.'
       write(6,'(A)')'This option retains its legacy short name. Related: --help spectrum.'
      case('build')
       write(6,'(A)')'Intel oneAPI Fortran and Intel MPI (load the site environment first):'
       write(6,'(A)')'  make ifx          -> build/vaspberry-ifx'
       write(6,'(A)')'  make ifx-mpi      -> build/vaspberry-ifx-mpi'
       write(6,'(A)')'  mpiexec -n 4 build/vaspberry-ifx-mpi --help task'
       write(6,'(A)')'GNU alternatives:'
       write(6,'(A)')'  make serial       -> build/vaspberry (vaspberry-gfortran alias)'
       write(6,'(A)')'  make mpi          -> build/vaspberry-mpi'
       write(6,'(A)')'Intel Classic: make ifort / make ifort-mpi when installed at the site.'
       write(6,'(A)')'Current serial and MPI binaries use the same vaspberry.f source.'
       write(6,'(A)')'Use matching compiler/MPI runtime, BLAS/LAPACK, and byte RECL.'
       write(6,'(A)')'Historical serial source:'
       write(6,'(A)')'  make serial SERIAL_SOURCE=vaspberry_gfortran_serial.f '//achar(92)
       write(6,'(A)')'    BUILD_DIR=build-legacy'
       write(6,'(A)')'That historical source has no modern --task interface or Kubo/spin tasks.'
       write(6,'(A)')'Reference: docs/BUILD.md; help requires no VASP input files.'
      case('outputs')
       write(6,'(A)')'VASPBERRY writes numerical files for independent analysis and plots.'
       write(6,'(A)')'--curvature-csv PATH: normalized Kubo omega_z_A2, gaps and k indices.'
       write(6,'(A)')'Individual-band CSV includes band energies; bundle CSV has no energy column.'
       write(6,'(A)')'Pair CSV: undivided numerators and energies for later occupation/gap sums.'
       write(6,'(A)')'Fukui-method DAT: plaquette flux/curvature with a Chern-number header.'
       write(6,'(A)')'Spin CSVs: sector summaries, k maps and projected-spin spectrum.'
       write(6,'(A)')'For Origin, spreadsheets or plotting scripts, skip # metadata lines to'
       write(6,'(A)')'read CSV column names, but preserve the metadata with your analysis.'
       write(6,'(A)')'Kubo and spin CSVs require new files; parent directories must exist.'
       write(6,'(A)')'Legacy DAT and wavefunction files may be overwritten: use a fresh cwd.'
       write(6,'(A)')'Require final PASS status where provided; incomplete files are not results.'
       write(6,'(A)')'Transport needs a full BZ mesh and k weights; a path is for plotting only.'
       write(6,'(A)')'Reusable pair workflow: python3 tools/vaspberry_kubo.py import-pairs/pair-hall;'
       write(6,'(A)')'plot the resulting scans with tools/plot_hall.py or your preferred tool.'
       write(6,'(A)')'Reference: docs/OUTPUT_FORMAT.md and docs/HANDS_ON.md.'
      case('examples')
       write(6,'(A)')'These commands use example band/mesh choices; match your own VASP data.'
       write(6,'(A)')'Run from a fresh result directory with the executable/input paths set.'
       write(6,'(A)')'Set the absolute executable path:'
       write(6,'(A)')'  VB=/absolute/path/to/build/vaspberry-ifx'
       write(6,'(A)')'  "$VB" --task chern --wavecar /data/WAVECAR --mesh 12,12 --bands 1:18'
       write(6,'(A)')'  "$VB" --task kubo --input-dir /data --bands 1:18 '//achar(92)
       write(6,'(A)')'    --curvature-csv KUBO.csv'
       write(6,'(A)')'  "$VB" --task kubo-pairs --kubo-source wavecar --wavecar /data/WAVECAR --pairs-csv PAIRS.csv'
       write(6,'(A)')'Backslashes continue commands onto the next line; options take one value.'
       write(6,'(A)')'For MPI, use the MPI executable with mpiexec -n N before the command.'
       write(6,'(A)')'For a copyable public-WAVECAR example see docs/NATIVE_COMMANDS.md.'
       write(6,'(A)')'More runnable workflows: docs/HANDS_ON.md and examples/features/.'
       write(6,'(A)')'Use --help TASK for task-specific inputs, units and output names.'
      end select
      end subroutine vaspberry_help_topic

      subroutine vaspberry_help_legacy_topic(topic)
      implicit none
      character*(*) topic
      select case(trim(topic))
      case('overview')
       write(6,'(A)')'Historical serial source: vaspberry_gfortran_serial.f'
       write(6,'(A)')'This build accepts legacy calculation flags, not modern --task options.'
       write(6,'(A)')'Help: --help [TOPIC] or -h [TOPIC], as a separate command.'
       write(6,'(A)')'  --help task       supported legacy calculations'
       write(6,'(A)')'  --help bands      legacy band selection'
       write(6,'(A)')'  --help options    supported option topics'
       write(6,'(A)')'  --help all        original complete legacy reference'
       write(6,'(A)')'Input defaults to WAVECAR. Help does not read input or write results.'
       write(6,'(A)')'Example: vaspberry -f WAVECAR -kx 12 -ky 12 -ii 1 -if 18'
       write(6,'(A)')'No Kubo or projected-spin tasks are implemented in this historical source.'
       write(6,'(A)')'For current functionality build vaspberry.f; see --help build.'
       write(6,'(A)')'Reference: docs/NATIVE_COMMANDS.md (current source) and --help all.'
      case('task')
       write(6,'(A)')'This historical build uses legacy flags, not --task NAME.'
       write(6,'(A)')'  chern         default Berry-flux/Chern-number calculation'
       write(6,'(A)')'  z2            -z2 1, full even Gamma-centered SOC mesh'
       write(6,'(A)')'  optical       -cd 1, selected-band circular selectivity'
       write(6,'(A)')'  wavefunction  -wf N, existing real-space reconstruction'
       write(6,'(A)')'  velocity      -vel 1, historical velocity diagnostic'
       write(6,'(A)')'No Kubo, projected-spin or optical-spectrum task is implemented here.'
       write(6,'(A)')'Use --help TOPIC for supported tasks; --help all gives the old reference.'
       write(6,'(A)')'Build the current vaspberry.f for the modern interface.'
      case('options')
       write(6,'(A)')'Supported help topics for this historical serial source:'
       write(6,'(A)')'  bands mesh input-dir wavecar spinor output'
       write(6,'(A)')'  wavefunction-band kpoint real-grid imaginary'
       write(6,'(A)')'Calculation syntax uses legacy flags: -f, -kx/-ky, -ii/-if/-is, -s, -o,'
       write(6,'(A)')'-wf, -k, -ng, -im. Long option names are help topics only.'
       write(6,'(A)')'More help: task build outputs examples all.'
       write(6,'(A)')'Current-source Kubo/spin/optical-spectrum options are unavailable here.'
      case('chern')
       write(6,'(A)')'Legacy Chern-number/Berry-flux calculation on a full periodic 2D mesh.'
       write(6,'(A)')'  vaspberry -f WAVECAR -kx 12 -ky 12 -ii 1 -if 18'
       write(6,'(A)')'Fukui-Hatsugai-Suzuki (FHS), also called the Fukui method.'
       write(6,'(A)')'Writes BERRYCURV.dat; select an isolated band/subspace and check convergence.'
       write(6,'(A)')'Modern --task/--mesh/--bands syntax is unavailable in this source.'
      case('z2')
       write(6,'(A)')'Legacy interface to the shared Fukui-Hatsugai (FH) n-field Z2 method.'
       write(6,'(A)')'  vaspberry -f WAVECAR -z2 1 -kx 12 -ky 12 -ii 1 -if 10'
       write(6,'(A)')'Requires nonmagnetic gapped time-reversal-symmetric SOC, two components,'
       write(6,'(A)')'full even Gamma-centered mesh Nx,Ny>=4, kz=0 modulo G, four TRIM once.'
       write(6,'(A)')'Select all occupied bands with even rank, and keep empty source bands.'
       write(6,'(A)')'Use ICHARG=11, ISYM=-1. Report only a final result_status=PASS CSV.'
       write(6,'(A)')'See --help all for full requirements. Modern --task syntax is unavailable.'
      case('optical')
       write(6,'(A)')'Historical selected-transition circular selectivity: -cd 1.'
       write(6,'(A)')'  vaspberry -f WAVECAR -cd 1 -kx 12 -ky 12 -ii 18 -if 19'
       write(6,'(A)')'The endpoints select initial/final bands. Default output: CIRC_DICHROISM.dat.'
       write(6,'(A)')'This source has no spectrum or angle-option implementation.'
       write(6,'(A)')'Use the current vaspberry.f for modern optical functionality.'
      case('wavefunction')
       write(6,'(A)')'Historical reconstruction: -wf BAND -k KPOINT -ng NX,NY,NZ -im 0|1.'
       write(6,'(A)')'Select the appropriate Gamma index; defaults are first k, inferred grid,'
       write(6,'(A)')'and no imaginary companion. Writes PARCHG-W-K###-E###-SPIN# files.'
       write(6,'(A)')'Use a fresh directory. Modern --task syntax is unavailable here.'
      case('velocity')
       write(6,'(A)')'Historical velocity diagnostic: -vel 1 -is BAND.'
       write(6,'(A)')'Writes VEL_EXPT.dat; choose a single band and fresh working directory.'
       write(6,'(A)')'This retained historical implementation is not the corrected current'
       write(6,'(A)')'velocity task. Prefer the current vaspberry.f for quantitative use.'
      case('bands')
       write(6,'(A)')'Legacy selection: -is N or -ii FIRST -if LAST (one-based inclusive).'
       write(6,'(A)')'Example: -ii 1 -if 18. Optical -cd 1 uses initial/final endpoints.'
       write(6,'(A)')'The default boundary is inferred from source occupations.'
       write(6,'(A)')'Modern --bands is not accepted by this historical calculation parser.'
      case('mesh')
       write(6,'(A)')'Legacy mesh dimensions: -kx NX -ky NY; defaults 12 and 12.'
       write(6,'(A)')'They must match the stored full periodic 2D mesh; no k points are generated.'
       write(6,'(A)')'Z2 needs an even Gamma-centered mesh with Nx,Ny>=4.'
       write(6,'(A)')'Modern --mesh is not accepted by this historical calculation parser.'
      case('input-dir')
       write(6,'(A)')'--input-dir DIR : directory for an omitted WAVECAR; default cwd.'
       write(6,'(A)')'DIR must exist. An explicit -f PATH overrides only WAVECAR, in any order.'
       write(6,'(A)')'Relative -f PATH remains relative to invocation cwd; outputs stay in cwd.'
       write(6,'(A)')'Example: vaspberry --input-dir /data/run -kx 12 -ky 12 -ii 1 -if 18'
       write(6,'(A)')'Historical path fields support 75 characters; overlong/empty paths fail.'
      case('wavecar')
       write(6,'(A)')'-f PATH selects the input WAVECAR; default WAVECAR.'
       write(6,'(A)')'Requires full-complex RTAG=45200 data and matching reader/writer RECL units.'
       write(6,'(A)')'GNU builds use byte RECL. This source stores most paths in 75 characters.'
       write(6,'(A)')'Modern --wavecar is not accepted by this historical calculation parser.'
      case('spinor')
       write(6,'(A)')'-s auto|1|2 : infer or assert the component count, default auto.'
       write(6,'(A)')'Checks every spin/k layout; not a spin degeneracy factor or SOC switch.'
       write(6,'(A)')'Modern --spinor is not accepted by this historical calculation parser.'
      case('output')
       write(6,'(A)')'-o LABEL changes legacy output naming, not the working directory.'
       write(6,'(A)')'Default BERRYCURV.dat, CIRC_DICHROISM.dat, VEL_EXPT.dat or NFIELD.dat'
       write(6,'(A)')'depends on the selected mode; wavefunction names use k and band indices.'
       write(6,'(A)')'Use a fresh directory because DAT files may be replaced.'
       write(6,'(A)')'Modern --output is not accepted by this historical calculation parser.'
      case('wavefunction-band')
       write(6,'(A)')'-wf N selects the one-based wavefunction band for reconstruction.'
       write(6,'(A)')'This historical source does not accept --wavefunction-band.'
      case('kpoint')
       write(6,'(A)')'-k N selects the one-based source k index for wavefunction reconstruction.'
       write(6,'(A)')'Default 1; use the Gamma point. This source does not accept --kpoint.'
      case('real-grid')
       write(6,'(A)')'-ng NX,NY,NZ sets the real-space reconstruction grid; no spaces.'
       write(6,'(A)')'Default is basis-derived. This source does not accept --real-grid.'
      case('imaginary')
       write(6,'(A)')'-im 1 writes the imaginary wavefunction companion; default 0.'
       write(6,'(A)')'This historical source does not accept --imaginary.'
      case('build')
       write(6,'(A)')'Intel oneAPI Fortran and Intel MPI (load the site environment first):'
       write(6,'(A)')'  make ifx          -> build/vaspberry-ifx'
       write(6,'(A)')'  make ifx-mpi      -> build/vaspberry-ifx-mpi'
       write(6,'(A)')'  mpiexec -n 4 build/vaspberry-ifx-mpi --help task'
       write(6,'(A)')'GNU alternatives:'
       write(6,'(A)')'  make serial       -> build/vaspberry (vaspberry-gfortran alias)'
       write(6,'(A)')'  make mpi          -> build/vaspberry-mpi'
       write(6,'(A)')'Intel Classic: make ifort / make ifort-mpi when installed at the site.'
       write(6,'(A)')'Current serial and MPI binaries use the same vaspberry.f source.'
       write(6,'(A)')'Use matching compiler/MPI runtime, BLAS/LAPACK, and byte RECL.'
       write(6,'(A)')'Historical serial source:'
       write(6,'(A)')'  make serial SERIAL_SOURCE=vaspberry_gfortran_serial.f '//achar(92)
       write(6,'(A)')'    BUILD_DIR=build-legacy'
       write(6,'(A)')'That historical source has no modern --task interface or Kubo/spin tasks.'
       write(6,'(A)')'Reference: docs/BUILD.md; help requires no VASP input files.'
      case('outputs')
       write(6,'(A)')'Historical output is legacy DAT/wavefunction data, not Kubo/spin CSV.'
       write(6,'(A)')'Use a fresh working directory; files may be overwritten.'
       write(6,'(A)')'For Z2, require the final PASS status in Z2_FIELD.csv.'
       write(6,'(A)')'Build the current vaspberry.f for reusable Kubo and projected-spin output.'
       write(6,'(A)')'See --help build and docs/OUTPUT_FORMAT.md.'
      case('examples')
       write(6,'(A)')'Historical serial source examples (use your actual bands and mesh):'
       write(6,'(A)')'  vaspberry -f WAVECAR -kx 12 -ky 12 -ii 1 -if 18'
       write(6,'(A)')'  vaspberry -f WAVECAR -z2 1 -kx 12 -ky 12 -ii 1 -if 10'
       write(6,'(A)')'  vaspberry -f WAVECAR -cd 1 -kx 12 -ky 12 -ii 18 -if 19'
       write(6,'(A)')'Use a fresh working directory. See --help TOPIC for task prerequisites.'
       write(6,'(A)')'Modern --task commands require building the current vaspberry.f.'
      end select
      end subroutine vaspberry_help_legacy_topic
