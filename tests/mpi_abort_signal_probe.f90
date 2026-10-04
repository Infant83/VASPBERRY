program mpi_abort_signal_probe
  use, intrinsic :: iso_fortran_env, only: error_unit
  use, intrinsic :: iso_c_binding, only: c_int
  implicit none
  include 'mpif.h'
  integer :: ierr, rank, code, pause_result
  character(len=32) :: mode
  logical :: abort_here
  interface
    integer(c_int) function c_sleep(seconds) bind(C, name='sleep')
      import c_int
      integer(c_int), value :: seconds
    end function c_sleep
  end interface
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call get_command_argument(1, mode)
  if (trim(mode) == 'success') then
    write(error_unit, '(A,I0)') 'PROBE_SUCCESS rank=', rank
    flush(error_unit)
    call MPI_FINALIZE(ierr)
    stop
  end if
  code = 1
  abort_here = .false.
  select case (trim(mode))
  case ('all')
    abort_here = .true.
  case ('rank0')
    abort_here = rank == 0
  case ('rank1')
    abort_here = rank == 1
  case ('abort7')
    abort_here = .true.
    code = 7
  case default
    write(error_unit, '(A)') 'INVALID_PROBE_MODE'
    flush(error_unit)
    call MPI_ABORT(MPI_COMM_WORLD, 19, ierr)
    stop 19
  end select
  if (abort_here) then
    write(error_unit, '(A,I0,A,I0)') 'PROBE_ABORT rank=', rank, ' code=', code
    flush(error_unit)
    call MPI_ABORT(MPI_COMM_WORLD, code, ierr)
    stop 23 ! MPI_ABORT returning would be an unexpected failure.
  end if
  ! An uninformed peer performs no collective. Abort must terminate it.
  do
    pause_result = c_sleep(1_c_int)
  end do
end program mpi_abort_signal_probe
