program test_mpi_runtime
  implicit none
  include 'mpif.h'

  integer :: ierr, rank, nranks, i
  integer :: field_ok
  real(8), allocatable :: local_field(:), reduced_field(:)

  call MPI_INIT(ierr)
  call require_mpi(ierr == MPI_SUCCESS, "MPI_INIT failed")
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call require_mpi(ierr == MPI_SUCCESS, "MPI_COMM_RANK failed")
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
  call require_mpi(ierr == MPI_SUCCESS, "MPI_COMM_SIZE failed")
  call require_mpi(nranks >= 2, "run this test with at least two ranks")

  allocate(local_field(nranks), reduced_field(nranks))
  local_field = 0.0_8
  reduced_field = 0.0_8
  local_field(rank + 1) = real(rank + 1, kind=8)

  call MPI_REDUCE(local_field, reduced_field, nranks, &
                  MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call require_mpi(ierr == MPI_SUCCESS, "MPI_DOUBLE_PRECISION reduce failed")

  field_ok = 0
  if (rank == 0) then
    do i = 1, nranks
      call require_mpi(abs(reduced_field(i) - real(i, kind=8)) < 1.0e-14_8, &
                       "disjoint field reduction changed a value")
    end do
    field_ok = 1
  end if

  call MPI_BCAST(field_ok, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call require_mpi(ierr == MPI_SUCCESS, "MPI_INTEGER broadcast failed")
  call require_mpi(field_ok == 1, "root field status was not broadcast")

  deallocate(local_field, reduced_field)
  call MPI_FINALIZE(ierr)
  if (ierr /= MPI_SUCCESS) stop 1

contains

  subroutine require_mpi(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (.not. condition) then
      write(*,'(A,I0,A,A)') "rank ", rank, ": ", trim(message)
      call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
  end subroutine require_mpi

end program test_mpi_runtime
