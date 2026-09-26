! Check the LP64 complex LAPACK interface used by the native Z2 code.
! diag(3i, 2) has singular values 3 and 2, including a complex-valued input.
program test_lapack_runtime
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, parameter :: dp = kind(0.0d0)
  integer :: info, lwork
  complex(dp) :: a(2, 2), u(1, 1), vt(1, 1), query(1)
  complex(dp), allocatable :: work(:)
  real(dp) :: singular(2), rwork(10)
  external :: zgesvd

  if (storage_size(info) /= 32) error stop 'LP64 requires 32-bit LAPACK integers'
  a = cmplx(0.0_dp, 0.0_dp, kind=dp)
  a(1, 1) = cmplx(0.0_dp, 3.0_dp, kind=dp)
  a(2, 2) = cmplx(2.0_dp, 0.0_dp, kind=dp)
  call zgesvd('N', 'N', 2, 2, a, 2, singular, u, 1, vt, 1, query, -1, rwork, info)
  if (info /= 0) error stop 'ZGESVD workspace query failed'
  if (.not. ieee_is_finite(real(query(1), dp))) error stop 'Invalid ZGESVD workspace'
  lwork = int(real(query(1), dp))
  if (lwork < 6 .or. lwork > 1000000) error stop 'Unexpected ZGESVD workspace'
  allocate(work(lwork))
  call zgesvd('N', 'N', 2, 2, a, 2, singular, u, 1, vt, 1, work, lwork, rwork, info)
  if (info /= 0) error stop 'ZGESVD calculation failed'
  if (.not. all(ieee_is_finite(singular))) error stop 'Non-finite singular values'
  if (maxval(abs(singular - [3.0_dp, 2.0_dp])) > 1.0d-12) &
    error stop 'Incorrect complex ZGESVD singular values'
  print *, 'PASS: LP64 complex ZGESVD singular values =', singular
end program test_lapack_runtime
