program test_z2_helpers
  implicit none

  integer :: i, ierr, info
  integer :: gs(3), gt(3), gtback(3), expected(3)
  real(8) :: ks(3), kt(3), trims(3,4), phase, smin, pi
  real(8) :: norm_value, norm_reference, norm_legacy
  complex(4) :: raw(4)
  complex(8) :: up, down, theta_up, theta_down
  complex(8) :: twice_up, twice_down
  complex(8) :: s1(1,1), s2(2,2), s18(18,18)

  pi = 4.0_8 * atan(1.0_8)
  trims = reshape((/ &
       0.0_8, 0.0_8, 0.0_8, &
       0.5_8, 0.0_8, 0.0_8, &
       0.0_8, 0.5_8, 0.0_8, &
       0.5_8, 0.5_8, 0.0_8 /), shape(trims))
  gs = (/ 2, -3, 1 /)

  do i = 1, 4
    ks = trims(:,i)
    kt = ks
    call z2_map_gvector(gt, gs, ks, kt, 1, ierr)
    call require(ierr == 0, "TRIM mapping returned an error")
    expected = -gs + nint(-ks-kt)
    call require(all(gt == expected), "TRIM reciprocal shift is wrong")
    call require(maxval(abs(dble(gt)+kt+dble(gs)+ks)) < 1.0e-13_8, &
                 "mapped physical plane wave is not time reversed")
    call z2_map_gvector(gtback, gt, kt, ks, 1, ierr)
    call require(ierr == 0 .and. all(gtback == gs), &
                 "reciprocal mapping is not an involution")
  end do

  ks = (/ 0.25_8, 0.125_8, 0.0_8 /)
  kt = ks + (/ 1.0_8, -1.0_8, 0.0_8 /)
  call z2_map_gvector(gt, gs, ks, kt, 0, ierr)
  call require(ierr == 0, "direct reciprocal mapping failed")
  call require(all(gt == gs + (/ -1, 1, 0 /)), &
               "direct reciprocal shift is wrong")

  kt = ks + (/ 0.1_8, 0.0_8, 0.0_8 /)
  call z2_map_gvector(gt, gs, ks, kt, 1, ierr)
  call require(ierr == 1, "noninteger reciprocal mapping was accepted")

  up = cmplx(0.2_8, -0.7_8, kind=8)
  down = cmplx(-0.3_8, 0.4_8, kind=8)
  call z2_apply_theta(up, down, theta_up, theta_down)
  call z2_apply_theta(theta_up, theta_down, twice_up, twice_down)
  call require(abs(twice_up+up) < 1.0e-14_8, "Theta^2 up != -1")
  call require(abs(twice_down+down) < 1.0e-14_8, "Theta^2 down != -1")

  raw(1) = cmplx(0.1, 0.2, kind=4)
  raw(2) = cmplx(-0.31, 0.47, kind=4)
  raw(3) = cmplx(0.12345, -0.54321, kind=4)
  raw(4) = cmplx(-0.77, -0.09, kind=4)
  call z2_raw_norm(raw, size(raw), norm_value)
  norm_reference = 0.0_8
  norm_legacy = 0.0_8
  do i = 1, size(raw)
    norm_reference = norm_reference + dble(real(raw(i)))**2 &
                    + dble(aimag(raw(i)))**2
    norm_legacy = norm_legacy + dble(abs(raw(i)))**2
  end do
  call require(abs(norm_value-norm_reference) < 1.0e-15_8, &
               "complex*8 norm was not accumulated in double precision")
  call require(abs(norm_legacy-norm_reference) > 1.0e-10_8, &
               "norm regression fixture does not expose old rounding")

  s1(1,1) = cmplx(3.0_8, 4.0_8, kind=8)
  call z2_link_svd_phase(s1, 1, phase, smin, info)
  call require(info == 0, "rank-one link factorization failed")
  call require(abs(smin-5.0_8) < 1.0e-14_8, "rank-one singular value")
  call require(abs(phase-atan2(4.0_8,3.0_8)) < 1.0e-14_8, &
               "rank-one determinant phase")

  s2 = cmplx(0.0_8, 0.0_8, kind=8)
  s2(1,1) = cmplx(1.0_8, 0.0_8, kind=8)
  s2(2,2) = cmplx(1.0e-8_8, 0.0_8, kind=8)
  call z2_link_svd_phase(s2, 2, phase, smin, info)
  call require(info == 0, "ill-conditioned link factorization failed")
  call require(abs(smin-1.0e-8_8) < 1.0e-18_8, &
               "hidden small singular value was not measured")

  s2 = cmplx(0.0_8, 0.0_8, kind=8)
  s2(1,2) = cmplx(1.0_8, 0.0_8, kind=8)
  s2(2,1) = cmplx(1.0_8, 0.0_8, kind=8)
  call z2_link_svd_phase(s2, 2, phase, smin, info)
  call require(info == 0, "pivoted link factorization failed")
  call require(abs(abs(phase)-pi) < 1.0e-13_8, &
               "LU pivot parity was omitted from determinant phase")

  s18 = cmplx(0.0_8, 0.0_8, kind=8)
  do i = 1, 18
    s18(i,i) = cmplx(0.1_8, 0.0_8, kind=8)
  end do
  call z2_link_svd_phase(s18, 18, phase, smin, info)
  call require(info == 0, "high-rank link factorization failed")
  call require(abs(smin-0.1_8) < 1.0e-13_8, &
               "well-conditioned high-rank link was misclassified")

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not. condition) then
      write(*,'(A)') "FAIL: "//trim(message)
      stop 1
    end if
  end subroutine require

end program test_z2_helpers
