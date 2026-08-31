program test_z2_helpers
  implicit none

  integer :: i, ierr, info, file_unit
  integer :: gs(3), gt(3), gtback(3), expected(3)
  real(8) :: ks(3), kt(3), trims(3,4), phase, smin, pi
  real(8) :: norm_value, norm_reference, norm_legacy
  complex(4) :: raw(4)
  complex(8) :: up, down, theta_up, theta_down
  complex(8) :: twice_up, twice_down
  complex(8) :: s1(1,1), s2(2,2), s18(18,18)
  complex(4) :: wave_record(4)
  complex(8) :: coeff_up(3), coeff_down(3)
  integer :: loop_k(5), nplane(5), trim_flag(5), tr_flag(5)
  integer :: basis_max(3)
  real(8) :: loop_coord(3,5), mesh_coord(3,1)
  real(8) :: rb1(3), rb2(3), rb3(3), fixture_ecut
  logical :: same_file, path_ok, deletable, exists
  character(4096) :: cwd

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

  loop_k = 1
  nplane = 4
  trim_flag = 0
  tr_flag = 0
  trim_flag(1) = 1
  basis_max = (/ 1, 0, 0 /)
  loop_coord = 0.0_8
  mesh_coord = 0.0_8
  loop_coord(1,1) = 0.5_8
  mesh_coord(1,1) = 0.5_8
  rb1 = (/ 1.0_8, 0.0_8, 0.0_8 /)
  rb2 = (/ 0.0_8, 1.0_8, 0.0_8 /)
  rb3 = (/ 0.0_8, 0.0_8, 1.0_8 /)
  fixture_ecut = 2.0_8
  wave_record(1) = cmplx(0.25, -0.125, kind=4)
  wave_record(2) = cmplx(-0.375, 0.5, kind=4)
  wave_record(3) = cmplx(0.125, 0.25, kind=4)
  wave_record(4) = cmplx(-0.5, -0.25, kind=4)
  open(unit=10, file="z2_get_state.tmp", status="replace", &
       access="direct", form="unformatted", recl=32, iostat=ierr)
  call require(ierr == 0, "could not create get_z2_state fixture")
  do i = 1, 4
    write(10, rec=i, iostat=ierr) wave_record
    call require(ierr == 0, "could not write get_z2_state fixture")
  end do
  call get_z2_state(coeff_up, coeff_down, 1, 2, loop_k, nplane, &
       2, 1, 1, 3, trim_flag, tr_flag, loop_coord, mesh_coord, &
       basis_max, 2, rb1, rb2, rb3, fixture_ecut)
  call require(abs(coeff_up(1)+conjg(dcmplx(wave_record(3)))) &
               < 1.0e-14_8, "M1 Theta up(G=0 -> G=-1) failed")
  call require(abs(coeff_down(1)-conjg(dcmplx(wave_record(1)))) &
               < 1.0e-14_8, "M1 Theta down(G=0 -> G=-1) failed")
  call require(abs(coeff_up(2)+conjg(dcmplx(wave_record(4)))) &
               < 1.0e-14_8, "M1 Theta up(G=-1 -> G=0) failed")
  call require(abs(coeff_down(2)-conjg(dcmplx(wave_record(2)))) &
               < 1.0e-14_8, "M1 Theta down(G=-1 -> G=0) failed")
  call require(abs(coeff_up(3))+abs(coeff_down(3)) < 1.0e-14_8, &
               "M1 reconstruction populated a plane wave outside cutoff")
  close(10, status="delete", iostat=ierr)
  call require(ierr == 0, "could not remove get_z2_state fixture")

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

  open(newunit=file_unit, file="z2_same_file_target.tmp", &
       status="replace", action="write")
  write(file_unit,'(A)') "identity fixture"
  close(file_unit)
  call z2_paths_same("z2_same_file_target.tmp", &
                     "./z2_same_file_target.tmp", same_file, path_ok)
  call require(path_ok .and. same_file, &
               "relative alias was not recognized as the same file")
  call get_environment_variable("PWD", cwd, status=ierr)
  call require(ierr == 0, "PWD is unavailable for absolute-path test")
  call z2_paths_same("z2_same_file_target.tmp", &
                     trim(cwd)//"/z2_same_file_target.tmp", &
                     same_file, path_ok)
  call require(path_ok .and. same_file, &
               "absolute alias was not recognized as the same file")
  call execute_command_line( &
       "ln -sf z2_same_file_target.tmp z2_same_file_link.tmp", &
       exitstat=ierr)
  call require(ierr == 0, "could not create same-file symlink fixture")
  call z2_paths_same("z2_same_file_target.tmp", &
                     "z2_same_file_link.tmp", same_file, path_ok)
  call require(path_ok .and. same_file, &
               "symlink alias was not recognized as the same file")
  call z2_paths_same("z2_same_file_target.tmp", &
                     "z2_missing_output.tmp", same_file, path_ok)
  call require(path_ok .and. .not. same_file, &
               "missing output path was not safely classified")
  call z2_paths_same("z2_missing_input.tmp", &
                     "z2_same_file_target.tmp", same_file, path_ok)
  call require(.not. path_ok, &
               "unresolvable input path did not fail closed")

  open(newunit=file_unit, file="z2_owned_output.tmp", &
       status="replace", action="write")
  write(file_unit,'(A)') "# schema=VASPBERRY_Z2_FIELD"
  close(file_unit)
  call z2_output_deletable("z2_owned_output.tmp", deletable)
  call require(deletable, "owned Z2 output was not recognized")
  call delete_z2_file("z2_owned_output.tmp")
  inquire(file="z2_owned_output.tmp", exist=exists)
  call require(.not. exists, "owned Z2 output was not deleted")

  open(newunit=file_unit, file="z2_unowned_output.tmp", &
       status="replace", action="write")
  write(file_unit,'(A)') "not a VASPBERRY Z2 output"
  close(file_unit)
  call z2_output_deletable("z2_unowned_output.tmp", deletable)
  call require(.not. deletable, "unowned output was marked deletable")
  inquire(file="z2_unowned_output.tmp", exist=exists)
  call require(exists, "unowned output was modified")

  call execute_command_line( &
       "rm -f z2_same_file_link.tmp z2_same_file_target.tmp " // &
       "z2_unowned_output.tmp", exitstat=ierr)
  call require(ierr == 0, "could not clean path-safety fixtures")

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
