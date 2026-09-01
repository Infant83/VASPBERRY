program test_z2_legacy_atomic
  implicit none

  integer :: unit_no
  integer :: isp, ispin, ispinor, irecl, nk, nkx, nky, nband
  integer :: kperiod, nini, nmax, kext, iz2, icd, iz, ivel, nprocs
  real(8) :: ecut, cell_area, rvari, rvari2
  real(8) :: b1(3), b2(3), b3(3)
  real(8) :: xrecivec(3,4), xrecilat(3,4), xrvari(4)
  real(8) :: rvari3(4), rvari4(4)
  character(75) :: fonameo, foname, filename
  logical :: exists, deletable

  isp = 1
  ispin = 1
  ispinor = 2
  irecl = 24
  nk = 1
  nkx = 1
  nky = 1
  nband = 2
  kperiod = 1
  nini = 1
  nmax = 2
  kext = 0
  iz2 = 1
  icd = 0
  iz = 1
  ivel = 0
  nprocs = 1
  ecut = 1.0_8
  cell_area = 1.0_8
  rvari = 0.0_8
  rvari2 = 0.0_8
  b1 = (/ 1.0_8, 0.0_8, 0.0_8 /)
  b2 = (/ 0.0_8, 1.0_8, 0.0_8 /)
  b3 = (/ 0.0_8, 0.0_8, 1.0_8 /)
  xrecivec = 0.0_8
  xrecilat = 0.0_8
  xrvari = 0.0_8
  rvari3 = 0.0_8
  rvari4 = 0.0_8
  foname = "NFIELD"
  fonameo = " "
  filename = "synthetic-WAVECAR"

  open(newunit=unit_no, file="Z2_FIELD.invalid.csv", &
       status="replace", action="write")
  write(unit_no,'(A)') "# schema=VASPBERRY_Z2_FIELD"
  write(unit_no,'(A)') "# result_status=INCOMPLETE"
  close(unit_no)
  open(newunit=unit_no, file="Z2_FIELD.tmp", &
       status="replace", action="write")
  write(unit_no,'(A)') "# schema=VASPBERRY_Z2_FIELD"
  write(unit_no,'(A)') "# result_status=PASS"
  close(unit_no)

  call write_result(isp,ispin,ispinor,fonameo,foname,filename, &
       irecl,ecut,nk,nkx,nky,nband,b1,b2,b3,kperiod,cell_area, &
       nini,nmax,xrecivec,xrecilat,kext,xrvari,rvari,rvari2, &
       rvari3,rvari4,iz2,icd,iz,ivel,nprocs)

  inquire(file="NFIELD.dat", exist=exists)
  call require(exists, "atomic legacy target was not published")
  inquire(file="NFIELD.tmp", exist=exists)
  call require(.not. exists, "legacy temporary file was not consumed")
  call z2_legacy_deletable("NFIELD.dat", deletable)
  call require(deletable, "published legacy output lacks ownership markers")
  inquire(file="Z2_FIELD.csv", exist=exists)
  call require(.not. exists, "PASS field was published before final commit")

  call finalize_z2_pass(foname)
  inquire(file="Z2_FIELD.csv", exist=exists)
  call require(exists, "final PASS field was not published")
  inquire(file="Z2_FIELD.invalid.csv", exist=exists)
  call require(.not. exists, "INCOMPLETE sentinel survived PASS commit")
  inquire(file="Z2_FIELD.tmp", exist=exists)
  call require(.not. exists, "staged PASS field survived final commit")

  call delete_z2_file("Z2_FIELD.csv")
  call delete_z2_legacy_file("NFIELD.dat")

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not. condition) then
      write(*,'(A)') "FAIL: "//trim(message)
      stop 1
    end if
  end subroutine require

end program test_z2_legacy_atomic
