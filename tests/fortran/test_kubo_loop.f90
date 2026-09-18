! Synthetic complex plane-wave eigenvectors; calls the real production loop.
program test_kubo_loop
  implicit none
  integer :: spinor, np, nband, nk, nbmax(3), nplist(1), i
  real(8) :: b1(3),b2(3),b3(3),wk(3,1),lower(1),upper(1),unit,expected
  real(8) :: hbar,c,em,metertoang
  complex(4) :: coeffn(6),coeffm(6),phase
  nband=2;nk=1;nbmax=1
  b1=[1d0,0d0,0d0];b2=[0d0,1d0,0d0];b3=[0d0,0d0,1d0];wk=0d0
  ! Same historical unit conversion as the production bare-momentum route.
  hbar=6.582119569e-16;c=2.99792458e8;em=.510998950e6;metertoang=1.e10
  unit=hbar**4/em**2*metertoang**4*c**4
  expected=sqrt(3d0)/36d0*unit
  phase=cmplx(-.5,sqrt(3.)/2.)
  do spinor=1,2
    np=3*spinor;nplist=np
    coeffn=0;coeffm=0
    coeffn(1:3)=1./sqrt(real(np))
    coeffm(1:3)=[cmplx(1.,0.),phase,conjg(phase)]/sqrt(real(np))
    if(spinor==2)then
      coeffn(4:6)=coeffn(1:3);coeffm(4:6)=coeffm(1:3)
    endif
    open(10,status='scratch',access='direct',form='unformatted',recl=64)
    write(10,rec=4)(coeffn(i),i=1,np)
    write(10,rec=5)(coeffm(i),i=1,np)
    call kubo_berry_curvature(lower,b1,b2,b3,wk,1,nband,10d0,spinor,nplist,nbmax,6,nk,1,1,0,0)
    call kubo_berry_curvature(upper,b1,b2,b3,wk,1,nband,10d0,spinor,nplist,nbmax,6,nk,2,1,0,0)
    if(abs(lower(1)-expected)>2d-6*abs(expected))stop 1
    if(abs(upper(1)+expected)>2d-6*abs(expected))stop 2
    if(abs(lower(1)+upper(1))>1d-12)stop 3
    if(spinor==2)then
      call write_kubo_band_csv('fixture.csv',1,1,1,nk,nband,wk,lower)
      call write_kubo_band_csv('fixture.csv',1,2,1,nk,nband,wk,upper)
    endif
    close(10)
  enddo
  print *, 'KUBO_LOOP_PASS'
end program

subroutine ener_read(ener,isp,ik,nk,nband)
  implicit none
  integer isp,ik,nk,nband
  real(8) ener(nband)
  ener=[0d0,2d0]
end subroutine

subroutine plindx(ig,ncnt,ispinor,wk,b1,b2,b3,nbmax,np,ecut,npmax)
  implicit none
  integer ncnt,ispinor,nbmax(3),np,npmax,ig(3,npmax)
  real(8) wk(3),b1(3),b2(3),b3(3),ecut
  ig=0;ncnt=3;ig(1,2)=1;ig(2,3)=1
end subroutine

subroutine vaspberry_fail
  stop 4
end subroutine
