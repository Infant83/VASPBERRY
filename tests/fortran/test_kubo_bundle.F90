! Orthogonal complex plane-wave states. Tests actual production WAVECAR loop.
program test_kubo_bundle
  implicit none
#ifdef MPI_USE
  include 'mpif.h'
#endif
  integer, parameter :: nk=3, nband=3, npmax=8
  integer :: spinor,np,nbmax(3),nplist(nk),i,ik,isp,nprocs,myrank,comm,ierr
  real(8) :: b1(3),b2(3),b3(3),wk(3,nk),omega(nk),gap(nk)
  real(8) :: lower(nk),middle(nk),upper(nk),total(nk),unit,expected
  real(8) :: hbar,c,em,metertoang
  complex(4) :: coeff(npmax,nband),saved(npmax,2)
  character(32) :: mode
  common /fixture_mode/ mode
  call get_command_argument(1,mode)
  nprocs=1;myrank=0;comm=0
#ifdef MPI_USE
  call MPI_INIT(ierr)
  comm=MPI_COMM_WORLD
  call MPI_COMM_SIZE(comm,nprocs,ierr)
  call MPI_COMM_RANK(comm,myrank,ierr)
#endif
  nbmax=1;b1=[1d0,0d0,0d0];b2=[0d0,1d0,0d0];b3=[0d0,0d0,1d0];wk=0d0
  hbar=6.582119569e-16;c=2.99792458e8;em=.510998950e6;metertoang=1.e10
  unit=hbar**4/em**2*metertoang**4*c**4
  call check_kubo_bundle_gaps(2,nk,nband,1,2)
  do spinor=1,2
    np=4*spinor;nplist=np
    coeff=0
    coeff(1:4,1)=[(1.,0.),(1.,0.),(1.,0.),(1.,0.)]/sqrt(real(np))
    coeff(1:4,2)=[(1.,0.),(0.,1.),(-1.,0.),(0.,-1.)]/sqrt(real(np))
    coeff(1:4,3)=[(1.,0.),(-1.,0.),(1.,0.),(-1.,0.)]/sqrt(real(np))
    if(spinor==2)coeff(5:8,:)=coeff(1:4,:)
    if(mode=='rotated')then
      saved=coeff(:,1:2)
      coeff(:,1)=(saved(:,1)+saved(:,2))/sqrt(2.)
      coeff(:,2)=(saved(:,1)-saved(:,2))/sqrt(2.)
    endif
    open(10,status='scratch',access='direct',form='unformatted',recl=128)
    do isp=1,2
      do ik=1,nk
        do i=1,nband
          write(10,rec=3+(ik-1)*(nband+1)+nk*(nband+1)*(isp-1)+i)coeff(1:np,i)
        enddo
      enddo
    enddo
    call kubo_bundle_curvature(omega,gap,b1,b2,b3,wk,1,nband,10d0,spinor, &
         nplist,nbmax,npmax,nk,1,2,nprocs,myrank,comm)
    do ik=1,nk
      expected=unit/(8d0*(ik+1d0)**2)
      if(abs(omega(ik)-expected)>3d-6*abs(expected))stop 10
      if(abs(gap(ik)-(ik+1d0))>1d-14)stop 11
    enddo
    if(mode=='normal')then
      call kubo_berry_curvature(lower,b1,b2,b3,wk,1,nband,10d0,spinor,nplist,nbmax,npmax,nk,1,nprocs,myrank,comm)
      call kubo_berry_curvature(middle,b1,b2,b3,wk,1,nband,10d0,spinor,nplist,nbmax,npmax,nk,2,nprocs,myrank,comm)
      call kubo_berry_curvature(upper,b1,b2,b3,wk,1,nband,10d0,spinor,nplist,nbmax,npmax,nk,3,nprocs,myrank,comm)
      if(maxval(abs(omega-lower-middle))>1d-12)stop 12
      if(maxval(abs(lower+middle+upper))>1d-12)stop 13
      call kubo_bundle_curvature(total,gap,b1,b2,b3,wk,1,nband,10d0,spinor, &
           nplist,nbmax,npmax,nk,2,2,nprocs,myrank,comm)
      if(maxval(abs(total-middle))>1d-12)stop 14
    endif
    ! Restore bundle gaps after the singleton check.
    call kubo_bundle_curvature(omega,gap,b1,b2,b3,wk,1,nband,10d0,spinor, &
         nplist,nbmax,npmax,nk,1,2,nprocs,myrank,comm)
    if(spinor==2.and.myrank==0)then
      call write_kubo_bundle_csv('bundle.csv',1,nk,nband,1,2,wk,omega,gap)
      call write_kubo_bundle_csv('bundle.csv',2,nk,nband,1,2,wk,omega,gap)
    endif
    call kubo_bundle_curvature(total,gap,b1,b2,b3,wk,1,nband,10d0,spinor, &
         nplist,nbmax,npmax,nk,1,3,nprocs,myrank,comm)
    if(any(total/=0d0))stop 15
    if(any(gap/=huge(1d0)))stop 16
    if(spinor==2.and.myrank==0)call write_kubo_bundle_csv('all.csv',1,nk,nband,1,3,wk,total,gap)
    close(10)
  enddo
  if(myrank==0)print *, 'KUBO_BUNDLE_PASS'
#ifdef MPI_USE
  call MPI_FINALIZE(ierr)
#endif
end program

subroutine ener_read(ener,isp,ik,nk,nband)
  use,intrinsic :: ieee_arithmetic,only: ieee_value,ieee_quiet_nan
  implicit none
  integer isp,ik,nk,nband
  real(8) ener(nband)
  character(32) mode
  common /fixture_mode/ mode
  ener=[0d0,0d0,ik+1d0]
  if(mode=='normal')ener=[0d0,1d0,ik+2d0]
  if(mode=='external-gap')ener(3)=0d0
  if(mode=='near-gap')ener(3)=5d-6
  if(mode=='later-spin'.and.isp==2)ener(3)=0d0
  if(mode=='nonfinite')ener(3)=ieee_value(ener(3),ieee_quiet_nan)
end subroutine

subroutine plindx(ig,ncnt,ispinor,wk,b1,b2,b3,nbmax,np,ecut,npmax)
  implicit none
  integer ncnt,ispinor,nbmax(3),np,npmax,ig(3,npmax)
  real(8) wk(3),b1(3),b2(3),b3(3),ecut
  ig=0;ncnt=4;ig(1,2)=1;ig(2,3)=1;ig(3,4)=1
end subroutine

subroutine vaspberry_fail
  implicit none
#ifdef MPI_USE
  include 'mpif.h'
  integer ierr
  call MPI_ABORT(MPI_COMM_WORLD,17,ierr)
#endif
  stop 17
end subroutine
