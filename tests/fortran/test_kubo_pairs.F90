! Exercise production export with complex plane-wave states and exact degeneracy.
program test_kubo_pairs
  use,intrinsic :: ieee_arithmetic,only: ieee_value,ieee_quiet_nan
  implicit none
#ifdef MPI_USE
  include 'mpif.h'
#endif
  integer,parameter :: nk=3,nband=3,npmax=8,npairs=3
  integer :: spinor,np,nbmax(3),nplist(nk),i,ik,isp,ispin,nprocs,myrank,comm,ierr
  real(8) :: a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),wk(3,nk)
  real(8) :: values(3,npairs),fallback(3,npairs),ener(nband),e2(nband),omega(nk),gap(nk),recovered
  complex(4) :: coeff(npmax,nband)
  character(32) :: mode
  character(64) :: path
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
  a1=2d0*acos(-1d0)*b1;a2=2d0*acos(-1d0)*b2;a3=2d0*acos(-1d0)*b3
  if(mode=='overflow')then
    nplist=4
    call write_kubo_pairs_csv('overflow.csv',1,1,1,nk,50000,a1,a2,a3,b1,b2,b3,wk,10d0, &
         nplist,nbmax,npmax,nprocs,myrank,comm)
    stop 99
  endif
  do spinor=1,2
    np=4*spinor;nplist=np
    coeff=0
    coeff(1:4,1)=[(1.,0.),(1.,0.),(1.,0.),(1.,0.)]/sqrt(real(np))
    coeff(1:4,2)=[(1.,0.),(0.,1.),(-1.,0.),(0.,-1.)]/sqrt(real(np))
    coeff(1:4,3)=[(1.,0.),(-1.,0.),(1.,0.),(-1.,0.)]/sqrt(real(np))
    if(spinor==2)coeff(5:8,:)=coeff(1:4,:)
    if(mode=='nonfinite-coeff')coeff(1,2)=cmplx(ieee_value(0.,ieee_quiet_nan),0.)
    ispin=3-spinor
    open(10,status='scratch',access='direct',form='unformatted',recl=128)
    do isp=1,ispin
      do ik=1,nk
        do i=1,nband
          write(10,rec=3+(ik-1)*(nband+1)+nk*(nband+1)*(isp-1)+i)coeff(1:np,i)
        enddo
      enddo
    enddo
    path='pairs-scalar.csv'
    if(spinor==2)path='pairs-spinor.csv'
    ! Export first so failures exercise absence of the completion footer.
    do isp=1,ispin
      call write_kubo_pairs_csv(path,isp,ispin,spinor,nk,nband,a1,a2,a3,b1,b2,b3,wk,10d0, &
           nplist,nbmax,npmax,nprocs,myrank,comm)
    enddo
    call kubo_bundle_curvature(omega,gap,b1,b2,b3,wk,1,nband,10d0,spinor, &
         nplist,nbmax,npmax,nk,1,2,nprocs,myrank,comm)
    do ik=1,nk
      call kubo_pair_numerators(values,ener,nband,npairs,b1,b2,b3,wk(:,ik),1,nk,ik,10d0, &
           spinor,np,nbmax,npmax,67108864_8)
      call kubo_pair_numerators(fallback,e2,nband,npairs,b1,b2,b3,wk(:,ik),1,nk,ik,10d0, &
           spinor,np,nbmax,npmax,0_8)
      if(any(values/=fallback).or.any(ener/=e2))stop 10
      ! Pairs 1:3,2:3 recover the isolated occupied-bundle xy trace.
      recovered=values(3,2)/(ener(1)-ener(3))**2+values(3,3)/(ener(2)-ener(3))**2
      if(abs(recovered-omega(ik))>1d-12)stop 11
    enddo
    close(10)
  enddo
  if(myrank==0)print *, 'KUBO_PAIRS_PASS'
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
  ener=[0d0,0d0,ik+1d0]+(isp-1)*4d0
  if(mode=='normal')ener=[0d0,1d0,ik+2d0]+(isp-1)*4d0
  if(mode=='nonfinite-energy')ener(3)=ieee_value(ener(3),ieee_quiet_nan)
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
