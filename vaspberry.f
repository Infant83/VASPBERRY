! PROGRAM VASPBERRY Version 1.0 (f77) for VASP
! Written by Hyun-Jung Kim (angpangmokjang@hanmail.net, Infant@kias.re.kr) 
!  Korea Institute for Advanced Study (KIAS)
!  Dep. of Phys., Hanyang Univ.
! Copyright 2015. Hyun-Jung Kim All rights reserved.
! Evaluate berry curvature omega(k) for a closed loop on a small patches in k-space
! version 0.1 rough version. not working at all           : 2015. Mar. 17. H.-J. Kim
! version 0.2 error fix for k-loop finding                : 2015. Mar. 18. H.-J. Kim
! version 0.3 ISPIN=1, ISPIN=2 available                  : 2015. Mar. 23. H.-J. Kim
! version 0.4 sign error fix for BZ boundary              : 2015. Mar. 25. H.-J. Kim
! version 0.5 bug fix for defining overlap matrix S(k,k') : 2015. Mar. 31. H.-J. Kim
! version 0.6 including the option for circular dichroism : 2015. Apr. 10. H.-J. Kim
! version 0.7 set coefficients type to be complex*16 and  : 2015. Apr. 12. H.-J. Kim
!             implemented simple routine for grid extending  
! version 0.8 routine for velocity expectation value      : 2015. Apr. 30. H.-J. Kim
! version 0.9 Z2 invariant via Fukui's method:not finished: 2016. Feb. 18. H.-J. Kim
! version 0.91 routine for wavefunction plot: -wf n -k k  : 2016. May. 10. H.-J. Kim
! version 0.92 MPI parallelization implemented (MPI_USE)  : 2016. Jun. 27. H.-J. Kim & Y.-K. Kang & 
!              (during the CAC workshop & KIAS)                            S.-B. Cho & S.-W. Kim  &
!              -for Berrycurvature & Chern number evaluation               Y.-J. Choi & S.-H. Lee
! version 0.93 MPI parallelization implemented (MPI_USE)  : 2016. Jun. 29. H.-J. Kim 
!              -for Z2 invariant evaluation (routines are modified..)
! routine change: routine for getting a determiant -> get_det : 2018, Jun. 28. H.-J. Kim

! last update and bug fixes : 2020. Jun. 10. by H.-J. Kim 

!#define MPI_USE
!#undef  MPI_USE

      PROGRAM VASPBERRY

      implicit real*8 (a-h, o-z)
      complex*8, allocatable :: coeff(:)
      complex*16, allocatable :: Siju(:,:),Sijd(:,:),Sijt(:,:)
      complex*16, allocatable :: coeff1u(:),coeff1d(:)
      complex*16, allocatable :: coeff2u(:),coeff2d(:)
      complex*16,allocatable :: cener(:)
      real*8,    allocatable :: berrycurv(:),recivec(:,:),wklp(:,:,:)
      real*8,    allocatable :: berrycurv_tot(:)
      real*8,    allocatable :: rnfield(:),rnfield_tot(:)
      real*8,    allocatable :: recivec_tot(:,:)
      real*8,    allocatable :: recilat_tot(:,:)
      real*8,    allocatable :: xrecivec(:,:),xberrycurv(:),wklist(:,:)
      real*8,    allocatable :: wnklist(:,:),xrnfield(:)
      real*8,    allocatable :: recilat(:,:),xrecilat(:,:),occ(:)
      real*8,    allocatable :: selectivity(:),xselectivity(:)
      real*16,   allocatable :: ener(:)
      integer,   allocatable :: ig(:,:),nplist(:)
      dimension selectivitymax(4),selectivitymin(4)
      dimension a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),a2xa3(3),sumkg(3)
      dimension wk(3),wkk(3,5),ikk(5),isgg(2,5),npl(5),itr(5),itrim(5)
      dimension nbmax(3),xb(3),berrymax(4),berrymin(4),ng(3),rs(3)
      complex*16 csum1,csum2
      complex*16  detS(4),detA,detLOOP
      integer k, n, nkx, nky,nini,nmax,ns,ne,icd,ivel
      character*75 filename,foname,fonameo,fbz,ver_tag,vdirec
      data c/0.262465831d0/ ! constant c = 2m/hbar**2 [1/eV Ang^2]
      real*8  rfield,rnnfield,rnnfield_bottom
      real*8  rnnfield_tot,rnnfield_bottom_tot
      real*8, allocatable:: w_half_klist(:,:)
      integer, allocatable:: i_half_klist(:),i_trim_klist(:)
      integer :: myrank, nprocs, ierr
#ifdef MPI_USE
      include 'mpif.h'
      INTEGER         status(MPI_STATUS_SIZE)
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      if(myrank == 0)then
       write(6,*)"THIS IS ROOT:",myrank
      endif
      if(myrank == 0)then
       time_1=MPI_WTIME()
      endif
#else
      nprocs=1
      myrank=0
#endif

      ver_tag="# VASPBERRY (Ver 1.0), by Hyun-Jung Kim. 2018. Aug. 23."
      pi=4.*atan(1.)

      !default settings
      kperiod=2
      nkx=12 ;nky=12 
      ispinor=2  ! 2 for soc, 1 for non-soc
      itr=0;itrim=0

!!$*  reading general informations
      call parse(filename,foname,nkx,nky,ispinor,icd,ixt,fbz,
     &   ivel,iz,ihf,nini,nmax,kperiod,it,iskp,ine,ver_tag,
     &   iwf,ikwf,ng,rs,imag)
      if(myrank == 0)call creditinfo(ver_tag)

      if (ixt .ne. 0) call extendingBZ(fbz,ixt)
      if (it .eq. 1) call test
      if(myrank == 0)write(6,'(A,A)')"# File reading... : ",filename
      call inforead(irecl,ispin,nk,nband,ecut,a1,a2,a3,filename)
      if (ispin .eq. 2) ispinor=1
      allocate(cener(nband),ener(nband),occ(nband))
      if(myrank == 0)write(6,'(A,I9)')"# TOTAL RECORD LENGTH = ",irecl
      allocate(wklist(3,nk),nplist(nk))
      if(iz == 0)then
       iz2=1 ! enhancing factor for variable size define in subroutines
       allocate(berrycurv(nk))
       allocate(berrycurv_tot(nk))
       allocate(xberrycurv(kperiod*2*kperiod*2*nk))
      elseif(iz == 1)then
       iz2=4 ! enhancing factor for variable size define in subroutines
       allocate(rnfield(nk*iz2),rnfield_tot(nk*iz2))
       allocate(xrnfield(kperiod*2*kperiod*2*nk*iz2))
       allocate(wnklist(3,nk*iz2))
       allocate(w_half_klist(3,nk),i_half_klist(nk),i_trim_klist(nk))
      endif
      do ik=1,nk 
       ne_temp0=ne
       ne=0
       irec=3+(ik-1)*(nband+1)
       read(10,rec=irec) xnplane, (wk(i),i=1,3),
     &                  (cener(nn),occ(nn),nn=1,nband)
       wklist(:,ik)=wk(:)
       nplist(ik)=nint(xnplane)
       do n=1,nband; ne=ne+nint(occ(n)); enddo
       if(ik .gt. 1 .and. ne_temp0 .ne. ne .and. ine .eq. 0) then
         write(6,*)"error. !!! ne(K) /= ne(K') !!!",ne,ne_temp0,ik ;stop
       endif
      enddo
      if(ine .ne. 0) ne=ine ! manually specified ne ; useful for the semimetal
      ! check whether multi or single band calculation is performed
      if((nini.eq.nmax))then
       nini=nmax
       else if(nmax .eq. 999999)then
        if(ispin .eq. 2 .and. icd .eq. 0) nmax=ne
        if(ispin .eq. 1 .and. ispinor .eq. 2 .and. icd .eq. 0) nmax=ne
        if(ispin .eq. 1 .and. ispinor .eq. 1 .and. icd .eq. 0) nmax=ne/2
        if(ispin .eq. 2 .and. icd .eq. 1)then;nini=ne;nmax=nini+1;endif
        if(ispin .eq. 1 .and. ispinor .eq. 1 .and. icd .eq. 1)then
         nini=ne/2
         nmax=nini+1
        endif
        if(ispinor .eq. 2 .and. icd .eq. 1)then
         nini=ne
         nmax=nini+1
        endif
      endif ! check multi or single ?
      ns=nmax-nini+1
      if(myrank == 0)then
       write(6,'(A,I6)')   "# NELECT     : ",ne*ispin
       if (ispinor .eq. 2)then
        write(6,'(A,I6,A)')"# ISPIN      : ",ispin," (LSORBIT =.TRUE.)"
       else
        write(6,'(A,I6,A)')"# ISPIN      : ",ispin," (LSORBIT =.FALSE.)"
       endif
       write(6,'(A,F11.4)')  "# ENCUT (eV) : ",ecut
       write(6,'(A,I6)')     "# NKPOINT    : ",nk
       write(6,'(A,I6,A,I4)')"#  K-GRID    : ",nkx,"   X",nky
       write(6,'(A,I6)')     "# NBANDS     : ",nband
       write(6,'(A,3F13.6)') "# LATTVEC A1 : ",(a1(i),i=1,3)
       write(6,'(A,3F13.6)') "# LATTVEC A2 : ",(a2(i),i=1,3)
       write(6,'(A,3F13.6)') "# LATTVEC A3 : ",(a3(i),i=1,3)
      endif
      call recilatt(b1,b2,b3,dSkxky, a1,a2,a3,nkx,nky)
      if(myrank == 0)then
       write(6,'(A,3F13.6)') "# RECIVEC B1 : ",(b1(i),i=1,3)
       write(6,'(A,3F13.6)') "# RECIVEC B2 : ",(b2(i),i=1,3)
       write(6,'(A,3F13.6)') "# RECIVEC B3 : ",(b3(i),i=1,3)
       if(icd.eq.0)write(6,'(A,F13.6)')  "#  dk^2 = |dk1xk2| = ",dSkxky
      endif
      call reciproperty(nbmax,npmax, b1,b2,b3,ecut,ispinor)
      if(myrank == 0)then
       write(6,'(A,I6)')     "# NPMAX      : ",npmax
       write(6,'(A,I6)')     "#  NB1MAX    : ",nbmax(1)
       write(6,'(A,I6)')     "#  NB2MAX    : ",nbmax(2)
       write(6,'(A,I6)')     "#  NB3MAX    : ",nbmax(3)
      endif
#ifdef MPI_USE
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
      allocate(ig(3,npmax))
      allocate(coeff(npmax))
      if(iz == 1)then
       allocate(wklp(3,5,nk*4))
       allocate(recivec(3,4*nk),recilat(3,4*nk))
       allocate(recivec_tot(3,4*nk),recilat_tot(3,4*nk))
       allocate(xrecivec(3,kperiod*2*kperiod*2*nk*4))
       allocate(xrecilat(3,kperiod*2*kperiod*2*nk*4))
      elseif(iz ==0)then
       allocate(wklp(3,5,nk))
       allocate(recivec(3,nk),recilat(3,nk))
       allocate(recivec_tot(3,nk),recilat_tot(3,nk))
       allocate(xrecivec(3,kperiod*2*kperiod*2*nk))
       allocate(xrecilat(3,kperiod*2*kperiod*2*nk))
      endif
      nbtot=(2*nbmax(1)+2)*(2*nbmax(2)+2)*(2*nbmax(3)+1)
      allocate(coeff1u(nbtot),coeff1d(nbtot))
      allocate(coeff2u(nbtot),coeff2d(nbtot))
      if(icd.eq.1) allocate(selectivity(nk))
      if(icd.eq.1) allocate(xselectivity(kperiod*2*kperiod*2*nk))
      if(ivel .eq. 1 .and. myrank==0) then
        ni=nini 
        nj=nini
       call vel_expectation(a1,a2,a3,b1,b2,b3,
     &        kperiod,nbmax,npmax,ecut,ispinor,ispin,nband,ne,nk,wklist,
     &        nplist,ni,nj,irecl,filename,foname)
      endif

      if(iwf .ge. 1 .and. myrank==0)then ! plotting wavefunction
       ilp=1
       ikk(ilp)=ikwf
       iband=iwf
       npl(ilp)=nplist(ikk(ilp))
       wk(:)=wklist(:,ikk(ilp))
       nx=ng(1)
       ny=ng(2)
       nz=ng(3)
       if(nx*ny*nz .eq. 0) then
        nx=nbmax(1)
        ny=nbmax(2)
        nz=nbmax(3)
       endif
       write(6,'(A)')" "
       write(6,'(A)')"# Evaluating wavefunction..."
       call wf_plot(a1,a2,a3,b1,b2,b3,wk,
     &     nbmax,npmax,ecut,ispinor,ispin,nband,nk,
     &     npl,nmax,irecl,filename,foname,nx,ny,nz,iband,ikk,1,rs,imag)

      endif !iwf end
      
      if(iz .eq. 1) then  !get Z2 invariant using Fukui's method
       del=1E-6 ! criterion for k-point find
       call set_BZ_fukui(nnk,wnklist,w_half_klist,i_half_klist,
     &                   i_trim_klist,nhf, ispin,wklist,del,nk,iz2)

       do isp=1,ispin !ispin start
      if(myrank == 0)then
       write(6,*)" "
       if(isp.eq.1.and.ispin.eq.2) write(6,'(A,I2)')"# SPIN : UP"
       if(isp.eq.2.and.ispin.eq.2) write(6,'(A,I2)')"# SPIN : DN"
       write(6,*)" "
      endif
#ifdef MPI_USE
       nk_rank=nnk/nprocs
       nk_rank_mod=mod(nnk,nprocs)

       if(myrank<nk_rank_mod)then
           ink=myrank*(nk_rank+1)+1
           ifk=(myrank+1)*(nk_rank+1)
       else
           ink=myrank*nk_rank+nk_rank_mod+1
           ifk=(myrank+1)*nk_rank+nk_rank_mod
       endif

       recivec=0.
       recilat=0.
       recivec_tot=0.
       recilat_tot=0.
       rnnfield=0.
       rnnfield_bottom=0.
       rnnfield_tot=0.
       rnnfield_bottom_tot=0.
       rnfield_tot=0.
       rnfield=0.
#else
       ink=1;ifk=nnk
#endif
       do ik=ink, ifk !ik - loop
        call get_nfield(rnfield(ik),wklp, isp,ik,nk,nband,nkx,nky,
     &                      wnklist,ihf,nplist,nnk,iz,ns,kperiod,
     &                      w_half_klist,i_half_klist,i_trim_klist,nhf,
     &                      nbmax,npmax,ecut,ispinor,ispin,ne,irecl,
     &                      nmax,nini,b1,b2,b3,iz2)
        if((wklp(2,1,ik)+wklp(2,3,ik))/2. .gt. -del .and.
     &    (wklp(2,1,ik)+wklp(2,3,ik))/2. .lt. 0.5+del ) then
         rnnfield=rnnfield+rnfield(ik)
        else
         rnnfield_bottom=rnnfield_bottom+rnfield(ik)
        endif

!       write(6,'(A)')"#__________________________
!    &________________________________ "
!       write(6,'(A)')"# Berry Curvature F (A^-2) :
!    & -Im[logPI_S(det(S(K_s,K_s+1)))]/dk^2"
!       write(6,'(A)')"# N-field strength       :
!    & Sum_s{Im[log(det(S(K_s,k_s+1)))]}/2pi - F/2pi "

        do j=1,3
         recivec(j,ik)=(wklp(1,1,ik)+wklp(1,3,ik))*b1(j)*0.5+
     &                 (wklp(2,1,ik)+wklp(2,3,ik))*b2(j)*0.5+
     &                 (wklp(3,1,ik)+wklp(3,3,ik))*b3(j)*0.5
         recilat(j,ik)=(wklp(j,1,ik)+wklp(j,3,ik))*0.5
        enddo

!       write(6,'(A)')"#      kx        ky        kz(A^-1)
!    &    n-field strength      ,IK"
!       write(6,'(A,3F11.6,A,F8.3,I4,A)')'#',(recivec(i,ik),i=1,3),
!    &                             "      ",rnfield(ik),ik,"th-K"

       enddo   ! ik    end

#ifdef MPI_USE
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(rnnfield,rnnfield_tot,1,MPI_REAL8,
     &                MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(rnnfield_bottom,rnnfield_bottom_tot,1,MPI_REAL8,
     &                MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(rnfield,rnfield_tot,iz2*nk,MPI_REAL8,
     &                MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(recivec,recivec_tot,3*iz2*nk,MPI_REAL8,
     &                MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(recilat,recilat_tot,3*iz2*nk,MPI_REAL8,
     &                MPI_SUM,0,MPI_COMM_WORLD,ierr)

      rnnfield=rnnfield_tot
      rnnfield_bottom=rnnfield_bottom_tot
      rnfield=rnfield_tot
      recivec=recivec_tot
      recilat=recilat_tot
#endif

        if(myrank == 0)then
         write(6,*)" "
         write(6,'(A)')"# DONE!"
         
         if(nini .eq. nmax) then
          write(6,'(A,I4)')"# Z2 invariant for the BAND : ",nini
         else
          write(6,'(A,I4,A,I4)')"# Z2 invariant for the BANDS : ",nini
     &                                                       ," - ",nmax
         endif
         write(6,'(A,I2)')"# Z2 Invariant =    ",
     &                    mod((nint(rnnfield)),2)
         write(6,'(A,I2)')"# Z2 Invariant(bottom) =    ",
     &                    mod((nint(rnnfield_bottom)),2)
        
         call get_ext_variable(xrecivec,xrecilat,xrnfield,kext,
     &                    recivec,recilat,rnfield,nnk,kperiod,nk,iz2,
     &                    b1,b2,b3) ! extending over exteded-BZ
         
         call get_sorted_xrvari(xrecivec,xrecilat,xrnfield, 
     &                          kext,kperiod,nk,iz2) ! sorting
         call write_result(isp,ispin,ispinor,fonameo,foname,filename,
     &                     irecl,ecut,nk,nkx,nky,nband,b1,b2,b3,kperiod,
     &                     dSkxky,nini,nmax,xrecivec,xrecilat,kext,
     &                     xrnfield,rnnfield,rnnfield_bottom,0,0,
     &                     iz2,icd,iz,ivel,nprocs)
        endif !myrank ==0
       enddo !ispin end

       goto 9999
      endif !if Z2 end

      do isp=1,ispin   ! ispin start
       if(icd+ivel+iz+iwf.eq.0)then
        chernnumber=0.
      if(myrank == 0)then
        write(6,*)" "
        if(isp.eq.1.and.ispin.eq.2) write(6,'(A,I2)')"# SPIN : UP"
        if(isp.eq.2.and.ispin.eq.2) write(6,'(A,I2)')"# SPIN : DN"
        write(6,*)" "
      endif

#ifdef MPI_USE
        nk_rank=nk/nprocs
        nk_rank_mod=mod(nk,nprocs)
        if(myrank<nk_rank_mod)then
            ink=myrank*(nk_rank+1)+1
            ifk=(myrank+1)*(nk_rank+1)
        else
            ink=myrank*nk_rank+nk_rank_mod+1
            ifk=(myrank+1)*nk_rank+nk_rank_mod
        endif
        if(myrank == 0)then
         time_2=MPI_Wtime()
        endif
        berrycurv=0.
        berrycurv_tot=0.
        recivec=0.
        recilat=0.
        recivec_tot=0.
        recilat_tot=0.
#else
        ink=1;ifk=nk
#endif
        do ik=ink, ifk     !ik loop.
         call klpfind(wkk, isp,ik,nk,nband,nkx,nky,0,0,iz,iz2) !k-loop(4pts) of ik-kpoint (K)
         write(6,*)" "
         write(6,490) ik,(wkk(i,1),i=1,3)
  490    format("#* Closed loop for KPOINT",I5," : (",3F9.5,"  )")
         call kindxfind(ikk,isgg,npl,itr,itrim, wklist,nplist,wkk,nk,
     &                  nband,
     &                  iz,w_half_klist,i_half_klist,0,0,iz2) !find K-loop for ik
         write(6,500)(ikk(i),i=1,5)
  500    format('# K1>K2>K3>K4>K5=K1 (k-index):',
     &          I5,'  >',I5,'  >',I5,'  >',I5,'  >',I5)
         do j=1,5
          write(6,510)j,(wkk(i,j),i=1,3),j,(isgg(i,j),i=1,2)
  510     format('# K',I1,' = (',3F9.5,'  )',
     &          ', G',I1,' = (n1,n2) = (',I1,',',I1,')')
         enddo
         wklp(:,:,ik) = wkk(:,:)  ! closed loop (C) for each K : RECIVEC

!!$* get overlap matrix S_ij(k,k+1) over the C and PI_s [det S(k_s,k_s+1)], s=1,4
         allocate(Siju(ns,ns),Sijd(ns,ns),Sijt(ns,ns))
         do ilp=1,4  ! berry curvature loop for ik
          Siju=(0.,0.)
          Sijd=(0.,0.)
          Sijt=(0.,0.)
          ! construct overlap matrix S_ij(ikk1,ikk2)
!         write(6,570)ilp,ilp+1,ilp,ilp+1
! 570     format('#  Constructing overlap matrix, S_ij(K',I1,',K',I1,
!    &          ') = <u_i,K',I1,'(r)|u_j,K',I1,'(r)>')
          do ni=nini, nmax   ! calculate upto valence band maximum
           ncnt=0;coeff1u=(0.,0.);coeff1d=(0.,0.);coeff=(0.,0.)
           read(10,rec=(3+(ikk(ilp)-1)*(nband+1)+
     &                  nk*(nband+1)*(isp-1)+ni))(coeff(i),i=1,npl(ilp))
           do ig3=0,2*nbmax(3);    ig3p=ig3
            if (ig3.gt.nbmax(3))   ig3p=ig3-2*nbmax(3)-1
            do ig2=0,2*nbmax(2);   ig2p=ig2
             if (ig2.gt.nbmax(2))  ig2p=ig2-2*nbmax(2)-1
             do ig1=0,2*nbmax(1);  ig1p=ig1
              if (ig1.gt.nbmax(1)) ig1p=ig1-2*nbmax(1)-1
               call get_etot(etot, ilp,ni,b1,b2,b3,isgg,ig1p,ig2p,ig3p,
     &                          ig1,ig2,ig3,wkk,itrim,itr)

              if (etot.lt.ecut) then; ncnt=ncnt+1
               call get_incnt(incnt, ilp,ni,ig1p,ig2p,ig3p,nbmax,
     &                               itr,itrim,isgg,wkk)
               if(ispinor .eq. 2)then 
                coeff1d(incnt)=coeff(ncnt+npl(ilp)/2)
               endif
               coeff1u(incnt)=coeff(ncnt)
              endif
             enddo  !loop for ig1                       
            enddo   !loop for ig2
           enddo    !loop for ig3
           if (ispinor*ncnt.ne.npl(ilp)) then
            write(0,*) '*ni error - computed NPL=',2*ncnt,
     &                 ' != input nplane=',npl(ilp);stop
           endif
        
           do nj=nini, nmax
            ncnt=0;coeff2u=(0.,0.);coeff2d=(0.,0.);coeff=(0.,0.)
            read(10,rec=(3+(ikk(ilp+1)-1)*(nband+1)+nk*(nband+1)*(isp-1)
     &                   +nj)) (coeff(i),i=1,npl(ilp+1))
            do ig3=0,2*nbmax(3); ig3p=ig3
             if (ig3.gt.nbmax(3)) ig3p=ig3-2*nbmax(3)-1
             do ig2=0,2*nbmax(2); ig2p=ig2
              if (ig2.gt.nbmax(2)) ig2p=ig2-2*nbmax(2)-1
              do ig1=0,2*nbmax(1); ig1p=ig1
               if (ig1.gt.nbmax(1)) ig1p=ig1-2*nbmax(1)-1
               call get_etot(etot,ilp+1,nj,b1,b2,b3,isgg,ig1p,ig2p,ig3p,
     &                          ig1,ig2,ig3,wkk,itrim,itr)
               if (etot.lt.ecut) then; ncnt=ncnt+1
               call get_incnt(incnt, ilp+1,nj,ig1p,ig2p,ig3p,nbmax,
     &                               itr,itrim,isgg,wkk)
                if(ispinor .eq. 2)then
                 coeff2d(incnt)=coeff(ncnt+npl(ilp+1)/2)
                endif
                coeff2u(incnt)=coeff(ncnt)
               endif
              enddo   !loop for ig1
             enddo    !loop for ig2
            enddo     !loop for ig3
            if (ispinor*ncnt.ne.npl(ilp+1)) then
             write(0,*) '*nj error - computed NPL=',2*ncnt,
     &                  ' != input nplane=',npl(ilp+1);stop
            endif
            if(ispinor .eq. 2)then
             Siju(ni-nmax+ns,nj-nmax+ns)=dot_product(coeff1u,coeff2u)
             Sijd(ni-nmax+ns,nj-nmax+ns)=dot_product(coeff1d,coeff2d)
             Sijt(ni-nmax+ns,nj-nmax+ns)=Siju(ni-nmax+ns,nj-nmax+ns)+
     &                                   Sijd(ni-nmax+ns,nj-nmax+ns)
             else if (ispinor .eq. 1)then
              Sijt(ni-nmax+ns,nj-nmax+ns)=dot_product(coeff1u,coeff2u)
            endif
           enddo  ! loop for nj
          enddo   ! loop for ni

          if(nini .eq. nmax)then   ! get determinant : det(S)
           detS(ilp)=Sijt(1,1)
           else
!           call getdetA(detA, Sijt,ns)
            call get_det(detA, Sijt,ns)
            detS(ilp)=detA;detA=(0.,0.)
          endif
         enddo    ! loop for ilp

         detLOOP=detS(1)*detS(2)*detS(3)*detS(4)
         write(6,'(A)')"# "
         write(6,600)detLOOP
  600    format('#===>PI_S[det(S(K_s,K_s+1))] =',F16.8,'  +',F16.8,' i')
        
         berrycurv(ik)=-1.*aimag(log(detLOOP))/dSkxky

         write(6,'(A)')"#__________________________
     &________________________________ "
         write(6,'(A)')"# Berry Curvature (A^-2) :
     & -Im[logPI_S(det(S(K_s,K_s+1)))]/dk^2"
         do j=1,3
          recivec(j,ik)=(wklp(1,1,ik)+wklp(1,3,ik))*b1(j)*0.5+
     &                  (wklp(2,1,ik)+wklp(2,3,ik))*b2(j)*0.5+
     &                  (wklp(3,1,ik)+wklp(3,3,ik))*b3(j)*0.5
          recilat(j,ik)=(wklp(j,1,ik)+wklp(j,3,ik))*0.5
         enddo
         write(6,'(A)')"#      kx        ky        kz(A^-1)
     &   Berry Curvature (A^-2)"
         write(6,'(A,3F11.6,A,F16.6)')'#',(recivec(i,ik),i=1,3),"     ",
     &                                    berrycurv(ik)

#ifdef MPI_USE
         chernnumber=chernnumber+berrycurv(ik)*dSkxky/(2.*pi)
#else 
         chernnumber=chernnumber+berrycurv(ik)*dSkxky/(2.*pi)
         if (ik .eq. 1)then
          berrymax(4)=berrycurv(ik)
          berrymax(1)=recilat(1,ik)
          berrymax(2)=recilat(2,ik)
          berrymax(3)=recilat(3,ik)
          berrymin(4)=berrycurv(ik)
          berrymin(1)=recilat(1,ik)
          berrymin(2)=recilat(2,ik)
          berrymin(3)=recilat(3,ik)
          else if (ik .ge. 2 .and. berrycurv(ik) .ge. berrymax(4))then
           berrymax(4)=berrycurv(ik)
           berrymax(1)=recilat(1,ik)
           berrymax(2)=recilat(2,ik)
           berrymax(3)=recilat(3,ik)
          else if (ik .ge.2 .and. berrycurv(ik) .le. berrymin(4)) then
           berrymin(4)=berrycurv(ik)
           berrymin(1)=recilat(1,ik)
           berrymin(2)=recilat(2,ik)
           berrymin(3)=recilat(3,ik)
         endif
#endif
        
         deallocate(Sijt)
         deallocate(Siju)
         deallocate(Sijd)
        enddo         ! ik loop over

#ifdef MPI_USE
        if(myrank == 0)then
         time_3=MPI_Wtime()
         write(6,*)" "
         write(6,'(A)')"# DONE!"
        endif
#else
        write(6,*)" "
        write(6,'(A)')"# DONE!"
#endif
        
!!$*  SUMMARIZTION and output the results..
      if(myrank == 0)then
        write(6,'(A)')"# Chern Number is sum of Berry Curvature on 1BZ"
        if(nini .eq. nmax) then
         write(6,'(A,I4)')"# Chern Number for the BAND : ",nini
         else
          write(6,'(A,I4,A,I4)')"# Chern Number for the BANDS : ",nini
     &                                                     ," - ",nmax
        endif
      endif

#ifdef MPI_USE
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(chernnumber,chernnumber_tot,1,MPI_REAL8,
     &                  MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(berrycurv,berrycurv_tot,nk,MPI_REAL8,
     &                  MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(recivec,recivec_tot,3*nk,MPI_REAL8,
     &                  MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(recilat,recilat_tot,3*nk,MPI_REAL8,
     &                  MPI_SUM,0,MPI_COMM_WORLD,ierr)
         chernnumber=chernnumber_tot
         berrycurv(:)=berrycurv_tot(:)
         recivec(:,:)=recivec_tot(:,:)
         recilat(:,:)=recilat_tot(:,:)
#endif
       if(myrank==0)then
        write(6,'(A,F16.6)')"# Chern Number =    ",chernnumber
       endif
!KKKKK NOTE : START-sorting

       if(myrank==0)then
       call get_ext_variable(xrecivec,xrecilat,xberrycurv,kext,
     &                  recivec,recilat,berrycurv,nk,kperiod,nk,iz2,
     &                  b1,b2,b3) ! extending over exteded-BZ

       call get_sorted_xrvari(xrecivec,xrecilat,xberrycurv,
     &                        kext,kperiod,nk,iz2) ! sorting
#ifdef MPI_USE
       call write_result(isp,ispin,ispinor,fonameo,foname,filename,
     &                   irecl,ecut,nk,nkx,nky,nband,b1,b2,b3,kperiod,
     &                   dSkxky,nini,nmax,xrecivec,xrecilat,kext,
     &                   xberrycurv,chernnumber,0.,0,0,
     &                   iz2,icd,iz,ivel,nprocs)
#else
       call write_result(isp,ispin,ispinor,fonameo,foname,filename,
     &                   irecl,ecut,nk,nkx,nky,nband,b1,b2,b3,kperiod,
     &                   dSkxky,nini,nmax,xrecivec,xrecilat,kext,
     &                   xberrycurv,chernnumber,0.,berrymax,berrymin,
     &                   iz2,icd,iz,ivel,nprocs)
#endif
       endif

!!! ######### LOOP for OPTICAL SELECTIVITY ###########################################
       elseif(icd.eq.1 .and. myrank==0)then   ! calculate optical selectivity
        call optical_selectivity(selectivity, 
     &                           selectivitymax,selectivitymin,
     &             b1,b2,b3,wklist,isp,
     &             nband,ecut,ispinor,nplist,nbmax,npmax,nk,nini,nmax)
        do ik=1,nk
         do j=1,3
          recivec(j,ik)=wklist(1,ik)*b1(j)+
     &                  wklist(2,ik)*b2(j)+
     &                  wklist(3,ik)*b3(j)
          recilat(j,ik)=wklist(j,ik)
         enddo
        enddo

       call get_ext_variable(xrecivec,xrecilat,xselectivity,kext,
     &                  recivec,recilat,selectivity,nk,kperiod,nk,iz2,
     &                  b1,b2,b3) ! extending over exteded-BZ

       call get_sorted_xrvari(xrecivec,xrecilat,xselectivity,
     &                        kext,kperiod,nk,iz2) ! sorting
       call write_result(isp,ispin,ispinor,fonameo,foname,filename,
     &                  irecl,ecut,nk,nkx,nky,nband,b1,b2,b3,kperiod,
     &                  dSkxky,nini,nmax,xrecivec,xrecilat,kext,
     &                  xselectivity,0,0.,selectivitymax,selectivitymin,
     &                  iz2,icd,iz,ivel,nprocs)

       endif !icd over
!!! ######### LOOP END for OPTICAL SELECTIVITY ###########################################

      enddo !ispin loop over

      if(iskp .eq. 1 .and. myrank==0)call write_special_kpoint(b1,b2,b3)
#ifdef MPI_USE
      if(myrank==0)then
#endif
      write(6,'(A)')"# DONE! "
      do isp=1, ispin
       if(isp .eq. 1 .and. ispinor .eq. 1 .and. ispin .eq. 2) then
        write(6,'(A,A,A)')"#  Results are summarized in ",TRIM(foname),
     &                    ".UP.dat for spin-1"
        else if (isp.eq.2 .and. ispinor.eq.1 .and. ispin.eq.2)then
         write(6,'(A,A,A)')"#  Results are summarized in ",
     &                     TRIM(foname),".DN.dat for spin-2"
        else if (isp .eq. 1 .and. ispinor .eq. 2) then
         write(6,'(A,A,A)')"#  Results are summarized in ",
     &                     TRIM(foname),".dat"
        else if (isp.eq.1 .and. ispinor.eq.1 .and. ispin.eq.1) then
         write(6,'(A,A,A)')"#  Results are summarized in ",
     &                     TRIM(foname),".dat"
       endif
      enddo
#ifdef MPI_USE
      endif
#endif
 9999 if(myrank==0) write(6,*)"end of program"

#ifdef MPI_USE
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

      deallocate(ig)
      deallocate(coeff)
      deallocate(coeff1u)
      deallocate(coeff1d)
      deallocate(coeff2u)
      deallocate(coeff2d)
      deallocate(recilat)
      deallocate(recivec)
      deallocate(xrecilat)
      deallocate(xrecivec)
      deallocate(wklp)
      deallocate(wklist)
      deallocate(cener)
      deallocate(ener)
      if(iz == 0)then
       deallocate(berrycurv)
       deallocate(berrycurv_tot)
       deallocate(xberrycurv)
      elseif(iz == 1)then
       deallocate(rnfield)
       deallocate(rnfield_tot)
       deallocate(wnklist)
       deallocate(xrnfield)
       deallocate(i_trim_klist)
       deallocate(i_half_klist)
       deallocate(w_half_klist)
      endif
      deallocate(occ)
      deallocate(nplist)
      if(icd.eq.1)deallocate(selectivity)
      if(icd.eq.1)deallocate(xselectivity)
#ifdef MPI_USE
      if(myrank == 0 .and. iz+ivel+icd+iwf .eq. 0)then
       time_4=MPI_Wtime()
       write(6,*) "Data reading      : ",time_2-time_1
       write(6,*) "Parallel sequence : ",time_3-time_1
       write(6,*) "End sequence      : ",time_4-time_3
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_FINALIZE(ierr)
#endif
      end program

!!$*  subroutine for writing results
      subroutine write_result(isp,ispin,ispinor,fonameo,foname,filename,
     &                  irecl,ecut,nk,nkx,nky,nband,b1,b2,b3,kperiod,
     &                  dSkxky,nini,nmax,xrecivec,xrecilat,kext,
     &                  xrvari,rvari,rvari2,rvari3,rvari4,
     &                  iz2,icd,iz,ivel,nprocs)
      implicit real*8 (a-h,o-z)
      dimension b1(3),b2(3),b3(3)
      real*8  xrecivec(3,kperiod*2*kperiod*2*nk*iz2)
      real*8  xrecilat(3,kperiod*2*kperiod*2*nk*iz2)
      real*8  xrvari(kperiod*2*kperiod*2*nk*iz2)
      real*8  rvari,rvari2
      real*8  xb(3),xtemp,rvari3(4),rvari4(4)
      character*75 filename,foname,fonameo

       if(isp .eq. 1 .and. ispinor .eq. 1 .and. ispin.eq.2) then
        write(fonameo,'(A,A)')TRIM(foname),'.UP.dat'
        else if (isp .eq. 2 .and. ispinor .eq. 1 .and. ispin.eq.2)then
         write(fonameo,'(A,A)')TRIM(foname),'.DN.dat'
        else if (isp .eq. 1 .and. ispinor .eq. 2) then
         write(fonameo,'(A,A)')TRIM(foname),'.dat'
        else if (isp .eq. 1 .and. ispinor .eq. 1 .and. ispin.eq.1)then
         write(fonameo,'(A,A)')TRIM(foname),'.dat'
       endif
       open(32,file=fonameo,status='unknown')
       write(32,'(A,I4,A)')"# Job running on ",nprocs," total cores"
       write(32,'(A,A)')   "# File reading...  : ",filename
       write(32,'(A,I9)')"# TOTAL RECORD LENGTH = ",irecl
       if (ispinor .eq. 2)then
        write(32,'(A,I6,A)')"# ISPIN            : ",ispin,
     &                      " (LSORBIT = .TRUE.)"
        else
         write(32,'(A,I6,A)')"# ISPIN            : ",ispin,
     &                       " (LSORBIT = .FALSE.)"
       endif
       write(32,'(A,F11.4)')  "# ENCUT (eV)       : ",ecut
       write(32,'(A,I6)')     "# NKPOINT          : ",nk
       write(32,'(A,I6,A,I4)')"#  K-GRID          : ",nkx,"   X",nky
       write(32,'(A,I6)')     "# NBANDS           : ",nband
       write(32,'(A,3F13.6)') "# RECIVEC B1 (A^-1): ",(b1(i),i=1,3)
       write(32,'(A,3F13.6)') "# RECIVEC B2       : ",(b2(i),i=1,3)
       write(32,'(A,3F13.6)') "# RECIVEC B3       : ",(b3(i),i=1,3)
       write(32,'(A,F13.6)')   "#  dk^2 = |dk1xk2| = ",dSkxky
       write(32,*)" "

       if(nini .eq. nmax) then
        if(iz == 1)then
         write(32,'(A,I4)')"# Z2 invariant for the BAND : ",nini
        elseif(iz+ivel+icd .eq. 0)then
         write(32,'(A,I4)')"# Chern Number for the BAND : ",nmax
        endif
       else
        if(iz == 1)then
         write(32,'(A,I4,A,I4)')"# Z2 invariant for the BANDS : ",nini,
     &                         " - ",nmax
        elseif(iz+ivel+icd .eq. 0)then
         write(32,'(A,I4,A,I4)')"# Chern Number for the BANDS : ",nini,
     &                         "    -  ",nmax
        endif
       endif
 
       if(iz == 1)then !Z2 INVARIANT
        write(32,'(A,I2)')"# Z2 Invariant (top) =    ",
     &                  mod((nint(rvari)),2)
        write(32,'(A,I2)')"# Z2 Invariant (bottom) =    ",
     &                  mod((nint(rvari2)),2)
        write(32,'(A)')"# Berry Curvature F (A^-2) :
     &  -Im[logPI_S(det(S(K_s,K_s+1)))]/dk^2"
        write(32,'(A)')"# N-field strength       :
     &  Sum_s{Im[log(det(S(K_s,k_s+1)))]}/2pi - F/2pi "
        write(32,'(A)')"# (cart) kx        ky        kz(A^-1)
     &    n-field strength      ,   (recip)kx        ky        kz"

       elseif(icd == 1)then !OPTICAL SELECTIVITY
        write(32,'(A,I4,A,I4)')"# OPTICAL SELECTIVITY BETWEEN BANDS: ",
     &                         nini,"    -  ",nmax
        write(32,'(A)')"# n(k,w_cv)= |P(k,s,cv,+)|^2 - |P(k,s,cv,-)|^2"
        write(32,'(A)')"#            ---------------------------------"
        write(32,'(A)')"#            |P(k,s,cv,+)|^2 + |P(k,s,cv,-)|^2"
        write(32,'(A)')"#  The TRANSITION MATRIX ELEMENT P ="
        write(32,'(A)')"#   P(k,s,cv,+ or -) = 1/sqrt(2)[p_x(k,cv,s) +
     &(or -) i*p_y(k,cv,s)]"
        write(32,'(A)')"#  THE INTERBAND TRANSITION MATRIX p_x,y ="
        write(32,'(A)')"#   p_x,y(k,cv,s)=<psi(k,c,s)|-i*hbar*1/dx(y)|
     &psi(k,v,s>"
        write(32,'(A,4F16.6)')"# MAXVAL of SELECTIVITY at kx,ky,kz 
     &(in reci)= ",(rvari3(i),i=1,4)
        write(32,'(A,4F16.6)')"# MINVAL of SELECTIVITY at kx,ky,kz 
     &(in reci)= ",(rvari4(i),i=1,4)
        write(32,'(A)')"# (cart) kx        ky        kz(A^-1)
     &   selectivity(n(k)),        (recip)kx        ky        kz"

       else !BERRYCURVATURE
        write(32,'(A)')"# Chern Number is sum of 
     &Berry Curvature over 1BZ"
        write(32,'(A,F16.6)')"# Chern Number =   ",rvari
        write(32,'(A,4F16.6)')"# MAXVAL of BERRYCURV at kx,ky,kz 
     &(in reci)= ",(rvari3(i),i=1,4)
        write(32,'(A,4F16.6)')"# MINVAL of BERRYCURV at kx,ky,kz 
     &(in reci)= ",(rvari4(i),i=1,4)
        write(32,'(A)')"# Berry Curvature (A^-2) :
     &-Im[logPI_S(det(S(K_s,K_s+1)))]/dk^2"
        write(32,'(A)')"# (cart) kx        ky        kz(A^-1)
     &   Berry Curvature (A^-2),   (recip)kx        ky        kz"
       endif

       do ik=1,kext
        write(32,'(3F11.6,A,F11.6,A,3F11.6)')(xrecivec(i,ik),i=1,3),
     &        "     ",xrvari(ik),"            ",
     &        (xrecilat(i,ik),i=1,3)
       enddo

       close(32)
       return
      end subroutine write_result

!!$*  subroutine for sorting
      subroutine get_sorted_xrvari(xrecivec,xrecilat,xrvari,
     &                        kext,kperiod,nk,iz2)
      implicit real*8 (a-h,o-z)
      real*8  xrecivec(3,kperiod*2*kperiod*2*nk*iz2)
      real*8  xrecilat(3,kperiod*2*kperiod*2*nk*iz2)
      real*8  xrvari(kperiod*2*kperiod*2*nk*iz2)
      real*8  xb(3),xtemp

      do k=kext-1,1,-1  ! sorting kx
       do j=1,k
        if(xrecivec(1,j+1) .gt. xrecivec(1,j))then
         xb(:)=xrecivec(:,j)
         xrecivec(:,j)=xrecivec(:,j+1)
         xrecivec(:,j+1)=xb(:)
         xb(:)=xrecilat(:,j)
         xrecilat(:,j)=xrecilat(:,j+1)
         xrecilat(:,j+1)=xb(:)
         xtemp=xrvari(j)
         xrvari(j)=xrvari(j+1)
         xrvari(j+1)=xtemp
        endif
       enddo
      enddo
      do k=kext-1,1,-1  ! sorting ky
       do j=1,k
        if(xrecivec(1,j+1) .eq. xrecivec(1,j))then
         if(xrecivec(2,j+1) .gt. xrecivec(2,j))then
         xb(:)=xrecivec(:,j)
         xrecivec(:,j)=xrecivec(:,j+1)
         xrecivec(:,j+1)=xb(:)
         xb(:)=xrecilat(:,j)
         xrecilat(:,j)=xrecilat(:,j+1)
         xrecilat(:,j+1)=xb(:)
         xtemp=xrvari(j)
         xrvari(j)=xrvari(j+1)
         xrvari(j+1)=xtemp
         endif
        endif
       enddo
      enddo
      return
      end subroutine get_sorted_xrvari

!!$*  subroutine for extending data_set distribution over extended BZ: 2D array
      subroutine get_ext_variable(xrecivec,xrecilat,xrvari,kext, 
     &                  recivec,recilat,rvari,nnk,kperiod,nk,iz2,
     &                  b1,b2,b3) !iz2=4 for z2, other 1
      implicit real*8 (a-h,o-z)
      real*8  recivec(3,nk*iz2),recilat(3,nk*iz2)
      real*8  xrecivec(3,kperiod*2*kperiod*2*nk*iz2)
      real*8  xrecilat(3,kperiod*2*kperiod*2*nk*iz2)
      real*8  rvari(nk*iz2),xrvari(kperiod*2*kperiod*2*nk*iz2)
      dimension b1(3),b2(3),b3(3)

       kk=0   ! extend variable distribution over extended BZ 
       do ib2=-1*(kperiod-1)+1,kperiod
        do ib1=-1*(kperiod-1)+1,kperiod
         do ik=1,nnk
          kk=kk+1
          xrecivec(1,kk)=recivec(1,ik) +(ib1-1)*b1(1)+(ib2-1)*b2(1)
          xrecivec(2,kk)=recivec(2,ik) +(ib1-1)*b1(2)+(ib2-1)*b2(2)
          xrecivec(3,kk)=recivec(3,ik) +(ib1-1)*b1(3)+(ib2-1)*b2(3)
          xrecilat(1,kk)=recilat(1,ik) +(ib1-1)
          xrecilat(2,kk)=recilat(2,ik) +(ib2-1)
          xrecilat(3,kk)=recilat(3,ik)
          xrvari(kk)=rvari(ik)
          kext=kk
         enddo
        enddo
       enddo

      return
      end subroutine 

!!$*  subroutine for getting planewave coefficient for given k & n
      subroutine wf_plot(a1,a2,a3,b1,b2,b3,wk,
     &    nbmax,npmax,ecut,ispinor,ispin,nband,nk,
     &    npl,nmax,irecl,filename,foname,nx,ny,nz,iband,ikk,ilp,rs,imag)
 
      implicit real*8 (a-h,o-z)
      complex*8 coeff(npmax),coeffi(npmax)
      complex*16 csum(ispinor,nx*ny*nz)
      real*8 wklist(3,nk),w_half_klist(3,nk),pi,pi2
      real*8  recivec(3),recilat(3,nk*4)
      dimension a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),a2xa3(3),sumkg(3)
      dimension wk(3),wkg(3),rr(3),rs(3)
      dimension wkgr(npmax)
      integer nplist(nk)
      integer ikk(5),isgg(2,5),npl(5),itrim(5),ig(3,npmax)
      integer nbmax(3)
      integer nk,nband,npmax,kperiod,ispin,irecl
      character*75 filename,foname,fonameo,fonameoi
      data c/0.262465831d0/ ! constant c = 2m/hbar**2 [1/eV Ang^2]
      pi=4.*atan(1.)
      pi2=pi*2.

      itrim=0     
      call vcross(a2xa3,a2,a3)
      vol=dot_product(a1,a2xa3)
      ngrid=nx*ny*nz
      nline=int(ngrid/5.)
      nresi=mod(ngrid,5)

       do isp=1,ispin
        write(6,'(A,I1)')" ##Processing for spin-",isp
        csum=(0.,0.);wkgr=0.
        write(fonameo,'(A,A,I1)')TRIM(foname),"-SPIN",isp
        call get_coeff(coeff,isgg, 1, iband,ikk,npl,nband,nk,isp,
     &                              npmax,itrim)
        call plindx(ig,ncnt,ispinor,wk,b1,b2,b3,nbmax,npl(1),ecut,npmax)

        open(unit=18,file=fonameo,status='unknown')
        call write_CHGCAR_head(18,fonameo,ikk,isp,iband,a1,a2,a3,
     &                         nx,ny,nz,rs)

        if(imag .eq. 1)then
         write(fonameoi,'(A,A,I1)')TRIM(foname),"-IM-SPIN",isp
         open(unit=19,file=fonameoi,status='unknown')
        call write_CHGCAR_head(19,fonameoi,ikk,isp,iband,a1,a2,a3,
     &                         nx,ny,nz,rs)
        endif

        ii=0
        do i3=1,nz
         do i2=1,ny
          do i1=1,nx
           ii=ii+1
           do iispinor=1,ispinor
            wkgr(:)=(wk(1)+ig(1,1:ncnt))*((i1-1)/dble(nx)+rs(1)) +
     &              (wk(2)+ig(2,1:ncnt))*((i2-1)/dble(ny)+rs(2)) +
     &              (wk(3)+ig(3,1:ncnt))*((i3-1)/dble(nz)+rs(3))
           
            csum(iispinor,ii)= sum(
     &       coeff(1+npl(1)/2*(iispinor-1):ncnt+npl(1)/2*(iispinor-1))*
     &       cdexp(pi2*cmplx(0.,1.)*wkgr(1:ncnt))*dsqrt(vol) )
           enddo !ispinor
          enddo !i1
         enddo !i2
         if(mod(i3,nint(nz/10.)) == 0 )then
          write(6,'(F5.1,A)',advance='yes')i3/dble(nz)*100,"%"
         endif
        enddo !i3

        do iline=1,nline
         write(18,'(5E15.7)')(real(csum(1,(iline-1)*5+i)),i=1,5)
        enddo
        if(nresi .ge. 1)then
         write(18,'(5E15.7)')(real(csum(1,(nline)*5+i)),i=1,nresi)
        endif

        if(ispinor .ge. 2)then
         write(18,*)nx,ny,nz
         do iline=1,nline
          write(18,'(5E15.7)')(real(csum(2,(iline-1)*5+i)),i=1,5)
         enddo
         if(nresi .ge. 1)then
          write(18,'(5E15.7)')(real(csum(2,(nline)*5+i)),i=1,nresi)
         endif
        endif

        if(imag .eq. 1)then
         do iline=1,nline
          write(19,'(5E15.7)')(dimag(csum(1,(iline-1)*5+i)),i=1,5)
         enddo
         if(nresi .ge. 1)then
          write(19,'(5E15.7)')(dimag(csum(1,(nline)*5+i)),i=1,nresi)
         endif

         if(ispinor .ge. 2)then
          write(19,*)nx,ny,nz
          do iline=1,nline
           write(19,'(5E15.7)')(dimag(csum(2,(iline-1)*5+i)),i=1,5)
          enddo
          if(nresi .ge. 1)then
           write(19,'(5E15.7)')(dimag(csum(2,(nline)*5+i)),i=1,nresi)
          endif
         endif
        endif

       enddo !ispin
      
      close(18)
      if(imag .eq. 1) close(19)
      stop
      end subroutine wf_plot

!!$*  subroutine for CHGCAR header (lattice info + atomic coordinates) writing
      subroutine write_CHGCAR_head(ID,fonameo,ikk,isp,iband,a1,a2,a3,
     &                             nx,ny,nz,rs)
        implicit real*8 (a-h,o-z)
        dimension a1(3),a2(3),a3(3),ikk(1),coord(3)
        dimension n_atom(10),rs(3),rss(3)
        character*4 at_name(10),const(3)
        character*75 fonameo
        character dummy

        !get total number of atoms : read EIGENVAL header
        open(unit=ID+10,file='EIGENVAL',status='old',iostat=IERR)
        if (IERR.ne.0) write(6,*) 'open error EIGENVAL - iostat =',IERR
        read(ID+10,*)natom,idummy,idummy,idummy
        close(ID+10)

        !get atomic info: read POSCAR
        open(unit=ID+10,file='POSCAR',status='old',iostat=IERR)
        if (IERR.ne.0) write(6,*) 'open error POSCAR - iostat =',IERR
        do i=1,6 ;read(ID+10,*) dummy ;enddo
         
        itot=0
        itype=1
        iselect=0 
        icont=0
        n_atom=0
        idirect=0

        do while (ifinish .ne. 1)
         read(ID+10,*) n_atom(1:itype)
         itot = sum(n_atom)
         if (itot .ne. natom ) then
          ifinish = 0
          itype = itype + 1
         elseif (itot .eq. natom ) then
          ifinish = 1
         elseif (itot .gt. natom) then
          write(6,'(A,I4)')"ERROR: total number of atom exeeds ",natom
         endif
         if(ifinish .ne. 1) backspace (unit=ID+10)
        enddo
        backspace (unit=ID+10)
        backspace (unit=ID+10)
        read(ID+10,*) at_name(1:itype)
        read(ID+10,*) dummy

        write(ID,'(A,I3,A,I3,A,I1)')"WAVEFUNCTION: BAND= ",iband,
     &                         " ,KP= ",ikk(1)," SPIN= ",isp
        write(ID,*)1.0000
        write(ID,'(3F20.16)')a1(:)
        write(ID,'(3F20.16)')a2(:)
        write(ID,'(3F20.16)')a3(:)
        write(ID,*)at_name(1:itype)
        write(ID,*)n_atom(1:itype)

        do while(icont .eq. 0) ! check selective or dynamics
         read(ID+10,*) dummy
         if (dummy == "S" .or. dummy == "s") then 
          write(ID,'(A)') "Selective dynamics"
          icont=0;iselect=1
         elseif(dummy == "D" .or. dummy == "d")then
          write(ID,'(A)') "Direct"
          icont=1;idirect=1
          rss(:)=rs(:)
         elseif(dummy == "C" .or. dummy == "c" .or. dummy == "k" .or. 
     &          dummy == "K")then
          write(ID,'(A)') "Cartesian"
          icont=1
          do j=1,3
           rss(j)=rs(1)*a1(j) + rs(2)*a2(j) + rs(3)*a3(j)
          enddo
         endif
        enddo

        do i=1,natom
         if(iselect .eq. 1)then
          read(ID+10,*) coord(1:3),const(1:3)
          if(idirect .eq. 1)then
           if(coord(1)+rss(1) .ge. 1.) rss(1)= rss(1) - 1.
           if(coord(2)+rss(2) .ge. 1.) rss(2)= rss(2) - 1.
           if(coord(3)+rss(3) .ge. 1.) rss(3)= rss(3) - 1.
          endif

          write(ID,'(3F20.16,3X, A, A, A)') coord(:)+rss(:),const(1:3)
         elseif(iselect .eq. 0)then
          read(ID+10,*) coord(1:3)
          write(ID,'(3F20.16)') coord(:)+rss(:)
         endif
        enddo

        write(ID,*)" "
        write(ID,*)nx, ny, nz

       return
      end subroutine write_CHGCAR_head

!!$*  subroutine for setting k-points in BZ
      subroutine set_BZ_fukui(nnk,wnklist,w_half_klist,i_half_klist,
     &                   i_trim_klist,nhf, ispin,wklist,del,nk,iz2)

      implicit real*8 (a-h,o-z)
      real*8  wklist(3,nk),wnklist(3,nk*iz2),w_half_klist(3,nk)
      integer i_half_klist(nk),i_trim_klist(nk)
      dimension wk(3)

! set new kpoints : TR partners
      i_trim_klist=0
      do isp=1,ispin ! ispin start
       nnk=0;nhf=0
       do ik=1, nk !ik loop.
        wk(:)=wklist(:,ik)
        if(abs(wk(1)) .le. del .and. abs(wk(2)) .le. del)then
         nnk=nnk+1
         nhf=nhf+1
         wnklist(:,nnk)=wk(:)  !gamma
         w_half_klist(:,nhf)=wk(:) !save only B+ k-point's index
         i_half_klist(nhf)=ik
         i_trim_klist(nhf)=4
        elseif(abs(wk(1)-.5) .le. del .and. abs(wk(2)-.0) .le. del)then
         nnk=nnk+1
         nhf=nhf+1
         wnklist(:,nnk)=wk(:)  !M1
         w_half_klist(:,nhf)=wk(:)
         i_half_klist(nhf)=ik
         i_trim_klist(nhf)=1
        elseif(abs(wk(1)-.0) .le. del .and. abs(wk(2)-.5) .le. del)then
         nnk=nnk+1
         nhf=nhf+1
         wnklist(:,nnk)=wk(:)  !M2
         w_half_klist(:,nhf)=wk(:)
         i_half_klist(nhf)=ik
         i_trim_klist(nhf)=2
        elseif(abs(wk(1)-.5) .le. del .and. abs(wk(2)-.5) .le. del)then
         nnk=nnk+1
         nhf=nhf+1
         wnklist(:,nnk)=wk(:)  !M3
         w_half_klist(:,nhf)=wk(:)
         i_half_klist(nhf)=ik
         i_trim_klist(nhf)=3
        elseif(wk(2) .gt. .0 .and. wk(2) .lt. .5
     &                       .and. (.5-wk(2)) .gt. del
     &                       .and. (.5-wk(1)) .gt. del
     &                       .and. (wk(2)-.0) .gt. del)then
         nnk=nnk+1
         nhf=nhf+1
         wnklist(:,nnk)=wk(:) ! B+ interior
         w_half_klist(:,nhf)=wk(:)
         i_half_klist(nhf)=ik

         !find TR partner
         nnk=nnk+1
         wnklist(:,nnk)=-wk(:)

!       elseif(wk(1).gt. del .and. (.5-wk(1)) .gt. del
        elseif( wk(1) .le. -del .and. (abs(wk(2) -.5) .le. del))then
!    &                          .and. (abs(wk(2) -.5) .le. del))then
         nnk=nnk+1
         nhf=nhf+1
         wnklist(:,nnk)=wk(:)  !B+ edge (upper)
         w_half_klist(:,nhf)=wk(:)
         i_half_klist(nhf)=ik

         !find TR + T(G2) partner 
         nnk=nnk+1
         wnklist(:,nnk)=-wk(:)
         wnklist(2,nnk)=wnklist(2,nnk)+1.

        elseif(wk(1).gt. del .and. (.5-wk(1)) .gt. del
     &                       .and. (abs(wk(2) -.0) .le. del)) then
         nnk=nnk+1
         nhf=nhf+1
         wnklist(:,nnk)=wk(:)  !B+ edge (lower)
         w_half_klist(:,nhf)=wk(:)
         i_half_klist(nhf)=ik

         !find TR partner
         nnk=nnk+1
         wnklist(:,nnk)=-wk(:)

        elseif(abs(.5-wk(1)) .le. del .and. (.5-wk(2)) .gt. del
     &                                .and. (wk(2)-.0) .gt. del)then
         nnk=nnk+1
         nhf=nhf+1
         wnklist(:,nnk)=wk(:)  !B+ edge (right)
         w_half_klist(:,nhf)=wk(:)
         i_half_klist(nhf)=ik

         !find TR + T(G1) partner 
         nnk=nnk+1
         wnklist(:,nnk)=-wk(:)
         wnklist(1,nnk)=wnklist(1,nnk)+1.

        endif
       enddo !ik
      enddo !ispin end

!     do i=1,nnk
!     write(6,*)(wnklist(jj,i),jj=1,3)
!     enddo
!     stop
      return
      end subroutine set_BZ_fukui

!!$*  subroutine for get nfield strength
      subroutine get_nfield(rfield,wklp, isp,ik,nk,nband,nkx,nky,
     &                      wnklist,ihf,nplist,nnk,iz,ns,kperiod,
     &                      w_half_klist,i_half_klist,i_trim_klist,nhf,
     &                      nbmax,npmax,ecut,ispinor,ispin,ne,irecl,
     &                      nmax,nini,b1,b2,b3,iz2)
      implicit real*8 (a-h,o-z)
      complex*8  coeff(npmax)
      complex*16, allocatable :: Siju(:,:),Sijd(:,:),Sijt(:,:)
      dimension nbmax(3),wkk(3,5),wklp(3,5,nk*iz2)
      complex*16 coeff1u((2*nbmax(1)+2)*(2*nbmax(2)+2)*(2*nbmax(3)+1))
      complex*16 coeff1d((2*nbmax(1)+2)*(2*nbmax(2)+2)*(2*nbmax(3)+1))
      complex*16 coeff2u((2*nbmax(1)+2)*(2*nbmax(2)+2)*(2*nbmax(3)+1))
      complex*16 coeff2d((2*nbmax(1)+2)*(2*nbmax(2)+2)*(2*nbmax(3)+1))
      complex*16 detS(4),detA,detLOOP
      real*8  wnklist(3,nk*iz2),w_half_klist(3,nk),rfield
      dimension b1(3),b2(3),b3(3)
      integer nplist(nk),i_half_klist(nk),i_trim_klist(nk)
      integer itr(5),itrim(5),ikk(5),isgg(2,5),npl(5)
      integer ni,nj,ne,nk,nband,np,npmax,kperiod,ispin,irecl
      data c/0.262465831d0/ ! constant c = 2m/hbar**2 [1/eV Ang^2]
      pi=4.*atan(1.)

      call klpfind(wkk, isp,ik,nk,nband,nkx,nky,wnklist,ihf,iz,iz2) !k-loop(4pts) of ik-kpoint (K),counter-clock
      call kindxfind(ikk,isgg,npl,itr,itrim, wnklist,nplist,wkk,nnk,
     &               nband,iz,w_half_klist,i_half_klist,i_trim_klist,
     &               nhf,iz2) !find index

      wklp(:,:,ik) = wkk(:,:)  ! closed loop (C) for each K : RECIVEC

!!$* get overlap matrix S_ij(k,k+1) over the C and PI_s [detS(k_s,k_s+1)], s=1,4
        allocate(Siju(ns,ns),Sijd(ns,ns),Sijt(ns,ns))
        do ilp=1,4  ! loop for ilp
         Siju=(0.,0.)  !initialize
         Sijd=(0.,0.)
         Sijt=(0.,0.)
        ! construct overlap matrix S_ij(ikk1,ikk2)
        do ni=nini, nmax   ! calculate upto valence band maximum
         ncnt=0;coeff1u=(0.,0.);coeff1d=(0.,0.);coeff=(0.,0.)
         call get_coeff(coeff,isgg, ilp,ni,ikk,npl,nband,nk,isp,
     &                              npmax,itrim)
         do ig3=0,2*nbmax(3);    ig3p=ig3
          if (ig3.gt.nbmax(3))   ig3p=ig3-2*nbmax(3)-1
          do ig2=0,2*nbmax(2);   ig2p=ig2
           if (ig2.gt.nbmax(2))  ig2p=ig2-2*nbmax(2)-1
           do ig1=0,2*nbmax(1);  ig1p=ig1
            if (ig1.gt.nbmax(1)) ig1p=ig1-2*nbmax(1)-1
             call get_etot(etot, ilp,ni,b1,b2,b3,isgg,ig1p,ig2p,ig3p,
     &                           ig1,ig2,ig3,wkk,itrim,itr)
            if (etot.lt.ecut) then; ncnt=ncnt+1
             call get_incnt(incnt, ilp,ni,ig1p,ig2p,ig3p,nbmax,
     &                             itr,itrim,isgg,wkk)
            if(ispinor .eq. 2)then !spinor-dn
             if(itr(ilp) .eq. 0)then
              if(itrim(ilp) .eq. 0)then
               coeff1d(incnt)=coeff(ncnt+npl(ilp)/2) !for |u(dn)>
              elseif(itrim(ilp) .ge. 1 .and. mod(ni,2) .eq. 1)then
               coeff1d(incnt)=coeff(ncnt+npl(ilp)/2) !for |u(dn)>
              elseif(itrim(ilp) .ge. 1 .and. mod(ni,2) .eq. 0)then
               coeff1d(incnt)=conjg(coeff(ncnt)) !TRIM for |u*(2n-1,up)>
              endif
             elseif(itr(ilp) .eq. 1)then
              coeff1d(incnt)=conjg(coeff(ncnt)) !TRS for |u*(up)>
             endif
            endif

            if(itr(ilp) .eq. 0)then  !spinor-up
             if(itrim(ilp) .eq. 0)then
              coeff1u(incnt)=coeff(ncnt) !for |u(up)>
             elseif(itrim(ilp) .ge. 1 .and. mod(ni,2) .eq. 1)then
              coeff1u(incnt)=coeff(ncnt) !for |u(up)>
             elseif(itrim(ilp) .ge. 1 .and. mod(ni,2) .eq. 0)then
              if(ispinor .eq. 1)then
               coeff1u(incnt)=conjg(coeff(ncnt)) !TRIM for |u*(up)>
              elseif(ispinor .eq. 2)then
               coeff1u(incnt)=-conjg(coeff(ncnt+npl(ilp)/2)) !TRIM for -|u*(dn,2n-1)>
              endif
             endif
            elseif(itr(ilp) .eq. 1)then
             if(ispinor .eq. 2)then
              coeff1u(incnt)=-conjg(coeff(ncnt+npl(ilp)/2)) !TRS for -|u*(dn)>
             elseif(ispinor .eq. 1)then
              coeff1u(incnt)=conjg(coeff(ncnt)) !TRS for |u*>
             endif
            endif

            endif
           enddo  !loop for ig1                       
          enddo   !loop for ig2
         enddo    !loop for ig3
         if (ispinor*ncnt.ne.npl(ilp)) then
          write(0,*) '*ni error - computed NPL=',2*ncnt,
     &               ' != input nplane=',npl(ilp);stop
         endif
         do nj=nini, nmax
          ncnt=0;coeff2u=(0.,0.);coeff2d=(0.,0.);coeff=(0.,0.)
          call get_coeff(coeff,isgg, ilp+1,nj,ikk,npl,nband,nk,isp,
     &                               npmax,itrim)
          do ig3=0,2*nbmax(3); ig3p=ig3
           if (ig3.gt.nbmax(3)) ig3p=ig3-2*nbmax(3)-1
           do ig2=0,2*nbmax(2); ig2p=ig2
            if (ig2.gt.nbmax(2)) ig2p=ig2-2*nbmax(2)-1
            do ig1=0,2*nbmax(1); ig1p=ig1
             if (ig1.gt.nbmax(1)) ig1p=ig1-2*nbmax(1)-1
             call get_etot(etot, ilp+1,nj,b1,b2,b3,isgg,ig1p,ig2p,ig3p,
     &                           ig1,ig2,ig3,wkk,itrim,itr)
             if (etot.lt.ecut) then; ncnt=ncnt+1
              call get_incnt(incnt, ilp+1,nj,ig1p,ig2p,ig3p,nbmax,
     &                              itr,itrim,isgg,wkk)

             if(ispinor .eq. 2)then  !spinor-dn
              if(itr(ilp+1) .eq. 0)then
               if(itrim(ilp+1) .ge. 0)then
                coeff2d(incnt)=coeff(ncnt+npl(ilp+1)/2) !c(dn,n)
               elseif(itrim(ilp+1) .ge. 1 .and. mod(nj,2) .eq. 1)then
                coeff2d(incnt)=coeff(ncnt+npl(ilp+1)/2) !c(dn,n)
               elseif(itrim(ilp+1) .ge. 1 .and. mod(nj,2) .eq. 0)then
                coeff2d(incnt)=conjg(coeff(ncnt)) !c*(up,n-1)
               endif
              elseif(itr(ilp+1) .eq. 1)then
               coeff2d(incnt)=conjg(coeff(ncnt)) !c*(up)
              endif
             endif

             if(itr(ilp+1) .eq. 0)then  !spinor-up 
              if(itrim(ilp+1) .eq. 0)then
               coeff2u(incnt)=coeff(ncnt)
              elseif(itrim(ilp+1) .ge. 1 .and. mod(nj,2) .eq. 1)then
               coeff2u(incnt)=coeff(ncnt)
              elseif(itrim(ilp+1) .ge. 1 .and. mod(nj,2) .eq. 0)then
               if(ispinor .eq. 1)then
                coeff2u(incnt)=conjg(coeff(ncnt))  !c(n-1)*
               elseif(ispinor .eq. 0)then
                coeff2u(incnt)=-conjg(coeff(ncnt+npl(ilp+1)/2)) !-c*(dn,n-1)
               endif
              endif
             elseif(itr(ilp+1) .eq. 1)then
              if(ispinor .eq. 2)then
               coeff2u(incnt)=-conjg(coeff(ncnt+npl(ilp+1)/2))
              elseif(ispinor .eq. 1)then
               coeff2u(incnt)=conjg(coeff(ncnt))
              endif
             endif

             endif
            enddo   !loop for ig1
           enddo    !loop for ig2
          enddo     !loop for ig3
          if (ispinor*ncnt.ne.npl(ilp+1)) then
           write(0,*) '*nj error - computed NPL=',2*ncnt,
     &                ' != input nplane=',npl(ilp+1);stop
          endif
          if(ispinor .eq. 2)then
           Siju(ni-nmax+ns,nj-nmax+ns)=dot_product(coeff1u,coeff2u)
           Sijd(ni-nmax+ns,nj-nmax+ns)=dot_product(coeff1d,coeff2d)
           Sijt(ni-nmax+ns,nj-nmax+ns)=Siju(ni-nmax+ns,nj-nmax+ns)+
     &                                 Sijd(ni-nmax+ns,nj-nmax+ns)
           else if (ispinor .eq. 1)then
            Sijt(ni-nmax+ns,nj-nmax+ns)=dot_product(coeff1u,coeff2u)
          endif
         enddo !nj
        enddo !ni

        if(nini .eq. nmax)then   ! get determinant : det(S)
         detS(ilp)=Sijt(1,1)
         else
!         call getdetA(detA, Sijt,ns)
          call get_det(detA, Sijt,ns)
          detS(ilp)=detA;detA=(0.,0.)
        endif
        enddo !ilp

        detLOOP=detS(1)*detS(2)*detS(3)*detS(4)
        rfield=+sum(aimag(log(detS(1:4))))/(2.*pi)
     &               - (aimag(log(detLOOP)))/(2.*pi)
!       write(6,'(A)')"# "
! 704   format('#===>PI_S[det(S(K_s,K_s+1))] =',F16.8,'  +',F16.8,' i'
!    &         ,I4,'th-K')
        if(sum(itrim(1:4)) .ge. 1)then 
         if(itrim(1)+itrim(3) .ge. 1)nparity=1
         if(itrim(2)+itrim(4) .ge. 1)nparity=-1
         itr_sum=sum(itrim(1:4))*nparity
        else
         itr_sum=0
        endif
        write(6,704)ik,rfield,log(detLOOP),ikk(:),wkk(:,ilp),
     &              (wkk(:,1)+wkk(:,3))/2.,itr_sum,sum(log(detS(:)))
  704   format('#IK=',I4,' NK=',F5.2,' P(detS)=',F16.8,'+',F16.8,'i',
     &         ' KLP=',5I4,' WK=',3F8.4, ' WKC=',3F8.4,' STRM=',I2,
     &         " S(detS)=",F16.8,'+',F16.8,'i')

      deallocate(Sijt)
      deallocate(Siju)
      deallocate(Sijd)

      return
      end subroutine get_nfield

!!$*  subroutine for getting planewave coefficient for given k & n
      subroutine get_coeff(coeff,isgg, iilp,nni,ikk,npl,nband,nk,isp,
     &                                npmax,itrim)

      implicit real*8 (a-h,o-z)  
      dimension itrim(5),isgg(2,5),ikk(5)
      complex*8  coeff(npmax)
      integer    npl(5)
      irecl=3+(ikk(iilp)-1)*(nband+1)+nk*(nband+1)*(isp-1)+nni
      if(itrim(iilp) .eq. 0)then
       read(10,rec=irecl)(coeff(i),i=1,npl(iilp))
      elseif(itrim(iilp) .ge. 1 .and. mod(nni,2) .eq. 1)then
       read(10,rec=irecl)(coeff(i),i=1,npl(iilp))
      elseif(itrim(iilp) .ge. 1 .and. mod(nni,2) .eq. 0)then
       read(10,rec=irecl-1)(coeff(i),i=1,npl(iilp))
      endif


      return
      end subroutine

!!$*  subroutine for getting incnt with given G vector 
      subroutine get_incnt(incnt, iilp,nni,ig1p,ig2p,ig3p,nbmax,
     &                            itr,itrim,isgg,wkk)
      implicit real*8 (a-h,o-z)
      dimension nbmax(3),itrim(5),itr(5),isgg(2,5)
      real*8    wkk(3,5)
      data c/0.262465831d0/ ! constant c = 2m/hbar**2 [1/eV Ang^2]
      if(itr(iilp) .eq. 0)then
       if(itrim(iilp) .eq. 0 )then
        incnt=(ig3p+nbmax(3))*(2*nbmax(2)+1)*(2*nbmax(1)+1)+
     &        (ig2p-isgg(2,iilp)+nbmax(2))*(2*nbmax(1)+1)+
     &        (ig1p-isgg(1,iilp)+nbmax(1))+1
       elseif(itrim(iilp) .ge. 1 .and. mod(nni,2) .eq. 1)then  ! TRIM & 2N-1 state
        incnt=(+ig3p+nbmax(3))*(2*nbmax(2)+1)*(2*nbmax(1)+1)+
     &        (+ig2p-isgg(2,iilp)*0+nbmax(2))*(2*nbmax(1)+1)+
     &        (+ig1p-isgg(1,iilp)*0+nbmax(1))+1
       elseif(itrim(iilp) .ge. 1 .and. mod(nni,2) .eq. 0)then  ! TRIM & 2N state
        incnt=(-ig3p+nbmax(3))*(2*nbmax(2)+1)*(2*nbmax(1)+1)+
     &        (-ig2p-isgg(2,iilp)+nbmax(2))*(2*nbmax(1)+1)+
     &        (-ig1p-isgg(1,iilp)+nbmax(1))+1
       endif
      elseif(itr(iilp) .eq. 1)then !TIME REVERSAL
       incnt=(-ig3p+nbmax(3))*(2*nbmax(2)+1)*(2*nbmax(1)+1)+
     &       (-ig2p-isgg(2,iilp)+nbmax(2))*(2*nbmax(1)+1)+
     &       (-ig1p-isgg(1,iilp)+nbmax(1))+1
      endif

      return
      end subroutine

!!$*  subroutine for getting etot within given K & G vector
      subroutine get_etot(etot, iilp,nni,b1,b2,b3,isgg,ig1p,ig2p,ig3p,
     &                          ig1,ig2,ig3,wkk,itrim,itr)

      implicit real*8 (a-h,o-z)
      dimension itrim(5),itr(5),isgg(2,5),sumkg(3)
      real*8    b1(3),b2(3),b3(3),wkk(3,5)
      data c/0.262465831d0/ ! constant c = 2m/hbar**2 [1/eV Ang^2]
      do j=1,3
       if(itr(iilp) .eq. 0) then
        if(itrim(iilp) .eq. 0)then
         sumkg(j)=(wkk(1,iilp)+ig1p-isgg(1,iilp))*b1(j)+
     &            (wkk(2,iilp)+ig2p-isgg(2,iilp))*b2(j)+
     &            (wkk(3,iilp)+ig3p)*b3(j)
        elseif(itrim(iilp) .ge. 1 .and. mod(nni,2) .eq. 1)then
         sumkg(j)=(wkk(1,iilp)+ig1p-isgg(1,iilp)*0)*b1(j)+
     &            (wkk(2,iilp)+ig2p-isgg(2,iilp)*0)*b2(j)+
     &            (wkk(3,iilp)+ig3p)*b3(j)
        elseif(itrim(iilp) .ge. 1 .and. mod(nni,2) .eq. 0)then
         sumkg(j)=(wkk(1,iilp)-ig1p-isgg(1,iilp)*0)*b1(j)+
     &            (wkk(2,iilp)-ig2p-isgg(2,iilp)*0)*b2(j)+
     &            (wkk(3,iilp)-ig3p)*b3(j)
        endif
       elseif(itr(iilp) .eq. 1)then !TIME REVERSAL
        sumkg(j)=(-wkk(1,iilp)-ig1p-isgg(1,iilp))*b1(j)+
     &           (-wkk(2,iilp)-ig2p-isgg(2,iilp))*b2(j)+
     &           (-wkk(3,iilp)-ig3p)*b3(j)
       endif
      enddo
      gtot=sqrt(dot_product(sumkg,sumkg))
      etot=gtot**2/c
      
      return
      end subroutine

!!$*  subroutine for computing optical selectivity in the given k-point
      subroutine optical_selectivity(selectivity,
     &                               selectivitymax,selectivitymin,
     &           b1,b2,b3,wklist,isp,
     &           nband,ecut,ispinor,nplist,nbmax,npmax,nk,nini,nmax)
      implicit real*8 (a-h,o-z)
      dimension nbmax(3),nplist(nk),wk(3)
      real*8    selectivity(nk),b1(3),b2(3),b3(3),wklist(3,nk)
      real*8    selectivitymax(4),selectivitymin(4)
      complex*16 ctrans_mtrx_left,ctrans_mtrx_right
      complex*16 cinter_mtrx_x,cinter_mtrx_y
      complex*16 coeffv(npmax),coeffc(npmax)
      complex*8  coeff(npmax)
      complex*16 coeffvu(npmax),coeffvd(npmax)
      complex*16 coeffcu(npmax),coeffcd(npmax)
      integer :: ig(3,npmax)
!     data hbar/6.58211928e-16/ 
      data hbar/1./ !Here, I will set hbar = 1. for the simplicity
      do ik=1,nk
       cinter_mtrx_x=(0.,0.)
       cinter_mtrx_y=(0.,0.)
       coeff=(0.,0.)
       coeffv=(0.,0.);coeffc=(0.,0.)
       coeffvu=(0.,0.);coeffcu=(0.,0.)
       coeffvd=(0.,0.);coeffcd=(0.,0.)
       wk(:)=wklist(:,ik)
       np=nplist(ik)
       call plindx(ig,ncnt, ispinor,wk,b1,b2,b3,nbmax,np,ecut,npmax)
       read(10,rec=(3+(ik-1)*(nband+1)+
     &                nk*(nband+1)*(isp-1)+nini))(coeff(i),i=1,np)
       coeffv=coeff;coeff=(0.,0.)
       read(10,rec=(3+(ik-1)*(nband+1)+
     &                nk*(nband+1)*(isp-1)+nmax))(coeff(i),i=1,np)
       coeffc=coeff;coeff=(0.,0.)
       do iplane=1,ncnt
        xkgx=(wk(1)+ig(1,iplane))*b1(1)+
     &       (wk(2)+ig(2,iplane))*b2(1)+
     &       (wk(3)+ig(3,iplane))*b3(1)
        xkgy=(wk(1)+ig(1,iplane))*b1(2)+
     &       (wk(2)+ig(2,iplane))*b2(2)+
     &       (wk(3)+ig(3,iplane))*b3(2)
        if(ispinor .eq. 2) then
         coeffvu(iplane)=coeffv(iplane)
         coeffvd(iplane)=coeffv(iplane+ncnt)
         coeffcu(iplane)=coeffc(iplane)
         coeffcd(iplane)=coeffc(iplane+ncnt)
         cinter_mtrx_x=cinter_mtrx_x+
     &            hbar*conjg(coeffcu(iplane))*xkgx*coeffvu(iplane)+
     &            hbar*conjg(coeffcd(iplane))*xkgx*coeffvd(iplane)
         cinter_mtrx_y=cinter_mtrx_y+
     &            hbar*conjg(coeffcu(iplane))*xkgy*coeffvu(iplane)+
     &            hbar*conjg(coeffcd(iplane))*xkgy*coeffvd(iplane)
         else if(ispinor .eq. 1) then
          cinter_mtrx_x=cinter_mtrx_x+
     &            hbar*conjg(coeffc(iplane))*xkgx*coeffv(iplane)
          cinter_mtrx_y=cinter_mtrx_y+
     &            hbar*conjg(coeffc(iplane))*xkgy*coeffv(iplane)
        endif
       enddo ! iplane loop end

       ctrans_mtrx_left  = cinter_mtrx_x + (0.,1.)*cinter_mtrx_y
       ctrans_mtrx_right = cinter_mtrx_x - (0.,1.)*cinter_mtrx_y
       selectivity(ik) = 
     &   ((abs(ctrans_mtrx_left))**2 - (abs(ctrans_mtrx_right))**2)/
     &   ((abs(ctrans_mtrx_left))**2 + (abs(ctrans_mtrx_right))**2)
       write(6,'(A,I4,4F11.6)')"# IK, K(reci), SELECTIVITY(n) : ",
     &                         ik,wk,selectivity(ik)
       if (ik. eq. 1)then
        selectivitymax(4)=selectivity(ik)
        selectivitymax(1)=wklist(1,ik)
        selectivitymax(2)=wklist(2,ik)
        selectivitymax(3)=wklist(3,ik)
        selectivitymin(4)=selectivity(ik)
        selectivitymin(1)=wklist(1,ik)
        selectivitymin(2)=wklist(2,ik)
        selectivitymin(3)=wklist(3,ik)
        else if (ik.ge.2.and.selectivity(ik).ge.selectivitymax(4))then
         selectivitymax(4)=selectivity(ik)
         selectivitymax(1)=wklist(1,ik)
         selectivitymax(2)=wklist(2,ik)
         selectivitymax(3)=wklist(3,ik)
        else if (ik.ge.2.and.selectivity(ik).le.selectivitymin(4))then
         selectivitymin(4)=selectivity(ik)
         selectivitymin(1)=wklist(1,ik)
         selectivitymin(2)=wklist(2,ik)
         selectivitymin(3)=wklist(3,ik)
       endif
      enddo !ik loop end

      return
      end subroutine optical_selectivity

!!$*  subroutine for computing velocity expectation value for state psi(n,k)
      subroutine vel_expectation(a1,a2,a3,b1,b2,b3,
     &   kperiod,nbmax,npmax,ecut,ispinor,ispin,nband,ne,nk,wklist,
     &   nplist,ni,nj,irecl,filename,foname)
      implicit real*8 (a-h, o-z)
      complex*8  coeff(npmax)
      complex*16 coeffi(npmax),coeffj(npmax)
      complex*16 coeffiu(npmax),coeffid(npmax)
      complex*16 coeffju(npmax),coeffjd(npmax)
      complex*16 vel_x,vel_y
      real*8  vel_expt(2,nk,ispin),wklist(3,nk) 
      real*8  recivec(3,nk),recilat(3,nk)
      real*8  xvel_expt(2,kperiod*2*kperiod*2*nk,ispin)
      real*8  xrecivec(3,kperiod*2*kperiod*2*nk)
      real*8  xrecilat(3,kperiod*2*kperiod*2*nk)
      real*8  vel_x_expt_max(4),vel_y_expt_max(4)
      real*8  vel_x_expt_min(4),vel_y_expt_min(4)
      dimension a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),a2xa3(3),sumkg(3)
      dimension wk(3),nbmax(3),xb(3)
      integer ig(3,npmax),nplist(nk)
      integer ni,nj,ne,nk,nband,np,npmax,kperiod,ispin,irecl
      character*75 filename,foname,fonameo
      data c/0.262465831d0/ ! constant c = 2m/hbar**2 [1/eV Ang^2]
      data hbar/6.58211928E-16/ !h/2pi [eV * s]
      data xm/0.510998910E+6/ ! electron mass (eV/c^2)
      data anginv/1.0E+10/     ! inverse angstrom (1/Ang)
      pi=4.*atan(1.)
      do isp=1,ispin
       do ik=1, nk
        coeff=(0.,0.)
        coeffi=(0.,0.);coeffj=(0.,0.)
        coeffiu=(0.,0.);coeffid=(0.,0.)
        coeffju=(0.,0.);coeffjd=(0.,0.)
        wk(:)=wklist(:,ik)
        np=nplist(ik)
        if(ni .ne. nj)then
         call plindx(ig,ncnt, ispinor,wk,b1,b2,b3,nbmax,np,ecut,npmax)
         read(10,rec=(3+(ik-1)*(nband+1)+
     &                nk*(nband+1)*(isp-1)+ni))(coeff(i),i=1,np)
         coeffi=coeff;coeff=(0.,0.)
         read(10,rec=(3+(ik-1)*(nband+1)+
     &                nk*(nband+1)*(isp-1)+nj))(coeff(i),i=1,np)
         coeffj=coeff;coeff=(0.,0.)
         else if(ni .eq. nj)then
          call plindx(ig,ncnt, ispinor,wk,b1,b2,b3,nbmax,np,ecut,npmax)
          read(10,rec=(3+(ik-1)*(nband+1)+
     &                nk*(nband+1)*(isp-1)+ni))(coeff(i),i=1,np)
          coeffi=coeff;coeff=(0.,0.)
        endif
        vel_x=(0.,0.);vel_y=(0.,0.)
        do iplane=1,ncnt
         xkgx=(wk(1)+ig(1,iplane))*b1(1)+
     &        (wk(2)+ig(2,iplane))*b2(1)+
     &        (wk(3)+ig(3,iplane))*b3(1)
         xkgy=(wk(1)+ig(1,iplane))*b1(2)+
     &        (wk(2)+ig(2,iplane))*b2(2)+
     &        (wk(3)+ig(3,iplane))*b3(2)
         if(ispinor .eq. 2) then
          coeffiu(iplane)=coeffi(iplane)
          coeffid(iplane)=coeffi(iplane+ncnt)
          if(ni .eq. nj)then
           coeffju(iplane)=coeffi(iplane)
           coeffjd(iplane)=coeffi(iplane+ncnt)
           else
           coeffju(iplane)=coeffj(iplane)
           coeffjd(iplane)=coeffj(iplane+ncnt)
          endif
          ! -i*hbar/m <psi|d/dx|psi>
          vel_x=vel_x+     
     &    anginv*hbar/xm*(conjg(coeffiu(iplane))*xkgx*coeffju(iplane)+
     &                    conjg(coeffid(iplane))*xkgx*coeffjd(iplane))
          ! -i*hbar/m <psi|d/dy|psi> 
          vel_y=vel_y+     
     &    anginv*hbar/xm*(conjg(coeffiu(iplane))*xkgy*coeffju(iplane)+
     &                    conjg(coeffid(iplane))*xkgy*coeffjd(iplane))
          else if (ispinor .eq. 1) then
           coeffiu(iplane)=coeffi(iplane)
           if(ni .eq. nj)then
            coeffju(iplane)=coeffi(iplane)
            else
             coeffju(iplane)=coeffj(iplane)
           endif
           vel_x=vel_x+
     &     anginv*hbar/xm*conjg(coeffiu(iplane))*xkgx*coeffju(iplane)
           vel_y=vel_y+
     &     anginv*hbar/xm*conjg(coeffiu(iplane))*xkgy*coeffju(iplane)
         endif ! ispinor
        enddo  ! iplane
        vel_expt(1,ik,isp)=real(vel_x,8)
        vel_expt(2,ik,isp)=real(vel_y,8)
        write(6,'(A,I4,5F11.6)')"# IK, K(reci), VEL_EXPT(x,y) : ",
     &                         ik,wk,(vel_expt(i,ik,isp),i=1,2)
       enddo  !ik

       if (ik .eq. 1)then
        vel_x_expt_max(4)=vel_expt(1,ik,isp)
        vel_x_expt_max(1)=wklist(1,ik)
        vel_x_expt_max(2)=wklist(2,ik)
        vel_x_expt_max(3)=wklist(3,ik)
        vel_x_expt_min(4)=vel_expt(1,ik,isp)
        vel_x_expt_min(1)=wklist(1,ik)
        vel_x_expt_min(2)=wklist(2,ik)
        vel_x_expt_min(3)=wklist(3,ik)
        else if(ik.ge.2.and.vel_expt(1,ik,isp).ge.vel_x_expt_max(4))then
         vel_x_expt_max(4)=vel_expt(1,ik,isp)
         vel_x_expt_max(1)=wklist(1,ik)
         vel_x_expt_max(2)=wklist(2,ik)
         vel_x_expt_max(3)=wklist(3,ik)
        else if(ik.ge.2.and.vel_expt(1,ik,isp).le.vel_x_expt_min(4))then
         vel_x_expt_min(4)=vel_expt(1,ik,isp)
         vel_x_expt_min(1)=wklist(1,ik)
         vel_x_expt_min(2)=wklist(2,ik)
         vel_x_expt_min(3)=wklist(3,ik)
       endif

       if (ik .eq. 1)then
        vel_y_expt_max(4)=vel_expt(2,ik,isp)
        vel_y_expt_max(1)=wklist(1,ik)
        vel_y_expt_max(2)=wklist(2,ik)
        vel_y_expt_max(3)=wklist(3,ik)
        vel_y_expt_min(4)=vel_expt(2,ik,isp)
        vel_y_expt_min(1)=wklist(1,ik)
        vel_y_expt_min(2)=wklist(2,ik)
        vel_y_expt_min(3)=wklist(3,ik)
        else if(ik.ge.2.and.vel_expt(2,ik,isp).ge.vel_y_expt_max(4))then
         vel_y_expt_max(4)=vel_expt(2,ik,isp)
         vel_y_expt_max(1)=wklist(1,ik)
         vel_y_expt_max(2)=wklist(2,ik)
         vel_y_expt_max(3)=wklist(3,ik)
        else if(ik.ge.2.and.vel_expt(2,ik,isp).le.vel_y_expt_min(4))then
         vel_y_expt_min(4)=vel_expt(2,ik,isp)
         vel_y_expt_min(1)=wklist(1,ik)
         vel_y_expt_min(2)=wklist(2,ik)
         vel_y_expt_min(3)=wklist(3,ik)
       endif
        do ik=1,nk
         do j=1,3
          recivec(j,ik)=wklist(1,ik)*b1(j)+
     &                  wklist(2,ik)*b2(j)+
     &                  wklist(3,ik)*b3(j)
          recilat(j,ik)=wklist(j,ik)
         enddo
        enddo

        kk=0   ! extend berry curvature distribution over extended BZ 
        do ib2=-1*(kperiod-1)+1,kperiod
         do ib1=-1*(kperiod-1)+1,kperiod  ! you may adjust these values as you wish..
          do ik=1,nk
           kk=kk+1
           xrecivec(1,kk)=recivec(1,ik) +(ib1-1)*b1(1)+(ib2-1)*b2(1)
           xrecivec(2,kk)=recivec(2,ik) +(ib1-1)*b1(2)+(ib2-1)*b2(2)
           xrecivec(3,kk)=recivec(3,ik) +(ib1-1)*b1(3)+(ib2-1)*b2(3)
           xrecilat(1,kk)=recilat(1,ik) +(ib1-1)
           xrecilat(2,kk)=recilat(2,ik) +(ib2-1)
           xrecilat(3,kk)=recilat(3,ik)
           xvel_expt(1,kk,isp)=vel_expt(1,ik,isp)
           xvel_expt(2,kk,isp)=vel_expt(2,ik,isp)
           kext=kk
          enddo
         enddo
        enddo

!$$*  sorting k-points and the corresponding optical selectivity on the
!     periodically repeated data grid
        write(6,*)" "
        write(6,'(A)')"# SORTING K-grids..."
        do k=kext-1,1,-1  ! sorting kx
         do j=1,k
          if(xrecivec(1,j+1) .gt. xrecivec(1,j))then
           xb(:)=xrecivec(:,j)
           xrecivec(:,j)=xrecivec(:,j+1)
           xrecivec(:,j+1)=xb(:)
           xb(:)=xrecilat(:,j)
           xrecilat(:,j)=xrecilat(:,j+1)
           xrecilat(:,j+1)=xb(:)
           xtemp=xvel_expt(1,j,isp)
           xvel_expt(1,j,isp)=xvel_expt(1,j+1,isp)
           xvel_expt(1,j+1,isp)=xtemp
           xtemp=xvel_expt(2,j,isp)
           xvel_expt(2,j,isp)=xvel_expt(2,j+1,isp)
           xvel_expt(2,j+1,isp)=xtemp
          endif
         enddo
        enddo
        do k=kext-1,1,-1  ! sorting ky
         do j=1,k
          if(xrecivec(1,j+1) .eq. xrecivec(1,j))then
           if(xrecivec(2,j+1) .gt. xrecivec(2,j))then
           xb(:)=xrecivec(:,j)
           xrecivec(:,j)=xrecivec(:,j+1)
           xrecivec(:,j+1)=xb(:)
           xb(:)=xrecilat(:,j)
           xrecilat(:,j)=xrecilat(:,j+1)
           xrecilat(:,j+1)=xb(:)
           xtemp=xvel_expt(1,j,isp)
           xvel_expt(1,j,isp)=xvel_expt(1,j+1,isp)
           xvel_expt(1,j+1,isp)=xtemp
           xtemp=xvel_expt(2,j,isp)
           xvel_expt(2,j,isp)=xvel_expt(2,j+1,isp)
           xvel_expt(2,j+1,isp)=xtemp
           endif
          endif
         enddo
        enddo

        if(isp .eq. 1 .and. ispinor .eq. 1 .and. ispin.eq.2) then
         write(fonameo,'(A,A)')TRIM(foname),'.UP.dat'
         else if (isp .eq. 2 .and. ispinor .eq. 1 .and. ispin.eq.2)then
          write(fonameo,'(A,A)')TRIM(foname),'.DN.dat'
         else if (isp .eq. 1 .and. ispinor .eq. 2) then
          write(fonameo,'(A,A)')TRIM(foname),'.dat'
         else if (isp .eq. 1 .and. ispinor .eq. 1 .and. ispin.eq.1)then
          write(fonameo,'(A,A)')TRIM(foname),'.dat'
        endif
        open(61,file=fonameo,status='unknown')
        write(61,'(A,A)')"# File reading... : ",filename
        write(61,'(A,I9)')"# TOTAL RECORD LENGTH = ",irecl
        if (ispinor .eq. 2)then
         write(61,'(A,I6,A)')"# ISPIN            : ",ispin,
     &                       " (LSORBIT = .TRUE.)"
         else
          write(61,'(A,I6,A)')"# ISPIN            : ",ispin,
     &                        " (LSORBIT = .FALSE.)"
        endif
        write(61,'(A,F11.4)')  "# ENCUT (eV)       : ",ecut
        write(61,'(A,I6)')     "# NKPOINT          : ",nk
        write(61,'(A,I6)')     "# NBANDS           : ",nband
        write(61,'(A,3F13.6)') "# RECIVEC B1 (A^-1): ",(b1(i),i=1,3)
        write(61,'(A,3F13.6)') "# RECIVEC B2       : ",(b2(i),i=1,3)
        write(61,'(A,3F13.6)') "# RECIVEC B3       : ",(b3(i),i=1,3)
        write(61,*)" "

        write(61,'(A,I4)')"# VEOLOCITY EXPECTATION VALUE of BAND:",ni
        write(61,'(A)')"# <v(n,k)>= 1/hbar dE(k)/dk =
     & 1/m_e<psi(n,k)|p|psi(n,k)>, p=-i*hbar*d/dx,m_e=elect_rest_mass"
        write(61,'(A,4F16.6)')"# MAXVAL of VEL_EXPT <v_x> 
     &(in reci)= ",(vel_x_expt_max(i),i=1,4)
        write(61,'(A,4F16.6)')"# MINVAL of VEL_EXPT <v_x> 
     &(in reci)= ",(vel_x_expt_min(i),i=1,4)
        write(61,'(A,4F16.6)')"# MAXVAL of VEL_EXPT <v_y> 
     &(in reci)= ",(vel_y_expt_max(i),i=1,4)
        write(61,'(A,4F16.6)')"# MINVAL of VEL_EXPT <v_y> 
     &(in reci)= ",(vel_y_expt_min(i),i=1,4)
        write(61,'(A)')"# (cart) kx     ky     kz(A^-1)
     &    vel_expt(vx(n,k), vy(n,k))(m/s)  (recip)kx      ky      kz"
        do ik=1,kext
         write(61,'(3F11.6,A,2F10.6,A,3F11.6)')(xrecivec(i,ik),i=1,3),
     &        "   ",(xvel_expt(i,ik,isp),i=1,2),"       ",
     &        (xrecilat(i,ik),i=1,3)
        enddo
        close(61)

      enddo   !isp

!     write(6,'(A)')"# DONE! "
!     do isp=1, ispin
!      if(isp .eq. 1 .and. ispinor .eq. 1 .and. ispin .eq. 2) then
!       write(6,'(A,A,A)')"#  Results are summarized in ",TRIM(foname),
!    &                    ".UP.dat for spin-1"
!       else if (isp.eq.2 .and. ispinor.eq.1 .and. ispin.eq.2)then
!        write(6,'(A,A,A)')"#  Results are summarized in ",
!    &                     TRIM(foname),".DN.dat for spin-2"
!       else if (isp .eq. 1 .and. ispinor .eq. 2) then
!        write(6,'(A,A,A)')"#  Results are summarized in ",
!    &                     TRIM(foname),".dat"
!       else if (isp.eq.1 .and. ispinor.eq.1 .and. ispin.eq.1) then
!        write(6,'(A,A,A)')"#  Results are summarized in ",
!    &                     TRIM(foname),".dat"
!      endif
!     enddo
      close(10)
      return
      end subroutine vel_expectation


!!$*  subroutine for computing reciprocal lattice vector
      subroutine recilatt(b1,b2,b3,dSkxky, a1,a2,a3,nkx,nky)
      implicit real*8 (a-h,o-z)
      dimension a1(3),a2(3),a3(3), b1(3),b2(3),b3(3),a2xa3(3)
      dimension b1xb2(3)
      pi=4.*atan(1.)
      call vcross(a2xa3,a2,a3)
      Vcell=dot_product(a1,a2xa3)
      a3mag=dsqrt(dot_product(a3,a3))
      call vcross(b1,a2,a3);call vcross(b2,a3,a1);call vcross(b3,a1,a2)
      b1=2.*pi*b1/Vcell ; b2=2.*pi*b2/Vcell ; b3=2.*pi*b3/Vcell

      call vcross(b1xb2,b1,b2)
      dSkxky=dsqrt(dot_product(b1xb2,b1xb2))/real(nkx)/real(nky)
      return
      end subroutine recilatt

!!$*  subroutine for computing vector cross-product
      subroutine vcross(a,b,c)
      implicit real*8(a-h,o-z)
      dimension a(3),b(3),c(3)
      a(1)=b(2)*c(3)-b(3)*c(2)
      a(2)=b(3)*c(1)-b(1)*c(3)
      a(3)=b(1)*c(2)-b(2)*c(1)
      return
      end subroutine vcross

!!$*  subroutine for computing reciprocal properties
      subroutine reciproperty(nbmax,npmax, b1,b2,b3,ecut,ispinor)
      implicit real*8(a-h,o-z)
      dimension b1(3),b2(3),b3(3),vtmp(3),nbmax(3)
      data c/0.262465831d0/
      pi=4.*atan(1.)

      b1mag=dsqrt(b1(1)**2+b1(2)**2+b1(3)**2)
      b2mag=dsqrt(b2(1)**2+b2(2)**2+b2(3)**2)
      b3mag=dsqrt(b3(1)**2+b3(2)**2+b3(3)**2)

      phi12=acos((b1(1)*b2(1)+b1(2)*b2(2)+b1(3)*b2(3))/(b1mag*b2mag))
      call vcross(vtmp,b1,b2)
      vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
      sinphi123=(b3(1)*vtmp(1)+b3(2)*vtmp(2)+b3(3)*vtmp(3))/(vmag*b3mag)
      nb1maxA=(dsqrt(ecut*c)/(b1mag*abs(sin(phi12))))+1
      nb2maxA=(dsqrt(ecut*c)/(b2mag*abs(sin(phi12))))+1
      nb3maxA=(dsqrt(ecut*c)/(b3mag*abs(sinphi123)))+1
      npmaxA=nint(4.*pi*nb1maxA*nb2maxA*nb3maxA/3.)

      phi13=acos((b1(1)*b3(1)+b1(2)*b3(2)+b1(3)*b3(3))/(b1mag*b3mag))
      call vcross(vtmp,b1,b3)
      vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
      sinphi123=(b2(1)*vtmp(1)+b2(2)*vtmp(2)+b2(3)*vtmp(3))/(vmag*b2mag)
      phi123=abs(asin(sinphi123))
      nb1maxB=(dsqrt(ecut*c)/(b1mag*abs(sin(phi13))))+1
      nb2maxB=(dsqrt(ecut*c)/(b2mag*abs(sinphi123)))+1
      nb3maxB=(dsqrt(ecut*c)/(b3mag*abs(sin(phi13))))+1
      npmaxB=nint(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)

      phi23=acos((b2(1)*b3(1)+b2(2)*b3(2)+b2(3)*b3(3))/(b2mag*b3mag))
      call vcross(vtmp,b2,b3)
      vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
      sinphi123=(b1(1)*vtmp(1)+b1(2)*vtmp(2)+b1(3)*vtmp(3))/(vmag*b1mag)
      phi123=abs(asin(sinphi123))
      nb1maxC=(dsqrt(ecut*c)/(b1mag*abs(sinphi123)))+1
      nb2maxC=(dsqrt(ecut*c)/(b2mag*abs(sin(phi23))))+1
      nb3maxC=(dsqrt(ecut*c)/(b3mag*abs(sin(phi23))))+1
      npmaxC=nint(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)

      nbmax(1)=max0(nb1maxA,nb1maxB,nb1maxC)       ! maximum 
      nbmax(2)=max0(nb2maxA,nb2maxB,nb2maxC)
      nbmax(3)=max0(nb3maxA,nb3maxB,nb3maxC)

      !! multiply 'ispinor' to handle two component spinors
      npmax=ispinor*min0(npmaxA,npmaxB,npmaxC)
      return
      end subroutine reciproperty
      
!!$*  subroutine for computing planewave G index
      subroutine plindx(ig,ncnt,
     &    ispinor,wk,b1,b2,b3,nbmax,np,ecut,npmax)
      implicit real*8(a-h,o-z)
      dimension wk(3),sumkg(3),b1(3),b2(3),b3(3),nbmax(3)
      integer :: ig(3,npmax)
      data c/0.262465831d0/
      ncnt=0
      do ig3=0,2*nbmax(3)
       ig3p=ig3
       if (ig3.gt.nbmax(3)) ig3p=ig3-2*nbmax(3)-1
       do ig2=0,2*nbmax(2)
       ig2p=ig2
        if (ig2.gt.nbmax(2)) ig2p=ig2-2*nbmax(2)-1
        do ig1=0,2*nbmax(1)
        ig1p=ig1
         if (ig1.gt.nbmax(1)) ig1p=ig1-2*nbmax(1)-1
         do j=1,3
          sumkg(j)=(wk(1)+ig1p)*b1(j)+
     &             (wk(2)+ig2p)*b2(j)+
     &             (wk(3)+ig3p)*b3(j)
         enddo
         gtot=sqrt(dot_product(sumkg,sumkg))
         etot=gtot**2/c
         if (etot.lt.ecut) then
          ncnt=ncnt+1
          ig(1,ncnt)=ig1p
          ig(2,ncnt)=ig2p
          ig(3,ncnt)=ig3p
         end if
        enddo
       enddo
      enddo
      if (ispinor*ncnt.ne.np) then
       write(0,*) '*** error - computed ispinor*ncnt=',ispinor*ncnt,
     &            ' != input nplane=',np;stop
      endif
      return
      end subroutine plindx

!!$*  subroutine for finding k-loop set for certain "k"-point
      subroutine klpfind(wkk, isp,ik,nk,nband,nkx,nky,wnklist,ihf,iz,
     &                        iz2)
      implicit real*8(a-h,o-z)
      dimension wkk(3,5),wk(3),dk1(3),dk2(3)
      real*8 occ(nband),wnklist(3,nk*iz2)
      real*16 ener(nband)
      dk1=0.;dk2=0. !;ne=0 
      dk1(1)=1./real(nkx,8)  ; dk2(2)=1./real(nky,8)  !vector for k-grid
! asign k-loop set wklp for each k-point
      if(iz .eq. 0)then
       if(ihf .eq. 0)then
        call kread(wk,nplane,ener,occ, isp,ik,nk,nband)
       elseif(ihf .eq. 1)then
        wk(:)=wnklist(:,ik)
       endif
       wkk(:,1)=wk(:)
       wkk(:,2)=wk(:)+dk1(:)
       wkk(:,3)=wk(:)+dk1(:)+dk2(:)
       wkk(:,4)=wk(:)       +dk2(:)
       wkk(:,5)=wk(:)
      elseif(iz .eq. 1)then
       if(ihf .eq. 0)then
        call kread(wk,nplane,ener,occ, isp,ik,nk,nband)
       elseif(ihf .eq. 1)then
        wk(:)=wnklist(:,ik)
       endif
       wkk(:,1)=wk(:)
       wkk(:,2)=wk(:)-dk1(:)
       wkk(:,3)=wk(:)-dk1(:)-dk2(:)
       wkk(:,4)=wk(:)       -dk2(:)
       wkk(:,5)=wk(:)
      endif
  
!     do i=1,5
!     write(6,'(i4,3F10.4)')ik, (wkk(jj,i),jj=1,3)
!     enddo
      return
      end subroutine klpfind

!!$*  subroutine for finding k-point index for certain wk
      subroutine kindxfind(ikk,isgg,npl,itr,itrim, wklist,nplist,wkk,
     &                     nk,nband,iz,w_half_klist,
     &                     i_half_klist,i_trim_klist,nhf,iz2)

      implicit real*8(a-h,o-z)
      dimension wkk(3,5),ikk(5),isgg(2,5),wklist(3,nk*iz2),itr(5)
      dimension npl(5),nplist(nk),itrim(5)
      real*8    w_half_klist(3,nk)
      integer   i_half_klist(nk),i_trim_klist(nk)
      del=1E-6 ! criterion for k-point find

      ikk=0;isgg=0;itr=0;itrim=0
!     ! find k-point index ikk for wkk set
      if(iz .eq. 0)then
       do ilp=1,5   !ilp
        do ik=1,nk  !ik
         d1=wklist(1,ik)-wkk(1,ilp)
         d2=wklist(2,ik)-wkk(2,ilp)
         d3=wklist(3,ik)-wkk(3,ilp)
         distk=dsqrt(d1**2 + d2**2 + d3**2)
         if (distk .lt. 1E-5)then
          ikk(ilp)=ik
          npl(ilp)=nplist(ik)
          else if(abs(1.-distk).lt.1E-5)then
           if (abs(1.-dsqrt(d1**2)).lt.1E-5)then
            ikk(ilp)=ik
            npl(ilp)=nplist(ik)
            isgg(1,ilp)=1
           else if(abs(1.-dsqrt(d2**2)).lt.1E-5)then
            ikk(ilp)=ik
            npl(ilp)=nplist(ik)
            isgg(2,ilp)=1
           endif
          else if(abs(sqrt(2.)-distk).lt.1E-5)then
           ikk(ilp)=ik
           isgg(1,ilp)=1
           isgg(2,ilp)=1
           npl(ilp)=nplist(ik)
         endif
        enddo !ik
       enddo  !ilp

      elseif(iz .eq. 1)then
       do ilp=1,5 !ilp
        do ik=1,nhf !ik
         d1=w_half_klist(1,ik)-wkk(1,ilp)
         d2=w_half_klist(2,ik)-wkk(2,ilp)
         d3=w_half_klist(3,ik)-wkk(3,ilp)
         distk=dsqrt(d1**2 + d2**2 + d3**2)      !B+ region including TRIM and edge (all nhf)

         d1t=-w_half_klist(1,ik)-wkk(1,ilp)
         d2t=-w_half_klist(2,ik)-wkk(2,ilp)
         d3t=-w_half_klist(3,ik)-wkk(3,ilp)
         distkt=dsqrt(d1t**2 + d2t**2 + d3t**2)  !time-reversal

         d1tx=-w_half_klist(1,ik)+1.-wkk(1,ilp)
         d2tx=-w_half_klist(2,ik)-wkk(2,ilp)
         d3tx=-w_half_klist(3,ik)-wkk(3,ilp)
         distktx=dsqrt(d1tx**2 + d2tx**2 + d3tx**2) !time-reversal+G1(B- right edge or B- right outside)

         d1ty=-w_half_klist(1,ik)-wkk(1,ilp)
         d2ty=-w_half_klist(2,ik)+1.-wkk(2,ilp)
         d3ty=-w_half_klist(3,ik)-wkk(3,ilp)
         distkty=dsqrt(d1ty**2 + d2ty**2 + d3ty**2) !time-reversal+G2(B+ upper edge or B+ upper ouside)

         d1txy=-w_half_klist(1,ik)+1.-wkk(1,ilp)
         d2txy=-w_half_klist(2,ik)+1.-wkk(2,ilp)
         d3txy=-w_half_klist(3,ik)-wkk(3,ilp)
         distktxy=dsqrt(d1txy**2 + d2txy**2 + d3txy**2) !time-reversal+G1+G2 (B+ upper left outside corner)

         d1x=w_half_klist(1,ik)+1.-wkk(1,ilp)
         d2x=w_half_klist(2,ik)-wkk(2,ilp)
         d3x=w_half_klist(3,ik)-wkk(3,ilp)
         distkx=dsqrt(d1x**2 + d2x**2 + d3x**2) !G1 (B+ right edge)

         d1x_=w_half_klist(1,ik)-1.-wkk(1,ilp)
         d2x_=w_half_klist(2,ik)-wkk(2,ilp)
         d3x_=w_half_klist(3,ik)-wkk(3,ilp)
         distkx_=dsqrt(d1x_**2 + d2x_**2 + d3x_**2) !G1 (B+ left edge)

         d1y=w_half_klist(1,ik)-wkk(1,ilp)
         d2y=w_half_klist(2,ik)-1.-wkk(2,ilp)
         d3y=w_half_klist(3,ik)-wkk(3,ilp)
         distky_=dsqrt(d1y**2 + d2y**2 + d3y**2) !-G2 (B- bottom edge)

         d1m1=-0.5-wkk(1,ilp)
         d2m1=0.0-wkk(2,ilp)
         d3m1=0.0-wkk(3,ilp)
         distkm1=dsqrt(d1m1**2 + d2m1**2 + d3m1**2) !-M1

         d1m2=0.0-wkk(1,ilp)
         d2m2=-0.5-wkk(2,ilp)
         d3m2=0.0-wkk(3,ilp)
         distkm2=dsqrt(d1m2**2 + d2m2**2 + d3m2**2) !-M2

         d1m3=-0.5-wkk(1,ilp)
         d2m3=-0.5-wkk(2,ilp)
         d3m3=0.0-wkk(3,ilp)
         distkm3=dsqrt(d1m3**2 + d2m3**2 + d3m3**2) !-M3

         d1m31= 0.5-wkk(1,ilp)
         d2m31=-0.5-wkk(2,ilp)
         d3m31=0.0-wkk(3,ilp)
         distkm31=dsqrt(d1m31**2 + d2m31**2 + d3m31**2) !-M3+G1

         d1m32=-0.5-wkk(1,ilp)
         d2m32= 0.5-wkk(2,ilp)
         d3m32=0.0-wkk(3,ilp)
         distkm32=dsqrt(d1m32**2 + d2m32**2 + d3m32**2) !-M3+G2

         if (distk .lt. 1E-5)then
          ikk(ilp)=i_half_klist(ik) !Original 
          npl(ilp)=nplist(i_half_klist(ik))
          if(i_trim_klist(ik) .ge. 1)then 
            itrim(ilp)=i_trim_klist(ik)
          endif
         elseif (distkt .lt. 1E-5 .and. distk .gt. 1E-5)then
          if(distkm1 .lt. 1E-5)then
           ikk(ilp)=i_half_klist(ik) ! -M1
           npl(ilp)=nplist(i_half_klist(ik))
           itrim(ilp)=i_trim_klist(ik)
           isgg(1,ilp)=-1
          elseif(distkm2 .lt. 1E-5)then
           ikk(ilp)=i_half_klist(ik) ! -M2
           npl(ilp)=nplist(i_half_klist(ik))
           itrim(ilp)=i_trim_klist(ik)
           isgg(2,ilp)=-1
          elseif(distkm3 .lt. 1E-5)then
           ikk(ilp)=i_half_klist(ik) ! -M3
           npl(ilp)=nplist(i_half_klist(ik))
           itrim(ilp)=i_trim_klist(ik)
           isgg(:,ilp)=-1
          else
           ikk(ilp)=i_half_klist(ik) ! TRS
           npl(ilp)=nplist(i_half_klist(ik))
           itr(ilp)=1
          endif
         elseif (distktx .lt. 1E-5 .and. distk .gt. 1E-5)then
          if(distkm31 .lt. 1E-5)then
           ikk(ilp)=i_half_klist(ik) ! M3-G2
           npl(ilp)=nplist(i_half_klist(ik))
           itrim(ilp)=i_trim_klist(ik)
           isgg(2,ilp)=-1
          else
           ikk(ilp)=i_half_klist(ik) !TRS + G1
           npl(ilp)=nplist(i_half_klist(ik))
           itr(ilp)=1
           isgg(1,ilp)=1
          endif
         elseif (distkty .lt. 1E-5 .and. distk .gt. 1E-5)then
          if(distkm32 .lt. 1E-5)then
           ikk(ilp)=i_half_klist(ik) ! M3-G1
           npl(ilp)=nplist(i_half_klist(ik))
           itrim(ilp)=i_trim_klist(ik)
           isgg(1,ilp)=-1
          else
           ikk(ilp)=i_half_klist(ik) ! TRS + G2
           npl(ilp)=nplist(i_half_klist(ik))
           itr(ilp)=1
           isgg(2,ilp)=1
          endif
         elseif (distkx_ .lt. 1E-5 .and. distk .gt. 1E-5
     &                   .and. wkk(2,ilp) .gt. del
     &                   .and. abs(0.5-wkk(2,ilp)) .gt. del )then
          ikk(ilp)=i_half_klist(ik) ! B+ left edge
          npl(ilp)=nplist(i_half_klist(ik))
          isgg(1,ilp)=-1
         elseif (distky_ .lt. 1E-5 .and. distk .gt. 1E-5
     &                             .and. wkk(1,ilp) .gt.  -0.5+del
     &                             .and. wkk(1,ilp) .lt. -del)then
          ikk(ilp)=i_half_klist(ik) ! B- bottom edge
          npl(ilp)=nplist(i_half_klist(ik))
          isgg(2,ilp)=-1
         endif

        enddo !ik
       enddo !ilp
      endif

      return
      end subroutine kindxfind

!!$*  subroutine for finding k-point vector in the given k index
      subroutine kread(wk,nplane,ener,occ, isp,k,nk,nband)
      implicit real*8(a-h,o-z)
      dimension wk(3)
      complex*16 cener(nband)
      real*8 occ(nband)
      real*16 ener(nband)
      irec=3+(k-1)*(nband+1)+nk*(nband+1)*(isp-1)  !record addres for "k"-point
      read(10,rec=irec) xnplane,(wk(i),i=1,3),
     &(cener(nn),occ(nn),nn=1,nband)
      nplane=nint(xnplane);ener=real(cener)
      if (kpoint.gt.nk) then
       write(0,*) '*** error - selected k=',k,' > max k=',nk;stop
      endif

      return
      end subroutine kread

!This function is to calculate determinant of the complex matrix
!The source is adoped from : https://dualm.wordpress.com/2012/01/06/computing-determinant-in-fortran/
      subroutine get_det(determinant, mat, N)
         implicit none
         integer*4, intent(in) :: N
         complex*16  mat(N,N)
         integer*4   i, info
         integer*4   ipiv(N)
         complex*16  determinant
         real*8      sgn
        
         ipiv = 0
         call zgetrf(N, N, mat, N, ipiv, info)
        
         determinant = (1d0,0d0)
         do i = 1, N
             determinant = determinant*mat(i, i)
         end do
        
         sgn = 1d0
         do i = 1, N
             if(ipiv(i) /= i) then
                 sgn = -sgn
             end if
         end do
         determinant = sgn*determinant

      endsubroutine

!$$*  subroutine to get det(A) where A is n x n complex matrix
!   This subroutine is adopted from below,
!   "http://computer-programming-forum.com/49-fortran/9e079d718158c944.htm"
!   The author is : Che-Ping Su
!c****************************************************
!    Description:     Calculate the determinant of
!                     general complex matrix
!    Input:           A - matrix to be calculated
!                     N - the order of matrix A
!    Output:          DT - Determinant of matrix A
!    Last Updated:      03/17/1998
!**************************************************** 
      subroutine getdetA(DT, S,N)
      implicit real*8(a-h,o-z)
      complex*16  S(N,N),DT,TM,TC,W(N,1)
      complex*16  A(N,N)
      A=(0.,0.)
      A=S
      L=1
      K=2
      flag=(0.,0.)
      DT=cmplx(1.)
   10 TM=A(L,L)
      if(A(L,L) .eq. (0.0,0.0)) then
       do I= L+1,N
        if(A(L,I) .ne. (0.0,0.0)) then
         do J=1,N
          w(J,1)=A(J,L) 
          A(J,L)=A(J,I)
          A(J,I)=W(J,1)
         enddo
         flag=flag+(1.,0.)
         goto 10  
        endif
       enddo
      endif

      do 20 J=L,N
       A(L,J)=A(L,J)/TM
   20 continue
      do 30 I=k,N
       TC=A(I,L)
       do 30 J=L,N
        A(I,J)=A(I,J)-A(L,J)*TC
   30 continue 
      L=L+1
      k=k+1
      DT=DT*TM
      if(L-N) 10,40,40
   40 DT=(-1.,0.)**flag*DT*A(N,N)
      return
      end subroutine getdetA
 
!!$   parse command line arguments
      subroutine parse(filename,foname,nkx,nky,ispinor,icd,ixt,fbz,
     &    ivel,iz,ihf,nini,nmax,kperiod,it,iskp,ine,ver_tag,
     &    iwf,ikwf,ng,rs,imag)
      implicit real*8(a-h,o-z)
      character*75 filename,foname,fbz,ver_tag,vdirec
      real*8 x,y
      character*20 option,value
      integer iarg,narg,ia,nkx,nky,ispinor,iskp,ine,ng(3)
      dimension rs(3)

      nini=1;it=0;iskp=0;ine=0;icd=0;ixt=0;ivel=0;iz=0;ihf=0
      iwf=0;ikwf=1;ng=0;imag=0;rs=0.
      nmax=999999
      iarg=iargc()
      nargs=iarg/2
      filename="WAVECAR"
      foname="BERRYCURV"
      fbz="BERRYCURV.tot.dat"
      if(iarg.ne.2*nargs) then
         call help(ver_tag)
      endif
      do ia=1,nargs
         call getarg(2*ia-1,option)
         call getarg(2*ia,value)
         if(option == "-f") then
            read(value,*) filename
           else if(option == "-o") then
            read(value,*) foname
           else if(option == "-kx") then
            read(value,*) nkx
           else if(option == "-ky") then
            read(value,*) nky
           else if(option == "-s") then
            read(value,*) ispinor
           else if(option == "-ii") then
            read(value,*) nini
           else if(option == "-if") then
            read(value,*) nmax
           else if(option == "-is") then
            read(value,*) nini
            nini=nini;nmax=nini
           else if(option == "-kp") then
            read(value,*) kperiod
           else if(option == "-t") then
            read(value,*) it
           else if(option == "-skp") then
            read(value,*) iskp
           else if(option == "-ne") then
            read(value,*) ine
           else if(option == "-ixt") then
            read(value,*) ixt
           else if(option == "-fbz") then
            read(value,*) fbz
           else if(option == "-cd") then
            read(value,*) icd
           else if(option == "-vel") then
            read(value,*) ivel
           else if(option == "-z2") then
            read(value,*) iz
           else if(option == "-hf") then
            read(value,*) ihf
           else if(option == "-wf") then
            read(value,*) iwf
           else if(option == "-k") then
            read(value,*) ikwf
           else if(option == "-ng") then
            read(value(:),*)ng(1:3)
           else if(option == "-ishift") then
            read(value(:),*)rs(1:3)
           else if(option == "-im") then
            read(value,*) imag
           else if(option =="-h") then
              call help(ver_tag)
           else
           call help(ver_tag)
          endif
      enddo
      if(icd.eq.1 .and. TRIM(foname) .ne. 'BERRYCURV' )then
       write(foname,'(A,A)')"CIRC_DICHROISM.",TRIM(foname)
      else if (icd .eq. 1 .and. TRIM(foname) .eq. 'BERRYCURV') then
        foname="CIRC_DICHROISM"
      else if (icd+ivel .eq. 0 .and. TRIM(foname) .ne. 'BERRYCURV')then
        write(foname,'(A,A)')"BERRYCURV.",TRIM(foname)
      else if (ivel .eq. 1 .and. TRIM(foname) .ne. 'BERRYCURV') then
        write(foname,'(A,A)')"VEL_EXPT.",TRIM(foname)
      else if (ivel .eq. 1 .and. TRIM(foname) .eq. 'BERRYCURV') then
        foname="VEL_EXPT"
      else if (iz .eq. 1 .and. TRIM(foname) .eq. 'BERRYCURV') then
        foname="NFIELD"
      endif

      if(iwf .ge. 1)then
       ine=iwf
       if(iwf .lt. 10) then
        if(ikwf .lt. 10) then
         write(foname,'(A,I1,A,I1)')"PARCHG-W-K00",ikwf,"-E00",iwf
        elseif(ikwf .ge. 10 .and. ikwf .lt. 100)then
         write(foname,'(A,I2,A,I1)')"PARCHG-W-K0",ikwf,"-E00",iwf
        elseif(ikwf .ge. 100 ) then
         write(foname,'(A,I3,A,I1)')"PARCHG-W-K",ikwf,"-E00",iwf
        endif
       elseif(iwf .ge. 10 .and. iwf .lt. 100)then
        if(ikwf .lt. 10) then
         write(foname,'(A,I1,A,I2)')"PARCHG-W-K00",ikwf,"-E0",iwf
        elseif(ikwf .ge. 10 .and. ikwf .lt. 100)then
         write(foname,'(A,I2,A,I2)')"PARCHG-W-K0",ikwf,"-E0",iwf
        elseif(ikwf .ge. 100 ) then
         write(foname,'(A,I3,A,I2)')"PARCHG-W-K",ikwf,"-E0",iwf
        endif
       elseif(iwf .ge. 100 ) then
        if(ikwf .lt. 10) then
         write(foname,'(A,I1,A,I3)')"PARCHG-W-K00",ikwf,"-E",iwf
        elseif(ikwf .ge. 10 .and. ikwf .lt. 100)then
         write(foname,'(A,I2,A,I3)')"PARCHG-W-K0",ikwf,"-E",iwf
        elseif(ikwf .ge. 100 ) then
         write(foname,'(A,I3,A,I3)')"PARCHG-W-K",ikwf,"-E",iwf
        endif
       endif
      endif

      return
      end subroutine parse

!!$*  subroutine for reading basic information
      subroutine inforead(irecl,ispin,nk,nband,ecut,a1,a2,a3,filename)
      implicit real*8(a-h,o-z)
      character*75 filename
      dimension a1(3),a2(3),a3(3)

      irecl=24
      open(unit=10,file=filename,access='direct',recl=irecl,
     & iostat=iost,status='old')
      if (iost .ne. 0) write(6,*) '0.open error - iostat =',iost

      read(10,rec=1)xirecl,xispin,xiprec !RDUM,RISPIN,RTAG(in real type)
      close(10)
      irecl=nint(xirecl);ispin=nint(xispin);iprec=nint(xiprec) ! set to integer
      if(iprec.eq.45210) then
       write(0,*) '*** error - WAVECAR_double requires complex*16';stop
      endif
      open(unit=10,file=filename,access='direct',recl=irecl,
     &iostat=iost,status='old')
      if (iost.ne.0) write(6,*) '1.open error - iostat =',iost
      read(10,rec=2) xnk,xnband,ecut,                 !RNKPTS,RNB_TOT,ENCUT
     &(a1(j),j=1,3),(a2(j),j=1,3),(a3(j),j=1,3)       !A1(3),A2(3),A3(3)
      nk=nint(xnk)
      nband=nint(xnband)

      return
      end subroutine inforead

      subroutine creditinfo(ver_tag)
      character*75 ver_tag

      write(6,*)ver_tag
      write(6,*)"#This program calculates (1)berry curvature omega(k) "
      write(6,*)"#for closed loop C on a small patches in k-space,"
      write(6,*)"#and (2) degree of optical selectivity between "
      write(6,*)"#two bands specified. "
      write(6,*)" "
      write(6,*)"#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(6,*)"#! Copyright 2015. Hyun-Jung Kim All rights reserved.!"
      write(6,*)"#!           (angpangmokjang@hanmail.net)            !"
      write(6,*)"#!             Hanyang Univ. 2015.Apr.05.            !"
      write(6,*)"#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(6,*)" "
      write(6,*)"#*Ref1: T. Fukui, Y. Hatsugai, and H. Suzuki, "
      write(6,*)"#       J. J. Phys. Soc. Jap. 74, 1674 (2005),"
      write(6,*)"#     'Chern Numbers in Discretized Brillouin Zone: "
      write(6,*)"#  Efficient Method of Computing (Spin) Hall 
     &Conductances'"
      write(6,*)"#*Ref2: R. Resta, J. Phys.: Condens. Matter. 12, R107 
     &(2000)"
      write(6,*)"#     'Menifestations of Berry's phase in molecules "
      write(6,*)"#      and condensed matter'"
      write(6,*)"#*Ref3: W. Yao, D. Xiao, and Q. Niu, PRB 77, 235406 
     &(2008)"
      write(6,*)"#     'Valley-dependent optoelectronics from inversion"
      write(6,*)"#      symmetry breaking'"
      write(6,*)"#*Ref4: http://www.andrew.cmu.edu/user/feenstra/
     &wavetrans/"
      write(6,*)"#       R. M. Feenstra and M. Widom                  "
      write(6,*)"#       some routines has been adopted from WAVETRANS"
      write(6,*)"#*      routines for Wavefunction reading and G "
      write(6,*)"#       matching "
      write(6,*)" "
      write(6,*)"#*Syntax:"
      write(6,*)"#       berry -f file -kx nx -ky ny -s 2 -ii ni -if nf"
      write(6,*)"#*  (or)berry -f file -kx nx -ky ny -s 2 -is n "
      write(6,*)" "
      write(6,*)"#*For the detailed help: ./berry -h"
      write(6,*)" "
      write(6,*)" "
      end subroutine creditinfo

      subroutine write_special_kpoint(b1,b2,b3)
      implicit real*8(a-h,o-z)
      dimension SKP(3),SP(3),b1(3),b2(3),b3(3)
      open(41,file='SKP.dat',status='unknown')
      do i=1,3;SKP(i)=( 0.                     * b1(i))+
     &                ( 0.                     * b2(i))+ 
     &                ( 0.                     * b3(i));enddo
      write(41,'(3F12.6,A)')(SKP(i),i=1,3),"  # G"

      do i=1,3;SKP(i)=( 2. / 3.                * b1(i))+
     &                ( 1. / 3.                * b2(i))+ 
     &                ( 0.                     * b3(i));enddo
      write(41,'(3F12.6,A)')(SKP(i),i=1,3),"  # K1"
      do i=1,3;SKP(i)=(-1. / 3.                * b1(i))+
     &                ( 1. / 3.                * b2(i))+ 
     &                ( 0.                     * b3(i));enddo
      write(41,'(3F12.6,A)')(SKP(i),i=1,3),"  # K1"

      do i=1,3;SKP(i)=( 1. / 3.                * b1(i)) +
     &                ( 2. / 3.                * b2(i)) + 
     &                ( 0.                     * b3(i));enddo
      write(41,'(3F12.6,A)')(SKP(i),i=1,3),"  # K2"
      do i=1,3;SKP(i)=( 1. / 3.                * b1(i))+
     &                (-1. / 3.                * b2(i))+
     &                ( 0.                     * b3(i));enddo
      write(41,'(3F12.6,A)')(SKP(i),i=1,3),"  # K2"

      do i=1,3;SKP(i)=( 1. / 2.                * b1(i))+
     &                ( 0.                     * b2(i))+ 
     &                ( 0.                     * b3(i));enddo
      write(41,'(3F12.6,A)')(SKP(i),i=1,3),"  # M1"

      do i=1,3;SKP(i)=( 1. / 2.                * b1(i))+
     &                ( 1. / 2.                * b2(i))+ 
     &                ( 0.                     * b3(i));enddo
      write(41,'(3F12.6,A)')(SKP(i),i=1,3),"  # M2"

      do i=1,3;SKP(i)=( 0.                     * b1(i))+
     &                ( 1. / 2.                * b2(i))+ 
     &                ( 0.                     * b3(i));enddo
      write(41,'(3F12.6,A)')(SKP(i),i=1,3),"  # M3"

      close(41)
      end subroutine write_special_kpoint

      subroutine extendingBZ(fbz,ixt)
      implicit real*8 (a-h, o-z)
      character*75  fbz
      character*200  A,S,P
      dimension b1(3),b2(3),b3(3),xb(3)
      real*8, allocatable :: xrecivec(:,:),xrecilat(:,:),xdata_(:)
      real*8, allocatable ::recivec(:,:),recilat(:,:),data_(:)
      integer IOstatus    

      IOstatus=0;II=1;ik=0
      open(51,file=fbz,status='unknown')
      open(61,file='EXT.dat',status='unknown')
      read(51,'(A)',IOSTAT=IOstatus)A
      S=A(1:1)
      do while (TRIM(S) .eq. '#' .or. II .eq. 35)
       write(6,'(A)')TRIM(A)
       write(61,'(A)')TRIM(A)
       read(51,'(A)',IOSTAT=IOstatus)A
       if (TRIM(A(3:9)) .eq. 'NKPOINT') then
        P=A(23:27);read(P,*)nk
        else if(TRIM(A(3:12)) .eq. 'RECIVEC B1') then
         P=A(26:34);read(P,*)b1(1)
         P=A(39:47);read(P,*)b1(2)
         P=A(52:60);read(P,*)b1(3)
        else if (TRIM(A(3:12)) .eq. 'RECIVEC B2') then
         P=A(26:34);read(P,*)b2(1)
         P=A(39:47);read(P,*)b2(2)
         P=A(52:60);read(P,*)b2(3)
        else if (TRIM(A(3:12)) .eq. 'RECIVEC B3') then
         P=A(26:34);read(P,*)b3(1)
         P=A(39:47);read(P,*)b3(2)
         P=A(52:60);read(P,*)b3(3)
       endif
       S=TRIM(A(1:1))
       II=ICHAR(TRIM(A))
      enddo
      allocate(recivec(3,nk))
      allocate(recilat(3,nk))
      allocate(data_(nk))
      allocate(xrecivec(3,ixt*2*ixt*2*nk))
      allocate(xrecilat(3,ixt*2*ixt*2*nk))
      allocate(xdata_(ixt*2*ixt*2*nk))

      backspace(51)
      do while(IOstatus .eq. 0)
       ik=ik+1
       read(51,*,IOSTAT=IOstatus)(recivec(i,ik),i=1,3),data_(ik),
     &                           (recilat(i,ik),i=1,3)
      enddo

      do ik=1,nk
       write(6,'(3F11.6,A,F16.6,A,3F11.6)')(recivec(i,ik),i=1,3),
     &      "     ",data_(ik),
     &      "                  ",(recilat(i,ik),i=1,3)
      enddo

      kk=0   ! extend 
      do ib2=-1*(ixt-1)+1,ixt
       do ib1=-1*(ixt-1)+1,ixt 
        do ik=1,nk
         kk=kk+1
         xrecivec(1,kk)=recivec(1,ik) +(ib1-1)*b1(1)+(ib2-1)*b2(1)
         xrecivec(2,kk)=recivec(2,ik) +(ib1-1)*b1(2)+(ib2-1)*b2(2)
         xrecivec(3,kk)=recivec(3,ik) +(ib1-1)*b1(3)+(ib2-1)*b2(3)
         xrecilat(1,kk)=recilat(1,ik) +(ib1-1)
         xrecilat(2,kk)=recilat(2,ik) +(ib2-1)
         xrecilat(3,kk)=recilat(3,ik)
         xdata_(kk)=data_(ik)
         kext=kk
        enddo
       enddo
      enddo

      write(6,*)" "
      write(6,'(A,I1,A,I1)')"# EXTENDING DATA GRID by : ",ixt,' x ',ixt
      write(6,'(A)')"# SORTING K-grids..."
      do k=ixt-1,1,-1  ! sorting kx
       do j=1,k
        if(xrecivec(1,j+1) .gt. xrecivec(1,j))then
         xb(:)=xrecivec(:,j)
         xrecivec(:,j)=xrecivec(:,j+1)
         xrecivec(:,j+1)=xb(:)
         xb(:)=xrecilat(:,j)
         xrecilat(:,j)=xrecilat(:,j+1)
         xrecilat(:,j+1)=xb(:)
         xtemp=xdata_(j)
         xdata_(j)=xdata_(j+1)
         xdata_(j+1)=xtemp
        endif
       enddo
      enddo
      do k=ixt-1,1,-1  ! sorting ky
       do j=1,k
        if(xrecivec(1,j+1) .eq. xrecivec(1,j))then
         if(xrecivec(2,j+1) .gt. xrecivec(2,j))then
         xb(:)=xrecivec(:,j)
         xrecivec(:,j)=xrecivec(:,j+1)
         xrecivec(:,j+1)=xb(:)
         xb(:)=xrecilat(:,j)
         xrecilat(:,j)=xrecilat(:,j+1)
         xrecilat(:,j+1)=xb(:)
         xtemp=xdata_(j)
         xdata_(j)=xdata_(j+1)
         xdata_(j+1)=xtemp
         endif
        endif
       enddo
      enddo

      do ik=1,kext
       write(61,'(3F11.6,A,F16.6,A,3F11.6)')(xrecivec(i,ik),i=1,3),
     &      "     ",xdata_(ik),"                  ",
     &      (xrecilat(i,ik),i=1,3)
      enddo
      write(6,'(A)')"# DONE! result is in 'EXT.dat' "

      stop
      end subroutine extendingBZ

      subroutine test
      implicit real*8 (a-h, o-z)
      complex*8   a,b,c,d
      character*75  foname
       a=(2.2,-1.3)
       b=(3.2,-5.4)
      
      c=a+(0.,1.)*b       
      d=a-(0.,1.)*b       

      write(6,*)"AAA",a
      write(6,*)"BBB",b
      write(6,*)"CCC",c,abs(c)
      write(6,*)"DDD",d,abs(d)

      stop
      end subroutine test

      subroutine help(ver_tag)
      character*75 ver_tag
      write(6,*)"          **** PROGRAM INSTRUCTION ***"
      write(6,*)" "
      write(6,*)ver_tag
      write(6,*)" "
      write(6,*)"*LIMITATION : -This program is ONLY for 2D system."
      write(6,*)"            : -It tries to find the VBM for first KPT,"
      write(6,*)"            :  it may result in some problem when you"
      write(6,*)"            :  dealing with metallic states. "
      write(6,*)"            :  Be careful!"
      write(6,*)"*NOTE1 The x,y,z components of each G value are given"
      write(6,*)"       in terms of the ig values and the components "
      write(6,*)"       of the recip. lattice vectors according to:"
      write(6,*)"        ig1*b1_x + ig2*b2_x + ig3*b3_x,"
      write(6,*)"        ig1*b1_y + ig2*b2_y + ig3*b3_y, and"
      write(6,*)"        ig1*b1_z + ig2*b2_z + ig3*b3_z, respectively,"
      write(6,*)"       with"
      write(6,*)"        ig1=ig(1,iplane),"
      write(6,*)"        ig2=ig(2,iplane), and"
      write(6,*)"        ig3=ig(3,iplane),"
      write(6,*)"       where iplane=1,2,...,nplane(k,ispin) is an "
      write(6,*)"       index incrementing the plane waves for specific"
      write(6,*)"       k and spin values"
      write(6,*)" "
      write(6,*)"*NOTE2 The energy eigenvalues are complex, as provided"
      write(6,*)"       in the WAVECAR file, but the imaginary part is "
      write(6,*)"       zero (at least for cases investigated thus far)"
      write(6,*)" "
      write(6,*)"*Syntax:berry -f file -kx nx -ky ny -s 2 -ii ni -if nf"
      write(6,*)"*   (or)berry -f file -kx nx -ky ny -s 2 -is n "
      write(6,*)" "
      write(6,*)"             ### POSSIBLE OPTIONS ###"
      write(6,*)" -f filename      : File name to be read"
      write(6,*)"                  : Default: WAVECAR"
      write(6,*)" -kx(ky) kx(ky)   : k-point grid of your system"
      write(6,*)" -s 2 or 1        : for the noncollinear case, -s 2"
      write(6,*)"                  : for the collinear or NM,   -s 1"
      write(6,*)"                  :  Default : 2 if ISPIN 1"
      write(6,*)"                  :          : 1 if ISPIN 2"
      write(6,*)" -ii(if) ni(nf)   : Once specified, berry curvature "
      write(6,*)"                  : for multiband (from ni-th to nf-th"
      write(6,*)"                  : state) states will be evaluated,"
      write(6,*)"                  : unless '-cd' is not 1."
      write(6,*)"                  :  Default -cd 0 -> -ii  1  -if VBM "
      write(6,*)"                  :          -cd 1 -> -ii VBM -if CBM"
      write(6,*)" -is n            : if specified, berry curvature for"
      write(6,*)"                  : single band (n-th band) will be "
      write(6,*)"                  : evaluated." 
      write(6,*)" -o filename      : Print Berry curvature distribution"
      write(6,*)"                  : to the 'filename.dat'             "
      write(6,*)"                  :  Default : BERRYCURV.dat          "
      write(6,*)" -kp np           : Print Berry curvature distribution"
      write(6,*)"                  : with extended BZ with 'np x np'-BZ"
      write(6,*)"                  :  Default : 2              "
      write(6,*)" -skp 1(or 0)     : Specify whether special K-points "
      write(6,*)"                  : will be printed in the SKP.dat "
      write(6,*)"                  : You may specify the points by hand"
      write(6,*)"                  : into the source code, the routine,"
      write(6,*)"                  : 'write_special_kpoint'"
      write(6,*)"                  :  Default : 0"
      write(6,*)" -ne  ne          : Specify total number of electrons"
      write(6,*)"                  : Usual insulating cases, you don't"
      write(6,*)"                  : need to specify it, program will "
      write(6,*)"                  : find total number of electrons,"
      write(6,*)"                  : but for semimetallic system, it "
      write(6,*)"                  : may be more clear to put by hand"
      write(6,*)" -cd  1(or 0)     : Calculate spin- and k-resolved"
      write(6,*)"                  : degree of optical polarization,"
      write(6,*)"                  : degree of circular polarization,"
      write(6,*)"                  : between valence & conduction band"
      write(6,*)"                  : If set to '1', Berry cuvature will"
      write(6,*)"                  : not be evaluated"
      write(6,*)"                  : You may set -ii and -if together,"
      write(6,*)"                  : defining VBM & CBM, respectively"
      write(6,*)"                  :  Default : 0"
      write(6,*)" -ixt np -fbz f   : **For the special purpose"
      write(6,*)"                  : Read file 'f' and extend data to"
      write(6,*)"                  : np x np periodic field. Note that"
      write(6,*)"                  : -ixt option should be used with "
      write(6,*)"                  :'-fbz filename_of_.dat_file' option"
      write(6,*)" -vel 1 -is n     : **For the special purpose"
      write(6,*)"                  : Calculate velocity expectation "
      write(6,*)"                  : vaule (v_x, v_y) of n-th state"
      write(6,*)" -z2  1           : Calculate Z2 invariant   "
      write(6,*)"                  : Using Fukui's method (JSPJ 76,"
      write(6,*)"                  : 053702 2007) which is lattice "
      write(6,*)"                  : version of Fu and Kane method"
      write(6,*)"                  : (PRB 74, 195312 2006) "
      write(6,*)" -hf  1           : Use half BZ, ISYM=1 or 2, "
      write(6,*)"                  : currently only works with -z2 1"
      write(6,*)" -wf nb -k nk     : **For the special purpose"
      write(6,*)"  -ng nx,ny,nz    : Calculate real-space wavefunction"
      write(6,*)"  -im 1           : output will be written in "
      write(6,*)"  -ishift rx,ry,rz: PARCHG-W-K$k-E$wf-(IM)-SPIN$ispin"
      write(6,*)"                  : nb and nk : index of band & kpoint"
      write(6,*)"                  : -ng option determins nx,ny,nz grid"
      write(6,*)"                  :  for the cube: default=nbmax"
      write(6,*)"                  : -im option: plot imaginary?"
      write(6,*)"                  : -ishift: shift origin with respect"
      write(6,*)"                  :  to the direct coord. rx,ry,rz"
      write(6,*)" "
      write(6,*)"* default: -f WAVECAR -kx 2 -ky 2 -s 2 -ii 1 -if VBM 
     &-kp 1"
      write(6,*)"* here, VBM is valence band maximum"
      write(6,*)" "
      write(6,*)"*Compilation:  gfortran or ifort. "
      write(6,*)" Flag '-assume byterecl' is required for ifort."
      write(6,*)" for OSX,  -Wl,-stack_size,0x80000000 may be required"
      write(6,*)" ex-noMPI)ifort -fpp -assume byterecl -mkl 
     &-o vaspberry vaspberry.f"
      write(6,*)" ex-MPI)mpif90 -DMPI_USE -mkl -fpp -assume byterecl 
     &-o vaspberry vaspberry.f"
      write(6,*)" ex-noMPI-gfortran) gfortran -I/opt/local/include 
     &-L/opt/local/lib/lapack/ -l lapack -o vaspberry vaspberry_gfortran
     &_serial.f"

!     write(6,*)" ex-MPI) mpif90 -DMPI_USE -mkl -fpp -assume byterecl -o vaspberry vaspberry.f "

      stop
      end subroutine help
