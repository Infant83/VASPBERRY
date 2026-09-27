! Native projected-spin Chern numbers in the WAVECAR pseudo metric.
! No time-reversal gauge replacement and no energy-denominator Kubo step.
      subroutine spinchern_run(wavepath,outpath,prefix,nx,ny,
     & first,last,ispin,components,nk,nb,ecut,a1,a2,a3,axis,
     & etol,stol,rank,nproc)
      use, intrinsic :: ieee_arithmetic
      implicit none
#ifdef MPI_USE
      include 'mpif.h'
#endif
      character*(*) wavepath,outpath,prefix
      integer nx,ny,first,last,ispin,components,nk,nb,rank,nproc
      integer r,k,i,j,m,n,ix,iy,jx,jy,k2,dir,ng,ng2,ios,rec
      integer npmax,nbmax(3),info,neg,pos,neg0,pos0,ierr,q
      integer, allocatable :: npl(:),grid(:,:),gx(:),gy(:)
      integer, allocatable :: g1(:,:),g2(:,:)
      real*8 ecut,a1(3),a2(3),a3(3),axis(3),etol,stol
      real*8 lat(3,3),b(3,3),cross(3),volume,area,pi
      real*8 rot(3,3),localaxis(3),xnp,offset(2),u,v,delta(3)
      real*8 gap,indirect,phase,smin,maxphase,chern(3),difference
      real*8, parameter :: ktol=1d-7,linktol=1d-6
      real*8, allocatable :: wk(:,:),en(:,:),occ(:),eval(:,:)
      real*8, allocatable :: egap(:),grammin(:),gramcond(:)
      real*8, allocatable :: link(:,:,:),lsv(:,:),flux(:,:)
      real*8, allocatable :: cov(:,:)
      complex*16, allocatable :: ec(:),c1(:,:),c2(:,:),vec(:,:,:)
      complex*16, allocatable :: mat(:,:),sector(:,:)
      logical occupied
      pi=4d0*atan(1d0)
      r=last-first+1
      if(ispin.ne.1.or.components.ne.2)then
       call spinchern_error('requires ISPIN=1 two-component spinors',0)
      endif
      if(first.lt.1.or.last.gt.nb.or.r.lt.1.or.r.eq.nb)then
       call spinchern_error('select a proper group within 1:NBANDS',0)
      endif
      if(nx.lt.2.or.ny.lt.2.or.int(nx,8)*ny.ne.nk)then
       call spinchern_error('requires a complete periodic NX by NY mesh',0)
      endif
! Bound the replicated sector frames and the MPI collective counts.
      if(int(r,8)*r*nk.gt.32000000_8.or.
     &   int(nb,8)*nk.gt.32000000_8)then
       call spinchern_error('selected frame exceeds bounded storage',0)
      endif
      lat(:,1)=a1;lat(:,2)=a2;lat(:,3)=a3
      cross=(/a2(2)*a3(3)-a2(3)*a3(2),
     &        a2(3)*a3(1)-a2(1)*a3(3),
     &        a2(1)*a3(2)-a2(2)*a3(1)/)
      volume=dot_product(a1,cross)
      b(:,1)=2d0*pi*cross/volume
      cross=(/a3(2)*a1(3)-a3(3)*a1(2),
     &        a3(3)*a1(1)-a3(1)*a1(3),
     &        a3(1)*a1(2)-a3(2)*a1(1)/)
      b(:,2)=2d0*pi*cross/volume
      cross=(/a1(2)*a2(3)-a1(3)*a2(2),
     &        a1(3)*a2(1)-a1(1)*a2(3),
     &        a1(1)*a2(2)-a1(2)*a2(1)/)
      b(:,3)=2d0*pi*cross/volume
      cross=(/b(2,1)*b(3,2)-b(3,1)*b(2,2),
     &        b(3,1)*b(1,2)-b(1,1)*b(3,2),
     &        b(1,1)*b(2,2)-b(2,1)*b(1,2)/)
      area=sqrt(sum(cross**2))/dble(nx)/dble(ny)
      call spinchern_frame(outpath,nk,nb,lat,rot)
      localaxis=matmul(transpose(rot),axis)
      localaxis=localaxis/sqrt(sum(localaxis**2))
      call reciproperty(nbmax,npmax,b(:,1),b(:,2),b(:,3),ecut,2)
      allocate(wk(3,nk),en(nb,nk),npl(nk),ec(nb),occ(nb))
      allocate(grid(nx,ny),gx(nk),gy(nk),egap(nk))
      allocate(vec(r,r,nk),eval(r,nk),grammin(nk),gramcond(nk))
      allocate(link(3,2,nk),lsv(2,nk),cov(2,nk),flux(3,nk))
      occupied=first.eq.1
      do k=1,nk
       rec=3+(k-1)*(nb+1)
       read(10,rec=rec,iostat=ios)xnp,wk(:,k),
     &     (ec(i),occ(i),i=1,nb)
       if(ios.ne.0)call spinchern_error('WAVECAR header read failed',k)
       if(.not.all(ieee_is_finite(wk(:,k))).or.
     &    .not.all(ieee_is_finite(real(ec,8))).or.
     &    .not.all(ieee_is_finite(occ)))then
        call spinchern_error('nonfinite WAVECAR header',k)
       endif
       npl(k)=nint(xnp)
       en(:,k)=real(ec,8)
       if(any(en(2:nb,k).lt.en(1:nb-1,k)-etol))then
        call spinchern_error('band energies are not ordered',k)
       endif
       gap=huge(1d0)
       do n=first,last
        do m=1,nb
         if(m.ge.first.and.m.le.last)cycle
         gap=min(gap,abs(en(n,k)-en(m,k)))
        enddo
       enddo
       egap(k)=gap
       if(gap.le.etol)then
        write(0,*)'spin-chern parent gap eV / threshold:',gap,etol
        call spinchern_error('parent subspace is not energy isolated',k)
       endif
       if(any(abs(occ(first:last)-1d0).gt.1d-3))occupied=.false.
       if(first.gt.1)then
        if(any(abs(occ(:first-1)).gt.1d-3))occupied=.false.
       endif
       if(last.lt.nb)then
        if(any(abs(occ(last+1:)).gt.1d-3))occupied=.false.
       endif
      enddo
      indirect=huge(1d0)
      if(first.eq.1.and.last.lt.nb)then
       indirect=minval(en(last+1,:))-maxval(en(last,:))
      endif
! Establish grid index independently of source order and signed k convention.
      offset(1)=minval(modulo(wk(1,:),1d0))
      offset(2)=minval(modulo(wk(2,:),1d0))
      grid=0
      do k=1,nk
       u=modulo(wk(1,k)-offset(1),1d0)*nx
       v=modulo(wk(2,k)-offset(2),1d0)*ny
       if(abs(u-nint(u)).gt.ktol.or.abs(v-nint(v)).gt.ktol.or.
     &    abs(wk(3,k)-wk(3,1)-nint(wk(3,k)-wk(3,1))).gt.ktol)then
        call spinchern_error('nonuniform or nonplanar periodic mesh',k)
       endif
       ix=modulo(nint(u),nx)+1;iy=modulo(nint(v),ny)+1
       if(grid(ix,iy).ne.0)call spinchern_error('duplicate mesh point',k)
       grid(ix,iy)=k;gx(k)=ix;gy(k)=iy
      enddo
      if(any(grid.eq.0))call spinchern_error('missing mesh point',0)
! Every rank obtains identical eigenspaces by gathering disjoint k frames.
      vec=(0d0,0d0);eval=0d0;grammin=0d0;gramcond=0d0
      do k=rank+1,nk,nproc
       allocate(c1(npl(k),r),g1(3,npmax))
       call spinchern_read(k,first,last,nk,nb,npl(k),wk(:,k),
     &      nbmax,npmax,ecut,b,c1,g1,ng)
       call spinchern_split(c1,ng,r,localaxis,vec(:,:,k),
     &      eval(:,k),grammin(k),gramcond(k),k)
       deallocate(c1,g1)
      enddo
#ifdef MPI_USE
      call MPI_ALLREDUCE(MPI_IN_PLACE,vec,r*r*nk,
     & MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,eval,r*nk,
     & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,grammin,nk,
     & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,gramcond,nk,
     & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif
      neg0=count(eval(:,1).lt.-stol);pos0=count(eval(:,1).gt.stol)
      do k=1,nk
       if(minval(abs(eval(:,k))).le.stol)then
        write(0,*)'spin-chern min |Pauli| / threshold:',
     &              minval(abs(eval(:,k))),stol
        call spinchern_error('projected-spin zero gap is closed',k)
       endif
       neg=count(eval(:,k).lt.-stol);pos=count(eval(:,k).gt.stol)
       if(neg.ne.neg0.or.pos.ne.pos0)then
        call spinchern_error('projected-spin sector rank changes',k)
       endif
      enddo
      link=0d0;lsv=0d0;cov=0d0
      do k=rank+1,nk,nproc
       allocate(c1(npl(k),r),g1(3,npmax))
       call spinchern_read(k,first,last,nk,nb,npl(k),wk(:,k),
     &      nbmax,npmax,ecut,b,c1,g1,ng)
       ix=gx(k);iy=gy(k)
       do dir=1,2
        jx=ix;jy=iy
        if(dir.eq.1)jx=modulo(ix,nx)+1
        if(dir.eq.2)jy=modulo(iy,ny)+1
        k2=grid(jx,jy)
        delta=wk(:,k)-wk(:,k2)
        if(dir.eq.1)delta(1)=delta(1)+1d0/nx
        if(dir.eq.2)delta(2)=delta(2)+1d0/ny
        if(maxval(abs(delta-nint(delta))).gt.ktol)then
         call spinchern_error('noninteger periodic sewing shift',k)
        endif
        allocate(c2(npl(k2),r),g2(3,npmax))
        call spinchern_read(k2,first,last,nk,nb,npl(k2),wk(:,k2),
     &       nbmax,npmax,ecut,b,c2,g2,ng2)
        allocate(mat(r,r))
        call spinchern_overlap(c1,c2,g1,g2,ng,ng2,r,nbmax,
     &       nint(delta),mat,cov(dir,k))
        deallocate(c2,g2)
! Full normalized frame and both sectors use the same pseudo link metric.
        mat=matmul(conjg(transpose(vec(:,:,k))),
     &                  matmul(mat,vec(:,:,k2)))
        lsv(dir,k)=huge(1d0)
        do q=1,3
         if(q.eq.1)then
          i=neg0+1;j=r
         else if(q.eq.2)then
          i=1;j=neg0
         else
          i=1;j=r
         endif
         if(j.lt.i)then
          link(q,dir,k)=0d0
          cycle
         endif
         allocate(sector(j-i+1,j-i+1))
         sector=mat(i:j,i:j)
         if(.not.all(ieee_is_finite(real(sector,8))).or.
     &      .not.all(ieee_is_finite(aimag(sector))))then
          call spinchern_error('nonfinite sector link',k)
         endif
         call z2_link_svd_phase(sector,j-i+1,phase,smin,info)
         if(info.ne.0.or.smin.le.linktol)then
          write(0,*)'spin-chern link axis,sector,smin:',dir,q,smin
          call spinchern_error('singular or poorly resolved link',k)
         endif
         link(q,dir,k)=phase
         lsv(dir,k)=min(lsv(dir,k),smin)
         deallocate(sector)
        enddo
        deallocate(mat)
       enddo
       deallocate(c1,g1)
      enddo
#ifdef MPI_USE
      call MPI_ALLREDUCE(MPI_IN_PLACE,link,6*nk,
     & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,lsv,2*nk,
     & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,cov,2*nk,
     & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif
      do k=1,nk
       ix=gx(k);iy=gy(k)
       jx=modulo(ix,nx)+1;jy=modulo(iy,ny)+1
       do q=1,3
        phase=link(q,1,k)+link(q,2,grid(jx,iy))-
     &        link(q,1,grid(ix,jy))-link(q,2,k)
        flux(q,k)=-atan2(sin(phase),cos(phase))
       enddo
      enddo
      maxphase=maxval(abs(flux))
      if(maxphase.ge.pi-1d-5)then
       call spinchern_error('plaquette phase touches branch cut; refine mesh',0)
      endif
      chern=sum(flux,dim=2)/(2d0*pi)
      difference=abs(chern(1)+chern(2)-chern(3))
      if(difference.gt.1d-6)then
       call spinchern_error('sector sum disagrees with parent; refine mesh',0)
      endif
      if(rank.eq.0)then
       call spinchern_write(wavepath,outpath,prefix,nx,ny,first,last,
     & nk,r,wk,axis,rot,etol,stol,egap,indirect,occupied,
     & grammin,gramcond,eval,neg0,pos0,flux,area,lsv,cov,chern)
       write(6,*)'# Spin Chern numbers (pseudo metric): C+, C-, Ccharge, Cspin'
       write(6,'(4ES24.16)')chern(1),chern(2),chern(3),
     &                       (chern(1)-chern(2))/2d0
       write(6,*)'# Spin Chern result_status=PASS'
      endif
      deallocate(wk,en,npl,ec,occ,grid,gx,gy,egap,vec,eval,
     & grammin,gramcond,link,lsv,cov,flux)
      end subroutine spinchern_run

      subroutine spinchern_error(message,k)
      implicit none
      character*(*) message
      integer k
      write(0,*)'*** error - spin-chern: ',trim(message),'; k=',k
      call vaspberry_fail
      end subroutine spinchern_error

! Read only the actual stored eigenvectors; never synthesize Kramers partners.
      subroutine spinchern_read(k,first,last,nk,nb,np,wk,
     & nbmax,npmax,ecut,b,c,g,ng)
      use, intrinsic :: ieee_arithmetic
      implicit none
      integer k,first,last,nk,nb,np,npmax,nbmax(3),ng,ios,n,r
      real*8 wk(3),ecut,b(3,3)
      complex*16 c(np,last-first+1)
      integer g(3,npmax)
      complex*8, allocatable :: raw(:)
      r=last-first+1
      allocate(raw(np))
      call plindx(g,ng,2,wk,b(:,1),b(:,2),b(:,3),
     & nbmax,np,ecut,npmax)
      if(ng*2.ne.np)call spinchern_error('coefficient/G layout mismatch',k)
      do n=first,last
       read(10,rec=3+(k-1)*(nb+1)+n,iostat=ios)raw
       if(ios.ne.0)call spinchern_error('coefficient read failed',k)
       c(:,n-first+1)=cmplx(raw,kind=8)
      enddo
      if(.not.all(ieee_is_finite(real(c,8))).or.
     &   .not.all(ieee_is_finite(aimag(c))))then
       call spinchern_error('nonfinite coefficients',k)
      endif
      deallocate(raw)
      end subroutine spinchern_read

! Generalized Hermitian eigenproblem in the raw pseudo-wavefunction metric.
      subroutine spinchern_split(c,ng,r,axis,v,e,gmin,gcond,k)
      use, intrinsic :: ieee_arithmetic
      implicit none
      integer ng,r,k,lwork,info,i
      real*8 axis(3),e(r),gmin,gcond
      real*8, allocatable :: ge(:),rw(:)
      complex*16 c(2*ng,r),v(r,r)
      complex*16, allocatable :: o(:,:),test(:,:),ud(:,:),s(:,:)
      complex*16, allocatable :: work(:)
      allocate(o(r,r),test(r,r),ud(r,r),s(r,r),ge(r))
      lwork=max(1,4*r);allocate(work(lwork),rw(max(1,3*r)))
      o=matmul(conjg(transpose(c)),c)
      test=o
      call zheev('N','U',r,test,r,ge,work,lwork,rw,info)
      if(info.ne.0.or..not.all(ieee_is_finite(ge)))then
       call spinchern_error('Gram eigensolver failed',k)
      endif
      gmin=minval(ge)
      if(gmin.le.1d-12.or.gmin/maxval(ge).le.1d-10)then
       call spinchern_error('singular or ill-conditioned pseudo Gram metric',k)
      endif
      gcond=maxval(ge)/gmin
      s=axis(3)*(matmul(conjg(transpose(c(:ng,:))),c(:ng,:))-
     &           matmul(conjg(transpose(c(ng+1:,:))),c(ng+1:,:)))
      ud=matmul(conjg(transpose(c(:ng,:))),c(ng+1:,:))
      s=s+cmplx(axis(1),-axis(2),8)*ud+
     &    cmplx(axis(1),axis(2),8)*conjg(transpose(ud))
      call zhegv(1,'V','U',r,s,r,o,r,e,work,lwork,rw,info)
      if(info.ne.0.or..not.all(ieee_is_finite(e)).or.
     & .not.all(ieee_is_finite(real(s,8))).or.
     & .not.all(ieee_is_finite(aimag(s))))then
       call spinchern_error('projected-spin generalized eigensolver failed',k)
      endif
      if(maxval(abs(e)).gt.1d0+1d-8)then
       call spinchern_error('projected Pauli spectrum violates metric bound',k)
      endif
      v=s
      test=matmul(conjg(transpose(v)),
     & matmul(matmul(conjg(transpose(c)),c),v))
      do i=1,r
       test(i,i)=test(i,i)-1d0
      enddo
      if(maxval(abs(test)).gt.1d-8)then
       call spinchern_error('sector Gram orthonormality residual',k)
      endif
      deallocate(o,test,ud,s,ge,work,rw)
      end subroutine spinchern_split

! Exact integer reciprocal sewing; missing plane waves are zero padded.
      subroutine spinchern_overlap(c1,c2,g1,g2,ng1,ng2,r,nbmax,
     & shift,m,coverage)
      implicit none
      integer ng1,ng2,r,nbmax(3),shift(3),g1(3,*),g2(3,*)
      integer i,j,idx,ngrid,gg(3),matched
      integer, allocatable :: lookup(:)
      complex*16 c1(2*ng1,r),c2(2*ng2,r),m(r,r)
      complex*16, allocatable :: aligned(:,:)
      real*8 coverage
      ngrid=(2*nbmax(1)+1)*(2*nbmax(2)+1)*(2*nbmax(3)+1)
      allocate(lookup(ngrid),aligned(2*ng1,r))
      lookup=0;aligned=(0d0,0d0);matched=0
      do i=1,ng1
       gg=g1(:,i)+nbmax
       idx=1+gg(1)+(2*nbmax(1)+1)*(gg(2)+(2*nbmax(2)+1)*gg(3))
       lookup(idx)=i
      enddo
      do j=1,ng2
       gg=g2(:,j)-shift
       if(any(abs(gg).gt.nbmax))cycle
       gg=gg+nbmax
       idx=1+gg(1)+(2*nbmax(1)+1)*(gg(2)+(2*nbmax(2)+1)*gg(3))
       i=lookup(idx)
       if(i.eq.0)cycle
       aligned(i,:)=c2(j,:)
       aligned(i+ng1,:)=c2(j+ng2,:)
       matched=matched+1
      enddo
      m=matmul(conjg(transpose(c1)),aligned)
      coverage=dble(matched)/max(ng1,ng2)
      deallocate(lookup,aligned)
      end subroutine spinchern_overlap

! Require actual OUTCAR spin rotation and matching source dimensions/lattice.
      subroutine spinchern_frame(path,nk,nb,lat,rot)
      use, intrinsic :: ieee_arithmetic
      implicit none
      character*(*) path
      integer nk,nb,ios,unit,i,j,p,nrot,nlat,on,onb,foundk,foundb
      real*8 lat(3,3),rot(3,3),ol(3,3),test(3,3),det
      character*1024 line
      character*16 lx,ly,lz
      unit=97;nrot=0;nlat=0;foundk=0;foundb=0
      open(unit,file=trim(path),status='old',action='read',iostat=ios)
      if(ios.ne.0)call spinchern_error('cannot read matching OUTCAR spin frame',0)
      do
       read(unit,'(A)',iostat=ios)line
       if(ios.lt.0)exit
       if(ios.ne.0)call spinchern_error('OUTCAR read failed',0)
       p=index(line,'NKPTS')
       if(p.gt.0)then
        j=index(line(p:),'=')
        if(j.gt.0)then
         read(line(p+j:),*,iostat=ios)on
         if(ios.eq.0)then
          foundk=foundk+1
          if(on.ne.nk)call spinchern_error('OUTCAR NKPTS mismatches WAVECAR',0)
         endif
        endif
       endif
       p=index(line,'NBANDS')
       if(p.gt.0)then
        j=index(line(p:),'=')
        if(j.gt.0)then
         read(line(p+j:),*,iostat=ios)onb
         if(ios.eq.0)then
          foundb=foundb+1
          if(onb.ne.nb)call spinchern_error('OUTCAR NBANDS mismatches WAVECAR',0)
         endif
        endif
       endif
       if(index(line,'direct lattice vectors').gt.0)then
        do i=1,3
         read(unit,*,iostat=ios)ol(:,i)
         if(ios.ne.0)call spinchern_error('malformed OUTCAR lattice',0)
        enddo
        if(.not.all(ieee_is_finite(ol)).or.
     &     maxval(abs(ol-lat)).gt.1d-5)then
         call spinchern_error('OUTCAR lattice mismatches WAVECAR',0)
        endif
        nlat=nlat+1
       endif
       if(index(line,'transformation matrix from SAXIS to cartesian').gt.0)then
        nrot=nrot+1
        read(unit,'(A)',iostat=ios)line
        if(ios.ne.0)call spinchern_error('missing OUTCAR spin-frame separator',0)
        do i=1,3
         read(unit,*,iostat=ios)rot(i,1),lx,rot(i,2),ly,rot(i,3),lz
         if(ios.ne.0)call spinchern_error('malformed OUTCAR spin frame',0)
         if(lx.ne.'m_x'.or.ly.ne.'m_y'.or.lz.ne.'m_z')then
          call spinchern_error('unknown OUTCAR spin-frame format',0)
         endif
        enddo
       endif
      enddo
      close(unit)
      if(nrot.ne.1.or.nlat.lt.1.or.foundk.lt.1.or.foundb.lt.1)then
       call spinchern_error('OUTCAR needs unambiguous spin frame and source metadata',0)
      endif
      if(.not.all(ieee_is_finite(rot)))then
       call spinchern_error('nonfinite OUTCAR spin frame',0)
      endif
      test=matmul(transpose(rot),rot)
      do i=1,3
       test(i,i)=test(i,i)-1d0
      enddo
      det=rot(1,1)*(rot(2,2)*rot(3,3)-rot(3,2)*rot(2,3))-
     &    rot(1,2)*(rot(2,1)*rot(3,3)-rot(3,1)*rot(2,3))+
     &    rot(1,3)*(rot(2,1)*rot(3,2)-rot(3,1)*rot(2,2))
      if(maxval(abs(test)).gt.2d-5.or.abs(det-1d0).gt.2d-5)then
       call spinchern_error('OUTCAR spin frame is not a proper rotation',0)
      endif
! Printed OUTCAR precision is finite; normalize the transformed axis later.
      end subroutine spinchern_frame

      subroutine spinchern_write(wavepath,outpath,prefix,nx,ny,first,last,
     & nk,r,wk,axis,rot,etol,stol,gap,indirect,occupied,
     & grammin,gramcond,eval,neg,pos,flux,area,lsv,cov,chern)
      implicit none
      character*(*) wavepath,outpath,prefix
      integer nx,ny,first,last,nk,r,neg,pos,k,j,ios,u
      real*8 wk(3,nk),axis(3),rot(3,3),etol,stol,gap(nk),indirect
      real*8 grammin(nk),gramcond(nk),eval(r,nk),flux(3,nk),area
      real*8 lsv(2,nk),cov(2,nk),chern(3),center(3)
      logical occupied,exists
      character*512 paths(3)
      paths(1)=trim(prefix)//'_CHERN.csv'
      paths(2)=trim(prefix)//'_BERRY.csv'
      paths(3)=trim(prefix)//'_SPECTRUM.csv'
      do j=1,3
       inquire(file=trim(paths(j)),exist=exists)
       if(exists)call spinchern_error('output exists; choose a new --output prefix',0)
      enddo
      do j=1,3
       u=90+j
       open(u,file=trim(paths(j)),status='new',action='write',iostat=ios)
       if(ios.ne.0)call spinchern_error('cannot create new output',0)
       write(u,'(A)',iostat=ios)'# schema=VASPBERRY_SPIN_CHERN_V1'
       if(ios.ne.0)call spinchern_error('output write failed',0)
       write(u,'(A)')'# producer=VASPBERRY 1.6.1'
       write(u,'(A)')'# build_status=RELEASE_1.6.1'
       write(u,'(A)')'# representation=WAVECAR_PSEUDO_GRAM'
       write(u,'(A)')'# physical_paw_validated=false'
       write(u,'(A)')'# scope=selected_subspace'
       write(u,'(A,L1)')'# source_occupations_match_integer_projector=',occupied
       write(u,'(A)')'# method=PROJECTED_PAULI_GENERALIZED_EIGEN_FUKUI'
       write(u,'(A)')'# spin_convention=C_spin=(C_plus-C_minus)/2'
       write(u,'(A)')'# spin_basis=Cartesian; local frame from OUTCAR'
       write(u,'(A)')'# flux_convention=-Arg(Ux Uy_next conj(Ux_next) conj(Uy))'
       write(u,'(A)')'# mesh_convergence=NOT_ESTABLISHED_BY_ONE_MESH'
       write(u,'(A,A)')'# wavecar=',trim(wavepath)
       write(u,'(A,A)')'# outcar=',trim(outpath)
       write(u,'(A,3ES24.16)')'# spin_axis_cartesian=',axis
       do k=1,3
        write(u,'(A,I0,A,3ES24.16)')'# spin_frame_row_',k,'=',rot(k,:)
       enddo
       write(u,'(A,I0,A,I0)')'# selected_bands=',first,':',last
       write(u,'(A,I0,A,I0)')'# mesh=',nx,',',ny
       write(u,'(A,ES24.16)')'# energy_gap_threshold_eV=',etol
       write(u,'(A,ES24.16)')'# spin_gap_threshold_dimensionless=',stol
       write(u,'(A)')'# link_min_singular_threshold=1e-6'
       write(u,'(A,ES24.16)')'# min_external_energy_gap_eV=',minval(gap)
       if(first.eq.1)then
        write(u,'(A,ES24.16)')'# global_indirect_gap_eV=',indirect
        write(u,'(A,L1)')'# has_global_energy_gap=',indirect.gt.etol
       else
        write(u,'(A)')'# has_global_energy_gap=NOT_APPLICABLE'
       endif
       write(u,'(A,ES24.16)')'# min_abs_spin_eigenvalue=',minval(abs(eval))
       write(u,'(A,ES24.16)')'# min_gram_eigenvalue=',minval(grammin)
       write(u,'(A,ES24.16)')'# max_gram_condition=',maxval(gramcond)
       write(u,'(A,ES24.16)')'# min_link_singular_value=',minval(lsv)
       write(u,'(A,ES24.16)')'# min_G_coverage_fraction=',minval(cov)
       write(u,'(A,ES24.16)')'# max_abs_plaquette_flux_rad=',maxval(abs(flux))
       write(u,'(A)')'# periodic_sewing=INTEGER_G_SHIFT_VALIDATED'
      enddo
      write(92,'(A)')'k_index,k1,k2,k3,flux_plus_rad,flux_minus_rad,'//
     & 'flux_parent_rad,plaquette_area_A-2,omega_plus_A2,omega_minus_A2'
      do k=1,nk
       center=wk(:,k)+(/0.5d0/nx,0.5d0/ny,0d0/)
       center=modulo(center,1d0)
       write(92,'(I0,9(",",ES24.16E3))',iostat=ios)k,center,
     &     flux(:,k),area,flux(1,k)/area,flux(2,k)/area
       if(ios.ne.0)call spinchern_error('Berry CSV write failed',k)
      enddo
      write(93,'(A)')'k_index,k1,k2,k3,spin_eigen_index,'//
     & 'projected_pauli_eigenvalue,rank_plus,rank_minus,'//
     & 'external_energy_gap_eV,min_gram_eigenvalue,gram_condition'
      do k=1,nk
       do j=1,r
        write(93,'(I0,3(",",ES24.16E3),",",I0,",",ES24.16E3,'//
     &       '2(",",I0),3(",",ES24.16E3))',iostat=ios)
     &    k,wk(:,k),j,eval(j,k),pos,neg,gap(k),grammin(k),gramcond(k)
        if(ios.ne.0)call spinchern_error('spectrum CSV write failed',k)
       enddo
      enddo
      do u=92,93
       write(u,'(A)',iostat=ios)'# result_status=PASS'
       if(ios.ne.0)call spinchern_error('CSV completion write failed',0)
       close(u,iostat=ios)
       if(ios.ne.0)call spinchern_error('CSV close failed',0)
      enddo
      write(91,'(A)')'C_plus,C_minus,C_charge,delta_C,C_spin,'//
     & 'C_parent,rank_plus,rank_minus'
      write(91,'(ES24.16E3,5(",",ES24.16E3),2(",",I0))',iostat=ios)
     & chern(1),chern(2),chern(1)+chern(2),chern(1)-chern(2),
     & (chern(1)-chern(2))/2d0,chern(3),pos,neg
      if(ios.ne.0)call spinchern_error('summary CSV write failed',0)
      write(91,'(A)',iostat=ios)'# result_status=PASS'
      if(ios.ne.0)call spinchern_error('summary completion failed',0)
      close(91,iostat=ios)
      if(ios.ne.0)call spinchern_error('summary close failed',0)
      end subroutine spinchern_write
