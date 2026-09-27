! Local spin-sector differential curvature from an explicitly approximate
! canonical-momentum response tangent. Not an exact PAW/projector derivative.
      subroutine spinkubo_run(wavepath,outpath,prefix,nx,ny,mesh_seen,
     & sum_max,first,last,ispin,components,nk,nb,ecut,a1,a2,a3,axis,
     & etol,stol,rank,nproc)
      use, intrinsic :: ieee_arithmetic
      implicit none
#ifdef MPI_USE
      include 'mpif.h'
#endif
      character*(*) wavepath,outpath,prefix
      integer nx,ny,mesh_seen,first,last,ispin,components,nk,nb
      integer rank,nproc,r,k,n,m,i,a,ng,npmax,nbmax(3),ios,ierr
      integer neg,pos,sum_max,nsum
      integer, allocatable :: npl(:),g(:,:)
      real*8 ecut,a1(3),a2(3),a3(3),axis(3),etol,stol
      real*8 lat(3,3),b(3,3),rot(3,3),localaxis(3),area(3),pi
      real*8 volume,cross(3),fac,np,kcart(3),gap
      real*8, allocatable :: wk(:,:),en(:,:),occ(:),eval(:,:)
      real*8, allocatable :: omega(:,:,:),base(:,:,:),mix(:,:,:)
      real*8, allocatable :: gmin(:),gcond(:),egap(:)
      complex*16, allocatable :: ec(:),c(:,:),dc(:,:,:),weighted(:,:)
      complex*16, allocatable :: d(:,:)
      r=last-first+1
      if(ispin.ne.1.or.components.ne.2)then
       call spinkubo_error('requires ISPIN=1 two-component spinors',0)
      endif
      if(first.lt.1.or.last.gt.nb.or.r.lt.1.or.r.eq.nb)then
       call spinkubo_error('select a proper group within 1:NBANDS',0)
      endif
      nsum=nb
      if(sum_max.gt.0)nsum=sum_max
      if(nsum.lt.last.or.nsum.gt.nb.or.nsum-r.lt.1)then
       call spinkubo_error('sum-bands must include selected bands and an external state, within NBANDS',0)
      endif
      if(int(nb,8)*nk.gt.32000000_8)then
       call spinkubo_error('source metadata exceeds bounded storage',0)
      endif
      pi=4d0*atan(1d0)
      lat(:,1)=a1;lat(:,2)=a2;lat(:,3)=a3
      cross=(/a2(2)*a3(3)-a2(3)*a3(2),a2(3)*a3(1)-a2(1)*a3(3),
     &        a2(1)*a3(2)-a2(2)*a3(1)/)
      volume=dot_product(a1,cross)
      if(.not.ieee_is_finite(volume).or.abs(volume).le.1d-12)then
       call spinkubo_error('invalid lattice volume',0)
      endif
      b(:,1)=2d0*pi*cross/volume
      cross=(/a3(2)*a1(3)-a3(3)*a1(2),a3(3)*a1(1)-a3(1)*a1(3),
     &        a3(1)*a1(2)-a3(2)*a1(1)/)
      b(:,2)=2d0*pi*cross/volume
      cross=(/a1(2)*a2(3)-a1(3)*a2(2),a1(3)*a2(1)-a1(1)*a2(3),
     &        a1(1)*a2(2)-a1(2)*a2(1)/)
      b(:,3)=2d0*pi*cross/volume
      cross=(/b(2,1)*b(3,2)-b(3,1)*b(2,2),
     &        b(3,1)*b(1,2)-b(1,1)*b(3,2),
     &        b(1,1)*b(2,2)-b(2,1)*b(1,2)/)
      if(.not.all(ieee_is_finite(b)).or.
     &   .not.all(ieee_is_finite(cross)))then
       call spinkubo_error('nonfinite reciprocal geometry',0)
      endif
      area=0d0
      if(mesh_seen.eq.1)then
       if(nx.lt.2.or.ny.lt.2.or.int(nx,8)*ny.ne.nk)then
        call spinkubo_error('explicit mesh needs complete NX by NY points',0)
       endif
       area=cross/dble(nx)/dble(ny)
      endif
      call spinchern_frame(outpath,nk,nb,lat,rot)
      localaxis=matmul(transpose(rot),axis)
      localaxis=localaxis/sqrt(sum(localaxis**2))
      call reciproperty(nbmax,npmax,b(:,1),b(:,2),b(:,3),ecut,2)
      allocate(wk(3,nk),en(nb,nk),npl(nk),ec(nb),occ(nb),egap(nk))
      allocate(eval(r,nk),omega(3,3,nk),base(2,3,nk),mix(2,3,nk))
      allocate(gmin(nk),gcond(nk))
      do k=1,nk
       read(10,rec=3+(k-1)*(nb+1),iostat=ios)np,wk(:,k),
     &      (ec(i),occ(i),i=1,nb)
       if(ios.ne.0)call spinkubo_error('WAVECAR header read failed',k)
       if(.not.all(ieee_is_finite(wk(:,k))).or.
     &    .not.all(ieee_is_finite(real(ec,8))))then
        call spinkubo_error('nonfinite k point or energies',k)
       endif
       npl(k)=nint(np);en(:,k)=real(ec,8)
       if(any(en(2:,k).lt.en(:nb-1,k)-etol))then
        call spinkubo_error('source energies are not ordered',k)
       endif
! Cutting an external degenerate multiplet would make the partial sum gauge dependent.
       if(nsum.lt.nb)then
        if(abs(en(nsum+1,k)-en(nsum,k)).le.etol)then
         call spinkubo_error('sum-bands cuts an unresolved external degeneracy',k)
        endif
       endif
! Parent isolation always uses ALL source bands, including omitted sum states.
       gap=huge(1d0)
       do n=first,last
        do m=1,nb
         if(m.ge.first.and.m.le.last)cycle
         gap=min(gap,abs(en(n,k)-en(m,k)))
        enddo
       enddo
       egap(k)=gap
       if(gap.le.etol)then
        write(0,*)'spin-kubo external gap eV / threshold:',gap,etol
        call spinkubo_error('parent subspace is not energy isolated',k)
       endif
      enddo
      if(mesh_seen.eq.1)call spinkubo_mesh(wk,nk,nx,ny)
      omega=0d0;base=0d0;mix=0d0;eval=0d0;gmin=0d0;gcond=0d0
! Same bare-momentum factor as the existing standard native Kubo formula.
! D_a=hbar*v_a=(hbar^2/m_e)*(k+G)_a, in eV Angstrom.
      fac=(6.582119569d-16*2.99792458d8*1d10)**2/0.510998950d6
      do k=rank+1,nk,nproc
       if(int(npl(k),8)*(nsum+5_8*r).gt.64000000_8)then
        call spinkubo_error('response tangent exceeds bounded storage',k)
       endif
       allocate(c(npl(k),nsum),g(3,npmax),dc(npl(k),r,3))
       allocate(weighted(npl(k),r),d(nsum,r))
       call spinchern_read(k,1,nsum,nk,nb,npl(k),wk(:,k),
     &        nbmax,npmax,ecut,b,c,g,ng)
       do a=1,3
        do i=1,ng
         kcart=matmul(b,wk(:,k)+g(:,i))
         weighted(i,:)=fac*kcart(a)*c(i,first:last)
         weighted(i+ng,:)=fac*kcart(a)*c(i+ng,first:last)
        enddo
        d=matmul(conjg(transpose(c)),weighted)
        do n=first,last
         do m=1,nsum
          if(m.ge.first.and.m.le.last)then
           d(m,n-first+1)=(0d0,0d0)
          else
           d(m,n-first+1)=d(m,n-first+1)/(en(n,k)-en(m,k))
          endif
         enddo
        enddo
        dc(:,:,a)=matmul(c,d)
       enddo
       call spinkubo_curvature(c(:,first:last),dc,ng,r,localaxis,
     &      stol,omega(:,:,k),base(:,:,k),mix(:,:,k),eval(:,k),
     &      gmin(k),gcond(k),k)
       deallocate(c,g,dc,weighted,d)
      enddo
#ifdef MPI_USE
      call MPI_ALLREDUCE(MPI_IN_PLACE,omega,9*nk,
     & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,base,6*nk,
     & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,mix,6*nk,
     & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,eval,r*nk,
     & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,gmin,nk,
     & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,gcond,nk,
     & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif
      neg=count(eval(:,1).lt.-stol);pos=r-neg
      do k=1,nk
       if(count(eval(:,k).lt.-stol).ne.neg)then
        call spinkubo_error('projected-spin sector rank changes',k)
       endif
      enddo
      if(rank.eq.0)then
       call spinkubo_write(wavepath,outpath,prefix,nk,nb,r,first,last,
     &  ecut,wk,en,b,axis,rot,etol,stol,egap,gmin,gcond,eval,
     &  omega,base,mix,neg,pos,mesh_seen,nx,ny,area,nsum)
       write(6,*)'# Spin Kubo result_status=PASS; canonical tangent proxy'
      endif
      deallocate(wk,en,npl,ec,occ,egap,eval,omega,base,mix,gmin,gcond)
      end subroutine spinkubo_run

! Public math kernel: dc must be the derivative/proxy of this SAME raw C.
! omega(first,second): first=plus,minus,parent; second=yz,zx,xy.
! base/mix first=plus,minus. No derivative of eigenvector gauge is needed.
      subroutine spinkubo_curvature(c,dc,ng,r,axis,stol,
     & omega,base,mix,eval,gmin,gcond,k)
      use, intrinsic :: ieee_arithmetic
      implicit none
      integer ng,r,k,a,j,n,m,s,aa(3),bb(3)
      real*8 axis(3),stol,omega(3,3),base(2,3),mix(2,3)
      real*8 eval(r),gmin,gcond,value,scale
      complex*16 c(2*ng,r),dc(2*ng,r,3)
      complex*16, allocatable :: v(:,:),u(:,:),tangent(:,:,:)
      complex*16, allocatable :: su(:,:),sb(:,:),ds(:,:,:)
      complex*16 z
      if(.not.all(ieee_is_finite(real(dc,8))).or.
     & .not.all(ieee_is_finite(aimag(dc))))then
       call spinkubo_error('nonfinite response tangent',k)
      endif
      allocate(v(r,r),u(2*ng,r),tangent(2*ng,r,3))
      allocate(su(2*ng,r),sb(2*ng,r),ds(r,r,3))
      call spinchern_split(c,ng,r,axis,v,eval,gmin,gcond,k)
      if(minval(abs(eval)).le.stol)then
       call spinkubo_error('projected-spin zero gap is closed',k)
      endif
      u=matmul(c,v)
      do a=1,3
       tangent(:,:,a)=matmul(dc(:,:,a),v)
       tangent(:,:,a)=tangent(:,:,a)-
     &       matmul(u,matmul(conjg(transpose(u)),tangent(:,:,a)))
      enddo
      call spinkubo_sigma(u,ng,r,axis,su)
      do a=1,3
       call spinkubo_sigma(tangent(:,:,a),ng,r,axis,sb)
       ds(:,:,a)=matmul(conjg(transpose(tangent(:,:,a))),su)+
     &          matmul(conjg(transpose(u)),sb)
      enddo
      aa=(/2,3,1/);bb=(/3,1,2/)
      omega=0d0;base=0d0;mix=0d0
      do j=1,3
       do n=1,r
        s=1
        if(eval(n).lt.0d0)s=2
        z=dot_product(tangent(:,n,aa(j)),tangent(:,n,bb(j)))
        base(s,j)=base(s,j)-2d0*aimag(z)
        do m=1,r
! No division within one sector, even for degenerate projected spin.
         if((eval(m).gt.0d0).eqv.(eval(n).gt.0d0))cycle
         z=conjg(ds(m,n,aa(j)))*ds(m,n,bb(j))/
     &                          (eval(n)-eval(m))**2
         mix(s,j)=mix(s,j)-2d0*aimag(z)
        enddo
       enddo
       omega(1:2,j)=base(:,j)+mix(:,j)
       omega(3,j)=sum(base(:,j))
      enddo
      if(.not.all(ieee_is_finite(omega)).or.
     &   .not.all(ieee_is_finite(base)).or.
     &   .not.all(ieee_is_finite(mix)))then
       call spinkubo_error('nonfinite spin-sector curvature',k)
      endif
      scale=max(1d0,maxval(abs(base)),maxval(abs(mix)))
      if(maxval(abs(sum(mix,dim=1))).gt.1d-10*scale.or.
     & maxval(abs(omega(1,:)+omega(2,:)-omega(3,:))).gt.1d-10*scale)then
       call spinkubo_error('sector curvature sum identity failed',k)
      endif
      deallocate(v,u,tangent,su,sb,ds)
      end subroutine spinkubo_curvature

      subroutine spinkubo_sigma(c,ng,r,axis,out)
      implicit none
      integer ng,r
      real*8 axis(3)
      complex*16 c(2*ng,r),out(2*ng,r)
      out(:ng,:)=axis(3)*c(:ng,:)+
     &               cmplx(axis(1),-axis(2),8)*c(ng+1:,:)
      out(ng+1:,:)=cmplx(axis(1),axis(2),8)*c(:ng,:)-
     &               axis(3)*c(ng+1:,:)
      end subroutine spinkubo_sigma

      subroutine spinkubo_error(message,k)
      implicit none
      character*(*) message
      integer k
      write(0,*)'*** error - spin-kubo: ',trim(message),'; k=',k
      call vaspberry_fail
      end subroutine spinkubo_error

      subroutine spinkubo_mesh(wk,nk,nx,ny)
      implicit none
      integer nk,nx,ny,k,ix,iy
      integer grid(nx,ny)
      real*8 wk(3,nk),offset(2),u,v
      offset(1)=minval(modulo(wk(1,:),1d0))
      offset(2)=minval(modulo(wk(2,:),1d0))
      grid=0
      do k=1,nk
       u=modulo(wk(1,k)-offset(1),1d0)*nx
       v=modulo(wk(2,k)-offset(2),1d0)*ny
       if(abs(u-nint(u)).gt.1d-7.or.abs(v-nint(v)).gt.1d-7.or.
     &    abs(wk(3,k)-wk(3,1)-nint(wk(3,k)-wk(3,1))).gt.1d-7)then
        call spinkubo_error('nonuniform or nonplanar explicit mesh',k)
       endif
       ix=modulo(nint(u),nx)+1;iy=modulo(nint(v),ny)+1
       if(grid(ix,iy).ne.0)call spinkubo_error('duplicate mesh point',k)
       grid(ix,iy)=k
      enddo
      if(any(grid.eq.0))call spinkubo_error('missing mesh point',0)
      end subroutine spinkubo_mesh

      subroutine spinkubo_write(wavepath,outpath,prefix,nk,nb,r,first,last,
     & ecut,wk,en,b,axis,rot,etol,stol,gap,gmin,gcond,eval,
     & omega,base,mix,neg,pos,mesh_seen,nx,ny,area,nsum)
      use, intrinsic :: ieee_arithmetic
      implicit none
      character*(*) wavepath,outpath,prefix
      integer nk,nb,r,first,last,neg,pos,mesh_seen,nx,ny,k,j,u,ios,nfiles
      integer nsum
      real*8 ecut,wk(3,nk),en(nb,nk),b(3,3),axis(3),rot(3,3)
      real*8 etol,stol,gap(nk),gmin(nk),gcond(nk),eval(r,nk)
      real*8 omega(3,3,nk),base(2,3,nk),mix(2,3,nk),area(3)
      real*8 kcart(3),distance,delta(3),estimate(3),pi,values(21)
      real*8 integral_values(6)
      logical exists
      character*512 paths(3)
! Validate accumulated/derived output values before creating any result file.
      pi=4d0*atan(1d0);estimate=0d0;distance=0d0
      if(.not.all(ieee_is_finite(area)))then
       call spinkubo_error('nonfinite oriented area',0)
      endif
      do k=1,nk
       kcart=matmul(b,wk(:,k))
       if(.not.all(ieee_is_finite(kcart)))then
        call spinkubo_error('nonfinite Cartesian k coordinate',k)
       endif
       if(k.gt.1)then
        delta=matmul(b,wk(:,k)-wk(:,k-1))
        distance=distance+sqrt(sum(delta**2))
        if(.not.ieee_is_finite(distance))then
         call spinkubo_error('nonfinite accumulated path distance',k)
        endif
       endif
       if(mesh_seen.eq.1)then
        estimate=estimate+matmul(omega(:,:,k),area)/(2d0*pi)
        if(.not.all(ieee_is_finite(estimate)))then
         call spinkubo_error('nonfinite accumulated curvature integral',k)
        endif
       endif
      enddo
      integral_values=(/estimate(1),estimate(2),
     & estimate(1)+estimate(2),estimate(1)-estimate(2),
     & (estimate(1)-estimate(2))/2d0,estimate(3)/)
      if(.not.all(ieee_is_finite(integral_values)))then
       call spinkubo_error('nonfinite derived sector integral',0)
      endif
      paths(1)=trim(prefix)//'_KUBO.csv'
      paths(2)=trim(prefix)//'_KUBO_SPECTRUM.csv'
      paths(3)=trim(prefix)//'_KUBO_INTEGRAL.csv'
      nfiles=2
      if(mesh_seen.eq.1)nfiles=3
      do j=1,3
! Reject stale integral even when the new run is only a path calculation.
       inquire(file=trim(paths(j)),exist=exists)
       if(exists)call spinkubo_error('output exists; choose a new --output prefix',0)
      enddo
      do j=1,nfiles
       u=90+j
       open(u,file=trim(paths(j)),status='new',action='write',iostat=ios)
       if(ios.ne.0)call spinkubo_error('cannot create output',0)
       write(u,'(A)')'# schema=VASPBERRY_SPIN_SECTOR_KUBO_V1'
       write(u,'(A)')'# producer=VASPBERRY 1.6.1'
       write(u,'(A)')'# build_status=RELEASE_1.6.1'
       write(u,'(A)')'# derivative=CANONICAL_MOMENTUM_DERIVATIVE_PROXY'
       write(u,'(A)')'# representation=PSEUDO_GRAM'
       write(u,'(A)')'# physical_paw_validated=false'
       write(u,'(A)')'# exact_pseudo_projector_derivative=false'
       write(u,'(A)')'# observable=PROJECTED_SPIN_SECTOR_GEOMETRIC_CURVATURE_PROXY'
       write(u,'(A)')'# conventional_spin_current_hall=false'
       write(u,'(A)')'# scope=selected_subspace'
       write(u,'(A)')'# source_occupations=NOT_USED_FIXED_RANK_PROJECTOR'
       write(u,'(A)')'# source_band_sum=FINITE_RETAINED_WAVECAR_BANDS'
       write(u,'(A)')'# missing_terms=PAW_AUGMENTATION_NONLOCAL_SOC_VERTEX_AND_BASIS_DERIVATIVE'
       write(u,'(A)')'# geometric_integer_certified=false'
       write(u,'(A)')'# integer_rounding_applied=false'
       write(u,'(A)')'# spin_convention=C_spin=(C_plus-C_minus)/2'
       write(u,'(A)')'# curvature_convention=-2Im(inner_product_of_derivatives)'
       write(u,'(A)')'# decomposition=parent_projected_plus_spin_mixing'
       write(u,'(A)')'# spin_mixing_includes_projected_spin_derivative=true'
       write(u,'(A)')'# tensor_components=yz,zx,xy; Cartesian derivatives'
       write(u,'(A,A)')'# wavecar=',trim(wavepath)
       write(u,'(A,A)')'# outcar=',trim(outpath)
       write(u,'(A,ES24.16)')'# source_encut_eV=',ecut
       write(u,'(A,I0)')'# source_nbands=',nb
       write(u,'(A,I0)')'# sum_band_max=',nsum
       write(u,'(A,I0)')'# summed_external_bands=',nsum-r
       write(u,'(A)')'# parent_isolation_checked_against=ALL_STORED_SOURCE_BANDS'
       if(nsum.lt.nb)then
        write(u,'(A,ES24.16)')'# sum_boundary_gap_eV=',
     &   minval(abs(en(nsum+1,:)-en(nsum,:)))
       endif
       write(u,'(A,I0,A,I0)')'# selected_bands=',first,':',last
       write(u,'(A,3ES24.16)')'# spin_axis_cartesian=',axis
       do k=1,3
        write(u,'(A,I0,A,3ES24.16)')'# spin_frame_row_',k,'=',rot(k,:)
       enddo
       write(u,'(A,ES24.16)')'# energy_gap_threshold_eV=',etol
       write(u,'(A,ES24.16)')'# spin_gap_threshold_dimensionless=',stol
       write(u,'(A,ES24.16)')'# min_external_energy_gap_eV=',minval(gap)
       write(u,'(A,ES24.16)')'# min_abs_spin_eigenvalue=',minval(abs(eval))
       if(mesh_seen.eq.1)then
        write(u,'(A)')'# integration=EXPLICIT_UNIFORM_FULL_2D_MESH_ESTIMATE'
        write(u,'(A,I0,A,I0)')'# mesh=',nx,',',ny
        write(u,'(A,3ES24.16)')'# oriented_cell_area_A-2=',area
       else
        write(u,'(A)')'# integration=NONE_ARBITRARY_K_LIST'
       endif
      enddo
      write(91,'(A)')'k_index,k1,k2,k3,kx_A-1,ky_A-1,kz_A-1,k_distance_A-1,'//
     & 'omega_plus_yz_A2,omega_minus_yz_A2,omega_parent_yz_A2,'//
     & 'omega_plus_zx_A2,omega_minus_zx_A2,omega_parent_zx_A2,'//
     & 'omega_plus_xy_A2,omega_minus_xy_A2,omega_parent_xy_A2,'//
     & 'parent_projected_plus_yz_A2,parent_projected_minus_yz_A2,'//
     & 'parent_projected_plus_zx_A2,parent_projected_minus_zx_A2,'//
     & 'parent_projected_plus_xy_A2,parent_projected_minus_xy_A2,'//
     & 'spin_mixing_plus_yz_A2,spin_mixing_minus_yz_A2,'//
     & 'spin_mixing_plus_zx_A2,spin_mixing_minus_zx_A2,'//
     & 'spin_mixing_plus_xy_A2,spin_mixing_minus_xy_A2,'//
     & 'rank_plus,rank_minus,external_energy_gap_eV,min_abs_spin_eigenvalue,'//
     & 'min_gram_eigenvalue,gram_condition,selected_energy_min_eV,selected_energy_max_eV'
      distance=0d0
      do k=1,nk
       kcart=matmul(b,wk(:,k))
       if(k.gt.1)then
        delta=matmul(b,wk(:,k)-wk(:,k-1))
        distance=distance+sqrt(sum(delta**2))
       endif
       values(:9)=reshape(omega(:,:,k),(/9/))
       values(10:15)=reshape(base(:,:,k),(/6/))
       values(16:21)=reshape(mix(:,:,k),(/6/))
       write(91,'(I0,28(",",ES24.16E3),2(",",I0),6(",",ES24.16E3))',
     & iostat=ios)k,wk(:,k),kcart,distance,values,pos,neg,gap(k),
     & minval(abs(eval(:,k))),gmin(k),gcond(k),en(first,k),en(last,k)
       if(ios.ne.0)call spinkubo_error('point CSV write failed',k)
      enddo
      write(92,'(A)')'k_index,k1,k2,k3,spin_eigen_index,'//
     & 'projected_pauli_eigenvalue,rank_plus,rank_minus'
      do k=1,nk
       do j=1,r
        write(92,'(I0,3(",",ES24.16E3),",",I0,",",ES24.16E3,2(",",I0))',
     &    iostat=ios)k,wk(:,k),j,eval(j,k),pos,neg
        if(ios.ne.0)call spinkubo_error('spectrum CSV write failed',k)
       enddo
      enddo
      if(mesh_seen.eq.1)then
       write(93,'(A)')'C_est_plus,C_est_minus,C_est_charge,delta_C_est,'//
     &                 'C_est_spin,C_est_parent,rank_plus,rank_minus'
       write(93,'(ES24.16E3,5(",",ES24.16E3),2(",",I0))',iostat=ios)
     &  integral_values,pos,neg
       if(ios.ne.0)call spinkubo_error('integral CSV write failed',0)
      endif
      do u=90+nfiles,91,-1
       write(u,'(A)',iostat=ios)'# result_status=PASS'
       if(ios.ne.0)call spinkubo_error('completion write failed',0)
       close(u,iostat=ios)
       if(ios.ne.0)call spinkubo_error('CSV close failed',0)
      enddo
      end subroutine spinkubo_write
