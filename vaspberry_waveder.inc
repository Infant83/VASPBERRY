! Native standard-WAVEDER charge Kubo route. Independent of legacy momentum.
! C retains the producer longitudinal PAW optical connection convention.
! Contract occupied or explicitly selected geometry from stored pairs.
! No square completion, Hermitian averaging or second energy denominator.
      subroutine waveder_kubo_run(wavepath,derpath,inpath,outpath,
     & csvpath,nx,ny,meshseen,first,last,ns,components,nk,nb,
     & a1,a2,a3,rank)
      use, intrinsic :: ieee_arithmetic
      implicit none
      interface
       subroutine wd_read_matrix(path,nb,nk,ns,nd,c)
        character*(*) path
        integer nb,nk,ns,nd
        complex*8, allocatable, intent(out) :: c(:,:,:,:,:)
       end subroutine wd_read_matrix
      end interface
      character*(*) wavepath,derpath,inpath,outpath,csvpath
      integer nx,ny,meshseen,first,last,ns,components,nk,nb,rank
      integer s,k,n,m,a,b,rec,ios,nd,r,count0,mult
      integer occupied(ns),kubo_per_band,kubo_outcar_seen
      integer kubo_selection_explicit
      character*256 kubo_band_spec
      common /kubo_source_ints/ kubo_per_band,kubo_outcar_seen
      common /kubo_selection_ints/ kubo_selection_explicit
      common /kubo_selection_chars/ kubo_band_spec
      real*8 a1(3),a2(3),a3(3),lat(3,3),np,nelect
      real*8 vbm,cbm,pi,vol,cross(3),bv(3,3),area(3),total
      real*8 chern(ns),value,term,correction,t,y
      real*8, allocatable :: q(:,:,:),en(:,:,:),occ(:,:,:)
      real*8, allocatable :: omega(:,:,:),gap(:,:),weights(:)
      complex*16, allocatable :: ec(:)
      complex*8, allocatable :: c(:,:,:,:,:)
      complex*16 cx,cy
      if(rank.ne.0)return
      pi=4d0*atan(1d0)
      lat(:,1)=a1;lat(:,2)=a2;lat(:,3)=a3
      if(ns.lt.1.or.ns.gt.2.or.nb.lt.2.or.nk.lt.1.or.
     &   int(ns,8)*nk*nb.gt.32000000_8)then
       call wd_error('unsupported source dimensions or storage bound')
      endif
      if(csvpath.eq.wavepath.or.csvpath.eq.derpath.or.
     &   csvpath.eq.inpath.or.csvpath.eq.outpath)then
       call wd_error('output CSV must differ from every source path')
      endif
      call wd_read_matrix(derpath,nb,nk,ns,nd,c)
      allocate(q(3,nk,ns),en(nb,nk,ns),occ(nb,nk,ns),ec(nb))
      allocate(omega(3,nk,ns),gap(nk,ns),weights(nk))
      do s=1,ns
       do k=1,nk
        rec=3+((s-1)*nk+k-1)*(nb+1)
        read(10,rec=rec,iostat=ios)np,q(:,k,s),
     &     (ec(n),occ(n,k,s),n=1,nb)
        if(ios.ne.0)call wd_error('cannot read WAVECAR state headers')
        en(:,k,s)=real(ec,8)
       enddo
      enddo
      if(.not.all(ieee_is_finite(q)).or.
     &   .not.all(ieee_is_finite(en)).or.
     &   .not.all(ieee_is_finite(occ)))then
       call wd_error('nonfinite WAVECAR state headers')
      endif
      do s=1,ns
       if(any(en(2:nb,:,s).lt.en(:nb-1,:,s)-1d-8))then
        call wd_error('WAVECAR band energies are not ordered')
       endif
      enddo
      mult=1
      if(ns.eq.1.and.components.eq.1)mult=2
      call wd_validate_run(inpath,outpath,ns,components,nk,nb,
     & lat,q,en,occ,mult,nelect,weights)
      if(kubo_selection_explicit.eq.1)then
       call wd_selected_run(wavepath,derpath,inpath,outpath,csvpath,
     &  nx,ny,meshseen,first,last,ns,components,nk,nb,nd,
     &  a1,a2,a3,q,en,weights,c,kubo_per_band,kubo_band_spec)
       return
      endif
! Source occupations may be smeared. They select a candidate leading rank;
! the integer NELECT and common global gap certify the T=0 filling below.
      vbm=-huge(1d0);cbm=huge(1d0)
      do s=1,ns
       count0=count(occ(:,1,s).gt.0.5d0)
       if(last.ne.999999)count0=last
       if(count0.lt.1.or.count0.ge.nb)then
        call wd_error('T=0 filling needs occupied and stored empty states')
       endif
       occupied(s)=count0
       do k=1,nk
        if(count(occ(:,k,s).gt.0.5d0).ne.count0)then
         call wd_error('source filling is not a fixed occupied bundle; '//
     &     'metallic/finite-temperature Kubo is unsupported')
        endif
        gap(k,s)=en(count0+1,k,s)-en(count0,k,s)
       enddo
       vbm=max(vbm,maxval(en(count0,:,s)))
       cbm=min(cbm,minval(en(count0+1,:,s)))
      enddo
      if(abs(dble(sum(occupied)*mult)-nelect).gt.5.1d-5)then
       call wd_error('T=0 occupied ranks disagree with OUTCAR NELECT')
      endif
      if(cbm.le.vbm)then
       call wd_error('no common global insulating gap; metallic Kubo '//
     &   'requires a separately validated matrix route')
      endif
      if(minval(gap).le.0.002d0)then
       call wd_error('occupied-empty gap must exceed standard producer '//
     &   'DEG_THRESHOLD=0.002 eV')
      endif
      if(any(occupied.gt.nd))then
       call wd_error('rectangular WAVEDER lacks an occupied ket column; '//
     &   'no missing elements are filled or inferred')
      endif
      omega=0d0
      do s=1,ns
       r=occupied(s)
       do k=1,nk
        do a=1,3
         b=mod(a,3)+1
         value=0d0;correction=0d0
         do n=1,r
          do m=r+1,nb
           cx=cmplx(c(m,n,k,s,b),kind=8)
           cy=cmplx(c(m,n,k,s,mod(b,3)+1),kind=8)
           term=-2d0*aimag(conjg(cx)*cy)
! Compensated accumulation, including valid zero optical elements.
           y=term-correction;t=value+y
           correction=(t-value)-y;value=t
          enddo
         enddo
         omega(a,k,s)=value
        enddo
       enddo
      enddo
      if(.not.all(ieee_is_finite(omega)))then
       call wd_error('nonfinite occupied-empty curvature contraction')
      endif
      chern=0d0;total=0d0
      if(meshseen.eq.1)then
       call wd_check_mesh(q(:,:,1),weights,nk,nx,ny)
       cross=(/a2(2)*a3(3)-a2(3)*a3(2),
     &         a2(3)*a3(1)-a2(1)*a3(3),
     &         a2(1)*a3(2)-a2(2)*a3(1)/)
       vol=dot_product(a1,cross)
       if(abs(vol).lt.1d-12)call wd_error('singular source lattice')
       bv(:,1)=2d0*pi*cross/vol
       cross=(/a3(2)*a1(3)-a3(3)*a1(2),
     &         a3(3)*a1(1)-a3(1)*a1(3),
     &         a3(1)*a1(2)-a3(2)*a1(1)/)
       bv(:,2)=2d0*pi*cross/vol
       area=(/bv(2,1)*bv(3,2)-bv(3,1)*bv(2,2),
     &        bv(3,1)*bv(1,2)-bv(1,1)*bv(3,2),
     &        bv(1,1)*bv(2,2)-bv(2,1)*bv(1,2)/)/dble(nk)
       do s=1,ns
        do k=1,nk
         chern(s)=chern(s)+dot_product(omega(:,k,s),area)/(2d0*pi)
        enddo
       enddo
       total=dble(mult)*sum(chern)
       if(.not.ieee_is_finite(total))call wd_error('nonfinite integral')
      endif
      call open_kubo_csv(csvpath,.true.)
      write(96,'(A)')'# schema=VASPBERRY_WAVEDER_KUBO_OCCUPIED_V1'
      write(96,'(A)')'# kubo_source=WAVEDER'
      write(96,'(A)')'# vaspberry_version=1.6.5'
      write(96,'(A)')'# source_operator=VASP_5.4.4_LONGITUDINAL_PAW_OPTICAL'
      write(96,'(A)')'# result_kind=T0_INSULATING_OCCUPIED_BUNDLE_TRACE'
      write(96,'(A)')'# complete_occupied_bundle=true'
      write(96,'(A)')'# representation=VASP longitudinal PAW optical connection C; Angstrom; source convention retained'
      write(96,'(A)')'# normalization=STANDARD_MINUS_TWO_IM'
      write(96,'(A)')'# curvature_formula=-2 Im sum_empty_m,occupied_n conj(C_mn_a)*C_mn_b'
      write(96,'(A)')'# energy_denominator=ALREADY_IN_CONNECTION; no extra division'
      write(96,'(A)')'# matrix_completion=NONE; internal occupied pairs excluded analytically'
      write(96,'(A)')'# source_precision=complex64; complex128 compensated accumulation'
      write(96,'(A)')'# source_association=USER_ASSERTED_SAME_RUN; consistency_checked'
      write(96,'(A)')'# association_limit=WAVEDER has no energies,kpoints,or producer hashes'
      write(96,'(A)')'# accuracy_scope=stored optical operator; full SOC/+U velocity not independently certified'
      write(96,'(A)')'# mesh_and_nbands_convergence=NOT_ESTABLISHED_BY_ONE_RUN'
      write(96,'(A)')'# producer_cluster_threshold_eV=0.002'
      write(96,'(A)')'# source_occupations=smearing_allowed; response uses certified integer T0 filling'
      write(96,'(A,A)')'# wavecar=',trim(wavepath)
      write(96,'(A,A)')'# waveder=',trim(derpath)
      write(96,'(A,A)')'# incar=',trim(inpath)
      write(96,'(A,A)')'# outcar=',trim(outpath)
      write(96,'(A,I0)')'# source_nbands=',nb
      write(96,'(A,I0)')'# source_ndbands=',nd
      write(96,'(A,I0)')'# source_nkpoints=',nk
      write(96,'(A,I0)')'# expected_rows=',ns*nk
      write(96,'(A,I0)')'# source_nspin=',ns
      write(96,'(A,I0)')'# spinor_components=',components
      write(96,'(A,I0)')'# physical_spin_multiplicity=',mult
      do s=1,ns
       write(96,'(A,I0,A,I0)')'# occupied_bands_spin_',s,'=',occupied(s)
      enddo
      write(96,'(A,ES24.16E3)')'# global_gap_eV=',cbm-vbm
      write(96,'(A,ES24.16E3)')'# vbm_eV_source_zero=',vbm
      write(96,'(A,ES24.16E3)')'# cbm_eV_source_zero=',cbm
      if(meshseen.eq.1)then
       write(96,'(A,I0,A,I0)')'# full_2d_mesh=',nx,',',ny
       write(96,'(A)')'# sheet_orientation=oriented reciprocal b1,b2 plane; no cell-height factor'
       write(96,'(A,ES24.16E3)')'# total_chern=',total
       write(96,'(A,ES24.16E3)')'# sheet_hall_e2_over_h=',-total
       do s=1,ns
        write(96,'(A,I0,A,ES24.16E3)')'# chern_spin_',s,'=',chern(s)
       enddo
      else
       write(96,'(A)')'# integration=NONE; source k-point trace only'
      endif
      write(96,'(A)',iostat=ios)'spin,k_index,kx_frac,ky_frac,kz_frac,'//
     & 'omega_z_A2,min_external_gap_eV'
      if(ios.ne.0)call wd_error('CSV header write failed')
      do s=1,ns
       do k=1,nk
        write(96,'(I0,",",I0,5(",",ES24.16E3))',
     &   iostat=ios)s,k,q(:,k,s),omega(3,k,s),gap(k,s)
        if(ios.ne.0)call wd_error('CSV data write failed')
       enddo
      enddo
      call finish_kubo_csv(csvpath,.true.)
      write(6,'(A,A)')'# Standard WAVEDER occupied trace: ',trim(csvpath)
      write(6,'(A,I0,A,I0)')'# Stored matrix: ',nb,' x ',nd
      if(meshseen.eq.1)then
       write(6,'(A,ES24.16E3)')'# Total Chern estimate: ',total
       write(6,'(A,ES24.16E3)')'# Sheet Hall [e^2/h]: ',-total
      endif
      end subroutine waveder_kubo_run

! Selected geometry uses unordered cross-pairs. The lower ket is stored
! whenever either orientation is available. Internal equal-weight pairs
! cancel without reconstructing missing entries in the high-high block.
      subroutine wd_selected_run(wavepath,derpath,inpath,outpath,
     & csvpath,nx,ny,meshseen,first,last,ns,components,nk,nb,nd,
     & a1,a2,a3,q,en,weights,c,perband,spec)
      use, intrinsic :: ieee_arithmetic
      implicit none
      character*(*) wavepath,derpath,inpath,outpath,csvpath,spec
      integer nx,ny,meshseen,first,last,ns,components,nk,nb,nd,perband
      integer i,j,k,s,a,b,d,r,nout,mult,ios,p,z,colon,lo,hi,n,id
      integer pos,bandrow(nb),ids(nb)
      character*256 token
      character*32 digits
      character*(12*nb) bandids
      logical selected(nb),needed,noexternal
      real*8 a1(3),a2(3),a3(3),q(3,nk,ns),en(nb,nk,ns),weights(nk)
      real*8 pi,vol,cross(3),bv(3,2),area(3),de,limit,term,sgn
      real*8 y,t,chern(ns),cluster_cutoff,denom_guard
      real*8, allocatable :: omega(:,:,:,:),corr(:,:,:,:),gap(:,:,:)
      complex*8 c(nb,nd,nk,ns,3)
      complex*16 forward(3),reverse(3)
      selected=.false.;bandrow=0
      if(len_trim(spec).eq.0)then
       hi=last
       if(hi.eq.999999)hi=nb
       if(first.lt.1.or.hi.lt.first.or.hi.gt.nb)
     &  call wd_error('selected bands are outside source NBANDS')
       selected(first:hi)=.true.
      else
       p=1
       do
        z=index(spec(p:),',')
        if(z.eq.0)then
         z=len_trim(spec)+1
        else
         z=p+z-1
        endif
        token=spec(p:z-1)
        colon=index(token,':')
        if(colon.eq.0)then
         read(token,*,iostat=ios)lo
         hi=lo
        else
         token(colon:colon)=','
         read(token,*,iostat=ios)lo,hi
        endif
        if(ios.ne.0.or.lo.lt.1.or.hi.lt.lo.or.hi.gt.nb)
     &   call wd_error('selected bands are outside source NBANDS')
        if(any(selected(lo:hi)))
     &   call wd_error('duplicate/overlapping selected bands')
        selected(lo:hi)=.true.
        if(z.gt.len_trim(spec))exit
        p=z+1
       enddo
      endif
      r=count(selected);n=0;bandids='';pos=1
      do i=1,nb
       if(.not.selected(i))cycle
       n=n+1;ids(n)=i;bandrow(i)=n
       write(digits,'(I0)')i
       if(n.gt.1)then
        bandids(pos:pos)=',';pos=pos+1
       endif
       bandids(pos:pos+len_trim(digits)-1)=trim(digits)
       pos=pos+len_trim(digits)
      enddo
      nout=1
      if(perband.eq.1)nout=r
      noexternal=r.eq.nb.and.perband.eq.0
      allocate(omega(3,nout,nk,ns),corr(3,nout,nk,ns),gap(nout,nk,ns))
      omega=0d0;corr=0d0;gap=huge(1d0)
      cluster_cutoff=0.002d0
      denom_guard=dble(real(1d-10,kind=4))
      do s=1,ns
       do k=1,nk
! Full-source transitive clusters: any membership change inside an
! adjacent-link cluster splits it, including through unselected bridges.
        do i=2,nb
         if(abs(en(i,k,s)-en(i-1,k,s)).le.cluster_cutoff)then
          if(perband.eq.1)then
           if(selected(i).or.selected(i-1))call wd_error(
     &      'per-band state belongs to an unresolved producer cluster; '//
     &      'select the complete cluster as a trace instead')
          else
           if(selected(i).neqv.selected(i-1))call wd_error(
     &      'selected trace splits a full-source transitive producer '//
     &      'cluster (DEG_THRESHOLD=0.002 eV)')
          endif
         endif
        enddo
! S selects the target states, not the virtual-state cutoff: visit
! every source pair. Equal unit weights cancel internal trace pairs;
! resolved band output instead needs those internal pairs as well.
        do i=1,nb-1
         do j=i+1,nb
          needed=selected(i).neqv.selected(j)
          if(perband.eq.1)needed=selected(i).or.selected(j)
          if(.not.needed)cycle
          if(i.gt.nd)call wd_error(
     &     'required selected-band pair is missing in the WAVEDER '//
     &     'high-high block; regenerate more ket columns, do not pad')
          de=abs(en(j,k,s)-en(i,k,s))
          if(de.le.denom_guard)call wd_error(
     &     'required pair was erased by the producer denominator guard')
! Optical C already contains the source-energy denominator. Its
! reverse orientation is fixed by conjugation; never divide again
! or fabricate the absent high-high block of a rectangular file.
          forward=cmplx(c(j,i,k,s,:),kind=8)
          if(j.le.nd)then
           reverse=conjg(cmplx(c(i,j,k,s,:),kind=8))
           do d=1,3
            limit=1d-10+1d-8*max(abs(forward(d)),abs(reverse(d)))*de
            if(abs(forward(d)-reverse(d))*de.gt.limit)then
             call wd_error('required pair violates strict D '//
     &        'Hermiticity/conjugation consistency')
            endif
           enddo
          endif
          if(perband.eq.0)then
           gap(1,k,s)=min(gap(1,k,s),de)
          else
           if(selected(i))gap(bandrow(i),k,s)=min(gap(bandrow(i),k,s),de)
           if(selected(j))gap(bandrow(j),k,s)=min(gap(bandrow(j),k,s),de)
          endif
          do a=1,3
           b=mod(a,3)+1
           term=-2d0*aimag(conjg(forward(b))*forward(mod(b,3)+1))
           if(perband.eq.0)then
            sgn=-1d0
            if(selected(i))sgn=1d0
            y=sgn*term-corr(a,1,k,s);t=omega(a,1,k,s)+y
            corr(a,1,k,s)=(t-omega(a,1,k,s))-y;omega(a,1,k,s)=t
           else
            if(selected(i))then
             id=bandrow(i)
             y=term-corr(a,id,k,s);t=omega(a,id,k,s)+y
             corr(a,id,k,s)=(t-omega(a,id,k,s))-y;omega(a,id,k,s)=t
            endif
            if(selected(j))then
             id=bandrow(j)
             y=-term-corr(a,id,k,s);t=omega(a,id,k,s)+y
             corr(a,id,k,s)=(t-omega(a,id,k,s))-y;omega(a,id,k,s)=t
            endif
           endif
          enddo
         enddo
        enddo
       enddo
      enddo
      if(.not.all(ieee_is_finite(omega)))call wd_error(
     & 'nonfinite selected-band curvature contraction')
      pi=4d0*atan(1d0);chern=0d0
      if(meshseen.eq.1)then
       call wd_check_mesh(q(:,:,1),weights,nk,nx,ny)
       cross=(/a2(2)*a3(3)-a2(3)*a3(2),
     &         a2(3)*a3(1)-a2(1)*a3(3),
     &         a2(1)*a3(2)-a2(2)*a3(1)/)
       vol=dot_product(a1,cross)
       if(abs(vol).lt.1d-12)call wd_error('singular source lattice')
       bv(:,1)=2d0*pi*cross/vol
       cross=(/a3(2)*a1(3)-a3(3)*a1(2),
     &         a3(3)*a1(1)-a3(1)*a1(3),
     &         a3(1)*a1(2)-a3(2)*a1(1)/)
       bv(:,2)=2d0*pi*cross/vol
       area=(/bv(2,1)*bv(3,2)-bv(3,1)*bv(2,2),
     &        bv(3,1)*bv(1,2)-bv(1,1)*bv(3,2),
     &        bv(1,1)*bv(2,2)-bv(2,1)*bv(1,2)/)/dble(nk)
       do s=1,ns
        do k=1,nk
         do n=1,nout
          chern(s)=chern(s)+dot_product(omega(:,n,k,s),area)/(2d0*pi)
         enddo
        enddo
       enddo
       if(.not.all(ieee_is_finite(chern)))call wd_error('nonfinite selected integral')
      endif
      mult=1
      if(ns.eq.1.and.components.eq.1)mult=2
      call open_kubo_csv(csvpath,.true.)
      if(perband.eq.1)then
       write(96,'(A)')'# schema=VASPBERRY_WAVEDER_KUBO_BAND_V1'
       write(96,'(A)')'# result_kind=ISOLATED_SELECTED_BAND_CURVATURE'
       write(96,'(A)')'# intermediate_bands=ALL_SOURCE_BANDS_EXCEPT_SELF'
      else
       write(96,'(A)')'# schema=VASPBERRY_WAVEDER_KUBO_BUNDLE_V1'
       write(96,'(A)')'# result_kind=ISOLATED_SELECTED_BUNDLE_TRACE'
       write(96,'(A)')'# intermediate_bands=EXTERNAL_TO_SELECTED_BUNDLE_WITHIN_SOURCE_NBANDS'
      endif
      write(96,'(A)')'# kubo_source=WAVEDER'
      write(96,'(A)')'# vaspberry_version=1.6.5'
      write(96,'(A)')'# source_operator=VASP_5.4.4_LONGITUDINAL_PAW_OPTICAL'
      write(96,'(A)')'# representation=VASP longitudinal PAW optical connection C; Angstrom; source convention retained'
      write(96,'(A)')'# normalization=STANDARD_MINUS_TWO_IM'
      write(96,'(A)')'# occupation_weighting=NONE'
      write(96,'(A)')'# required_pair_coverage=PASS'
      write(96,'(A)')'# required_pair_hermiticity=PASS; D_atol=1e-10; D_rtol=1e-8'
      write(96,'(A)')'# curvature_formula=unordered pair K_ij=-2 Im conj(C_ji_a)*C_ji_b; selection difference'
      write(96,'(A)')'# energy_denominator=ALREADY_IN_CONNECTION; no extra division'
      write(96,'(A)')'# matrix_completion=NONE; internal equal-weight pairs cancel analytically'
      write(96,'(A)')'# source_precision=complex64; complex128 compensated accumulation'
      write(96,'(A)')'# source_association=USER_ASSERTED_SAME_RUN; consistency_checked'
      write(96,'(A)')'# association_limit=WAVEDER has no energies,kpoints,or producer hashes'
      write(96,'(A)')'# accuracy_scope=stored optical operator; full SOC/+U velocity not independently certified'
      write(96,'(A)')'# mesh_and_nbands_convergence=NOT_ESTABLISHED_BY_ONE_RUN'
      write(96,'(A)')'# producer_cluster_threshold_eV=0.002'
      write(96,'(A)')'# producer_cluster_scope=FULL_SOURCE_TRANSITIVE_ADJACENCY'
      write(96,'(A)')'# source_occupations=not used for selected geometric curvature'
      write(96,'(A,A)')'# wavecar=',trim(wavepath)
      write(96,'(A,A)')'# waveder=',trim(derpath)
      write(96,'(A,A)')'# incar=',trim(inpath)
      write(96,'(A,A)')'# outcar=',trim(outpath)
      write(96,'(A,A)')'# band_ids=',trim(bandids)
      write(96,'(A,I0)')'# band_min=',ids(1)
      write(96,'(A,I0)')'# band_max=',ids(r)
      write(96,'(A,I0)')'# band_rank=',r
      write(96,'(A,I0)')'# source_nbands=',nb
      write(96,'(A,I0)')'# source_ndbands=',nd
      write(96,'(A,I0)')'# source_nkpoints=',nk
      write(96,'(A,I0)')'# expected_rows=',ns*nk*nout
      write(96,'(A,I0)')'# source_nspin=',ns
      write(96,'(A,I0)')'# spinor_components=',components
      write(96,'(A,I0)')'# physical_spin_multiplicity=',mult
      if(noexternal)then
       write(96,'(A)')'# no_external_states=true'
       write(96,'(A)')'# zero_trace_scope=TRUNCATED_WAVEDER_BASIS'
      else
       write(96,'(A)')'# no_external_states=false'
      endif
      if(meshseen.eq.1)then
       write(96,'(A,I0,A,I0)')'# full_2d_mesh=',nx,',',ny
       write(96,'(A)')'# sheet_orientation=oriented reciprocal b1,b2 plane; no cell-height factor'
       write(96,'(A)')'# integration=UNIT_WEIGHT_SELECTED_GEOMETRY; not total physical Hall'
       write(96,'(A,ES24.16E3)')'# selected_chern_sum=',sum(chern)
       do s=1,ns
        write(96,'(A,I0,A,ES24.16E3)')'# selected_chern_spin_',s,'=',chern(s)
       enddo
      else
       write(96,'(A)')'# integration=NONE; source k-point selected geometry only'
      endif
      if(perband.eq.1)then
       write(96,'(A)',iostat=ios)'spin,k_index,band,kx_frac,ky_frac,kz_frac,'//
     &  'energy_eV,omega_z_A2,min_gap_eV'
      else
       write(96,'(A)',iostat=ios)'spin,k_index,kx_frac,ky_frac,kz_frac,'//
     &  'omega_z_A2,min_external_gap_eV'
      endif
      if(ios.ne.0)call wd_error('CSV header write failed')
      do s=1,ns
       do k=1,nk
        if(perband.eq.1)then
         do n=1,r
          write(96,'(I0,",",I0,",",I0,6(",",ES24.16E3))',iostat=ios)
     &     s,k,ids(n),q(:,k,s),en(ids(n),k,s),omega(3,n,k,s),gap(n,k,s)
          if(ios.ne.0)call wd_error('CSV data write failed')
         enddo
        else if(noexternal)then
         write(96,'(I0,",",I0,4(",",ES24.16E3),",NA")',iostat=ios)
     &    s,k,q(:,k,s),omega(3,1,k,s)
        else
         write(96,'(I0,",",I0,5(",",ES24.16E3))',iostat=ios)
     &    s,k,q(:,k,s),omega(3,1,k,s),gap(1,k,s)
        endif
        if(ios.ne.0)call wd_error('CSV data write failed')
       enddo
      enddo
      call finish_kubo_csv(csvpath,.true.)
      write(6,'(A,A)')'# Standard WAVEDER selected geometry: ',trim(csvpath)
      write(6,'(A,A)')'# Selected band IDs: ',trim(bandids)
      write(6,'(A,I0,A,I0)')'# Stored matrix: ',nb,' x ',nd
      if(meshseen.eq.1)write(6,'(A,ES24.16E3)')
     & '# Selected unit-weight Chern contribution (not total Hall): ',sum(chern)
      end subroutine wd_selected_run

      subroutine wd_read_matrix(path,nb,nk,ns,nd,c)
      use, intrinsic :: ieee_arithmetic
      implicit none
      character*(*) path
      integer nb,nk,ns,nd,u,ios
      integer*4 marker,trailer,dims(4)
      integer*8 bytes,elements
      real*8 node,plasma(3,3)
      complex*8, allocatable, intent(out) :: c(:,:,:,:,:)
      logical exists
      inquire(file=trim(path),exist=exists,size=bytes,iostat=ios)
      if(ios.ne.0.or..not.exists)then
       call wd_error('WAVEDER is required for standard Kubo. Generate '//
     &   'it with LOPTICS, or explicitly select --kubo-source wavecar '//
     &   'to acknowledge the bare-momentum approximation')
      endif
      open(newunit=u,file=trim(path),status='old',action='read',
     & access='stream',form='unformatted',convert='little_endian',iostat=ios)
      if(ios.ne.0)call wd_error('cannot open WAVEDER')
      read(u,iostat=ios)marker,dims,trailer
      if(ios.ne.0.or.marker.ne.16.or.trailer.ne.16)then
       call wd_error('unsupported WAVEDER header/record markers; '//
     &   'only little-endian complex64 four-record format is supported')
      endif
      nd=dims(2)
      if(dims(1).ne.nb.or.dims(3).ne.nk.or.dims(4).ne.ns.or.
     &   nd.lt.1.or.nd.gt.nb)then
       call wd_error('WAVEDER dimensions disagree with WAVECAR')
      endif
      elements=int(nb,8)*nd*nk*ns*3
      if(elements.gt.32000000_8.or.bytes.ne.128_8+8_8*elements)then
       call wd_error('WAVEDER length/storage bound invalid; split records '//
     &   'and trailing data are unsupported')
      endif
      read(u,iostat=ios)marker,node,trailer
      if(ios.ne.0.or.marker.ne.8.or.trailer.ne.8)then
       call wd_error('invalid WAVEDER scalar record')
      endif
      read(u,iostat=ios)marker,plasma,trailer
      if(ios.ne.0.or.marker.ne.72.or.trailer.ne.72)then
       call wd_error('invalid WAVEDER plasmon record')
      endif
      allocate(c(nb,nd,nk,ns,3),stat=ios)
      if(ios.ne.0)call wd_error('cannot allocate WAVEDER storage')
      read(u,iostat=ios)marker,c,trailer
      if(ios.ne.0.or.int(marker,8).ne.elements*8.or.
     & int(trailer,8).ne.elements*8)then
       call wd_error('invalid/truncated WAVEDER matrix record')
      endif
      close(u,iostat=ios)
      if(ios.ne.0)call wd_error('cannot close WAVEDER input')
      if(.not.ieee_is_finite(node).or.
     & .not.all(ieee_is_finite(plasma)).or.
     & .not.all(ieee_is_finite(real(c))).or.
     & .not.all(ieee_is_finite(aimag(c))))then
       call wd_error('nonfinite WAVEDER input')
      endif
      end subroutine wd_read_matrix

      subroutine wd_check_mesh(q,weights,nk,nx,ny)
      implicit none
      integer nk,nx,ny,k,ix,iy,ios
      integer, allocatable :: grid(:,:)
      real*8 q(3,nk),weights(nk),offset(2),u,v,w
      if(nx.lt.1.or.ny.lt.1.or.int(nx,8)*ny.ne.nk)then
       call wd_error('explicit mesh must cover all NKPTS exactly')
      endif
      if(any(abs(weights-1d0/dble(nk)).gt.0.00051d0))then
       call wd_error('integration requires uniform full-grid OUTCAR weights')
      endif
      allocate(grid(nx,ny),stat=ios)
      if(ios.ne.0)call wd_error('cannot allocate mesh validation grid')
      grid=0
      offset(1)=minval(modulo(q(1,:),1d0))
      offset(2)=minval(modulo(q(2,:),1d0))
      do k=1,nk
       u=modulo(q(1,k)-offset(1),1d0)*nx
       v=modulo(q(2,k)-offset(2),1d0)*ny
       w=q(3,k)-q(3,1)
       if(abs(u-nint(u)).gt.1d-7.or.abs(v-nint(v)).gt.1d-7.or.
     &    abs(w-nint(w)).gt.1d-7)then
        call wd_error('integration requires a uniform periodic 2D mesh')
       endif
       ix=modulo(nint(u),nx)+1;iy=modulo(nint(v),ny)+1
       if(grid(ix,iy).ne.0)call wd_error('duplicate periodic mesh point')
       grid(ix,iy)=k
      enddo
      if(any(grid.eq.0))call wd_error('incomplete periodic mesh')
      end subroutine wd_check_mesh

! Validate the producer branch and compare effective OUTCAR state tables with
! WAVECAR. WAVEDER itself cannot authenticate its same-run file association.
      subroutine wd_validate_run(inpath,outpath,ns,components,nk,nb,
     & lat,q,en,occ,mult,nelect,weights)
      use, intrinsic :: ieee_arithmetic
      implicit none
      character*(*) inpath,outpath
      integer ns,components,nk,nb,mult,u,ios,i,j,k,s,idx,p,nassign
      integer tablecount,spinseen,nlat,ngrid,na,ib,val
      integer, parameter :: ntag=15
      character*32 tags(ntag),values(ntag),token
      character*64 names(512),key
      character*256 ivals(512)
      character*2048 line,work
      character*16 word
      logical found,complete,converged,optics,dielectric,version
      logical seen(ntag),tables(nk,ns),wd_bool,expected_bool
      real*8 latout(3,3),grid(4,nk),coord(3),energy,occupation
      real*8 nelect,x,old,wd_real
      real*8 lat(3,3),q(3,nk,ns),en(nb,nk,ns),occ(nb,nk,ns)
      real*8 weights(nk),record_occ
      tags=(/'ISPIN                           ',
     &       'NKPTS                           ',
     &       'NBANDS                          ',
     &       'LNONCOLLINEAR                   ',
     &       'LSORBIT                         ',
     &       'LNABLA                          ',
     &       'LREAL                           ',
     &       'LHFCALC                         ',
     &       'METAGGA                         ',
     &       'LEPSILON                        ',
     &       'LVEL                            ',
     &       'ISYM                            ',
     &       'NSW                             ',
     &       'NELECT                          ',
     &       'DEG_THRESHOLD                   '/)
      call wd_read_incar(inpath,nassign,names,ivals)
      do i=1,3
       select case(i)
       case(1)
        key='LOPTICS';expected_bool=.true.
       case(2)
        key='LPEAD';expected_bool=.false.
       case(3)
        key='LNABLA';expected_bool=.false.
       end select
       found=.false.
       do j=1,nassign
        if(names(j).eq.key)then
         found=.true.
         if(wd_bool(ivals(j)).neqv.expected_bool)then
          call wd_error('unsupported explicit INCAR '//trim(key))
         endif
        endif
       enddo
       if(.not.found)call wd_error('INCAR requires explicit '//trim(key))
      enddo
      do j=1,nassign
       if(names(j).eq.'LBERRY_EXPORT')then
        if(wd_bool(ivals(j)))call wd_error('standard WAVEDER requires '//
     &    'LBERRY_EXPORT false; custom optical producers are separate')
       endif
      enddo
      open(newunit=u,file=trim(outpath),status='old',action='read',iostat=ios)
      if(ios.ne.0)call wd_error('same-run completed OUTCAR is required')
      seen=.false.;tables=.false.
      complete=.false.;converged=.false.;optics=.false.
      dielectric=.false.;version=.false.
      s=1;spinseen=0;tablecount=0;nlat=0;ngrid=0
      do
       read(u,'(A)',iostat=ios)line
       if(ios.lt.0)exit
       if(ios.ne.0)call wd_error('OUTCAR text read failed')
       work=adjustl(line)
       if(index(work,'vasp.').eq.1)then
        if(index(work,'vasp.5.4.4.').ne.1.and.
     &     index(work,'vasp.5.4.4 ').ne.1)then
         call wd_error('only audited standard VASP 5.4.4 is supported')
        endif
        version=.true.
       endif
       if(index(line,'General timing and accounting informations for this job:').gt.0)
     &   complete=.true.
       if(index(line,'aborting loop because EDIFF is reached').gt.0)
     &   converged=.true.
       if(index(line,'OPTICS:').gt.0.and.index(line,'cpu time').gt.0)
     &   optics=.true.
       if(index(line,'frequency dependent IMAGINARY DIELECTRIC FUNCTION').gt.0)
     &   dielectric=.true.
       do i=1,ntag
        call wd_out_token(work,tags(i),token,found)
        if(.not.found)cycle
        if(seen(i))then
         if(i.ge.4.and.i.le.11)then
          if(wd_bool(values(i)).neqv.wd_bool(token))then
           call wd_error('ambiguous effective OUTCAR '//trim(tags(i)))
          endif
         else
          if(wd_real(values(i)).ne.wd_real(token))then
           call wd_error('ambiguous effective OUTCAR '//trim(tags(i)))
          endif
         endif
        endif
        values(i)=token;seen(i)=.true.
       enddo
       if(index(line,'direct lattice vectors').gt.0.and.
     &    index(line,'reciprocal lattice vectors').gt.0)then
        do i=1,3
         read(u,*,iostat=ios)latout(:,i)
         if(ios.ne.0)call wd_error('malformed OUTCAR lattice')
        enddo
        if(.not.all(ieee_is_finite(latout)).or.
     &     maxval(abs(lat-latout)).gt.6d-9)then
         call wd_error('OUTCAR lattice mismatches WAVECAR')
        endif
        nlat=nlat+1
       endif
       if(index(line,'k-points in reciprocal lattice and weights:').gt.0)then
        do k=1,nk
         read(u,*,iostat=ios)grid(:,k)
         if(ios.ne.0)call wd_error('malformed OUTCAR k-point grid')
        enddo
        if(.not.all(ieee_is_finite(grid)).or.any(grid(4,:).lt.0d0))then
         call wd_error('invalid OUTCAR k-point weights/coordinates')
        endif
        do s=1,ns
         if(maxval(abs(grid(:3,:)-q(:,:,s))).gt.6d-9)then
          call wd_error('OUTCAR k-point list mismatches WAVECAR')
         endif
        enddo
        s=1;ngrid=ngrid+1
       endif
       if(index(work,'spin component ').eq.1)then
        read(work(16:),*,iostat=ios)s
        if(ios.ne.0.or.s.lt.1.or.s.gt.ns)then
         call wd_error('invalid OUTCAR spin-component table')
        endif
        spinseen=spinseen+1
       endif
       if(index(work,'k-point ').eq.1.and.index(work,':').gt.0.and.
     &    index(work,'plane waves:').eq.0)then
        p=index(work,':')
        read(work(9:p-1),*,iostat=ios)k
        if(ios.ne.0.or.k.lt.1.or.k.gt.nk)then
         call wd_error('invalid OUTCAR eigenvalue k-point index')
        endif
        read(work(p+1:),*,iostat=ios)coord
        if(ios.ne.0.or..not.all(ieee_is_finite(coord)))then
         call wd_error('invalid OUTCAR eigenvalue coordinates')
        endif
        read(u,'(A)',iostat=ios)line
        if(ios.ne.0.or.index(line,'band No.').eq.0.or.
     &     index(line,'occupation').eq.0)then
         call wd_error('missing OUTCAR eigenvalue header')
        endif
        if(maxval(abs(coord-q(:,k,s))).gt.5.1d-5)then
         call wd_error('OUTCAR eigenvalue k-points mismatch WAVECAR')
        endif
! Earlier SCF eigenvalue tables may differ. Retain final tables below by
! recording mismatch status for each source (spin,k) and overwriting it.
        if(k.ne.mod(tablecount,nk)+1.or.
     &     s.ne.mod(tablecount/nk,ns)+1)then
         call wd_error('OUTCAR final eigenvalue spin/k order incomplete or ambiguous')
        endif
        if(k.eq.1.and.s.eq.1)tables=.false.
        tables(k,s)=.true.
        do i=1,nb
         read(u,*,iostat=ios)ib,energy,occupation
         if(ios.ne.0.or.ib.ne.i.or..not.ieee_is_finite(energy).or.
     &      .not.ieee_is_finite(occupation))then
          call wd_error('malformed OUTCAR eigenvalue/occupation table')
         endif
         if(abs(energy-en(i,k,s)).gt.5.1d-5.or.
     &      abs(occupation-dble(mult)*occ(i,k,s)).gt.5.1d-6)then
          tables(k,s)=.false.
         endif
        enddo
        tablecount=tablecount+1
       endif
      enddo
      close(u)
      if(.not.(complete.and.converged.and.optics.and.dielectric.and.version))then
       call wd_error('OUTCAR must show completed, electronically converged '//
     &   'standard VASP 5.4.4 static optical calculation')
      endif
      if(.not.all(seen))then
       do i=1,ntag
        if(.not.seen(i))call wd_error('OUTCAR lacks effective '//trim(tags(i)))
       enddo
      endif
      if(nlat.lt.1.or.ngrid.lt.1.or.tablecount.lt.nk*ns.or.
     & mod(tablecount,nk*ns).ne.0.or..not.all(tables).or.
     & (ns.eq.2.and.spinseen.lt.2))then
       call wd_error('OUTCAR final state tables incomplete or mismatch WAVECAR')
      endif
      if(wd_real(values(1)).ne.dble(ns).or.
     & wd_real(values(2)).ne.dble(nk).or.
     & wd_real(values(3)).ne.dble(nb))then
       call wd_error('OUTCAR dimensions mismatch WAVECAR')
      endif
      if(wd_bool(values(4)).neqv.(components.eq.2))then
       call wd_error('OUTCAR noncollinear flag mismatches WAVECAR spinors')
      endif
      if(wd_bool(values(5)).and.components.ne.2)then
       call wd_error('OUTCAR LSORBIT requires spinor WAVECAR')
      endif
      do i=6,11
       if(wd_bool(values(i)))call wd_error('unsupported optical branch: '//trim(tags(i)))
      enddo
      if(wd_real(values(12)).ne.-1d0.or.wd_real(values(13)).ne.0d0)then
       call wd_error('standard route requires effective ISYM=-1 and NSW=0')
      endif
      nelect=wd_real(values(14))
      if(nelect.le.0d0.or.wd_real(values(15)).ne.0.002d0)then
       call wd_error('invalid NELECT or nonstandard DEG_THRESHOLD; require 0.002 eV')
      endif
! Explicit INCAR values must agree with the final effective settings.
      do j=1,nassign
       do i=1,ntag
        if(names(j).ne.tags(i))cycle
        if(i.ge.4.and.i.le.11)then
         if(wd_bool(ivals(j)).neqv.wd_bool(values(i)))then
          call wd_error('INCAR/OUTCAR setting mismatch: '//trim(tags(i)))
         endif
        else
         x=wd_real(ivals(j));old=wd_real(values(i))
         if(i.eq.3)then
          if(x.lt.1d0.or.x.gt.old.or.x.ne.dble(nint(x)))then
           call wd_error('INCAR NBANDS exceeds effective OUTCAR band count')
          endif
         elseif(x.ne.old)then
          call wd_error('INCAR/OUTCAR setting mismatch: '//trim(tags(i)))
         endif
        endif
       enddo
      enddo
      weights=grid(4,:)
      end subroutine wd_validate_run

      subroutine wd_read_incar(path,nassign,names,values)
      implicit none
      character*(*) path
      integer nassign,u,ios,i,j,p,q,start,last
      character*64 names(512),name
      character*256 values(512),value
      character*2048 line,part
      logical duplicate,system_text
      nassign=0
      open(newunit=u,file=trim(path),status='old',action='read',iostat=ios)
      if(ios.ne.0)call wd_error('same-run INCAR is required')
      do
       read(u,'(A)',iostat=ios)line
       if(ios.lt.0)exit
       if(ios.ne.0)call wd_error('INCAR read failed')
       p=index(line,'!');q=index(line,'#')
       if(p.gt.0)line(p:)=' '
       if(q.gt.0)line(q:)=' '
       start=1;last=len_trim(line);system_text=.false.
       do while(start.le.last)
        p=index(line(start:last),';')
        if(p.eq.0)p=last-start+2
        part=adjustl(line(start:start+p-2))
        start=start+p
        if(len_trim(part).eq.0)cycle
        q=index(part,'=')
! SYSTEM is descriptive text and may itself contain semicolons.
! Later actual assignments still parse; other malformed fragments fail.
        if(q.eq.0.and.system_text)cycle
        if(q.lt.2.or.q.eq.len_trim(part).or.q.gt.64.or.
     &    len_trim(part)-q.gt.256)then
         call wd_error('unsupported/ambiguous INCAR assignment')
        endif
        name=trim(part(:q-1));value=trim(adjustl(part(q+1:)))
        call wd_upper(name)
        system_text=name.eq.'SYSTEM'
        duplicate=.false.
        do j=1,nassign
         if(names(j).eq.name)then
          if(values(j).ne.value)call wd_error('conflicting INCAR tag '//trim(name))
          duplicate=.true.
         endif
        enddo
        if(.not.duplicate)then
         nassign=nassign+1
         if(nassign.gt.512)call wd_error('INCAR exceeds supported tag count')
         names(nassign)=name;values(nassign)=value
        endif
       enddo
      enddo
      close(u)
      end subroutine wd_read_incar

      subroutine wd_out_token(line,name,token,found)
      implicit none
      character*(*) line,name,token
      character*2048 part
      integer p,q,j,n
      logical found
      found=.false.;token=' '
      part=adjustl(line)
      if(name.eq.'NKPTS'.or.name.eq.'NBANDS')then
       if(index(part,'k-points ').ne.1)return
       p=index(part,trim(name))
       if(p.eq.0)return
       part=part(p:)
      endif
      n=len_trim(name)
      if(part(:n).ne.trim(name))return
      part=adjustl(part(n+1:))
      if(part(:1).ne.'=')return
      part=adjustl(part(2:));q=len_trim(part)
      do j=1,len_trim(part)
       if(part(j:j).eq.' '.or.part(j:j).eq.';'.or.
     &    iachar(part(j:j)).eq.9)then
        q=j-1;exit
       endif
      enddo
      if(q.lt.1.or.q.gt.len(token))call wd_error('malformed OUTCAR '//trim(name))
      token=part(:q);found=.true.
      end subroutine wd_out_token

      logical function wd_bool(token)
      implicit none
      character*(*) token
      character*256 value
      value=adjustl(token)
      call wd_upper(value)
      select case(trim(value))
      case('T','TRUE','.T.','.TRUE.')
       wd_bool=.true.
      case('F','FALSE','.F.','.FALSE.')
       wd_bool=.false.
      case default
       call wd_error('invalid or unsupported VASP logical '//trim(token))
      end select
      end function wd_bool

      real*8 function wd_real(token)
      use, intrinsic :: ieee_arithmetic
      implicit none
      character*(*) token
      integer ios
      read(token,*,iostat=ios)wd_real
      if(ios.ne.0)call wd_error('invalid VASP number '//trim(token))
      if(.not.ieee_is_finite(wd_real))call wd_error('nonfinite VASP metadata')
      end function wd_real

      subroutine wd_upper(value)
      implicit none
      character*(*) value
      integer i,k
      do i=1,len_trim(value)
       k=iachar(value(i:i))
       if(k.ge.97.and.k.le.122)value(i:i)=achar(k-32)
      enddo
      end subroutine wd_upper

      subroutine wd_error(message)
      implicit none
      character*(*) message
      write(0,'(A)')'*** error - WAVEDER Kubo: '//trim(message)
      call vaspberry_fail
      end subroutine wd_error

! Absolute display paths retain explicit cwd-relative semantics, without chdir.
      subroutine wd_print_inputs(wavepath,derpath,inpath,outpath)
      use, intrinsic :: iso_c_binding, only: c_char,c_ptr,c_size_t,
     & c_associated,c_null_char
      implicit none
! The C ABI avoids incompatible GNU-subroutine/Intel-function GETCWD APIs.
      interface
       function wd_getcwd(buffer,size) bind(C,name='getcwd') result(pointer)
        import c_char,c_ptr,c_size_t
        character(kind=c_char),intent(out) :: buffer(*)
        integer(c_size_t),value :: size
        type(c_ptr) :: pointer
       end function wd_getcwd
      end interface
      character*(*) wavepath,derpath,inpath,outpath
      character*4096 cwd
      character(kind=c_char) cwd_buffer(4096)
      type(c_ptr) cwd_pointer
      character*256 paths(4)
      character*7 names(4)
      integer i
      paths(1)=wavepath;paths(2)=derpath
      paths(3)=inpath;paths(4)=outpath
      names=(/'WAVECAR','WAVEDER','INCAR  ','OUTCAR '/)
      cwd_pointer=wd_getcwd(cwd_buffer,int(size(cwd_buffer),c_size_t))
      if(.not.c_associated(cwd_pointer))then
       call wd_error('cannot determine invocation working directory')
      endif
      cwd=''
      do i=1,size(cwd_buffer)
       if(cwd_buffer(i).eq.c_null_char)exit
       cwd(i:i)=cwd_buffer(i)
      enddo
      write(6,'(A,A)')'# Input path base (invocation cwd): ',trim(cwd)
      do i=1,4
       if(paths(i)(1:1).eq.'/')then
        write(6,'(A,A,A,A)')'# Input ',trim(names(i)),': ',trim(paths(i))
       else
        write(6,'(A,A,A,A)')'# Input ',trim(names(i)),': ',
     &    trim(cwd)//'/'//trim(paths(i))
       endif
      enddo
      flush(6)
      end subroutine wd_print_inputs
