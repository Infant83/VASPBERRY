program test_kubo_parser
  implicit none
  character(256) filename,foname,fbz,ver_tag,klist_fname,sw_fname,atlist_fname,kubo_csv
  character(256) kubo_pairs
  integer nkx,nky,ispinor,icd,ixt,ivel,ikubo,iz,ihf,nini,nmax,nn,kperiod,it,iskp,ine
  integer iwf,ikwf,ng(3),imag,nediv,ikubo_bundle
  real(8) rs(3),init_e,fina_e,sigma,theta,phi
  logical flag_atom_project
  nkx=2;nky=2;ispinor=2;kperiod=1;ver_tag='parser regression'
  call parse(filename,foname,nkx,nky,ispinor,icd,ixt,fbz, &
             ivel,ikubo,iz,ihf,nini,nmax,nn,kperiod,it,iskp,ine,ver_tag, &
             iwf,ikwf,ng,rs,imag,init_e,fina_e,nediv,sigma, &
             klist_fname,sw_fname,atlist_fname,flag_atom_project,theta,phi,kubo_csv,ikubo_bundle,kubo_pairs)
  write(*,'(A)')trim(foname)
end program

subroutine help(ver_tag)
  character(*) ver_tag
  stop 1
end subroutine

subroutine vaspberry_fail
  stop 2
end subroutine
