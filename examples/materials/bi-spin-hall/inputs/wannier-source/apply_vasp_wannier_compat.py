"""Apply the verified VASP 5.4.4 Wannier spin-column edits to a licensed source.

Only original insertion text and short API-call replacements are distributed.
The licensed source module is supplied by the user and is never bundled.
This program refuses other source revisions and never overwrites a file.
"""
from pathlib import Path
import argparse,hashlib,json
BEFORE = '73fca10313523732609932f3d8cd5f5efaf8e90d3399728ce66675a99667b234'
AFTER = '30cf701fa9b7f29f515e905a1fb3255531f894a4b5f489ea2aa0ff8ec701b4f4'
EDITS = [[336, 336, '      LOGICAL, PRIVATE, SAVE :: REUSE_MMN=.FALSE.,CACHE_WANPROJ=.FALSE.\n'], [10174, 10174, "         REUSE_MMN=.FALSE.\n         CALL RDATAB(.FALSE.,INCAR,IU5,'LREUSE_MMN','=','#',';','L', &\n        &            IDUM,RDUM,CDUM,REUSE_MMN,CHARAC,N,1,IERR)\n         IF (IERR/=0.AND.IERR/=3) STOP 'LREUSE_MMN parse error'\n         IF (REUSE_MMN.AND.LWANNIER90_RUN) STOP 'LREUSE_MMN requires export-only mode'\n         CACHE_WANPROJ=.FALSE.\n         CALL RDATAB(.FALSE.,INCAR,IU5,'LCACHE_WANPROJ','=','#',';','L', &\n        &            IDUM,RDUM,CDUM,CACHE_WANPROJ,CHARAC,N,1,IERR)\n         IF (IERR/=0.AND.IERR/=3) STOP 'LCACHE_WANPROJ parse error'\n\n"], [12813, 12813, '      INTEGER, ALLOCATABLE :: proj_s(:),proj_cache_key(:)\n      REAL(q), ALLOCATABLE :: proj_s_qaxis(:,:)\n      GDEF, ALLOCATABLE :: A_cache(:,:,:,:)\n      LOGICAL, ALLOCATABLE :: A_cache_valid(:,:)\n      INTEGER :: NACTIVE,NKEYS,JCNTR,IKEY\n      LOGICAL :: LBUILD_RADIAL,LMATCH\n\n'], [17438, 17438, '      ALLOCATE(proj_s(num_bands_tot),proj_s_qaxis(3,num_bands_tot),proj_cache_key(num_bands_tot))\n      proj_s=0;proj_s_qaxis=0;proj_cache_key=0\n'], [17923, 17963, '     &                   exclude_bands,proj_s,proj_s_qaxis)\n'], [18770, 18770, "      CALLMPI( M_sum_i(WDES%COMM,proj_s,num_bands_tot) )\n      CALLMPI( M_sum_d(WDES%COMM,proj_s_qaxis,3*num_bands_tot) )\n      NACTIVE=0;NKEYS=0\n      IF (spinors) THEN\n         IF (MAXVAL(ABS(WDES%SAXIS-(/0._q,0._q,1._q/)))>1.E-10_q) &\n        &   STOP 'Private spinor compatibility export supports SAXIS=+z only'\n      ENDIF\n      DO ICNTR=1,num_bands_tot\n         IF (proj_l(ICNTR)==0.AND.proj_m(ICNTR)==0.AND.proj_radial(ICNTR)==0) CYCLE\n         NACTIVE=NACTIVE+1\n         IF (spinors) THEN\n            IF (ABS(proj_s(ICNTR))/=1) STOP 'Invalid Wannier spin label'\n            IF (MAXVAL(ABS(proj_s_qaxis(:,ICNTR)-(/0._q,0._q,1._q/)))>1.E-10_q) &\n           &   STOP 'Private spinor compatibility export supports projection spin axis +z only'\n         ELSE\n            proj_s(ICNTR)=1\n         ENDIF\n         DO JCNTR=1,ICNTR-1\n            IF (proj_cache_key(JCNTR)==0) CYCLE\n            LMATCH=ALL(proj_site(:,ICNTR)==proj_site(:,JCNTR)).AND. &\n           & proj_radial(ICNTR)==proj_radial(JCNTR).AND.proj_zona(ICNTR)==proj_zona(JCNTR).AND. &\n           & proj_s(ICNTR)==proj_s(JCNTR)\n            IF (LMATCH) THEN\n               proj_cache_key(ICNTR)=proj_cache_key(JCNTR)\n               EXIT\n            ENDIF\n         ENDDO\n         IF (proj_cache_key(ICNTR)==0) THEN\n            NKEYS=NKEYS+1;proj_cache_key(ICNTR)=NKEYS\n         ENDIF\n      ENDDO\n      IF (NACTIVE/=num_wann) STOP 'Active spin-labelled projections must equal num_wann'\n      IF (IO%IU0>=0) WRITE(IO%IU0,*) 'Audited spin-labelled projections/cache groups:',NACTIVE,NKEYS\n\n"], [20598, 20598, "      IF (REUSE_MMN) THEN\n         INQUIRE(FILE=seed_name//'.mmn',EXIST=LDUM)\n         IF (.NOT.LDUM.OR.WDES%ISPIN/=1) STOP 'LREUSE_MMN requires an existing audited spinor mmn'\n         IF (IO%IU0>=0) WRITE(IO%IU0,*) 'Reusing closed externally audited mmn; export-only stage'\n      ELSE\n\n"], [22296, 22296, '      ENDIF ! REUSE_MMN\n'], [22414, 22448, ''], [22586, 22598, ''], [24745, 24760, '      NPROJ=0\n      IF (CACHE_WANPROJ) THEN\n         ALLOCATE(A_cache(WDES%NB_TOT,(LMAX+1)**2,num_kpts,NKEYS),A_cache_valid(num_kpts,NKEYS))\n         A_cache_valid=.FALSE.\n      ENDIF\n'], [24828, 24870, ''], [25022, 25022, "         IF (NPROJ>SIZE(A_matrix,2)) STOP 'Projection column exceeds A_matrix allocation'\n         ISPINOR=1\n         IF (spinors.AND.proj_s(ICNTR)==-1) ISPINOR=2\n         IKEY=proj_cache_key(ICNTR)\n         LBUILD_RADIAL=.TRUE.\n         IF (CACHE_WANPROJ) LBUILD_RADIAL=.NOT.ALL(A_cache_valid(:,IKEY))\n"], [25642, 25642, '         IF (LBUILD_RADIAL) THEN\n'], [26133, 26133, '         ENDIF ! LBUILD_RADIAL\n'], [26239, 26373, '            IF (CACHE_WANPROJ) THEN\n               IF (A_cache_valid(NKI,IKEY)) THEN\n                  A=A_cache(:,:,NKI,IKEY)\n               ELSE\n                  CALL CALC_OVERLAP_GN( &\n                 & LMAX,FG,proj_site(:,ICNTR),W,kpt_latt(:,NKI),ISP,ISPINOR,P,CQIJ,LATT_CUR,T_INFO,A)\n                  A_cache(:,:,NKI,IKEY)=A;A_cache_valid(NKI,IKEY)=.TRUE.\n               ENDIF\n            ELSE\n               CALL CALC_OVERLAP_GN( &\n              & LMAX,FG,proj_site(:,ICNTR),W,kpt_latt(:,NKI),ISP,ISPINOR,P,CQIJ,LATT_CUR,T_INFO,A)\n            ENDIF\n'], [27469, 27507, '         IF (LBUILD_RADIAL) DEALLOCATE(FTMP,FR,FGTMP,FG)\n'], [27700, 27719, '      IF (CACHE_WANPROJ) DEALLOCATE(A_cache,A_cache_valid)\n']]

def digest(data):return hashlib.sha256(data).hexdigest()
def transform(data):
    if digest(data)!=BEFORE:
        raise ValueError('Unsupported mlwf.F revision or already modified source; original VASP 5.4.4 input required')
    previous=len(data)
    for start,end,text in reversed(EDITS):
        if not (0<=start<=end<=previous):raise ValueError('Invalid edit bounds')
        data=data[:start]+text.encode('utf-8')+data[end:];previous=start
    if digest(data)!=AFTER:raise ValueError('Modified source identity check failed')
    return data

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True,help='User-licensed original src/mlwf.F')
    p.add_argument('--output',type=Path,required=True,help='New output file; never overwrites the source')
    a=p.parse_args()
    if a.output.exists() or a.output.with_name(a.output.name+'.json').exists():
        raise ValueError('Output or its audit file already exists')
    original=a.source.read_bytes();modified=transform(original)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('xb') as f:f.write(modified)
    record=dict(status='PASS_SOURCE_TRANSFORM',source_sha256=BEFORE,output_sha256=AFTER,
      edits=len(EDITS),source_bytes=len(original),output_bytes=len(modified),
      scope='VASP5.4.4 v2 spin-labelled projection mapping; serial +z spin axes; exact optional raw overlap cache. Compilation and actual export checks remain required.')
    a.output.with_name(a.output.name+'.json').write_text(json.dumps(record,indent=2)+'\n')
    print(a.output)
if __name__=='__main__':main()
