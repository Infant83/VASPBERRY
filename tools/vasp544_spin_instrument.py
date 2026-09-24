"""Add original, opt-in VASPBERRY exporters to an isolated licensed VASP copy.

No VASP implementation is distributed here. Two exact source checksums and
unique structural insertion points limit compatibility to the audited 5.4.4
sources. Existing VASP calculations are unchanged unless export is enabled.
Use a serial complex ncl build. Never apply this to your only working copy.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re

BASE_HASHES = {
    'linear_optics.F': '1eac9ddc5ad014572f797675c6a25c50d01aa234d5b6035562783626e6302cd5',
    'elinear_response.F': '93aa633c326e2b75c5de4eae5946c3aa799e945d315265fadedf0eea8d009434',
}
TEMPLATES = Path(__file__).with_name('vasp544_spin_bridge')


def sha(value):
    return hashlib.sha256(value).hexdigest()


def one(text, pattern, replacement):
    out, count = re.subn(pattern, replacement, text, flags=re.M)
    if count != 1:
        raise ValueError(f'expected one compatible insertion point, found {count}: {pattern}')
    return out


def in_routine(text, name, transform):
    match = re.search(r'(?mi)^\s*SUBROUTINE\s+'+name+r'\s*\(', text)
    if not match:
        raise ValueError('missing routine '+name)
    end = re.search(r'(?mi)^\s*END\s+SUBROUTINE\s+'+name+r'\b[^\n]*', text[match.start():])
    if not end:
        raise ValueError('missing routine end '+name)
    stop = match.start()+end.end()
    return text[:match.start()]+transform(text[match.start():stop])+text[stop:]


def after_call(text, name, addition):
    matches = list(re.finditer(r'(?mi)^\s*CALL\s+'+name+r'\s*\(', text))
    if len(matches) != 1:
        raise ValueError('ambiguous call '+name)
    pos = matches[0].end(); depth = 1
    while pos < len(text) and depth:
        depth += (text[pos] == '(')-(text[pos] == ')'); pos += 1
    if depth:
        raise ValueError('unterminated call '+name)
    return text[:pos]+'\n'+addition+text[pos:]


def patched_sources(original):
    spin = (TEMPLATES/'spin_routines.inc').read_text()
    optical = (TEMPLATES/'optical_routines.inc').read_text()
    response = original['elinear_response.F']
    response = one(response, r'(^  IMPLICIT NONE\s*$)', lambda m:m[1]+'''
  LOGICAL, SAVE :: LSPIN_EXPORT=.FALSE.
  COMPLEX(q), ALLOCATABLE, SAVE :: BSP_PSEUDO(:,:,:,:,:),BSP_DELTA(:,:,:,:,:)
  COMPLEX(q), ALLOCATABLE, SAVE :: BSP_XI(:,:,:,:,:),BSP_T(:,:,:,:,:),BSP_V(:,:,:,:,:)
''')
    response = one(response, r'(^CONTAINS\s*$)', lambda m:m[1]+'\n'+spin)

    def response_hooks(body):
        body = after_call(body, 'MRG_CEL', '    IF (LSPIN_EXPORT) CALL BSP_RAW(W0,WXI,IDIR)')
        body = one(body, r'(^[ \t]*CALL\s+OVERL_ALL\([^\n]*CDIJ1\)[ \t]*$)',
                   lambda m:m[1]+'\n       IF (LSPIN_EXPORT) CALL BSP_CORRECTION(W0,WXI,LMDIM,CQIJ,CDIJ1,IDIR)')
        return body
    response = in_routine(response, 'LRF_RPHI0', response_hooks)
    optic = original['linear_optics.F']
    optic = one(optic, r'(^  GDEFS,[^\n]*:: CDER_BETWEEN_STATES\([^\n]*$)', lambda m:m[1]+'''
  LOGICAL, SAVE :: LBERRY_EXPORT=.FALSE.
  REAL(q), PARAMETER :: BERRY_DEG_THRESHOLD=1E-10_q
  GDEF, ALLOCATABLE, SAVE :: BERRY_CDER_DP(:,:,:,:,:)
''')

    def optical_hooks(body):
        body = one(body, r'(^    USE pead\s*$)', lambda m:m[1]+'\n    USE fock, ONLY : LHFCALC')
        body = one(body, r'(^    INTEGER IDIR,[^\n]*ISTEP\s*$)', lambda m:m[1]+'''
    INTEGER BERRY_IDUM, BERRY_N, BERRY_IERR, BERRY_IU
    REAL(q) BERRY_RDUM
    COMPLEX(q) BERRY_CDUM
    CHARACTER(40) BERRY_CHAR
''')
        body = one(body, r'(^    EMAX=MAX_ENERGY_UNOCCUPIED[^\n]*$)', lambda m:'''    LBERRY_EXPORT=.FALSE.
    BERRY_IU=-1
    CALL RDATAB(.TRUE.,'INCAR',BERRY_IU,'LBERRY_EXPORT','=','#',';','L', &
         BERRY_IDUM,BERRY_RDUM,BERRY_CDUM,LBERRY_EXPORT,BERRY_CHAR,BERRY_N,1,BERRY_IERR)
    IF ((BERRY_IERR/=0.AND.BERRY_IERR/=3).OR.(BERRY_IERR==0.AND.BERRY_N<1)) STOP 'Invalid LBERRY_EXPORT'
    IF (LBERRY_EXPORT) THEN
#if defined(MPI) || defined(gammareal)
       STOP 'VASPBERRY exporter requires serial complex ncl build'
#endif
       IF (INFO%LREAL.OR.LUSEPEAD().OR.LNABLA.OR.LHFCALC.OR.LDO_METAGGA().OR.SYMM%ISYM/=-1) &
          STOP 'Unsupported VASPBERRY PAW optical branch'
    ENDIF
    LSPIN_EXPORT=.FALSE.
    CALL RDATAB(.TRUE.,'INCAR',BERRY_IU,'LSPIN_EXPORT','=','#',';','L', &
         BERRY_IDUM,BERRY_RDUM,BERRY_CDUM,LSPIN_EXPORT,BERRY_CHAR,BERRY_N,1,BERRY_IERR)
    IF ((BERRY_IERR/=0.AND.BERRY_IERR/=3).OR.(BERRY_IERR==0.AND.BERRY_N<1)) STOP 'Invalid LSPIN_EXPORT'
    IF (LSPIN_EXPORT) THEN
      IF (.NOT.LBERRY_EXPORT) STOP 'LSPIN_EXPORT requires LBERRY_EXPORT'
      CALL BSP_INIT(W,LMDIM,CQIJ)
    ENDIF
'''+m[1])
        body = one(body, r'(^[ \t]*CALL FIND_DEG_CLUSTERS\(([^\n]+)\)[ \t]*$)',
                   lambda m:'    IF (LBERRY_EXPORT) THEN\n       CALL FIND_DEG_CLUSTERS('+m[2]+', THRESHOLD=BERRY_DEG_THRESHOLD)\n    ELSE\n'+m[1]+'\n    ENDIF')
        body = one(body, r'(^[ \t]*IF\s*\(LVEL\)\s*THEN[ \t]*\n)(?=[ \t]*NBANDS_CDER[ \t]*=[ \t]*WDES%NB_TOT)',
                   '    IF (LVEL.OR.LBERRY_EXPORT) THEN\n')
        body = one(body, r'(^    ENERGY_DER=0[ \t]*$)', lambda m:m[1]+'''
    IF (LBERRY_EXPORT) THEN
       IF (ALLOCATED(BERRY_CDER_DP)) DEALLOCATE(BERRY_CDER_DP)
       ALLOCATE(BERRY_CDER_DP(WDES%NB_TOT,NBANDS_CDER,WDES%NKPTS,WDES%ISPIN,3))
       BERRY_CDER_DP=0
    ENDIF
''')
        body = after_call(body, 'LRF_EPSILON', '''    IF (LBERRY_EXPORT) THEN
       CALL BERRY_WRITE_CONNECTION(WDES,W,LATT_CUR,EFERMI)
       IF (LSPIN_EXPORT) CALL BSP_WRITE(W,LATT_CUR,EFERMI)
       DEALLOCATE(BERRY_CDER_DP)
    ENDIF
''')
        return body
    optic = in_routine(optic, 'LR_OPTIC', optical_hooks)
    optic = in_routine(optic, 'LRF_EPSILON', lambda body:one(body, r'(^[ \t]*CALL INPROD_W\([^\n]*$)',
                          lambda m:'       IF (LBERRY_EXPORT) CALL BERRY_INPROD_DP(W1,W0,BERRY_CDER_DP(:,:,:,:,IDIR))\n'+m[1]))
    optic = one(optic, r'(^  SUBROUTINE INPROD_W\()', lambda m:optical+'\n\n'+m[1])
    return {'elinear_response.F': response, 'linear_optics.F': optic}


def instrument(directory):
    source = Path(directory)/'src'
    manifest = Path(directory)/'vaspberry-spin-producer.json'
    if manifest.exists():
        raise ValueError('producer manifest already exists; use a fresh isolated copy')
    original = {}
    for name, digest in BASE_HASHES.items():
        raw = (source/name).read_bytes()
        if sha(raw) != digest:
            raise ValueError('unsupported or already modified VASP source: '+name)
        original[name] = raw.decode()
    modified = patched_sources(original)  # Validate all anchors before any write.
    records = {n:{'before_sha256': BASE_HASHES[n], 'after_sha256': sha(t.encode())} for n,t in modified.items()}
    for name, text in modified.items():
        (source/name).write_text(text)
    record = {'schema': 'vaspberry.vasp544-spin-producer', 'version': 1,
              'source_files': records, 'instrumenter_sha256': sha(Path(__file__).read_bytes()),
              'template_sha256': {n:sha((TEMPLATES/n).read_bytes()) for n in ('spin_routines.inc', 'optical_routines.inc')},
              'supported_scope': 'Serial complex VASP5.4.4 ncl; audited source hashes; PAW spin-independent overlap; Cartesian +z SAXIS.',
              'new_files_are_original_instrumentation': True}
    manifest.write_text(json.dumps(record, indent=2)+'\n')
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('isolated_vasp_directory', type=Path)
    print(json.dumps(instrument(parser.parse_args().isolated_vasp_directory), indent=2))
