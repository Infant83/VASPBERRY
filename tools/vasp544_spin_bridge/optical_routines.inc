! Original VASPBERRY optical stream instrumentation; no VASP implementation body.
  SUBROUTINE BERRY_INPROD_DP(W1,W0,CHAM)
      USE wave_high
      USE dfast
      IMPLICIT NONE
      TYPE(wavespin) :: W1,W0
      GDEF :: CHAM(:,:,:,:)
      TYPE(wavedes1) :: WD0,WD1
      TYPE(wavefuna) :: WA
      GDEF, ALLOCATABLE :: BLOCK(:,:)
      INTEGER ISP,NK,NPOS,NS,NA,NB0,NB1

      NB0=W0%WDES%NB_TOT
      NB1=W1%WDES%NB_TOT
      IF (W0%WDES%DO_REDIS.OR.W1%WDES%DO_REDIS) THEN
         WRITE(0,*) 'BERRY_INPROD_DP: redistribution is not supported in serial v1'
         STOP 84
      ENDIF
      IF (SIZE(CHAM,1)/=NB0.OR.SIZE(CHAM,2)/=NB1.OR.NB0/=NB1) STOP 85
      NS=MIN(NSTRIP_STANDARD_GLOBAL,NB1)
      ALLOCATE(BLOCK(NB0,NS))
      DO ISP=1,W0%WDES%ISPIN
         DO NK=1,W0%WDES%NKPTS
            CALL SETWDES(W0%WDES,WD0,NK)
            CALL SETWDES(W1%WDES,WD1,NK)
            WA=ELEMENTS(W1,WD1,ISP)
            DO NPOS=1,NB1,NS
               NA=MIN(NS,NB1+1-NPOS)
               BLOCK=0
               CALL ORTH2_NOSUBINDEX( &
                    W0%CW(1,1,NK,ISP),WA%CW_RED(1,NPOS),W0%CPROJ(1,1,NK,ISP), &
                    WA%CPROJ_RED(1,NPOS),NB0,NPOS,NA,WD0%NPL_RED,WD0%NPRO_O_RED, &
                    WD0%NRPLWV_RED,WD0%NPROD_RED,BLOCK(1,1))
               CHAM(:,NPOS:NPOS+NA-1,NK,ISP)=BLOCK(:,1:NA)
            ENDDO
         ENDDO
      ENDDO
! No conjugate triangle fill or Hermitian averaging: retain computed evidence.
      DEALLOCATE(BLOCK)
  END SUBROUTINE BERRY_INPROD_DP

  SUBROUTINE BERRY_WRITE_CONNECTION(WDES,W,LATT_CUR,EFERMI)
      USE wave
      USE lattice
      USE constant
      IMPLICIT NONE
      TYPE(wavedes) :: WDES
      TYPE(wavespin) :: W
      TYPE(latt) :: LATT_CUR
      REAL(q) :: EFERMI
      INTEGER :: IU,IERR
      INTEGER(KIND=4) :: HEADER(7)
      CHARACTER(LEN=32) :: MAGIC,FOOTER
      MAGIC='VBERRY_PAW_CONNECTION_V1'
      FOOTER='VBERRY_CONNECTION_COMPLETE_V1'
      HEADER=(/1,WDES%NKPTS,WDES%NB_TOT,WDES%ISPIN,3,8,16/)
      OPEN(NEWUNIT=IU,FILE='BERRY_CONNECTION.bin',ACCESS='STREAM',FORM='UNFORMATTED', &
           STATUS='REPLACE',ACTION='WRITE',IOSTAT=IERR)
      IF (IERR/=0) STOP 86
! Native little-endian stream for this tested Intel build; reader verifies magic,
! scalar precision, dimensions, exact byte length and complete footer.
! Fortran axes: band_bra,band_ket,k,spin,Cartesian. No optical/occupation weights.
      WRITE(IU,IOSTAT=IERR) MAGIC,HEADER,BERRY_DEG_THRESHOLD,EFERMI, &
           LATT_CUR%A,LATT_CUR%B*TPI,WDES%VKPT,WDES%WTKPT, &
           REAL(W%CELTOT,q),W%FERTOT,ENERGY_DER,BERRY_CDER_DP,FOOTER
      IF (IERR/=0) STOP 87
      CLOSE(IU,IOSTAT=IERR)
      IF (IERR/=0) STOP 88
      WRITE(*,*) 'VASPBERRY: completed BERRY_CONNECTION.bin (experimental PAW optical connection)'
  END SUBROUTINE BERRY_WRITE_CONNECTION
