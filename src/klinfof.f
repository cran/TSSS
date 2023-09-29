C     PROGRAM 4.2  KLINFO
      SUBROUTINE KLINFOF( IDISTG,PG,IDISTF,PF,XMIN,XMAX,NINT,DX,
     *                   FKLI,GINT )
C
      INCLUDE 'TSSS.h'
C
C  ...  Driver program of the subroutine KLINFO  ...
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  PARAMG(3), PARAMF(3)
cxx      DIMENSION PG(2), PF(2)
cxx      DIMENSION NINT(4), DX(4), FKLI(4), GINT(4)
C
      INTEGER IDISTG, IDISTF, NINT(4)
      DOUBLE PRECISION PG(2), PF(2), XMIN, XMAX, DX(4), FKLI(4),
     1                 GINT(4), GAUSS, CAUCHY
c local
      INTEGER II
C
      EXTERNAL  GAUSS
      EXTERNAL  CAUCHY
cc      DATA  PARAMG /0.0D0, 1.0D0, 0.0D0/
cc      DATA  PARAMF /0.1D0, 1.5D0, 0.0D0/
C
cc      XMIN = -8.0D0
cc      XMAX =  8.0D0
cc      WRITE(6,600)  XMIN, XMAX
cc      WRITE(6,610)
      DO 10 II=1,4
cc      NINT = (XMAX-XMIN+1.0D-5)*2**(II-1)
cc   10 CALL KLINFO(GAUSS,GAUSS ,PAG,PAF,XMIN,XMAX,NINT,FKLI,GINT)
cxx         NINT(II) = (XMAX-XMIN+1.0D-5)*2**(II-1)
         NINT(II) = INT((XMAX-XMIN+1.0D-5)*2**(II-1))
         IF( IDISTG .EQ. 1 ) THEN
            IF( IDISTF .EQ. 1 ) THEN
               CALL KLINFO( GAUSS,GAUSS ,PG,PF,XMIN,XMAX,
     *                            NINT(II),DX(II),FKLI(II),GINT(II) )
            ELSE
               CALL KLINFO( GAUSS, CAUCHY,PG,PF,XMIN,XMAX,
     *                            NINT(II),DX(II),FKLI(II),GINT(II) )
            END IF
         ELSE
            IF( IDISTF .EQ. 1 ) THEN
               CALL KLINFO( CAUCHY,GAUSS ,PG,PF,XMIN,XMAX,
     *                            NINT(II),DX(II),FKLI(II),GINT(II) )
            ELSE
               CALL KLINFO( CAUCHY,CAUCHY,PG,PF,XMIN,XMAX,
     *                            NINT(II),DX(II),FKLI(II),GINT(II) )
            END IF
         END IF
   10 CONTINUE
cc      STOP
      RETURN
cxx  600 FORMAT( 1H ,'XMIN =',F5.1,3X,'XMAX =',F5.1 )
cxx  610 FORMAT( 1H ,'  XMIN  XMAX   NINT',3X,'DX',9X,'FKLI',11X,'GINT' )
      E N D
      SUBROUTINE  KLINFO( DISTG,DISTF,PARAMG,PARAMF,XMIN,XMAX,NINT,
cc     *                    FKLI,GINT )
     *                     DX,FKLI,GINT )  
C
C  ...  This subroutine computes Kullback-Leibler information  ...
C
C     Inputs:
C        DISTG:    function name for the true density
C        DISTF:    function name for the model density
C        PARAMG:   parameter vector of the true density
C        PARAMF:   parameter vector of the model density
C        XMIN:     lower limit of integration
C        XMAX:     upper limit of integration
C        NINT:     number of function evaluation
C     Outputs:
C        FKLI:     Kullback-Leibler information number, I(g;f)
C        GINT:     integration of g(y) over [XMIN,XMAX]
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cx      DIMENSION  PARAMG(*), PARAMF(*)
cxx      DIMENSION  PARAMG(2), PARAMF(2)
C
      INTEGER NINT
      DOUBLE PRECISION DISTG, DISTF, PARAMG(2), PARAMF(2), XMIN, XMAX,
     1                 DX, FKLI, GINT
c local 
      INTEGER I
      DOUBLE PRECISION XX, GX, FX 
C
cxxx      DX = (XMAX-XMIN)/DFLOAT(NINT)
      DX = (XMAX-XMIN)/DBLE(NINT)
      FKLI = 0.0D0
      GINT = 0.0D0
C
      DO 10  I=0,NINT
      XX = XMIN + DX*I
      GX = DISTG( XX,PARAMG )
      FX = DISTF( XX,PARAMF )
      IF( I.GE.1 .AND. I.LT.NINT )  THEN
         FKLI = FKLI + DLOG( GX/FX )*GX
         GINT = GINT + GX
      ELSE
         FKLI = FKLI + DLOG( GX/FX )*GX/2.0D0
         GINT = GINT + GX/2.0D0
      END IF
   10 CONTINUE
C
      FKLI = FKLI*DX
      GINT = GINT*DX
cc      WRITE(6,600) XMIN, XMAX, NINT, DX, FKLI, GINT
C
      RETURN
cxx  600 FORMAT( 1H ,2F6.2,I5,F9.4,D17.8,F12.8 )
      E N D
