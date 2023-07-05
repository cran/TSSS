C     PROGRAM 6.2  MARSPC
      SUBROUTINE MARSPCF( M,L,A,E,NF,P,AMP,ANG,COH,FNC,FRNC )
C
      INCLUDE 'TSSS.h'
C
C  ...  This program computes cross spectra and noise contribution  ...
C
C     Inputs:
C        IDEV:    Input device
C        M:       AR order
C        L:       Dimension of the time series
C        E(I,J):  Innovation covariance matrix
C        A(II,I,J):  AR coefficient matrices
C     Parametrs:
C        NF:      Number of frequencies
C        MJ:      Adjustable dimension
C        MJ1:     Adjustable dimension
C     Written by G.K.
cc      PARAMETER( NF=200,MJ=3,MJ1=7,IDEV=1 )
cxx      IMPLICIT  REAL*8(A-H)
cxx      IMPLICIT  COMPLEX*16(O-Z)
cc      DIMENSION  A(MJ,MJ,MJ1), E(MJ,MJ)
cc      DIMENSION  P(0:NF,MJ,MJ), FNC(0:NF,MJ,MJ), COH(0:NF,MJ,MJ)
cc      DIMENSION  AMP(0:NF,MJ,MJ), ANG(0:NF,MJ,MJ)
cxx      DIMENSION  A(L,L,M), E(L,L), FRNC(0:NF,L,L)
cxx      DIMENSION  P(0:NF,L,L), FNC(0:NF,L,L), COH(0:NF,L,L)
cxx      DIMENSION  AMP(0:NF,L,L), ANG(0:NF,L,L)
C
      INTEGER M, L, NF
      DOUBLE PRECISION A(L,L,M), E(L,L), AMP(0:NF,L,L), ANG(0:NF,L,L),
     1                 COH(0:NF,L,L), FNC(0:NF,L,L), FRNC(0:NF,L,L)
c local
      COMPLEX(KIND(0d0)) P(0:NF,L,L)
c
      AMP(:,:,:) = 0.0d0
      ANG(:,:,:) = 0.0d0
      COH(:,:,:) = 0.0d0
C
cc      open( IDEV,file='ar2.dat' )
cc      READ(IDEV,*)  M, L
cc      DO 10 I=1,L 
cc   10 READ(IDEV,*)  (E(I,J),J=1,L)
cc      DO 20 II=1,M
cc      DO 20  I=1,L
cc   20 READ(IDEV,*)  (A(I,J,II),J=1,L)
cc      close( IDEV )
C
C  ...  Cross spectrum  ...
C
cc      CALL  MARSPC( M,L,A,E,NF,MJ,MJ1,P,FNC,AMP,ANG,COH )
      CALL  MARSPC( M,L,A,E,NF,P,FNC,AMP,ANG,COH )
C
C  ...  Print out plot cross spectra, coherency and noise
C       contribution  ...
C
cc      CALL  PRMSPC( M,L,NF,MJ,P,FNC,COH )
cxx      CALL  PRMSPC( M,L,NF,P,FNC,COH,FRNC )
      CALL  PRMSPC( L,NF,FNC,FRNC )
C      CALL  PTCSPC( M,L,NF,MJ,P,FNC,AMP,ANG,COH )
C
cc      STOP
      RETURN
      E N D
cc      SUBROUTINE  PRMSPC( M,L,NF,MJ,P,FNC,COH )
cxx      SUBROUTINE  PRMSPC( M,L,NF,P,FNC,COH,FRNC )
      SUBROUTINE  PRMSPC( L,NF,FNC,FRNC )
C
C  ...  Print out cross spectra and relative noise contribution  ...
C
C     Inputs:
C        M:    AR order
C        L:    Dimension of the time series
C        NF:   Number of frequencies
C        MJ:   Adjustable dimension of P and FNC
C        P(II,I,J):   Cross spectrum of Y(I) and Y(J) at II-th frequency
C        FNC(II,I,J): Noise contribution of Y(J) to Y(I)
C        COH(II,I,J): Simple coherency
C
cxx      IMPLICIT COMPLEX*16(O-Z)
cxx      IMPLICIT REAL*8(A-H)
cc      DIMENSION  P(0:NF,MJ,MJ), FNC(0:NF,MJ,MJ), FRNC(0:200)
cc      DIMENSION  COH(0:NF,MJ,MJ)
cxx      DIMENSION  P(0:NF,L,L), FNC(0:NF,L,L), FRNC(0:NF,L,L)
cxx      DIMENSION  COH(0:NF,L,L)
C
      INTEGER L, NF
      DOUBLE PRECISION FNC(0:NF,L,L), FRNC(0:NF,L,L)
C
cc      WRITE(6,600)
cc      WRITE(6,610)  M, L, NF
cc      WRITE(6,620)
cc      DO 10 I=1,L
cc      DO 10 J=1,L
cc      WRITE(6,630)  I, J
cc   10 WRITE(6,640)  (P(II,I,J),II=0,NF)
cc      WRITE(6,650)
cc      DO 20 I=1,L-1
cc      DO 20 J=I+1,L
cc      WRITE(6,630)  I, J
cc   20 WRITE(6,660)  (COH(II,I,J),II=0,NF)
cc      WRITE(6,670)
cxx      DO 40 I=1,L
      DO 41 I=1,L
      DO 40 J=1,L
cc      WRITE(6,630)  I, J
      DO 30 II=0,NF
cc      IF(J.EQ.1)  FRNC(II) = FNC(II,I,J)/FNC(II,I,L)
cc   30 IF(J.NE.1)  FRNC(II) = (FNC(II,I,J)-FNC(II,I,J-1))/FNC(II,I,L)
      IF(J.EQ.1)  FRNC(II,I,J) = FNC(II,I,J)/FNC(II,I,L)
cxx   30 IF(J.NE.1)  FRNC(II,I,J) = (FNC(II,I,J)-FNC(II,I,J-1))/FNC(II,I,L)
      IF(J.NE.1)  FRNC(II,I,J) = (FNC(II,I,J)-FNC(II,I,J-1))/FNC(II,I,L)
   30 CONTINUE
cc   40 WRITE(6,660)  (FRNC(II),II=0,NF)
   40 CONTINUE
   41 CONTINUE
C
cxx  600 FORMAT( 1H ,'PROGRAM 6.2:  SPECTRUM ANALYSIS BY MAR MODEL' )
cxx  610 FORMAT( 1H ,'M(ORDER) =',I3,5X,'L(DIMENSION) =',I2,5X,
cxx     *            'NF(NUMBER OF FREQUENCIES =',I3 )
cxx  620 FORMAT( 1H ,'CROSS SPECTRUM' )
cxx  630 FORMAT( 1H ,'*** I =',I2,2X,'J =',I2,' ***' )
cxx  640 FORMAT( 2(' (',D13.5,',',D13.5,')',4X) )
cxx  650 FORMAT( 1H ,'COHERENCY' )
cxx  660 FORMAT( 1H ,10F7.4 )
cxx  670 FORMAT( 1H ,'RELATIVE NOISE CONTRIBUTION' )
      RETURN
      E N D
cc      SUBROUTINE  MARSPC( M,L,A,E,NF,MJ,MJ1,P,FNC,AMP,ANG,COH )
      SUBROUTINE  MARSPC( M,L,A,E,NF,P,FNC,AMP,ANG,COH )
C
C  ...  Cross spectra, amplitude, phase and noise contribution  ...
C
C     Inputs:
C        M:       AR order
C        L:       Dimension of the time series
C        E(I,J):  Innovation covariance matrix
C        A(I,J,II):  AR coefficient matrices
C        NF:      Number of frequencies
C        MJ:      Adjustable dimension
C        MJ1:     Adjustable dimension
C     Outputs:
C        P(II,I,J):    Cross spectra
C        FNC(II,I,J):  Noise contribution of Y(J) to Y(I)
C        COH(II,I,J):  Simple coherency
C        AMP(II,I,J):  Amplitude spectra
C        ANG(II,I,J):  Phase spectra
C
cxx      IMPLICIT  REAL*8(A-H)
cxx      IMPLICIT  COMPLEX*16(O-Z)
cc      DIMENSION  A(MJ,MJ,MJ1), E(MJ,MJ), C(0:20)
cc      DIMENSION  P(0:NF,MJ,MJ), FNC(0:NF,MJ,MJ), COH(0:NF,MJ,MJ)
cc      DIMENSION  AMP(0:NF,MJ,MJ), ANG(0:NF,MJ,MJ)
cc      DIMENSION  ZA(0:200,10,10), ZB(10,10)
cc      DIMENSION  FC(0:200), FS(0:200), WRK(10,10)
cxx      DIMENSION  A(L,L,M), E(L,L), C(0:M)
cxx      DIMENSION  P(0:NF,L,L), FNC(0:NF,L,L), COH(0:NF,L,L)
cxx      DIMENSION  AMP(0:NF,L,L), ANG(0:NF,L,L)
cxx      DIMENSION  ZA(0:NF,L,L), ZB(L,L)
cxx      DIMENSION  FC(0:NF), FS(0:NF), WRK(L,L)
C
      INTEGER M, L, NF
      DOUBLE PRECISION A(L,L,M), E(L,L), FNC(0:NF,L,L), AMP(0:NF,L,L),
     1                 ANG(0:NF,L,L), COH(0:NF,L,L)
      COMPLEX(KIND(0d0)) P(0:NF,L,L)
c local
      DOUBLE PRECISION C(0:M), FC(0:NF), FS(0:NF), FSUM
      COMPLEX(KIND(0d0)) ZA(0:NF,L,L), ZB(L,L), WRK(L,L), ZDET, SUM
C
cxx      DO 30 I=1,L
      DO 31 I=1,L
      DO 30 J=1,L
      C(0) = 0.0D0
      IF( I.EQ.J )  C(0) = 1.0D0
      DO 10 II=1,M
cxx   10 C(II) = -A(I,J,II)
      C(II) = -A(I,J,II)
   10 CONTINUE
C
      CALL  FOURIE( C(0),M+1,NF+1,FC(0),FS(0) )
C
      DO 20 II=0,NF
cxx   20 ZA(II,I,J) = DCMPLX( FC(II),FS(II) )
cxxx      ZA(II,I,J) = DCMPLX( FC(II),FS(II) )
      ZA(II,I,J) = CMPLX( FC(II),FS(II), KIND=8 )
   20 CONTINUE
   30 CONTINUE
   31 CONTINUE
C
      DO 200  II=0,NF
cxx      DO 40 I=1,L
      DO 41 I=1,L
      DO 40 J=1,L
cxx   40 ZB(I,J) = ZA(II,I,J)
      ZB(I,J) = ZA(II,I,J)
   40 CONTINUE
   41 CONTINUE
C
cc      CALL  CINV( ZB,ZDET,L,10 )
      CALL  CINV( ZB,ZDET,L )
C
cxx      DO 60 I=1,L
      DO 61 I=1,L
      DO 60 J=1,L
      SUM = (0.0D0,0.0D0)
      DO 50 IJ=1,L
cxx   50 SUM = SUM + ZB(I,IJ)*E(IJ,J)
      SUM = SUM + ZB(I,IJ)*E(IJ,J)
   50 CONTINUE
cxx   60 WRK(I,J) = SUM
      WRK(I,J) = SUM
   60 CONTINUE
   61 CONTINUE
cxx      DO 80 I=1,L
      DO 81 I=1,L
      DO 80 J=1,L
      SUM = 0.0D0
      DO 70 IJ=1,L
cxx   70 SUM = SUM + WRK(I,IJ)*DCONJG( ZB(J,IJ) )
cxx      SUM = SUM + WRK(I,IJ)*DCONJG( ZB(J,IJ) )
      SUM = SUM + WRK(I,IJ)*CONJG( ZB(J,IJ) )
   70 CONTINUE
cxx   80 P(II,I,J) = SUM
      P(II,I,J) = SUM
   80 CONTINUE
   81 CONTINUE
C
C  ...  Amplitude and phase  ...
C
cxx      DO 90 I=1,L-1
      DO 91 I=1,L-1
      DO 90 J=I+1,L
cxx        AMP(II,I,J) = DSQRT( DREAL(P(II,I,J))**2 + DIMAG(P(II,I,J))**2)
cxx        ANG(II,I,J) = DATAN( DIMAG(P(II,I,J))/DREAL(P(II,I,J)) )
        AMP(II,I,J) = DSQRT( REAL(P(II,I,J))**2 + AIMAG(P(II,I,J))**2)
        ANG(II,I,J) = DATAN( AIMAG(P(II,I,J))/REAL(P(II,I,J)) )
cxx        IF( DIMAG(P(II,I,J)).GT.0.0D0.AND.DREAL(P(II,I,J)).LT.0.0D0 )
        IF( AIMAG(P(II,I,J)).GT.0.0D0 .AND. REAL(P(II,I,J)).LT.0.0D0 )
     *        ANG(II,I,J) = ANG(II,I,J) + 3.1415926535D0
cxx        IF( DIMAG(P(II,I,J)).LT.0.0D0.AND.DREAL(P(II,I,J)).LT.0.0D0 )
        IF( AIMAG(P(II,I,J)).LT.0.0D0 .AND. REAL(P(II,I,J)).LT.0.0D0 )
     *        ANG(II,I,J) = ANG(II,I,J) - 3.1415926535D0
   90 CONTINUE
   91 CONTINUE
C
C  ...  Simple coherency  ...
C
cxx      DO 100 I=1,L-1
      DO 101 I=1,L-1
      DO 100 J=I+1,L
cxx        COH(II,I,J) = (DREAL(P(II,I,J))**2 + DIMAG(P(II,I,J))**2)/
cxx     *                (P(II,I,I)*P(II,J,J))
cxx        COH(II,I,J) = DREAL((DREAL(P(II,I,J))**2 + DIMAG(P(II,I,J))**2)/
cxx     *                (P(II,I,I)*P(II,J,J)))
        COH(II,I,J) = DBLE((REAL(P(II,I,J))**2 + AIMAG(P(II,I,J))**2)/
     *                (P(II,I,I)*P(II,J,J)))
  100 CONTINUE
  101 CONTINUE
C
C  ...  Power contribution  ...
C
      DO 120 I=1,L
      FSUM = 0.0D0
      DO 110 J=1,L
cxx      FSUM = FSUM + ZB(I,J)*DCONJG( ZB(I,J) )*E(J,J)
cxx      FSUM = FSUM + DREAL(ZB(I,J)*DCONJG( ZB(I,J) )*E(J,J))
      FSUM = FSUM + REAL(ZB(I,J)*CONJG( ZB(I,J) )*E(J,J))
cxx  110 FNC(II,I,J) = FSUM
      FNC(II,I,J) = FSUM
  110 CONTINUE
  120 CONTINUE
C
  200 CONTINUE
C
      RETURN
      E N D
cc      SUBROUTINE CINV( X,DET,M,MJ )
      SUBROUTINE CINV( X,DET,M )
C
C  ...  Inverse and determinant of a complex matrix X  ...
C
C       INPUTS:
C          X:     M*M SQUARE MATRIX
C          M:     DIMENSION OF X
C          MJ:    ABSOLUTE DIMENSION OF X IN THE MAIN PROGRAM
C
C       OUTPUTS:
C          X:     INVERSE OF X
C          DET:   DETERMINANT OF X
C
cxx      IMPLICIT  COMPLEX*16 (A-H,O-Z)
cc      DIMENSION  X(MJ,MJ), IND(100)
cxx      DIMENSION  X(M,M), IND(M)
C
      INTEGER M
      COMPLEX(KIND(0d0)) DET, X(M,M)
c local
      INTEGER IND(M)
      COMPLEX(KIND(0d0)) XMAX, XTEMP
C
      DET = 1.0D0
      DO 60 L=1,M
      XMAX = 0.1D-10
      IMAX = 0
      DO 10 I=L,M
cxx      IF( CDABS( X(I,L) ).GT.CDABS( XMAX ) )  THEN
      IF( ABS( X(I,L) ).GT.ABS( XMAX ) )  THEN
         XMAX = X(I,L)
         IMAX = I
      END IF
   10 CONTINUE
      IND(L) = IMAX
      IF( IMAX.NE.L )  THEN
         IF( IMAX.LE.0 )  THEN
            DET = 0.0D0
            RETURN
         END IF
C
C     ROW INTERCHANGE
         DO 20 J=1,M
         XTEMP = X(IMAX,J)
         X(IMAX,J) = X(L,J)
cxx   20    X(L,J) = XTEMP
         X(L,J) = XTEMP
   20    CONTINUE
         DET = -DET
      END IF
      DET = DET*XMAX
      X(L,L) = 1.0D0
      DO 30 J=1,M
cxx   30 X(L,J) = X(L,J)/XMAX
      X(L,J) = X(L,J)/XMAX
   30 CONTINUE
      DO 50 I=1,M
      IF(I.NE.L)  THEN
         XTEMP  =X(I,L)
         X(I,L) = 0.0D0
         DO 40 J=1,M
cxx   40    X(I,J) = X(I,J) - XTEMP*X(L,J)
         X(I,J) = X(I,J) - XTEMP*X(L,J)
   40    CONTINUE
      END IF
   50 CONTINUE
   60 CONTINUE
C
C     COLUMN INTERCHANGE
      IF( M.GT.1 )  THEN
         DO 110 J=1,M-1
         JJ = IND(M-J)
         IF( JJ.NE.M-J )  THEN
            DO 100 I=1,M
            XTEMP    = X(I,JJ)
            X(I,JJ)  = X(I,M-J)
cxx  100       X(I,M-J) = XTEMP
            X(I,M-J) = XTEMP
  100       CONTINUE
         END IF
  110    CONTINUE
      END IF
      RETURN
      E N D
