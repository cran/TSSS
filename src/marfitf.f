C     PROGRAM 7.2  MARFIT
      SUBROUTINE MARFIT( Y,N,L,LAG,AMIN,VMIN,AIC,MMIN )
C
      INCLUDE 'TSSS.h'
C
C  ...  This program fits multivariate AR model  ...
C
C     The following inputs are required in the subroutine READMD.
C        TITLE:   title of the data
C        N:       data length
C        ID:      dimension of the observation
C        IFM:     = 1   ((Y(I,J),J=1,ID),I=1,N)
C                 = 2   ((Y(I,J),I=1,N),J=1,ID)
C        Y(I,J):  multi-variate time series
C     Parameters:
C        MJ:      adjustable dimension of Y; (MJ.GE.N)
C        MJ1:     adjustable dimension of Y, C, R; (MJ1.GE.ID)
C        LAG:     maximum lag of the cross-covariance function
C
cc      PARAMETER( MJ0=500,MJ=4,MJ1=10,LAG=MJ1,IDEV=1 )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  Y(MJ0,MJ), OUTMIN(10), OUTMAX(10)
cc      DIMENSION  C(0:MJ1,MJ,MJ), R(0:MJ1,MJ,MJ), YMEAN(10)
cc      DIMENSION  AMIN(MJ1,MJ,MJ), VMIN(MJ,MJ), AIC(0:MJ1)
cc      DATA  OUTMIN/10*-1.0D30/, OUTMAX/10*1.0D30/
cxx      DIMENSION  Y(N,L), OUTMIN(L), OUTMAX(L)
cxx      DIMENSION  C(0:LAG,L,L), R(0:LAG,L,L), YMEAN(L)
cxx      DIMENSION  AMIN(LAG,L,L), VMIN(L,L), AIC(0:LAG)
C
      INTEGER N, L, LAG, MMIN
      DOUBLE PRECISION Y(N,L), AMIN(LAG,L,L), VMIN(L,L), AIC(0:LAG)
c local
      DOUBLE PRECISION OUTMIN(L), OUTMAX(L), C(0:LAG,L,L), R(0:LAG,L,L),
     1                 YMEAN(L)
C
cxx      DO 100 I = 1,L
cxx         OUTMIN(I) = -1.0D30
cxx         OUTMAX(I) = 1.0D30
cxx  100 CONTINUE
      OUTMIN(1:L) = -1.0D30
      OUTMAX(1:L) = 1.0D30
C
C  ...  read in multivariate time series  ...
C
cc      CALL  READMD( IDEV,MJ0,Y,N,L )
C
C  ...  cross-covariance function  ...
C
cc      CALL  CRSCOR( Y,N,L,LAG,MJ0,MJ,OUTMIN,OUTMAX,C,R,YMEAN )
      CALL  CRSCOR( Y,N,L,LAG,OUTMIN,OUTMAX,C,R,YMEAN )
C
C  ...  Yule-Walker method for MAR model  ...
C
cc      CALL  MYULE( L,LAG,N,C,MJ,MJ1,AMIN,VMIN,MMIN,AIC )
      CALL  MYULE( L,LAG,N,C,AMIN,VMIN,MMIN,AIC )
C
C  ...  print out estimated model  ...
C
cc      CALL  PRMAR( L,N,LAG,MJ,MJ1,AMIN,VMIN,MMIN,AIC )
cc      STOP
      RETURN
      E N D
cc      SUBROUTINE  MYULE( L,LAG,N,C,MJ,MJ1,AMIN,VMIN,MMIN,AIC )
      SUBROUTINE  MYULE( L,LAG,N,C,AMIN,VMIN,MMIN,AIC )
C
C  ...  Yule-Walker method for fitting multivariate AR model  ...
C
C     Inputs:
C        L:         Dimension of the time series
C        N:         Data length
C        LAG:       Highest order of fitted AR models
C        C(K,I,J):  Crosscovariance function
C        MJ, MJ1:   Adjustable dimensions
C     Outputs:
C        AMIN(K,I,J):  AR coefficient of the AIC best model
C        VMIN(I,J):    Innovation covariance matrix of AIC best model
C        MMIN:         MAICE order
C        AIC(I):       AIC's of the AR models
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  C(0:MJ1,MJ,MJ)
cc      DIMENSION  AMIN(MJ1,MJ,MJ), VMIN(MJ,MJ), AIC(0:MJ1)
cc      DIMENSION  A(20,10,10), A0(20,10,10), B(20,10,10), B0(20,10,10)
cc      DIMENSION  V(10,10), U(10,10), W(10,10)
cxx      DIMENSION  C(0:LAG,L,L)
cxx      DIMENSION  AMIN(LAG,L,L), VMIN(L,L), AIC(0:LAG)
cxx      DIMENSION  A(LAG,L,L), A0(LAG,L,L), B(LAG,L,L), B0(LAG,L,L)
cxx      DIMENSION  V(L,L), U(L,L), W(L,L)
C
      INTEGER L, LAG, N, MMIN
      DOUBLE PRECISION C(0:LAG,L,L), AMIN(LAG,L,L), VMIN(L,L),
     1                 AIC(0:LAG)
c local
      INTEGER I, IJ, J, K, M, MJ2
      DOUBLE PRECISION A(LAG,L,L), A0(LAG,L,L), B(LAG,L,L), B0(LAG,L,L),
     1                 V(L,L), U(L,L), W(L,L), UDET, VDET, AICMIN, SUM,
     2                 SUM1, SUM2, PI2
C
      DATA  PI2 /6.28318531D0/
cc      MJ2 = 10
      MJ2 = L
C
cxx      DO 10 I=1,L
      DO 11 I=1,L
      DO 10 J=1,L
      V(I,J) = C(0,I,J)
      U(I,J) = C(0,I,J)
cxx   10 VMIN(I,J) = V(I,J)
      VMIN(I,J) = V(I,J)
   10 CONTINUE
   11 CONTINUE
      CALL  INVDET( U,UDET,L,MJ2 )
      CALL  INVDET( V,VDET,L,MJ2 )
      AIC(0) = N*L*(DLOG(PI2) + 1) + N*DLOG(VDET) + L*(L+1)
      AICMIN = AIC(0)
      MMIN = 0
C
      DO 300 M=1,LAG
C
cxx      DO 110 I=1,L
      DO 111 I=1,L
      DO 110 J=1,L
      SUM = C(M,I,J)
cxx      DO 100 K=1,M-1
      DO 101 K=1,M-1
      DO 100 IJ=1,L
cxx  100 SUM = SUM - A0(K,I,IJ)*C(M-K,IJ,J)
      SUM = SUM - A0(K,I,IJ)*C(M-K,IJ,J)
  100 CONTINUE
  101 CONTINUE
cxx  110 W(I,J) = SUM
      W(I,J) = SUM
  110 CONTINUE
  111 CONTINUE
C
cxx      DO 130 I=1,L
      DO 131 I=1,L
      DO 130 J=1,L
      SUM1 = 0.0D0
      SUM2 = 0.0D0
      DO 120 IJ=1,L
      SUM1 = SUM1 + W(I,IJ)*U(IJ,J)
cxx  120 SUM2 = SUM2 + W(IJ,I)*V(IJ,J)
      SUM2 = SUM2 + W(IJ,I)*V(IJ,J)
  120 CONTINUE
      A(M,I,J) = SUM1
cxx  130 B(M,I,J) = SUM2
      B(M,I,J) = SUM2
  130 CONTINUE
  131 CONTINUE
C
cxx      DO 150 K=1,M-1
cxx      DO 150 I=1,L
      DO 152 K=1,M-1
      DO 151 I=1,L
      DO 150 J=1,L
      SUM1 = A0(K,I,J)
      SUM2 = B0(K,I,J)
      DO 140 IJ=1,L
      SUM1 = SUM1 - A(M,I,IJ)*B0(M-K,IJ,J)
cxx  140 SUM2 = SUM2 - B(M,I,IJ)*A0(M-K,IJ,J)
      SUM2 = SUM2 - B(M,I,IJ)*A0(M-K,IJ,J)
  140 CONTINUE
      A(K,I,J) = SUM1
cxx  150 B(K,I,J) = SUM2
      B(K,I,J) = SUM2
  150 CONTINUE
  151 CONTINUE
  152 CONTINUE
C
cxx      DO 170 I=1,L
      DO 171 I=1,L
      DO 170 J=1,L
      SUM1 = C(0,I,J)
      SUM2 = C(0,I,J)
cxx      DO 160 K=1,M
      DO 161 K=1,M
      DO 160 IJ=1,L
      SUM1 = SUM1 - A(K,I,IJ)*C(K,J,IJ)
cxx  160 SUM2 = SUM2 - B(K,I,IJ)*C(K,IJ,J)
      SUM2 = SUM2 - B(K,I,IJ)*C(K,IJ,J)
  160 CONTINUE
  161 CONTINUE
      V(I,J) = SUM1
      W(I,J) = SUM1
cxx  170 U(I,J) = SUM2
      U(I,J) = SUM2
  170 CONTINUE
  171 CONTINUE
C
cxx      DO 180 K=1,M
cxx      DO 180 J=1,L
      DO 182 K=1,M
      DO 181 J=1,L
      DO 180 I=1,L
      A0(K,I,J) = A(K,I,J)
cxx  180 B0(K,I,J) = B(K,I,J)
      B0(K,I,J) = B(K,I,J)
  180 CONTINUE
  181 CONTINUE
  182 CONTINUE

C
      CALL  INVDET( U,UDET,L,MJ2 )
      CALL  INVDET( V,VDET,L,MJ2 )
      AIC(M) = N*L*(DLOG(PI2) + 1) + N*DLOG(VDET) + L*(L+1) + 2*L*L*M
C
      IF( AIC(M).LT.AICMIN )  THEN
      AICMIN = AIC(M)
      MMIN   = M
cxx      DO 200 I=1,L
cxx      DO 200 J=1,L
      DO 202 I=1,L
      DO 201 J=1,L
      VMIN(I,J) = W(I,J)
      DO 200 K=1,M
cxx  200 AMIN(K,I,J) = A(K,I,J)
      AMIN(K,I,J) = A(K,I,J)
  200 CONTINUE
  201 CONTINUE
  202 CONTINUE
      END IF
  300 CONTINUE
C
      RETURN
      E N D
      SUBROUTINE INVDET( X,DET,M,MJ )
C
C  ...  Inverse and determinant of matrix X  ...
C
C       Inputs:
C          X:     M*M square matrix
C          M:     Dimension of X
C          MJ:    Adjustable dimension of X
C       Outputs:
C          X:     Inverse of X
C          DET:   Determinant of X
C
cxx      IMPLICIT  REAL*8 (A-H,O-Z)
cc      DIMENSION  X(MJ,MJ), IND(100)
cxx      DIMENSION  X(MJ,MJ), IND(M)
C
      INTEGER M, MJ
      DOUBLE PRECISION X(MJ,MJ), DET
c local
      INTEGER I, IMAX, J, JJ, L, IND(M)
      DOUBLE PRECISION XMAX, XTEMP
C
      DET = 1.0D0
      DO 60 L=1,M
      XMAX = 0.1D-10
      IMAX = 0
      DO 10 I=L,M
      IF( DABS( X(I,L) ).GT.DABS( XMAX ) )  THEN
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
C     Row interchange
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
   40 CONTINUE
      END IF
   50 CONTINUE
   60 CONTINUE
C
C     Column interchange
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
