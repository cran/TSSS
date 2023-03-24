C     PROGRAM  11.1  POLREG
      SUBROUTINE POLREG(Y,N,K,A,SIG2,AIC,IMIN,DATA)
C
      INCLUDE 'TSSS.h'
C
C  ...  Polynomial regression model  ...
C
C     INPUT:
C        K:       Order of polynomial regression
C     PARATERS:
C        MJ:      Adjustable dimension of Y (MJ.GE.N)
C        MJ1:     Adjustable dimension of X and D
C        MJ2:     Adjustable dimension of A, SIG2 and AIC
C        IDEV:    Input device
C        LAG:     Number of sine and cosine terms
C     @TEST.PN51:
C
cc      PARAMETER( MJ=1000,MJ1=200,MJ2=21,ISW=2,IDEV=1 )
      PARAMETER( MJ1=200 )
cxx      IMPLICIT   REAL*8 (A-H,O-Z)
cc      DIMENSION  Y(MJ), AIC(0:MJ2)
cc      DIMENSION  X(MJ1,MJ2), D(MJ1), A(MJ2,MJ2), SIG2(0:MJ2)
cxx      DIMENSION  Y(N), AIC(0:K)
cxx      DIMENSION  X(MJ1,K+1), A(K,K), SIG2(0:K), DATA(N)
C
      INTEGER N, K, IMIN
      DOUBLE PRECISION Y(N), A(K,K), SIG2(0:K), AIC(0:K), DATA(N)
c local
      DOUBLE PRECISION X(MJ1,K+1), SUM, XX
C
      EXTERNAL   SETXPL
C
cc      READ( 5,* )  K
cc      CALL  READTS( IDEV,Y,N )
C
cc      CALL  REDUCT( SETXPL,Y,D,N,0,K+1,MJ1,X )
cx      CALL  REDUCT( SETXPL,Y,N, 0,K1,MJ1,X )
      CALL  REDUCT1( SETXPL,Y,N, 0,K,MJ1,X )
C
cc      CALL  REGRES( X,K+1,N,MJ1,MJ2,A,SIG2,AIC,IMIN )
      CALL  REGRES( X,K,N,MJ1,A,SIG2,AIC,IMIN )
C
cc      CALL  PRPOL( N,K+1,MJ2,A,SIG2,AIC )
C      CALL  PTPOL( Y,N,A(1,IMIN),IMIN )
c-----    PTPOL   ----
      DO 20 I=1,N
      SUM = 0.0D0
      XX  = 1.0D0
      DO 10 J=1,IMIN
      SUM = SUM + A(J,IMIN)*XX
cxx   10 XX = XX*I
      XX = XX*I
   10 CONTINUE
cxx   20 DATA(I) = SUM
      DATA(I) = SUM
   20 CONTINUE
c--------------------
cc      STOP
      RETURN
      E N D
cc      SUBROUTINE  REDUCT( SETX,Z,D,NMK,N0,K,MJ1,X )
      SUBROUTINE  REDUCT1( SETX,Z,NMK,N0,K,MJ1,X )
C
C  ...  Successive Householder reduction  ...
C
C     Inputs:
C        SETX:    Name of the subroutine for making X(I,J)
C        Z(I):    Data vector
C        D(I):    Working area
C        NMK:     Number of actually used observations
C        N0:      Time point of the previous set ofobservations
C        K:       Heighest order of the model
C        MJ1:     Adjustable dimension of X
C     Output:
C        X(I,J):  data matrix
C
cxx      IMPLICIT  REAL*8( A-H,O-Z )
cc      DIMENSION  X(MJ1,1) , D(1), Z(1)
cx      DIMENSION  X(MJ1,1) , Z(1)
cxx      DIMENSION  X(MJ1,K+1) , Z(N0+NMK)
      INTEGER NMK, N0, K, MJ1
      DOUBLE PRECISION Z(N0+NMK), X(MJ1,K+1)
C
      L = MIN0( NMK,MJ1 )
      K1 = K + 1
      N1 = L
C
      CALL  SETX( Z,N0,L,K,MJ1,0,X )
cc      CALL  HUSHLD( X,D,MJ1,L,K1 )
      CALL  HUSHLD( X,MJ1,L,K1 )
      IF( N1 .GE. NMK )  RETURN
C
   10 L = MIN0( NMK-N1,MJ1-K1 )
C
      LK = L + K1
      N2 = N0 + N1
      CALL  SETX( Z,N2,L,K,MJ1,1,X )
cc      CALL  HUSHLD( X,D,MJ1,LK,K1 )
      CALL  HUSHLD( X,MJ1,LK,K1 )
      N1 = N1 + L
      IF( N1.LT.NMK )  GO TO 10
C
      RETURN
C
      E N D
C
C
      SUBROUTINE  SETXPL( Z,N0,L,K,MJ1,JSW,X )
C
C  ...  Data matrix for polynomial regression  ...
C
C     Inputs:
C        Z(I):    Data vector
C        N0:      Origin of the current observations
C        L:       Number of current observations
C        K:       Number of regressors
C                 (Order of polynomial = K-1)
C        MJ1:     Adjustable dimension of X
C        JSW=0:   Make initial data matrix
C           =1:   Apend L*(K+1) data matrix below the triangular one
C     Output:
C        X(I,J):  Data matrix
C
cc      REAL*8  X(MJ1,1), Z(1)
cxx      REAL*8  X(MJ1,K+1), Z(N0+L)
C
      INTEGER N0, L, K, MJ1, JSW
      DOUBLE PRECISION Z(N0+L), X(MJ1,K+1)
c local
      DOUBLE PRECISION SUM
C
      I0 = 0
      IF( JSW .EQ. 1 )     I0 = K+1
cxx      DO 10  I=1,L
      DO 11  I=1,L
        II = I + I0
        X(II,K+1) = Z(N0+I)
        X(II,1) = 1.0D0
        SUM = 1.0D0
      DO 10  J=1,K-1
        SUM = SUM*(N0+I)
cxx   10 X(II,J+1) = SUM
      X(II,J+1) = SUM
   10 CONTINUE
   11 CONTINUE
C
      RETURN
C
      E N D
