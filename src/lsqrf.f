C     PROGRAM  5.1  LSQR
cxx      SUBROUTINE LSQRF(Y,N,K,AIC,SIG2,IMIN,A,DATA)
cxxx      SUBROUTINE LSQRF(Y,N,K,MJ1,AIC,SIG2,IMIN,A,DATA)
      SUBROUTINE LSQR(Y,N,K,PERIOD,MJ1,AIC,SIG2,IMIN,A,DATA)
C
      INCLUDE 'TSSS.h'
C
C  ...  The least squares method via Householder transformation  ...
C
C     The following inputs are required in the subroutine  READTS:
C        TITLE:   title of the data
C        FORMAT:  reading format of the data
C        Y(I):    time series
C     PARATERS:
C        MJ:      Adjustable dimension of Y (MJ.GE.N)
C        MJ1:     Adjustable dimension of X and D
C        MJ2:     Adjustable dimension of A, SIG2 and AIC
C        IDEV:    Input device
C        LAG:     Number of sine and cosine terms
C     @TEST.PN51:
C
cc      PARAMETER (MJ=1000,MJ1=100,MJ2=22,ISW=2,IDEV=1,LAG=10,K=LAG*2+1)
cxx      PARAMETER (MJ1=100)
cxx      IMPLICIT   REAL*8 (A-H,O-Z)
cc      CHARACTER  TITLE*72
cc      DIMENSION  Y(MJ), AIC(0:MJ2)
cc      DIMENSION  X(MJ1,MJ2), D(MJ1), A(MJ2,MJ2), SIG2(0:MJ2)
cc      COMMON     /CMDATA/  TITLE
cxx      DIMENSION  Y(N), DATA(N), AIC(0:K)
cxx      DIMENSION  X(MJ1,K+1), A(K,K), SIG2(0:K)
C
      INTEGER N, K, PERIOD, IMIN
      DOUBLE PRECISION Y(N), AIC(0:K), SIG2(0:K), A(K,K), DATA(N)
c local
      DOUBLE PRECISION X(MJ1,K+1)
C
      EXTERNAL   SETXTP
C
cc      MJ2=K
C
cc      CALL  READTS( IDEV,Y,N )
C
cc      CALL  REDUCT( SETXTP,Y,D,N,0,K,MJ1,X )
cx      CALL  REDUCT( SETXTP,Y,N,0,K,MJ1,X )
cxxx      CALL  REDUCT1( SETXTP,Y,N,0,K,MJ1,X )
      CALL  REDUCT2( SETXTP,Y,N,0,K,PERIOD,MJ1,X )
C
cc      CALL  REGRES( X,K,N,MJ1,MJ2,A,SIG2,AIC,IMIN )
      CALL  REGRES( X,K,N,MJ1,A,SIG2,AIC,IMIN )
C
cc      CALL  PRREG( N,K,MJ2,A,SIG2,AIC )
cxx      CALL  PTTPL( Y,N,A(1,IMIN),IMIN,DATA )
cxxx      CALL  PTTPL( N,A(1,IMIN),IMIN,DATA )
      CALL  PTTPL( N,A(1,IMIN),IMIN,PERIOD,DATA )
C
cc      STOP
      RETURN
      E N D
C
cc      SUBROUTINE  REDUCT( SETX,Z,D,NMK,N0,K,MJ1,X )
      SUBROUTINE  REDUCT2( SETX,Z,NMK,N0,K,PERIOD,MJ1,X )
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
C          PERIOD:       Period of one cycle
C        MJ1:     Adjustable dimension of X
C     Output:
C        X(I,J):  data matrix
C
cxx      IMPLICIT  REAL*8( A-H,O-Z )
cc      DIMENSION  X(MJ1,1) , D(1), Z(1)
cx      DIMENSION  X(MJ1,1) , Z(1)
cxx      DIMENSION  X(MJ1,K+1) , Z(N0+NMK)
      INTEGER NMK, N0, K, MJ1, PERIOD
      DOUBLE PRECISION Z(N0+NMK), X(MJ1,K+1)
C
      L = MIN0( NMK,MJ1 )
      K1 = K + 1
      N1 = L
C
cxxx      CALL  SETX( Z,N0,L,K,MJ1,0,X )
      CALL  SETX( Z,N0,L,K,PERIOD,MJ1,0,X )
cc      CALL  HUSHLD( X,D,MJ1,L,K1 )
      CALL  HUSHLD( X,MJ1,L,K1 )
      IF( N1 .GE. NMK )  RETURN
C
   10 L = MIN0( NMK-N1,MJ1-K1 )
C
      LK = L + K1
      N2 = N0 + N1
cxxx      CALL  SETX( Z,N2,L,K,MJ1,1,X )
      CALL  SETX( Z,N2,L,K,PERIOD,MJ1,1,X )
cc      CALL  HUSHLD( X,D,MJ1,LK,K1 )
      CALL  HUSHLD( X,MJ1,LK,K1 )
      N1 = N1 + L
      IF( N1.LT.NMK )  GO TO 10
C
      RETURN
C
      E N D
C
cxxx      SUBROUTINE  SETXTP( Z,N0,L,K,MJ1,JSW,X )
      SUBROUTINE  SETXTP( Z,N0,L,K,PERIOD,MJ1,JSW,X )
C
C  ...  Data matrix for trigonometric polynomial regression  ...
C
C     Inputs:
C        Z(I):    Data vector
C        N0:      Origin of the current observations
C        L:       Number of current observations
C        K:       Number of regressors
C          PERIOD:       Period of one cycle
C        MJ1:     Adjustable dimension of X
C        JSW=0:   Make initial data matrix
C           =1:   Apend L*(K+1) data matrix below the triangular one
C     Output:
C        X(I,J):  Data matrix
C
cc      REAL*8  X(MJ1,1), Z(1), W
cxx      REAL*8  X(MJ1,K+1), Z(N0+L), W
      INTEGER N0, L, K, MJ1, JSW, PERIOD
      DOUBLE PRECISION Z(N0+L), X(MJ1,K+1)
c local
      DOUBLE PRECISION W
C
cxxx      W = 2*3.1415926536D0/365.0D0
      W = 2*3.1415926536D0/DBLE(period)
      I0 = 0
      IF( JSW .EQ. 1 )     I0 = K+1
cxx      DO 10  I=1,L
      DO 11  I=1,L
        II = I + I0
        X(II,K+1) = Z(N0+I)
        X(II,1) = 1.0D0
      DO 10  J=1,(K-1)/2
      X(II,J*2)   = DSIN( W*(I+N0)*J )
cxx   10 X(II,J*2+1) = DCOS( W*(I+N0)*J )
      X(II,J*2+1) = DCOS( W*(I+N0)*J )
   10 CONTINUE
   11 CONTINUE
C
      RETURN
C
      E N D

cc      SUBROUTINE PTTPL( Y,N,A,M )
cxx      SUBROUTINE PTTPL( Y,N,A,M,DATA )
cxxx      SUBROUTINE PTTPL( N,A,M,DATA )
      SUBROUTINE PTTPL( N,A,M,PERIOD,DATA )
C
C  ...  This subroutine draws fitted trigonometric polynomial  ...
C
C     Inputs:
C        Y(I):   Original data
C        N:      Data length
C        A(I):   Regression coefficients
C        M:      Order of regression model
C          PERIOD:       Period of one cycle
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      CHARACTER  VNAME*8
cc      DIMENSION  Y(N), A(M), DATA(1000)
cxx      DIMENSION  Y(N), A(M), DATA(N)
cc      DIMENSION  VNAME(5), VALUE(5)
C
      INTEGER N, M, PERIOD
      DOUBLE PRECISION A(M), DATA(N)
c local
      DOUBLE PRECISION W, SUM
C
cc      WX = 20.0
cc      WY = 6.0
cc      IPOS = 1
cc      CALL  PLOTS
C     call  plots( 1,0,0,1,0 )
C     call  form( 1 )
C     call  factor( 10.0 )
cc      CALL  HEADER( 'TRIGONOMETRIC REGRESSION',24,0,VNAME,VALUE )
C
cc      CALL  DELX( 0.0,DBLE(N),DX )
cc      CALL  MAXMIN( Y,N,YMIN,YMAX,DY )
cc      CALL  AXISXY( 3.0D0,15.0-WY,WX,WY,0.0D0,DBLE(N),YMIN,YMAX,
cc     *              DX,DY,0.2D0,1,10,2 )
cc      CALL  SYMBOL( 0.0,SNGL(WY)+0.1,0.2,'ORIGINAL AND TREND',0.0,18 )
cc      CALL  NEWPEN( 1 )
cc      CALL  PLOTY ( Y,N,YMIN,YMAX,WX,WY,IPOS,1 )
C
cxxx      W = 2*3.1415926536D0/365.0D0
      W = 2*3.1415926536D0/DBLE(period)
      DO 20 I=1,N
      SUM = A(1)
      DO 10 J=1,10
      IF( 2*J.LE.M )   SUM = SUM + A(2*J)*DSIN(W*I*J)
cxx   10 IF( 2*J+1.LE.M ) SUM = SUM + A(2*J+1)*DCOS(W*I*J)
      IF( 2*J+1.LE.M ) SUM = SUM + A(2*J+1)*DCOS(W*I*J)
   10 CONTINUE
cxx   20 DATA(I) = SUM
      DATA(I) = SUM
   20 CONTINUE
cc      CALL  NEWPEN( 2 )
cc      CALL  PLOTY( DATA,N,YMIN,YMAX,WX,WY,IPOS,1 )
cc      CALL  PLOTE
C     call  plot( 0.0,0.0,999 )
C
      RETURN
      E N D
