C     PROGRAM 7.1  ARFIT
ccx      SUBROUTINE ARFITF( Y,N,LAG,NF,ISW,SIG2,AIC,MAR,A,PAR,SP )
      SUBROUTINE ARFIT( Y,N,LAG,NF,MJ2,ISW,SIG2,AIC,MAR,A,PAR,SP )
C
      INCLUDE 'TSSS.h'
C
C  ...  AR model fitting  ...
C
C     Inputs:
C        IDEV:    Input device for time series
C        LAG:     Highest order of AR model
C        ISW:     Estimation procedure
C                 = 1:  Yule-Walker method
C                 = 2:  Least squares (Householder) method
C                 = 3:  Partial autoregression method
C                 = 4:  PARCOR method
C                 = 5:  Burg's algorithm (MEM)
C     The following inputs are required in the subroutine READTS.
C        TITLE:   Caption of the data set
C        N:       Data length
C        FORMAT:  Reading format
C        Y(I):    Time series, (i=1,...,N)
C     Parameters:
C        NMAX:    Adjustable dimension of Y (NMAX.GE.N)
C        MJ,MJ2:  Adjustable dimensions
C        NF:      Number of frequencies for computing spectrum
C     @TEST.PN71:  12/12/90,1/7/91 Y.I. and G.K. 9/4/91
C
cc      PARAMETER( NMAX=200,MJ=20,MJ2=100,NF=200,IDEV=1 )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION Y(NMAX), COV(0:MJ,4), SIG2(0:MJ), SP(0:NF)
cc      DIMENSION D(MJ2), X(MJ2,MJ+1)
cc      DIMENSION AIC(0:MJ), A(MJ,MJ), B(MJ), PAR(MJ)
cc      DIMENSION FE(NMAX), BE(NMAX)
cxx      DIMENSION Y(N), COV(0:LAG,4), SIG2(0:LAG), SP(0:NF)
cxx      DIMENSION X(MJ2,LAG+1)
cxx      DIMENSION AIC(0:LAG), A(LAG,LAG), B(LAG+1), PAR(LAG)
cxx      DIMENSION FE(N), BE(N)
c
      INTEGER N, LAG, NF, MJ2, ISW, MAR
      DOUBLE PRECISION Y(N), SIG2(0:LAG), AIC(0:LAG), A(LAG,LAG), 
     1                 PAR(LAG), SP(0:NF)
c local
      INTEGER I, NSUM
      DOUBLE PRECISION COV(0:LAG,4), X(MJ2,LAG+1), B(LAG+1), FE(N),
     1                 BE(N), OUTMIN, OUTMAX, YMEAN
c
      EXTERNAL  SETXAR
      DATA  OUTMIN/-1.0D30/, OUTMAX/1.0D30/
c
      X(1:MJ2,1:LAG+1) = 0.0D0
      PAR(1:LAG) = 0.0D0
cc      READ(5,*)  LAG, ISW
cc      CALL  READTS( IDEV,Y,N )
      CALL  MEAN( Y,N,-1.0D30,1.0D30,NSUM,YMEAN )
C
C  ...  Yule-Walker method  ...
C
      IF( ISW.EQ.1 )  THEN
cc         CALL  UNICOR( Y,N,LAG,OUTMIN,OUTMAX,COV )
         CALL  UNICOR( Y,N,LAG,OUTMIN,OUTMAX,COV,YMEAN )
         CALL  ARYULE( COV,N,LAG,SIG2,AIC,PAR,A,MAR)
      END IF
C
C  ...  Householder method  ...
C
      IF( ISW.EQ.2 )  THEN
cc         CALL  REDUCT( SETXAR,Y,D,N-LAG,0,LAG,MJ2,X )
cc         CALL  REGRES( X,LAG,N-LAG,MJ2,MJ,A,SIG2,AIC,MAR )
         CALL  REDUCT( SETXAR,Y,N-LAG,0,LAG,MJ2,X )
         CALL  REGRES( X,LAG,N-LAG,MJ2,A,SIG2,AIC,MAR )
cxx         CALL  PARCOR( A(1,MAR),MAR,PAR )
         CALL  PARCOR( A(1,LAG),LAG,PAR )
      END IF
C
C  ...  PARCOR method  ...
C
      IF( ISW.GE.3 )  THEN
         CALL  ARPCOR( Y,FE,BE,SIG2,AIC,LAG,N,PAR,ISW-2,MAR )
         DO 10 I=1,LAG
cxx   10    CALL  ARCOEF( PAR,I,A(1,I) )
            CALL  ARCOEF( PAR,I,A(1,I) )
   10    CONTINUE
      END IF
C
C  ...  Power spectrum  ...
C
      CALL  ARMASP( A(1,MAR),MAR,B,0,SIG2(MAR),NF,SP )
C
cc      CALL  PRARSP( N,LAG,A,MAR,PAR,SIG2,AIC,SP,MJ,NF,ISW )
C      CALL  PTAR( PAR,AIC,SP,LAG,NF,MAR,SIG2,ISW )
cc      STOP
      RETURN
      E N D
      SUBROUTINE  ARPCOR( Y,FE,BE,SIG2,AIC,K,N,PARCOR,ISW,MAR )
C
C  ...  PARCOR method  ...
C
C     Inputs:
C        Y(I):    Time series
C        N:       Data length
C        K:       Highest AR order
C        ISW:     Estimation procedure
C                 = 1:   Partial autoregression
C                 = 2:   PARCOR
C                 = 3:   Burg's algorithm (MEM)
C        FE,BE:   Working area
C     Outputs:
C        SIG2(I): Innovation variance
C        AIC(I):  AIC
C        PARCOR(I):  PARCOR
C        MAR:     MAICE order
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  Y(N), FE(N), BE(N)
cxx      DIMENSION  SIG2(0:K), AIC(0:K), PARCOR(K)
cc      DIMENSION  A(50), B(50), FA(50), BA(50)
cxx      DIMENSION  A(K), B(K), FA(K), BA(K)
C
      INTEGER K, N, ISW, MAR
      DOUBLE PRECISION Y(N), FE(N), BE(N), SIG2(0:K), AIC(0:K),
     1                 PARCOR(K)
c local
      INTEGER I, L, M
      DOUBLE PRECISION A(K), B(K), FA(K), BA(K), PI, SUM, AICM, FB, FF,
     1                 BB, FE0, X
C
      PI = 3.1415926535D0
C
      SUM = 0.0D0
      DO 5 I=K+1,N
cxx    5 SUM = SUM + Y(I)**2
      SUM = SUM + Y(I)**2
    5 CONTINUE
      SIG2(0) = SUM/(N-K)
      AIC(0) = (N-K)*( DLOG(2*PI) + 1+ DLOG(SIG2(0))) + 2
      MAR = 0
      AICM = AIC(0)
C
      DO 10 I=1,N
      FE(I) = Y(I)
cxx   10 BE(I) = Y(I)
      BE(I) = Y(I)
   10 CONTINUE
C
      DO 100 M=1,K
      FB = 0.0D0
      FF = 0.0D0
      BB = 0.0D0
      DO 20 I=M+1,N
      FB = FB + FE(I)*BE(I-M)
      FF = FF + BE(I-M)**2
cxx   20 BB = BB + FE(I)**2
      BB = BB + FE(I)**2
   20 CONTINUE
C
      IF( ISW.EQ.1) THEN
        A(M) = FB/FF
        B(M) = FB/BB
      END IF
      IF( ISW.EQ.2 ) THEN
        A(M) = FB/(DSQRT(FF*BB))
        B(M) = FB/(DSQRT(FF*BB))
      END IF
      IF( ISW.EQ.3 ) THEN
        A(M) = FB/((FF+BB)/2)
        B(M) = FB/((FF+BB)/2)
      END IF
C
      DO 50 L=1,M-1
      A(L) = FA(L)-A(M)*BA(M-L)
cxx   50 B(L) = BA(L)-B(M)*FA(M-L)
      B(L) = BA(L)-B(M)*FA(M-L)
   50 CONTINUE
      DO 60 L=1,M
      FA(L) = A(L)
cxx   60 BA(L) = B(L)
      BA(L) = B(L)
   60 CONTINUE
      DO 70 I=M+1,N
      FE0 = FE(I)
      FE(I) = FE(I)-A(M)*BE(I-M)
cxx   70 BE(I-M) = BE(I-M)-B(M)*FE0
      BE(I-M) = BE(I-M)-B(M)*FE0
   70 CONTINUE
C
      SUM = 0.0D0
      DO 80 I=K+1,N
cxx   80 SUM = SUM + FE(I)**2
      SUM = SUM + FE(I)**2
   80 CONTINUE
      X = SUM
      PARCOR(M) = A(M)
      SIG2(M) = X/(N-K)
      AIC(M) = (N-K)*( DLOG(2*PI) + 1+ DLOG(SIG2(M)))+ 2*(M+1)
      IF( AIC(M).LT.AICM )  THEN
         MAR  = M
         AICM = AIC(M)
      END IF
C
  100 CONTINUE
C
      RETURN
      E N D
