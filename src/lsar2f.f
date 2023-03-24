C     PROGRAM  8.2  LSAR2
      SUBROUTINE LSAR2( Y,N,K,N0,N1,N2,NE, AICS,AICMIN,MMIN )
C
      INCLUDE 'TSSS.h'
C
C  ...  Estimation of the change point  ...
C
C     Inputs:
C        K:       Highest order od AR model
C        [N0,NE]: Time interval used for model fitting
C        [N1,N2]: Candidate for change point
C     The following inputs are required in the subroutine READTS.
C        TITLE:   Caption of the data set
C        N:       Data length
C        FORMAT:  Reading format
C        Y(I):    Time series, (i=1,...,N)
C     Parameters:
C        IDEV:    Input device for time series
C        NMAX:    Adjustable dimension of Y (NMAX.GE.N)
C        MJ,K,NSPAN:  Adjustable dimensions
C     MODIFIED  2/15/93
C
cc      PARAMETER(NMAX=3000,MJ=100,NSPAN=1000,KMAX=20,K1=KMAX+1,IDEV=1)
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  Y(NMAX)
cc      DIMENSION  X(MJ,K1), AIC1(NSPAN), AIC2(NSPAN), AICS(NSPAN)
cxx      DIMENSION  Y(N)
cxx      DIMENSION  AIC1(N2-N1), AIC2(N2-N1), AICS(N2-N1)
C     DIMENSION  DATA(200)
C
      INTEGER N, K, N0, N1, N2, NE, MMIN
      DOUBLE PRECISION Y(N), AICS(N2-N1), AICMIN
c local
      DOUBLE PRECISION AIC1(N2-N1), AIC2(N2-N1)
C
cc      READ( 5,* )  K, N0, N1, N2, NE
      NS = 1
      M = N2 - N1
c-----
cc      MJ = NE-N0
       MJ = NS + K + 1
c-----
C
C  ...  Read time series  ...
C
cc      CALL  READTS( IDEV,Y,N )
C
C  ...   First models and AIC1's  ...
C
cc      CALL  UPDATE( X,Y,N0,N1,M,NS,K,MJ,AIC1 )
      CALL  UPDATE( Y,N,N0,N1,M,NS,K,MJ,AIC1 )
C
C  ...   Second models and AIC2's  ...
C
cc      CALL  BUPDAT( X,Y,N2,NE,M,NS,K,MJ,AIC2 )
      CALL  BUPDAT( Y,N2,NE,M,NS,K,MJ,AIC2 )
C
C  ...   AICs of the locally stationary AR models  ...
C
      DO 20 I=1,M
cxx   20 AICS(I) = AIC1(I) + AIC2(I)
      AICS(I) = AIC1(I) + AIC2(I)
   20 CONTINUE
C
      AICMIN = 1.0D30
c-----
      NMIN = 1
c-----
      DO 60 I=1,M
      IF( AICS(I) .GT. AICMIN )  GO TO 60
         AICMIN = AICS(I)
         MMIN = I
   60 CONTINUE
C
C  ...   Print out and plot results  ...
C
cc      CALL  PRLSAR( AICS,AICMIN,MMIN,M,N0,N1,N2,NE,K )
C      CALL  PTLSAR( Y,AICS,AICMIN,MMIN,N,M,N0,N1,N2,NE,K )
C
cc      STOP
      RETURN
      E N D
cc      SUBROUTINE  UPDATE( X,Z,N0,N1,M,NS,K,MJ,AIC )
      SUBROUTINE  UPDATE( Z,N,N0,N1,M,NS,K,MJ,AIC )
C
C  ...  Fit AR models to first part of data  ...
C
C     Inputs:
C        X:     Working area
C        Z:     Time series
C        N0,N1,M,NS:  Fit AR models on [N0,N1],[N0,N1+NS],...,
C                                                    [N0,N1+(M-1)*NS]
C        K:     Highest order of AR models
C        MJ:    Adjustable dimension
C     Output:
C        AIC(I):  AIC of the AR model fitted on [N0,N1+(I-1)*NS]
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  X(MJ,1), AIC(1), D(200), A(20,20), AICS(0:20)
cc      DIMENSION  Z(N1), SIG2(0:20)
cxx      DIMENSION  X(MJ,K+1), AIC(M), A(K,K), AICS(0:K)
cxx      DIMENSION  Z(N), SIG2(0:K)
C
C
      INTEGER N, N0, N1, M, NS, K, MJ
      DOUBLE PRECISION Z(N), AIC(M)
c local
      DOUBLE PRECISION SIG2(0:K), X(MJ,K+1), A(K,K), AICS(0:K)
      EXTERNAL  SETXAR
C
      NMK = N1 - K - N0
C
cc      CALL  REDUCT( SETXAR,Z,D,NMK,N0,K,MJ,X )
      CALL  REDUCT( SETXAR,Z,NMK,N0,K,MJ,X )   
C
      DO 100  I=1,M
      II = N1 + (I-1)*NS
c      CALL  REGRES( X,K,II-K-N0,MJ,20,A,SIG2,AICS,IMIN )
      CALL  REGRES( X,K,II-K-N0,MJ,A,SIG2,AICS,IMIN )
      AIC(I) = AICS(IMIN)
C
      CALL  SETXAR( Z,II-K,NS,K,MJ,1,X )
cc      CALL  HUSHL2( X,D,MJ,K+1+NS,K+1 )
      CALL  HUSHL2( X,MJ,K+1+NS,K+1 )
C
  100 CONTINUE
      RETURN
      E N D
cc      SUBROUTINE  BUPDAT( X,Z,N2,N,M,NS,K,MJ,AIC )
      SUBROUTINE  BUPDAT( Z,N2,N,M,NS,K,MJ,AIC )
C
C  ...  Fit AR models to the second part of data  ...
C
C     Inputs:
C        X:     Working area
C        Z:     Time series
C        N0,N2,N,M,NS:  Fit AR models on [N2-(M-1)*NS,N],...,
C                                           [N2-NS,S],[N2,N]
C        K:     Highest order of AR models
C        MJ:    Adjustable dimension
C     Output:
C        AIC(I):  AIC of the AR model fitted on [N0,N1+(I-1)*NS]
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  X(MJ,1), AIC(1), D(200), A(20,20), AICS(0:20)
cc      DIMENSION  Z(1), SIG2(0:20)
cxx      DIMENSION  X(MJ,K+1), AIC(M), A(K,K), AICS(0:K)
cx      DIMENSION  Z(1), SIG2(0:K)
cxx      DIMENSION  Z(N), SIG2(0:K)
C
      INTEGER N2, N, M, NS, K, MJ
      DOUBLE PRECISION Z(N), AIC(M)
c local
      DOUBLE PRECISION SIG2(0:K), X(MJ,K+1), A(K,K), AICS(0:K)
      EXTERNAL  SETXAR
C
      NMK = N - N2
      J = M
C
cc      CALL  REDUCT( SETXAR,Z,D,NMK,N2-K-NS,K,MJ,X )
      CALL  REDUCT( SETXAR,Z,NMK,N2-K-NS,K,MJ,X )
C
      DO 100  I=1,M
      II = N2 - (I-2)*NS
cc      CALL  REGRES( X,K,N-II,MJ,20,A,SIG2,AICS,IMIN )
      CALL  REGRES( X,K,N-II,MJ,A,SIG2,AICS,IMIN )
      AIC(J) = AICS(IMIN)
      J = J - 1
C
      CALL  SETXAR( Z,II-K-NS,NS,K,MJ,1,X )
cc      CALL  HUSHL2( X,D,MJ,K+1+NS,K+1 )
      CALL  HUSHL2( X,MJ,K+1+NS,K+1 )
C
  100 CONTINUE
      RETURN
      E N D
cc      SUBROUTINE  HUSHL2( X,D,MJ1,N,K )
      SUBROUTINE  HUSHL2( X,MJ1,N,K )
C
C  ...  Householder transformation for adding new observations  ...
C
C     Inputs:
C        X:     Data matrix (Householder reduced form augmented with
C                            new observations)
C        D:     Working area
C        MJ1:   Adjustable dimension
C        N:     Number of rows of X
C        K:     Number of columns of X
C     Output:
C        X:     Householder reduced form
C
cxx      IMPLICIT  REAL*8 (A-H,O-Z)
cxx      DIMENSION  X(MJ1,K) , D(MJ1)
C
      INTEGER MJ1, N, K
      DOUBLE PRECISION X(MJ1,K)
c local
      DOUBLE PRECISION D(MJ1), TOL, D1, H, G, S
C
           TOL = 1.0D-30
C
cxx      DO 100  II=1,K
      DO 101  II=1,K
         D1 = X(II,II)
         H  = D1**2
         DO 10  I=K+1,N
            D(I) = X(I,II)
cxx   10       H = H + D(I)**2
            H = H + D(I)**2
   10    CONTINUE
         IF( H .GT. TOL )  GO TO 20
         G = 0.0D00
         GO TO 100
   20    G = DSQRT( H )
         IF( D1 .GE. 0.0D00 )   G = -G
         H  = H - D1*G
         D1 = D1 - G
C
         IF( II .EQ. K )  GO TO 100
         II1 = II+1
         DO 60  J=II1,K
            S = D1*X(II,J)
            DO 40  I=K+1,N
cxx   40       S = S + D(I)*X(I,J)
               S = S + D(I)*X(I,J)
   40       CONTINUE
            S = S/H
            X(II,J) = X(II,J) - D1*S
            DO 50  I=K+1,N
cxx   50       X(I,J) = X(I,J) - D(I)*S
               X(I,J) = X(I,J) - D(I)*S
   50       CONTINUE
   60    CONTINUE
  100 X(II,II) = G
  101 CONTINUE
C
      RETURN
      E N D
