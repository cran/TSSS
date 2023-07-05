C     PROGRAM  3.2  FFTPER
cx      SUBROUTINE FFTPERF( Y,N,IWINDW,PE,SPE,NP )
      SUBROUTINE FFTPERF( Y,N,IWINDW,PE,SPE,NP,IFG )
C
      INCLUDE 'TSSS.h'
C
C  ...  This program computes periodogram via FFT  ...
C
C     The following inputs are required in READTS.
C        TITLE:   title of the data
C        N:       data length
C        Y(I):    time series
C     Parameters:
C        NMAX:    adjustable dimension of Y (NMAX.GE.N)
C        NF:      adjustable dimension of PE and SPE
C        MJ:      adjustable dimension of COV
C     @TEST.PN31:  5/16/89, 12/5/90, 12/21/90, 8/5/91, 7/19/92
C
cc      PARAMETER( NMAX=3000,NF=512,IDEV=1,JDEV=6,IWINDW=1 )
      PARAMETER( NF=1024 )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  Y(NMAX), PE(0:NF), SPE(0:NF)
cxx      DIMENSION  Y(N), PE(0:NF), SPE(0:NF)
C
      INTEGER N, IWINDW, NP, IFG
      DOUBLE PRECISION Y(N), PE(0:NF), SPE(0:NF)
C
cc      CALL  READTS( IDEV,Y,N )
C
      IF( IWINDW.EQ.0 )  LAG = N-1
cxx      IF( IWINDW.GT.0 )  LAG = 2*DSQRT( DFLOAT(N) )
cxxx      IF( IWINDW.GT.0 )  LAG = INT(2*DSQRT( DFLOAT(N) ))
      IF( IWINDW.GT.0 )  LAG = INT(2*DSQRT( DBLE(N) ))
      IF( IWINDW.EQ.0 )  NP  = (LAG+1)/2
      IF( IWINDW.GT.0 )  NP  = LAG
C
c----------
      NN = 0
c----------
      CALL  FFTPER( Y,N,NN,PE,NP )
cx      CALL  WINDOW( PE,NP,IWINDW,SPE )
      CALL  WINDOW( PE,NP,IWINDW,SPE,IFG )
cc      CALL  PRPER( JDEV,PE,SPE,N,NP,IWINDW )
C
cc      STOP
      RETURN
      E N D
      SUBROUTINE  FFTPER( Y,N,NN,P,N2 )
C
C ... periodogram computation by FFT ...
C
C     Inputs:
C        Y(I):   data
C        N:      data length
C        NN:     basic span for computing peiodogram (if NN>0)
C     Outputs:
C        P(J):   periodogram
C        NN:     basic span for computing peiodogram
C        N2:     number of frequencies
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION Y(N),P(0:1024),X(1024),FY(1024),WK(1024)
C
      INTEGER N, NN, N2
      DOUBLE PRECISION Y(N), P(0:1024)
c local
      DOUBLE PRECISION X(1024), FY(1024), WK(1024)
C
      IF(NN.LE.0) THEN
         IF( N.LE.1024 ) THEN
cxx           NP = ALOG( N-0.01 )/ALOG(2.0) + 1
           NP = INT(ALOG( N-0.01 )/ALOG(2.0) + 1)
           NN = 2**NP
           NPOOL = 1
         ELSE
           NN = 1024
           NPOOL = (N-1)/NN + 1
         END IF
      ELSE
cxx         NP = ALOG(NN-0.01)/ALOG(2.0) +1
         NP = INT(ALOG(NN-0.01)/ALOG(2.0) +1)
         NN = 2**NP
      IF( NN.GT.1024) NN=1024
         NPOOL = (N-1)/NN +1
      END IF
C
      N2 = NN/2
cxx      DO 10 I=0,N2
cxx   10 P(I) = 0.0D0
      P(0:N2) = 0.0D0
C
      DO 100 II=1,NPOOL
      ND = MIN(N,II*NN) - (II-1)*NN
      DO 20 I=1,ND
cxx   20 X(I) = Y((II-1)*NN+I)
      X(I) = Y((II-1)*NN+I)
   20 CONTINUE
      DO 30 I=ND+1,NN
cxx   30 X(I) = 0.0D0
      X(I) = 0.0D0
   30 CONTINUE
C
      CALL  FFTR2 ( X,NN,1,FY,WK )
C
      P(0)  = P(0) + FY(1)**2
      P(N2) = P(N2) + FY(N2+1)**2
      DO 40 I=0,N2-1
cxx   40 P(I) = P(I) + FY(I+1)**2 + FY(N2+1)**2
      P(I) = P(I) + FY(I+1)**2 + FY(N2+1)**2
   40 CONTINUE
  100 CONTINUE
C
      DO 50 I=0,N2
cxx   50 P(I) = P(I)/N 
      P(I) = P(I)/N
   50 CONTINUE

C
      RETURN
      E N D
      SUBROUTINE  FFTR2( X,N,ISW,FX,WRK )
C
C  ...  Fast Fourier Transform  ...
C
C     Inputs:
C        X(I):   data
C        N:      data length
C        ISW:
C        WRK:    working area
C     Output:
C        FX(J):  Fourier transform
C                  J=1,...,N/2+1   cosine transform
C                  J=N/2+2,...,N   sine transform
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  X(N), FX(N), WRK(N/4)
C
      INTEGER N, ISW
      DOUBLE PRECISION X(N), FX(N), WRK(N/4)
c local
      DOUBLE PRECISION PI2
C
      DATA  PI2 /6.2831853071796D0/
C
cxx      NP = DLOG( DFLOAT(N) )/DLOG( 2.0D0 ) + 1.0D-5
cxxx      NP = INT(DLOG( DFLOAT(N) )/DLOG( 2.0D0 ) + 1.0D-5)
      NP = INT(DLOG( DBLE(N) )/DLOG( 2.0D0 ) + 1.0D-5)
      N2 = N/2
      N4 = N/4
C
      IF( ISW.EQ.1 )  THEN
        DO 10 I=2,N4
cxx   10   WRK(I) = DSIN( (I-1)*PI2/N )
        WRK(I) = DSIN( (I-1)*PI2/N )
   10   CONTINUE
      END IF
C
      DO 20  I=1,N2
      FX(I)    = X(I) + X(I+N2)
cxx   20 FX(I+N2) = X(I) - X(I+N2)
      FX(I+N2) = X(I) - X(I+N2)
   20 CONTINUE
C
      IFG = 1
      JFG = 1
      K =  1
      M = N4
C
      DO 30  I=1,NP-1
      M2 = M*2
      K2 = K*2
      IF( M.GE.K )  THEN
        IF( IFG.LT.0 )  THEN
          CALL  FFTSB1( X,WRK,K,M,M2,K2,FX )
        ELSE
          CALL  FFTSB1( FX,WRK,K,M,M2,K2,X )
        END IF
      ELSE
        IF( JFG.EQ.1 )  THEN
          IF( IFG.LT.0 )  THEN
            CALL  FFTSB2( X,K2,M2,FX )
          ELSE
            CALL  FFTSB2( FX,K2,M2,X )
          END IF
          JFG = 2
        END IF
        IF( IFG.LT.0 )  THEN
          CALL  FFTSB3( FX,WRK,K,M,X )
        ELSE
          CALL  FFTSB3( X,WRK,K,M,FX )
        END IF
      END IF
      K = K*2
      M = M/2
cxx   30 IFG = -IFG
      IFG = -IFG
   30 CONTINUE
C
      IF( IFG.GT.0 )  THEN
         DO 40  I=1,N
cxx   40    FX(I) = X(I)
         FX(I) = X(I)
   40    CONTINUE
      END IF
C
      RETURN
      E N D
      SUBROUTINE  FFTSB2( X,M,L,Y )
C
C  ...  slave subroutine for FFTR2  ...
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  X(L,M),  Y(M,L)
C
      INTEGER M, L
      DOUBLE PRECISION X(L,M), Y(M,L)
C
      IF( M.GE.L )  THEN
cxx      DO 10  J=1,L
      DO 11  J=1,L
      DO 10  I=1,M
cxx   10 Y(I,J) = X(J,I)
      Y(I,J) = X(J,I)
   10 CONTINUE
   11 CONTINUE
      ELSE
cxx      DO 20  I=1,M
      DO 21  I=1,M
      DO 20  J=1,L
cxx   20 Y(I,J) = X(J,I)
      Y(I,J) = X(J,I)
   20 CONTINUE
   21 CONTINUE
      END IF
C
      RETURN
      E N D
      SUBROUTINE  FFTSB1( X,SINE,K,M,MJ1,MJ2,Y )
C
C  ...  slave subroutine for FFTR2  ...
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  X(MJ1,MJ2), Y(M,K,4), SINE(M,K)
C
      INTEGER K, M, MJ1, MJ2
      DOUBLE PRECISION X(MJ1,MJ2), SINE(M,K), Y(M,K,4)
c local
      DOUBLE PRECISION SUM1, SUM2
C
      DO 10  I=1,M
      Y(I,1,1) = X(I,1) + X(I+M,1)
      Y(I,1,3) = X(I,1) - X(I+M,1)
      Y(I,1,2) = X(I,K+1)
cxx   10 Y(I,1,4) = X(I+M,K+1)
      Y(I,1,4) = X(I+M,K+1)
   10 CONTINUE
C
cxx      DO 20  J=2,K
      DO 21  J=2,K
      DO 20  I=1,M
      SUM1 = SINE(1,K-J+2)*X(I+M,J) - SINE(1,J)    *X(I+M,J+K)
      SUM2 = SINE(1,J)    *X(I+M,J) + SINE(1,K-J+2)*X(I+M,J+K)
      Y(I,J,1)     = X(I,J) + SUM1
      Y(I,K-J+2,2) = X(I,J) - SUM1
      Y(I,J,3)     = SUM2 + X(I,J+K)
cxx   20 Y(I,K-J+2,4) = SUM2 - X(I,J+K)
      Y(I,K-J+2,4) = SUM2 - X(I,J+K)
   20 CONTINUE
   21 CONTINUE
C
      RETURN
      E N D
      SUBROUTINE  FFTSB3( X,SINE,K,M,Y )
C
C  ...  slave subroutine for FFTR2  ...
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  X(K,2,M,2),  SINE(M,K),  Y(K,4,M)
C
      INTEGER K, M
      DOUBLE PRECISION X(K,2,M,2), SINE(M,K), Y(K,4,M)
c local
      DOUBLE PRECISION SUM1, SUM2
C
cxx      DO 10  J=1,M
      DO 11  J=1,M
      Y(1,1,J) = X(1,1,J,1) + X(1,1,J,2)
      Y(1,3,J) = X(1,1,J,1) - X(1,1,J,2)
      Y(1,2,J) = X(1,2,J,1)
      Y(1,4,J) = X(1,2,J,2)
      DO 10  I=2,K
      SUM1 = SINE(1,K-I+2)*X(I,1,J,2) - SINE(1,I)    *X(I,2,J,2)
      SUM2 = SINE(1,I)    *X(I,1,J,2) + SINE(1,K-I+2)*X(I,2,J,2)
      Y(I,1,J)     = X(I,1,J,1) + SUM1
      Y(K-I+2,2,J) = X(I,1,J,1) - SUM1
      Y(I,3,J)     = SUM2 + X(I,2,J,1)
cxx   10 Y(K-I+2,4,J) = SUM2 - X(I,2,J,1)
      Y(K-I+2,4,J) = SUM2 - X(I,2,J,1)
   10 CONTINUE
   11 CONTINUE
C
      RETURN
      E N D
