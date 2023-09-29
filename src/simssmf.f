C     PROGRAM 15.1  SIMSSM
      SUBROUTINE SIMSSMF( M1,M2,M3,M,K,N,INI,SIG2,PERIOD,TAU1,TAU2,TAU3,
     &                    AR,X,Y )
C
      INCLUDE 'TSSS.h'
C
C  ...  SIMULATION BY GAUSSIAN STATE SPACE MODEL  ...
C
C     INPUTS:
C        M1:   TREND ORDER (0,1,2)
C        M2:   SEASONAL ORDER (0 OR 1)
C        M3:   AR ORDER (0,...,20)
C        N:    SIMULATION INTERVAL
C        NS:   Y(NS+1),...,Y(N) ARE ACTUALLY USED
C        INI:  INITIAL VALUE FOR RANDOM NUMBER GENERATOR
C     PARAMETERS:
C        MJ:   MAXIMUM STATE DIMENSION
C        NN:   MAXIMUM DATA LENGTH
C     MODIFED  2/16/93
C
cc      PARAMETER( MJ=20,NN=1000 )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      INTEGER    PERIOD
cc      DIMENSION  X(MJ), F(MJ,MJ), G(MJ,MJ), H(MJ)
cc      DIMENSION  Q(MJ,MJ), Y(NN), AR(20), R(1,1)
cc      INI = 1992092521
cxx      DIMENSION  X(M), F(M,M), G(M,K), H(M)
cxx      DIMENSION  Q(K,K), Y(N), AR(M3), R(1,1)
C
      INTEGER M1, M2, M3, M, K, N, INI, PERIOD
      DOUBLE PRECISION SIG2, TAU1, TAU2, TAU3, AR(M3), X(M), Y(N)
c local
      INTEGER INI0
      DOUBLE PRECISION F(M,M), G(M,K), H(M), Q(K,K), R(1,1)
C
cc      MM = 0
C
cc      READ( 5,* )  M1, M2, M3, N, NS, INI
cc      READ( 5,* )  SIG2
cc      IF( M1.GT.0 )  THEN
cc         READ( 5,* )  TAU1
cc         READ( 5,* )  (X(I),I=1,M1)
cc         MM = M1
cc      END IF
cc      IF( M2.GT.0 )  THEN
cc         READ( 5,* )  TAU2, PERIOD
cc         READ( 5,* )  (X(MM+I),I=1,PERIOD-1)
cc         MM = MM + PERIOD-1
cc      END IF
cc      IF( M3.GT.0 )  THEN
cc         READ( 5,* )  TAU3, (AR(I),I=1,M3)
cc         READ( 5,* )  (X(MM+I),I=1,M3)
cc      END IF
      INI0  = INI
C
C  ...  SET STATE SPACE MODEL  ...
C
      CALL  SETSEA( M1,M2,M3,PERIOD,AR,TAU1,TAU2,TAU3,SIG2,F,G,H,
cc     *              Q,R,M,K,MJ )
     *              Q,R,M,K )
C
C  ...  SIMULATION  ...
C
cc      CALL  SIMSSM( F,G,H,Q,R,X,N,M,1,K,INI,MJ,NN,Y )
      CALL  SIMSSM( F,G,H,Q,R,X,N,M,1,K,INI,N,Y )
C
C  ...  PLOT AND PRINT OUT SIMULATED TIME SERIES  ...
C
C      CALL  PTSIM( Y,N,NS,K,M1,M2,M3,INI0 )
cc      CALL  PRSIM( Y,N,NS,M1,M2,M3,INI0 )
C
cc      STOP
      RETURN
      E N D
cc      SUBROUTINE SIMSSM( F,G,H,Q,R,X,N,M,L,K,IX,MJ,NN,Y )
      SUBROUTINE SIMSSM( F,G,H,Q,R,X,N,M,L,K,IX,NN,Y )
C
C  ...  SIMULATION BY STATE SPACE MODEL  ...
C
C     INPUTS:
C        F,G,H:   COEFFICIENT MATRICES OF STATE SPACE MODEL
C        Q,R:     SYSTEM NOISE AND OBS. NOISE VARIANCES
C        X:       INITIAL STATE VECTOR
C        M,L,K:   DIMENSIONS OF THE STATE, SYSTEM NOISE AND OBS.
C        IX:      INITAL SEED FOR RANDOM NUMBER GENERATOR
C        MJ:      ADJUSTABLE DIMENSION
C        NN:      NUMBER OF DATA
C     OUTPUT:
C        Y(I):    SIMULATED TIME SERIES
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION X(MJ), F(MJ,MJ), G(MJ,MJ), H(L,MJ)
cc      DIMENSION Q(MJ,MJ), R(L,L), SQ(10,10), SR(10,10)
cc      DIMENSION XT(40), T(20), V(20), W(10), Y(NN,L)
cxx      DIMENSION X(M), F(M,M), G(M,K), H(L,M)
cxx      DIMENSION Q(K,K), R(L,L), SQ(K,K), SR(L,L)
cxx      DIMENSION XT(M), T(L), V(K), W(L), Y(NN,L)
C
      INTEGER N, M, L, K, IX, NN
      DOUBLE PRECISION F(M,M), G(M,K), H(L,M), Q(K,K), R(L,L), X(M),
     1                 Y(NN,L)
c local
      INTEGER I, II, J
      DOUBLE PRECISION SQ(K,K), SR(L,L), XT(M), T(L), V(K), W(L)
C
C  ...  CHOLESKY DECOMPOSITION OF Q AND R  ...
C
cc      CALL  CHOLES( Q,MJ,K,SQ,10 )
cc      CALL  CHOLES( R,L,L,SR,10 )
      CALL  CHOLES( Q,K,K,SQ,K )
      CALL  CHOLES( R,L,L,SR,L )
C
C----- Initialize Mersenne Twister with given seed value
      call init(IX)
c-----
      DO 100 II= 1,N
C
C  ...  X = F*X + V  ...
C
ccxx      CALL  WHITE( K,SQ,IX,V )
      CALL  WHITE( K,SQ,V )
C
cxx      DO 10 I=1,M
cxx   10 XT(I) = 0.0D0
      XT(1:M) = 0.0D0
cxx      DO 20 J=1,M
      DO 21 J=1,M
      DO 20 I=1,M
cxx   20 XT(I) = XT(I) + F(I,J)*X(J)
      XT(I) = XT(I) + F(I,J)*X(J)
   20 CONTINUE
   21 CONTINUE
cxx      DO 30 J=1,K
      DO 31 J=1,K
      DO 30 I=1,M
cxx   30 XT(I) = XT(I) + G(I,J)*V(J)
      XT(I) = XT(I) + G(I,J)*V(J)
   30 CONTINUE
   31 CONTINUE
      DO 40 I=1,M
cxx   40 X(I) = XT(I)
      X(I) = XT(I)
   40 CONTINUE
C
C  ...  T = H*X + W  ...
C
ccxx      CALL  WHITE( L,SR,IX,W )
      CALL  WHITE( L,SR,W )
cxx      DO 50 I=1,L
cxx   50 T(I) = 0.0D0
      T(1:L) = 0.0D0
cxx      DO 60 I=1,L
      DO 61 I=1,L
      DO 60 J=1,M
cxx   60 T(I) = T(I) + H(I,J)*X(J)
      T(I) = T(I) + H(I,J)*X(J)
   60 CONTINUE
   61 CONTINUE
      DO 80 I=1,L
cxx   80 Y(II,I) = T(I) + W(I)
      Y(II,I) = T(I) + W(I)
   80 CONTINUE
  100 CONTINUE
C
      RETURN
      E N D
ccxx      SUBROUTINE  WHITE( K,Q,IX,V )
      SUBROUTINE  WHITE( K,Q,V )
C
C  ...  GENERATE K-DIMENSIONAL GAUSSIAN RANDOM NUMBER  ...
C
C     INPUTS:
C        K:   DIMENSION OF THE RANDOM NUMBER
C        Q:   COVARIANCE
C        IX:  INITIAL SEED FOR RANDOM NUMBER GENERATOR
C     OUTPUT:
C        V:   K-DIMENSIONAL RANDOM NUMBER WITH COVARIANCE Q
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION Q(K,K),W(100),V(K)
cxx      DIMENSION Q(K,K),W(K),V(K)
C
      INTEGER K
      DOUBLE PRECISION Q(K,K), V(K)
c local
      INTEGER I, IC, J
      DOUBLE PRECISION RGAUSS, W(K), V2, S, SUM
C
      IC = 0
      DO 10 I=1,K
cc   10 W(I) = RGAUSS( IX )
cxx   10 W(I) = RGAUSS( IX,IC,V2,S )
ccxx      W(I) = RGAUSS( IX,IC,V2,S )
      W(I) = RGAUSS( IC,V2,S )
   10 CONTINUE
      DO 30 I=1,K
      SUM = 0.0D0
      DO 20 J=1,K
cxx   20 SUM = SUM + Q(J,I)*W(J)
      SUM = SUM + Q(J,I)*W(J)
   20 CONTINUE
cxx   30 V(I) = SUM
      V(I) = SUM
   30 CONTINUE
      RETURN
      E N D
cc      DOUBLE PRECISION FUNCTION  RGAUSS( IX )
ccxx      DOUBLE PRECISION FUNCTION  RGAUSS( IX,ICOUNT,V2,S )
      DOUBLE PRECISION FUNCTION  RGAUSS( ICOUNT,V2,S )
C
C  ...  GAUSSIAN RANDON NUMBER GENERATOR (MARSAGLIA'S ALGORITHM)  ...
C
cxx      IMPLICIT REAL*8( A-H,O-Z)
      INTEGER ICOUNT
      DOUBLE PRECISION V2, S
c local
cxx      DOUBLE PRECISION RUNI, U1, U2, V1
      DOUBLE PRECISION U1, U2, V1, random

cc      DATA  ICOUNT /0/
C
cxx      IF( IX.EQ.0 )  IX = 1990103011
      IF( ICOUNT.EQ.0 )  THEN
cxx   10 U1 = RUNI(IX)
cxx      U2 = RUNI(IX)
   10 U1 = random()
      U2 = random()
      V1 = 2*U1-1
      V2 = 2*U2-1
      S = V1**2 + V2**2
      IF(S.GE.1) GO TO 10
      RGAUSS = V1*SQRT( -2*DLOG(S)/S )
      ICOUNT = 1
C
      ELSE
        RGAUSS = V2*SQRT( -2*DLOG(S)/S)
        ICOUNT = 0
      END IF
      RETURN
      E N D
