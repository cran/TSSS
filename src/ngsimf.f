C     PROGRAM 15.2  NGSIM
      SUBROUTINE NGSIMF( M1,M2,M3,M,K,N,INI0,NOISEW,WMIN,WMAX,
     &                   PARAMW,NOISEV,VMIN,VMAX,PARAMV,PERIOD,AR,X,Y )
C
      INCLUDE 'TSSS.h'
C
C  ...  SIMULATION BY NON-GAUSSIAN STATE SPACE MODEL  ...
C
C     INPUTS:
C        M1:   TREND ORDER (0,1,2)
C        M2:   SEASONAL ORDER (0 OR 1)
C        M3:   AR ORDER (0,...,20)
C        N:    SIMULATION INTERVAL
C        NS:   Y(NS+1),...,Y(N) ARE ACTUALLY USED
C        INI:  INITIAL VALUE FOR RANDOM NUMBER GENERATOR
C        NOISEW:    TYPE OF OBSERVATIONAL NOISE
C        NOISEV:    TYPE OF SYSTEM NOISE
C        WMIN,WMAX: LOWER AND UPPER BOUND OF OBS. NOISE
C        VMIN,VMAX: LOWER AND UPPER BOUND OF SYSTEM NOISE
C        X(I): INITIAL STATE
C     PARAMETERS:
C        MJ:   MAXIMUM STATE DIMENSION
C        NN:   MAXIMUM DATA LENGTH
C     MODIFIED  2/16/93
C
cc      PARAMETER( MJ=20,NN=1000 )
      PARAMETER( L=1 )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      INTEGER    PERIOD
cc      DIMENSION  X(MJ), F(MJ,MJ), G(MJ,MJ), H(MJ)
cc      DIMENSION  Q(MJ,MJ), Y(NN), AR(20), R(1,1)
cxx      DIMENSION  X(M), F(M,M), G(M,K), H(M)
cxx      DIMENSION  Q(K,K), Y(N), AR(M3), R(L,L)
cxx      DIMENSION  PARAMV(3), PARAMW(3)
C
      INTEGER M1, M2, M3, M, K, N, INI0, NOISEW, NOISEV, PERIOD
      DOUBLE PRECISION WMIN, WMAX, PARAMW(3), VMIN, VMAX, PARAMV(3),
     1                 AR(M3), X(M), Y(N)
c local
      DOUBLE PRECISION F(M,M), G(M,K), H(M), Q(K,K), R(L,L), TAU1, TAU2,
     1                 TAU3, SIG2
C
      DATA   TAU1/1.0D0/, TAU2/1.0D0/, TAU3/1.0D0/, SIG2/1.0D0/
cc      INI = 1992092521
C
cc      READ( 5,* )  M1, M2, M3, N, NS, INI
C
cc      READ( 5,* )  NOISEW, WMIN, WMAX
      PARAMW(1) = 0.0D0
cc      IF( NOISEW.NE.2 )  READ(5,*)  PARAMW(2)
cc      IF( NOISEW.EQ.2 )  READ(5,*)  PARAMW(2), PARAMW(3)
C
cc      READ( 5,* )  NOISEV, VMIN, VMAX
      PARAMV(1) = 0.0D0
cc      IF( NOISEV.NE.2 )  READ(5,*)  PARAMV(2)
cc      IF( NOISEV.EQ.2 )  READ(5,*)  PARAMV(2), PARAMV(3)
C
cc      MM = 0
cc      IF( M1.GT.0 )  THEN
cc         READ( 5,* )  (X(I),I=1,M1)
cc         MM = M1
cc      END IF
cc      IF( M2.GT.0 )  THEN
cc         READ( 5,* )  PERIOD
cc         READ( 5,* )  (X(MM+I),I=1,PERIOD-1)
cc         MM = MM + PERIOD-1
cc      END IF
cc      IF( M3.GT.0 )  THEN
cc         READ( 5,* )  (AR(I),I=1,M3)
cc         READ( 5,* )  (X(MM+I),I=1,M3)
cc      END IF
      INI  = INI0
C
C  ...  SET STATE SPACE MODEL  ...
C
      CALL  SETSEA( M1,M2,M3,PERIOD,AR,TAU1,TAU2,TAU3,SIG2,F,G,H,
cc     *              Q,R,M,K,MJ )
     *              Q,R,M,K )
C
C  ...  SIMULATION  ...
C
      CALL  NGSIM( NOISEV,NOISEW,PARAMV,PARAMW,VMIN,VMAX,WMIN,WMAX,
cxx     *             F,G,H,Q,R,X,N,M,1,K,INI,MJ,NN,Y )
     *             F,G,H,Q,R,X,N,M,1,K,INI,Y )
C
C  ...  PLOT AND PRINT OUT SIMULATED TIME SERIES  ...
C
C      CALL  PTSIM( Y,N,NS,K,M1,M2,M3,INI0 )
cc      CALL  PRSIM( Y,N,NS,M1,M2,M3,INI0 )
C
cc      STOP
      RETURN
      E N D
      SUBROUTINE NGSIM( NOISEV,NOISEW,PV,PW,X0,X1,Y0,Y1,F,G,H,Q,R,
cxx     *                  X,N,M,L,K,IX,MJ,NN,Y )
     *                  X,N,M,L,K,IX,Y )
C
C  ...  SIMULATION BY NON-GAUSSIAN STATE SPACE MODEL  ...
C
C     INPUTS:
C        NOISEV:  TYPE OF THE SYSTEM NOISE
C        NOISEW:  TYPE OF THE OBSERVATIONAL NOISE
C        PV:      PARAMETER OF THE SYSTEM NOISE DENSITY
C        PW:      PARAMETER OF THE OBSER. NOISE DENSITY
C        X0,X1:   DOMAIN OF SYSTEM NOISE DENSITY
C        Y0,Y1:   DOMAIN OF OBSER. NOISE DENSITY
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
cxx      DIMENSION XT(M), T(L), V(K), W(L), Y(N,L)
cxx      DIMENSION FV(0:400), FW(0:400), XV(0:400), XW(0:400)
cxx      DIMENSION PV(3), PW(3)
C
      INTEGER NOISEV, NOISEW, N, M, L, K, IX
      DOUBLE PRECISION PV(3), PW(3), X0, X1, Y0, Y1, F(M,M), G(M,K),
     1                 H(L,M), Q(K,K), R(L,L), X(M), Y(N,L)
c local
      DOUBLE PRECISION GAUSS, PEARSN, DBLEXP, USERV1, USERW1, SQ(K,K),
     1                 SR(L,L), XT(M), T(L), V(K), W(L), FV(0:400),
     2                 FW(0:400), XV(0:400), XW(0:400), DXV, DXW
C
      EXTERNAL  GAUSS
      EXTERNAL  PEARSN
      EXTERNAL  DBLEXP
cc      EXTERNAL  USERV
cc      EXTERNAL  USERW
      EXTERNAL  USERV1
      EXTERNAL  USERW1
C
C  ...  DISTRIBUTION FUNCTION OF SYSTEM NOISE AND OBS. NOISE  ...
C
cc      IF(NOISEV.EQ.0)  CALL  DISTRI( USERV,PV,X0,X1,FV,XV,DXV )
cxx      IF(NOISEV.EQ.0)  CALL  DISTRI( USERV1,PV,X0,X1,FV,XV,DXV )
      IF(NOISEV.EQ.0)  CALL  DISTRI0( USERV1,X0,X1,FV,XV,DXV )
      IF(NOISEV.EQ.1)  CALL  DISTRI( GAUSS ,PV,X0,X1,FV,XV,DXV )
      IF(NOISEV.EQ.2)  CALL  DISTRI( PEARSN,PV,X0,X1,FV,XV,DXV )
      IF(NOISEV.EQ.3)  CALL  DISTRI( DBLEXP,PV,X0,X1,FV,XV,DXV )
C
cc      IF(NOISEW.EQ.0)  CALL  DISTRI( USERW,PW,Y0,Y1,FW,XW,DXW )
cxx      IF(NOISEW.EQ.0)  CALL  DISTRI( USERW1,PW,Y0,Y1,FW,XW,DXW )
      IF(NOISEW.EQ.0)  CALL  DISTRI0( USERW1,Y0,Y1,FW,XW,DXW )
      IF(NOISEW.EQ.1)  CALL  DISTRI( GAUSS ,PW,Y0,Y1,FW,XW,DXW )
      IF(NOISEW.EQ.2)  CALL  DISTRI( PEARSN,PW,Y0,Y1,FW,XW,DXW )
      IF(NOISEW.EQ.3)  CALL  DISTRI( DBLEXP,PW,Y0,Y1,FW,XW,DXW )
C
cc      CALL  CHOLES( Q,MJ,K,SQ,10 )
cc      CALL  CHOLES( R,L,L,SR,10 )
      CALL  CHOLES( Q,K,K,SQ,K )
      CALL  CHOLES( R,L,L,SR,L )
C
C----- Initialize Mersenne Twister with given seed value
      call init(IX)
c-----
C
      DO 100 II= 1,N
C
C  ...  X = F*X + V  ...
C
ccxx      CALL  NGNOIS( NOISEV,FV,XV,DXV,IX,SQ,K,V )
      CALL  NGNOIS( NOISEV,FV,XV,DXV,SQ,K,V )
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
ccxx      CALL  NGNOIS( NOISEW,FW,XW,DXW,IX,SR,L,W )
      CALL  NGNOIS( NOISEW,FW,XW,DXW,SR,L,W )
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
ccxx      SUBROUTINE  NGNOIS( NOISE,F,X,DX,IX,Q,K,V )
      SUBROUTINE  NGNOIS( NOISE,F,X,DX,Q,K,V )
C
C  ...  K-DIMENSIONAL NON-GAUSSIAN RANDOM NUMBER  ...
C
C     INPUTS:
C        NOISE:NOISE TYPE
C        F:    DISTRIBUTION FUNCTION OF NON-GAUSSIAN DISTRIBUTION
C        IX:   INITIAL SEED FOR RANDOM NUMBER GENERATOR
C        Q:    COVARIANCE
C        K:    DIMENSION OF THE RANDOM NUMBER
C     OUTPUT:
C        V:    K-DIMENSIONAL RANDOM NUMBER WITH COVARIANCE Q
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  Q(10,10), W(10), V(K), F(0:400), X(0:400)
cxx      DIMENSION  Q(K,K), W(K), V(K), F(0:400), X(0:400)
C
      INTEGER NOISE, K
      DOUBLE PRECISION F(0:400), X(0:400), DX, Q(K,K), V(K)
c local
      DOUBLE PRECISION RNG, RNG2, W(K), SUM
C
      DO 10 I=1,K
ccxx      IF( NOISE.LT.0 )  W(I) = RNG2( IX,NOISE )
      IF( NOISE.LT.0 )  W(I) = RNG2( NOISE )
cxx   10 IF( NOISE.GE.0 )  W(I) = RNG( IX,F,X,DX )
ccxx      IF( NOISE.GE.0 )  W(I) = RNG( IX,F,X,DX )
      IF( NOISE.GE.0 )  W(I) = RNG( F,X,DX )
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
      SUBROUTINE  DISTRI( FUNCT,PARAM,XMIN,XMAX,F,X,DX )
C
C  ...  DISTRIBUTION FUNCTION  ...
C
C     INPUTS:
C        FUNCT:  FUNCTION NAME
C        PARAM:  PARAMETER OF THE DENSITY
C        XMIN,XMAX:  DOMAIN OF THE  DENSITY
C     OUTPUTS:
C        F:      DISTRIBUTION FUNCTION
C        X:      LOCATION OF NODE
C        DX:     SIZE OF EACH BIN
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  F(0:400), P(0:400), X(0:400), PARAM(3)
C
      DOUBLE PRECISION FUNCT, PARAM(3), XMIN, XMAX, F(0:400), X(0:400),
     1                 DX
c local
      DOUBLE PRECISION P(0:400)
      EXTERNAL FUNCT
C
      K = 400
      DX = (XMAX-XMIN)/K
C
      DO 10 I=0,K
      X(I) = XMIN + I*DX
cxx   10 P(I) = FUNCT( X(I),PARAM )
      P(I) = FUNCT( X(I),PARAM )
   10 CONTINUE
cxx      DO 20 I=0,K
cxx   20 F(I) = 0.0D0
      F(0:K) = 0.0D0
      DO 30 I=1,K
cxx   30 F(I) = F(I-1) + ( P(I-1) + P(I) )*DX/2.0D0
      F(I) = F(I-1) + ( P(I-1) + P(I) )*DX/2.0D0
   30 CONTINUE
      DO 40 I=1,K
cxx   40 F(I) = F(I)/F(K)
      F(I) = F(I)/F(K)
   40 CONTINUE
C
      RETURN
      E N D
cxx      SUBROUTINE  DISTRI( FUNCT,PARAM,XMIN,XMAX,F,X,DX )
      SUBROUTINE  DISTRI0( FUNCT,XMIN,XMAX,F,X,DX )
C
C  ...  DISTRIBUTION FUNCTION  ...
C
C     INPUTS:
C        FUNCT:  FUNCTION NAME
C        PARAM:  PARAMETER OF THE DENSITY
C        XMIN,XMAX:  DOMAIN OF THE  DENSITY
C     OUTPUTS:
C        F:      DISTRIBUTION FUNCTION
C        X:      LOCATION OF NODE
C        DX:     SIZE OF EACH BIN
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  F(0:400), P(0:400), X(0:400), PARAM(3)
C
      DOUBLE PRECISION FUNCT, XMIN, XMAX, F(0:400), X(0:400), DX
c local
      DOUBLE PRECISION P(0:400)
      EXTERNAL FUNCT
C
      K = 400
      DX = (XMAX-XMIN)/K
C
      DO 10 I=0,K
      X(I) = XMIN + I*DX
cxx   10 P(I) = FUNCT( X(I),PARAM )
      P(I) = FUNCT( X(I) )
   10 CONTINUE
cxx      DO 20 I=0,K
cxx   20 F(I) = 0.0D0
      F(0:K) = 0.0D0
      DO 30 I=1,K
cxx   30 F(I) = F(I-1) + ( P(I-1) + P(I) )*DX/2.0D0
      F(I) = F(I-1) + ( P(I-1) + P(I) )*DX/2.0D0
   30 CONTINUE
      DO 40 I=1,K
cxx   40 F(I) = F(I)/F(K)
      F(I) = F(I)/F(K)
   40 CONTINUE
C
      RETURN
      E N D
ccxx      DOUBLE PRECISION FUNCTION  RNG( IX,F,X,DX )
      DOUBLE PRECISION FUNCTION  RNG( F,X,DX )
C
C  ...  NON-GAUSSIAN RANDON NUMBER GENERATOR  ...
C
C     INPUTS:
C        IX:     SEED
C        F:      DISTRIBUTION FUNCTION
C        X:      LOCATION
C        DX:     SIZE OF EACH BIN
cxx      IMPLICIT REAL*8( A-H,O-Z)
cxx      DIMENSION  F(0:400), X(0:400)
C
ccxx      INTEGER IX
      DOUBLE PRECISION F(0:400), X(0:400), DX
c local
cxx      DOUBLE PRECISION RUNI, U, V
      DOUBLE PRECISION U, V, random

C
cxx      IF( IX.EQ.0 )  IX = 1990103011
cxx      U = RUNI(IX)
      U = random()
C
      I = 0
   10 I = I+1
      IF( U.GT.F(I) )  GO TO 10
      IF( U.EQ.F(I) )  THEN
        V = X(I)
      ELSE
        V =((U-F(I-1))/( F(I)-F(I-1)) )*DX+X(I-1)
      END IF
      RNG = V
C
      RETURN
      E N D
ccxx      DOUBLE PRECISION FUNCTION  RNG2( IX,NOISE )
      DOUBLE PRECISION FUNCTION RNG2( NOISE )
C
C  ...  NON-GAUSSIAN RANDON NUMBER GENERATOR  ...
C
C     NOISE = -1:  CAUCHY RANDOM NUMBER
C           = -2:  EXPONENTIAL DISTRIBUTION
C           = -3:  DOUBLE EXPONENTIAL DISTRIBUTION
C
cxx      IMPLICIT REAL*8( A-H,O-Z)
C
ccxx      INTEGER IX, NOISE
cxx      DOUBLE PRECISION RUNI, PI, U
      DOUBLE PRECISION PI, U, random
C
      DATA  PI /3.1415926535D0/
C
cxx      IF( IX.EQ.0 )  IX = 1990103011
cxx   10 U = RUNI()
      U = random()
C
C  ...  CAUCHY DISTRIBUTION  ...
C
ccxx      IF( NOISE.EQ.-1 )  THEN
ccxx         IF( U.EQ.0.5D0 )  GO TO 10
         RNG2 = DTAN( PI*U )
ccxx      END IF
C
C  ...  EXPONENTIAL DISTRIBUTION  ...
C
      IF( NOISE.EQ.-2 )  THEN
ccxx         IF( U.EQ.0.0D0 )  GO TO 10
         RNG2 = -DLOG( U )
      END IF
C
C  ...  DOUBLE EXPONENTIAL DISTRIBUTION  ...
C
      IF( NOISE.EQ.-3 )  THEN
ccxx         IF( U.EQ.0.0D0 )  GO TO 10
         RNG2 = DEXP( -DEXP( U ) )
      END IF
C
      RETURN
      E N D
cc      DOUBLE PRECISION FUNCTION  USERV( X,PARAM )
cxx      DOUBLE PRECISION FUNCTION  USERV1( X,PARAM )
      DOUBLE PRECISION FUNCTION  USERV1( X )
C
C  ...  User supplied density function  ...
C       (The double exp. dist. with mean 0 )
C     Inputs:
C        X:
C        PARAM(1):  lambda
C     Output:
C        USERV:     density at X
C
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cxx      DIMENSION  PARAM(2)
C
      DOUBLE PRECISION X, C
      DATA  C/0.5771D0/
C
cc      USERV = DEXP( X-C - DEXP(X-C) )
      USERV1 = DEXP( X-C - DEXP(X-C) )
      RETURN
      E N D
cc      DOUBLE PRECISION FUNCTION  USERW( X,PARAM )
cxx      DOUBLE PRECISION FUNCTION  USERW1( X,PARAM )
      DOUBLE PRECISION FUNCTION  USERW1( X )
C
C  ...  User supplied density function  ...
C       (The double exp. dist. with mean 0 )
C     Inputs:
C        X:
C        PARAM(1):  lambda
C     Output:
C        USERW:     density at X
C
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cxx      DIMENSION  PARAM(2)
C
      DOUBLE PRECISION X, C
      DATA  C/0.5771D0/
C
cc      USERW = DEXP( X-C - DEXP(X-C) )
      USERW1 = DEXP( X-C - DEXP(X-C) )
      RETURN
      E N D
