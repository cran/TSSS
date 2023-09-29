C     PROGRAM 14.1  NGSMTH
      SUBROUTINE  NGSMTHF( Y,N,NOISEV,TAU2,BV,NOISEW,SIG2,BW,INITD,
     *                    TREND,SS,FF,NS,NFE,NPE,K )
C
      INCLUDE 'TSSS.h'
C
C  ...  NON-GAUSSIAN SMOOTHING  ...
C
C     INPUTS:
C        NOISEV:   TYPE OF SYSTEM NOISE DENSITY (0,1,2,3)
C        NOISEW:   TYPE OF OBSER. NOISE DENSITY (0,1,2,3,4)
C        TAU2:     VARIANCE OR DISPERSION OF SYSTEM NOISE
C        SIG2:     VARIANCE OR DISPERSION OF OBSERVATION NOISE
C        BV:       SHAPE PARAMETER OF SYSTEM NOISE (FOR NOISEV = 2)
C        BW:       SHAPE PARAMETER OF OBSER. NOISE (FOR NOISEW = 2)
C     PARAMETERS:
C        MJ:       ADJUSTABLE DIMENSION OF Y, PS, FS, SS (MJ >= N)
C        K:        NUMBER OF INTERVALS + 1
C     @NFILTER.F:   5/14/85, 6/11/87, 5/19/88, 2/14/91
C     MODIFIED  2/16/93
C
cc      PARAMETER( MJ=500, K=201 )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      CHARACTER  TITLE*72
cc      DIMENSION  Y(MJ), P(K), F(K), S(K), T(K), Q(-K:K)
cc      DIMENSION  YM(10), ST(7), TREND(MJ,7), LOC(MJ)
cc      REAL*4    PS(K,MJ), FS(K,MJ), SS(K,MJ)
cxx      DIMENSION  Y(N), F(K), Q(-K:K)
cxx      DIMENSION  YM(2), ST(7), TREND(NPE,7), LOC(NPE)
cxx      REAL*4     SS(K,NPE)
cxx      REAL*8     SS(K,NPE)
C
      INTEGER N, NOISEV, NOISEW, INITD, NS, NFE, NPE, K
      DOUBLE PRECISION Y(N), TAU2, BV, SIG2, BW, TREND(NPE,7),
     1                 SS(K,NPE), FF
c local
      INTEGER I, J, LOC(NPE)
      DOUBLE PRECISION F(K), Q(-K:K), YM(2), ST(7), DX, XMIN, XMAX,
     1                 OUTMIN, OUTMAX, SSUM
C
cc      COMMON   /COMAAA/  Y
cc      COMMON   /C91214/  XMIN, XMAX, FIGMIN, FIGMAX, YM
cc      COMMON   /C91215/  NOISEV, NOISEW, INITD, ITF, ITH
cc      COMMON   /C91216/  B, OUTMIN, OUTMAX
cc      COMMON   /C91218/  SIG2, TAU2, BW, BV
cc      COMMON   /C91219/  N, KK, NS, NFE, NPE
cc      COMMON   /CMDATA/  TITLE
cc      COMMON    / DDD /  FF , AIC , SD
C     NAMELIST  /PARAM/  XMIN, XMAX, TAU20, SIG20, BC, IOPT, B,
C    *                   OUTMIN, OUTMAX, NS, NFE, NPE, FIGMIN, FIGMAX,
C    *                   NOISEV, NOISEW, ITF, ITH
C
cc      EXTERNAL  GAUSS
cc      EXTERNAL  PEARSN
cc      EXTERNAL  TWOEXP
cc      EXTERNAL  USERV
C
C  ...  READ DATA  ...
C 
cc      CALL  READTS( 1,Y,N )
      CALL  MOMENT( Y,N,YM(1),YM(2) )
C
      LOC(1) = 0
cc      KK = K
cc      CALL  DEFALT
      CALL  DEFALT( Y,N,XMIN,XMAX,OUTMIN,OUTMAX )
cc      READ( 5,* )  NOISEV, TAU2, BV
cc      READ( 5,* )  NOISEW, SIG2, BW
C
      DX = (XMAX-XMIN)/(K-1)
cc      CALL  IDIST( F,K,YM(1),YM(2),XMIN,DX )
      CALL  IDIST( F,K,YM(1),YM(2),XMIN,DX,INITD )
      CALL  NORMLZ( F,K,DX,SSUM )
cc      IF( NOISEV.EQ.0 )  CALL  TRANS( USERV ,K,DX,TAU2,BV,Q )
cc      IF( NOISEV.EQ.1 )  CALL  TRANS( GAUSS ,K,DX,TAU2,BV,Q )
cc      IF( NOISEV.EQ.2 )  CALL  TRANS( PEARSN,K,DX,TAU2,BV,Q )
cc      IF( NOISEV.EQ.3 )  CALL  TRANS( TWOEXP,K,DX,TAU2,BV,Q )
      IF( NOISEV.EQ.0 )  CALL  TRANS1( K,DX,TAU2,BV,Q )
      IF( NOISEV.EQ.1 )  CALL  TRANS2( K,DX,TAU2,BV,Q )
      IF( NOISEV.EQ.2 )  CALL  TRANS3( K,DX,TAU2,BV,Q )
      IF( NOISEV.EQ.3 )  CALL  TRANS4( K,DX,TAU2,BV,Q )
cc      CALL  NGSMTH( Y,P,F,S,T,N,K,DX,XMIN,Q,FF,PS,FS,SS,LOC )
      CALL  NGSMTH( NOISEW,SIG2,BW,Y,F,N,K,DX,XMIN,Q,FF,SS,LOC,
     *              OUTMIN,OUTMAX,NS,NFE,NPE )
C
      DO 30 I=1,NPE
      DO 10 J=1,K
cc   10 S(J) = SS(J,I)
cc      CALL  PINTVL( S,K,XMIN,DX,ST )
cxx   10 F(J) = DBLE(SS(J,I))
cxx   10 F(J) = SS(J,I)
      F(J) = SS(J,I)
   10 CONTINUE
      CALL  PINTVL( F,K,XMIN,DX,ST )
      DO 20 J=1,7
cxx   20 TREND(I,J) = ST(J) + DX*LOC(I)
      TREND(I,J) = ST(J) + DX*LOC(I)
   20 CONTINUE
   30 CONTINUE
C
cc      CALL  PRNGSM( NOISEV,NOISEW,TAU2,SIG2,FF,TREND,MJ,N )
C
cc      CALL  PTNGSM( N,NPE,NOISEV,NOISEW,TAU2,SIG2,FF,B,Y,TREND,MJ,
cc     *              FIGMIN,FIGMAX,SS,K,LOC )
      CALL  POST3D( SS,LOC,K,NPE )
      RETURN
      E N D
cc      SUBROUTINE  DEFALT
      SUBROUTINE  DEFALT( Y,N,XMIN,XMAX,OUTMIN,OUTMAX )
C
C  ...  THIS SUBROUTINE SETS DEFAULT VALUES OF PARAMETERS  ...
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  Y(500)
cxx      DIMENSION  Y(N)
C
      INTEGER N
      DOUBLE PRECISION Y(N), XMIN, XMAX, OUTMIN, OUTMAX
c local
      DOUBLE PRECISION DY
C
cc      COMMON   /COMAAA/  Y
cc      COMMON   /C91214/  XMIN, XMAX, FIGMIN, FIGMAX, YM(10)
cc      COMMON   /C91215/  NOISEV, NOISEW, INITD, ITF, ITH
cc      COMMON   /C91216/  B, OUTMIN, OUTMAX
cc      COMMON   /C91219/  N, KK, NS, NFE, NPE
cc      NOISEV = 2
cc      NOISEW = 1
cc      INITD  = 1
cc      ITF    = 1
cc      ITH    = 1
cc      B = 1.00D0
      OUTMIN = -1.0D30
      OUTMAX =  1.0D30
cc      CALL  MAXMIN( Y,N,XMIN,XMAX,DY )
      CALL  MAXMINK( Y,N,XMIN,XMAX,DY )
cc      NS = 1
cc      NFE = N
cc      NPE = N
cc      FIGMIN = XMIN
cc      FIGMAX = XMAX
C
      RETURN
      E N D
cc      SUBROUTINE  NGSMTH( Y,P,F,S,T,N,K,DX,XMIN,Q,FF,PS,FS,SS,LOC )
      SUBROUTINE  NGSMTH( NOISEW,SIG2,BW,Y,F,N,K,DX,XMIN,Q,FF,SS,LOC,
     *                     OUTMIN,OUTMAX,NS,NFE,NPE )
C
C  ...  NON-GAUSSIAN SMOOTHER  ...
C
C     INPUTS:
C       Y(I):   TIME SERIES
C       P(I):   INITIAL DENSITY
C       N:      DATA LENGTH
C       K:      NUMBER OF INTERVALS IN STEP FUNCTION APPROXIMATION
C       DX:     WIDTH OF INTERVAL
C       XMIN:   MINIMUM OF THE INTERVAL
C       Q:      SYSTEM NOISE DENSITY
C     OUTPUTS:
C       FF:     LOG-LIKELIHOOD
C       LOC(I): LOCATION OF THE CENTER OF THE INTERVAL AT STEP I
C       SS:     SMOOTHED DENSITY
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  P(K), F(K), S(K), T(K), Y(N), Q(-K:K), LOC(N)
cc      REAL*4     FS(K,N), PS(K,N), SS(K,N)
cxx      DIMENSION  P(K), F(K), S(K), T(K), Y(N), Q(-K:K), LOC(NPE)
cxx      REAL*4     PS(K,NPE), SS(K,NPE)
cxx      DIMENSION     PS(K,NPE), SS(K,NPE)
C
      INTEGER NOISEW, N, K, NS, NFE, NPE, LOC(NPE)
      DOUBLE PRECISION SIG2, BW, Y(N), F(K), DX, XMIN, Q(-K:K), FF,
     1                 SS(K,NPE), OUTMIN, OUTMAX
c local
      INTEGER I, II, J
      DOUBLE PRECISION P(K), S(K), T(K), PS(K,NPE), PSUM, FINT, TSUM
C
cc      COMMON   /C91215/  NOISEV, NOISEW, INITD, ITF, ITH
cc      COMMON   /C91216/  B, OUTMIN, OUTMAX
cc      COMMON   /C91219/  NN, KK, NS, NFE, NPE
C
      FF = 0.0D0
C
      DO 200 II=NS,NPE
C
C  ...  CONVOLUTION (SYSTEM NOISE)  ...
C
      CALL  CONVOL( Q,F,K,P )
      CALL  NORMLZ( P,K,DX,PSUM )
C
C  ...  BAYES FORMULA  ...
C
      IF( Y(II).LE.OUTMIN .OR. Y(II).GE.OUTMAX .OR. II.GT.NFE ) THEN
      DO 110 I=1,K
cxx  110 F(I) = P(I)
      F(I) = P(I)
  110 CONTINUE
      ELSE
cc      CALL  BAYES( P,K,XMIN,DX,Y(II),F,LOC(II) )
      CALL  BAYES( NOISEW,SIG2,BW,P,K,XMIN,DX,Y(II),F,LOC(II) )
      CALL  NORMLZ( F,K,DX,FINT )
C
C  ...  LIKELIHOOD COMPUTATION  ...
C
      FF = FF + DLOG( FINT )
cc      IF( MOD(II,10).EQ.0)  WRITE(6,*) II,FF
      END IF
C
C  ...  SAVE FOR SMOOTHING  ...
C
      DO 130 I=1,K
      PS(I,II) = P(I)
cc  130 FS(I,II) = F(I)
cxx  130 SS(I,II) = F(I)
      SS(I,II) = F(I)
  130 CONTINUE
C
C  ...  SHIFT ORIGIN  ...
C
cc      CALL  SHIFT( F,K,T,II,N,LOC )
      CALL  SSHIFT( F,K,T,II,N,LOC )
C
  200 CONTINUE
C
C  ...  SMOOTHING  ...
C
cc      DO 190 J=NFE,NPE
cc      DO 190 I=1,K
cc  190 SS(I,J) = FS(I,J)
      DO 195 I=1,K
cc  195 S(I) = FS(I,NFE)
cxx  195 S(I) = SS(I,NFE)
      S(I) = SS(I,NFE)
  195 CONTINUE
C
      DO 300 II=NFE-1,NS,-1
cc      IF(MOD(II,10).EQ.0)  WRITE(6,*) II
      DO 210 I=1,K
      T(I) = 0.0D0
      P(I) = 0.0D0
cc  210 F(I) = FS(I,II)
cxx  210 F(I) = SS(I,II)
      F(I) = SS(I,II)
  210 CONTINUE
      DO 220 I=1,K
      J = I - (LOC(II+1)-LOC(II))
      IF( J.GE.1.AND.J.LE.K )  P(I) = PS(J,II+1)
cxx  220 IF( J.GE.1.AND.J.LE.K )  T(I) = S(J)
      IF( J.GE.1.AND.J.LE.K )  T(I) = S(J)
  220 CONTINUE
      DO 230 I=1,K
cxx  230 S(I) = T(I)
      S(I) = T(I)
  230 CONTINUE
C
cc      CALL  SCONVL( Q,S,P,F,K,T )
      CALL  SCONVLK( Q,S,P,F,K,T )
      CALL  NORMLZ( T,K,DX,TSUM )
C
      DO 240 I=1,K
      S(I) = T(I)
cxx  240 SS(I,II) = S(I)
      SS(I,II) = S(I)
  240 CONTINUE
  300 CONTINUE
C
      RETURN
      E N D
      SUBROUTINE  PINTVL( P,K,XMIN,DX,Y )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  P(K), Y(7), PROB(7), P1(401)
cxx      DIMENSION  P(K), Y(7), PROB(7), P1(K)
C
      INTEGER K
      DOUBLE PRECISION P(K), XMIN, DX, Y(7)
c local
      INTEGER I, II, J
      DOUBLE PRECISION PROB(7), P1(K), PP
C
      DATA  PROB /0.0013D0, 0.0227D0, 0.1587D0, 0.5000D0, 0.8413D0,
     *            0.9773D0,0.9987D0/
C
      P1(1) = 0.0
      DO 10 I=2,K
cxx   10 P1(I) = P1(I-1) + (P(I-1) + P(I))*DX/2.0
      P1(I) = P1(I-1) + (P(I-1) + P(I))*DX/2.0
   10 CONTINUE
C
      DO 31 J=1,7
cxx     DO 30 J=1,7
      PP = PROB(J)
cxx      DO 20 I=2,K
      DO 20 II=2,K
         I = II
      IF(P1(I-1).LE.PP .AND. P1(I).GT.PP)  GO TO 30
   20 CONTINUE
   30 Y(J) = XMIN + (I-2)*DX + DX*(PP - P1(I-1))/(P1(I) - P1(I-1))
   31 CONTINUE
C
      RETURN
      E N D
      SUBROUTINE  CONVOL( Q,S,K,P )
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cxx      DIMENSION  S(K), P(K), Q(-K:K)
C
      INTEGER K
      DOUBLE PRECISION Q(-K:K), S(K), P(K)
c local
      INTEGER I, J, J1, J2
      DOUBLE PRECISION SUM
C
      DO 20 I=1,K
      J1 = 1-I
      J2 = K-I
      SUM = 0.0
      DO 10 J=J1,J2
cxx   10 SUM = SUM + S(I+J)*Q(J)
      SUM = SUM + S(I+J)*Q(J)
   10 CONTINUE
cxx   20 P(I) = SUM
      P(I) = SUM
   20 CONTINUE
C
      RETURN
      E N D
cc      SUBROUTINE  SCONVL( Q,P,R,S,K,T )
      SUBROUTINE  SCONVLK( Q,P,R,S,K,T )
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cxx      DIMENSION  S(K), P(K), R(K), T(K), Q(-K:K)
C
      INTEGER K
      DOUBLE PRECISION Q(-K:K), P(K), R(K), S(K), T(K)
c local
      INTEGER I, J, J1, J2
      DOUBLE PRECISION SUM
C
      DO 20 I=1,K
      J1 = 1-I
      J2 = K-I
      SUM = 0.0D0
      DO 10 J=J1,J2
cxx   10 IF(P(I+J).GT.0.0D0)  SUM = SUM + P(I+J)/R(I+J)*Q(J)
      IF(P(I+J).GT.0.0D0)  SUM = SUM + P(I+J)/R(I+J)*Q(J)
   10 CONTINUE
cxx   20 T(I) = S(I)*SUM
      T(I) = S(I)*SUM
   20 CONTINUE
      RETURN
      E N D
      SUBROUTINE  TRANS1(K,DX,TAU2,BV,Q )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  Q(-K:K), PARAM(3)
C
      INTEGER K
      DOUBLE PRECISION DX, TAU2, BV, Q(-K:K)
c local
      INTEGER I, J
      DOUBLE PRECISION USERV, PARAM(3), X0, X, SUM
C
cc      EXTERNAL   FUNCT
C
      PARAM(1) = 0.0D0
      PARAM(2) = TAU2
      PARAM(3) = BV
C
      Q = 0
c
      DO 20 I=1-K,K-1
      X0 = -DX*I - DX/2
      SUM = (USERV(X0,PARAM) + USERV(X0+DX,PARAM))/2
      DO 10 J=1,49
      X = X0 + (DX*J)/50
cxx   10 SUM = SUM + USERV( X,PARAM )
      SUM = SUM + USERV( X,PARAM )
   10 CONTINUE
cxx   20 Q(I) = SUM*DX/50
      Q(I) = SUM*DX/50
   20 CONTINUE
C
      RETURN
      E N D
      SUBROUTINE  TRANS2( K,DX,TAU2,BV,Q )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  Q(-K:K), PARAM(3)
C
      INTEGER K
      DOUBLE PRECISION DX, TAU2, BV, Q(-K:K)
c local
      INTEGER I, J
      DOUBLE PRECISION GAUSS, PARAM(3), X0, X, SUM
C
cc      EXTERNAL   FUNCT
C
      PARAM(1) = 0.0D0
      PARAM(2) = TAU2
      PARAM(3) = BV
C
      Q = 0
c
      DO 20 I=1-K,K-1
      X0 = -DX*I - DX/2
      SUM = (GAUSS(X0,PARAM) + GAUSS(X0+DX,PARAM))/2
      DO 10 J=1,49
      X = X0 + (DX*J)/50
cxx   10 SUM = SUM + GAUSS( X,PARAM )
      SUM = SUM + GAUSS( X,PARAM )
   10 CONTINUE
cxx   20 Q(I) = SUM*DX/50
      Q(I) = SUM*DX/50
   20 CONTINUE
C
      RETURN
      E N D
      SUBROUTINE  TRANS3( K,DX,TAU2,BV,Q )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  Q(-K:K), PARAM(3)
C
      INTEGER K
      DOUBLE PRECISION DX, TAU2, BV, Q(-K:K)
c local
      INTEGER I, J
      DOUBLE PRECISION PEARSN, PARAM(3), X0, X, SUM
C
cc      EXTERNAL   FUNCT
C
      PARAM(1) = 0.0D0
      PARAM(2) = TAU2
      PARAM(3) = BV
C
      Q = 0
c
      DO 20 I=1-K,K-1
      X0 = -DX*I - DX/2
      SUM = (PEARSN(X0,PARAM) + PEARSN(X0+DX,PARAM))/2
      DO 10 J=1,49
      X = X0 + (DX*J)/50
cxx   10 SUM = SUM + PEARSN( X,PARAM )
      SUM = SUM + PEARSN( X,PARAM )
   10 CONTINUE
cxx   20 Q(I) = SUM*DX/50
      Q(I) = SUM*DX/50
   20 CONTINUE
C
      RETURN
      E N D
      SUBROUTINE  TRANS4( K,DX,TAU2,BV,Q )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  Q(-K:K), PARAM(3)
C
      INTEGER K
      DOUBLE PRECISION DX, TAU2, BV, Q(-K:K)
c local
      INTEGER I, J
      DOUBLE PRECISION TWOEXP, PARAM(3), X0, X, SUM
C
cc      EXTERNAL   FUNCT
C
      PARAM(1) = 0.0D0
      PARAM(2) = TAU2
      PARAM(3) = BV
C
      Q = 0
c
      DO 20 I=1-K,K-1
      X0 = -DX*I - DX/2
      SUM = (TWOEXP(X0,PARAM) + TWOEXP(X0+DX,PARAM))/2
      DO 10 J=1,49
      X = X0 + (DX*J)/50
cxx   10 SUM = SUM + TWOEXP( X,PARAM )
      SUM = SUM + TWOEXP( X,PARAM )
   10 CONTINUE
cxx   20 Q(I) = SUM*DX/50
      Q(I) = SUM*DX/50
   20 CONTINUE
C
      RETURN
      E N D

      SUBROUTINE  NORMLZ( P,K,DX,SUM )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  P(K)
      INTEGER K
      DOUBLE PRECISION P(K), DX, SUM
c local
      INTEGER I
C
      SUM = 0.0D0
      DO 10 I=1,K
cxx   10 SUM = SUM + P(I)
      SUM = SUM + P(I)
   10 CONTINUE
      SUM = SUM*DX
      DO 30 I=1,K
cxx   30 P(I) = P(I)/SUM
      P(I) = P(I)/SUM
   30 CONTINUE
C
      RETURN
      E N D
cc      SUBROUTINE  IDIST( P,K,P1,P2,XMIN,DX )
      SUBROUTINE  IDIST( P,K,P1,P2,XMIN,DX,INITD )
cxx      IMPLICIT REAL*8( A-H,O-Z )
cxx      DIMENSION  P(K), PARAM(3)
cc      COMMON   /C91215/  NOISEV, NOISEW, INITD, ITF, ITH
C
      INTEGER K, INITD 
      DOUBLE PRECISION P(K), P1, P2, XMIN, DX
c local
      INTEGER I
      DOUBLE PRECISION USERI, GAUSS, PARAM(3), X
C
      PARAM(1) = P1
      PARAM(2) = P2
C
      DO 10 I=1,K
      X = XMIN + DX*(I-1)
      IF( INITD.EQ.0 )  P(I) = USERI( X,PARAM )
      IF( INITD.EQ.1 )  P(I) = GAUSS( X,PARAM )
cc      IF( INITD.EQ.2 )  P(I) = UNIF ( X,PARAM )
cxx      IF( INITD.EQ.2 )  P(I) = UNIF ( X )
      IF( INITD.EQ.2 )  P(I) = 1.0D0
   10 CONTINUE
      RETURN
      END
cc      SUBROUTINE  BAYES( P,K,XMIN,DX,Y,F,LSHIFT )
      SUBROUTINE  BAYES( NOISEW,SIG2,BW,P,K,XMIN,DX,Y,F,LSHIFT )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  P(K), F(K), PARAM(3)
cc      COMMON   /C91215/  NOISEV, NOISEW, INITD, ITF, ITH
cc      COMMON   /C91218/  SIG2, TAU2, BW, BV
cc      EXTERNAL  USERW
cc      EXTERNAL  GAUSS
cc      EXTERNAL  PEARSN
cc      EXTERNAL  TWOEXP
cc      EXTERNAL  DBLEXP
C
      INTEGER NOISEW, K, LSHIFT
      DOUBLE PRECISION SIG2, BW, P(K), XMIN, DX, Y, F(K)
c local
      INTEGER I
      DOUBLE PRECISION USERW, GAUSS, PEARSN, TWOEXP, DBLEXP, PARAM(3)
C
      PARAM(2) = SIG2
      PARAM(3) = BW
C
      DO 10 I=1,K
      PARAM(1) = XMIN + DX*(I-1+LSHIFT)
      IF( NOISEW.EQ.0 )  F(I) = P(I)*USERW ( Y,PARAM )
      IF( NOISEW.EQ.1 )  F(I) = P(I)*GAUSS ( Y,PARAM )
      IF( NOISEW.EQ.2 )  F(I) = P(I)*PEARSN( Y,PARAM )
      IF( NOISEW.EQ.3 )  F(I) = P(I)*TWOEXP( Y,PARAM )
      IF( NOISEW.EQ.4 )  F(I) = P(I)*DBLEXP( Y,PARAM )
   10 CONTINUE
      RETURN
      E N D
      DOUBLE PRECISION FUNCTION  USERW( Y,PARAM )
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cxx      DIMENSION  PARAM(3)
C
      DOUBLE PRECISION Y, PARAM(3)
c local
      DOUBLE PRECISION C1, YMEAN, VAR
      DATA  C1  /2.506628275D0/
C
      YMEAN = 0.0D0
      VAR   = PARAM(1)
      USERW = DEXP( -(Y-YMEAN)**2/(2*VAR) )/(C1*DSQRT( VAR ))
      RETURN
      E N D
      DOUBLE PRECISION FUNCTION  USERV( X,PARAM )
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cxx      DIMENSION  PARAM(3)
C
      DOUBLE PRECISION X, PARAM(3)
c local
      DOUBLE PRECISION C1
      DATA  C1  /2.506628275D0/
C
      USERV = DEXP( -(X-PARAM(1))**2/(2*PARAM(2)) )
     *                /(C1*DSQRT( PARAM(2) ))
      RETURN
      E N D
      DOUBLE PRECISION FUNCTION  TWOEXP( X,PARAM )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  PARAM(2)
      DOUBLE PRECISION X, PARAM(2)
C
      TWOEXP = DEXP( -DABS(X-PARAM(1))*PARAM(2) )*PARAM(2)/2.0D0
      RETURN
C
      E N D
cc      SUBROUTINE  SHIFT( F,K,T,II,N,LOC )
      SUBROUTINE  SSHIFT( F,K,T,II,N,LOC )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cx      DIMENSION  F(K), T(K), LOC(*)
cxx      DIMENSION  F(K), T(K), LOC(N)
C
      INTEGER K, II, N, LOC(N)
      DOUBLE PRECISION F(K), T(K)
c local
      INTEGER I, IMAX, J
      DOUBLE PRECISION PMAX
C
C  ...  FIND THE POSTERIOR MODE AND SHIFT ORIGIN  ...
C
      PMAX = 0.0D0
c-----------
      IMAX = 1
c-----------
      DO 10  I=1,K
      IF( F(I).LE.PMAX )  GO TO 10
         PMAX = F(I)
         IMAX = I
   10 CONTINUE
      IF(II.LT.N)  LOC(II+1) = LOC(II) + IMAX - (K+1)/2
      DO 20 I=1,K
      J = I + IMAX - (K+1)/2
      T(I) = 0.0D0
cxx   20 IF(J.GE.1.AND.J.LE.K )  T(I) = F(J)
      IF(J.GE.1 .AND. J.LE.K )  T(I) = F(J)
   20 CONTINUE
      DO 30 I=1,K
cxx   30 F(I) = T(I)
      F(I) = T(I)
   30 CONTINUE
C
      RETURN
      E N D
      SUBROUTINE  POST3D( F,LOC,K,N )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      REAL*8     F(K,N)
cc      DIMENSION  YH(2000), FF(801), LOC(N)
cxx      REAL*8     FF(-K:2*K)
cxx      DIMENSION  LOC(N)
C
      INTEGER K, N, LOC(N)
      DOUBLE PRECISION F(K,N)
c local
      INTEGER I, II, I1, I2, J, NDIF, NN, N0
      DOUBLE PRECISION FF(-K:2*K)
C
cc      COMMON   /C91214/  X0, X1, F0, F1, YM(10)
cc      ANGLE = 70.0
C
cc      WX = 8.0
cc      WY = 2.0
cc      WZ =12.0
cc      CALL  DELX( X0,X1,DX )
cc      Z0 =  0.0
cc      Z1 =  N
cc      KK = 1
      NDIF = 1
cc      IF( N.GE.100 )  NDIF = 2
cc      IF( N.GE.200 )  NDIF = 4
cc      IF( N.GE.300 )  NDIF = 6
cc      IF( N.GE.500 )  NDIF = N/50
cc      DZ = NDIF*10
      NN = N/NDIF
      N0 = NDIF/2 + 1
cc      DO 10 I=1,K
cc   10 FF(I) = F(I,1)
cc      CALL  MAXMIN( FF,K,Y0,Y1,DY )
C
cc      CALL  PLOT3A( K,NN,WX,WY,WZ,ANGLE,X0,X1,DX,Y0,Y1,DY,Z0,Z1,DZ,
cc     *              KN,KK,YH )
c
      DO 100 J=N0,N,NDIF
cc      DO 80 I=1,K
      DO 80 I=-K,2*K
cxx   80 FF(I) = 0.0D0
      FF(I) = 0.0D0
   80 CONTINUE
      II = LOC(J)
      I1 = MAX0( 1,II )
      I2 = MIN0( K,K+II )
      DO 90 I=I1,I2
cxx   90 FF(I+II) = F(I,J)
      FF(I+II) = F(I,J)
   90 CONTINUE
cc      CALL  PLOT3B( FF,K,NN,WX,WY,WZ,ANGLE,Y0,Y1,KK,YH,ZS,ZC )
      DO 95 I=1,K
cxx   95 F(I,J) = FF(I)
      F(I,J) = FF(I)
   95 CONTINUE
  100 CONTINUE
cc      CALL  PLOT( -SNGL(ZC),-SNGL(ZS),-3 )
C
      RETURN
      E N D
      DOUBLE PRECISION FUNCTION  USERI( X,PARAM )
C
C  ...  User supplied density function  ...
C       (The following is an example of two-sided exponential dist.)
C
C     Inputs:
C        X:
C        PARAM(1):  mean
C        PARAM(2):  lambda
C     Output:
C        USERI:     density at X
C
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cxx      DIMENSION  PARAM(2)
C
      DOUBLE PRECISION X, PARAM(2)
c local
      DOUBLE PRECISION SIGMA
C
      SIGMA = DSQRT( PARAM(2) )
      USERI = SIGMA*DEXP( -SIGMA*DABS(X-PARAM(1)) )/2
      RETURN
      E N D
cc      SUBROUTINE  MAXMIN( X,N,XMIN0,XMAX0,DXL )
      SUBROUTINE  MAXMINK( X,N,XMIN0,XMAX0,DXL )
C
C  ...  This subroutine determines the minimum, the maximum and
C       the step width  ...
C
C     Inputs:
C        X(I):    data
C        N:       data length
C     Outputs:
C        XMIN0:   the minimum bound for the figure
C        XMAX0:   the maximum bound for the figure
C        DXL:     step width
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION X(N)
C
      INTEGER N
      DOUBLE PRECISION X(N), XMIN0, XMAX0, DXL
c local
      INTEGER I
      DOUBLE PRECISION XMIN, XMAX, DX, DIF
C
      XMIN = 1.0D30
      XMAX =-1.0D30
C
      DO 10 I=1,N
      IF( X(I) .LT. XMIN )   XMIN = X(I)
cxx   10 IF( X(I) .GT. XMAX )  XMAX = X(I)
      IF( X(I) .GT. XMAX )  XMAX = X(I)
   10 CONTINUE
      DX = XMAX-XMIN
      IF( DLOG10(DX) .GE. 0.0D0 )  DXL = INT( DLOG10(DX) )
      IF( DLOG10(DX) .LT. 0.0D0 )  DXL = INT( DLOG10(DX) )-1.0
      DXL = 10.0**DXL
      IF( DX/DXL.GT.6.0D0 )  DXL = DXL*2.0
      DIF = INT( DX/DXL )
      XMIN0 = INT( XMIN/DXL )*DXL
      XMAX0 = XMIN0 + DIF*DXL
      IF( XMIN0 .GT. XMIN )  XMIN0 = XMIN0 - DXL
   30 IF( XMAX0 .GE. XMAX )  GO TO 40
        XMAX0 = XMAX0 + DXL
        GO TO 30
   40 CONTINUE
C
      RETURN
      E N D
