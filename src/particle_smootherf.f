*     particle_smoother.f       for R function  JUl.12,2019
*     Monte-S.f   8/13/2002     Small memory version   
*     monte.f     12/06/2001    
C     @MONTE.JCGS.FIG32:                    JUN.  ,1995
C     @NFILTER.SIM4:    FORWARD SMOOTHING   NOV.20,1992
C     @NFILTER.SIM3:    SMOOTHING           NOV.19,1992
C     @NFILTER.SIM2:    MONTE CARLO METHOD  NOV.02,1992
C
      SUBROUTINE PFILTERF(Y,N,M,MODEL,LAG,INID,SIG2,TAU2,ALPHA,BIGTAU2,
     *                    SIG2I,XMIN,XMAX,IX,T,FF)
cc      PARAMETER( NMAX=500,IDEV=1,MM=1000000, LAG=20 )
cc      IMPLICIT REAL*8(A-H,O-Z)
cc	  real*4  PS, PST
cc      DIMENSION  Y(NMAX), T(NMAX,8), PS(MM,0:LAG), PST(MM,0:LAG)
cc      CHARACTER   DATE*8, TIME2*10
      INTEGER :: N, M, MODEL, LAG, INID, IX
      REAL(8) :: Y(N), SIG2, TAU2, ALPHA, BIGTAU2, SIG2I, XMIN, XMAX,
     *           T(N,8), FF
      INTEGER :: L
      REAL(8) :: PS(M,0:LAG), PST(M,0:LAG)
cc      COMMON  /INIT/ IX
cc      common  /COMIST/  IST

cc      IX = 1992110113
cc      IST = 0
C
cc      MODEL = 1
*           M is the number of particles. 1000 - 100000
*           LAG:  Fixed lag smoother LAG=0 => filter
*           IX:  initial seed for random number generator
cc      M = 1000000
cc      M = 100000
      L = 1
cc	  INID = 0
*
cc      IX = 2019071117
C     IX = 2002081119
*      IX = 1997102221
cc      IX0= IX
*
*     sig2  observation noise variance
*     tau2: system noise variance
cc      SIG2 = 1.045D0
cc      TAU2 = 3.5300D-5
cc      IF( MODEL.EQ.0 )  THEN
cc        SIG2 = 1.043D0
cc        TAU2 = 1.2200D-2
cc        SIG2 = 1.048D0
cc        TAU2 = 1.4000D-2
cc      END IF
cc      IF( MODEL.EQ.2 )  THEN
cc        SIG2 = 1.03D0
cc        TAU2 = 0.00013d0
cc      END IF
*
*      READTS reads in data: change this subroutine
*
cc      CALL  READTS( IDEV,Y,N )
cc      WRITE(6,600)  MODEL, M, L, IX, INID
cc  600 FORMAT( 1H ,'MODEL =',I2,3X,'M =',I6,3X,'L =',I3,3X,'IX =',I12,
cc     *            3X,'INID =',I3 )
C
cc      CALL  DATE_AND_TIME( DATE,TIME2 )     
cc      WRITE(6,650)  DATE, TIME2      
cc  650 FORMAT( 'Date',2X,A8,3X,'Time:',2X,A10 )
cc      CALL  FILTER( Y,N,M,L,T,PS,PST,MODEL,SIG2,TAU2,LAG,NMAX,FF,MM,
cc     *              INID )
      CALL  CFILTER( Y,N,M,L,T,PS,PST,MODEL,SIG2,TAU2,LAG,ALPHA,BIGTAU2,
     *               SIG2I,XMIN,XMAX,FF,INID,IX )
cc      CALL  DATE_AND_TIME( DATE,TIME2 )     
cc      WRITE(6,650)  DATE, TIME2      
C
cc      open( 6,file='results.dat' )
cc      CALL  PTDATA( Y,T,N,M,L,MODEL,SIG2,TAU2,IX0,LAG,NMAX,FF )
cc	close( 6 )
C
cc      STOP
      return
      E N D
cc      SUBROUTINE  FILTER( Y,N,M,L,T,PS,PST,MODEL,SIG2,TAU2,LAG,NMAX,FF,
cc     *                    MM, INID )
      SUBROUTINE  CFILTER( Y,N,M,L,T,PS,PST,MODEL,SIG2,TAU2,LAG,ALPHA,
     *                     BIGTAU2,SIG2I,XMIN,XMAX,FF,INID,IX)
cc      IMPLICIT REAL*8(A-H,O-Z)
cc	  real*4  PS, PST
cc      DIMENSION  Y(N), F(1000000), P(1000000), S(1000000)
cc      DIMENSION  PS(MM,0:LAG), T(NMAX,8), PST(MM,0:LAG)
cc      DIMENSION  TT(7)
      INTEGER :: N, M, L, MODEL, LAG, INID, IX
      REAL(8) :: Y(N), T(N,8), PS(M,0:LAG), PST(M,0:LAG), SIG2, TAU2,
     *           ALPHA, BIGTAU2, SIG2I, XMIN, XMAX, FF
      REAL(8) :: PI, TAU, F(M), P(M), S(M), TT(7), RNOR, RR, SUM, GAUS2, 
ccx     1           RUNIFT
     1           random, BIGTAU, GB, ymin, ymax
cc      COMMON  /INIT/ IX
      DATA  PI /3.1415926535D0/
c
c----- Initialize Mersenne Twister with given seed value
      call init(IX)
c-----
*
      TAU = DSQRT( TAU2 )
      BIGTAU = DSQRT( BIGTAU2 )
      GB = DSQRT( SIG2I )
      ymin = MINVAL(Y)
      ymax = MAXVAL(Y)
      FF = 0.0D0

*       Specify initial distribution : you can change these
      DO 10 I=1,M
      IF( INID.EQ.3 )  F(I) = 0.0D0
cc      IF( INID.EQ.2 )  F(I) = RCAU( 1.0D0 )
      IF( INID.EQ.2 )  F(I) = GB * DTAN( random()*PI )
cc      IF( INID.EQ.1 )  F(I) = -4.0D0 + 8*RUNIFT(IX)
      IF( INID.EQ.1 )  F(I) = ymin + (ymax - ymin)*random()
cc   10 IF( INID.EQ.0 )  F(I) = RNOR( 1.0D0 )
      IF( INID.EQ.0 )  F(I) = RNOR( GB )
   10 CONTINUE
*
*    Monte Carlo filtering   N is data length
*
      DO 100 II=1,N
*
*     Prediction
*
      I0 = min0( ii,lag )
cc
      DO 20 I=1,M
      IF( MODEL.EQ.0 ) P(I) = F(I) + RNOR( 1.0D0 )*TAU
cc      IF( MODEL.EQ.1 ) P(I) = F(I) + RCAU( 1.0D0 )*TAU
ccx      IF( MODEL.EQ.1 ) P(I) = F(I) + DTAN(RUNIFT(IX)*PI )*TAU
      IF( MODEL.EQ.1 ) P(I) = F(I) + DTAN( random()*PI )*TAU
      IF( MODEL.EQ.2 ) then
ccx      RR = RUNIFT(IX)
      RR = random()
cc      IF( RR.le.0.991d0 )  P(I) = F(I) + RNOR( 1.0D0 )*TAU
cc      IF( RR.GT.0.991d0 )  P(I) = F(I) + RNOR( 1.0d0 )*2.0d0
      IF( RR.le.ALPHA )  P(I) = F(I) + RNOR( 1.0D0 )*TAU
      IF( RR.GT.ALPHA )  P(I) = F(I) + RNOR( 1.0d0 )*BIGTAU
      end if
   20 CONTINUE
C
*     Filtering
*
*   GAUS2:  Gaussian density, change this if you want other distribution
      SUM = 0.0D0
      DO 30 I=1,M*L
      S(I) = GAUS2( Y(II)-P(I),SIG2 )
cc   30 SUM = SUM + S(I)
      SUM = SUM + S(I)
   30 CONTINUE
*
*     FF: log-likelihood
*
      FF = FF + DLOG( SUM/M )
cc   35 CONTINUE
*
*     Re-sampling
*
      S(1) = S(1)/SUM
      DO 40 I=2,M
cc   40 S(I) = S(I-1) + S(I)/SUM
      S(I) = S(I-1) + S(I)/SUM
   40 CONTINUE
C  
      JJ = 1
      DO 50 I=1,M
ccx      RR = (I-1+RUNIFT(IX))/M
      RR = (I-1+random())/M
      DO 60 J=JJ,M*L
      IF( S(J).GE.RR )  GO TO 55
   60 CONTINUE
cc   55 IF( J.GE.M )  J = M
   55 CONTINUE
      IF( J.GE.M )  J = M
      F(I) = P(J)
      JJ = J
      I0 = min0( ii,lag )
      DO 56 IJ=1,I0
cc   56 PST(I,IJ) = PS(J,IJ-1)
      PST(I,IJ) = PS(J,IJ-1)
   56 CONTINUE
   50 CONTINUE
*
cc      DO 65 IJ=1,I0
      DO 66 IJ=1,I0
      DO 65 I=1,M
cc   65 PS(I,IJ) = PST(I,IJ)
      PS(I,IJ) = PST(I,IJ)
   65 CONTINUE
   66 CONTINUE
C
      DO 70 I=1,M
cc   70 PS(I,0) = F(I)
      PS(I,0) = F(I)
   70 CONTINUE
cc   90 CONTINUE

      IF( II.LE.LAG )  GO TO 100
      do 320 I=1,8
cc  320 T(II-LAG,I) = 0.0d0
      T(II-LAG,I) = 0.0d0
  320 CONTINUE

cc      CALL  DENSTY( PS(1,LAG),M,TT )
      CALL  DENSTY1( PS(1,LAG),M,TT,XMIN,XMAX )

      DO 310 I=1,7
cc  310 T(II-LAG,I) = TT(I)
      T(II-LAG,I) = TT(I)
  310 CONTINUE

  100 CONTINUE
*
*     Compute Posterior Distributions
*     T(II,J), J=4: MEdian, other: 1,2,3 sigma interval
*              J=8  MEAN
*
cc      WRITE(6,660)  TAU2, SIG2, FF
cc  660 FORMAT( 1H ,F15.8,F10.4,F13.4 )
      do 200 II=0,LAG-1
cc      CALL  DENSTY( PS(1,II),M,TT )
      CALL  DENSTY1( PS(1,II),M,TT,XMIN,XMAX )
      DO 210 I=1,7
cc  210 T(N-II,I) = TT(I)
      T(N-II,I) = TT(I)
  210 CONTINUE
*      DO 110 I=1,M
*  110 SUM = SUM + PS(I,II)
*      T(II,8) = SUM/M
  200 CONTINUE
C
      RETURN
      E N D

      double precision FUNCTION  RNOR( GB )
cc      IMPLICIT REAL*8(A-H,O-Z)
      REAL(8) :: GB
      INTEGER :: I
ccx      REAL(8) :: X, Y, S, RUNIFT
      REAL(8) :: X, Y, S, random
cc      COMMON  /INIT/ IX
      DATA   I /0/
      data   y/0.0d0/
      data   s/0.0d0/
ccx      IF( IX.EQ.0 )  IX = 1991041911
      IF( I.GT.0 )  GO TO 30
ccx   10 X = 2.0*RUNIFT(IX) - 1.0
ccx      Y = 2.0*RUNIFT(IX) - 1.0
   10 X = 2.0*random() - 1.0
      Y = 2.0*random() - 1.0
      S = X**2 + Y**2
      IF( S.GE.1.0 )  GO TO 10
      S = DSQRT( -2.0*DLOG(S)/S )
      RNOR = X*S
      I = 1
      GO TO 40
   30 RNOR = Y*S
      I = 0
   40 CONTINUE
      RNOR = RNOR*GB
C      WRITE(6,*)  RNOR
      RETURN
      E N D
      DOUBLE PRECISION FUNCTION  GAUS2( X,SIG2 )
C
cc      IMPLICIT  REAL*8(A-H,O-Z)
      REAL(8) :: X, SIG2, C1
      DATA  C1  /2.506628275D0/
C
      GAUS2 = DEXP( -X**2/(2*SIG2) )/(C1*DSQRT(SIG2))
      RETURN
      E N D
      
      DOUBLE PRECISION FUNCTION  CAUC2( X,SIG2 )
C
cc      IMPLICIT  REAL*8(A-H,O-Z)
      REAL(8) :: X, SIG2, PI
      DATA  PI  /3.141592653D0/
C
      CAUC2 = DSQRT(SIG2)/PI/(X**2 + SIG2)
      RETURN
      E N D
      
      SUBROUTINE  SORT( Y,N )
cc      IMPLICIT REAL*8 (A-H,O-Z)
cc      DIMENSION  Y(N)
      INTEGER :: N
      REAL(8) :: Y(N), YY
C
      IF( Y(2).LT.Y(1) )  THEN
         YY   = Y(1)
         Y(1) = Y(2)
         Y(2) = YY
      END IF
      DO 100  I=3,N
      JJ = I
      DO 20  J=I-1,1,-1
         IF( Y(I).LT.Y(J) )  JJ = J
   20 CONTINUE
      IF( JJ.LT.I )  THEN
            YY = Y(I)
            DO 10 L=I-1,JJ,-1
cc   10          Y(L+1) = Y(L)
               Y(L+1) = Y(L)
   10       CONTINUE
            Y(JJ) = YY
      END IF
  100 CONTINUE
      RETURN
      E N D
      
cc      SUBROUTINE  DENSTY( P,M,T )
cc      IMPLICIT REAL*8(A-H,O-Z)
cc	  real*4  P
cc      REAL*4  PROB
cc      DIMENSION  P(M), Q(1000), QQ(0:1000), PROB(7), T(7)
      SUBROUTINE  DENSTY1( P,M,T,XMIN,XMAX )
      INTEGER :: M
      REAL(8) :: P(M), T(7),XMIN, XMAX
      PARAMETER(K = 1000)
      REAL(8) :: Q(K), QQ(0:K), PROB(7), PP, DX
cc      DATA  PROB /0.0013, 0.0227, 0.1587, 0.5000, 0.8413, 0.9773,
cc     *            0.9987/
      DATA  PROB /0.0013D0, 0.0227D0, 0.1587D0, 0.5D0, 0.8413D0,
     *             0.9773D0, 0.9987D0/
C
cc      XMIN = -4.0D0
cc      XMAX =  4.0D0
cc      K = 1000
      DX = (XMAX-XMIN)/K
C
cc      DO 10 I=1,K
cc   10 Q(I) = 0.0D0
      Q(1:k) = 0.0D0
      DO 20 I=1,M
cc      J = (P(I) - XMIN)/DX
      J = INT((P(I) - XMIN)/DX)
      IF( J.Le.0 )  J = 1
      IF( J.GT.K )  J = K
cc   20 Q(J) = Q(J) + 1
      Q(J) = Q(J) + 1
   20 CONTINUE
      DO 30 I=1,K
cc   30 Q(I) = Q(I)/M
      Q(I) = Q(I)/M
   30 CONTINUE
C
      QQ(0) = 0.0D0
      DO 40 I=1,K
cc   40 QQ(I) = QQ(I-1) + Q(I)
      QQ(I) = QQ(I-1) + Q(I)
   40 CONTINUE
C     WRITE(6,*)  (P(I),I=1,K)
C     WRITE(6,*)  (QQ(I),I=1,K)
C
cc      DO 60 J=1,7
      DO 61 J=1,7
      PP = PROB(J)
      DO 50 I=1,K
      IF(QQ(I-1).LE.PP .AND. QQ(I).GT.PP)  GO TO 60
   50 CONTINUE
   60 T(J) = XMIN + (I-1)*DX + DX*(PP - QQ(I-1))/(QQ(I) - QQ(I-1))
   61 CONTINUE
C
      RETURN
      E N D
