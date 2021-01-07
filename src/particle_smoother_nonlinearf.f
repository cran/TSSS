*     particle_smoother_nonlinear.f         JUl.12,2019   For nonlinear smoothing
*     particle_smoother.f       for R function  JUl.12,2019
*     Monte-S.f   8/13/2002     Small memory version   
*     monte.f     12/06/2001    
C     @MONTE.JCGS.FIG32:                    JUN.  ,1995
C     @NFILTER.SIM4:    FORWARD SMOOTHING   NOV.20,1992
C     @NFILTER.SIM3:    SMOOTHING           NOV.19,1992
C     @NFILTER.SIM2:    MONTE CARLO METHOD  NOV.02,1992
C
      SUBROUTINE PFILTERNF(Y,N,M,LAG,SIG2,TAU2,XMIN,XMAX,IX,T,FF)
cc      PARAMETER( NMAX=500,IDEV=1,MM=1000000, LAG=20 )
cc      IMPLICIT REAL*8(A-H,O-Z)
cc      real*4  PS, PST
cc      DIMENSION  Y(NMAX), T(NMAX,8), PS(MM,0:LAG), PST(MM,0:LAG)
cc      CHARACTER   DATE*8, TIME2*10
      INTEGER :: N, M, LAG, IX
      REAL(8) :: Y(N), SIG2, TAU2, XMIN, XMAX, T(N,8), FF
      INTEGER :: IST, L
      REAL(8) :: PS(M,0:LAG), PST(M,0:LAG)
cc      COMMON  /INIT/ IX
      common  /COMIST/  IST
cc      IX = 1992110113
      IST = 0
C
cc      MODEL = 1
*           M is the number of particles. 1000 - 100000
*           LAG:  Fixed lag smoother LAG=0 => filter
*           IX:  initial seed for random number generator
cc      M = 1000000
cc      M = 100000
      L = 1
cc      INID = 0
*      
cc      IX = 2019071117
C      IX = 2002081119
*      IX = 1997102221
cc      IX0= IX
*
*     sig2  observation noise variance
*     tau2: system noise variance
cc      SIG2 = 10.0D0
cc      TAU2 =  1.0D0
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
cc     *              INID,IX )
      CALL  FILTERNL( Y,N,M,L,T,PS,PST,SIG2,TAU2,LAG,XMIN,XMAX,IX,FF)
cc      CALL  DATE_AND_TIME( DATE,TIME2 )     
cc      WRITE(6,650)  DATE, TIME2      
C
cc      open( 6,file='results.dat' )
cc      CALL  PTDATA( Y,T,N,M,L,MODEL,SIG2,TAU2,IX0,LAG,NMAX,FF )
cc      close( 6 )
C
cc      STOP
      return
      E N D
cc      SUBROUTINE  FILTER( Y,N,M,L,T,PS,PST,MODEL,SIG2,TAU2,LAG,NMAX,FF,
cc     *                    MM, INID )
      SUBROUTINE  FILTERNL(Y,N,M,L,T,PS,PST,SIG2,TAU2,LAG,XMIN,XMAX,IX,
     *                     FF)
cc      IMPLICIT REAL*8(A-H,O-Z)
cc      real*4  PS, PST
cc      DIMENSION  Y(N), F(1000000), P(1000000), S(1000000)
cc      DIMENSION  PS(MM,0:LAG), T(NMAX,8), PST(MM,0:LAG)
cc      DIMENSION  TT(7)
      INTEGER :: N, M, L, LAG, IX
      REAL(8) :: Y(N), T(N,8), PS(M,0:LAG), PST(M,0:LAG), SIG2, TAU2,
     *            XMIN, XMAX, FF
      REAL(8) :: TAU, F(M), P(M), S(M), TT(7), RNOR, RR, SUM, GAUS2,
ccx     *           RUNIFT
     *           random 
cc      COMMON  /INIT/ IX
c----- Initialize Mersenne Twister with given seed value
      call init(IX)
c-----
*
      TAU = DSQRT( TAU2 )
      FF = 0.0D0
*       Specify initial distribution : you can change these
      DO 10 I=1,M
cc   10 F(I) = RNOR( 1.0D0 )*DSQRT(5.0D0)
      F(I) = RNOR( 1.0D0 )*DSQRT(5.0D0)
   10 CONTINUE
*
*    Monte Carlo filtering   N is data length
*
      DO 100 II=1,N
*
*     Prediction
*
      I0 = min0( ii,lag )

      DO 20 I=1,M
      P(I) = F(I)/2 + 25*F(I)/(1+F(I)**2) + 8*DCOS(1.2D0*II)
cc     *       + RNOR( 1.0D0 )*TAU
     *       + RNOR( 1.0D0 )*TAU
   20 CONTINUE
C
*     Filtering
*
*   GAUS2:  Gaussian density, change this if you want other distribution
      SUM = 0.0D0
      DO 30 I=1,M*L
      S(I) = GAUS2( Y(II)-P(I)**2/20.0D0,SIG2 )
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

