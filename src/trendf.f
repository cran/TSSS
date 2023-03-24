C     PROGRAM 11.2  TREND
      SUBROUTINE TREND(Y,N,M,IOPT,TAU20,DELTA,TAUMAX,SIG2M,FFMAX,AIC,
ccc     &                 XSS,RS)
     &                 XSS,VSS,RS)
C
      INCLUDE 'TSSS.h'
C
C  ...  State space model for trend estimation  ...
C
C     Inputs:
C        M:       Order of the trend model
C        IOPT:    Search method
C        TAU20:   Initial estimate of TAU2
C        DELTA:   Search width
C     Parameters:
C        NMAX:    Adjustable dimension of Y
C        MJ:      Adjustable dimension of XF, VF, etc.
C        MAXM:    Adjustable dimension of A, B, C
C        NC:      Number of components
C     @TEST.FILTER2:  SEP.08,1990
C     MODIFIED  2/15/93
C
cc      PARAMETER( NMAX=1000,MJ=2,K=1,NDIM=NMAX,IDEV=1 )
      PARAMETER( K=1 )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  Y(NMAX)
cc      DIMENSION  F(MJ,MJ), G(MJ), H(MJ), Q(K,K)
cc      DIMENSION  XPS(MJ,NMAX), XFS(MJ,NMAX), XSS(MJ,NMAX)
cc      DIMENSION  VPS(MJ,MJ,NMAX), VFS(MJ,MJ,NMAX), VSS(MJ,MJ,NMAX)
cc      DIMENSION  XF(MJ), VF(MJ,MJ)
cxx      DIMENSION  Y(N)
cxx      DIMENSION  F(M,M), G(M), H(M), Q(K,K)
cxx      DIMENSION  XPS(M,N), XFS(M,N), XSS(M,N)
cxx      DIMENSION  VPS(M,M,N), VFS(M,M,N), VSS(M,M,N)
cxx      DIMENSION  XF(M), VF(M,M)
cxx      DIMENSION  RS(N)
C
      INTEGER N, M, IOPT
      DOUBLE PRECISION Y(N), TAU20, DELTA, TAUMAX, SIG2M, FFMAX, AIC,
     1                 XSS(M,N), RS(N)
c local
      DOUBLE PRECISION F(M,M), G(M), H(M), Q(K,K), XPS(M,N), XFS(M,N),
     1                 VPS(M,M,N), VFS(M,M,N), VSS(M,M,N), XF(M),
     2                 VF(M,M), OUTMIN, OUTMAX, SIG2, YMEAN, VAR, TAU2,
     3                 R, FF
C
      DATA  OUTMIN, OUTMAX /-1.0D30, 1.0D30/
C
      NDIM = N
      MJ = M
CCCCCC  TREND MODEL  CCCCC
      SIG2 =1.0D0
CCCCCCCCCCCCCCCCCCCCCCCCCC
cc      READ( 5,* )  M, IOPT
cc      READ( 5,* )  TAU20, DELTA
C
C  ...  Read Time Series  ...
C
cc      CALL  READTS( IDEV,Y,N )
C
      NS = 1
      NFE = N
      NPE = N
      CALL  MOMENT( Y,N/10,YMEAN,VAR )
C
      FFMAX  = -1.0D30
      TAU2 = 0.0D0
C
      DO 100  II=1,19
      IF( IOPT.NE.0 )  TAU2 = TAU20 + DELTA*(II-9)
      IF( IOPT.EQ.0 .AND. M.EQ.1 )  TAU2 = 2.0D0**(-II)
      IF( IOPT.EQ.0 .AND. M.GT.1 )  TAU2 = 2.0D0**(-II-5)
C
C  ...  Set Up State Space Model for Trend Estimation  ...
C
cc      CALL  SETTRN( M,MJ,F,G,H,R )
      CALL  SETTRN( M,F,G,H,R )
C
C  ...  Set initial values  ...
C
cc      CALL  ISTATE( M,MJ,YMEAN,VAR,XF,VF )
      CALL  ISTATE( M,YMEAN,VAR,XF,VF )
C
C  ...  Kalman Filter  ...
C
      Q(1,1) = TAU2
cc      CALL  FILTER( Y,XF,VF,F,G,H,Q,R,M,K,0,NS,NFE,NPE,MJ,NMAX,
      CALL  FILTER( Y,XF,VF,F,G,H,Q,R,M,K,0,NS,NFE,NPE,
     *              NDIM,OUTMIN,OUTMAX,VFS,VPS,XFS,XPS,FF,SIG2 )
      IF( FF.GT.FFMAX )  THEN
         FFMAX  = FF
         TAUMAX = TAU2
         SIG2M  = SIG2
      END IF
cc  100 WRITE(6,600)  TAU2, SIG2, FF
  100 CONTINUE
      AIC = -2*FFMAX + 2*(M+2)
cc      WRITE(6,610)  TAUMAX, SIG2M, FFMAX, AIC
C
C  ...  Filtering with the best parameter  ...
C
cc      CALL  ISTATE( M,MJ,YMEAN,VAR,XF,VF )
      CALL  ISTATE( M,YMEAN,VAR,XF,VF )
      Q(1,1) = TAUMAX
cc      CALL  FILTER( Y,XF,VF,F,G,H,Q,R,M,K,0,NS,NFE,NPE,MJ,NMAX,
      CALL  FILTER( Y,XF,VF,F,G,H,Q,R,M,K,0,NS,NFE,NPE,
     *              NDIM,OUTMIN,OUTMAX,VFS,VPS,XFS,XPS,FF,SIG2 )
C
C  ... Fixed Interval Smoother  ...
C
cc      CALL  SMOOTH( F,M,MJ,NDIM,NS,NFE,NPE,VFS,VPS,XFS,XPS,
cc     *              VSS,XSS )
      CALL  SMOOTH( F,M,NDIM,NS,NFE,NPE,VFS,VPS,XFS,XPS,VSS,XSS )
C
C  ...  Plot and print data, trend and noise  ...
C
C      CALL  PTTRND( Y,XSS,VSS,N,MJ,M,TAUMAX,SIG2,FF,AIC )
cc      CALL  PRTRND( Y,M,TAUMAX,SIG2,FF,AIC,XSS,MJ,N )
      CALL  PRTRND( Y,XSS,MJ,N,RS )
C
cc      STOP
      RETURN
cxx  600 FORMAT( 1H ,5X,F15.10,F12.6,F13.5 )
cxx  610 FORMAT( 1H ,'TAUMAX =',F15.10,3X,'SIG2 =',F15.10,3X,
cxx     *            'FF =',F13.5,3X,'AIC =',F13.4 )
      E N D
cc      SUBROUTINE  PRTRND( Y,M,TAU2,SIG2,FF,AIC,XSS,MJ,N )
      SUBROUTINE  PRTRND( Y,XSS,MJ,N,DATA )
C
C  ...  Print out Estimates  ...
C
C     Inputs:
C        Y:      Time Series
C        M:      Trend order
C        TAU2:   Variance of the system noise
C        SIG2:   Variance of the observational noise
C        FF:     Log-Likelihood of the mode
C        AIC:    AIC of the model
C        XSS:    Smoothed state
C        MJ:     Adjustable dimension of F and G
C        N:      Data LENGTH
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      CHARACTER  TITLE*72
cc      DIMENSION  XSS(MJ,*), Y(*), DATA(400)
cx      DIMENSION  XSS(MJ,*), Y(N), DATA(N)
cxx      DIMENSION  XSS(MJ,N), Y(N), DATA(N)
cc      COMMON  /CMDATA/  TITLE
C
      INTEGER MJ, N
      DOUBLE PRECISION Y(N), XSS(MJ,N), DATA(N)
C
cc      WRITE(6,600)
cc      WRITE(6,610)  TITLE
cc      WRITE(6,620)  M, TAU2, SIG2, FF, AIC
cc      WRITE(6,630)
cc      WRITE(6,640)  (XSS(1,I),I=1,N)
      DO 10 I=1,N
cxx   10 DATA(I) = Y(I) - XSS(1,I)
      DATA(I) = Y(I) - XSS(1,I)
   10 CONTINUE
cc      WRITE(6,650)
cc      WRITE(6,640)  (DATA(I),I=1,N)
      RETURN
C
cxx  600 FORMAT( 1H1,'PROGRAM 9.1:  PREDICTION AND INTERPOLATION' )
cxx  610 FORMAT( 1H ,'DATA: ',A72 )
cxx  620 FORMAT( 1H ,'M =',I3,3X,'TAU2 =',D12.5,3X,'SIG2 =',D12.5,
cxx     *        3X,'LOG-L =',F12.4,3X,'AIC =',F12.4 )
cxx  630 FORMAT( 1H ,'***  TREND COMPONENT  ***' )
cxx  640 FORMAT( 1H ,5D13.6 )
cxx  650 FORMAT( 1H ,'***  RESIDUALS  ***' )
      E N D
