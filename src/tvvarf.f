C     PROGRAM 13.1  TVVAR
      SUBROUTINE TVVARF(Y,N0,M,TAU20,IOPT,DELTA,TVVAR,NORDAT,Y2,N,TREND,
     *                  NOISE,TAUMAX,SIG2M,FFMAX,AIC)
C
      INCLUDE 'TSSS_f.h'
C
C  ...  TIME-VARYING VARIANCE  ...
C
C     Inputs:
C        M:       Trend Order
C        IOPT:    Search Method
C        TAU20:   Initial Estimate of TAU2
C        DELTA:   Search Width
C     Parameters:
C        NMAX:    Adjustable dimension of Y
C        MJ:      Adjustable dimension of XF, VF, etc.
C        MAXM:    Adjustable dimension of A, B, C
C        NC:      Number of components
C     @TEST.FILTER2:  SEP.08,1990
C     MODIFIED  2/15/93
C
cc      !DEC$ ATTRIBUTES DLLEXPORT::TVVARF
C
cc      PARAMETER( MJ1=3000,NMAX=MJ1/2,MJ=2,K=1,NDIM=NMAX )
      PARAMETER( K=1 )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  Y(MJ1), Y2(NMAX)
cc      DIMENSION  F(MJ,MJ), G(MJ), H(MJ), Q(K,K)
cc      DIMENSION  XPS(MJ,NMAX), XFS(MJ,NMAX), XSS(MJ,NMAX)
cc      DIMENSION  VPS(MJ,MJ,NMAX), VFS(MJ,MJ,NMAX), VSS(MJ,MJ,NMAX)
cc      DIMENSION  XF(MJ), VF(MJ,MJ)
cxx      DIMENSION  Y(N0), Y2(N0/2)
cxx      DIMENSION  F(M,M), G(M,K), H(M), Q(K,K)
cxx      DIMENSION  XPS(M,N0/2), XFS(M,N0/2), XSS(M,N0/2)
cxx      DIMENSION  VPS(M,M,N0/2), VFS(M,M,N0/2), VSS(M,M,N0/2)
cxx      DIMENSION  XF(M), VF(M,M)
cxx      DIMENSION  TVVAR(N0/2), TREND(N0/2,3)
cxx      REAL*8     NORDAT(N0), NOISE(N0/2)
C
      INTEGER :: N0, M, IOPT, N
      REAL(8) :: Y(N0), TAU20, DELTA, TVVAR(N0/2), NORDAT(N0), Y2(N0/2),
     1           TREND(N0/2,3), NOISE(N0/2), TAUMAX, SIG2M, FFMAX, AIC
      REAL(8) :: F(M,M), G(M,K), H(M), Q(K,K), XPS(M,N0/2), XFS(M,N0/2),
     1           XSS(M,N0/2), VPS(M,M,N0/2), VFS(M,M,N0/2),
     2           VSS(M,M,N0/2), XF(M), VF(M,M), OUTMIN, OUTMAX, SIG2,
     3           YMIN, YMEAN, VAR, TAU2, R, FF
C
      DATA  OUTMIN, OUTMAX /-1.0D30, 1.0D30/
C
cc      IDEV = 1
      SIG2 =1.0D0
cc      READ( 5,* )  M, IOPT
cc      IF( IOPT.EQ.1 )  READ( 5,* )  TAU20, DELTA
C
C  ...  Read Time Series  ...
C
cc      CALL READTS ( IDEV,Y,N0 )
      N = N0/2
      YMIN = 1.0D30
      DO 10 I=1,N
      J = 2*I
      Y2(I) = (Y(J-1)**2 + Y(J)**2)/2
cxx   10 IF( Y2(I).GT.0.0D0 .AND. Y2(I).LT.YMIN )  YMIN = Y2(I)
      IF( Y2(I).GT.0.0D0 .AND. Y2(I).LT.YMIN )  YMIN = Y2(I)
   10 CONTINUE
      DO 20 I=1,N
cxx   20 Y2(I) = DLOG( DMAX1(Y2(I),YMIN/2) )
      Y2(I) = DLOG( DMAX1(Y2(I),YMIN/2) )
   20 CONTINUE
C
      NS = 1
      NFE = N
      NPE = N
      NDIM = N
      CALL  MOMENT( Y2,N/10,YMEAN,VAR )      
C
      FFMAX  = -1.0D30
C
      DO 100  II=1,19
cc      IF( IOPT.NE.0 )  TAU2 = TAU20 + DELTA*(II-9)
      TAU2 = TAU20 + DELTA*(II-9)
      IF( IOPT.EQ.0 .AND. M.EQ.1 )  TAU2 = 2.0D0**(-II)
      IF( IOPT.EQ.0 .AND. M.GT.1 )  TAU2 = 2.0D0**(-II-5)
C
C  ...  Set Up State Space Model for Trend Estimation  ...
C
cc      CALL  SETTRN( M,MJ,F,G,H,R )
      CALL  SETTRN( M,F,G,H,R )
      R = (3.14159265D0)**2/6
C
C  ...  Set initial values  ...
C
cc      CALL  ISTATE( M,MJ,YMEAN,VAR,XF,VF )
      CALL  ISTATE( M,YMEAN,VAR,XF,VF )
C
C  ...  Kalman Filter  ...
C
      Q(1,1) = TAU2
cc      CALL  FILTER( Y2,XF,VF,F,G,H,Q,R,M,K,1,NS,NFE,NPE,MJ,NMAX,
      CALL  FILTER( Y2,XF,VF,F,G,H,Q,R,M,K,1,NS,NFE,NPE,
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
cc      CALL  FILTER( Y2,XF,VF,F,G,H,Q,R,M,K,1,NS,NFE,NPE,MJ,NMAX,
      CALL  FILTER( Y2,XF,VF,F,G,H,Q,R,M,K,1,NS,NFE,NPE,
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
cc      CALL  PTTRND( Y2,XSS,VSS,N,MJ,M,TAUMAX,SIG2,FF,AIC )
cc      CALL  PRVAR( Y,M,TAUMAX,SIG2,FF,AIC,XSS,MJ,N,N0 )
cxx      CALL  PTTRND( Y2,XSS,VSS,N,MJ,M,TAUMAX,SIG2,FF,AIC,TREND,NOISE )
cxx      CALL  PRVAR( Y,M,TAUMAX,SIG2,FF,AIC,XSS,MJ,N,N0,TVVAR,NORDAT )
      CALL  PTTRND( Y2,XSS,VSS,N,M,SIG2,TREND,NOISE )
      CALL  PRVAR( Y,M,XSS,N,N0,TVVAR,NORDAT )
C
      RETURN
cxx  600 FORMAT( 1H ,5X,F15.10,F12.6,F13.5 )
cxx  610 FORMAT( 1H ,'TAUMAX =',F15.10,3X,'SIG2 =',F15.10,3X,
cxx     *            'FF =',F13.5,3X,'AIC =',F13.4 )
      E N D
cc      SUBROUTINE  PRVAR( Y,M,TAU2,SIG2,FF,AIC,XSS,MJ,N,N0 )
cxx      SUBROUTINE  PRVAR( Y,M,TAU2,SIG2,FF,AIC,XSS,MJ,N,N0,DATA,DATA1 )
      SUBROUTINE  PRVAR( Y,M,XSS,N,N0,DATA,DATA1 )
C
C  ...  Print out changing variance and normalized data  ...
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
cc      CHARACTER  TITLE*72
cc      DIMENSION  XSS(MJ,*), Y(*), DATA(1500)
cx      DIMENSION  XSS(M,*), Y(*), DATA(N), DATA1(N0)
cxx      DIMENSION  XSS(M,N), Y(N0), DATA(N), DATA1(N0)
cc      COMMON  /CMDATA/  TITLE
C
      INTEGER :: M, N, N0
      REAL(8) :: Y(N0), XSS(M,N), DATA(N), DATA1(N0)
C
cc      WRITE(6,600)
cc      WRITE(6,610)  TITLE
cc      WRITE(6,620)  M, TAU2, SIG2, FF, AIC
      DO 10 I=1,N
cxx   10 DATA(I) = DEXP(XSS(1,I)+0.57721D0)
      DATA(I) = DEXP(XSS(1,I)+0.57721D0)
   10 CONTINUE
cc      WRITE(6,630)
cc      WRITE(6,640)  (DATA(I),I=1,N)
      DO 20 I=1,N0
      J = (I+1)/2
cc   20 Y(I) = Y(I)/DSQRT( DATA(J) )
cxx   20 DATA1(I) = Y(I)/DSQRT( DATA(J) )
      DATA1(I) = Y(I)/DSQRT( DATA(J) )
   20 CONTINUE
cc      WRITE(6,650)
cc      WRITE(6,640)  (Y(I),I=1,N0)
      RETURN
C
cxx  600 FORMAT( 1H1,'PROGRAM 12.3:  CHANGING VARIANCE ESTIMATION' )
cxx  610 FORMAT( 1H ,'DATA: ',A72 )
cxx  620 FORMAT( 1H ,'M =',I3,3X,'TAU2 =',D12.5,3X,'SIG2 =',D12.5,
cxx     *        3X,'LOG-L =',F12.4,3X,'AIC =',F12.4 )
cxx  630 FORMAT( 1H ,'***  TIME VARYING VARIANCE  ***' )
cc  640 FORMAT( 1H ,5D13.6 )
cxx  640 FORMAT( 5(f13.6,', '))
cxx  650 FORMAT( 1H ,'***  NORMALIZED DATA  ***' )
      E N D
cc      SUBROUTINE  PTTRND( Y,XSS,VSS,N,MJ,M,TAU2,SIG2,FF,AIC )
cxx      SUBROUTINE  PTTRND( Y,XSS,VSS,N,MJ,M,TAU2,SIG2,FF,AIC,TREND,DATA)
      SUBROUTINE  PTTRND( Y,XSS,VSS,N,M,SIG2,TREND,DATA)
C
C  ...  Plot Original time series, trend ans residuals  ...
C
C     Inputs:
C        Y(I):   Time Series
C        XSS:    Smoothed state
C        VSS:    Covariance matrices
C        N:      Data length
C        MJ:     Adjustable dimension of F and G
C        M:      Model order
C        TAU2:   Variance of the model
C        FF:     Log-Likelihood of the mode
C        AIC:    AIC of the model
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      CHARACTER  VNAME*8
cxx      DIMENSION  Y(N)
cc      DIMENSION  XSS(MJ,N), VSS(MJ,MJ,N)
cc      DIMENSION  DATA(1500), VNAME(20), VALUE(20)
cxx      DIMENSION  XSS(M,N), VSS(M,M,N)
cxx      DIMENSION  DATA(N), TREND(N,-1:1)
C
      INTEGER :: N, M
      REAL(8) :: Y(N), XSS(M,N), VSS(M,M,N), SIG2, TREND(N,-1:1),
     1           DATA(N)
C
cc      VNAME(1) = 'M1    = '
cc      VNAME(2) = 'TAU2  = '
cc      VNAME(3) = 'SIG2  = '
cc      VNAME(4) = 'FF    = '
cc      VNAME(5) = 'AIC   = '
cc      VALUE(1) = M
cc      VALUE(2) = TAU2
cc      VALUE(3) = SIG2
cc      VALUE(4) = FF
cc      VALUE(5) = AIC
cc      CALL  PLOTS
C     call  plots( 1,0,0,1,0 )
C     call  form( 1 )
C     call  factor( 10.0 )
cc      CALL  HEADER( 'PROGRAM 10.2: TREND ESTIMATION ',38,5,
cc     *               VNAME,VALUE )
cc      WX = 15.0
cc      WY =  4.5
cc      DX =200.0
cc      FN = N
cc      IY = 10
C
C  ...  Plot original data  ...
C
cc      CALL  MAXMIN( Y,N,YMIN1,YMAX1,DY )
cc      CALL  NEWPEN( 2 )
cc      CALL  PLOT( 3.0,16.0-SNGL(WY),-3 )
cc      CALL  SYMBOL( 0.25,SNGL(WY)+0.25,0.25,'TIME SERIES',0.0,11 )
cc      CALL  AXISXY( 0.0D0,0.0D0,WX,WY,0.0D0,FN,YMIN1,YMAX1,DX,DY,
cc     *              0.2D0,1,IY,2 )
cc      CALL  PLOTY( Y,N,YMIN1,YMAX1,WX,WY,1,1 )
C
C  ...   PLOT  TREND  ...
C
cc      CALL  PLOT( 0.0,-SNGL(WY+1.0),-3 )
cc      CALL  NEWPEN( 2 )
cc      CALL  SYMBOL( 0.25,SNGL(WY)+0.25,0.25,'TREND',0.0,5 )
cc      CALL  AXISXY( 0.0D0,0.0D0,WX,WY,0.0D0,FN,YMIN1,YMAX1,DX,DY,
cc     *              0.2D0,1,IY,2)
C
      DO 120 J=-1,1
      DO 100 I=1,N
cc  100 DATA(I) = XSS(1,I) + DSQRT( VSS(1,1,I)*SIG2 )*J
cxx  100 TREND(I,J) = XSS(1,I) + DSQRT( VSS(1,1,I)*SIG2 )*J
      TREND(I,J) = XSS(1,I) + DSQRT( VSS(1,1,I)*SIG2 )*J
  100 CONTINUE
cc      CALL  NEWPEN( 1 )
cc      IF(J.EQ.0)  CALL  NEWPEN( 2 )
cc      CALL  PLOTY( DATA,N,YMIN1,YMAX1,WX,WY,1,1 )
  120 CONTINUE
C
C  ...  Plot residuals  ...
C
      DO 180 J=4,4,3
      DO 160 I=1,N
cxx  160 DATA(I) = Y(I) - XSS(1,I)
      DATA(I) = Y(I) - XSS(1,I)
  160 CONTINUE
cc      CALL  MAXMIN( DATA,N,DMIN,DMAX,DDY )
cc      YMIN3 = DY*INT(DMIN/DY)
cc      YMAX3 = DY*INT(DMAX/DY)
cc      IF( DMIN.LT.YMIN3 )  YMIN3 = YMIN3-DY
cc      IF( DMAX.GT.YMAX3 )  YMAX3 = YMAX3+DY
cc      WY3 = WY*(YMAX3-YMIN3)/(YMAX1-YMIN1)
C
cc      CALL  PLOT( 0.0,-SNGL(WY3+1.0),-3 )
cc      CALL  NEWPEN( 2 )
cc      CALL  SYMBOL( 0.25,SNGL(WY3)+0.25,0.25,'RESIDUAL',0.0,8 )
cc      CALL  AXISXY( 0.0D0,0.0D0,WX,WY3,0.0D0,FN,YMIN3,YMAX3,DX,DY,
cc     *              0.2D0,1,IY,2)
cc      CALL  NEWPEN( 1 )
cc      IF(J.EQ.4)  CALL  NEWPEN( 2 )
cc      CALL  PLOTY2( DATA,N,YMIN3,YMAX3,WX,WY3,1,1,0.0D0 )
  180 CONTINUE
C
cc      CALL  PLOTE
C     call  plot( 0.0,0.0,999 )
C
      RETURN
      E N D
cc      SUBROUTINE  SETTRN( M,MJ,F,G,H,R )
      SUBROUTINE  SETTRN( M,F,G,H,R )
C
C  ...  State space model for trend estimation  ...
C
C     Input:
C       M:   Order of the trend
C       MJ:  Adjustable Dimension of F
C     Outputs:
C       F:   M*M Transition Matrix
C       G,H: M Vector
C       R:   Observational Noise (=1)
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  F(MJ,MJ), G(MJ), H(MJ)
cxx      DIMENSION  F(M,M), G(M), H(M)
C
      INTEGER :: M
      REAL(8) :: F(M,M), G(M), H(M), R
C
cc      DO 10  J=1,MJ
cxx      DO 10  J=1,M
cxx      G(J) = 0.0D0
cxx      H(J) = 0.0D0
cc      DO 10  I=1,MJ
cxx      DO 10  I=1,M
cxx   10 F(I,J) = 0.0D0
      F(1:M,1:M) = 0.0D0
      G(1:M) = 0.0D0
      H(1:M) = 0.0D0
C
      IF( M.EQ.1 )  THEN
        F(1,1) = 1.0D0
      END IF
      IF( M.EQ.2 )  THEN
        F(1,1) = 2.0D0
        F(1,2) =-1.0D0
        F(2,1) = 1.0D0
      END IF
      IF( M.EQ.3 )  THEN
        F(1,1) = 3.0D0
        F(1,2) =-3.0D0
        F(1,3) = 1.0D0
        F(2,1) = 1.0D0
        F(3,2) = 1.0D0
      END IF
      G(1) = 1.0D0
      H(1) = 1.0D0
      R    = 1.0D0
C
      RETURN
      E N D
cc      SUBROUTINE  FILTER( Y,XF,VF,F,G,H,Q,R,M,K,ISW,NS,NFE,NPE,MJ,NMAX,
      SUBROUTINE  FILTER( Y,XF,VF,F,G,H,Q,R,M,K,ISW,NS,NFE,NPE,
     *                    NDIM,OUTMIN,OUTMAX,VFS,VPS,XFS,XPS,FF,SIG2 )
C
C  ...  Kalman Filter (General Form, L=1)  ...
C
C     Inputs:
C        Y:      time series
C        NS:     Start position of filtering
C        NFE:    End position of filtering
C        NPE:    End position of prediction
C        XF:     Initial state vector
C        VF:     Initial covariance matrix
C        M:      Dimension of the state vector
C        K:      Dimension of the system noise
C        F:      M*M matrix
C        G:      M*K matrix
C        H:      M vector
C        Q:      K*K matrix, system noise covariance
C        R:      observation variance
C        ISW:    = 0;   R will be estimated
C                = 1;   R is given
C        MJ:     Adjustable dimension of XF, VF
C        NDIM:   Adjustable dimension of XFS, XPS, VFS, VPS
C                = 0   XF, XP, VF, VP are not stored
C                > 0   They are stored for smoothing
C        NMAX    Adjustable dimension of Y
C        OUTMIN: Lower limit for detecting outliers
C        OUTMAX: Upper limit for detecting outliers
C     Outputs:
C        VFS:    Covariance matrices of the filter
C        VPS:    Covariance matrices of the predictor
C        XFS:    Mean vectors of the filter
C        XPS:    Mean vectors of the predictor
C        FF:     Log likelihood
C        SIG2:   Estimated observational noise variance
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  Y(NMAX)
cc      DIMENSION  F(MJ,MJ), G(MJ,K), H(MJ), Q(K,K)
cc      DIMENSION  XF(MJ), VF(MJ,MJ), XP(40), VP(40,40)
cc      DIMENSION  XFS(MJ,NDIM),    XPS(MJ,NDIM)
cc      DIMENSION  VFS(MJ,MJ,NDIM), VPS(MJ,MJ,NDIM)
cc      DIMENSION  WRK(40,40), VH(40), GAIN(40)
cxx      DIMENSION  Y(NPE)
cxx      DIMENSION  F(M,M), G(M,K), H(M), Q(K,K)
cxx      DIMENSION  XF(M), VF(M,M), XP(M), VP(M,M)
cxx      DIMENSION  XFS(M,NDIM),    XPS(M,NDIM)
cxx      DIMENSION  VFS(M,M,NDIM), VPS(M,M,NDIM)
cxx      DIMENSION  WRK(M,M), WRK1(M,K), VH(M), GAIN(M)
C
      INTEGER :: M, K, ISW, NS, NFE, NPE, NDIM
      REAL(8) :: Y(NPE), XF(M), VF(M,M), F(M,M), G(M,K), H(M), Q(K,K),
     1           R, OUTMIN, OUTMAX, VFS(M,M,NDIM), VPS(M,M,NDIM),
     2           XFS(M,NDIM), XPS(M,NDIM), FF, SIG2
      REAL(8) :: XP(M), VP(M,M), WRK(M,M), WRK1(M,K), VH(M), GAIN(M),
     1           PI, SDET, SUM, PERR, PVAR
C
      DATA   PI  /3.1415926535D0/
C
      SIG2 = 0.0D0
      SDET = 0.0D0
      NSUM = 0
C
      DO 500  II=NS,NPE
C
C  ...  ONE STEP AHEAD PREDICTION  ...
C
      DO 20  I=1,M
      SUM = 0.0D0
      DO 10  J=1,M
cxx   10 SUM = SUM + F(I,J)*XF(J)
      SUM = SUM + F(I,J)*XF(J)
   10 CONTINUE
cxx   20 XP(I) = SUM
      XP(I) = SUM
   20 CONTINUE
C
cxx   DO 40  I=1,M
      DO 41 I=1,M
      DO 40  J=1,M
      SUM = 0.0D0
      DO 30  JJ=1,M
cxx   30 SUM = SUM + F(I,JJ)*VF(JJ,J)
      SUM = SUM + F(I,JJ)*VF(JJ,J)
   30 CONTINUE
cxx   40 WRK(I,J) = SUM
      WRK(I,J) = SUM
   40 CONTINUE
   41 CONTINUE
C
cxx      DO 60  I=1,M
      DO 61  I=1,M
      DO 60  J=1,M
      SUM = 0.0D0
      DO 50 JJ=1,M
cxx   50 SUM = SUM + WRK(I,JJ)*F(J,JJ)
      SUM = SUM + WRK(I,JJ)*F(J,JJ)
   50 CONTINUE
cxx   60 VP(I,J) = SUM
      VP(I,J) = SUM
   60 CONTINUE
   61 CONTINUE
C
cxx      DO 80  I=1,M
      DO 81  I=1,M
      DO 80  J=1,K
      SUM = 0.0D0
      DO 70 JJ=1,K
cxx   70 SUM = SUM + G(I,JJ)*Q(JJ,J)
      SUM = SUM + G(I,JJ)*Q(JJ,J)
   70 CONTINUE
cc   80 WRK(I,J) = SUM
cxx   80 WRK1(I,J) = SUM
      WRK1(I,J) = SUM
   80 CONTINUE
   81 CONTINUE
C
cxx      DO 100  I=1,M
      DO 101  I=1,M
      DO 100  J=1,M
      SUM = VP(I,J)
      DO 90  JJ=1,K
cc   90 SUM = SUM + WRK(I,JJ)*G(J,JJ)
cxx   90 SUM = SUM + WRK1(I,JJ)*G(J,JJ)
      SUM = SUM + WRK1(I,JJ)*G(J,JJ)
   90 CONTINUE
cxx  100 VP(I,J) = SUM
      VP(I,J) = SUM
  100 CONTINUE
  101 CONTINUE
C
C  ...  FILTERING  ...
C
      IF( Y(II).GT.OUTMIN.AND.Y(II).LT.OUTMAX.AND. II.LE.NFE ) THEN
C
      DO 210  I=1,M
      SUM = 0.0D0
      DO 200  J=1,M
cxx  200 SUM = SUM + VP(I,J)*H(J)
      SUM = SUM + VP(I,J)*H(J)
  200 CONTINUE
cxx  210 VH(I) = SUM
      VH(I) = SUM
  210 CONTINUE
C
      PERR = Y(II)
      PVAR = R
      DO 220  I=1,M
      PERR = PERR - H(I)*XP(I)
cxx  220 PVAR = PVAR + H(I)*VH(I)
      PVAR = PVAR + H(I)*VH(I)
  220 CONTINUE
C
      DO 250  I=1,M
cxx  250 GAIN(I) = VH(I)/PVAR
      GAIN(I) = VH(I)/PVAR
  250 CONTINUE
C
      DO 290  I=1,M
cxx  290 XF(I) = XP(I) + GAIN(I)*PERR
      XF(I) = XP(I) + GAIN(I)*PERR
  290 CONTINUE
C
cxx      DO 310  I=1,M
      DO 311  I=1,M
      DO 310  J=1,M
cxx  310 VF(I,J) = VP(I,J) - GAIN(I)*VH(J)
      VF(I,J) = VP(I,J) - GAIN(I)*VH(J)
  310 CONTINUE
  311 CONTINUE
C
      SIG2 = SIG2 + PERR**2/PVAR
      SDET = SDET + DLOG(PVAR)
      NSUM = NSUM + 1
C
C  ...  MISSING OBSERVATION  ...
C
      ELSE
cxx      DO 350  I=1,M
      DO 351  I=1,M
      XF(I) = XP(I)
      DO 350  J=1,M
cxx  350 VF(I,J) = VP(I,J)
      VF(I,J) = VP(I,J)
  350 CONTINUE
  351 CONTINUE
      END IF
C
C  ...  SAVE MEAN AND COVARIANCE  ...
C
      IF( NDIM.GT.1 )  THEN
cxx      DO 360  I=1,M
      DO 361  I=1,M
      XPS(I,II) = XP(I)
      XFS(I,II) = XF(I)
      DO 360  J=1,M
      VPS(I,J,II) = VP(I,J)
cxx  360 VFS(I,J,II) = VF(I,J)
      VFS(I,J,II) = VF(I,J)
  360 CONTINUE
  361 CONTINUE
      END IF
C
  500 CONTINUE
      SIG2 = SIG2/NSUM
      IF(ISW.EQ.0)  FF = -0.5D0*(NSUM*(DLOG(PI*2*SIG2) + 1) + SDET)
      IF(ISW.EQ.1)  FF = -0.5D0*(NSUM*(DLOG(PI*2) + SIG2) + SDET)
C
      RETURN
      E N D
cc      SUBROUTINE  ISTATE( M,MJ,XMEAN,XVAR,XF,VF )
      SUBROUTINE  ISTATE( M,XMEAN,XVAR,XF,VF )
C
C  ...  Initial state ...
C
C     Inputs:
C        M:     Dimension of the state vector
C        MJ:    Adjustable dimension of F
C        XMEAN: Mean value
C        XVAR:  Varaince
C     Outputs:
C         XF:   State vector, X(0|0)
C         VF:   State covarance matrix, V(0|0)
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  XF(MJ), VF(MJ,MJ)
cxx      DIMENSION  XF(M), VF(M,M)
C
      INTEGER :: M
      REAL(8) :: XMEAN, XVAR, XF(M), VF(M,M)
C
cc      DO 10  I=1,MJ
cc      DO 10  J=1,MJ
cxx      DO 10  I=1,M
cxx      DO 10  J=1,M
cxx   10 VF(I,J) = 0.0D0
      VF(1:M,1:M) = 0.0D0

C
      DO 20  I=1,M
cxx   20 XF(I) = XMEAN
      XF(I) = XMEAN
   20 CONTINUE
C
      DO 30  I=1,M
cxx   30 VF(I,I) = XVAR
      VF(I,I) = XVAR
   30 CONTINUE
C
      RETURN
      E N D
