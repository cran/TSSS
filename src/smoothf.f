C     PROGRAM 9.1  SMOOTH
      SUBROUTINE SMOOTHF(Y,N,M,K,F,G,H,Q,R,XF,VF,NFE,NPE,OUTMIN,
     *                   OUTMAX,NMISS,N0,NN,XSS,VSS,FLK,AIC)
C
      INCLUDE 'TSSS.h'
C
C  ...  Prediction and Interpolation of Time Series  ...
C
C     Inputs:
C        NFE:      End point of filtering
C        NPE:      End point of prediction
C        M:        Order of the AR model
C        TAU2:     Innovation variance of the AR model
C        AR(I):    AR coefficients (II=1,M)
C        OUTMIN,OUTMAX:  Lower and upper limits of observations
C        NMISS:    Number of missed intervals
C        N0(I):    Start position of missed intervals (I=1,NMISS)
C        NN(I):    Number of missed observations (I=1,NMISS)
C     Inputs required in subroutine READTS
C        TITLE:    Title of the data set
C        N:        Data length
C        Y(I):     Tiem series (I=1,N)
C     Parameters:
C        IDEV:     Input device specification
C        NMAX:     Adjustable dimension of Y, YMISS (NMAX >= N)
C        MJ:       Highest dimension of the state (MJ >= M)
C        K:        Dimension of the system noise
C        ISW:      =1 (R is specified)
C     @TEST.FILTER2:  SEP.08,1990, SEP.02,1992
C     MODIFIED  2/15/93
C
cc      PARAMETER( NMAX=1000,MJ=20,K=1,ISW=1,IDEV=1 )
cxxx      PARAMETER( ISW=1 )
      INTEGER, PARAMETER :: ISW=1
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  Y(NMAX), YMISS(NMAX)
cc      DIMENSION  F(MJ,MJ), G(MJ), H(MJ), Q(K,K)
cc      DIMENSION  XPS(MJ,NMAX), XFS(MJ,NMAX), XSS(MJ,NMAX)
cc      DIMENSION  VPS(MJ,MJ,NMAX), VFS(MJ,MJ,NMAX), VSS(MJ,MJ,NMAX)
cc      DIMENSION  XF(MJ), VF(MJ,MJ)
cc      DIMENSION  N0(10), NN(10), AR(MJ)
cxx      DIMENSION  Y(N,L), YMISS(N,L)
cxx      DIMENSION  F(M,M), G(M,K), H(L,M), Q(K,K), R(L,L)
cxx      DIMENSION  XPS(M,NPE), XFS(M,NPE), XSS(M,NPE)
cxx      DIMENSION  VPS(M,M,NPE), VFS(M,M,NPE), VSS(M,M,NPE)
cxx      DIMENSION  XF(M), VF(M,M)
cxx      DIMENSION  N0(NMISS), NN(NMISS)
cxx      DIMENSION  YMEAN(L), YVAR(L)
C
      INTEGER N, M, K, NFE, NPE, NMISS, N0(NMISS), NN(NMISS)
      DOUBLE PRECISION Y(N), F(M,M), G(M,K), H(M), Q(K,K), R, XF(M),
     1                 VF(M,M), OUTMIN, OUTMAX, XSS(M,NPE),
     2                 VSS(M,M,NPE), FLK, AIC
c local
      INTEGER I, J, MJ, NMAX, NDIM, NS
      DOUBLE PRECISION YMISS(N), XPS(M,NPE), XFS(M,NPE), VPS(M,M,NPE),
     1                 VFS(M,M,NPE), YMEAN, YVAR, FF, SIG2
C
      NMAX = N
      MJ = M
      NDIM = NPE
C
cc      READ( 5,* )  NFE, NPE
cc      READ( 5,* )  M
cc      READ( 5,* )  TAU2
cc      READ( 5,* )  (AR(I),I=1,M)
cc      READ( 5,* )  OUTMIN, OUTMAX
cc      READ( 5,* )  NMISS
cc      DO 10 I=1,NMISS
cc   10 READ( 5,* )  N0(I), NN(I)
      NS = 1
C
C  ...  Read Time Series  ...
C
cc      CALL  READTS( IDEV,Y,N )
      CALL  MOMENT( Y,N,YMEAN,YVAR )
C
C  ...  SET MISSING OBSERVATIONS  ...
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 20 I=1,N
cc   20 YMISS(I) = Y(I) - YMEAN
      YMISS(I) = Y(I) - YMEAN
   20 CONTINUE
cxx      DO 30 J=1,NMISS
      DO 31 J=1,NMISS
      DO 30 I=1,NN(J)
cc   30 YMISS(N0(J)+I-1) = OUTMIN
      YMISS(N0(J)+I-1) = OUTMIN
   30 CONTINUE
   31 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  ...  Set State Space Model  ...
C
cc      CALL  SETFGH( M,MJ,K,AR,TAU2,F,G,H,Q,R )
C
C  ...  Set Initial State  ...
C
cc      CALL  ISTATE( M,MJ,0.0D0,YVAR,XF,VF )
C
C  ...  Kalman Filtering  ...
C
cc      CALL  FILTER( YMISS,XF,VF,F,G,H,Q,R,M,K,ISW,NS,NFE,NPE,MJ,NMAX,
cc     *              NMAX,OUTMIN,OUTMAX,VFS,VPS,XFS,XPS,FF,SIG2 )
      CALL  MFILTER( YMISS,N,XF,VF,F,G,H,Q,R,M,K,ISW,NS,NFE,NPE,
     *               NDIM,OUTMIN,OUTMAX,VFS,VPS,XFS,XPS,FF,SIG2 )
C
C  ...  Fixed Interval Smoothing  ...
C
cc      CALL  SMOOTH( F,M,MJ,NDIM,NS,NFE,NPE,VFS,VPS,XFS,XPS,
cc     *              VSS,XSS )
      CALL SMOOTH( F,M,NDIM,NS,NFE,NPE,VFS,VPS,XFS,XPS,VSS,XSS )
      FLK = -FF
      AIC = -2*FF + 2*(M+1)
C
C  ...  Plot and Print  ...
C
C      CALL  PTSTAT( Y,M,TAU2,FLK,AIC,XSS,VSS,N,NFE,NPE,MJ,
C     *              YMEAN,NMISS,N0,NN )
cc      CALL  PRSTAT( M,TAU2,FLK,AIC,XSS,VSS,YMEAN,MJ,NFE,NPE,
cc     *              NMISS,N0,NN )
C
      RETURN
      E N D
cc      SUBROUTINE  FILTER( Y,XF,VF,F,G,H,Q,R,M,K,ISW,NS,NFE,NPE,MJ,NMAX,
      SUBROUTINE  MFILTER( Y,N,XF,VF,F,G,H,Q,R,M,K,ISW,NS,NFE,NPE,NDIM,
     *                     OUTMIN,OUTMAX,VFS,VPS,XFS,XPS,FF,SIG2 )
C
C  ...  Kalman Filter (General Form, L>1)  ...
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
cxx      DIMENSION  Y(N,L)
cxx      DIMENSION  F(M,M), G(M,K), H(L,M), Q(K,K)
cxx      DIMENSION  XF(M), VF(M,M), XP(M), VP(M,M)
cxx      DIMENSION  XFS(M,NDIM), XPS(M,NDIM)
cxx      DIMENSION  VFS(M,M,NDIM), VPS(M,M,NDIM)
cxx      DIMENSION  WRK(M,M), WRK1(M,K), WRK2(L), VH(M,L), GAIN(M,L)
cxx      DIMENSION  R(L,L), PVAR(L,L), PERR(L)
C
      INTEGER N, M, K, ISW, NS, NFE, NPE, NDIM
      DOUBLE PRECISION Y(N), XF(M), VF(M,M), F(M,M), G(M,K), H(M),
     1                 Q(K,K),R, OUTMIN, OUTMAX, VFS(M,M,NDIM),
     2                 VPS(M,M,NDIM),XFS(M,NDIM), XPS(M,NDIM), FF, SIG2
c local
      INTEGER I, II, J, JJ, NSUM
      DOUBLE PRECISION XP(M), VP(M,M), WRK(M,M), WRK1(M,K), VH(M),
     1                 GAIN(M), PVAR, PERR, PI, SDET, SUM
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
cxx      DO 40  I=1,M
      DO 41  I=1,M
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
cc      IF( Y(II).GT.OUTMIN.AND.Y(II).LT.OUTMAX.AND. II.LE.NFE ) THEN
      IF( II.LE.NFE ) THEN
         IF  ( Y(II).GT.OUTMIN .AND. Y(II).LT.OUTMAX )  THEN
C
      DO 210  I=1,M
      SUM = 0.0D0
cc      DO 200  J=1,M
cc  200 SUM = SUM + VP(I,J)*H(J)
cc  210 VH(I) = SUM
      DO 200  JJ=1,M
cxx  200 SUM = SUM + VP(I,JJ)*H(J,JJ)
      SUM = SUM + VP(I,JJ)*H(JJ)
  200 CONTINUE
cxx  210 VH(I,J) = SUM
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
cc  250 GAIN(I) = VH(I)/PVAR
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
cc  310 VF(I,J) = VP(I,J) - GAIN(I)*VH(J)
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
      DO 352  I=1,M
      XF(I) = XP(I)
      DO 350  J=1,M
cxx  350 VF(I,J) = VP(I,J)
      VF(I,J) = VP(I,J)
  350 CONTINUE
  352 CONTINUE
      END IF
C
      ELSE
cxx      DO 351  I=1,M
      DO 353  I=1,M
      XF(I) = XP(I)
      DO 351  J=1,M
cxx  351 VF(I,J) = VP(I,J)
      VF(I,J) = VP(I,J)
  351 CONTINUE
  353 CONTINUE
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
cc      SIG2 = SIG2/NSUM
cxxx      SIG2 = SIG2/DFLOAT(NSUM)
      SIG2 = SIG2/DBLE(NSUM)
      IF(ISW.EQ.0)  FF = -0.5D0*(NSUM*(DLOG(PI*2*SIG2) + 1) + SDET)
      IF(ISW.EQ.1)  FF = -0.5D0*(NSUM*(DLOG(PI*2) + SIG2) + SDET)
C
      RETURN
      E N D
