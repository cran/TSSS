C     PROGRAM 12.1  SEASON
      SUBROUTINE SEASONF( Y0,N,M1,M2,M3,M4,PERIOD,JYEAR,MONTH,
     & TAU2,NS,NFE,NPE,AR,LOGT,IOPT,OUTMIN,OUTMAX,NMAX,MJ,
     & FF,OVAR,AIC,XSS,VSS,DEFF,IER1,IER2 ) 
C
      INCLUDE 'TSSS_f.h'
C
C  ...  SEASONAL ADJUSTMENT BY STATE SPACE MODELING  ...
C
C     Inputs:
C        M1:      Trend Order (M1=0,1,2,3)
C        M2:      Seasonal Order (M1=0,1,2)
C        M3:      AR Order
C        M4:      Trading day effect (M4=0 or 6)
C        PERIOD:  Number of seasons in one period
C                 (=12 for monthly data, =4 for quarterly data)
C        IPARAM:  =0    Use defalt initail values
C                 =1    Read intial values
C        OUTMIN:  Lower limit of the observations
C        OUTMAX:  Upper limit of the observations
C        JYEAR:   Starting year of the data (e.g. 1967)
C        MONTH:   Starting month of the data (e.g. 1)
C        TAU2(I): System noise variances   (I=1,4)
C        AR(I):   AR coefficients (I=1,M3)
C        NS:      Start position of filtering
C        NFE:     End point of filtering
C        NPE:     End point of prediction
C        IOPT:    = 1     To get MLE
C     Parameters:
C        NMAX:    Adjustable dimension of Y
C        MJ:      Adjustable dimension of XF, VF, etc.
C        MAXM:    Adjustable dimension of A, B, C
C        NC:      Number of components
C     @TEST.FILTER2O    NOV.29,1990, SEP.02,1992
C
cc      PARAMETER( NMAX=160,MJ=20,MAXM=22,NC=10)
      PARAMETER( MAXM=22,NC=9 )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      INTEGER*4  PERIOD
cc      DIMENSION  M(10), TAU2(10), TAU20(10)
cc      DIMENSION  XMEAN(10), XVAR(MAXM,10), AA(20), AR(10), PAR(10)
cc      DIMENSION  A(MAXM,10), B(MAXM,10), C(MAXM,10), Q(10,10)
cxx      DIMENSION  M(NC), TAU2(4)
cxx      DIMENSION  XMEAN(NC), XVAR(MAXM,NC), AA(M3+3), AR(M3), PAR(M3)
cxx      DIMENSION  A(MAXM,NC), B(MAXM,NC), C(MAXM,NC), Q(NC,NC)
cxx      DIMENSION  XPS(MJ,NMAX), XFS(MJ,NMAX), XSS(MJ,NMAX)
cxx      DIMENSION  VPS(MJ,MJ,NMAX), VFS(MJ,MJ,NMAX), VSS(MJ,MJ,NMAX)
cc      DIMENSION  XF(MJ), VF(MJ,MJ), MTYPE(10)
cxx      DIMENSION  XF(MJ), VF(MJ,MJ), MTYPE(NC)
cc      COMMON  /C92827/  M1, M2, M3, M4, PERIOD
cc      COMMON  /C92825/  OUTMIN, OUTMAX
cc      COMMON  /C92826/  Y(NMAX), REG(NMAX,7)
cxx      DIMENSION  Y(N), REG(NMAX,7)
cc      COMMON  /C92907/  ALIMIT
cc      COMMON  /CMMODL/  M, XMEAN, XVAR, N, NS, NFE, NPE
cc      COMMON  / CCC /  ISW, IPR, ISMT, IDIF
cxx      DIMENSION  Y0(N), DEFF(NMAX)
C
      INTEGER :: N, M1, M2, M3, M4, PERIOD, JYEAR, MONTH, NS, NFE, NPE,
     1           LOGT, IOPT, NMAX, MJ, IER1, IER2
      REAL(8) :: Y0(N), TAU2(4), AR(M3), OUTMIN, OUTMAX, FF, OVAR, AIC,
     1           XSS(MJ,NMAX), VSS(MJ,MJ,NMAX), DEFF(NMAX)
      INTEGER :: M(NC), MTYPE(NC)
      REAL(8) :: XMEAN(NC), XVAR(MAXM,NC), AA(M3+3), PAR(M3),
     1           A(MAXM,NC), B(MAXM,NC), C(MAXM,NC), Q(NC,NC),
     2           XPS(MJ,NMAX), XFS(MJ,NMAX), VPS(MJ,MJ,NMAX),
     3           VFS(MJ,MJ,NMAX), XF(MJ), VF(MJ,MJ), Y(N), REG(NMAX,7),
     4           SIG2, ALIMIT
C
      EXTERNAL  FFSEAS
C
      XSS = 0.0D0
      VSS = 0.0D0
      DEFF = 0.0D0
C
C  ...  Read Time Series  ...
C
cc      CALL  READTS( 1,Y,N )
      DO 100 I = 1,N
cxx  100 Y(I) = Y0(I)
      Y(I) = Y0(I)
  100 CONTINUE
C
C  ...  Read Model Orders  ...
C
cc      READ( 5,* )  M1, M2, M3, M4, PERIOD
C
C  ...  Set Defalt Parameters  ...
C
cc      IPR  = 7
      SIG2 = 1.0D0
      ALIMIT = 0.9D0
cc      CALL  SPARAM( M1,M2,M3,M4,N,TAU2,AR,OUTMIN,OUTMAX,NS,NFE,NPE,
cc     *              LOGT,IOPT )
C
cc      READ( 5,* )  IPARAM
cc      IF( IPARAM.EQ.1 )  THEN
cc         READ( 5,* )  (TAU2(I),I=1,3)
cc         READ( 5,* )  NS, NFE, NPE
cc         READ( 5,* )  LOGT, IOPT
cc         READ( 5,* )  OUTMIN, OUTMAX
cc         IF(M3.GT.0)  READ( 5,* )  (AR(I),I=1,M3)
cc      END IF
cc      WRITE(6,*) M1,M2,M3,M4,PERIOD
cc      WRITE(6,*) (TAU2(I),I=1,3)
cc      WRITE(6,*) (AR(I),I=1,M3)
C
cc      IF( M4.GT.0 )  READ( 5,* )  JYEAR, MONTH
C
C  ...  Log Transformation  ...
C
      IF( LOGT.EQ.1 )  THEN
         IER1 = -1
         DO 5 I=1,N
cc    5    Y(I) = DLOG10( Y(I) )
            IF( Y(I).LE.0 ) RETURN
cxx    5       Y(I) = DLOG10( Y(I) )
            Y(I) = DLOG10( Y(I) )
    5    CONTINUE
      END IF
      IER1 = 0
C
C  ...  Trading day factor  ...
C
cc      IF( M4.GT.0 )  CALL  TRADE( JYEAR,MONTH,N,NMAX,REG )
      IF( M4.GT.0 )  THEN
         CALL  TRADE( JYEAR,MONTH,NMAX,NMAX,REG )
cc      DO 10 I=1,N
cxx      DO 10 I=1,NMAX
      DO 11 I=1,NMAX
      DO 10 J=1,6
cxx   10 REG(I,J) = REG(I,J) - REG(I,7)
      REG(I,J) = REG(I,J) - REG(I,7)
   10 CONTINUE
   11 CONTINUE
      END IF
C
      CALL  PARCOR( AR,M3,PAR )
      NP = ID(M1) + ID(M2) + ID(M3)
      DO 20 I=1,NP
cxx   20 AA(I) = DLOG( TAU2(I)/(1.0D0-TAU2(I)) )
      AA(I) = DLOG( TAU2(I)/(1.0D0-TAU2(I)) )
   20 CONTINUE
c----------    2013/07/02
      ier2 = 0
      DO 30 I=1,M3
         if( ALIMIT .le. PAR(I) ) ier2 = -1
cxx   30 AA(NP+I) = DLOG( (ALIMIT+PAR(I))/(ALIMIT-PAR(I)) )
      AA(NP+I) = DLOG( (ALIMIT+PAR(I))/(ALIMIT-PAR(I)) )
   30 CONTINUE
      if( ier2 .eq. -1 ) return
c----------
cc      DO 40 I=1,NC
cc   40 TAU20(I) = TAU2(I)
C
C  ...  Maximum Likelihood Method  ...
C
cc      IF( IOPT.EQ.1 )  CALL  DAVIDN( FFSEAS,AA,NP+M3,2 )
      IF( IOPT.EQ.1 )  THEN
         CALL  DAVIDN1( FFSEAS,AA,NP+M3,2,
     *   Y,N,M1,M2,M3,M4,PERIOD,REG,OUTMIN,OUTMAX,ALIMIT,ISW,IDIF,
     *   M,XMEAN,XVAR,NS,NFE,NPE,NMAX,MJ,MAXM,NC,IER2 )
         if( ier2.ne.0 ) return
      END IF
C
      DO 50 I=1,NP
cxx   50 TAU2(I) = DEXP( AA(I) )/(1.0D0 + DEXP( AA(I) ) )
      TAU2(I) = DEXP( AA(I) )/(1.0D0 + DEXP( AA(I) ) )
   50 CONTINUE
      DO 60 I=1,M3
cxx   60 PAR(I)  = ALIMIT*(DEXP(AA(NP+I))-1.0D0)/(DEXP(AA(NP+I))+1.0D0)
      PAR(I)  = ALIMIT*(DEXP(AA(NP+I))-1.0D0)/(DEXP(AA(NP+I))+1.0D0)
   60 CONTINUE
      CALL  ARCOEF( PAR,M3,AR )
C
C  ...  The Maxmum Likelihood Model  ...
C
cc      CALL  SETABC( M,L,NC,MTYPE,TAU2,MAXM,MM,AR,Y,N,A,B,C,Q,XMEAN,
cc     *              XVAR )
      CALL  SETABC1( M1,M2,M3,M4,PERIOD,M,L,NC,MTYPE,TAU2,MAXM,MM,AR,
     *         Y,N,A,B,C,Q,XMEAN,XVAR,ier2 )
      if( ier2.ne.0 ) return
      CALL  ISTAT1( L,M,MJ,MAXM,A,XMEAN,XVAR,XF,VF )
C
C  ...  The Kalman Filter  ... 
C
cc      CALL  FILTR1( Y,XF,VF,A,B,C,Q,SIG2,REG,MTYPE,NMAX,M,MAXM,NC,L,
cc     *      NS,NFE,NPE,MJ,N, OUTMIN,OUTMAX,VFS,VPS,XFS,XPS,FF,OVAR )
      CALL  FILTR1( Y,XF,VF,A,B,C,Q,SIG2,REG,MTYPE,NMAX,M,MAXM,NC,L,
     *      NS,NFE,NPE,MJ,N, OUTMIN,OUTMAX,VFS,VPS,XFS,XPS,FF,OVAR )
C
C  ...  Fixed Interval Smoothing  ...
C
cc      CALL  SMOTH1( A,M,MAXM,L,NS,NFE,NPE,MJ,
      CALL  SMOTH1( A,M,MAXM,L,NS,NFE,NPE,NMAX,MJ,
     *              VFS,VPS,VSS,XFS,XPS,XSS )
      AIC = -2*FF + 2*(NP+M3+MM)
C
C  ...  Print out and draw estimates  ...
C
cc      CALL  PRSEAS( M1,M2,M3,M4,PERIOD,TAU2,TAU20,OVAR,FF,AIC,AR,
cc     *              XSS,N,REG,NMAX,MJ )
cxx      CALL  PRSEAS( M1,M2,M3,M4,PERIOD,XSS,DEFF,N,REG,NMAX,MJ )
      CALL  PRSEAS( M1,M2,M3,M4,PERIOD,XSS,DEFF,REG,NMAX,MJ )
C      CALL  PTSEAS( Y,XSS,VSS,N,NPE,MJ,M1,M2,M3,M4,
C     *              PERIOD,TAU2,OVAR,FF,AIC,REG,NMAX )
C
cc      STOP
      RETURN
      E N D
cc      SUBROUTINE  PRSEAS( M1,M2,M3,M4,LPER,TAU2,TAU20,OVAR,FF,AIC,
cc     *                    AR,XSS,N,REG,NMAX,MJ )
cxx      SUBROUTINE  PRSEAS( M1,M2,M3,M4,LPER,XSS,DEFF,N,REG,NMAX,MJ )
      SUBROUTINE  PRSEAS( M1,M2,M3,M4,LPER,XSS,DEFF,REG,NMAX,MJ )
C
C  ...  Print out Estimates  ...
C
C     Inputs:
C        M1-M4:  Model orders
C        LPER:   Period
C        TAU2:   Variance of the system noise
C        TAU20:  Variance of the system noise (initial estimat)
C        OVAR:   Variance of the observational noise
C        FF:     Log-Likelihood of the mode
C        AIC:    AIC of the model
C        AR:     AR coefficients
C        XSS:    Smoothed state vector
C        N:      Data length
C        REG:    Number of trading days
C        NMAX:   Adjustable dimension of REG
C        MJ:     Adjustable dimension of XSS
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      CHARACTER  TITLE*72
cc      DIMENSION  TAU2(*), TAU20(*), AR(*)
cxx      DIMENSION  XSS(MJ,NMAX), REG(NMAX,7)
cc      COMMON     /CMDATA/  TITLE
cxx      DIMENSION  DEFF(NMAX)
C
      INTEGER :: M1, M2, M3, M4, LPER, NMAX, MJ 
      REAL(8) :: XSS(MJ,NMAX), DEFF(NMAX), REG(NMAX,7)
      REAL(8) :: SUM
C
      NP = ID(M1) + ID(M2) + ID(M3)
      M12  = M1 + M2*(LPER-1)
      M123 = M12 + M3
cc      WRITE(6,600)
cc      WRITE(6,610) TITLE
cc      WRITE(6,620) M1, M2, M3, M4, LPER
cc      WRITE(6,630)
cc      WRITE(6,640)
cc      WRITE(6,650) (TAU20(I),I=1,NP)
cc      WRITE(6,660)
cc      WRITE(6,640)
cc      WRITE(6,650) (TAU2(I),I=1,NP)
cc      IF( M3.GT.0 )  THEN
cc      WRITE(6,690)
cc      WRITE(6,700) (AR(I),I=1,M3)
cc      END IF
cc      WRITE(6,680) OVAR, FF, AIC
cc      IF( M1.GT.0 )  THEN
cc        WRITE(6,720)
cc        WRITE(6,710)  (XSS(1,I),I=1,N)
cc      END IF
cc      IF( M2.GT.0 )  THEN
cc        WRITE(6,730)
cc        WRITE(6,710)  (XSS(M1+1,I),I=1,N)
cc      END IF
cc      IF( M3.GT.0 )  THEN
cc        WRITE(6,740)
cc        WRITE(6,710)  (XSS(M12+1,I),I=1,N)
cc      END IF
      IF( M4.GT.0 )  THEN
cc        DO 210 I=1,N
        DO 210 I=1,NMAX
        SUM = 0.0D0
        DO 200 J=1,6
cxx  200   SUM = SUM + XSS(M123+J,I)*REG(I,J)
        SUM = SUM + XSS(M123+J,I)*REG(I,J)
  200   CONTINUE
cc  210   REG(I,7) = SUM
cxx  210 DEFF(I) = SUM
      DEFF(I) = SUM
  210 CONTINUE
cc        WRITE(6,750)
cc        WRITE(6,710)  (REG(I,7),I=1,N)
      END IF
      RETURN
cxx  600 FORMAT( 1H1,'PROGRAM 10.1:  SEASONAL ADJUSTMENT' )
cxx  610 FORMAT( 1H ,A72 )
cxx  620 FORMAT( 1H ,'M1 =',I2,3X,'M2 =',I2,3X,'M3 =',I2,3X,'M4 =',I2,
cxx     *        'PERIOD =',I3 )
cxx  630 FORMAT( 1H ,5X,'***  INITIAL ESTIMATE  ***' )
cc  640 FORMAT( 1H ,'TAU2(I),I=1,NC:' )
cxx  640 FORMAT( 1H ,'TAU2(I),I=1,NP:' )
cxx  650 FORMAT( 1H ,5D13.6 )
cxx  660 FORMAT( 1H ,5X,'***  FINAL ESTIMATE  ***' )
cxx  680 FORMAT( 1H ,'SIG2 =',D13.6,3X,'LOG-LIKELIHOOD =',F13.4,3X,
cxx     *            'AIC =',F13.4 )
cxx  690 FORMAT( 1H ,'AR COEFFICIENTS' )
cxx  700 FORMAT( 1H ,5F10.6 )
cxx  710 FORMAT( 1H ,6F12.5 )
cxx  720 FORMAT( 1H ,'TREND COMPONENT' )
cxx  730 FORMAT( 1H ,'SEASONAL COMPONENT' )
cxx  740 FORMAT( 1H ,'AR COMPONENT' )
cxx  750 FORMAT( 1H ,'TRADING DAY EFFECT' )
      E N D
      SUBROUTINE  FILTR1( Y,XF,VF,A,B,C,Q,SIG2,REG,MTYPE,NMAX,M,MMAX,
     *                    NCM,NC,NS,N,NPE,MJ,NDIM,OUTMIN,OUTMAX,
     *                    VFS,VPS,XFS,XPS,FF,OVAR )
C
C  ...  Kalman filter  ...
C
C     Inputs:
C        Y:      time series
C        N:      data length
C        NS:     Start position of filtering
C        NE:     End position of prediction
C        XF:     Initial state vector
C        VF:     Initial covariance matrix
C        A:      Parameters of the matrix F
C        B:      Parameters of the matrix G
C        C:      Parameters of the matrix H
C        Q:      Covariance matrix of the system noise
C        M:      Dimension of the state vector
C        MMAX:   Upper limit of the dimension of each component
C        NC:     Number of components
C        MJ:     Adjustable dimension of XF, VF
C        NDIM:   Adjustable dimension of XFS, XPS, VFS, VPS
C                = 0   XF, XP, VF, VP are not stored
C                > 0   They are stored for smoothing
C        OUTMIN: Lower limit for detecting outliers
C        OUTMAX: Upper limit for detecting outliers
C        SIG2:   Variance of the observational noise
C     Outputs:
C        VFS:    Covariance matrices of the filter
C        VPS:    Covariance matrices of the predictor
C        XFS:    Mean vectors of the filter
C        XPS:    Mean vectors of the predictor
C        FF:     Log likelihood
C        OVAR:   Estimated variance
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  Y(NMAX), REG(NMAX,*)
cxx      DIMENSION  Y(NDIM), REG(NMAX,7)
cxx      DIMENSION  A(MMAX,NCM), B(MMAX,NCM), C(MMAX,NCM), M(NCM)
cc      DIMENSION  XF(MJ), VF(MJ,MJ), XP(40), VP(40,40)
cc      DIMENSION  XFS(MJ,N), XPS(MJ,N), Q(NCM,NCM)
cc      DIMENSION  VFS(MJ,MJ,NDIM), VPS(MJ,MJ,NDIM)
cc      DIMENSION  WRK(40,40), VH(40), GAIN(40)
cc      DIMENSION  I0(10), MTYPE(NCM)
cxx      DIMENSION  XF(MJ), VF(MJ,MJ), XP(MJ), VP(MJ,MJ)
cxx      DIMENSION  XFS(MJ,NMAX), XPS(MJ,NMAX), Q(NCM,NCM)
cxx      DIMENSION  VFS(MJ,MJ,NMAX), VPS(MJ,MJ,NMAX)
cxx      DIMENSION  WRK(MJ,MJ), VH(MJ), GAIN(MJ)
cxx      DIMENSION  I0(NCM), MTYPE(NCM)
C
      INTEGER :: NMAX, MMAX, NCM, NC, NS, N, NPE, MJ, NDIM, MTYPE(NCM),
     1           M(NCM)
      REAL(8) :: Y(NDIM), XF(MJ), VF(MJ,MJ), A(MMAX,NCM), B(MMAX,NCM),
     1           C(MMAX,NCM), Q(NCM,NCM), SIG2, REG(NMAX,7), OUTMIN,
     2           OUTMAX, VFS(MJ,MJ,NMAX), VPS(MJ,MJ,NMAX), XFS(MJ,NMAX),
     3           XPS(MJ,NMAX), FF, OVAR
      INTEGER :: I0(NCM)
      REAL(8) :: XP(MJ), VP(MJ,MJ), WRK(MJ,MJ), VH(MJ), GAIN(MJ), PI,
     1           SDET, SUM, PVAR, PERR
C
      DATA   PI  /3.1415926535D0/
C
      OVAR = 0.0D0
      SDET = 0.0D0
      NSUM = 0
      I0(1) = 0
      DO 10 I=2,NC
cxx   10 I0(I) = I0(I-1) + M(I-1)
      I0(I) = I0(I-1) + M(I-1)
   10 CONTINUE
      MM = I0(NC) + M(NC)
C
      DO 300  II=NS,NPE
C
C  ...  ONE STEP AHEAD PREDICTION  ...
C
cxx      DO 100  L=1,NC
      DO 101  L=1,NC
      XP(I0(L)+M(L)) = A(M(L),L)*XF(I0(L)+1)
      DO 100  I=1,M(L)-1
cxx  100 XP(I0(L)+I) = A(I,L)*XF(I0(L)+1) + XF(I0(L)+I+1)
      XP(I0(L)+I) = A(I,L)*XF(I0(L)+1) + XF(I0(L)+I+1)
  100 CONTINUE
  101 CONTINUE
C
cxx      DO 110  J=1,MM
cxx      DO 110  L=1,NC
      DO 112  J=1,MM
      DO 111  L=1,NC
      WRK(I0(L)+M(L),J) = A(M(L),L)*VF(I0(L)+1,J)
      DO 110  I=1,M(L)-1
cxx  110 WRK(I0(L)+I,J) = A(I,L)*VF(I0(L)+1,J) + VF(I0(L)+I+1,J)
      WRK(I0(L)+I,J) = A(I,L)*VF(I0(L)+1,J) + VF(I0(L)+I+1,J)
  110 CONTINUE
  111 CONTINUE
  112 CONTINUE
C
cxx      DO 120  I=1,MM
cxx      DO 120  L=1,NC
      DO 122  I=1,MM
      DO 121  L=1,NC
      VP(I,I0(L)+M(L)) = WRK(I,I0(L)+1)*A(M(L),L)
      DO 120  J=1,M(L)-1
cxx  120 VP(I,I0(L)+J) = WRK(I,I0(L)+1)*A(J,L) + WRK(I,I0(L)+J+1)
      VP(I,I0(L)+J) = WRK(I,I0(L)+1)*A(J,L) + WRK(I,I0(L)+J+1)
  120 CONTINUE
  121 CONTINUE
  122 CONTINUE
C
cxx      DO 130  J=1,NC
cxx      DO 130  L=1,NC
      DO 132  J=1,NC
      DO 131  L=1,NC
      DO 130  I=1,M(L)
cxx  130 WRK(I0(L)+I,J) = B(I,L)*Q(L,J)
      WRK(I0(L)+I,J) = B(I,L)*Q(L,J)
  130 CONTINUE
  131 CONTINUE
  132 CONTINUE
C
cxx      DO 140  I=1,MM
cxx      DO 140  L=1,NC
      DO 142  I=1,MM
      DO 141  L=1,NC
      DO 140  J=1,M(L)
cxx  140 VP(I,I0(L)+J) = VP(I,I0(L)+J) + WRK(I,L)*B(J,L)
      VP(I,I0(L)+J) = VP(I,I0(L)+J) + WRK(I,L)*B(J,L)
  140 CONTINUE
  141 CONTINUE
  142 CONTINUE
C
C  ...  FILTERING  ...
C
cc      IF( Y(II).GT.OUTMIN .AND. Y(II).LT.OUTMAX .AND. II.LE.N )  THEN
      IF( II.LE.N .AND. (Y(II).GT.OUTMIN .AND. Y(II).LT.OUTMAX) )  THEN
C
      DO 210  I=1,MM
      SUM = 0.0D0
      DO 205  L=1,NC
      IF( MTYPE(L).EQ.0 )  THEN
         DO 200  J=1,M(L)
cxx  200    SUM = SUM + VP(I,I0(L)+J)*C(J,L)
         SUM = SUM + VP(I,I0(L)+J)*C(J,L)
  200    CONTINUE
      ELSE
         SUM = SUM + VP(I,I0(L)+1)*REG(II,MTYPE(L))
      END IF
  205 CONTINUE
cxx  210 VH(I) = SUM
      VH(I) = SUM
  210 CONTINUE
C
      PVAR = SIG2
      DO 225 L=1,NC
      IF( MTYPE(L).EQ.0 )  THEN
         DO 220 J=1,M(L)
cxx  220    PVAR = PVAR + VH(I0(L)+J)*C(J,L)
         PVAR = PVAR + VH(I0(L)+J)*C(J,L)
  220    CONTINUE
      ELSE
         PVAR = PVAR + VH(I0(L)+1)*REG(II,MTYPE(L))
      END IF
  225 CONTINUE
c------------------   2013/07/01
      if( PVAR.LE.0 ) GO TO 400
c------------------
C
      DO 230  I=1,MM
cxx  230 GAIN(I) = VH(I)/PVAR
      GAIN(I) = VH(I)/PVAR
  230 CONTINUE
C
      PERR = Y(II)
      DO 245 L=1,NC
      IF( MTYPE(L).EQ.0 )  THEN
         DO 240 J=1,M(L)
cxx  240    PERR = PERR - C(J,L)*XP(I0(L)+J)
         PERR = PERR - C(J,L)*XP(I0(L)+J)
  240    CONTINUE
      ELSE
         PERR = PERR - REG(II,MTYPE(L))*XP(I0(L)+1)
      END IF
  245 CONTINUE
C
      DO 250  I=1,MM
cxx  250 XF(I) = XP(I) + GAIN(I)*PERR
      XF(I) = XP(I) + GAIN(I)*PERR
  250 CONTINUE
C
cxx      DO 260  J=1,MM
      DO 261  J=1,MM
      DO 260  I=1,MM
cxx  260 VF(I,J) = VP(I,J) - GAIN(I)*VH(J)
      VF(I,J) = VP(I,J) - GAIN(I)*VH(J)
  260 CONTINUE
  261 CONTINUE
C
      OVAR = OVAR + PERR**2/PVAR
      SDET = SDET + DLOG(PVAR)
      NSUM = NSUM + 1
C
C  ...  MISSING OBSERVATION  ...
C
      ELSE
cxx      DO 270  I=1,MM
      DO 271  I=1,MM
      XF(I) = XP(I)
      DO 270  J=1,MM
cxx  270 VF(I,J) = VP(I,J)
      VF(I,J) = VP(I,J)
  270 CONTINUE
  271 CONTINUE
      END IF
C
C  ...  SAVE MEAN AND COVARIANCE  ...
C
      IF( NDIM.GT.1 )  THEN
cxx      DO 280  I=1,MM
      DO 281  I=1,MM
      XPS(I,II) = XP(I)
      XFS(I,II) = XF(I)
      DO 280  J=1,MM
      VPS(I,J,II) = VP(I,J)
cxx  280 VFS(I,J,II) = VF(I,J)
      VFS(I,J,II) = VF(I,J)
  280 CONTINUE
  281 CONTINUE
      END IF
C
  300 CONTINUE
      OVAR = OVAR/NSUM
      FF = -0.5D0*(NSUM*DLOG(PI*2*OVAR) + SDET + NSUM)
C     WRITE(6,*)  'FF =',FF, OVAR, SDET, NSUM
C
      RETURN
c------------------   2013/07/01
  400 FF = -1.0D30
      RETURN
c------------------
      E N D
      SUBROUTINE  ISTAT1( NC,M,MJ,MAXM,A,XMEAN,XVAR,XF,VF )
C
C  ...  Initial state ...
C
C     Inputs:
C        NC:    Number of components
C        M:     Dimension of the state vector
C        MJ:    Adjustable dimension of F
C        MAXM:  Adjustable dimension of A
C        A:     Parameter of matrix F
C        XMEAN: Mean value
C        XVAR:  Varaince
C     Outputs:
C         XF:   State vector, X(0|0)
C         VF:   State covarance matrix, V(0|0)
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  M(NC), A(MAXM,NC), XMEAN(NC), XVAR(MAXM,NC)
cc      DIMENSION  XF(MJ), VF(MJ,MJ), I0(10)
cxx      DIMENSION  XF(MJ), VF(MJ,MJ), I0(NC)
C
      INTEGER :: NC, M(NC), MJ, MAXM
      REAL(8) :: A(MAXM,NC), XMEAN(NC), XVAR(MAXM,NC), XF(MJ), VF(MJ,MJ)
      INTEGER :: I0(NC)
      REAL(8) :: SUM
C
      I0(1) = 0
      DO 5 I=2,NC
cxx    5 I0(I) = I0(I-1) + M(I-1)
      I0(I) = I0(I-1) + M(I-1)
    5 CONTINUE
cxx      DO 10  I=1,MJ
cxx      DO 10  J=1,MJ
cxx   10 VF(I,J) = 0.0D0
      VF(1:MJ,1:MJ) = 0.0D0
C
      DO 70  L=1,NC
C
      XF(I0(L)+1) = XMEAN(L)
      SUM = 0.0D0
      DO 20  I=M(L),2,-1
      SUM = SUM + A(I,L)
cxx   20 XF(I0(L)+I) = SUM*XMEAN(L)
      XF(I0(L)+I) = SUM*XMEAN(L)
   20 CONTINUE
C
      VF(I0(L)+1,I0(L)+1) = XVAR(1,L)
      DO 40 I=2,M(L)
      SUM = 0.0D0
      DO 30 J=I,M(L)
cxx   30 SUM = SUM + A(J,L)*XVAR(J,L)
      SUM = SUM + A(J,L)*XVAR(J,L)
   30 CONTINUE
      VF(I0(L)+1,I0(L)+I) = SUM
cxx   40 VF(I0(L)+I,I0(L)+1) = SUM
      VF(I0(L)+I,I0(L)+1) = SUM
   40 CONTINUE
cxx      DO 60 I=2,M(L)
      DO 61 I=2,M(L)
      DO 60 J=I,M(L)
      SUM = 0.0D0
cxx      DO 50 I1=I,M(L)
      DO 51 I1=I,M(L)
      DO 50 J1=J,M(L)
cxx   50 SUM = SUM + A(I1,L)*A(J1,L)*XVAR(IABS(J1-J-I1+I)+1,L)
      SUM = SUM + A(I1,L)*A(J1,L)*XVAR(IABS(J1-J-I1+I)+1,L)
   50 CONTINUE
   51 CONTINUE
      VF(I0(L)+I,I0(L)+J) = SUM
cxx   60 VF(I0(L)+J,I0(L)+I) = SUM
      VF(I0(L)+J,I0(L)+I) = SUM
   60 CONTINUE
   61 CONTINUE
   70 CONTINUE
C
      RETURN
      E N D
cc      SUBROUTINE  SETABC( M,L,NC,MTYPE,TAU2,MAXM,MM,AR,Y,N,A,B,C,Q,
cc     *                    XMEAN,XVAR )
      SUBROUTINE  SETABC1( M1,M2,M3,M4,MPER,M,L,NC,MTYPE,TAU2,MAXM,
     *                                MM,AR,Y,N,A,B,C,Q,XMEAN,XVAR,ier )
C
C  ...  State space model for Seasonal Adjustment  ...
C
C     Input:
C       M:     Order of the trend
C       L:     Number of components
C       NC:    Adjustable dimension
C       MTYPE: Model type
C       TAU2:  System noise variance
C       MAXM:  Adjustable dimension of A, B, C
C       MM:    State dimension
C       AR:    AR coefficient
C       Y(I):  Time series
C       N:     Data length
C     Outputs:
C       A,B,C: Parameters of F, G, H
C       Q:     System noise covarance matrix
C       XMEAN: Initial value of each component
C       XVAR:  Initial variance of each component
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  A(MAXM,NC), B(MAXM,NC), C(MAXM,NC), Q(NC,NC)
cc      DIMENSION  M(NC), TAU2(NC), AR(NC), MTYPE(NC)
cxx      DIMENSION  M(NC), TAU2(4), AR(M3), MTYPE(NC)
cc      DIMENSION  Y(N), XMEAN(NC), XVAR(MAXM,NC), DUM(1), COV(0:10)
cxx      DIMENSION  Y(N), XMEAN(NC), XVAR(MAXM,NC), DUM(1), COV(0:M3)
cc      COMMON  /C92827/  M1, M2, M3, M4, MPER
C
      INTEGER :: M1, M2, M3, M4, MPER, M(NC), L, NC, MTYPE(NC), MAXM,
     1           MM, N, ier
      REAL(8) :: TAU2(4), AR(M3), Y(N), A(MAXM,NC), B(MAXM,NC),
     1           C(MAXM,NC), Q(NC,NC), XMEAN(NC), XVAR(MAXM,NC)
      REAL(8) :: DUM(1), COV(0:M3), YMEAN, YVAR
C
      ier = 0
      L = 0
      CALL  MOMENT( Y,N/4,YMEAN,YVAR )
      YVAR = 1.0D4
C
cxx      DO 10  J=1,NC
cxx      DO 10  I=1,MAXM
cxx      XVAR(I,J) = 0.0D0
cxx      C(I,J) = 0.0D0
cxx      A(I,J) = 0.0D0
cxx   10 B(I,J) = 0.0D0
cxx      DO 20  I=1,NC
cxx      DO 20  J=1,NC
cxx   20 Q(I,J) = 0.0D0
      XVAR(1:MAXM,1:NC) = 0.0D0
      C(1:MAXM,1:NC) = 0.0D0
      A(1:MAXM,1:NC) = 0.0D0
      B(1:MAXM,1:NC) = 0.0D0
      Q(1:NC,1:NC) = 0.0D0
C
C  ...  TREND MODEL  ...
C
      IF( M1.GT.0 )  THEN
        L = 1
        M(L) = M1
        MTYPE(L) = 0
        IF( M1.EQ.1 )  THEN
          A(1,1) = 1.0D0
        END IF
        IF( M1.EQ.2 )  THEN
          A(1,1) =  2.0D0
          A(2,1) = -1.0D0
        END IF
        IF( M1.EQ.3 )  THEN
          A(1,1) =  3.0D0
          A(2,1) = -3.0D0
          A(3,1) =  1.0D0
        END IF
        B(1,1) = 1.0D0
        C(1,1) = 1.0D0
        Q(1,1) = TAU2(1)
        XMEAN(1) = YMEAN
        XVAR(1,1) = YVAR
      END IF
C
C  ...  SEASONAL MODEL  ...
C
      IF( M2.GT.0 )  THEN
        L = L+1
        M(L) = M2*(MPER-1)
        MTYPE(L) = 0
        IF( M2.NE.2 )  THEN
          DO 30  I=1,M(L)
cxx   30     A(I,L) = -1.0D0
          A(I,L) = -1.0D0
   30     CONTINUE
        END IF
        IF( M2.EQ.2 )  THEN
          DO 40 I=1,MPER-1
          A(I,L) = -(I+1)
cxx   40     A(M(L)-I+1,L) = -I
          A(M(L)-I+1,L) = -I
   40     CONTINUE
        END IF
        B(1,L) =  1.0D0
        C(1,L) =  1.0D0
        Q(L,L) =  TAU2(L)
        XMEAN(L) = 0.0D0
        XVAR(1,L) = YVAR
      END IF
C
C  ...  AR MODEL  ...
C
      IF( M3.GT.0 )  THEN
        L = L+1
        M(L) = M3
        MTYPE(L) = 0
        DO 50  I=1,M3
cxx   50   A(I,L) = AR(I)
        A(I,L) = AR(I)
   50   CONTINUE
        B(1,L) = 1.0D0
        C(1,L) = 1.0D0
        Q(L,L) = TAU2(L)
cc        CALL  ARMCOV( M3,0,AR,DUM,1.0D0,M3,COV )
        CALL  ARMCOV( M3,0,AR,DUM,1.0D0,M3,COV,M3,IER )
        if( ier.ne.0 ) return
        XMEAN(L) = 0.0D0
        DO 60 I=1,M3
cxx   60   XVAR(I,L) = COV(I-1)
        XVAR(I,L) = COV(I-1)
   60   CONTINUE
      END IF
C
C  ...  TRADING DAY EFFECT MODEL  ...
C
      IF( M4.GT.0 )  THEN
      DO 70  I=L+1,L+6
      M(I) = 1
      MTYPE(I) = I-L
      XMEAN(I) = 0.0D0
      XVAR(1,I) = YVAR
      A(1,I) = 1.0D0
      B(1,I) = 1.0D0
      C(1,I) = 1.0D0
cxx   70 Q(I,I) = 0.0D0
      Q(I,I) = 0.0D0
   70 CONTINUE
      L = L + 6
      END IF
C
C  ...  STATE DIMENSION  ...
C
      MSUM = 0
      DO 80  I=1,L
cxx   80 MSUM = MSUM + M(I)
      MSUM = MSUM + M(I)
   80 CONTINUE
      MM = MSUM
C
      RETURN
      E N D
cc      SUBROUTINE  FFSEAS( K,AA,FF,IFG )
      SUBROUTINE  FFSEAS( Y,N,M1,M2,M3,M4,NPER,REG,K,AA,
     * OUTMIN,OUTMAX,ALIMIT,M,XMEAN,XVAR,NS,NFE,NPE,NMAX,MJ,MAXM,NC,
     * FF,IFG,ier )
C
C  ... Function for Maximum likelihood estimation (seasonal model) ...
C
C     Inputs:
C        K:      Number of parameters
C        AA(I):  Parameter vector
C     Outputs:
C        FF:     -(Log likelihood)
C        IFG:    =1 if some conditions are violated
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  M(10), TAU2(10), PAR(10)
cc      DIMENSION  XMEAN(10), XVAR(22,10), AA(20), AR(10)
cc      DIMENSION  A(22,9), B(22,9), C(22,9), Q(9,9)
cc      DIMENSION  XPS(25,300), XFS(25,300)
cc      DIMENSION  VPS(25,25,300), VFS(25,25,300)
cc      DIMENSION  XF(40), VF(40,40), MTYPE(10)
cc      DIMENSION  Y(NMAX), REG(N,7)
cxx      DIMENSION  M(NC), TAU2(4), PAR(M3)
cxx      DIMENSION  XMEAN(NC), XVAR(MAXM,NC), AA(K), AR(M3)
cxx      DIMENSION  A(MAXM,NC), B(MAXM,NC), C(MAXM,NC), Q(NC,NC)
cxx      DIMENSION  XPS(MJ,NMAX), XFS(MJ,NMAX)
cxx      DIMENSION  VPS(MJ,MJ,NMAX), VFS(MJ,MJ,NMAX)
cxx      DIMENSION  XF(MJ), VF(MJ,MJ), MTYPE(NC)
cxx      DIMENSION  Y(N), REG(NMAX,7)
C
      INTEGER :: N, M1, M2, M3, M4, NPER, K, NS, NFE, NPE, NMAX,
     1           MJ, MAXM, NC, IFG, ier, M(NC)
      REAL(8) :: Y(N), REG(NMAX,7), AA(K), OUTMIN, OUTMAX, ALIMIT,
     1           XMEAN(NC), XVAR(MAXM,NC), FF
      INTEGER :: MTYPE(NC)
      REAL(8) :: TAU2(4), PAR(M3), AR(M3), A(MAXM,NC), B(MAXM,NC),
     1           C(MAXM,NC), Q(NC,NC), XPS(MJ,NMAX), XFS(MJ,NMAX),
     2           VPS(MJ,MJ,NMAX), VFS(MJ,MJ,NMAX), XF(MJ), VF(MJ,MJ),
     3           SIG2, OVAR
C
cc      COMMON  /C92827/  M1, M2, M3, M4, NPER
cc      COMMON  /C92826/   Y(160), REG(160,7)
cc      COMMON  /C92825/   OUTMIN, OUTMAX
cc      COMMON  /C92907/  ALIMIT
cc      COMMON  /CMMODL/  M, XMEAN, XVAR, N, NS, NFE, NPE
cc      NMAX= 300
cc      MJ  = 25
cc      MAXM= 22
cc      NC  = 9
      NP  = ID(M1) + ID(M2) + ID(M3)
C
      DO 10 I=1,K
cxx   10 IF( DABS(AA(I)).GT.20.0D0 )  GO TO 100
      IF( DABS(AA(I)).GT.20.0D0 )  GO TO 100
   10 CONTINUE
      DO 20 I=1,NP
cxx   20 TAU2(I) = DEXP( AA(I) )/(1.0D0 + DEXP( AA(I) ))
      TAU2(I) = DEXP( AA(I) )/(1.0D0 + DEXP( AA(I) ))
   20 CONTINUE
      DO 30 I=1,M3
cxx   30 PAR(I)  = ALIMIT*(DEXP(AA(NP+I))-1.0D0)/(DEXP(AA(NP+I))+1.0D0)
      PAR(I)  = ALIMIT*(DEXP(AA(NP+I))-1.0D0)/(DEXP(AA(NP+I))+1.0D0)
   30 CONTINUE
      CALL  ARCOEF( PAR,M3,AR )
      SIG2 = 1.0D0
      IFG = 0
C
cc      CALL  SETABC( M,L,NC,MTYPE,TAU2,MAXM,MM,AR,Y,N,A,B,C,Q,
cc     *              XMEAN,XVAR )
      CALL  SETABC1( M1,M2,M3,M4,NPER,M,L,NC,MTYPE,TAU2,MAXM,MM,AR,
     *                     Y,N,A,B,C,Q,XMEAN,XVAR,ier )
      if( ier.ne.0 ) return
      CALL  ISTAT1( L,M,MJ,MAXM,A,XMEAN,XVAR,XF,VF )
cc      CALL  FILTR1( Y,XF,VF,A,B,C,Q,SIG2,REG,MTYPE,NMAX,M,MAXM,NC,L,
      CALL  FILTR1( Y,XF,VF,A,B,C,Q,SIG2,REG,MTYPE,NMAX,M,MAXM,NC,L,
     *      NS,NFE,NPE,MJ,N,OUTMIN,OUTMAX,VFS,VPS,XFS,XPS,FF,OVAR )
cc      WRITE(6,*) FF
      FF = -FF
      RETURN
C
  100 IFG = 1
      FF = 1.0D20
      RETURN
      E N D
      SUBROUTINE  TRADE( JYEAR,MONTH,N,MJ,TDAY )
C
C  ...  This subroutine computes the number of days of the week
C       in each month,
C
C     Inputs:
C        JYEAR:  Starting year of the data set
C        MONTH:  Starting month of the data set
C        N:      Data length
C        MJ:     Adjustable dimension of TDAY
C     Output:
C        TDAY(I,J):  Number of J-th days of the week in I-th month
C     Written by F.K.  Nov.1981
C
cxx      REAL*8  TDAY
cxx      DIMENSION  TDAY(MJ,7), IX(12)
C
      INTEGER :: JYEAR, MONTH, N, MJ
      REAL(8) :: TDAY(MJ,7)
      INTEGER :: IX(12)
C
      DATA   IX  /3,0,3,2,3,2,3,3,2,3,2,3/
C
      JS = JYEAR - 1900
      I0 = MOD( JS+(JS-1)/4,7 ) + 1
         I2 = 0
      JJ = 2 - MONTH
      II = 0
    5 II = II + 1
      I1 = II + JS - 1
      IX(2) = 0
      IF( MOD(I1,4) .EQ. 0 )  IX(2) = 1
      IF( MOD(I1,100).EQ.0 )  IX(2) = 0
      IF( MOD(I1,400).EQ.0 )  IX(2) = 1
cxx      DO 30  J=1,12
      DO 31  J=1,12
      DO 10  I=1,7
cc   10 IF( JJ.GT.0 )  TDAY(JJ,I) = 4
cxx   10 IF( JJ.GT.0 .AND. JJ.LE.MJ )  TDAY(JJ,I) = 4
      IF( JJ.GT.0 .AND. JJ.LE.MJ )  TDAY(JJ,I) = 4
   10 CONTINUE
C
      IE = IX(J)
      IF( IE .EQ. 0 )  GO TO 30
      DO 20  I=1,IE
      I2 = I0 + I
      IF( I2 .GT. 7 ) I2 = I2 - 7
cc   20 IF( JJ.GT.0 )  TDAY(JJ,I2) = 5
cxx   20 IF( JJ.GT.0 .AND. JJ.LE.MJ )  TDAY(JJ,I2) = 5
      IF( JJ.GT.0 .AND. JJ.LE.MJ )  TDAY(JJ,I2) = 5
   20 CONTINUE
      I0 = I2
   30 JJ = JJ + 1
   31 CONTINUE
      IF( JJ .GT.N ) RETURN
      GO TO 5
C
      E N D
cc      SUBROUTINE  DAVIDN( FUNCT,X,N,NDIF )
      SUBROUTINE  DAVIDN1( FUNCT,X,N,NDIF,
     * YY,NN,M1,M2,M3,M4,NPER,REG,OUTMIN,OUTMAX,ALIMIT,ISW,IDIF,
     * M,XMEAN,XVAR,NS,NFE,NPE,NMAX,MJ,MAXM,NC,IER )
C
C  ...  6/20/83, 12/19/92
C
cxx      IMPLICIT  REAL*8( A-H,O-Z )
cc      DIMENSION  X(40), DX(40), G(40), G0(40), Y(40)
cc      DIMENSION  H(40,40), WRK(40), S(40)
cxx      DIMENSION  X(N), DX(N), G(N), G0(N), Y(N)
cxx      DIMENSION  H(N,N), WRK(N), S(N)
cxx      DIMENSION  YY(NN), REG(NMAX,7)
cxx      DIMENSION  M(NC), XMEAN(NC), XVAR(MAXM,NC)
C
      INTEGER :: N, NDIF, NN, M1, M2, M3, M4, NPER, ISW, IDIF, NS, NFE,
     1           NPE, NMAX, MJ, MAXM, NC, IER, M(NC)
      REAL(8) :: X(N), YY(NN), REG(NMAX,7), OUTMIN, OUTMAX, ALIMIT,
     1           XMEAN(NC), XVAR(MAXM,NC)
      REAL(8) :: DX(N), G(N), G0(N), Y(N), H(N,N), WRK(N), S(N), TAU2,
     1           EPS1, EPS2, RAMDA, CONST1, XM, SUM, S1, S2, SS, STEM,
     2           ED, XMB
C
cc      COMMON     / CCC /  ISW, IPR, ISMT, IDIF
cc      COMMON     / DDD /  XM , AIC , SD
      EXTERNAL  FUNCT
      DATA        TAU2  /          1.0D-6  /
      DATA  EPS1, EPS2  / 1.0D-6 , 1.0D-6  /
      RAMDA  = 0.5D0
      CONST1 = 1.0D-30
      IPR = 0
c-----
      IDIF = NDIF
c-----
C
C          INITIAL ESTIMATE OF INVERSE OF HESSIAN
C
      ICOUNT = 0
cxx 1000 DO 20  I=1,N
cxx      DO 10  J=1,N
cxx   10 H(I,J) = 0.0D0
cxx      S(I)   = 0.0D0
cxx      DX(I)  = 0.0D0
cxx   20 H(I,I) = 1.0D0
 1000 CONTINUE
      H(1:N,1:N) = 0.0D0
      S(1:N)   = 0.0D0
      DX(1:N)  = 0.0D0
      DO 20  I=1,N
      H(I,I) = 1.0D0
   20 CONTINUE

      ISW = 0
C
cc      IF( NDIF.EQ.0 )  CALL  FUNCT( N,X,XM,G,IG )
cc      IF( NDIF.GE.1 )  CALL  FUNCND( FUNCT,N,X,XM,G,IG )
      IF( NDIF.EQ.0 )  CALL  FUNCT( YY,NN,M1,M2,M3,M4,NPER,REG,N,X,
     * OUTMIN,OUTMAX,ALIMIT,M,XMEAN,XVAR,NS,NFE,NPE,NMAX,MJ,
     * MAXM,NC,XM,IG,ier )
      IF( NDIF.GE.1 )  CALL  FUNCND1( FUNCT,N,X,XM,G,IG,IDIF,
     * YY,NN,M1,M2,M3,M4,NPER,REG,OUTMIN,OUTMAX,ALIMIT,M,XMEAN,XVAR,
     * NS,NFE,NPE,ISW,NMAX,MJ,MAXM,NC,ier )
      if( ier.ne.0 ) return
C
cc      IF( IPR .GE. 2 )   WRITE( 6,640 )   XM, SD, AIC
cc      WRITE(6,650) (X(I),I=1,N), XM, (G(I),I=1,N)
C
      ICC = 0
C      ITERATION
 2000 CONTINUE
      ICC = ICC + 1
      IF( ICC .EQ. 1 )  GO TO 120
C
      DO 40  I=1,N
cxx   40 Y(I) = G(I) - G0(I)
      Y(I) = G(I) - G0(I)
   40 CONTINUE
      DO 60  I=1,N
      SUM = 0.0D0
      DO 50  J=1,N
cxx   50 SUM = SUM + Y(J)*H(I,J)
      SUM = SUM + Y(J)*H(I,J)
   50 CONTINUE
cxx   60 WRK(I) = SUM
      WRK(I) = SUM
   60 CONTINUE
      S1 = 0.0D0
      S2 = 0.0D0
      DO 70  I=1,N
      S1 = S1 + WRK(I)*Y(I)
cxx   70 S2 = S2 + DX(I) *Y(I)
      S2 = S2 + DX(I) *Y(I)
   70 CONTINUE
      IF( S1.LE.CONST1 .OR. S2.LE.CONST1 )  GO TO 900
C     IF( S1 .LE. S2 )   GO TO 100
C
C          UPDATE THE INVERSE OF HESSIAN MATRIX
C
C               ---  BROYDEN-FLETCHER-GOLDFARB-SHANNO TYPE CORRECTION  -
C
cxx  100 CONTINUE
      STEM = S1 / S2 + 1.0D0
cxx      DO 110  I=1,N
      DO 111  I=1,N
      DO 110  J=I,N
      H(I,J) = H(I,J)-(DX(I)*WRK(J)+WRK(I)*DX(J)-DX(I)*DX(J)*STEM)/S2
cxx  110 H(J,I) = H(I,J)
      H(J,I) = H(I,J)
  110 CONTINUE
  111 CONTINUE
C
  120 CONTINUE
      SS = 0.0D0
      DO 140  I=1,N
      SUM = 0.0D0
      DO 130  J=1,N
cxx  130 SUM = SUM + H(I,J)*G(J)
      SUM = SUM + H(I,J)*G(J)
  130 CONTINUE
      SS  = SS + SUM**2
cxx  140 S(I) = -SUM
      S(I) = -SUM
  140 CONTINUE
C
      S1 = 0.0D0
      S2 = 0.0D0
      DO 150  I=1,N
      S1 = S1 + S(I)*G(I)
cxx  150 S2 = S2 + G(I)*G(I)
      S2 = S2 + G(I)*G(I)
  150 CONTINUE
C     DS2 = DSQRT(S2)
C     GTEM = DABS(S1)/DS2
C     IF( GTEM.LE.TAU1 .AND. DS2.LE.TAU2 )  GO TO 900
      IF( S1.LT.0.0D0 )  GO TO 200
      DO 170  I=1,N
      DO 160  J=1,N
cxx  160 H(I,J) = 0.0D0
      H(I,J) = 0.0D0
  160 CONTINUE
      H(I,I) = 1.0D0
cxx  170 S(I) = -S(I)
      S(I) = -S(I)
  170 CONTINUE
  200 CONTINUE
C
      ED = XM
C
C          LINEAR  SEARCH
C
cc      CALL  LINEAR( FUNCT,X,S,RAMDA,ED,N,IG )
      CALL  LINEAR1( FUNCT,X,S,RAMDA,ED,N,IG,
     * YY,NN,M1,M2,M3,M4,NPER,REG,OUTMIN,OUTMAX,ALIMIT,M,XMEAN,XVAR,
     * NS,NFE,NPE,ISW,NMAX,MJ,MAXM,NC,ier )
      if( ier.ne.0 ) return
C
cc      IF( IPR .GE. 2 )  WRITE( 6,630 )  RAMDA, ED, SD, AIC
C
      S1 = 0.0D0
      DO 210  I=1,N
      DX(I) = S(I)*RAMDA
      S1 = S1 + DX(I)*DX(I)
      G0(I) = G(I)
cxx  210 X(I) = X(I) + DX(I)
      X(I) = X(I) + DX(I)
  210 CONTINUE
      XMB = XM
      ISW = 0
C
cc      IF( NDIF.EQ.0 )  CALL  FUNCT( N,X,XM,G,IG )
cc      IF( NDIF.GE.1 )  CALL  FUNCND( FUNCT,N,X,XM,G,IG )
      IF( NDIF.EQ.0 )  CALL  FUNCT( YY,NN,M1,M2,M3,M4,NPER,REG,N,X, 
     * OUTMIN,OUTMAX,ALIMIT,M,XMEAN,XVAR,NS,NFE,NPE,NMAX,MJ,
     * MAXM,NC,XM,IG,ier )
      IF( NDIF.GE.1 )  CALL  FUNCND1( FUNCT,N,X,XM,G,IG,IDIF,
     * YY,NN,M1,M2,M3,M4,NPER,REG,OUTMIN,OUTMAX,ALIMIT,M,XMEAN,XVAR,
     * NS,NFE,NPE,ISW,NMAX,MJ,MAXM,NC,ier )
      if( ier.ne.0 ) return
cc      WRITE(6,650) (X(I),I=1,N), XM, (G(I),I=1,N)
C
      S2 = 0.D0
      DO 220  I=1,N
cxx  220 S2 = S2 + G(I)*G(I)
      S2 = S2 + G(I)*G(I)
  220 CONTINUE
      IF( DSQRT(S2) .LT. TAU2 )  GO TO 900
      IF( XMB/XM-1.D0.LT.EPS1 .AND. DSQRT(S1).LT.EPS2 )  GO TO 900
C     IF( ICC .GE. 5 )  GO TO 900
      GO TO 2000
  900 CONTINUE
      IF( IPR .LE. 0 )  RETURN
cc      WRITE( 6,600 )
cc      WRITE( 6,610 )  (X(I),I=1,N)
cc      WRITE( 6,620 )
cc      WRITE( 6,610 )  (G(I),I=1,N)
C
      ICOUNT  = ICOUNT + 1
      S2 = 0.D0
      DO 910  I=1,N
cxx  910 S2 = S2 + G(I)*G(I)
      S2 = S2 + G(I)*G(I)
  910 CONTINUE
      IF( S2.GT.1.0D0.AND.ICOUNT.LE.1 )   GO TO 1000
      RETURN
cxx  600 FORMAT( 1H0,'-----  X  -----' )
cxx  610 FORMAT( 1H ,10D13.5 )
cxx  620 FORMAT( 1H0,'***  GRADIENT  ***' )
cxx  630 FORMAT( 1H ,'LAMBDA =',D15.7,3X,'(-1)LOG LIKELIHOOD =',D23.15,
cxx     *        3X,'SD =',D22.15,5X,'AIC =',D23.15 )
cxx  640 FORMAT( 1H ,26X,'(-1)LOG-LIKELIHOOD =',D23.15,3X,'SD =',D22.15,
cxx     *        5X,'AIC =',D23.15 )
cxx  650 FORMAT( 10X,5F12.7 )
      E N D
cc      SUBROUTINE  FUNCND( FUNCT,M,A,F,G,IFG )
      SUBROUTINE  FUNCND1( FUNCT,M,A,F,G,IFG,IDIF,Y,N,M1,M2,M3,M4,
     * NPER,REG,OUTMIN,OUTMAX,ALIMIT,MM,XMEAN,XVAR,NS,NFE,NPE,ISW,
     * NMAX,MJ,MAXM,NC,ier )
C
C  ...  FUNCTION EVALUATION AND NUMERICAL DIFFERENCING  ...
C
cxx      IMPLICIT   REAL*8( A-H,O-Z )
cc      DIMENSION  A(M) , G(M) , B(20),GD(5)
cxx      DIMENSION  A(M) , G(M) , B(M)
cxx      DIMENSION  Y(N), REG(NMAX,7), MM(NC), XMEAN(NC), XVAR(MAXM,NC)
C
      INTEGER :: M, IFG, IDIF, N, M1, M2, M3, M4, NPER, NS, NFE, NPE,
     1           ISW, NMAX, MJ, MAXM, NC, ier, MM(NC)
      REAL(8) :: A(M), F, G(M), Y(N), REG(NMAX,7), OUTMIN, OUTMAX,
     1           ALIMIT, XMEAN(NC), XVAR(MAXM,NC)
      REAL(8) :: B(M), CONST, FF, FB
C
cc      COMMON  / CCC /  ISW , IPR, ISMT, IDIF
cc      COMMON  /CMFUNC/ DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(20),GC(20)
C     DATA       ICNT /0/
      CONST = 0.00001D0
C
cc      CALL  FUNCT( M,A,F,GD,IFG )
      CALL FUNCT( Y,N,M1,M2,M3,M4,NPER,REG,M,A,OUTMIN,OUTMAX,
     * ALIMIT,MM,XMEAN,XVAR,NS,NFE,NPE,NMAX,MJ,MAXM,NC,F,IFG,ier )
      if( ier.ne.0 ) return
      FB = F
      IF( ISW .GE. 1 )   RETURN
C
C     WRITE( 6,600 )   (A(I),I=1,M)
      DO 10  I=1,M
cxx   10 B(I) = A(I)
      B(I) = A(I)
   10 CONTINUE
C
      DO 30  II=1,M
      B(II) = A(II) + CONST
cc      CALL  FUNCT( M,B,FF,GD,IFG )
      CALL FUNCT( Y,N,M1,M2,M3,M4,NPER,REG,M,B,OUTMIN,OUTMAX,
     * ALIMIT,MM,XMEAN,XVAR,NS,NFE,NPE,NMAX,MJ,MAXM,NC,FF,IFG,ier )
      if( ier.ne.0 ) return
      IF( IDIF .EQ. 1 )  GO TO 20
      B(II) = A(II) - CONST
cc      CALL  FUNCT( M,B,FB,GD,IFG )
      CALL FUNCT( Y,N,M1,M2,M3,M4,NPER,REG,M,B,OUTMIN,OUTMAX,
     * ALIMIT,MM,XMEAN,XVAR,NS,NFE,NPE,NMAX,MJ,MAXM,NC,FB,IFG,ier )
      if( ier.ne.0 ) return
   20 G(II) = (FF-FB)/(CONST*IDIF)
      IF( G(II) .GT. 1.0D20 )  G(II) = (F-FB)/CONST
      IF( G(II) .LT.-1.0D20 )  G(II) = (FF-F)/CONST
      IF( FB.GT.F .AND. FF.GT.F )  G(II) = 0.0D0
cxx   30 B(II) = A(II)
      B(II) = A(II)
   30 CONTINUE
C
      RETURN
C
cxx  600 FORMAT( 3X,'---  PARAMETER  ---',(/,3X,5D13.5) )
cxx  610 FORMAT( 3X,'---  GRADIENT  ---',(/,3X,5D13.5) )
      E N D
cc      SUBROUTINE  LINEAR( FUNCT,X,H,RAM,EE,K,IG )
      SUBROUTINE  LINEAR1( FUNCT,X,H,RAM,EE,K,IG, Y,N,M1,M2,M3,M4,
     * NPER,REG,OUTMIN,OUTMAX,ALIMIT,M,XMEAN,XVAR,NS,NFE,NPE,ISW,
     * NMAX,MJ,MAXM,NC,ier )
C
C  ...  LINEAR SEARCH  ...
C
cxx      IMPLICIT  REAL*8( A-H,O-Z )
cc      INTEGER  RETURN,SUB
cc      DIMENSION  X(1), H(1), X1(40)
cc      DIMENSION  G(40)
cc      COMMON     / CCC /  ISW , IPR, ISMT, IDIF
cxx      DIMENSION  X(K), H(K), X1(K)
cxx      DIMENSION  Y(N), REG(NMAX,7)
cxx      DIMENSION  M(NC), XMEAN(NC), XVAR(MAXM,NC)
C
      INTEGER :: K, IG, N, M1, M2, M3, M4, NPER, NS, NFE, NPE, ISW,
     1           NMAX, MJ, MAXM, NC, ier, M(NC) 
      REAL(8) :: X(K), H(K), RAM, EE, Y(N), REG(NMAX,7), OUTMIN, OUTMAX,
     1           ALIMIT, XMEAN(NC), XVAR(MAXM,NC)
      REAL(8) :: X1(K), CONST2, HNORM, RAM1, RAM2, RAM3, E1, E2, E3,
     1           A1, A2, A3, B1, B2
C
C     IPR = 7
C     C1 = 1.0D-7
C     C2 = 1.0D-30
C
      ISW = 1
      IG  = 0
      RAM = 0.1D0
      CONST2 = 1.0D-60
   14 DO 15 I=1,K
cxx   15 IF( DABS(H(I)) .GT. 1.0D10 )  GO TO 16
      IF( DABS(H(I)) .GT. 1.0D10 )  GO TO 16
   15 CONTINUE
      GO TO 18
   16 DO 17 I=1,K
cxx   17 H(I) = H(I)*1.0D-10
      H(I) = H(I)*1.0D-10
   17 CONTINUE
      GO TO 14
   18 CONTINUE
      HNORM = 0.D0
      DO 10  I=1,K
cxx   10 HNORM = HNORM + H(I)**2
      HNORM = HNORM + H(I)**2
   10 CONTINUE
      HNORM = DSQRT( HNORM )
C
      RAM2 = RAM
      E1 =EE
      RAM1 = 0.D0
C
      DO 20  I=1,K
cxx   20 X1(I) = X(I) + RAM2*H(I)
      X1(I) = X(I) + RAM2*H(I)
   20 CONTINUE
cc      CALL  FUNCT( K,X1,E2,G,IG )
      CALL  FUNCT( Y,N,M1,M2,M3,M4,NPER,REG,K,X1,OUTMIN,OUTMAX,
     * ALIMIT,M,XMEAN,XVAR,NS,NFE,NPE,NMAX,MJ,MAXM,NC,E2,IG,ier )
      if( ier.ne.0 ) return
cc      IF(IPR.GE.7)  WRITE(6,2)  RAM2,E2
C
      IF( IG .EQ. 1  )  GO TO 50
      IF( E2 .GT. E1 )  GO TO 50
   30 RAM3 = RAM2*4.D0
      DO 40  I=1,K
cxx   40 X1(I) = X(I) + RAM3*H(I)
      X1(I) = X(I) + RAM3*H(I)
   40 CONTINUE
cc      CALL  FUNCT( K,X1,E3,G,IG )
      CALL  FUNCT( Y,N,M1,M2,M3,M4,NPER,REG,K,X1,OUTMIN,OUTMAX,
     * ALIMIT,M,XMEAN,XVAR,NS,NFE,NPE,NMAX,MJ,MAXM,NC,E3,IG,ier )
      if( ier.ne.0 ) return
      IF( IG.EQ.1 )  GO TO  500
cc      IF( IPR.GE.7 )  WRITE(6,3)  RAM3,E3
      IF( E3 .GT. E2 )  GO TO 70
      IF(RAM3.GT.1.0D10 .AND. E3.LT.E1)  GO TO 45
      IF(RAM3.GT.1.0D10 .AND. E3.GE.E1)  GO TO 46
      RAM1 = RAM2
      RAM2 = RAM3
      E1 = E2
      E2 = E3
      GO TO 30
C
   45 RAM = RAM3
      EE = E3
      RETURN
C
   46 RAM = 0.0D0
      RETURN
C
   50 RAM3 = RAM2

      E3 = E2
      RAM2 = RAM3*0.1D0
      IF( RAM2*HNORM .LT. CONST2 )  GO TO  400
      DO 60  I=1,K
cxx   60 X1(I) = X(I) + RAM2*H(I)
      X1(I) = X(I) + RAM2*H(I)
   60 CONTINUE
cc      CALL  FUNCT( K,X1,E2,G,IG )
      CALL  FUNCT( Y,N,M1,M2,M3,M4,NPER,REG,K,X1,OUTMIN,OUTMAX,
     * ALIMIT,M,XMEAN,XVAR,NS,NFE,NPE,NMAX,MJ,MAXM,NC,E2,IG,ier )
      if( ier.ne.0 ) return
cc      IF(IPR.GE.7)  WRITE(6,4)  RAM2,E2
      IF( E2.GT.E1 )  GO TO 50
C
cc   70 ASSIGN 80 TO RETURN
   70 CONTINUE
      IRET = 80
      GO TO 200
C
   80 DO 90  I=1,K
cxx   90 X1(I) = X(I) + RAM*H(I)
      X1(I) = X(I) + RAM*H(I)
   90 CONTINUE
cc      CALL  FUNCT( K,X1,EE,G,IG )
      CALL  FUNCT( Y,N,M1,M2,M3,M4,NPER,REG,K,X1,OUTMIN,OUTMAX,
     * ALIMIT,M,XMEAN,XVAR,NS,NFE,NPE,NMAX,MJ,MAXM,NC,EE,IG,ier )
      if( ier.ne.0 ) return
cc      IF(IPR.GE.7)  WRITE(6,5)  RAM,EE
C
      IFG = 0
cc      ASSIGN  300 TO  SUB
cc      ASSIGN 200 TO SUB
cc   95 ASSIGN 130 TO RETURN
      ISUB = 200
   95 CONTINUE
      IRET = 130
      IF( RAM .GT. RAM2 )  GO TO 110
      IF( EE .GE. E2 )  GO TO 100
      RAM3 = RAM2
      RAM2 = RAM
      E3 =E2
      E2 =EE
cc      GO TO  SUB,( 200,300 )
      IF( ISUB.EQ.200 ) GO TO 200
      IF( ISUB.EQ.300 ) GO TO 300
C
  100 RAM1 = RAM
      E1 = EE
cc      GO TO  SUB,( 200,300 )
      IF( ISUB.EQ.200 ) GO TO 200
      IF( ISUB.EQ.300 ) GO TO 300
C
  110 IF( EE .LE. E2 )  GO TO 120
      RAM3 = RAM
      E3 = EE
cc      GO TO  SUB,( 200,300 )
      IF( ISUB.EQ.200 ) GO TO 200
      IF( ISUB.EQ.300 ) GO TO 300
C
  120 RAM1 = RAM2
      RAM2 = RAM
      E1 = E2
      E2 = EE
cc      GO TO  SUB,( 200,300 )
      IF( ISUB.EQ.200 ) GO TO 200
      IF( ISUB.EQ.300 ) GO TO 300
C
  130 DO 140  I=1,K
cxx  140 X1(I) = X(I) + RAM*H(I)
      X1(I) = X(I) + RAM*H(I)
  140 CONTINUE
cc      CALL  FUNCT( K,X1,EE,G,IG )
      CALL  FUNCT( Y,N,M1,M2,M3,M4,NPER,REG,K,X1,OUTMIN,OUTMAX,
     * ALIMIT,M,XMEAN,XVAR,NS,NFE,NPE,NMAX,MJ,MAXM,NC,EE,IG,ier )
      if( ier.ne.0 ) return
cc      IF( IPR.GE.7 )  WRITE(6,6)  RAM,EE
cc      ASSIGN 200 TO SUB
      ISUB = 200
      IFG = IFG+1
cxx  600 FORMAT( 1H ,6D20.13 )
C     IF( HNORM*(RAM3-RAM1) .GT. 1.0D-12 )  GO TO 95
      IF( RAM3-RAM1 .GT. RAM2*1.0D-1 )  GO TO 95
C     IF( E1-EE .GT. 1.0D-15 )  GO TO 95
C     IF( E3-EE .GT. 1.0D-15 )  GO TO 95
C     IF(DX.LE.C1*XN+C2 .AND. DF.LE.C1*DABS(EE)+C2)  RETURN
C     IF( IFG .LE. 2 )  GO TO 95
C
      IF( E2 .LT. EE )  RAM = RAM2
      RETURN
C
C      -------  INTERNAL SUBROUTINE SUB1  -------
  200 IF( RAM3-RAM2 .GT. 5.0D0*(RAM2-RAM1) )  GO TO 202
      IF( RAM2-RAM1 .GT. 5.0D0*(RAM3-RAM2) )  GO TO 204
      A1 = (RAM3-RAM2)*E1
      A2 = (RAM1-RAM3)*E2
      A3 = (RAM2-RAM1)*E3
      B2 = (A1+A2+A3)*2.D0
      B1 = A1*(RAM3+RAM2) + A2*(RAM1+RAM3) + A3*(RAM2+RAM1)
      IF( B2 .EQ. 0.D0 )  GO TO 210
      RAM = B1 /B2
C
      IF( RAM .LE. RAM1 )  RAM = (RAM1 + RAM2)/2.0D0
      IF( RAM .GE. RAM3 )  RAM = (RAM2 + RAM3)/2.0D0
      IF( DABS(RAM-RAM2) .LE. 1.0D-15 )  RAM = (RAM2*4.D0 + RAM3)/5.0D0
cc      GO TO RETURN ,( 80,130 )
      IF( IRET.EQ.80 ) GO TO 80
      IF( IRET.EQ.130 ) GO TO 130
  202 RAM = (4.0D0*RAM2 + RAM3)/5.0D0
cc      GO TO RETURN, (80,130)
      IF( IRET.EQ.80 ) GO TO 80
      IF( IRET.EQ.130 ) GO TO 130
  204 RAM = (RAM1 + 4.0D0*RAM2)/5.0D0
cc      GO TO RETURN, (80,130)
      IF( IRET.EQ.80 ) GO TO 80
      IF( IRET.EQ.130 ) GO TO 130
C
  210 IG = 1
      RAM = RAM2
      RETURN
C
C      -------  INTERNAL SUBROUTINE SUB2  -------
C
  300 IF( RAM3-RAM2 .GT. RAM2-RAM1 )  GO TO 310
      RAM = (RAM1+RAM2)*0.5D0
cc      GO TO RETURN ,( 80,130 )
      IF( IRET.EQ.80 ) GO TO 80
      IF( IRET.EQ.130 ) GO TO 130
C
  310 RAM = (RAM2+RAM3)*0.5D0
cc      GO TO RETURN ,( 80,130 )
      IF( IRET.EQ.80 ) GO TO 80
      IF( IRET.EQ.130 ) GO TO 130
C
  400 RAM = 0.D0
      RETURN
C
  500 RAM = (RAM2+RAM3)*0.5D0
  510 DO 520  I=1,K
cxx  520 X1(I) = X(I) + RAM*H(I)
      X1(I) = X(I) + RAM*H(I)
  520 CONTINUE
cc      CALL  FUNCT( K,X1,E3,G,IG )
      CALL  FUNCT( Y,N,M1,M2,M3,M4,NPER,REG,K,X1,OUTMIN,OUTMAX,
     * ALIMIT,M,XMEAN,XVAR,NS,NFE,NPE,NMAX,MJ,MAXM,NC,E3,IG,ier )
      if( ier.ne.0 ) return
cc      IF( IPR.GE.7 )  WRITE(6,7)  RAM,E3
      IF( IG.EQ.1 )  GO TO 540
      IF( E3.GT.E2 )  GO TO 530
      RAM1 = RAM2
      RAM2 = RAM
      E1 = E2
      E2 = E3
      GO TO 500
C
  530 RAM3 = RAM
      GO TO 70
C
  540 RAM = (RAM2+RAM)*0.5D0
      GO TO 510
C
C ------------------------------------------------------------
cxx    1 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E1 =',D25.17 )
cxx    2 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E2 =',D25.17 )
cxx    3 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E3 =',D25.17 )
cxx    4 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E4 =',D25.17 )
cxx    5 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E5 =',D25.17 )
cxx    6 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E6 =',D25.17 )
cxx    7 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E7 =',D25.17 )
      E N D
