C     PROGRAM 10.1  ARMAFT
      SUBROUTINE ARMAFT( Y0,N,M,L,MLMAX,IPARAM,AR0,CMA0,SIG2,FLK,AIC,
     *                   AR,CMA,IER )
C
      INCLUDE 'TSSS.h'
C
C  ...  ARMA MODEL FITTING  ...
C
C     Inputs:
C        L:       AR Order (M <= 20)
C        M:       MA Order (M <= 20)
C        IPARAM:  =0    Use defalt initail values
C                 =1    Read intial values
C        AR(I):   AR coefficients (I=1,M3)
C     Parameters:
C        NMAX:    Adjustable dimension of Y
C        MJ:      Adjustable dimension of XF, VF, etc.
C     @TEST.FILTER2O    NOV.29,1990, SEP.02,1992
C
cc      PARAMETER( NMAX=1000,MJ=20 )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  AA(20), AR(20), PAR(20), CMA(20)
cxx      DIMENSION  AR0(M), CMA0(L)
cxx      DIMENSION  AA(M+L), AR(M), PAR(MLMAX), CMA(L)
cxx      DIMENSION  Y0(N), Y(N)
cc      COMMON  /C92825/  OUTMIN, OUTMAX
cc      COMMON  /C92826/  Y(NMAX)
cc      COMMON  /C92907/  ALIMIT
cc      COMMON  /C92908/  M, L, N
cc      COMMON  /C92909/  FLK, SIG2
cc      COMMON  / CCC /  ISW, IPR, ISMT, IDIF
C
      INTEGER N, M, L, MLMAX, IPARAM, IER 
      DOUBLE PRECISION Y0(N), AR0(M), CMA0(L), SIG2, FLK, AIC, AR(M),
     1                 CMA(L)
c local
      INTEGER I, IOPT, NSUM
      DOUBLE PRECISION Y(N), AA(M+L), PAR(MLMAX), ALIMIT, OUTMIN,
     1                 OUTMAX, YMEAN
      EXTERNAL  FFARMA
C
C  ...  Read Model Orders  ...
C
cc      READ( 5,* )  M, L, IPARAM
C
C  ...  Set Defalt Parameters  ...
C
cc      IPR  = 7
      ALIMIT = 0.95D0
cc      CALL  SPARA1( M,L,AR,CMA,OUTMIN,OUTMAX,IOPT )
      CALL  SPARA1( M,L,MLMAX,AR,CMA,OUTMIN,OUTMAX,IOPT )
      IF( IPARAM.EQ.1 )  THEN
         DO 5 I=1,M
cxx    5    AR(I) = AR0(I)
         AR(I) = AR0(I)
    5    CONTINUE
         DO 6 I=1,L
cxx    6    CMA(I) = CMA0(I)
         CMA(I) = CMA0(I)
    6    CONTINUE
cc         READ( 5,* )  (AR(I),I=1,M)
cc         READ( 5,* )  (CMA(I),I=1,L)
      END IF
C
C  ...  Read Time Series  ...
C
cc      CALL  READTS( 1,Y,N )
      DO 7 I = 1,N
cxx    7 Y(I) = Y0(I)
      Y(I) = Y0(I)
    7 CONTINUE
C
C  ...  Subtrac Mean Value  ...
C
cc      SUM = 0.0D0
cc      DO 10 I=1,N
cccxx   10 SUM = SUM + Y(I)
cc      SUM = SUM + Y(I)
cc   10 CONTINUE
cc      YMEAN = SUM/N
cxx      CALL  MEAN( Y,N,-1.0D30,1.0D30,SUM,YMEAN )
      CALL  MEAN( Y,N,-1.0D30,1.0D30,NSUM,YMEAN )
      DO 20 I=1,N
cxx   20 Y(I) = Y(I) - YMEAN
      Y(I) = Y(I) - YMEAN
   20 CONTINUE
C
cc      WRITE(6,*) M, L
cc      WRITE(6,*) (AR(I),I=1,M)
C
      CALL  PARCOR( AR,M,PAR )
      DO 30 I=1,M
cxx   30 AA(I) = DLOG( (ALIMIT+PAR(I))/(ALIMIT-PAR(I)) )
      AA(I) = DLOG( (ALIMIT+PAR(I))/(ALIMIT-PAR(I)) )
   30 CONTINUE
      CALL  PARCOR( CMA,L,PAR )
      DO 40 I=1,L
cxx   40 AA(M+I) = DLOG( (ALIMIT+PAR(I))/(ALIMIT-PAR(I)) )
      AA(M+I) = DLOG( (ALIMIT+PAR(I))/(ALIMIT-PAR(I)) )
   40 CONTINUE
C
C  ...  Maximum Likelihood Method  ...
C
cc      IF( IOPT.EQ.1 )  CALL  DAVIDN( FFARMA,AA,M+L,2 )
      ier = 0
      IF( IOPT.EQ.1 )  CALL  DAVIDN( FFARMA,AA,M+L,2,
     *    Y,N,M,L,MLMAX,OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
      if( ier.ne.0 ) return
C
      DO 50 I=1,M
cxx   50 PAR(I)  = ALIMIT*(DEXP(AA(I))-1.0D0)/(DEXP(AA(I))+1.0D0)
      PAR(I)  = ALIMIT*(DEXP(AA(I))-1.0D0)/(DEXP(AA(I))+1.0D0)
   50 CONTINUE
      CALL  ARCOEF( PAR,M,AR )
      DO 60 I=1,L
cxx   60 PAR(I)  = ALIMIT*(DEXP(AA(M+I))-1.0D0)/(DEXP(AA(M+I))+1.0D0)
      PAR(I)  = ALIMIT*(DEXP(AA(M+I))-1.0D0)/(DEXP(AA(M+I))+1.0D0)
   60 CONTINUE
      CALL  ARCOEF( PAR,L,CMA )
C
      AIC = -2*FLK + 2*(M+L+1)
C
C  ...  Print out the Maximum Likelihood Estimates  ...
C
cc      CALL  PRARMA( M,L,AR,CMA,SIG2,FLK,AIC )
C
cc      STOP
      RETURN
      E N D
C
C------------------------------
C  armafit, armafit2 commoon 
C------------------------------
C
cc      SUBROUTINE  SPARA1( M,L,AR,CMA,OUTMIN,OUTMAX,IOPT )
      SUBROUTINE  SPARA1( M,L,MLMAX,AR,CMA,OUTMIN,OUTMAX,IOPT )
C
C  ...  Set Default Parameters  ...
C
C     Inputs:
C       M:       AR order
C       L:       MA order
C     Outputs:
C       TAU2:    System noise variance
C       AR:      AR coefficients
C       OUTMIN:  Lower bound for outliers
C       OUTMAX:  Upper bound for outliers
C       IOPT:    (=1  MLE by numerical optimization)
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  AR(20), PAR(20), CMA(20)
cxx      DIMENSION  AR(M), PAR(MLMAX), CMA(L)
C
      INTEGER M, L, MLMAX, IOPT
      DOUBLE PRECISION AR(M), CMA(L), OUTMIN, OUTMAX
c local
      INTEGER I
      DOUBLE PRECISION PAR(MLMAX)
C
      DO 10 I=1,M
cxx   10 PAR(I) = -(-0.6D0)**I
      PAR(I) = -(-0.6D0)**I
   10 CONTINUE
      CALL  ARCOEF( PAR,M,AR )
      DO 20 I=1,L
cxx   20 PAR(I) = -(-0.6D0)**I
ccc      PAR(I) = -(-0.6D0)**I
      PAR(I) = -(-0.5D0)**I
   20 CONTINUE
      CALL  ARCOEF( PAR,L,CMA )
C
      OUTMIN = -1.0D30
      OUTMAX =  1.0D30
      IOPT = 1
C
      RETURN
      E N D
cc      SUBROUTINE  FILTR3( Y,XF,VF,A,B,C,M,MJ,NS,N,OUTMIN,OUTMAX,
cxx      SUBROUTINE  FILTR3( Y,XF,VF,A,B,C,M,NS,N,OUTMIN,OUTMAX,
      SUBROUTINE  FILTR3( Y,XF,VF,A,B,M,NS,N,OUTMIN,OUTMAX,
     *                    FF,OVAR )
C
C  ...  Kalman filter  ...
C
C     Inputs:
C        Y:      time series
C        N:      data length
C        NS:     Start position of filtering
C        XF:     Initial state vector
C        VF:     Initial covariance matrix
C        A:      Parameters of the matrix F
C        B:      Parameters of the matrix G
C        C:      Parameters of the matrix H
C        K:      Dimension of the state vector
C        MJ:     Adjustable dimension of XF, VF
C        OUTMIN: Lower limit for detecting outliers
C        OUTMAX: Upper limit for detecting outliers
C     Outputs:
C        FF:     Log likelihood
C        SIG2:   Estimated variance
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  A(MJ), B(MJ), C(MJ), Y(N)
cc      DIMENSION  XF(MJ), VF(MJ,MJ), XP(40), VP(40,40)
cc      DIMENSION  WRK(40,40), VH(40), GAIN(40)
cxx      DIMENSION  A(M), B(M), C(M), Y(N)
cxx      DIMENSION  XF(M), VF(M,M), XP(M), VP(M,M)
cxx      DIMENSION  WRK(M,M), VH(M), GAIN(M)
C
      INTEGER M, NS, N
      DOUBLE PRECISION Y(N), VF(M,M), A(M), B(M), OUTMIN, OUTMAX, FF,
     1                 OVAR
c local
      INTEGER I, II, J, NSUM
      DOUBLE PRECISION XF(M), XP(M), VP(M,M), WRK(M,M), VH(M), GAIN(M),
     1                 PI, SDET, PVAR, PERR
C
      DATA   PI  /3.1415926535D0/
C
      OVAR = 0.0D0
      SDET = 0.0D0
      NSUM = 0
C
      DO 300  II=NS,N
C
C  ...  ONE STEP AHEAD PREDICTION  ...
C
      XP(M) = A(M)*XF(1)
      DO 100  I=1,M-1
cxx  100 XP(I) = A(I)*XF(1) + XF(I+1)
      XP(I) = A(I)*XF(1) + XF(I+1)
  100 CONTINUE
C
cxx      DO 110  J=1,M
      DO 111  J=1,M
      WRK(M,J) = A(M)*VF(1,J)
      DO 110  I=1,M-1
cxx  110 WRK(I,J) = A(I)*VF(1,J) + VF(I+1,J)
      WRK(I,J) = A(I)*VF(1,J) + VF(I+1,J)
  110 CONTINUE
  111 CONTINUE
C
cxx      DO 120  I=1,M
      DO 121  I=1,M
      VP(I,M) = WRK(I,1)*A(M)
      DO 120  J=1,M-1
cxx  120 VP(I,J) = WRK(I,1)*A(J) + WRK(I,J+1)
      VP(I,J) = WRK(I,1)*A(J) + WRK(I,J+1)
  120 CONTINUE
  121 CONTINUE
C
cxx      DO 140  I=1,M
      DO 141  I=1,M
      DO 140  J=1,M
cxx  140 VP(I,J) = VP(I,J) + B(I)*B(J)
      VP(I,J) = VP(I,J) + B(I)*B(J)
  140 CONTINUE
  141 CONTINUE
C
C  ...  FILTERING  ...
C
      IF( Y(II).GT.OUTMIN .AND. Y(II).LT.OUTMAX .AND. II.LE.N )  THEN
C
      DO 210  I=1,M
cxx  210 VH(I) = VP(I,1)
      VH(I) = VP(I,1)
  210 CONTINUE
C
      PVAR = VH(1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF( PVAR.LE.1.0D-30 )  THEN
         FF = -1.0D20
         RETURN
      END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PERR = Y(II) - XP(1)
C
      DO 230  I=1,M
cxx  230 GAIN(I) = VH(I)/PVAR
      GAIN(I) = VH(I)/PVAR
  230 CONTINUE
C
      DO 250  I=1,M
cxx  250 XF(I) = XP(I) + GAIN(I)*PERR
      XF(I) = XP(I) + GAIN(I)*PERR
  250 CONTINUE
C
cxx      DO 260  J=1,M
      DO 261  J=1,M
      DO 260  I=1,M
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
cxx      DO 270  I=1,M
      DO 271  I=1,M
      XF(I) = XP(I)
      DO 270  J=1,M
cxx  270 VF(I,J) = VP(I,J)
      VF(I,J) = VP(I,J)
  270 CONTINUE
  271 CONTINUE
      END IF
C
  300 CONTINUE
      OVAR = OVAR/NSUM
      FF = -0.5D0*(NSUM*DLOG(PI*2*OVAR) + SDET + NSUM)
C
      RETURN
      E N D
cc      SUBROUTINE  ISTAT3( M,L,AR,CMA,MJ,XF,VF )
      SUBROUTINE  ISTAT3( M,L,MM,AR,CMA,XF,VF,IER )
C
C  ...  Initial state ...
C
C     Inputs:
C        M:     AR order
C        L:     MA order
C        AR:    AR coefficient
C        CMA:   MA coefficient
C        MJ:    Adjustable dimension of F
C     Outputs:
C         XF:   State vector, X(0|0)
C         VF:   State covarance matrix, V(0|0)
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  AR(*), CMA(*)
cc      DIMENSION  XF(MJ), VF(MJ,MJ), COV(0:20), G(0:20)
cxx      DIMENSION  AR(M), CMA(L)
cxx      DIMENSION  XF(MM), VF(MM,MM), COV(0:MM), G(0:MM)
C
      INTEGER M, L, MM, IER 
      DOUBLE PRECISION AR(M), CMA(L), XF(MM), VF(MM,MM)
c local
      INTEGER I, I1, J, J1, JMIN
      DOUBLE PRECISION COV(0:MM), G(0:MM), SUM
C
cc      MM = MAX0( M,L+1 )
cc      DO 10  I=1,MJ
cxx      DO 10  I=1,MM
      DO 11  I=1,MM
      XF(I) = 0.0D0
cc      DO 10  J=1,MJ
      DO 10  J=1,MM
cxx   10 VF(I,J) = 0.0D0
      VF(I,J) = 0.0D0
   10 CONTINUE
   11 CONTINUE
C
cc      CALL  ARMCOV( M,L,AR,CMA,1.0D0,MM,COV )
      CALL  ARMCOV( M,L,AR,CMA,1.0D0,MM,COV,MM,IER )
      if( ier.ne.0 ) return
      CALL  IMPULS( M,L,AR,CMA,MM,G )
      VF(1,1) = COV(0)
      DO 50 I=2,MM
      SUM = 0.0D0
      DO 30 J=I,M
cxx   30 SUM = SUM + AR(J)*COV(J-I+1)
      SUM = SUM + AR(J)*COV(J-I+1)
   30 CONTINUE
      DO 40 J=I-1,L
cxx   40 SUM = SUM - CMA(J)*G(J-I+1)
      SUM = SUM - CMA(J)*G(J-I+1)
   40 CONTINUE
      VF(1,I) = SUM
cxx   50 VF(I,1) = SUM
      VF(I,1) = SUM
   50 CONTINUE
cxx      DO 100 I=2,MM
      DO 101 I=2,MM
      DO 100 J=I,MM
      SUM = 0.0D0
cxx      DO 60 I1=I,M
      DO 61 I1=I,M
      DO 60 J1=J,M
cxx   60 SUM = SUM + AR(I1)*AR(J1)*COV(IABS(J1-J-I1+I))
      SUM = SUM + AR(I1)*AR(J1)*COV(IABS(J1-J-I1+I))
   60 CONTINUE
   61 CONTINUE
cxx      DO 70 I1=I,M
      DO 71 I1=I,M
c  modified 2019/12/09 =============
ccc   DO 70 J1=J-I+I1,L
      JMIN = MAX(J-1,J-I+I1)
      DO 70 J1=JMIN,L
c ==================================
cxx   70 SUM = SUM - AR(I1)*CMA(J1)*G(IABS(J1-J-I1+I))
      SUM = SUM - AR(I1)*CMA(J1)*G(IABS(J1-J-I1+I))
   70 CONTINUE
   71 CONTINUE
cxx      DO 80 I1=J,M
      DO 81 I1=J,M
c  modified 2019/12/09 =============
ccc   DO 80 J1=I-J+I1,L
      JMIN = MAX( I-1,I-J+I1)
      DO 80 J1=JMIN,L
c ==================================
cxx   80 SUM = SUM - AR(I1)*CMA(J1)*G(IABS(J1-I-I1+J))
      SUM = SUM - AR(I1)*CMA(J1)*G(IABS(J1-I-I1+J))
   80 CONTINUE
   81 CONTINUE
c  modified 2019/12/09 =============
ccc      DO 90 I1=I-1,L
      DO 90 I1=I-1,L+I-J
ccc      SUM = SUM + CMA(I1)*CMA(J1)
      SUM = SUM + CMA(I1)*CMA(J-I+I1)
   90 CONTINUE
c ==================================
      VF(I,J) = SUM
cxx  100 VF(J,I) = SUM
      VF(J,I) = SUM
  100 CONTINUE
  101 CONTINUE
C
      RETURN
      E N D
cc      SUBROUTINE  SETABC( M,L,AR,CMA,A,B,C )
      SUBROUTINE  SETABC( M,L,AR,CMA,A,B,C,MM )
C
C  ...  State space model for Seasonal Adjustment  ...
C
C     Input:
C       M:     AR Order
C       L:     MA Order
C       AR:    AR coefficient
C       CMA:   MA coefficient
C     Outputs:
C       A,B,C: Parameters of F, G, H
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  A(*), B(*), C(*)
cc      DIMENSION  AR(*), CMA(*)
cxx      DIMENSION  A(MM), B(MM), C(MM)
cxx      DIMENSION  AR(M), CMA(L)
C
      INTEGER M, L, MM
      DOUBLE PRECISION AR(M), CMA(L), A(MM), B(MM), C(MM)
c local
      INTEGER I
C
cc      MM = MAX0( M,L+1 )
cxx      DO 10  I=1,MM
cxx      C(I) = 0.0D0
cxx      A(I) = 0.0D0
cxx   10 B(I) = 0.0D0
      C(1:MM) = 0.0D0
      A(1:MM) = 0.0D0
      B(1:MM) = 0.0D0

C
      DO 20  I=1,M
cxx   20 A(I) = AR(I)
      A(I) = AR(I)
   20 CONTINUE
      B(1) = 1.0D0
      DO 30 I=1,L
cxx   30 B(I+1) =-CMA(I)
      B(I+1) =-CMA(I)
   30 CONTINUE
      C(1) = 1.0D0
C
      RETURN
      E N D
cc      SUBROUTINE  FFARMA( K,AA,FF,GDUMMY,IFG )
cxx      SUBROUTINE  FFARMA( K,AA,FF,GDUMMY,IFG,
      SUBROUTINE  FFARMA( K,AA,FF,IFG,
     * Y,N,M,L,MM,OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
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
cc      DIMENSION  PAR(20), GDUMMY(20)
cc      DIMENSION  AA(20), AR(20), CMA(20)
cc      DIMENSION  A(20), B(20), C(20)
cc      DIMENSION  XF(40), VF(40,40)
cc      COMMON  /C92826/   Y(1000)
cc      COMMON  /C92825/   OUTMIN, OUTMAX
cc      COMMON  /C92907/  ALIMIT
cc      COMMON  /C92908/   M, L, N
cc      COMMON  /C92909/   FLK, SIG2
cxx      DIMENSION  PAR(MM)
cxx      DIMENSION  AA(K), AR(M), CMA(L)
cxx      DIMENSION  A(MM), B(MM), C(MM)
cxx      DIMENSION  XF(MM), VF(MM,MM)
cxx      DIMENSION  Y(N)
C
      INTEGER K, IFG, N, M, L, MM, IER
      DOUBLE PRECISION AA(K), FF, Y(N), OUTMIN, OUTMAX, ALIMIT, FLK,
     1                 SIG2
c local
      INTEGER I
      DOUBLE PRECISION PAR(MM), AR(M), CMA(L), A(MM), B(MM), C(MM),
     1                 XF(MM), VF(MM,MM)
C avoid floating-point exceptions
      INTEGER IFPLIM
cc      IFPLIM = 709
      IFPLIM = 87
C
cc      MJ  = 20
C
cc      MM = MAX0( M,L+1 )
c
      IER = 0
      DO 101 I=1,M
      IF( DABS(AA(I)).GT.IFPLIM ) IER = -1
cxx      IF( DABS(AA(I)).GT.20.0D0 )  GO TO 100
      IF( DABS(AA(I)).GT.30.0D0 )  GO TO 100
  101 CONTINUE
c
      DO 10 I=1,M
cxx   10 PAR(I)  = ALIMIT*(DEXP(AA(I))-1.0D0)/(DEXP(AA(I))+1.0D0)
      PAR(I)  = ALIMIT*(DEXP(AA(I))-1.0D0)/(DEXP(AA(I))+1.0D0)
   10 CONTINUE
      CALL  ARCOEF( PAR,M,AR )
      DO 20 I=1,L
cxx   20 PAR(I)  = ALIMIT*(DEXP(AA(M+I))-1.0D0)/(DEXP(AA(M+I))+1.0D0)
      PAR(I)  = ALIMIT*(DEXP(AA(M+I))-1.0D0)/(DEXP(AA(M+I))+1.0D0)
   20 CONTINUE
      CALL  ARCOEF( PAR,L,CMA )
      IFG = 0
C
cc      CALL  SETABC( M,L,AR,CMA,A,B,C )
cc      CALL  ISTAT3( M,L,AR,CMA,MJ,XF,VF )
cc      CALL  FILTR3( Y,XF,VF,A,B,C,MM,MJ,1,N,OUTMIN,OUTMAX,FLK,SIG2 )
      CALL  SETABC( M,L,AR,CMA,A,B,C,MM )
      CALL  ISTAT3( M,L,MM,AR,CMA,XF,VF,IER )
      if( ier.ne.0 ) return
cxx      CALL  FILTR3( Y,XF,VF,A,B,C,MM,1,N,OUTMIN,OUTMAX,FLK,SIG2 )
      CALL  FILTR3( Y,XF,VF,A,B,MM,1,N,OUTMIN,OUTMAX,FLK,SIG2 )
      FF = -FLK
      RETURN
C
  100 IFG = 1
      FF = 1.0D20
      RETURN
      E N D
cc      SUBROUTINE  DAVIDN( FUNCT, X, N, NDIF )
      SUBROUTINE  DAVIDN( FUNCT,X,N,NDIF, YY,NN,M,L,MLMAX,
     *                   OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
C
C  ...  6/20/83, 12/19/92
C
cxx      IMPLICIT  REAL*8( A-H,O-Z )
cc      DIMENSION  X(40), DX(40), G(40), G0(40), Y(40)
cc      DIMENSION  H(40,40), WRK(40), S(40)
cc      COMMON     / CCC /  ISW, IPR, ISMT, IDIF
cc      COMMON     / DDD /  XM , AIC , SD
cxx      DIMENSION  X(N), DX(N), G(N), G0(N), Y(N)
cxx      DIMENSION  H(N,N), WRK(N), S(N)
cxx      DIMENSION  YY(NN)
C
      INTEGER N, NDIF, NN, M, L, MLMAX, IER
      DOUBLE PRECISION X(N), YY(NN), OUTMIN, OUTMAX, ALIMIT, FLK, SIG2
c local
      INTEGER I, IDIF, IPR, ICOUNT, ICC, IG, ISW, J
      DOUBLE PRECISION DX(N), G(N), G0(N), Y(N), H(N,N), WRK(N), S(N),
     1                 TAU2, EPS1, EPS2, RAMDA, CONST1, SUM, S1, S2,
     2                 STEM, SS, ED, XM, XMB
C
      EXTERNAL  FUNCT
      DATA        TAU2  /          1.0D-6  /
      DATA  EPS1, EPS2  / 1.0D-6 , 1.0D-6  /
      RAMDA  = 0.5D0
      CONST1 = 1.0D-30
      IPR = 0
c-------   
      IDIF = NDIF
c-------
C
C          INITIAL ESTIMATE OF INVERSE OF HESSIAN
C
      ICOUNT = 0
 1000 CONTINUE
cxx 1000 DO 20  I=1,N
cxx      DO 10  J=1,N
cxx   10 H(I,J) = 0.0D0
cxx      S(I)   = 0.0D0
cxx      DX(I)  = 0.0D0
cxx   20 H(I,I) = 1.0D0
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
cxx      IF( NDIF.EQ.0 ) CALL  FUNCT ( N,X,XM,G,IG, YY,NN,M,L,MLMAX,
      IF( NDIF.EQ.0 ) CALL  FUNCT ( N,X,XM,IG, YY,NN,M,L,MLMAX,
     *                    OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
      IF( NDIF.GE.1 ) CALL  FUNCND( FUNCT,N,X,XM,G,IG,YY,NN,M,L,MLMAX,
     *                    OUTMIN,OUTMAX,ALIMIT,ISW,IDIF,FLK,SIG2,IER )
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
      CALL  LINEAR( FUNCT,X,S,RAMDA,ED,N,IG, YY,NN,M,L,
cxx     *         MLMAX,OUTMIN,OUTMAX,ALIMIT,ISW,IDIF,FLK,SIG2,IER )
     *         MLMAX,OUTMIN,OUTMAX,ALIMIT,ISW,FLK,SIG2,IER )
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
cxx      IF( NDIF.EQ.0 )  CALL  FUNCT ( N,X,XM,G,IG, YY,NN,M,L,MLMAX,
      IF( NDIF.EQ.0 )  CALL  FUNCT ( N,X,XM,IG, YY,NN,M,L,MLMAX,
     *                     OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
      IF( NDIF.GE.1 )  CALL  FUNCND( FUNCT,N,X,XM,G,IG,YY,NN,M,L,MLMAX,
     *                     OUTMIN,OUTMAX,ALIMIT,ISW,IDIF,FLK,SIG2,IER )
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
      SUBROUTINE  FUNCND( FUNCT,M,A,F,G,IFG, Y,N,MM,L,MLMAX,
     *                    OUTMIN,OUTMAX,ALIMIT,ISW,IDIF,FLK,SIG2,IER )
C
C  ...  FUNCTION EVALUATION AND NUMERICAL DIFFERENCING  ...
C
cxx      IMPLICIT   REAL*8( A-H,O-Z )
cc      DIMENSION  A(M) , G(M) , B(20),GD(5)
cc      COMMON  / CCC /  ISW , IPR, ISMT, IDIF
cc      COMMON  /CMFUNC/ DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(20),GC(20)
cxx      DIMENSION  A(M) , G(M) , B(M)
cxx      DIMENSION  Y(N)
C
      INTEGER M, IFG, N, MM, L, MLMAX, ISW, IDIF, IER
      DOUBLE PRECISION A(M), F, G(M), Y(N), OUTMIN, OUTMAX, ALIMIT,
     1                 FLK, SIG2
c local
      INTEGER I, II
      DOUBLE PRECISION B(M), CONST, FB, FF
      EXTERNAL FUNCT
C
C     DATA       ICNT /0/
      CONST = 0.00001D0
C
cc      CALL  FUNCT( M,A,F,GD,IFG )
cxx      CALL  FUNCT( M,A,F,GD,IFG, Y,N,MM,L,MLMAX,
      CALL  FUNCT( M,A,F,IFG, Y,N,MM,L,MLMAX,
     *                    OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
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
cxx      CALL  FUNCT( M,B,FF,GD,IFG, Y,N,MM,L,MLMAX,
      CALL  FUNCT( M,B,FF,IFG, Y,N,MM,L,MLMAX,
     *                    OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
      if( ier.ne.0 ) return
      IF( IDIF .EQ. 1 )  GO TO 20
      B(II) = A(II) - CONST
cc      CALL  FUNCT( M,B,FB,GD,IFG )
cxx      CALL  FUNCT( M,B,FB,GD,IFG, Y,N,MM,L,MLMAX,
      CALL  FUNCT( M,B,FB,IFG, Y,N,MM,L,MLMAX,
     *                    OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
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
      SUBROUTINE  LINEAR( FUNCT,X,H,RAM,EE,K,IG, Y,N,M,L,
cxx     * MLMAX,OUTMIN,OUTMAX,ALIMIT,ISW,IDIF,FLK,SIG2,IER )
     * MLMAX,OUTMIN,OUTMAX,ALIMIT,ISW,FLK,SIG2,IER )

C
C  ...  LINEAR SEARCH  ...
C
cxx      IMPLICIT  REAL*8( A-H,O-Z )
cc      INTEGER  RETURN,SUB
cc      DIMENSION  X(1), H(1), X1(K)
cc      DIMENSION  G(40)
cc      COMMON     / CCC /  ISW , IPR, ISMT, IDIF
cxx      DIMENSION  X(K), H(K), X1(K)
cxx      DIMENSION Y(N)
C
      INTEGER K, IG, N, M, L, MLMAX, ISW, IER
      DOUBLE PRECISION X(K), H(K), RAM, EE, Y(N), OUTMIN, OUTMAX,
     1                 ALIMIT, FLK, SIG2
c local
      INTEGER I, IRET, ISUB, IFG, ire510
      DOUBLE PRECISION X1(K), CONST2, HNORM, RAM1, RAM2, RAM3, E1, E2,
     1                 E3, A1, A2, A3, B1, B2
      EXTERNAL FUNCT
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
cxx      CALL  FUNCT( K,X1,E2,G,IG, Y,N,M,L,MLMAX,
      CALL  FUNCT( K,X1,E2,IG, Y,N,M,L,MLMAX,
     *                    OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
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
cxx      CALL  FUNCT( K,X1,E3,G,IG, Y,N,M,L,MLMAX,
      CALL  FUNCT( K,X1,E3,IG, Y,N,M,L,MLMAX,
     *                    OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
      if( ier.ne.0 ) return
       IF( IG.EQ.1 )  GO TO  500
cc      IF( IPR.GE.7 )  WRITE(6,3)  RAM3,E3
      IF( E3 .GT. E2 )  GO TO 70
c----- for debug
c         IF( E3 .EQ. E2 )  GO TO 70
c-----
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
cxx      CALL  FUNCT( K,X1,E2,G,IG, Y,N,M,L,MLMAX,
      CALL  FUNCT( K,X1,E2,IG, Y,N,M,L,MLMAX,
     *                    OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
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
cxx      CALL  FUNCT( K,X1,EE,G,IG, Y,N,M,L,MLMAX,
      CALL  FUNCT( K,X1,EE,IG, Y,N,M,L,MLMAX,
     *                    OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
      if( ier.ne.0 ) return
cc      IF(IPR.GE.7)  WRITE(6,5)  RAM,EE
C
      IFG = 0
cc      ASSIGN  300 TO  SUB
cc      ASSIGN 200 TO SUB
cc   95 ASSIGN 130 TO RETURN
      ISUB = 300
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
cxx      CALL  FUNCT( K,X1,EE,G,IG, Y,N,M,L,MLMAX,
      CALL  FUNCT( K,X1,EE,IG, Y,N,M,L,MLMAX,
     *                    OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
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
cxx 19/12/07
         ire510 = 0
cxx
  510 DO 520  I=1,K
cxx  520 X1(I) = X(I) + RAM*H(I)
      X1(I) = X(I) + RAM*H(I)
  520 CONTINUE
cc      CALL  FUNCT( K,X1,E3,G,IG )
cxx      CALL  FUNCT( K,X1,E3,G,IG, Y,N,M,L,MLMAX,
      CALL  FUNCT( K,X1,E3,IG, Y,N,M,L,MLMAX,
     *                    OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
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
cxx      GO TO 510
cxx 19/12/07
         if (RAM .EQ. RAM2) return
         ire510 = ire510 + 1
         if (ire510 .le. 100) GO TO 510
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
