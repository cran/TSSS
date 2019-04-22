      INCLUDE 'TSSS_f.h'
C------------------------------------------------- 
C     ARCOEF  ---  arfit, armaft, season, tvar
C     FOUGER  ---  spgrh, tvar
C     PARCOR  ---  arfit, armaft, season, tvar
C     REGRES  ---  lsqr, arfit, lsar1, lsar2, polreg
C     COMAIC  ---  lsqr, arfit, lsar1, lsar2, polreg
C     HUSHLD  ---  lsqr, arfit, marlsq, lsar1, lsar2 
C     RECOEF  ---  lsqr, arfit, lsar1, lsar2, polreg
C     REDUCT  ---  lsqr, arfit, lsar1, lsar2, polreg
C     SETXAR  ---  arfit, lsar1, lsar2
C     ARYULE  ---  arma, arfit
C     ARMASP  ---  arma, arfit, lsar1, tvspc
C     FOURIE  ---  period, arma, marspc, arfit, tvspc
C     MOMENT ---  trend, season, tvvar, ngsmth
C     SMOOTH  ---  tvvar, smooth
C     GINVRS  ---  tvvar, smooth, season
C     SETSEA  ---  simssm. ngsim
C     CHOLES  ---  simssm. ngsim
C     CRSCOR  ---  crscor, marfit
C     MEAN  ---  unicor, crscor, period, arfit, marfit
C     AUTCOR  --- unicor, arfit
C     AUTCOV  ---  unicor, period, arfit
C     UNICOR  ---  unicor, arfit
C     ERRACF  ---  unicor, arfit
C     WINDOW  ---  period, fftper
C     ARMCOV  ---  arma, armaft, season
C     DECOM  ---  arma, armaft, season
C     IMPULS  ---  arma, armaft, season
C     SOLVE  ---  arma, armaft, season
C     SMOTH1  ---  seasom, tvar
C
C     ID  ---  simssm, ngsim, season
C     INIT  ---  simssm, ngsim
C     GAUSS  ---  densty, klinfo, ngsmth, ngsim
C     PEARSN  ---  densty, ngsmth, ngsim
C     DBLEXP  ---  densty, ngsim, ngsmth
C     CAUCHY  ---  densty, klinfo
C-------------------------------------------------
C
      SUBROUTINE  ARCOEF( PAR,K,A )
C
C  ...  This subroutine computes AR coefficients from PARCOR  ...
C
C     Inputs:
C        PAR(I):   PARCOR
C        K:        Order of AR model
C     Output:
C        A(I):     AR coefficient
C
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cc      DIMENSION  PAR(K), A(K), AA(50)
cxx      DIMENSION  PAR(K), A(K), AA(K)
      INTEGER :: K
      REAL(8) :: PAR(K), A(K), AA(K)
C
      DO 30  II=1,K
      A(II)  = PAR(II)
      AA(II) = PAR(II)
      IF( II-1.LE.0 )  GO TO 30
      DO 10  J=1,II-1
cxx   10 A(J) = AA(J) - PAR(II)*AA(II-J)
      A(J) = AA(J) - PAR(II)*AA(II-J)
   10 CONTINUE
      IF( II.LT.K )  THEN
        DO 20  J=1,II-1
cxx   20   AA(J) = A(J)
        AA(J) = A(J)
   20   CONTINUE
      END IF
   30 CONTINUE
      RETURN
      E N D
c
C
      SUBROUTINE FOUGER(G,LGP1,FC,FS,LF1)
C
C     FOURIER TRANSFORM (GOERTZEL METHOD)
C     THIS SUBROUTINE COMPUTES FOURIER TRANSFORM OF G(I),I=0,1,...,LG AT
C     FREQUENCIES K/(2*LF), K=0,1,...,LF AND RETURNS COSINE TRANSFORM IN
C     FC(K) AND SINE TRANSFORM IN FS(K).
C
C       INPUTS:
C          G(I):   ORIGINAL DATA (I=0,1,...,LG)
C          LG1:    = LG+1
C          LF1:    = LF+1
C
C       OUTPUTS:
C          FC(I):  COSINE TRANSFORM OF G  (I=0,1,...,LF)
C          FB(I):  SINE TRANSFORM OF G  (I=0,1,...,LF)
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION G(LGP1),FC(LF1),FS(LF1)
      INTEGER :: LGP1, LF1
      REAL(8) :: G(LGP1),FC(LF1),FS(LF1)
      REAL(8) :: PI, ALF, T, AK, TK, CK, SK, CK2, UM0, UM1, UM2
      LG=LGP1-1
      LF=LF1-1
C     REVERSAL OF G(I),I=1,...,LGP1 INTO G(LG3-I)   LG3=LGP1+1
      IF(LGP1.LE.1) GO TO 110
      LG3=LGP1+1
      LG4=LGP1/2
      DO 100 I=1,LG4
      I2=LG3-I
      T=G(I)
      G(I)=G(I2)
cxx  100 G(I2)=T
      G(I2)=T
  100 CONTINUE
  110 PI=3.1415926536D0
      ALF=LF
      T=PI/ALF
      DO 10 K=1,LF1
      AK=K-1
      TK=T*AK
      CK=DCOS(TK)
      SK=DSIN(TK)
      CK2=CK+CK
      UM1=0.0D0
      UM2=0.0D0
      IF(LG.EQ.0) GO TO 12
      DO 11 I=1,LG
      UM0=CK2*UM1-UM2+G(I)
      UM2=UM1
cxx   11 UM1=UM0
      UM1=UM0
   11 CONTINUE
   12 FC(K)=CK*UM1-UM2+G(LGP1)
      FS(K)=SK*UM1
   10 CONTINUE
      RETURN
      END
C
C
      SUBROUTINE  PARCOR( A,K,PAR )
C
C  ...  This subriutine computes PARCOR for AR coefficients  ...
C
C       Inputs:
C          A:   Vector of AR coefficients
C          K:   Order of the model
C
C       Output:
C          PAR: Vector of partial autocorrelations (PARCOR)
C
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cc      DIMENSION  A(K), PAR(K), G(50)
cxx      DIMENSION  A(K), PAR(K), G(K)
      INTEGER :: K
      REAL(8) :: A(K), PAR(K)
      REAL(8) :: G(K), S
C
      DO 10  I=1,K
cxx   10 PAR(I) = A(I)
      PAR(I) = A(I)
   10 CONTINUE
C
c-----
      IF( K .EQ. 1 )   RETURN                                           
c-----
      DO 40  II=K-1,1,-1
      S = 1.D0 - PAR(II+1)**2
      DO 20  I=1,II
cxx   20 G(I) = (PAR(I) + PAR(II+1)*PAR(II-I+1))/S
      G(I) = (PAR(I) + PAR(II+1)*PAR(II-I+1))/S
   20 CONTINUE
      I2 = (II+1)/2
      IF( MOD( II,2 ).EQ.1 )  G(I2) = PAR(I2)/(1.D0 - PAR(II+1))
      DO 30  I=1,II
cxx   30 PAR(I) = G(I)
      PAR(I) = G(I)
   30 CONTINUE
   40 CONTINUE
C
      RETURN
      E N D
C
C
cc      SUBROUTINE  REGRES( X,K,N,MJ1,MJ2,A,SIG2,AIC,IMIN )
      SUBROUTINE  REGRES( X,K,N,MJ1,A,SIG2,AIC,IMIN )
C
C  ...  Regression model fitting  ...
C  ...  Order of the model is selected by AIC  ...
C
C     Inputs:
C        X:      Householder reduced form (upper triangular matrix)
C        K:      Maximum number of regressors
C        N:      Number of data
C        MJ1:    Adjustable dimension of X
C        MJ2:    Adjustable dimension of A, ISG2 and AIC
C     Outputs:
C        A(I,M): Regression coefficients of the model with order M
C        SIG2:   Residual variances
C        AIC:    AIC's
C        IMIN:   MAICE order
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  X(MJ1,1), A(MJ2,MJ2), SIG2(0:MJ2), AIC(0:MJ2)
cxx      DIMENSION  X(MJ1,K+1), A(K,K), SIG2(0:K), AIC(0:K)
      INTEGER :: K, N, MJ1, IMIN
      REAL(8) :: X(MJ1,K+1), A(K,K), SIG2(0:K), AIC(0:K)
      REAL(8) :: AICM
C
      CALL  COMAIC( X,N,K,MJ1,SIG2,AIC )
C
      IMIN = 0
      AICM = AIC(0)
      DO 10  M=1,K
      IF( AIC(M).LT.AICM )  THEN
         IMIN = M
         AICM = AIC(M)
      END IF
      CALL  RECOEF( X,M,K,MJ1,A(1,M) )
   10 CONTINUE
C
      RETURN
      E N D
C
C
      SUBROUTINE  COMAIC( X,N,K,MJ1,SIG2,AIC )
C
C  ...  This subroutine computes residual variances and AIC's  ...
C
C     Inputs:
C        X(I,J):  Householder reduced form
C        N:       Data length
C        K:       Highest order
C        MJ1:     Adjustable dimension of X
C     Outputs:
C        SIG2(I): residual variance of the model with order I
C        AIC(I):  AIC of the model with order I
C
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cc      DIMENSION  X(MJ1,1), AIC(0:K), SIG2(0:K)
cxx      DIMENSION  X(MJ1,K+1), AIC(0:K), SIG2(0:K)
      INTEGER :: N, K, MJ1
      REAL(8) :: X(MJ1,K+1), SIG2(0:K), AIC(0:K)
      REAL(8) :: PI2, PVAR
      DATA  PI2/6.28318531D0/
C
      PVAR = 0.0D0
      DO 10 I=K,0,-1
      PVAR = PVAR + X(I+1,K+1)**2
      SIG2(I) = PVAR / N
cxx   10 AIC(I)  = N*DLOG( PI2*SIG2(I) ) + N + 2*(I+1)
      AIC(I)  = N*DLOG( PI2*SIG2(I) ) + N + 2*(I+1)
   10 CONTINUE
C
      RETURN
      E N D
C
C
cc      SUBROUTINE  HUSHLD( X,D,MJ1,N,K )
      SUBROUTINE  HUSHLD( X,MJ1,N,K )
C
C  ...  Householder transformation  ...
C
C     Inputs:
C        X(I,J):  Original matrix
C        D(I):    Working area
C        MJ1:     Adjustable dimension of X
C        N:       Number of rows of X
C        K:       Number of columns of X
C     Output:
C        X(I,J):  Householder reduced form (upper triangular form)
C
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cc      DIMENSION  X(MJ1,1), D(MJ1)
cxx      DIMENSION  X(MJ1,K), D(MJ1)
      INTEGER :: MJ1, N, K
      REAL(8) :: X(MJ1,K)
      REAL(8) :: D(MJ1), TOL, H, G, F, S
C
      TOL = 1.0D-60
C
cxx      DO 100  II=1,K
      DO 101  II=1,K
        H = 0.0D0
        DO 10  I=II,N
          D(I) = X(I,II)
cxx   10     H = H + D(I)**2
          H = H + D(I)**2
   10 CONTINUE
        IF( H .GT. TOL )  GO TO 20
        G = 0.0D0
        GO TO 100
   20   G = DSQRT( H )
        F = X(II,II)
        IF( F .GE. 0.0D0 )   G = -G
        D(II) = F - G
        H = H - F*G
C
        DO 30  I=II+1,N
cxx   30   X(I,II) = 0.0D0
        X(I,II) = 0.0D0
   30   CONTINUE
        DO 60  J=II+1,K
          S = 0.0D0
          DO 40  I=II,N
cxx   40     S = S + D(I)*X(I,J)
            S = S + D(I)*X(I,J)
   40     CONTINUE
          S = S/H
          DO 50  I=II,N
cxx   50     X(I,J) = X(I,J) - D(I)*S
            X(I,J) = X(I,J) - D(I)*S
   50     CONTINUE
   60   CONTINUE
  100 X(II,II) = G
  101 CONTINUE
C
      RETURN
      E N D
C
C
      SUBROUTINE  RECOEF( X,M,K,MJ,A )
C
C  ...  Regression coefficients  ...
C
C     Inputs:
C        X(I,J):  Householder reduced form
C        M:       Number of actually used regressors
C        K:       Heighest order
C        MJ:      Adjustable dimension of X
C     Output:
C        A(I):    Vector of regression coefficients
C
cxx      IMPLICIT REAL*8 (A-H,O-Z)
cc      DIMENSION  X(MJ,1), A(1)
cxx      DIMENSION  X(MJ,K+1), A(M)
      INTEGER :: M, K, MJ
      REAL(8) :: X(MJ,K+1), A(M)
      REAL(8) :: SUM
C
      A(M) = X(M,K+1)/X(M,M)
c-----
      IF( M .EQ. 1 ) RETURN
c-----
      DO 20 I=M-1,1,-1
      SUM = X(I,K+1)
      DO 10 J=I+1,M
cxx   10 SUM  = SUM - A(J)*X(I,J)
      SUM  = SUM - A(J)*X(I,J)
   10 CONTINUE
cxx   20 A(I) = SUM/X(I,I)
      A(I) = SUM/X(I,I)
   20 CONTINUE
C
      RETURN
      E N D
C
C
cc      SUBROUTINE  REDUCT( SETX,Z,D,NMK,N0,K,MJ1,X )
      SUBROUTINE  REDUCT1( SETX,Z,NMK,N0,K,MJ1,X )
C
C  ...  Successive Householder reduction  ...
C
C     Inputs:
C        SETX:    Name of the subroutine for making X(I,J)
C        Z(I):    Data vector
C        D(I):    Working area
C        NMK:     Number of actually used observations
C        N0:      Time point of the previous set ofobservations
C        K:       Heighest order of the model
C        MJ1:     Adjustable dimension of X
C     Output:
C        X(I,J):  data matrix
C
cxx      IMPLICIT  REAL*8( A-H,O-Z )
cc      DIMENSION  X(MJ1,1) , D(1), Z(1)
cx      DIMENSION  X(MJ1,1) , Z(1)
cxx      DIMENSION  X(MJ1,K+1) , Z(N0+NMK)
      INTEGER :: NMK, N0, K, MJ1
      REAL(8) :: Z(N0+NMK), X(MJ1,K+1)
C
      L = MIN0( NMK,MJ1 )
      K1 = K + 1
      N1 = L
C
      CALL  SETX( Z,N0,L,K,MJ1,0,X )
cc      CALL  HUSHLD( X,D,MJ1,L,K1 )
      CALL  HUSHLD( X,MJ1,L,K1 )
      IF( N1 .GE. NMK )  RETURN
C
   10 L = MIN0( NMK-N1,MJ1-K1 )
C
      LK = L + K1
      N2 = N0 + N1
      CALL  SETX( Z,N2,L,K,MJ1,1,X )
cc      CALL  HUSHLD( X,D,MJ1,LK,K1 )
      CALL  HUSHLD( X,MJ1,LK,K1 )
      N1 = N1 + L
      IF( N1.LT.NMK )  GO TO 10
C
      RETURN
C
      E N D
C
C
      SUBROUTINE  REDUCT( SETX,Z,NMK,N0,K,MJ1,X )
C
C  ...  Successive Householder reduction  ...
C
C     Inputs:
C        SETX:    Name of the subroutine for making X(I,J)
C        Z(I):    Data vector
C        D(I):    Working area
C        NMK:     Number of actually used observations
C        N0:      Time point of the previous set ofobservations
C        K:       Heighest order of the model
C        MJ1:     Adjustable dimension of X
C     Output:
C        X(I,J):  data matrix
C
cxx      IMPLICIT  REAL*8( A-H,O-Z )
cc      DIMENSION  X(MJ1,1) , D(1), Z(1)
cx      DIMENSION  X(MJ1,1) , Z(1)
cxx      DIMENSION  X(MJ1,K+1) , Z(N0+NMK+K)
      INTEGER :: NMK, N0, K, MJ1
      REAL(8) :: Z(N0+NMK+K), X(MJ1,K+1)
C
      L = MIN0( NMK,MJ1 )
      K1 = K + 1
      N1 = L
C
cdd      CALL  SETX( Z,N0,L,K,MJ1,0,X )
      CALL SETXAR( Z,N0,L,K,MJ1,0,X )
cc      CALL  HUSHLD( X,D,MJ1,L,K1 )
      CALL  HUSHLD( X,MJ1,L,K1 )
      IF( N1 .GE. NMK )  RETURN
C
   10 L = MIN0( NMK-N1,MJ1-K1 )
C
      LK = L + K1
      N2 = N0 + N1
      CALL  SETX( Z,N2,L,K,MJ1,1,X )
cc      CALL  HUSHLD( X,D,MJ1,LK,K1 )
      CALL  HUSHLD( X,MJ1,LK,K1 )
      N1 = N1 + L
      IF( N1.LT.NMK )  GO TO 10
C
      RETURN
C
      E N D
C
C
      SUBROUTINE  SETXAR( Z,N0,L,K,MJ1,JSW,X )
C
C  ...  Data matrix for AR model  ...
C
C     Inputs:
C        Z(I):    Time series
C        N0:      Origin of the current observations
C        L:       Number of current observations
C        K:       Number of regressors
C        MJ1:     Adjustable dimension of X
C        JSW=0:   Make initial data matrix
C           =1:   Apend L*(K+1) data matrix below the triangular one
C     Output:
C        X(I,J):  Data matrix
C
cc      REAL*8  X(MJ1,1), Z(1)
cxx      REAL*8  X(MJ1,K+1), Z(N0+L+K)
      INTEGER :: N0, L, K, MJ1, JSW
      REAL(8) :: Z(N0+L+K), X(MJ1,K+1)
C
      I0 = 0
      IF( JSW .EQ. 1 )     I0 = K+1
cxx      DO 10  I=1,L
      DO 11  I=1,L
        II = I + I0
        JJ = N0 + K + I
        X(II,K+1) = Z(JJ)
      DO 10  J=1,K
        JJ = JJ - 1
cxx   10 X(II,J) = Z(JJ)
      X(II,J) = Z(JJ)
   10 CONTINUE
   11 CONTINUE
C
      RETURN
      E N D
C
C
      SUBROUTINE ARYULE( C,N,MAXM,SIG2,AIC,PARCOR,A,MAR )
C
C  ...  Yule-Walker method  ...
C
C     Inputs:
C        C(I):    Autocovariance function
C        N:       Data length
C        MAXM:    Highest AR order
C     Outputs:
C        SIG2(I): Innovation variance
C        AIC(I):  AIC
C        PARCOR(I):  PARCOR
C        AMIN:     AR coefficients of the best model
C        MAR:      Selected order of the model
C
cxx      IMPLICIT REAL*8 (A-H,O-Z)
cxx      DIMENSION  C(0:MAXM), SIG2(0:MAXM), AIC(0:MAXM)
cxx      DIMENSION  PARCOR(MAXM), A(MAXM,MAXM)
      INTEGER :: N, MAXM, MAR
      REAL(8) :: C(0:MAXM), SIG2(0:MAXM), AIC(0:MAXM), PARCOR(MAXM),
     1           A(MAXM,MAXM)
      REAL(8) :: CONST, AICMIN, SUM
C
      CONST = N*(DLOG(2*3.1415926535D0) + 1)
C
      SIG2(0) = C(0)
      AIC(0) = CONST + N*DLOG(SIG2(0)) + 2
      AICMIN = AIC(0)
      MAR = 0
C
      DO 50 M=1,MAXM
      SUM  = C(M)
      IF (M.GE.2) THEN
      DO 10 I=1,M-1
cxx   10 SUM  = SUM - A(I,M-1)*C(M-I)
      SUM  = SUM - A(I,M-1)*C(M-I)
   10 CONTINUE
      END IF
      A(M,M) = SUM /SIG2(M-1)
      IF (M.GE.2) THEN
      DO 20 J=1,M-1
cxx   20 A(J,M) = A(J,M-1)-A(M,M)*A(M-J,M-1)
      A(J,M) = A(J,M-1)-A(M,M)*A(M-J,M-1)
   20 CONTINUE
      END IF
      SIG2(M) = SIG2(M-1)*(1.0D0-A(M,M)**2)
      AIC(M) = CONST + N*DLOG(SIG2(M)) + 2*(M+1)
      PARCOR(M) = A(M,M)
      IF( AIC(M).LT.AICMIN )  THEN
         AICMIN = AIC(M)
         MAR = M
      END IF
   50 CONTINUE
      RETURN
      E N D
C
C
      SUBROUTINE ARMASP( A,M,B,L,SIG2,NF,SP )
C
C  ...  Logarithm of the power spectrum of the ARMA model  ...
C
C     Inputs:
C        M:     AR order
C        L:     MA order
C        A(I):  AR coefficient
C        B(I):  MA coefficient
C        SIG2:  Innovation variance
C        NF:    Number of frequencies
C     Output:
C        SP(I): Power spectrum (in log scale)
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION A(M), B(L)
cc      DIMENSION SP(0:NF), H(0:500), FR(0:500), FI(0:500)
cxx      DIMENSION SP(0:NF), H(0:M+L), FR(0:NF), FI(0:NF)
      INTEGER :: M, L, NF
      REAL(8) :: A(M), B(L), SIG2, SP(0:NF)
      REAL(8) :: H(0:M+L), FR(0:NF), FI(0:NF)
C
      H(0) = 1.0D0
      DO 10 I=1,M
cxx   10 H(I) = -A(I)
      H(I) = -A(I)
   10 CONTINUE
C
      CALL  FOURIE( H,M+1,NF+1,FR,FI )
C
      DO 20 I=0,NF
cxx   20 SP(I) = SIG2/( FR(I)**2 + FI(I)**2 )
      SP(I) = SIG2/( FR(I)**2 + FI(I)**2 )
   20 CONTINUE
C
      IF (L .EQ. 0) GO TO 41
      H(0) = 1.0D0
      DO 30 I=1,L
cxx   30 H(I) = -B(I)
      H(I) = -B(I)
   30 CONTINUE
      CALL  FOURIE( H,L+1,NF+1,FR,FI )
      DO 40 I=0,NF
cxx   40 SP(I) = SP(I)*( FR(I)**2 + FI(I)**2 )
      SP(I) = SP(I)*( FR(I)**2 + FI(I)**2 )
   40 CONTINUE
   41 CONTINUE
C
      DO 50 I=0,NF
cxx   50 SP(I) = DLOG10( SP(I) )
      SP(I) = DLOG10( SP(I) )
   50 CONTINUE
C
      RETURN
      E N D
C
C
      SUBROUTINE FOURIE( X,N,M,FC,FS )
C
C  ...  Discrete Fourier transformation by Goertzel method  ...
C
C     Inputs:
C        X(I):   data (I=1,N)
C        N:      data length
C        M:      number of Fourier components
C        FC(J):  Fourier cosine transform (J=1,M)
C        FS(J):  Fourier sine transform   (J=1,M)
C
cxx      IMPLICIT REAL*8 (A-H,O-Z)
cxx      DIMENSION  X(N), FC(M), FS(M)
      INTEGER :: N, M
      REAL(8) :: X(N), FC(M), FS(M)
      REAL(8) :: PI, W, CI, SI, T0, T1, T2
      DATA  PI/3.14159265358979D0/
C
      W = PI/(M-1)
      DO 20 I=1,M
      CI = DCOS(W*(I-1))
      SI = DSIN(W*(I-1))
      T1 = 0.0
      T2 = 0.0
      DO 10 J=N,2,-1
        T0 = 2*CI*T1 - T2 + X(J)
        T2 = T1
cxx   10   T1 = T0
      T1 = T0
   10 CONTINUE
      FC(I) = CI*T1 - T2 + X(1)
cxx   20 FS(I) = SI*T1
      FS(I) = SI*T1
   20 CONTINUE
C
      RETURN
      E N D
      SUBROUTINE  MOMENT( Y,N,YMEAN,VAR )
C
C  ...  Mean and variance of the data  ...
C
C     Inputs:
C       Y(I):   data
C       N:      data length
C     Outputs:
C       YMEAN:  mean
C       VAR:    variance
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  Y(N)
      INTEGER :: N
      REAL(8) :: Y(N), YMEAN, VAR, SUM
C
      SUM = 0.0D0
      DO 10 I=1,N
cxx   10 SUM = SUM + Y(I)
      SUM = SUM + Y(I)
   10 CONTINUE
      YMEAN = SUM/N
      SUM = 0.0D0
      DO 20 I=1,N
cxx   20 SUM = SUM + (Y(I)-YMEAN)**2
      SUM = SUM + (Y(I)-YMEAN)**2
   20 CONTINUE
      VAR = SUM/N
C
      RETURN
      E N D
cc      SUBROUTINE  SMOOTH( F,M,MJ,NDIM,NS,NFE,NPE,VFS,VPS,XFS,XPS,
cc     *                    VSS,XSS )
      SUBROUTINE  SMOOTH( F,M,NDIM,NS,NFE,NPE,VFS,VPS,XFS,XPS,VSS,XSS )
C
C  ...  Fixed Interval Smoother (General Form)  ...
C
C     Inputs:
C        NS:     Start position of filtering
C        NFE:    End position of filtering
C        NPE:    End position of prediction
C        M:      Dimension of the state vector
C        F:      M*M matrix
C        MJ:     Adjustable dimension of XF, VF
C        NDIM:   Adjustable dimension of XFS, XPS, VFS, VPS
C        NMAX    Adjustable dimension of Y
C        VFS:    Covariance matrices of the filter
C        VPS:    Covariance matrices of the predictor
C        XFS:    Mean vectors of the filter
C        XPS:    Mean vectors of the predictor
C     Outputs:
C        VSS:    Covariance matrices of the smoother
C        XSS:    Mean vectors of the smoother
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  F(MJ,MJ)
cc      DIMENSION  XS(40), VS(40,40), VP(40,40)
cc      DIMENSION  XFS(MJ,NDIM), XPS(MJ,NDIM), XSS(MJ,NDIM)
cc      DIMENSION  VFS(MJ,MJ,NDIM), VPS(MJ,MJ,NDIM), VSS(MJ,MJ,NDIM)
cc      DIMENSION  WRK(40,40), SGAIN(40,40)
cxx      DIMENSION  F(M,M)
cxx      DIMENSION  XS(M), VS(M,M), VP(M,M)
cxx      DIMENSION  XFS(M,NDIM), XPS(M,NDIM), XSS(M,NDIM)
cxx      DIMENSION  VFS(M,M,NDIM), VPS(M,M,NDIM), VSS(M,M,NDIM)
cxx      DIMENSION  WRK(M,M), SGAIN(M,M)
      INTEGER :: M, NDIM, NS, NFE, NPE
      REAL(8) :: F(M,M), VFS(M,M,NDIM), VPS(M,M,NDIM), XFS(M,NDIM),
     1           XPS(M,NDIM), VSS(M,M,NDIM), XSS(M,NDIM)
      REAL(8) :: XS(M), VS(M,M), VP(M,M), WRK(M,M), SGAIN(M,M), VDET,
     1           SUM
C
C
C  ...  SMOOTHING  ...
C
cxx      DO 10 II=NFE,NPE
cxx      DO 10  I=1,M
      DO 12 II=NFE,NPE
      DO 11  I=1,M
      XSS(I,II)   = XFS(I,II)
      DO 10  J=1,M
cxx   10 VSS(I,J,II) = VFS(I,J,II)
      VSS(I,J,II) = VFS(I,J,II)
   10 CONTINUE
   11 CONTINUE
   12 CONTINUE
cxx      DO 20  I=1,M
      DO 21  I=1,M
      XS(I)   = XFS(I,NFE)
      DO 20  J=1,M
cxx   20 VS(I,J) = VFS(I,J,NFE)
      VS(I,J) = VFS(I,J,NFE)
   20 CONTINUE
   21 CONTINUE
C
      DO 500  II=NFE-1,NS,-1
C
      NZERO = 0
      DO 100 I=1,M
cxx  100 IF( VFS(I,I,II).GT.1.0D-12 )  NZERO = NZERO + 
      IF( VFS(I,I,II).GT.1.0D-12 )  NZERO = NZERO + 1
  100 CONTINUE
C
      IF( NZERO.EQ.0 )  THEN
cxx         DO 110  I=1,M
         DO 111  I=1,M
         XS(I)     = XFS(I,II)
         XSS(I,II) = XFS(I,II)
         DO 110  J=1,M
         VS(I,J)     = VFS(I,J,II)
cxx  110    VSS(I,J,II) = VFS(I,J,II)
         VSS(I,J,II) = VFS(I,J,II)
  110    CONTINUE
  111    CONTINUE
C
      ELSE
cxx      DO 410  I=1,M
      DO 411  I=1,M
      DO 410  J=1,M
cxx  410 VP(I,J) = VPS(I,J,II+1)
      VP(I,J) = VPS(I,J,II+1)
  410 CONTINUE
  411 CONTINUE
C
cc      CALL  GINVRS( VP,VDET,M,40 )
      CALL  GINVRS( VP,VDET,M )
C
cxx      DO 425  I=1,M
      DO 426  I=1,M
      DO 425  J=1,M
      SUM = 0.0D0
      DO 420  IJ=1,M
cxx  420 SUM = SUM + VFS(I,IJ,II)*F(J,IJ)
      SUM = SUM + VFS(I,IJ,II)*F(J,IJ)
  420 CONTINUE
cxx  425 WRK(I,J) = SUM
      WRK(I,J) = SUM
  425 CONTINUE
  426 CONTINUE
C
cxx      DO 440  I=1,M
      DO 441  I=1,M
      DO 440  J=1,M
      SUM = 0.0D0
      DO 430 IJ=1,M
cxx  430 SUM = SUM + WRK(I,IJ)*VP(IJ,J)
      SUM = SUM + WRK(I,IJ)*VP(IJ,J)
  430 CONTINUE
cxx  440 SGAIN(I,J) = SUM
      SGAIN(I,J) = SUM
  440 CONTINUE
  441 CONTINUE
C
cxx      DO 450  I=1,M
      DO 451  I=1,M
      XS(I) = XFS(I,II)
      DO 450  J=1,M
      WRK(I,J) = 0.0D0
cxx  450 VS (I,J) = VFS(I,J,II)
      VS (I,J) = VFS(I,J,II)
  450 CONTINUE
  451 CONTINUE
C
cxx      DO 460  J=1,M
      DO 461  J=1,M
      DO 460  I=1,M
cxx  460 XS(I) = XS(I) + SGAIN(I,J)*(XSS(J,II+1) - XPS(J,II+1))
      XS(I) = XS(I) + SGAIN(I,J)*(XSS(J,II+1) - XPS(J,II+1))
  460 CONTINUE
  461 CONTINUE
C
cxx      DO 470  J=1,M
cxx      DO 470 IJ=1,M
      DO 472  J=1,M
      DO 471 IJ=1,M
      DO 470  I=1,M
cxx  470 WRK(I,J)=WRK(I,J) + SGAIN(I,IJ)*(VSS(IJ,J,II+1)-VPS(IJ,J,II+1))
      WRK(I,J)=WRK(I,J) + SGAIN(I,IJ)*(VSS(IJ,J,II+1)-VPS(IJ,J,II+1))
  470 CONTINUE
  471 CONTINUE
  472 CONTINUE
C
cxx      DO 480  J=1,M
cxx      DO 480 IJ=1,M
      DO 482  J=1,M
      DO 481 IJ=1,M
      DO 480  I=1,M
cxx  480 VS(I,J) = VS(I,J) + WRK(I,IJ)*SGAIN(J,IJ)
      VS(I,J) = VS(I,J) + WRK(I,IJ)*SGAIN(J,IJ)
  480 CONTINUE
  481 CONTINUE
  482 CONTINUE
      DO 485 I=1,M
cxx  485 IF( VS(I,I).LT.0.0D0 )  VS(I,I) = 0.0D0
      IF( VS(I,I).LT.0.0D0 )  VS(I,I) = 0.0D0
  485 CONTINUE
C
cxx      DO 490  I=1,M
      DO 491  I=1,M
      XSS(I,II) = XS(I)
      DO 490  J=1,M
cxx  490 VSS(I,J,II) = VS(I,J)
      VSS(I,J,II) = VS(I,J)
  490 CONTINUE
  491 CONTINUE
      END IF
C
  500 CONTINUE
C
      RETURN
      E N D
cc      SUBROUTINE  GINVRS( A,DET,M,MJ )
      SUBROUTINE  GINVRS( A,DET,M )
C
C  ...  Generalized inverse of a square matrix A  ...
C
C     Inputs:
C        A:     M*M matrix
C        M:     Dimension of A
C        MJ:    Adjustable dimension of A
C     Outputs:
C        A:     Generalize inverse of A
C        DET:   Determinant of A
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  A(MJ,MJ), IND(50)
cxx      DIMENSION  A(M,M), IND(M+1)
      INTEGER :: M, IND(M+1)
      REAL(8) :: A(M,M), DET, EPS, AMAX, SUM
C
      EPS = 1.0D-10
C
      DO 10  I=1,M
cxx   10 IND(I) = I
       IND(I) = I
   10 CONTINUE
cc------------
      I0 = 0
      LMAX = 0
cc------------
      DO 60  L=1,M
      AMAX = 0.0D0
      DO 20  I=L,M
      IF( A(IND(I),IND(I)).GT.AMAX )  THEN
        AMAX = A(IND(I),IND(I))
        I0 = I
      END IF
   20 CONTINUE
      IF( AMAX.GT.EPS*A(IND(1),IND(1)) )  THEN
         IMAX = IND(I0)
         DO 30  I=I0,L+1,-1
cxx   30    IND(I) = IND(I-1)
         IND(I) = IND(I-1)
   30    CONTINUE
         IND(L) = IMAX
         LMAX   = L
cxx         DO 40  I=L+1,M
         DO 41  I=L+1,M
         A(IND(I),IMAX) = -A(IND(I),IMAX)/A(IMAX,IMAX)
         DO 40  J=L+1,M
cxx   40    A(IND(I),IND(J)) = A(IND(I),IND(J))
         A(IND(I),IND(J)) = A(IND(I),IND(J))
     *                    + A(IND(I),IMAX)*A(IMAX,IND(J))
   40    CONTINUE
   41    CONTINUE
      ELSE
cxx         DO 50  I=L,M
         DO 51  I=L,M
         DO 50  J=L,M
cxx   50    A(IND(I),IND(J)) = 0.0D0
         A(IND(I),IND(J)) = 0.0D0
   50    CONTINUE
   51    CONTINUE
         GO TO 70
      END IF
   60 CONTINUE
C
   70 DET = 1.0D0
      DO 80  I=1,M
C     DET = DET*A(IND(I),IND(I))
C     IF( A(IND(I),IND(I)).GT.EPS*AMAX )  THEN
      IF( A(IND(I),IND(I)).GT.0.0D0 )  THEN
        A(IND(I),IND(I)) = 1.0D0/A(IND(I),IND(I))
      ELSE
        A(IND(I),IND(I)) = 0.0D0
      END IF
   80 CONTINUE
C
      MS = MIN0( M-1,LMAX )
      DO 200  L=MS,1,-1
      DO 100 J=L+1,M
      SUM = 0.0D0
      DO 90  I=L+1,M
cxx   90 SUM = SUM + A(IND(I),IND(L))*A(IND(I),IND(J))
      SUM = SUM + A(IND(I),IND(L))*A(IND(I),IND(J))
   90 CONTINUE
cxx  100 A(IND(L),IND(J)) = SUM
      A(IND(L),IND(J)) = SUM
  100 CONTINUE
      SUM = A(IND(L),IND(L))
      DO 110  I=L+1,M
cxx  110 SUM = SUM + A(IND(L),IND(I))*A(IND(I),IND(L))
      SUM = SUM + A(IND(L),IND(I))*A(IND(I),IND(L))
  110 CONTINUE
      A(IND(L),IND(L)) = SUM
      DO 120  I=L+1,M
cxx  120 A(IND(I),IND(L)) = A(IND(L),IND(I))
      A(IND(I),IND(L)) = A(IND(L),IND(I))
  120 CONTINUE
      IMAX = IND(L)
      DO 130 I=L+1,M
      IF( IMAX.GT.IND(I) )  THEN
         IND(I-1) = IND(I)
         IND(I)   = IMAX
      END IF
  130 CONTINUE
  200 CONTINUE
C
      RETURN
      E N D

cc      SUBROUTINE  CRSCOR( Y,N,ID,LAG,MJ,MJ1,OUTMIN,OUTMAX,C,R,YMEAN )
       SUBROUTINE  CRSCOR( Y,N,ID,LAG,OUTMIN,OUTMAX,C,R,YMEAN )
C
C ... cross correlation function computation ...
C
C     Inputs:
C        Y(I,J):   multi-variate time series
C        N:        data length
C        ID:       dimension of the observation
C        LAG:      maximum lag of cross-covariance
C        MJ:       adjustable dimension (MJ.GE.N)
C        MJ1:      adjustable dimension (MJ1.GE.ID)
C        OUTMIN:   bound for outliers in low side
C        OUTMAX:   bound for outliers in high side
C     Outputs:
C        C:        sample cross-covariance function
C        R:        sample cross-correlation function
C        YMEAN:    sample mean vector
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  Y(MJ,MJ1), OUTMIN(MJ1), OUTMAX(MJ1)
cc      DIMENSION  C(0:LAG,MJ1,MJ1), R(0:LAG,MJ1,MJ1)
cc      DIMENSION  YMEAN(MJ1), NSUM(10)
cxx      DIMENSION  Y(N,ID), OUTMIN(ID), OUTMAX(ID)
cxx      DIMENSION  C(0:LAG,ID,ID), R(0:LAG,ID,ID)
cxx      DIMENSION  YMEAN(ID), NSUM(ID)
      INTEGER :: N, ID, LAG, NSUM(ID)
      REAL(8) :: Y(N,ID), OUTMIN(ID), OUTMAX(ID), C(0:LAG,ID,ID),
     1           R(0:LAG,ID,ID), YMEAN(ID), SUM
C
      DO 10 J=1,ID
cxx   10 CALL  MEAN( Y(1,J),N,OUTMIN(J),OUTMAX(J),NSUM(J),YMEAN(J) )
      CALL  MEAN( Y(1,J),N,OUTMIN(J),OUTMAX(J),NSUM(J),YMEAN(J) )
   10 CONTINUE
C
cxx      DO 30 I=1,ID
cxx      DO 30 J=1,ID
      DO 32 I=1,ID
      DO 31 J=1,ID
cc      WRITE(6,*)  I,J
      DO 30 L=0,LAG
      SUM = 0.0D0
      NNN = 0
      DO 20 II=L+1,N
      IF( Y(II,I).GT.OUTMIN(I).AND.Y(II,I).LT.OUTMAX(I) )  THEN
        IF( Y(II-L,J).GT.OUTMIN(J).AND.Y(II-L,J).LT.OUTMAX(J) )  THEN
          SUM = SUM + (Y(II,I)-YMEAN(I))*(Y(II-L,J)-YMEAN(J))
          NNN = NNN + 1
        END IF
      END IF
   20 CONTINUE
cxx   30 C(L,I,J)=SUM/DSQRT( DBLE( NSUM(I)*NSUM(J) ) )
      C(L,I,J)=SUM/DSQRT( DBLE( NSUM(I)*NSUM(J) ) )
   30 CONTINUE
   31 CONTINUE
   32 CONTINUE
C
cxx      DO 40 I=1,ID
cxx      DO 40 J=1,ID
      DO 42 I=1,ID
      DO 41 J=1,ID
      DO 40 L=0,LAG
cxx   40 R(L,I,J) = C(L,I,J)/DSQRT(C(0,I,I)*C(0,J,J))
      R(L,I,J) = C(L,I,J)/DSQRT(C(0,I,I)*C(0,J,J))
   40 CONTINUE
   41 CONTINUE
   42 CONTINUE
      RETURN
C
      E N D

      SUBROUTINE  MEAN( Y,N,OUTMIN,OUTMAX,NSUM,YMEAN )
C
C  ...  This subroutine computes sample mean  ...
C
C     Inputs:
C        Y(I):    time series
C        N:       data length
C        OUTMIN:  bound for outliers in low side
C        OUTMAX:  bound for outliers in high side
C     Outputs:
C        NSUM:    number of non-outlier observations
C        YMEAN:   sample mean
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION Y(N)
      INTEGER :: N, NSUM
      REAL(8) :: Y(N), OUTMIN, OUTMAX, YMEAN, YSUM
C
      NSUM = 0
      YSUM = 0.0D0
      DO 10 I=1,N
      IF( Y(I) .GT.OUTMIN .AND. Y(I).LT.OUTMAX ) THEN
         NSUM = NSUM + 1
         YSUM = YSUM + Y(I)
      END IF
   10 CONTINUE
      YMEAN = YSUM/NSUM
C
      RETURN
      E N D
      SUBROUTINE AUTCOR( C,MAXLAG,R )
C
C  ...  This subroutine compute sample autocorrelation  ...
C
C     Inputs:
C        C(I):    sample autocovariance
C        MAXLAG:  maximum lag of autocovariance
C     Output:
C        R(I):    sample autocorrelation
C
cxx      IMPLICIT REAL*8 (A-H,O-Z)
cxx      DIMENSION  C(0:MAXLAG), R(0:MAXLAG)
      INTEGER :: MAXLAG
      REAL(8) :: C(0:MAXLAG), R(0:MAXLAG)
C
      DO 10 LAG=0,MAXLAG
cxx   10 R(LAG) = C(LAG)/C(0)
      R(LAG) = C(LAG)/C(0)
   10 CONTINUE
      RETURN
      E N D
cc      SUBROUTINE AUTCOV( Y,N,MAXLAG,OUTMIN,OUTMAX,C )
      SUBROUTINE AUTCOV( Y,N,MAXLAG,OUTMIN,OUTMAX,C,YMEAN )
C
C  ...  This subroutine computes sample autocovariance  ...
C
C     Inputs:
C        Y(I):    time series
C        N:       data length
C        MAXLAG:  maximum lag of autocovariance
C        OUTMIN:  bound for outliers in low side
C        OUTMAX:  bound for outliers in high side
C     OUTPUT:
C        C(I):    autocovariance function
C
cxx      IMPLICIT REAL*8( A-H,O-Z )
cxx      DIMENSION Y(N), C(0:MAXLAG )
      INTEGER :: N, MAXLAG
      REAL(8) :: Y(N), OUTMIN, OUTMAX, C(0:MAXLAG ), YMEAN, SUM
C
C  ...  sample mean  ...
C
      CALL  MEAN( Y,N,OUTMIN,OUTMAX,NSUM,YMEAN )
C
C  ...  sample autocovariance  ...
C
      DO 20 LAG = 0,MAXLAG
      SUM = 0.0D0
      DO 10 I=LAG+1,N
      IF( Y(I).GT.OUTMIN .AND. Y(I).LT.OUTMAX )   THEN
         IF( Y(I-LAG).GT.OUTMIN .AND. Y(I-LAG).LT.OUTMAX ) THEN
         SUM = SUM + ( Y(I)-YMEAN)*( Y(I-LAG) - YMEAN )
         END IF
      END IF
   10 CONTINUE
cxx   20 C(LAG) = SUM/NSUM
      C(LAG) = SUM/NSUM
   20 CONTINUE
C
      RETURN
      E N D
cc      SUBROUTINE  UNICOR( Y,N,MAXLAG,OUTMIN,OUTMAX,COV )
      SUBROUTINE  UNICOR( Y,N,MAXLAG,OUTMIN,OUTMAX,COV,YMEAN )
C
C  ...  This subroutine computes sample autocovariance function,
C       sample autocorrelationfunction and their error bounds  ...
C
C     Inputs:
C        Y(I):    time series
C        N:       data length
C        MAXLAG:  maximum lag of autocovariance
C        OUTMIN:  bound for outliers in low side
C        OUTMAX:  bound for outliers in high side
C     Output:
C        COV(I,J):  I=0,...,MAXLAG
C                   J = 1   sample autocovariance
C                     = 2   sample autocorrelation
C                     = 3   error bound for sample autocovariance
C                     = 4   error bound for sample autocorrelation
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  Y(N), COV(0:MAXLAG,4)
      INTEGER :: N, MAXLAG
      REAL(8) :: Y(N), OUTMIN, OUTMAX, COV(0:MAXLAG,4), YMEAN
C
C  ...  sample autocovariance function  ...
C
cc      CALL  AUTCOV( Y,N,MAXLAG,OUTMIN,OUTMAX,COV )
      CALL  AUTCOV( Y,N,MAXLAG,OUTMIN,OUTMAX,COV,YMEAN )
C
C  ...  sample autocorrelation function  ...
C
      CALL  AUTCOR( COV,MAXLAG,COV(0,2) )
C
C  ...  error bounds  ...
C
      CALL  ERRACF( COV,N,MAXLAG,COV(0,3),COV(0,4) )
C
      RETURN
      E N D
      SUBROUTINE ERRACF( C,N,MAXLAG,CERR,RERR )
C
C  ...  This subroutine computes error bounds for sample autocovariance
C       and autocorrelation function  ...
C
C     Inputs:
C        C(I):     sample autocovariance
C        N:        data length
C        MAXLAG:   maximum lag of autocovariance
C     Outputs:
C        CERR(I):  error bound for I-th autocovariance
C        RERR(I):  error bound for I-th autocorrelation
C
cxx      IMPLICIT REAL*8 (A-H,O-Z)
cxx      DIMENSION  C(0:MAXLAG), CERR(0:MAXLAG), RERR(0:MAXLAG)
      INTEGER :: N, MAXLAG
      REAL(8) :: C(0:MAXLAG), CERR(0:MAXLAG), RERR(0:MAXLAG), SUM
C
      SUM = C(0)**2
      CERR(0)   = DSQRT( 2*SUM/N )
      DO 10 LAG=1,MAXLAG
      IF(LAG.GE.2)  SUM = SUM + 2*C(LAG-1)**2
cxx   10 CERR(LAG) = DSQRT( SUM/N )
      CERR(LAG) = DSQRT( SUM/N )
   10 CONTINUE
C
      RERR(0) = 0.0D0
      DO 20 LAG=1,MAXLAG
cxx   20 RERR(LAG) = CERR(LAG)/C(0)
      RERR(LAG) = CERR(LAG)/C(0)
   20 CONTINUE
      RETURN
      E N D
C
C
cc      SUBROUTINE  WINDOW( PE,NP,IWINDW,SPE )
      SUBROUTINE  WINDOW( PE,NP,IWINDW,SPE,IFG )
C
C  ...  Smoothing by spectral window and log-transformation  ...
C
C     Inputs:
C        PE(I):    raw specrum
C        NP:       number of frequencies
C        IWINDW:   window type (0: box-car, 1: Hanning, 2: Hamming)
C     Outputs:
C        SPE(I):   logarithm of smoothed periodogram
C
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cxx      DIMENSION  PE(0:NP), SPE(0:NP)
cxx      DIMENSION  W(0:1,2)
      INTEGER :: NP, IWINDW, IFG
      REAL(8) ::PE(0:NP), SPE(0:NP), W(0:1,2), PMIN
      DATA  W/0.5D0, 0.25D0, 0.54D0, 0.23D0/
C
      IF( IWINDW.EQ.0 )  THEN
        PMIN = 1.0D30
        DO 10 I=0,NP
cxx   10   IF( PE(I).GT.0.0D0 .AND. PE(I).LT.PMIN )  PMIN = PE(I)
        IF( PE(I).GT.0.0D0 .AND. PE(I).LT.PMIN )  PMIN = PE(I)
   10   CONTINUE
        DO 20 I=0,NP
cxx   20   SPE(I) = DMAX1( PE(I),PMIN )
        SPE(I) = DMAX1( PE(I),PMIN )
   20   CONTINUE
      ELSE
        SPE(0) = W(0,IWINDW)*PE(0) + W(1,IWINDW)*PE(1)*2
        SPE(NP)= W(0,IWINDW)*PE(NP)+ W(1,IWINDW)*PE(NP-1)*2
        DO 30  I=1,NP-1
cxx   30   SPE(I) = W(0,IWINDW)*PE(I) + W(1,IWINDW)*(PE(I-1) + PE(I+1))
        SPE(I) = W(0,IWINDW)*PE(I) + W(1,IWINDW)*(PE(I-1) + PE(I+1))
   30   CONTINUE
      END IF
c---------- 2013/07/03
      IFG = 0
      DO 40 I=0,NP
cxx   40 IF( SPE(I) .LE. 0 ) IFG = -1
      IF( SPE(I) .LE. 0 ) IFG = -1
   40 CONTINUE
      IF( IFG .NE. 0 ) GO TO 60
c----------
      DO 50 I=0,NP
c   50 SPE(I) = DLOG10( SPE(I) )
c      DO 50 I=0,NP
c		 if (SPE(I) .gt. 0) then
           SPE(I) = DLOG10( SPE(I) )
c		 else
c		   SPE(I) = -10
c         end if
   50 CONTINUE
C
   60 CONTINUE
      RETURN
      E N D
C
C
      SUBROUTINE  SETSEA( M1,M2,M3,IPER,AR,TAU1,TAU2,TAU3,SIG2,
cc     *                    F,G,H,Q,R,M,K,MJ )
     *                    F,G,H,Q,R,M,K )
C
C  ...  SET STATE SPACE MODEL FOR SEASONAL ADJUSTMENT  ...
C
C     INPUTS:
C        M1,M2,M3:   MODEL ORDERS
C        IPER:       NUMBER OF SEASONS IN ONE PERIOD
C        AR(I):      AR COEFFICIENTS
C        TAU1-TAU3:  SYSTEM NOISE VARIANCES
C        MJ:         ADJUSTABLE DIMENSION
C     OUTPUTS:
C        F,G,H:      MATRICES OF STATE SPACE MODEL
C        Q,R:        SYSTEM NOISE AND OBS NOISE VARIANCES
C        M:          STATE DIMENSION
C        K:          DIMENSION OF SYSTEM NOISE
C     MODIFIED  2/16/93
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  F(MJ,MJ), G(MJ,MJ), H(MJ), Q(MJ,MJ), AR(*), R(1,1)
cxx      DIMENSION  F(M,M), G(M,K), H(M), Q(K,K), AR(M3), R(1,1)
      INTEGER :: M1, M2, M3, IPER, M, K 
      REAL(8) :: AR(M3), TAU1, TAU2, TAU3, SIG2, F(M,M), G(M,K), H(M),
     1           Q(K,K), R(1,1)
C
cc      M = M1 + M2*(IPER-1) + M3
cc      K = ID(M1) + ID(M2) + ID(M3)
cxx      DO 10 I=1,M
cxx      H(I) = 0.0D0
cxx      DO 10 J=1,M
cxx   10 F(I,J) = 0.0D0
cxx      DO 20 I=1,M
cxx      DO 20 J=1,K
cxx   20 G(I,J) = 0.0D0
cxx      DO 30 I=1,K
cc   30 Q(I,I) = 0.0D0
cxx      DO 30 J=1,K
cxx   30 Q(I,J) = 0.0D0
      H(1:M) = 0.0D0
      F(1:M,1:M) = 0.0D0
      G(1:M,1:K) = 0.0D0
      Q(1:K,1:K) = 0.0D0
C
      IF( M1.GT.0 )  THEN
        IF( M1.EQ.1 )  F(1,1) = 1.0D0
        IF( M1.EQ.2 )  THEN
          F(1,1) = 2.0D0
          F(1,2) =-1.0D0
          F(2,1) = 1.0D0
        END IF
        G(1,1) = 1.0D0
        H(1)   = 1.0D0
        Q(1,1) = TAU1
      END IF
C
C  ...  SEASONAL COMPONENT  ...
C
      IF( M2.GT.0 )  THEN
        L1 = ID(M1) + 1
        DO 40 I=1,IPER-1
cxx   40   F(M1+1,M1+I) = -1.0D0
        F(M1+1,M1+I) = -1.0D0
   40   CONTINUE
        DO 50 I=2,IPER-1
cxx   50   F(M1+I,M1+I-1) = 1.0D0
        F(M1+I,M1+I-1) = 1.0D0
   50   CONTINUE
        G(M1+1,L1) = 1.0D0
        H(M1+1)    = 1.0D0
        Q(L1,L1)   = TAU2
      END IF
C
C  ...  AR COMPONENT  ...
C
      IF( M3.GT.0 )  THEN
        M12 = M1 + M2*(IPER-1)
        L12 = ID(M1) + ID(M2) + 1
        DO 60 I=1,M3
cxx   60   F(M12+1,M12+I) = AR(I)
        F(M12+1,M12+I) = AR(I)
   60   CONTINUE
        DO 70 I=2,M3
cxx   70   F(M12+I,M12+I-1) = 1.0D0
        F(M12+I,M12+I-1) = 1.0D0
   70   CONTINUE
        G(M12+1,L12) = 1.0D0
        H(M12+1)     = 1.0D0
        Q(L12,L12)   = TAU3
      END IF
C
      R(1,1) = SIG2
C
      RETURN
      E N D
      SUBROUTINE  CHOLES ( X,MJ,N,Y,NJ )
C
C  ...  CHOLESKY DECOMPOSITION  ...
C
C     INPUTS:
C        X(I,J):   SYMMETRIC MATRIX
C        N:        DIMENSION OF Z
C        MJ:       ADJUSTABLE DIMENSION OF X
C        NJ:       ADJUSTABLE DIMENSION OF Y
C     OUTPUT:
C        Y(I,J):   LOWER TRIANGULAR MATRIX; X = Y*Y'
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  X(MJ,MJ), Y(NJ,NJ)
      INTEGER :: MJ, N, NJ
      REAL(8) :: X(MJ,MJ), Y(NJ,NJ), SUM1, SUM2
C
cxx      DO 10 I=1,N
cxx      DO 10 J=1,N
cxx   10 Y(I,J) = 0.0D0
      Y(1:N,1:N) = 0.0D0
C
      DO 100 J = 1,N
      SUM1 = X(J,J)
      DO 20 K=1,J-1
cxx   20 SUM1 = SUM1 - Y(J,K)**2
      SUM1 = SUM1 - Y(J,K)**2
   20 CONTINUE
      IF( SUM1.GT.0.0D0 ) Y(J,J) = DSQRT( SUM1 )
      IF( SUM1.EQ.0.0D0 ) Y(J,J) = 0.0D0
      DO 40 I=J+1,N
      SUM2 = 0.0D0
      DO 30 K=1,J-1
cxx   30 SUM2 = SUM2 + Y(I,K)*Y(J,K)
      SUM2 = SUM2 + Y(I,K)*Y(J,K)
   30 CONTINUE
cxx   40 Y(I,J) = ( X(I,J) - SUM2 ) / Y(J,J)
      Y(I,J) = ( X(I,J) - SUM2 ) / Y(J,J)
   40 CONTINUE
  100 CONTINUE
C
      RETURN
      E N D
cc      SUBROUTINE ARMCOV( M,L,A,B,SIG2,K,COV )
      SUBROUTINE ARMCOV( M,L,A,B,SIG2,K,COV,KMAX,IER )
C
C ...  Autocovariance Function of ARMA model  ...
C
C     Inputs:
C        M:     AR order
C        L:     MA order
C        A(I):  AR coefficient
C        B(I):  MA coefficient
C        SIG2:  innovation variance
C        K:     Required maximum lag of autocovariance
C     Output:
C        COV(I):  Autocovariance function
C     Y.I.
cxx      IMPLICIT REAL*8( A-H,O-Z )
cc      DIMENSION  A(*), B(*), COV(0:K), G(0:100), X(30,30)
cc      DIMENSION  Z(100), UL(30,30), IPS(100)
cxx      DIMENSION  A(M), B(L), COV(0:K), G(0:KMAX), X(M+1,M+1)
cxx      DIMENSION  Z(M+1), UL(M+1,M+1), IPS(M+1)
      INTEGER :: M, L, K, KMAX, IER, IPS(M+1)
      REAL(8) :: A(M), B(L), SIG2, COV(0:K)
      REAL(8) :: G(0:KMAX), X(M+1,M+1), Z(M+1), UL(M+1,M+1), SUM
C
cc      KMAX = MAX(M,L,K)
      CALL  IMPULS( M,L,A,B,KMAX,G )
C
cxx      DO 10 I=1,M+1
cxx      DO 10 J=1,M+1
cxx   10 X(I,J) = 0.0D0
      X(1:M+1,1:M+1) = 0.0D0
      DO 20 I=1,M+1
cxx   20 X(I,I) = 1.0D0
      X(I,I) = 1.0D0
   20 CONTINUE
c
      IF( M.GT.0 ) THEN

cxx      DO 30 I=1,M
      DO 31 I=1,M
      DO 30 J=2,M-I+2
cxx   30 X(I,J) = X(I,J) - A(I+J-2)
      X(I,J) = X(I,J) - A(I+J-2)
   30 CONTINUE
   31 CONTINUE
cxx      DO 40 I=2,M+1
      DO 41 I=2,M+1
      DO 40 J=1,I-1
cxx   40 X(I,J) = X(I,J) - A(I-J)
      X(I,J) = X(I,J) - A(I-J)
   40 CONTINUE
   41 CONTINUE

      END IF
c
C
cc      CALL  DECOM( M+1,X,30,UL,IPS )
      CALL  DECOM( M+1,X,UL,IPS,IER )
      if( ier.ne.0 ) return
C
      SUM = 1.0D0
c
      IF( L.GT.0 ) THEN

      DO 50 J=1,L
cxx   50 SUM = SUM - B(J)*G(J)
      SUM = SUM - B(J)*G(J)
   50 CONTINUE

      END IF
c
      Z(1)= SIG2*SUM
c
      IF( M.GT.0 ) THEN

      DO 70 I=2,M+1
      SUM = 0.0D0
      DO 60 J=I-1,L
cc   60 SUM = SUM - B(J)*G(J-I+1)
cxx   60 IF( L.GT. 0) SUM = SUM - B(J)*G(J-I+1)
      IF( L.GT. 0) SUM = SUM - B(J)*G(J-I+1)
   60 CONTINUE
cxx   70 Z(I) = SIG2*SUM
      Z(I) = SIG2*SUM
   70 CONTINUE

      END IF
c
C
cc      CALL  SOLVE( M+1,UL,30,Z,COV,IPS)
      CALL  SOLVE( M+1,UL,Z,COV,IPS )
C
      DO 100 J=M+1,K
      SUM = 0.0D0
      DO 80 I=1,M
cc   80 SUM = SUM + A(I)*COV(J-I)
cxx   80 IF( M.GT.0 ) SUM = SUM + A(I)*COV(J-I)
      IF( M.GT.0 ) SUM = SUM + A(I)*COV(J-I)
   80 CONTINUE
      DO 90 I=J,L
cc   90 SUM = SUM - B(I)*G(I-J)*SIG2
cxx   90 IF( L.GT.0 ) SUM = SUM - B(I)*G(I-J)*SIG2
      IF( L.GT.0 ) SUM = SUM - B(I)*G(I-J)*SIG2
   90 CONTINUE
cxx  100 COV(J) = SUM
      COV(J) = SUM
  100 CONTINUE
C
      RETURN
      E N D
C
cc      SUBROUTINE DECOM( N,A,MJ,UL,IPS )
      SUBROUTINE DECOM( N,A,UL,IPS,IER )
C
C  ...  UL decomposition:  A = L*U  ...
C
C     Inputs:
C        N:      Dimension of the matrix A
C        A(I,J): N*N positive definite matrix
C        MJ:     Adjustable dimension of A and UL
C     Outputs:
C        UL(I,J):  L-I and U
C        IPS:    Index vector
C     Y.I.
C        IER:    Error code
C
cxx      IMPLICIT REAL*8(A-H,O-Z )
cc      DIMENSION A(MJ,*),UL(MJ,*),SCALES(100),IPS(100)
cxx      DIMENSION A(N,N),UL(N,N),SCALES(N),IPS(N)
      INTEGER :: N, IPS(N), IER
      REAL(8) :: A(N,N), UL(N,N), SCALES(N), RNORM, BIG, SIZE, PIVOT, TM
C
      IER = 0
C
      DO 20 I=1,N
      IPS(I) = I
      RNORM = 0.0D0
      DO 10 J=1,N
      UL(I,J) = A(I,J)
      IF(RNORM.LT.ABS(UL(I,J)) ) RNORM = ABS( UL(I,J) )
   10 CONTINUE
      IF( RNORM .NE. 0.0D0 )  THEN
          SCALES(I) = 1/RNORM
      ELSE
          SCALES(I) = 0.0D0
cc          CALL  SING(0)         
          IER = 1
      END IF
   20 CONTINUE
      if( ier.ne.0 ) return
C
cc-------------
      INDEX = 0
cc-------------
      DO 60 K=1,N-1
      BIG = 0.0D0
      DO 30 I=K,N
      SIZE = ABS( UL(IPS(I),K) )*SCALES( IPS(I) )
      IF( BIG.LT.SIZE ) THEN
          BIG = SIZE
          INDEX = I
      END IF
   30 CONTINUE
      IF( BIG.EQ. 0.0D0 )  THEN
cc          CALL  SING(1)
          IER = 2
          GO TO 60
      END IF
      IF( INDEX.NE.K ) THEN
      J = IPS(K)
      IPS(K) = IPS(INDEX)
      IPS(INDEX) = J
      END IF
C
      PIVOT = UL(IPS(K),K)
      DO 50 I=K+1,N
      TM = UL( IPS(I),K)/PIVOT
      UL( IPS(I),K) = TM
      IF( TM.NE. 0.0D0 )  THEN
      DO 40 J = K+1,N
cxx   40 UL( IPS(I),J ) = UL( IPS(I),J)-TM*UL( IPS(K),J)
      UL( IPS(I),J ) = UL( IPS(I),J)-TM*UL( IPS(K),J)
   40 CONTINUE
C     WRITE(6,*) (UL(IPS(I),J),J=1,N)
      END IF
   50 CONTINUE
   60 CONTINUE
      if( ier.ne.0 ) return
C
cc      IF( UL(IPS(N),N) .EQ. 0.0D0 )   CALL  SING(2)
      IF( UL(IPS(N),N) .EQ. 0.0D0 )   IER = 3
      RETURN
      E N D
C
      SUBROUTINE IMPULS( M,L,A,B,K,G )
C
C ...  Impulse Response Function  ...
C
C     Inputs:
C        M:     AR order
C        L:     MA order
C        A(I):  AR coefficient
C        B(I):  MA coefficient
C        K:     Required maximum lag of impulse respose
C     Output:
C        G(I):  Impulse response function
C     Y.I.
cxx      IMPLICIT REAL*8( A-H,O-Z )
cc      DIMENSION A(*), B(*), G(0:K)
cxx      DIMENSION A(M), B(L), G(0:K)
      INTEGER :: M, L, K
      REAL(8) :: A(M), B(L), G(0:K), SUM
C
      G(0) = 1.0
      DO  20 I=1,K
      SUM = 0.0D0
      IF(I.LE.L) SUM = -B(I)
      DO  10 J=1,I
cxx   10 IF(J.LE.M) SUM = SUM + A(J)*G(I-J)
      IF(J.LE.M) SUM = SUM + A(J)*G(I-J)
   10 CONTINUE
cxx   20 G(I) = SUM
      G(I) = SUM
   20 CONTINUE
C
      RETURN
      E N D
C
cc      SUBROUTINE  SOLVE( N,UL,MJ,B,X,IPS )
      SUBROUTINE  SOLVE( N,UL,B,X,IPS )
C
C  ...  Solve Ax=b using UL obtained by DECOM  ...
C
C     Inputs:
C        N:     Dimension of UL and B
C        UL:    LU decomposition of A
C        MJ:    Adjustable dimension of A
C        B:
C        IPS:   index vector
C     Output:
C        X:     Solution
C     Y.I.
cxx      IMPLICIT REAL*8( A-H,O-Z )
cc      DIMENSION UL(MJ,*),B(*),X(*),IPS(100)
cxx      DIMENSION UL(N,N),B(N),X(N),IPS(N)
      INTEGER :: N, IPS(N)
      REAL(8) :: UL(N,N), B(N),X(N), SUM
C
      DO 20 I=1,N
      SUM = 0.0D0
      DO 10 J=1,I-1
cxx   10 SUM = SUM + UL(IPS(I),J)*X(J)
      SUM = SUM + UL(IPS(I),J)*X(J)
   10 CONTINUE
cxx   20 X(I) = B(IPS(I)) - SUM
      X(I) = B(IPS(I)) - SUM
   20 CONTINUE
C
      DO 40 I=N,1,-1
      SUM = 0.0D0
      DO 30 J=I+1,N
cxx   30 SUM = SUM + UL(IPS(I),J)*X(J)
      SUM = SUM + UL(IPS(I),J)*X(J)
   30 CONTINUE
cxx   40 X(I) = ( X(I)-SUM )/UL(IPS(I),I)
      X(I) = ( X(I)-SUM )/UL(IPS(I),I)
   40 CONTINUE
      RETURN
cxx  600 FORMAT(1H ,'N=',I10,/,(5X,'IPS=',I10 ) )
      E N D
C
C
cc      SUBROUTINE  SMOTH1( A,M,MMAX,NC,NS,N,NE,MJ,
      SUBROUTINE  SMOTH1( A,M,MMAX,NC,NS,N,NE,NMAX,MJ,
     *                    VFS,VPS,VSS,XFS,XPS,XSS )
C
C  ...  Fixed Interval Smoother (General Form)  ...
C
C     Inputs:
C        A:      Parameter of matrix F
C        M:      Dimension of the state vector
C        MMAX:   Adjustable dimension of A
C        NC:     Number of components
C        N:      Data length
C        NS:     Start position of filtering
C        NE:     End position of filtering
C        MJ:     Adjustable dimension of XF, VF
C        VFS:    Covariance matrices of the filter
C        VPS:    Covariance matrices of the predictor
C        XFS:    Mean vectors of the filter
C        XPS:    Mean vectors of the predictor
C     Outputs:
C        VSS:    Covariance matrices of the smoother
C        XSS:    Mean vectors of the smoother
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  A(MMAX,NC), M(NC)
cc      DIMENSION  XS(40), VS(40,40), VP(40,40)
cc      DIMENSION  XFS(MJ,N), XPS(MJ,N), XSS(MJ,N)
cc      DIMENSION  VFS(MJ,MJ,N), VPS(MJ,MJ,N), VSS(MJ,MJ,N)
cc      DIMENSION  WRK(40,40), SGAIN(40,40)
cc      DIMENSION  I0(10)
cxx      DIMENSION  XS(MJ), VS(MJ,MJ), VP(MJ,MJ)
cxx      DIMENSION  XFS(MJ,NMAX), XPS(MJ,NMAX), XSS(MJ,NMAX)
cxx      DIMENSION  VFS(MJ,MJ,NMAX), VPS(MJ,MJ,NMAX), VSS(MJ,MJ,NMAX)
cxx      DIMENSION  WRK(MJ,MJ), SGAIN(MJ,MJ)
cxx      DIMENSION  I0(NC)
      INTEGER :: MMAX, NC, N S, N, NE, NMAX, MJ, M(NC), I0(NC)
      REAL(8) :: A(MMAX,NC), VFS(MJ,MJ,NMAX), VPS(MJ,MJ,NMAX),
     1           VSS(MJ,MJ,NMAX), XFS(MJ,NMAX), XPS(MJ,NMAX),
     2           XSS(MJ,NMAX)
      REAL(8) :: XS(MJ), VS(MJ,MJ), VP(MJ,MJ), WRK(MJ,MJ), SGAIN(MJ,MJ),
     1           SUM, VDET
C
      I0(1) = 0
      DO 10 I=2,NC
cxx   10 I0(I) = I0(I-1) + M(I-1)
      I0(I) = I0(I-1) + M(I-1)
   10 CONTINUE
      MM = I0(NC) + M(NC)
C
C  ...  SMOOTHING  ...
C
      NSS = MIN0( N,NE )
cxx      DO 20  I=1,MM
      DO 21  I=1,MM
      XS(I)      = XFS(I,NSS)
      XSS(I,NSS) = XFS(I,NSS)
      DO 20  J=1,MM
      VS(I,J)      = VFS(I,J,NSS)
cxx   20 VSS(I,J,NSS) = VFS(I,J,NSS)
      VSS(I,J,NSS) = VFS(I,J,NSS)
   20 CONTINUE
   21 CONTINUE
cxx      DO 30 II=NSS+1,NE
cxx      DO 30 I=1,MM
      DO 32 II=NSS+1,NE
      DO 31 I=1,MM
      XSS(I,II) = XFS(I,II)
      DO 30 J=1,MM
cxx   30 VSS(I,J,II) = VFS(I,J,II)
      VSS(I,J,II) = VFS(I,J,II)
   30 CONTINUE
   31 CONTINUE
   32 CONTINUE
C
      DO 500  II=NSS-1,NS,-1
C
      NZERO = 0
      DO 100 I=1,MM
cxx  100 IF( VFS(I,I,II).GT.1.0D-12 )  NZERO = NZERO + 1
      IF( VFS(I,I,II).GT.1.0D-12 )  NZERO = NZERO + 1
  100 CONTINUE
C
      IF( NZERO.EQ.0 )  THEN
cxx         DO 110  I=1,MM
         DO 111  I=1,MM
         XS(I)     = XFS(I,II)
         XSS(I,II) = XFS(I,II)
         DO 110  J=1,MM
         VS(I,J)     = VFS(I,J,II)
cxx  110    VSS(I,J,II) = VFS(I,J,II)
         VSS(I,J,II) = VFS(I,J,II)
  110 CONTINUE
  111 CONTINUE
C
      ELSE
cxx      DO 410  I=1,MM
      DO 411  I=1,MM
      DO 410  J=1,MM
cxx  410 VP(I,J) = VPS(I,J,II+1)
      VP(I,J) = VPS(I,J,II+1)
  410 CONTINUE
  411 CONTINUE
C
cc      CALL  GINVRS( VP,VDET,MM,40 )
      CALL  GINVRS( VP,VDET,MM )
C
cxx      DO 420  I=1,MM
cxx      DO 420  L=1,NC
      DO 422  I=1,MM
      DO 421  L=1,NC
      WRK(I,I0(L)+M(L)) = VFS(I,I0(L)+1,II)*A(M(L),L)
      DO 420  J=1,M(L)-1
cxx  420 WRK(I,I0(L)+J) = VFS(I,I0(L)+1,II)*A(J,L) + VFS(I,I0(L)+J+1,II)
      WRK(I,I0(L)+J) = VFS(I,I0(L)+1,II)*A(J,L) + VFS(I,I0(L)+J+1,II)
  420 CONTINUE
  421 CONTINUE
  422 CONTINUE
C
cxx      DO 440  I=1,MM
      DO 441  I=1,MM
      DO 440  J=1,MM
      SUM = 0.0D0
      DO 430 IJ=1,MM
cxx  430 SUM = SUM + WRK(I,IJ)*VP(IJ,J)
      SUM = SUM + WRK(I,IJ)*VP(IJ,J)
  430 CONTINUE
cxx  440 SGAIN(I,J) = SUM
      SGAIN(I,J) = SUM
  440 CONTINUE
  441 CONTINUE
C
cxx      DO 450  I=1,MM
      DO 451  I=1,MM
      XS(I) = XFS(I,II)
      DO 450  J=1,MM
      WRK(I,J) = 0.0D0
cxx  450 VS (I,J) = VFS(I,J,II)
      VS (I,J) = VFS(I,J,II)
  450 CONTINUE
  451 CONTINUE
C
cxx      DO 460  J=1,MM
      DO 461  J=1,MM
      DO 460  I=1,MM
cxx  460 XS(I) = XS(I) + SGAIN(I,J)*(XSS(J,II+1) - XPS(J,II+1))
      XS(I) = XS(I) + SGAIN(I,J)*(XSS(J,II+1) - XPS(J,II+1))
  460 CONTINUE
  461 CONTINUE
C
cxx      DO 470  J=1,MM
cxx      DO 470 IJ=1,MM
      DO 472  J=1,MM
      DO 471 IJ=1,MM
      DO 470  I=1,MM
cxx  470 WRK(I,J) = WRK(I,J) + SGAIN(I,IJ)*(VSS(IJ,J,II+1)-VPS(IJ,J,II+1))
      WRK(I,J) = WRK(I,J) + SGAIN(I,IJ)*(VSS(IJ,J,II+1)-VPS(IJ,J,II+1))
  470 CONTINUE
  471 CONTINUE
  472 CONTINUE
C
cxx      DO 480  J=1,MM
cxx      DO 480 IJ=1,MM
      DO 482  J=1,MM
      DO 481 IJ=1,MM
      DO 480  I=1,MM
cxx  480 VS(I,J) = VS(I,J) + WRK(I,IJ)*SGAIN(J,IJ)
      VS(I,J) = VS(I,J) + WRK(I,IJ)*SGAIN(J,IJ)
  480 CONTINUE
  481 CONTINUE
  482 CONTINUE
      DO 485 I=1,MM
cxx  485 IF( VS(I,I).LT.0.0D0 )  VS(I,I) = 0.0D0
      IF( VS(I,I).LT.0.0D0 )  VS(I,I) = 0.0D0
  485 CONTINUE
C
cxx      DO 490  I=1,MM
      DO 491  I=1,MM
      XSS(I,II) = XS(I)
      DO 490  J=1,MM
cxx  490 VSS(I,J,II) = VS(I,J)
      VSS(I,J,II) = VS(I,J)
  490 CONTINUE
  491 CONTINUE
      END IF
C
  500 CONTINUE
C
      RETURN
      E N D
C------------------------------------------------------------------
C
      FUNCTION  ID( K )                                                 
C                                                                       
C  ...  ID = 1:    IF K > 0                                             
C       ID = 0:    OTHERWISE                                            
C                                                                       
      ID = 0                                                            
      IF( K .GT. 0 )  ID = 1                                            
      RETURN                                                            
      E N D                                                             

      SUBROUTINE INIT(IX)
      INTEGER :: IX, v(8)
C
      if ( IX .ge. 0 ) then
          call init_genrand64(IX)
      else
         call date_and_time( values=v )
         call init_genrand64( sum(v) )
      endif
      RETURN
      END

      DOUBLE PRECISION FUNCTION  GAUSS( X,PARAM )
C
C  ...  Gaussian (normal) distribution  ...
C
C     Inputs:
C        X:
C        PARAM(1):  mean
C        PARAM(2):  variance
C     Output:
C        GAUSS:     density at X
C
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cxx      DIMENSION  PARAM(2)
      REAL(8) :: X, PARAM(2), C1
      DATA  C1  /2.506628275D0/
C
      GAUSS = DEXP( -(X-PARAM(1))**2/(2*PARAM(2)) )/(C1*DSQRT(PARAM(2)))
      RETURN
      E N D
      DOUBLE PRECISION FUNCTION  PEARSN( X,PARAM )
C
C  ...  Pearson family of  distributions  ...
C
C     Inputs:
C        X:
C        PARAM(1):  location parameter, mu
C        PARAM(2):  dispersion parameter, tau square
C        PARAM(3):  shape parameter
C     Output:
C        PEARSN:    density at X
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  PARAM(3)
      REAL(8) :: X, PARAM(3), PI, dgammafn
      DATA  PI/3.1415926535D0/
C
CXX      PEARSN = DGAMMA(PARAM(3))/DGAMMA(PARAM(3)-0.5D0)
      PEARSN = dgammafn(PARAM(3))/dgammafn(PARAM(3)-0.5D0)
     *                  /DSQRT(PI)*PARAM(2)**(PARAM(3)-0.5D0)
     *                  /((X-PARAM(1))**2 + PARAM(2))**PARAM(3)
      RETURN
C
      END

      DOUBLE PRECISION FUNCTION  DBLEXP( X,PARAM )
C
C  ...  double exponential distribution  f(x) = exp(x - exp(x))  ...
C
C     Inputs:
C        X:
C     Output:
C        DBLEXP:    density at X
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      dimension PARAM(3)
      REAL(8) :: X, PARAM(3)
C
cxx 2018/07/02      DBLEXP = DEXP( X-DEXP(X) )
      DBLEXP = DEXP( (X-PARAM(1))-DEXP(X-PARAM(1)) )
      RETURN
C
      E N D

      DOUBLE PRECISION FUNCTION  CAUCHY( X,PARAM )
C
C  ...  Cauchy distribution  ...
C
C     Inputs:
C        X:
C        PARAM(1):  location parameter, mu
C        PARAM(2):  dispersion parameter, tau square
C     Output:
C        CAUCHY:    density at X
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  PARAM(2)
      REAL(8) :: X, PARAM(2), PI
      DATA  PI /3.1415926535D0/
C
      CAUCHY = DSQRT( PARAM(2) )/(PARAM(2) + (X-PARAM(1))**2)/PI
      RETURN
C
      E N D
