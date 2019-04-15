C     PROGRAM 7.3  MARLSQ
      SUBROUTINE MARLSQF( Y,N,L,LAG,A,V,LMAX,AIC )
C
      INCLUDE 'TSSS_f.h'
C
C  ...  MAR model fitting (least squares method)  ...
C
C     Inputs:
C        MT:    Input device
C        LAG:   Heighest AR order
C     Inputs for READMD
C        TITLE: Tile of the data set
C        N:     Data length
C        L:     Dimension
C        IFM:   Control variable for reading order
C        FORM:  Reading format
C        Y(N,J):   ID-dimensional time series
C
cc      PARAMETER( MJ=1000,MJ1=200,MJ2=5,MJ3=20,MJ4=MJ2*(MJ3+1),IDEV=1 )
ccxx      PARAMETER( MJ1=200 )
cxx      IMPLICIT  REAL*8 (A-H,O-Z)
cc      DIMENSION  Y(MJ,MJ2), X(MJ1,MJ4), D(MJ1) 
cc      DIMENSION  A(MJ2,MJ2,MJ3) , V(MJ2,MJ2)
cc      DIMENSION  Y(N,L), X(MJ1,L*(LAG+1))
cxx      DIMENSION  Y(N,L), X((L+1)*(LAG+1),L*(LAG+1))
cxx      DIMENSION  A(L,L,LAG) , V(L,L)
C
      INTEGER :: N, L, LAG, LMAX
      REAL(8) :: Y(N,L), A(L,L,LAG), V(L,L), AIC
      REAL(8) :: X((L+1)*(LAG+1),L*(LAG+1))
C
      EXTERNAL   SETMAR
C
      MJ1 = (L+1)*(LAG+1)
C
      MJ = N
      MJ2 = L
      MJ3 = LAG
      MJ4 = L*(LAG+1)
C
cc      IPR = 3
cc      READ( 5,* )  LAG
C
C  ...  Read multi-variate time series  ...
C
cc      CALL  READMD( IDEV,MJ,Y,N,L )
      N0 = 0
      NMK = N - LAG
C
C  ...  Housholder reduction  ...
C
cc      CALL  MREDCT( SETMAR,Y,D,NMK,N0,LAG,L,MJ,MJ1,X )
      CALL  MREDCT( SETMAR,Y,NMK,N0,LAG,L,MJ,MJ1,X )
C
C  ...  Least squares method for MAR model  ...
C
cc      CALL  MARFIT( X,D,NMK,L,LAG,MJ1,MJ2,MJ3,IPR,A,V,LMAX,AIC )
      CALL  MARFIT2( X,NMK,L,LAG,MJ1,IPR,A,V,LMAX,AIC )
C
C  ...  Print out fitted model  ...
C
cc      CALL  PRMAR2( A,V,L,LMAX,LAG,N,MJ2,MJ3 )
C     WRITE(2,*)  LMAX
C     DO 97 I=1,L
C  97 WRITE(2,600) (V(I,J),J=1,L)
C     DO 98 II=1,LMAX
C     DO 98 I=1,L
C  98 WRITE(2,600) (A(I,J,II),J=1,L)
C 600 FORMAT( 3D15.7 )
C
cc      STOP
      RETURN
      E N D
cc      SUBROUTINE  MARFIT( X,D,N,ID,M,MJ1,MJ2,MJ3,IPR,B,E,LMAX,AICS )
      SUBROUTINE  MARFIT2( X,N,ID,M,MJ1,IPR,B,E,LMAX,AICS )
C
C  ...  Multivariate AR model fitting (Householder method)  ...
C
C       Inputs:
C          X:     Householder reduced form
C          D:     Working area
C          N:     Data length
C          ID:    Dimension of the series
C          M:     Highest AR order
C          IPR:   Print out control
C          MJ1, MJ2, MJ3:  Adjustable dimensions
C       Outputs:
C          B:     AR-coefficient matrices
C          E:     innovation covariance matrix
C          LMAX:  Order of the MAICE model
C          AICS:  Total AIC of the model
C
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cc      DIMENSION  X(MJ1,1), D(1), E(MJ2,1), B(MJ2,MJ2,MJ3)
cc      DIMENSION  A(100), AIC(51), SD(51), EX(10)
cc      DIMENSION  IND(100), JND(100)
cxx      DIMENSION  X(MJ1,(M+1)*ID), E(ID,ID), B(ID,ID,M)
cxx      DIMENSION  A((M+1)*ID), AIC(M+1), SD(M+1), EX(ID)
cxx      DIMENSION  IND((M+1)*ID), JND((M+1)*ID)
C
      INTEGER :: N, ID, M, MJ1, IPR, LMAX
      REAL(8) :: X(MJ1,(M+1)*ID), B(ID,ID,M), E(ID,ID), AICS
      INTEGER :: IND((M+1)*ID), JND((M+1)*ID)
      REAL(8) :: A((M+1)*ID), AIC(M+1), SD(M+1), EX(ID), PI2, SUM,
     1           AICSUM, AICMIN, SDMIN, AICX
C
      DATA  PI2/6.28318531D0/
C
      MJ2 = ID
      MJ3 = M
C
      MD = M*ID
      MD2 = (M+1)*ID
      AICSUM = 0.D0
      LMAX = 0
cxx      DO 20  I=1,MJ2
cxx      DO 20  J=1,MJ2
cxx      DO 10  K=1,MJ3
cxx   10 B(I,J,K) = 0.D0
cxx   20 E(I,J) = 0.D0
      B(1:MJ2,1:MJ2,1:MJ3) = 0.D0
      E(1:MJ2,1:MJ2) = 0.D0
      DO 30  I=1,MD2
cxx   30 JND(I) = I
      JND(I) = I
   30 CONTINUE
C
      DO 200     II=1,ID
C
      JJ = II - 1
      KK = MD + JJ
      IF( II .GT. 1 )  THEN
C
C  ...  Addition of instantaneous variable  ...
C
      DO 110  I=1,MD2
      J = JND(I)
cxx  110 IND(J) = I
      IND(J) = I
  110 CONTINUE
      JND(1) = 1
      DO 120  I=1,JJ
cxx  120 JND(I) = MD+I
      JND(I) = MD+I
  120 CONTINUE
      DO 130  I=1,MD
cxx  130 JND(I+JJ) = I
      JND(I+JJ) = I
  130 CONTINUE
C     I1 = JJ + 1
      DO 140  I=JJ+1,ID
      J = MD + I
cxx  140 JND(J) = J
      JND(J) = J
  140 CONTINUE
C
cc      CALL  HUSHL1( X,D,MJ1,MD2,MD2,1,IND,JND )
      CALL  HUSHL1( X,MJ1,MD2,MD2,1,IND,JND )
      END IF
C
C  ...   AIC's of the models for II-th variable  ...
C
      DO 160  I=1,M+1
      K = (I-1)*ID + JJ
      SUM = 0.0D0
      DO 150 IJ=K+1,KK+1
cxx  150 SUM = SUM + X(IJ,KK+1)**2
      SUM = SUM + X(IJ,KK+1)**2
  150 CONTINUE
      SD(I) = SUM/N
cxx  160 AIC(I) = N*DLOG( PI2*SD(I) ) + N + 2*(K+1)
      AIC(I) = N*DLOG( PI2*SD(I) ) + N + 2*(K+1)
  160 CONTINUE
C
C  ...  Order determination by AIC  ...
C
      CALL  MAICE( AIC,SD,M,IPR-2,AICMIN,SDMIN,IMIN )
C
      K0 = IMIN*ID + JJ
      CALL  SRCOEF( X,K0,KK,N,MJ1,JND,A,EX(II),AICX )
      DO 170 I=1,II-1
cxx  170 E(II,I)   = -A(I)
      E(II,I)   = -A(I)
  170 CONTINUE
      E(II,II)  = 1.0D0
cxx      DO 180 IJ=1,IMIN
      DO 181 IJ=1,IMIN
      DO 180 J=1,ID
cxx  180 B(II,J,IJ) = A((IJ-1)*ID+J+II-1)
      B(II,J,IJ) = A((IJ-1)*ID+J+II-1)
  180 CONTINUE
  181 CONTINUE
      AICSUM = AICSUM + AICMIN
      IF( IMIN.GT.LMAX )  LMAX = IMIN
C
  200 CONTINUE
C
cc      IF(IPR.GE.1)  WRITE( 6,3 )
cxx      CALL  MCOEF( B,E,EX,ID,LMAX,MJ2,MJ3 )
      CALL  MCOEF( B,E,EX,ID,LMAX,MJ3 )
cc      IF(IPR.GE.1)  WRITE( 6,611 )  AICSUM
      AICS = AICSUM
C
      RETURN
cxx    3 FORMAT\( 1H0,132(1H-) )
cxx   15 FORMAT( 1H0,'*****  FINAL ESTIMATE  *****' )
cxx  611 FORMAT( 1H0,'AIC =',F15.3 )
      END
cxx      SUBROUTINE  MCOEF( B,E,EX,ID,LMAX,MJ2,MJ3 )
      SUBROUTINE  MCOEF( B,E,EX,ID,LMAX,MJ3 )
C
C  ...  AR coefficient matrices and innovation covariance matrix  ...
C
C       Inputs:
C          B:     Regression coefficient matrix
C          E:     B(0)
C          EX:    Residual variances of orthogonarized model
C          ID:    Dimension
C          LMAX:  AR-order
C          MJ2:   Adjustable dimension of B and E
C          MJ3:   djustable dimension of B
C       Outputs:
C          B:     AR-coefficient matrix
C          E:     Innovation covaraince matrix
C
cxx      IMPLICIT  REAL*8 ( A-H,O-Z )
cc      DIMENSION  B(MJ2,MJ2,MJ3), E(MJ2,1), EX(1), EE(24,24)
cc      MJ5 = 24
cxx      DIMENSION  B(ID,ID,MJ3), E(ID,ID), EX(ID), EE(ID,ID)
C
      INTEGER :: ID, LMAX, MJ3 
      REAL(8) :: B(ID,ID,MJ3), E(ID,ID), EX(ID)
      REAL(8) :: EE(ID,ID), SUM
C
cc      CALL  TRIINV( E,ID,MJ2,MJ5,EE )
      CALL  TRIINV( E,ID,EE )
C
cxx      DO 30  II=1,LMAX
cxx      DO 20  I=1,ID
      DO 32  II=1,LMAX
      DO 21  I=1,ID
      DO 20  J=1,ID
      SUM = 0.D0
      DO 10  JJ=1,I
cxx   10 SUM = SUM + EE(I,JJ)*B(JJ,J,II)
      SUM = SUM + EE(I,JJ)*B(JJ,J,II)
   10 CONTINUE
cxx   20 E(I,J) = SUM
      E(I,J) = SUM
   20 CONTINUE
   21 CONTINUE
cxx      DO 30  I=1,ID
      DO 31  I=1,ID
      DO 30  J=1,ID
cxx   30 B(I,J,II) = E(I,J)
      B(I,J,II) = E(I,J)
   30 CONTINUE
   31 CONTINUE
   32 CONTINUE
cxx      DO 50  I=1,ID
      DO 51  I=1,ID
      DO 50  J=1,I
      SUM = 0.D0
      DO 40  II=1,J
cxx   40 SUM = SUM + EE(I,II)*EE(J,II)*EX(II)
      SUM = SUM + EE(I,II)*EE(J,II)*EX(II)
   40 CONTINUE
      E(I,J) = SUM
cxx   50 E(J,I) = SUM
      E(J,I) = SUM
   50 CONTINUE
   51 CONTINUE
cc      DO 99 I=1,ID
cc   99 WRITE(6,*) (E(I,J),J=1,ID)
C
      RETURN
      E N D
cc      SUBROUTINE  MREDCT( MSETX,Z,D,NMK,N0,LAG,ID,MJ,MJ1,X )
      SUBROUTINE  MREDCT( MSETX,Z,NMK,N0,LAG,ID,MJ,MJ1,X )
C
C  ...  This subroutine computes Householder reduced form of the
C       data matrix for fitting MAR model  ...
C
C       Inputs:
C          MSETX:  Function name
C          Z:      Multivariate time series
C          D:      Working area
C          NMK:    Actual data length
C          N0:     Origin of new data
C          LAG:    Highest order of the AR model
C          ID:     Dimension of observations
C          MJ,MJ1: Adjustable dimension
C         OUTPUT:
C          X:      Householder reduced form (upper triangular form)
C
cxx      IMPLICIT  REAL*8 (A-H,O-Z)
cc      DIMENSION  X(MJ1,1) , D(1), Z(MJ,1)
cx      DIMENSION  X(MJ1,1) , Z(MJ,1)
cxx      DIMENSION  X(MJ1,(LAG+1)*ID) , Z(MJ,ID)
C
      INTEGER :: NMK, N0, LAG, ID, MJ, MJ1
      REAL(8) :: Z(MJ,ID), X(MJ1,(LAG+1)*ID)
C
      L = MIN0( NMK,MJ1 )
      KD1 = (LAG+1)*ID
      N1 = L
C
      CALL  MSETX( Z,N0,L,LAG,ID,MJ,MJ1,0,X )
cc      CALL  HUSHLD( X,D,MJ1,L,KD1 )
      CALL  HUSHLD( X,MJ1,L,KD1 )
C
      IF( N1 .GE. NMK )     RETURN
  100 L = MIN0( NMK-N1,MJ1-KD1 )
      LK = L + KD1
      N2 = N0 + N1
C
      CALL  MSETX( Z,N2,L,LAG,ID,MJ,MJ1,1,X )
cc      CALL  HUSHLD( X,D,MJ1,LK,KD1 )
      CALL  HUSHLD( X,MJ1,LK,KD1 )
      N1 = N1 + L
      IF( N1 .LT. NMK )     GO TO 100
C
      RETURN
      E N D
      SUBROUTINE  SETMAR( Z,N0,L,LAG,ID,MJ,MJ1,JSW,X )
C
C  ...  Make matrix X for multi-variate AR model  ...
C
C       Inputs:
C          Z:      Multivariate time series
C          N0:     Origin of new data
C          L:      Number of new observations
C          LAG:    Highest order of AR model
C          ID:     Dimension of time series
C          MJ,MJ1: Adjustable dimension
C          JSW:    =0   To make initail L*(LAG+1) data matrix
C                  =1   To augment original (LAG+1)*(LAG+1) matrix X
C                       by L*(LAG+1) data matrix of additional obs.
C       Outputs
C          X:      L*(LAG+1) matrix            if  JSW = 0
C                  (LAG+1+L)*(LAG+1) matrix    if  JSW = 1
C
cc      REAL*8  X(MJ1,1), Z(MJ,1)
cxx      REAL*8  X(MJ1,(LAG+1)*ID), Z(MJ,ID)
C
      INTEGER :: N0, L, LAG, ID, MJ, MJ1, JSW
      REAL(8) :: Z(MJ,ID), X(MJ1,(LAG+1)*ID)
C
      KD = LAG*ID
      KD1 = (LAG+1)*ID
      I0 = 0
      IF( JSW.EQ.1 )  I0 = KD1
C
cxx      DO 30  II=1,L
      DO 31  II=1,L
      I1 = N0 + LAG + II
      I2 = I0 + II
      DO 10  J=1,ID
      J2 = KD + J
cxx   10 X(I2,J2) = Z(I1,J)
      X(I2,J2) = Z(I1,J)
   10 CONTINUE
      DO 30  JJ=1,LAG
      I1 = I1 - 1
      J1 = (JJ-1)*ID
      DO 20 J=1,ID
      J2 = J1 + J
cxx   20 X(I2,J2) = Z(I1,J)
      X(I2,J2) = Z(I1,J)
   20 CONTINUE
   30 CONTINUE
   31 CONTINUE
C
      RETURN
      E N D
cc      SUBROUTINE  TRIINV( X,M,MJ,MJ1,Y )
      SUBROUTINE  TRIINV( X,M,Y )
C
C  ...  Lower triangular matrix inversion  ...
C
C       Inputs:
C          X:    Triangular matrix with unit diagonal elements
C          M:    Dimension of matrix X
C          MJ:   Adjustable dimension of X
C          MJ1:  Adjustable dimension of Y
C       Output:
C          Y:    Inverse of X
C
cxx      IMPLICIT  REAL * 8  ( A-H , O-Z )
cc      DIMENSION  X(MJ,1) , Y(MJ1,1)
cxx      DIMENSION  X(M,M) , Y(M,M)
C
      INTEGER :: M
      REAL(8) :: X(M,M) , Y(M,M)
      REAL(8) :: SUM
C
cxx      DO 10  I=1,M-1
      DO 11  I=1,M-1
      DO 10  J=I,M
cxx   10 Y(I,J) = 0.0D0
      Y(I,J) = 0.0D0
   10 CONTINUE
   11 CONTINUE
      DO 20  I=1,M
cxx   20 Y(I,I) = 1.0D0
      Y(I,I) = 1.0D0
   20 CONTINUE
cxx      DO 40  J=1,M-1
      DO 41  J=1,M-1
      DO 40  I=J+1,M
      SUM = 0.0D0
      DO 30  II=1,I-J
      JJ = II + J - 1
cxx   30 SUM = SUM + X(I,JJ) * Y(JJ,J)
      SUM = SUM + X(I,JJ) * Y(JJ,J)
   30 CONTINUE
cxx   40 Y(I,J) = -SUM
      Y(I,J) = -SUM
   40 CONTINUE
   41 CONTINUE
      RETURN
      E N D
      SUBROUTINE  MAICE( AIC,SD,K,ISW,AICM,SDM,IMIN )
C
C       INPUTS:
C          AIC:   VECTOR OF AIC'S
C          SD:    VECTOR OF INNOVATION VARIANCES
C          K:     UPPER LIMIT OF THE ORDER
C          ISW:   =0   OUTPUTS ARE SUPPRESSED
C                 >0   AIC'S ARE DIPLAIED
C       OUTPUTS:
C          AICM:  MINIMUM AIC
C          SDM:   MAICE INNOVATION VARIANCE
C          IMIN:  MAICE ORDER
C
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cc      DIMENSION  AIC(1), SD(1)
cxx      DIMENSION  AIC(K+1), SD(K+1)
C
      INTEGER :: K, ISW, IMIN 
      REAL(8) :: AIC(K+1), SD(K+1), AICM, SDM
C
C       SEARCH FOR THE MINIMUM OF AIC(I)
C
      IMIN = 0
      SDM  = SD(1)
      AICM = AIC(1)
      DO 20 I=1,K
      IF( AIC(I+1).LT.AICM )  THEN
         IMIN = I
         SDM  = SD(I+1)
         AICM = AIC(I+1)
      END IF
   20 CONTINUE
      IF( ISW .LE. 0 )     RETURN
C
C       DISPLAY OF AIC'S
C
cc      WRITE( 6,5 )
cc      DO 30 I=1,K+1
cc      II = I - 1
cc      DIC = AIC(I) - AICM
cc   30 WRITE( 6,6 )  II, SD(I), AIC(I), DIC
cc      WRITE( 6,7 )  IMIN, AICM, SDM
C
      RETURN
cxx    5 FORMAT( 1H ,4X,'M',9X,'SIG2(M)',13X,'AIC(M)',7X,'AIC(M)-AICMIN' )
cxx    6 FORMAT( 1H ,I5,D20.10,2F16.3 )
cxx    7 FORMAT( 1H0,'MAICE ORDER =',I3,5X,'AIC =',F15.3,5X,'SIG2 =',
cxx     *             D17.10 )
      E N D
cc      SUBROUTINE  HUSHL1( X,D,MJ,K,L,M,IND,JND )
      SUBROUTINE  HUSHL1( X,MJ,K,L,M,IND,JND )
C
C  ...  Householder trasnformation of the matrix X  ...
C
C       Inputs:
C         X(I,J): (K+1)*(K+1) matrix with the form specified by IND(I)
C         D(I):   Working area
C         MJ:     Adjustable dimension of X
C         K:      Number of columns of X
C         M:      Starting position of the Householder transformation
C         L:      End position of the Householder transformation
C         IND:    Specification of the present form of X
C                 << IND(J) = I means J-th regressor is variable I>>
C         JND:    Specification of the required form of X
C                 << JND(I) = J means I-th regressor is variable J >>
C       Output:
C         X(I,J)  Transformed matrix with the form specified by JND(I)
C
cxx      IMPLICIT  REAL*8 (A-H,O-Z)
cc      DIMENSION  X(MJ,1), D(40), IND(40), JND(40)
cxx      DIMENSION  X(MJ,K), D(MJ), IND(K), JND(K)
C
      INTEGER :: MJ, K, L, M, IND(K), JND(K)
      REAL(8) :: X(MJ,K)
      REAL(8) :: D(MJ), TOL, F, G, H, S
C
      TOL = 1.0D-60
C
      NN = 0
      DO 100  II=M,L
      JJ = JND(II)
      NN = MAX0( NN,IND(JJ) )
      H = 0.0D0
      DO 10  I=II,NN
      D(I) = X(I,JJ)
cxx   10 H = H + D(I)**2
      H = H + D(I)**2
   10 CONTINUE
      IF( H.GT.TOL )  GO TO 20
      G = 0.0D0
      GO TO 100
   20 G = DSQRT( H )
      F = X(II,JJ)
      IF( F.GE.0.0D0 )  G = -G
      D(II) = F - G
      H = H - F * G
C
      DO 30  I=II+1,NN
cxx   30 X(I,JJ) = 0.D0
      X(I,JJ) = 0.D0
   30 CONTINUE
      IF( II.EQ.K )  GO TO 100
      DO 60  J1=II+1,K
      J = JND(J1)
      S = 0.0D0
      DO 40  I=II,NN
cxx   40 S = S + D(I)*X(I,J)
      S = S + D(I)*X(I,J)
   40 CONTINUE
      S = S / H
      DO 50  I=II,NN
cxx   50 X(I,J) = X(I,J) - D(I)*S
      X(I,J) = X(I,J) - D(I)*S
   50 CONTINUE
   60 CONTINUE
cxx  100 X(II,JJ) = G
      X(II,JJ) = G
  100 CONTINUE
      RETURN
      E N D
      SUBROUTINE  SRCOEF( X,M,K,N,MJ,JND,A,SIG2,AIC )
C
C  ...  Subset regression coefficients, residual variance and AIC  ...
C
C       Iuputs:
C         X:     Householder reduced form with index vector JND(I)
C         M:     Number of regressors
C         K:     Heighest order of the models
C         N:     Data length
C         MJ:    Adjustable dimension of X
C         JND(I): Specification of the I-th regressor
C       Outputs:
C         A(I):  Regression coefficients
C         SIG2:  Innovation variance
C         AIC:   AIC
C
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cc      DIMENSION  X(MJ,1), A(1), JND(1)
cxx      DIMENSION  X(MJ,K+1), A(M), JND(M)
C
      INTEGER :: M, K, N, MJ, JND(M)
      REAL(8) :: X(MJ,K+1), A(M), SIG2, AIC
      REAL(8) :: PI2, SUM, SD
C
      DATA  PI2/6.28318531D0/
C
C  ...  REGRESSION COEFFICIENTS...
C
      L = JND(M)
      A(M) = X(M,K+1)/X(M,L)
      DO 10  II=1,M-1
      I = M - II
      SUM = X(I,K+1)
      DO 20  J=I+1,M
      L = JND(J)
cxx   20 SUM = SUM - A(J)*X(I,L)
      SUM = SUM - A(J)*X(I,L)
   20 CONTINUE
      L = JND(I)
cxx   10 A(I) = SUM/X(I,L)
      A(I) = SUM/X(I,L)
   10 CONTINUE
C
C  ...  RESIDUAL VARIANCE AND AIC  ...
C
      SD = 0.0D0
      DO 30  I=M+1,K+1
cxx   30 SD = SD + X(I,K+1)**2
      SD = SD + X(I,K+1)**2
   30 CONTINUE
      SIG2 = SD/N
      AIC  = N*DLOG( PI2*SIG2 ) + N + 2.0D0*(M+1)
C
      RETURN
      E N D

