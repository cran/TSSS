      SUBROUTINE ARMA(M,L,A,B,SIG2,K,KMAX,NF,G,COV,PAR,SP,
     * ROOTA,ROOTB,IER,JER)
C
      INCLUDE 'TSSS.h'
C
C     PROGRAM 6.1  ARMA
C
C  ...  This program analyses time series via ARMA modeling  ...
C
C     Inputs:
C        IDEV:  Input devie (=5: CONSOLE)
C        M:     AR order
C        L:     MA order
C        SIG2:  Innovation variance
C        A(I):  AR coefficients
C        B(I):  MA coefficients
C        K:     Maximum lag of autocovariance function (K<MJ+1)
C     Parameter:
C        MJ:    Adjustable dimension of A, B, etc..
C        NF:    Number of frequencies in evaluating spectrum
C     Programmed by Y.I and G.K.
C
C     Outputs:
C         IER:  =1 : MATRIX WITH ZERO ROW IN DECOMPOSE
C               =2 : SINGULAR MATRIX IN DECOMPOSE.ZERO DIDIVIDE IN SOLVE
C               =3 : CONVERGENCE IN IMPRUV.MATRIX IS NEARLY SINGULAR
C         JER:  =1 : NON-CONVERGENCE AT POLYRT
C
cc      PARAMETER( NF=200,MJ=100,IDEV=1 )
cxx      IMPLICIT REAL*8( A-H,O-Z )
cc      DIMENSION  A(MJ), B(MJ), PAR(MJ), ROOTA(MJ,2), ROOTB(MJ,2)
cc      DIMENSION  G(0:MJ), COV(0:MJ), SP(0:NF)
cc      DIMENSION  WRK1(0:MJ), WRK2(0:MJ), WRK3(MJ)
cc      DATA  PAR/MJ*0.0D0/
cxx      DIMENSION  A(M), B(L), PAR(K), ROOTA(M,2), ROOTB(L,2)
cxx      DIMENSION  G(0:KMAX), COV(0:K), SP(0:NF)
cxx      DIMENSION  WRK1(0:K), WRK2(0:K), WRK3(K,K)
C
      INTEGER M, L, K, KMAX, NF, IER, JER
      DOUBLE PRECISION A(M), B(L), SIG2, G(0:KMAX), COV(0:K), PAR(K),
     1                 SP(0:NF), ROOTA(M,2), ROOTB(L,2)
c local
      INTEGER MAR, JER1, JER2
      DOUBLE PRECISION WRK1(0:K), WRK2(0:K), WRK3(K,K)
C
cc      WRITE( 6,* )  'K = ?'
cc      READ( 5,* )   K
cc      OPEN( IDEV,FILE='arma.dat' )
cc      READ(IDEV,*)  M, L
cc      READ(IDEV,*)  SIG2
cc      READ(IDEV,*)  (A(I),I=1,M)
cc      READ(IDEV,*)  (B(I),I=1,L)
cc      CLOSE( IDEV )
C
cc      KMAX = MAX(M,L,K)
      CALL  IMPULS( M,L,A,B,K,G )
cc      CALL  ARMCOV( M,L,A,B,SIG2,K,COV )
      CALL  ARMCOV( M,L,A,B,SIG2,K,COV,KMAX,IER )
      if( ier.ne.0 ) return
c------------
      PAR(1:K) = 0.0D0
c------------
      CALL  PARCOR( A,M,PAR )
      CALL  ARCOEF( PAR,M,A )
cc      IF( L.GT.0 )  CALL  ARYULE( COV,1000,K,WRK1,WRK2,PAR,WRK3,MAR )
      IF( L.GT.0 )  CALL  ARYULE( COV,0,K,WRK1,WRK2,PAR,WRK3,MAR )
      CALL  ARMASP( A,M,B,L,SIG2,NF,SP )
cc      CALL  CHROOT( A,M,ROOTA,MJ )
cc      CALL  CHROOT( B,L,ROOTB,MJ )
      CALL  CHROOT( A,M,ROOTA,M,JER1 )
      CALL  CHROOT( B,L,ROOTB,L,JER2 )
      JER = JER1
      IF( JER2 .NE. 0 ) JER = JER + JER2 + 1
cc      CALL  PRARMA( M,L,A,B,G,K,COV,K,PAR,SP,NF,ROOTA,ROOTB,MJ )
C      CALL  PTARMA( G,K,COV,K,PAR,SP,NF,ROOTA,M,ROOTB,L,MJ )
      RETURN
      E N D
cc      SUBROUTINE  CHROOT( A,M,ROOT,MJ )
      SUBROUTINE  CHROOT( A,M,ROOT,MJ,IER )
C
C  ...  Characteristic roots of the AR or MA operator  ...
C
C     Inputs:
C        A:     Vector of AR or MA coefficients
C        M:     Order of the model
C        MJ:    Adjustable dimension of ROOT
C     Output:
C        ROOT:  Characteristic roots (real part,imaginary part)
C
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cxx      DIMENSION  A(M), ROOT(MJ,2)
cc      DIMENSION  C(50), CW(50)
cxx      DIMENSION  C(M+1)
C
      INTEGER M, MJ, IER
      DOUBLE PRECISION A(M), ROOT(MJ,2)
c local
      INTEGER I, MMAX
      DOUBLE PRECISION C(M+1)
C
      IER = 0
C
      IF( M.EQ.0 )  RETURN
      DO 10  I=1,M
cxx   10 C(I) = -A(M-I+1)
      C(I) = -A(M-I+1)
   10 CONTINUE
      C(M+1) = 1.0D0
      MMAX = M
C
C  ... Characteristic roots of operator 1-A(1)*S- ... -A(M)*S**M=0  ...
C
cc      CALL  POLYRT( C,CW,MMAX,ROOT(1,1),ROOT(1,2),IER )
      CALL  POLYRT( C,MMAX,ROOT(1,1),ROOT(1,2),IER )
cc      IF( IER.NE.0 )   WRITE(6,600)
C
      RETURN
cxx  600 FORMAT( 1H0,'*****  NON-CONVERGENCE AT POLYRT  *****' )
      E N D
cc      SUBROUTINE  POLYRT( A,B,M,ROOTR,ROOTI,IER )
      SUBROUTINE  POLYRT( A,M,ROOTR,ROOTI,IER )
C
C  ...  This subroutine finds the roots of the equation
C            A(1) + A(2)*S + ... + A(M)*(S**M) = 0
C       by Newton-Raphson method
C
C     Inputs:
C        A:     Coefficients of the equation
C        B:     Working area
C        M:     Degree of the polynomial
C     Outputs:
C        ROOTR:   Real parts of the roots
C        ROOTI:   Imaginary parts of the roots
C        IER:     Error code to indicate non-convergence
C
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cc      DIMENSION A(1), B(1), ROOTR(1), ROOTI(1)
cxx      DIMENSION A(M+1), B(M+3), ROOTR(M), ROOTI(M)
C
      INTEGER M, IER 
      DOUBLE PRECISION A(M+1), ROOTR(M), ROOTI(M)
c local
      INTEGER I, IFIT, ISW, ICT, ITEM, IN, J, JSW, K, KX, KXX, KJI,
     1        K2, L
      DOUBLE PRECISION B(M+3), XPR, YPR, X, XO, Y, YO, UX, UY, U, V,
     1                 XT, YT, XT2, YT2, FI, SUM, DX, DY, TEM, ALPH
C
      IFIT = 0
      ISW = 0
      JSW = 0
      K = M
      IER = 0
      KX = K
      KXX = K+1
      KJI = KXX
      K2 = 1
c------------
      XPR = 0.0D0
      YPR = 0.0D0
c------------
      DO 10  I=1,KXX
      A(I) = A(I)/A(K+1)
      J = KXX - I+1
cxx   10 B(J) = A(I)
      B(J) = A(I)
   10 CONTINUE
   20 XO = 0.5D-02
      YO = 0.1D-01
      IN = 0
   30 X = XO
C
      XO = -10.D0*YO
      YO = -10.D0*X
      X = XO
      Y = YO
      ICT = 0
      IN = IN+1
      GO TO 50
C
   40 IFIT = 1
      XPR = X
      YPR = Y
C
   50 UX = 0.D0
      UY = 0.D0
      V  = 0.D0
      YT = 0.D0
      XT = 1.D0
      U = B(K+1)
      IF( U .EQ. 0.D0 )  GO TO 140
      DO 60  I=1,K
      L = K-I+1
      XT2 = X*XT - Y*YT
      YT2 = X*YT+Y*XT
      U = U + B(L)*XT2
      V = V + B(L)*YT2
      FI = I
      UX = UX + FI*XT*B(L)
      UY = UY - FI*YT*B(L)
      XT = XT2
      YT = YT2
   60 CONTINUE
      SUM = UX**2 + UY**2
      IF( SUM .EQ. 0.D0 )  GO TO 100
      DX = (V*UY - U*UX)/SUM
      X = X + DX
      DY = -(U*UY + V*UX)/SUM
      Y = Y + DY
      IF( DABS(DY)+DABS(DX) .LT. 1.0D-10 )  GO TO 80
C
      ICT = ICT+1
      IF( ICT .LT. 500 )  GO TO 50
      ISW = 1
      IF( IN .GE. 5 )  GO TO 70
      IF( IFIT .NE. 0 )  GO TO 80
   65 ISW = 0
      IFIT = 0
      GO TO 30
C
   70 IF( IFIT .EQ. 0 )  GO TO 300
      JSW = 1
   80 DO 90  L=1,KXX
      J = KJI-L+1
      TEM = A(J)
      A(J) = B(L)
cxx   90 B(L) = TEM
      B(L) = TEM
   90 CONTINUE
      ITEM = K
      K = KX
      KX = ITEM
  100 IF( IFIT .EQ. 0 )  GO TO 40
      IF( JSW .EQ. 1 )   GO TO 110
cc      IF( ISW-1 )  120,65,120
      IF( ISW-1 .LT. 0 )  GO TO 120
      IF( ISW-1 .EQ. 0 )  GO TO 65
      IF( ISW-1 .GT. 0 )  GO TO 120
  110 X = XPR
      Y = YPR
      ISW = 0
      JSW = 0
  120 IFIT = 0
      IF( X .EQ. 0.D0 )  GO TO 130
      IF( DABS( Y/X ) .LT. 1.0D-08 )  GO TO 150
  130 ALPH = 2*X
      SUM = X*X + Y*Y
      K = K-2
      GO TO 160
C
  140 X = 0.D0
      KX = KX-1
      KXX = KXX-1
  150 Y = 0.D0
      SUM = 0.D0
      ALPH = X
      K = K-1
  160 B(2) = B(2) + ALPH*B(1)
cxx  170 DO 180  L=2,K
      DO 180  L=2,K
cxx  180 B(L+1) = B(L+1) + ALPH*B(L) - SUM*B(L-1)
      B(L+1) = B(L+1) + ALPH*B(L) - SUM*B(L-1)
  180 CONTINUE
cc  190 ROOTI(K2) = Y
cc      ROOTR(K2) = X
  190 IF( K2.LE.M ) ROOTI(K2) = Y
      IF( K2.LE.M ) ROOTR(K2) = X
      K2 = K2+1
      IF( SUM .EQ. 0.D0 )  GO TO 200
      Y = -Y
      SUM = 0.D0
      GO TO 190
  200 IF (K .LE. 0 )  RETURN
      GO TO 20
C
  300 IER = 1
      RETURN
      E N D
