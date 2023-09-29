C     PROGRAM  3.1  PERIOD
cx      SUBROUTINE PERIODF(Y,N,IWINDW,OUTMIN,OUTMAX,PE,SPE,NP)
ccb      SUBROUTINE PERIODF(Y,N,IWINDW,OUTMIN,OUTMAX,PE,SPE,NP,IFG)
ccxx      SUBROUTINE PERIODF(Y,N,NP,IWINDW,OUTMIN,OUTMAX,PE,SPE,IFG)
      SUBROUTINE PERIODF(Y,N,NP,IWINDW,LAG,OUTMIN,OUTMAX,PE,SPE,IFG)
C
      INCLUDE 'TSSS.h'
C
C  ...  This program computes periodogram  ...
C
C     The following inputs are required in READTS.
C        TITLE:   title of the data
C        N:       data length
C        Y(I):    time series
C     Parameters:
C        NMAX:    adjustable dimension of Y (NMAX.GE.N)
C        MJ:      adjustable dimension of COV
C        NF:      adjustable dimension of PE and SPE
C        IDEV:    input device
C        JDEV:    output device
C        IWINDW:  window type (0: box-car, 1: Hanning, 2: Hamming)
C     @TEST.PN31:  5/16/89, 12/5/90, 12/21/90, 8/5/91
C
cc      PARAMETER( NMAX=200,MJ=NMAX,NF=512,IDEV=1,JDEV=6,IWINDW=1 )
ccb      PARAMETER( NF=512 )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  Y(NMAX), PE(0:NF), SPE(0:NF), COV(0:MJ)
ccb      DIMENSION  Y(N), PE(0:NF), SPE(0:NF), COV(0:1024)
cxx      DIMENSION  Y(N), PE(0:NP), SPE(0:NP), COV(0:1024)
C
      INTEGER N, NP, IWINDW, LAG, IFG
      DOUBLE PRECISION Y(N), OUTMIN, OUTMAX, PE(0:NP), SPE(0:NP)
c local
ccxx      DOUBLE PRECISION COV(0:1024)
      DOUBLE PRECISION COV(0:LAG)
C
cc      DATA  OUTMIN/-1.0D30/, OUTMAX/ 1.0D30/
C
cc      CALL  READTS( IDEV,Y,N )
C
cc      IF( IWINDW.EQ.0 )  LAG = MIN0( N-1,MJ )
ccxx      IF( IWINDW.EQ.0 )  LAG = N-1
cxx      IF( IWINDW.GT.0 )  LAG = 2*DSQRT( DBLE(N) )
ccxx      IF( IWINDW.GT.0 )  LAG = INT(2*DSQRT( DBLE(N) ))
ccxx      IF( IWINDW.EQ.0 )  NP  = (LAG+1)/2
ccxx      IF( IWINDW.GT.0 )  NP  = LAG
C
      CALL  PERIOD( Y,N,LAG,OUTMIN,OUTMAX,NP,0,COV,PE )
C
cx      CALL  WINDOW( PE,NP,IWINDW,SPE )
      CALL  WINDOW( PE,NP,IWINDW,SPE,IFG )
C
cc      CALL  PRPER( JEV,PE,SPE,N,NP,IWINDW )
C
cc      STOP
      RETURN
      E N D
      SUBROUTINE  PERIOD( Y,N,LAG,OUTMIN,OUTMAX,NN,ISW,COV,P )
C
C ... Periodogram computation by Fourier transf. of autocovariance  ...
C
C     Inputs:
C        Y(I):    time series
C        N:       data length
C        LAG:     maximum lag of autocovariance
C        NN:      number of frequencies
C        OUTMIN:  lower bound for detecting outliers
C        OUTMAX:  upper bound for detecting outliers
C        ISW = 0: autocovariance is not given
C            = 1:                   given
C        COV:     autocovariance (if ISW=1)
C     Outputs:
C        COV(I):  autocovariance (if ISW=0)
C        P(I):    periodogram (raw spectrum)
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
ccb      DIMENSION  Y(N), COV(0:1024), P(0:512), FC(513), FS(513)
cxx      DIMENSION  Y(N), COV(0:1024), P(0:NN), FC(NN+1), FS(NN+1)
C
      INTEGER N, LAG, NN, ISW
ccxx      DOUBLE PRECISION Y(N), OUTMIN, OUTMAX, COV(0:1024), P(0:NN)
      DOUBLE PRECISION Y(N), OUTMIN, OUTMAX, COV(0:LAG), P(0:NN)
c local
      INTEGER I
      DOUBLE PRECISION FC(NN+1), FS(NN+1), YMEAN
C
ccxx      IF( LAG.GE.1023 )  LAG = 1023
      IF( ISW.EQ.0 )  CALL  AUTCOV( Y,N,LAG,OUTMIN,OUTMAX,COV,YMEAN )
C
C  ...  Fourier transformation by Goertzel method  ...
C
      CALL  FOURIE( COV(0),LAG+1,NN+1,FC,FS )
C
      DO 10 I=0,NN
cxx   10 P(I) =  2.0D0*FC(I+1) - COV(0)
      P(I) =  2.0D0*FC(I+1) - COV(0)
   10 CONTINUE
C
      RETURN
      E N D

