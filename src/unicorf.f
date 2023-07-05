C     PROGRAM 2.1  UNICOR
      SUBROUTINE UNICORF(Y,N,LAG,OUTMIN,OUTMAX,COV,YMEAN)
C
      INCLUDE 'TSSS.h'
C
C  ...  This program computes sample autocovariance and sample
C       autocorrelation function  ...
C
C     The following inputs are required in the subroutine READTS.
C        TITLE:   caption of the data set
C        N:       data length
C        FORMAT:  reading format
C        Y(I):    time series, (i=1,...,N)
C     Parameters:
C        NMAX:    adjustable dimension of Y (NMAX.GE.N)
C        LAG:     maximum lag of autocovariance function
C        IDEV:    input file specification
C        JDEV:    output file specification
C     5/16,1989 modified 12/05/90, 12/06/90, 7/18/92
C
cc      PARAMETER( NMAX=1000,LAG=30,IDEV=1,JDEV=6 )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  Y(NMAX), COV(0:LAG,4)
cxx      DIMENSION  Y(N), COV(0:LAG,4)
cc      DATA  OUTMIN/-1.0D30/, OUTMAX/ 1.0D30/
C
      INTEGER N, LAG
      DOUBLE PRECISION Y(N), OUTMIN, OUTMAX, COV(0:LAG,4), YMEAN
C
C  ...  read in data  ...
C
cc      CALL  READTS( IDEV,Y,N )
C
C  ...  compute autocovariance, autocorrelation and their error bound  .
C
cc      CALL  UNICOR( Y,N,LAG,OUTMIN,OUTMAX,COV )
      CALL  UNICOR( Y,N,LAG,OUTMIN,OUTMAX,COV,YMEAN )
C
C  ...  print out results  ...
C
cc      CALL  PRACOR( JDEV,COV,N,LAG )
cc      STOP
      RETURN
      E N D
