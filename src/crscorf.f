C     PROGRAM 2.3  CRSCOR
      SUBROUTINE CRSCORF( Y,N,L,LAG,OUTMIN,OUTMAX,C,R,YMEAN )
C
      INCLUDE 'TSSS.h'
C
C  ...  This program computes sample cross-covariance and
C       sample cross-correlation functions  ...
C
C     The following inputs are required in the subroutine READMD.
C        TITLE:   title of the data
C        N:       data length
C        L:       dimension of the observation
C        IFM:     = 1   ((Y(I,J),J=1,L),I=1,N)
C                 = 2   ((Y(I,J),I=1,N),J=1,L)
C        Y(I,J):  multi-variate time series
C     Parameters:
C        MJ:      adjustable dimension of Y; (MJ.GE.N)
C        MJ1:     adjustable dimension of Y, C, R; (MJ1.GE.L)
C        LAG:     maximum lag of the cross-covariance function
C        IDEV:    input device specification
C     MAR.24,1989:  modified  12/20/90, 7/18/92
C
cc      PARAMETER( MJ=1000,MJ1=7,LAG=50,IDEV=1,JDEV=6 )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  Y(MJ,MJ1), OUTMIN(10), OUTMAX(10)
cc      DIMENSION  C(0:LAG,MJ1,MJ1), R(0:LAG,MJ1,MJ1), YMEAN(10)
cc      DATA  OUTMIN/10*-1.0D30/, OUTMAX/10*1.0D30/
cxx      DIMENSION  Y(N,L), OUTMIN(L), OUTMAX(L)
cxx      DIMENSION  C(0:LAG,L,L), R(0:LAG,L,L), YMEAN(L)
C
      INTEGER N, L, LAG
      DOUBLE PRECISION Y(N,L), OUTMIN(L), OUTMAX(L), C(0:LAG,L,L),
     1                 R(0:LAG,L,L), YMEAN(L)
C
C  ...  read in multivariate time series  ...
C
cc      CALL  READMD( IDEV,MJ,Y,N,L )
C
C  ...  cross-covariance function  ...
C
cc      CALL  CRSCOR( Y,N,L,LAG,MJ,MJ1,OUTMIN,OUTMAX,C,R,YMEAN )
      CALL  CRSCOR( Y,N,L,LAG,OUTMIN,OUTMAX,C,R,YMEAN )
C
C  ...  print out and plot cross-correlation function  ...
C
cc      CALL  PRMCOR( JDEV,C,R,L,LAG,MJ1 )
C
cc      STOP
      RETURN
      E N D
