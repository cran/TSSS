C     PROGRAM 4.1  DENSTY
      SUBROUTINE DENSTYF( MODEL, PARAM, XMIN, XMAX, K, F ) 
C
      INCLUDE 'TSSS.h'
C
C  ...  This program draws probability density function  ...
C
C     The following inputs are required in the main program:
C        MODEL:   function number (1: GAUSS, 2: CAUCHY, etc.)
C        XMIN:    lower bound of the interval
C        XMAX:    upper bound of the interval
C        PARAM:   parameter vector
C     @TEST.PN41: DEC.26,1990, 8/6/91
C
cc      PARAMETER( K=201 )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      CHARACTER*24   TITLE1(0:7)
cc      CHARACTER*72   TITLE
cxx      DIMENSION  F(K), PARAM(3), NP(0:7)
cc      COMMON  /CMDATA/  TITLE
C
      INTEGER MODEL, K
      DOUBLE PRECISION PARAM(3), XMIN, XMAX, F(K)
c local
      INTEGER NP(0:7)
      DOUBLE PRECISION USERF, GAUSS, CAUCHY, PEARSN, EXPNTL, CHISQR,
     1                 DBLEXP, UNIFRM
C
      EXTERNAL  USERF
      EXTERNAL  GAUSS
      EXTERNAL  CAUCHY
      EXTERNAL  PEARSN
      EXTERNAL  EXPNTL
      EXTERNAL  CHISQR
      EXTERNAL  DBLEXP
      EXTERNAL  UNIFRM
C
cc     DATA TITLE1/'USER SUPPLIED FUNCTION  ','NORMAL DISTRIBUTION     '
cc    *          ,'CAUCHY DISTRIBUTION     ','PEARSON DISTRIBUTION    '
cc    *          ,'EXPONENTIAL DISTRIBUTION','CHI-SQUARE DISTRIBUTION '
cc    *          ,'DOUBLE EXPONENTIAL DIST.','UNIFORM DISTRIBUTION    '/
      DATA  NP/0,2,2,3,1,1,0,2/
C
cc      WRITE(6,*)  'INPUT MODEL NUMBER ?'
cc      READ(5,*)   MODEL
cc      WRITE(6,*)  'INPUT XMIN AND XMAX ?'
cc      READ(5,*)   XMIN, XMAX
cc      IF( MODEL.EQ.0 )  THEN
cc         WRITE(6,*)  'INPUT THE NUMBER OF PARAMETERS'
cc         READ(5,*)   NP(0)
cc      END IF
cc      WRITE(6,*)  'INPUT PARAM(I),I=1,',NP(MODEL),' ?'
cc      READ(5,*)   (PARAM(I),I=1,NP(MODEL))
cc      TITLE = TITLE1(MODEL)
cc      DO 10 I=1,6
cc   10 TITLE = TITLE//'        '
C
      IF(MODEL.EQ.1)  CALL  DENSTY( GAUSS ,F,K,PARAM,XMIN,XMAX )
      IF(MODEL.EQ.2)  CALL  DENSTY( CAUCHY,F,K,PARAM,XMIN,XMAX )
      IF(MODEL.EQ.3)  CALL  DENSTY( PEARSN,F,K,PARAM,XMIN,XMAX )
      IF(MODEL.EQ.4)  CALL  DENSTY( EXPNTL,F,K,PARAM,XMIN,XMAX )
      IF(MODEL.EQ.5)  CALL  DENSTY( CHISQR,F,K,PARAM,XMIN,XMAX )
      IF(MODEL.EQ.6)  CALL  DENSTY( DBLEXP,F,K,PARAM,XMIN,XMAX )
      IF(MODEL.EQ.7)  CALL  DENSTY( UNIFRM,F,K,PARAM,XMIN,XMAX )
      IF(MODEL.EQ.0)  CALL  DENSTY( USERF ,F,K,PARAM,XMIN,XMAX )
C      CALL  PTDENS( F,K,XMIN,XMAX,PARAM,NP(MODEL) )
cc      CALL  PRDENS( F,K )
cc      STOP
      RETURN
cxx  600 FORMAT( 1H ,10F8.4 )
      E N D
      SUBROUTINE  DENSTY( DIST,P,K,PARAM,XMIN,XMAX )
C
C  ...  This subroutine evaluates values of density  ...
C               DIST(X), X=XMIN,...,XMAX
C     Inputs:
C        DIST:    name of function
C        PARAM:   parameters of the density
C        XMIN:    minimum of X
C        XMAX:    maximum of X
C        K:       number of location, I-th location is XMIN + I*DX
C                 where  DX = (I-1)*(XMAX-XMIN)/(K-1)
C     OUTPUT:
C        P(I):    density of DIST at I-th location
C
cxx      IMPLICIT REAL*8( A-H,O-Z )
cx      DIMENSION  P(K), PARAM(*)
cxx      DIMENSION  P(K), PARAM(3)
C
      INTEGER K
      DOUBLE PRECISION P(K), PARAM(3), XMIN, XMAX
c local
      INTEGER I
      DOUBLE PRECISION DX, X, DIST
C
      EXTERNAL  DIST
C
      DX = (XMAX-XMIN)/(K-1)
      DO 10 I=1,K
      X = XMIN + DX*(I-1)
cxx   10 P(I) = DIST( X,PARAM )
      P(I) = DIST( X,PARAM )
   10 CONTINUE
      RETURN
      E N D
      DOUBLE PRECISION FUNCTION  EXPNTL( X,PARAM )
C
C  ...  Exponential  distribution  ...
C
C     Inputs:
C        X:
C        PARAM(1):  lambda
C     Output:
C        EXPNTL:    density at X
C
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cxx      DIMENSION  PARAM(1)
      DOUBLE PRECISION X, PARAM(1)
C
      IF( X.GE.0.0D0 )  EXPNTL = PARAM(1)*DEXP( -PARAM(1)*X )
      IF( X.LT.0.0D0 )  EXPNTL = 0.0D0
      RETURN
      E N D
      DOUBLE PRECISION FUNCTION  CHISQR( X,PARAM )
C
C  ...  Chi-square  distributions  ...
C
C     Inputs:
C        X:
C        PARAM(1):  degree of freedoms, k
C     Output:
C        CHISQR:    density at X
C
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cx      DIMENSION  PARAM(*)
cxx      DIMENSION  PARAM(1)
      DOUBLE PRECISION X, PARAM(1), dgammafn
C
      CHISQR = 0.0D0
      IF( X.GT.0.0D0 ) CHISQR = DEXP( -X/2 )*(X/2)**(PARAM(1)/2-1.D0)
CXX     *                           /(2*DGAMMA(PARAM(1)/2))
     *                           /(2*dgammafn(PARAM(1)/2))
cxx      IF( X.LE.0.0D0 ) CHISQR = 0.0D0
      RETURN
      E N D
      DOUBLE PRECISION FUNCTION  UNIFRM( X,PARAM )
C
C  ...  Uniform distribution on [a,b]  ...
C
C     Inputs:
C        X:
C        PARAM(1):  a
C        PARAM(2):  b
C     Output:
C        UNIFRM:    density at X
C
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cxx      DIMENSION  PARAM(2)
      DOUBLE PRECISION X, PARAM(2)
C
      IF( X.GT.PARAM(1) .AND. X.LE.PARAM(2) )  THEN
         UNIFRM = 1.0D0/(PARAM(2)-PARAM(1))
      ELSE
         UNIFRM = 0.0D0
      END IF
      RETURN
      E N D
      DOUBLE PRECISION FUNCTION  USERF( X,PARAM )
C
C  ...  User supplied density function  ...
C       (The following is an example of two-sided exponential dist.)
C
C     Inputs:
C        X:
C        PARAM(1):  lambda
C     Output:
C        USERF:     density at X
C
cxx      IMPLICIT  REAL*8(A-H,O-Z)
cx      DIMENSION  PARAM(2)
cxx      DIMENSION  PARAM(1)
      DOUBLE PRECISION X, PARAM(1)
C
      IF( X.GE.0.0D0 )  THEN
         USERF = PARAM(1)*DEXP( -PARAM(1)*X )/2
      ELSE
         USERF = PARAM(1)*DEXP(  PARAM(1)*X )/2
      END IF
      RETURN
      E N D

