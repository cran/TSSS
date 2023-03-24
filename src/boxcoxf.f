C     PROGRAM  4.4   BOXCOX
      SUBROUTINE BOXCOXF(Y, N, AICZT, FFZT, AICZ, FFZ, ZMEAN, ZVAR, ZZ)
C
      INCLUDE 'TSSS.h'
C
C  ...  Box-Cox transformation
C
C     The following inputs are required in READTS
C        TITLE:   title of the data set
C        N:       data length
C        Y(I):    data
C     12/27/90 Y.I. 8/7/91 G.K.
C
cc      PARAMETER(MJ=1000)
cxx      IMPLICIT REAL*8(A-H,O-Z )
cc      CHARACTER  TITLE*72
cc      DIMENSION  Y(MJ), Z(MJ)
cc      COMMON  /CMDATA/  TITLE
cxx      DIMENSION Y(N), Z(N), ZZ(N)
cxx      DIMENSION AICZT(21), FFZT(21), AICZ(21), FFZ(21)
cxx      DIMENSION ZMEAN(21), ZVAR(21)
C
      INTEGER N
      DOUBLE PRECISION Y(N), AICZT(21), FFZT(21), AICZ(21), FFZ(21),
     1                 ZMEAN(21), ZVAR(21), ZZ(N)
c local
      DOUBLE PRECISION Z(N), YMEAN, YVAR, FFY, AICY, A, ZJACOB, AICM
C
cc      CALL  READTS( 1,Y,N )
      CALL  GAUSSM( Y,N,YMEAN,YVAR,FFY,AICY )
cc      WRITE(6,600)
cc      WRITE(6,610)  TITLE
cc      WRITE(6,620)
      I = 0
      AICM=AICZT(1)
      DO 200 II=10,-10,-1
         I = I+1
      A = II/10.0D0
      CALL  BOXCOX( Y,N,A,Z,ZJACOB )
C
cc      CALL  GAUSSM( Z,N,ZMEAN,ZVAR,FFZ,AICZ )
      CALL  GAUSSM( Z,N,ZMEAN(I),ZVAR(I),FFZ(I),AICZ(I) )
C
cc      FFZT = FFZ + ZJACOB
cc      AICZT = AICZ-2*ZJACOB
      FFZT(I) = FFZ(I) + ZJACOB
      AICZT(I) = AICZ(I)-2*ZJACOB
cc      WRITE(6,630)  A, AICZT, FFZT, AICZ, FFZ, ZMEAN, ZVAR
c-----
      IF( I.EQ.1 ) AICM=AICZT(1)
      IF( AICZT(I).LE.AICM ) THEN
         DO 100 J=1,N
cxx  100    ZZ(J) = Z(J)
         ZZ(J) = Z(J)
  100    CONTINUE
         AICM = AICZT(I)
      END IF
c-----
  200 CONTINUE
cc      STOP
      RETURN
cxx  600 FORMAT( 1H ,'PROGRAM 4.4:   BOX-COX TRANSFORMATION' )
cxx  610 FORMAT( 1H ,A72 )
cxx  620 FORMAT( 1H ,'LAMBDA',5X,'AIC''',8X,'LL''',7X,'AIC',9X,'LL',
cxx     *            9X,'MEAN',9X,'VARIANCE' )
cxx  630 FORMAT( 1H ,F5.2,4F11.2,2D15.6 )
      E N D
      SUBROUTINE BOXCOX( Y,N,A,Z,ZJACOB )
C
C  ...  Box-Cox transformation:  Z = (Y**A - 1)/A  ...
C
C     Inputs:
C        Y(I):   data
C        N:      data length
C        A:      lambda
C     Outputs:
C        Z(I):   transformed data
C        ZJACOB: log(Jacobian), log |dz/dy|
C
cxx      IMPLICIT REAL*8(A-H,O-Z )
cxx      DIMENSION  Y(N), Z(N)
C
      INTEGER N
      DOUBLE PRECISION Y(N), A, Z(N), ZJACOB
c local
      DOUBLE PRECISION SUM
C
      SUM = 0.0D0
      DO 10 I=1,N
      IF( A .NE. 0.0D0 ) THEN
         Z(I) = (Y(I)**A - 1)/A
         SUM = SUM + (A-1)*DLOG( DABS( Y(I) ) )
      ELSE
         Z(I) = DLOG( Y(I) )
         SUM = SUM - DLOG( DABS( Y(I) ) )
      END IF
   10 CONTINUE
      ZJACOB = SUM
C
      RETURN
      E N D
      SUBROUTINE GAUSSM( Y,N,YMEAN,YVAR,FF,AIC )
C
C  ...  This subroutine fits Gaussian distribution model  ...
C
C     Inputs:
C        Y(I):   data
C        N:      data length
C     Outputs:
C        YMEAN:  mean
C        YVAR:   variance
C        FF:     log-likelihood
C        AIC:    AIC = -2FF + 4
C
cxx      IMPLICIT REAL*8( A-H,O-Z )
cxx      DIMENSION Y(N)
C
      INTEGER N
      DOUBLE PRECISION Y(N), YMEAN, YVAR, FF, AIC
c local
      DOUBLE PRECISION PI, SUM
C
      DATA  PI/3.1415926535D0/
C
      SUM = 0.0D0
      DO 10 I=1,N
cxx   10 SUM = SUM + Y(I)
      SUM = SUM + Y(I)
   10 CONTINUE
      YMEAN = SUM/N
C
      SUM = 0.0D0
      DO 20 I=1,N
cxx   20 SUM = SUM + (Y(I)-YMEAN)**2
      SUM = SUM + (Y(I)-YMEAN)**2
   20 CONTINUE
      YVAR = SUM/N
      FF  = -0.5D0*N*(DLOG( 2*PI*YVAR ) + 1)
      AIC = -2*FF + 4
C
      RETURN
      E N D
