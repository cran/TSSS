C     PROGRAM 13.3  TVSPC
      SUBROUTINE TVSPC( N, M, NOBS, NF, IVAR, SIG2, AR, VAR, SP )
C
      INCLUDE 'TSSS.h'
C
C  ...  CHANGING SPECTRUM ...
C
C     INPUTS:
C        N:    NUMBER OF POINTS
C        M:    AR ORDER
C        K:    TREND ORDER
C        NOBS: SAMPLING INTERVAL
C        IVAR=1:  FOR VARIANCE CORRECTION
C        AR(I,J): TIME VARYING AR COEFFICIENT
C        VAR(I):  TIME VARYING VARIANCE
C     PARAMETRS:
C        NMAX: ADJUSTABLE DIMENSION OF AR (NMAX >= N)
C        MJ:   ADJUSTABLE DIMENSION OF AR (MJ   >= M)
C        NF:   NUMBER OF FREQUENCIES
C     @TEST.UNI:  5/16/89  1/07/91  9/09/92
C
cc      PARAMETER( NMAX=200,NF=200,MJ=10,IDEV1=2,IDEV2=3 )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      CHARACTER  TITLE*72,  VNAME*8
cc      DIMENSION  AR(MJ,NMAX), VNAME(5), VALUE(5), VAR(3000)
cxx      DIMENSION  AR(M,N), VAR(N*NOBS), SP(0:NF,N)
cc      COMMON  /CMDATA/  TITLE
cc      DATA    VAR /3000*1.0D0/
C
      INTEGER N, M, NOBS, NF, IVAR
      DOUBLE PRECISION SIG2, AR(M,N), VAR(N*NOBS), SP(0:NF,N)
C
cc      OPEN( IDEV1,FILE='A:$LASERLIB$TEMP3.DAT' )
cc      READ(IDEV1,1) TITLE
cc      READ(IDEV1,*)  N, M, K, NOBS, IVAR
cc      READ(IDEV1,*)  SIG2
cc      DO 10 J=1,N
cc   10 READ(IDEV1,*) (AR(I,J),I=1,M)
cc      CLOSE( IDEV1 )
      NN = N*NOBS
cc      IF(IVAR.EQ.1 )  OPEN( IDEV2,FILE='A:$LASERLIB$TEMP4.DAT' )
cc      IF(IVAR.EQ.1 )  READ( IDEV2,*) (VAR(I),I=1,NN)
cc      IF(IVAR.EQ.1 )  CLOSE( IDEV2 )
cc    1 FORMAT( A72 )
      IF( IVAR.NE.1 )  THEN
cxx         DO 100 I = 1,NN
cxx  100    VAR(I) = 1.0D0
         VAR(1:NN) = 1.0D0
      END IF
C
C  ...  PLOT CHANGING SPECTRUM  ...
C
C
cc      VNAME(1) = 'N      ='
cc      VNAME(2) = 'M      ='
cc      VNAME(3) = 'K      ='
cc      VNAME(4) = 'NOBS   ='
cc      VNAME(5) = 'SIG2   ='
cc      VALUE(1) = N
cc      VALUE(2) = M
cc      VALUE(3) = K
cc      VALUE(4) = NOBS
cc      VALUE(5) = SIG2
cc      CALL  PLOTS
C     call  plots( 1,0,0,1,0 )
C     call  form( 1 )
C     call  factor( 10.0 )
cc      CALL  HEADER( 'PROGRAM 12.2:  CHANGING SPECTRUM',32,5,VNAME,
cc     *              VALUE )
cc      CALL  PLOT( 3.0,2.0,-3 )
cc      CALL  PT3DSP( AR,SIG2,MJ,M,N,NOBS,VAR )
      CALL  PT3DSP( AR,SIG2,M,N,NOBS,NF,VAR,SP )

cc      CALL  PLOTE
C     call  plot( 0.0,0.0,999 )
      RETURN
      E N D
cc      SUBROUTINE  PT3DSP( A,SIG2,MJ,M,N,NOBS,VAR )
      SUBROUTINE  PT3DSP( A,SIG2,M,N,NOBS,NF,VAR,SP )
C
C  ...  COMPUTE AND DRAW TIME VARYING SPECTRUM  ...
C
C     INPUTS:
C        A,SIG2,M:  AR MODELS
C        N:         DETA LENGTH
C        NOBS:      LOCAL STATIONARY SPAN
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION A(MJ,N), SP(0:200), YH(2000), VAR(3000), B(20)
cxx      DIMENSION A(M,N), SP(0:NF,N), VAR(N*NOBS)
C
      INTEGER M, N, NOBS, NF
      DOUBLE PRECISION A(M,N), SIG2, VAR(N*NOBS), SP(0:NF,N)
c local
      DOUBLE PRECISION B
C
cc      ANGLE = 70.0
cc      NF = 200
C
cc      WX = 5.0
cc      WY = 3.0
cc      WZ =10.0
cc      X0 = 0.0
cc      X1 = 0.5
cc      DX = 0.5
cc      Z0 =  0.0
cc      Z1 =  N*NOBS
cc      DZ =  500.0
cc      K = 2
cc      CALL  ARMASP( A,M,B,0,SIG2,NF,SP )
      CALL  ARMASP( A,M,B,0,SIG2,NF,SP(0,1) )
cc      CALL  MAXMIN( SP(0),NF+1,Y0,Y1,DY )
C
C  ...  DRAW 3D-AXIS  ...
C
cc      CALL  PLOT3A( NF,M,WX,WY,WZ,ANGLE,X0,X1,DX,Y0,Y1,DY,Z0,Z1,
cc     *              DZ,KN,K,YH )
      DO 100 J=1,N
C
C  ...  SPECTRUM OF AR MODEL AT TIME J*NOBS   ...
C
cc      CALL  ARMASP( A(1,J),M,B,0,SIG2,NF,SP )
      CALL  ARMASP( A(1,J),M,B,0,SIG2,NF,SP(0,J) )
C
C  ...  DRAW SPECTRUM AND SHIFT THE ORIGIN  ...
C
      DO 10 I=0,NF
cc   10 SP(I) = SP(I) + DLOG10( VAR(J*NOBS-NOBS/2) )
cxx   10 SP(I,J) = SP(I,J) + DLOG10( VAR(J*NOBS-NOBS/2) )
      SP(I,J) = SP(I,J) + DLOG10( VAR(J*NOBS-NOBS/2) )
   10 CONTINUE
cc      CALL  PLOT3B( SP(0),NF+1,N,WX,WY,WZ,ANGLE,Y0,Y1,K,YH,ZS,ZC )
  100 CONTINUE
C      CALL  PLOT( -SNGL(ZC),-SNGL(ZS),-3 )
C
      RETURN
      E N D
     
