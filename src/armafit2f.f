C     PROGRAM 10.1  ARMAFT
      SUBROUTINE ARMAFT2( Y0,N,MMAX,LMAX,MLMAX,TSIG2,TFF,TAIC,AR,CMA,
     1                    IER )
C
      INCLUDE 'TSSS.h'
C
C  ...  ARMA MODEL FITTING  ...
C
C     Inputs:
C        M:       AR Order (M <= 20)
C        L:       MA Order (M <= 20)
C        IPARAM:  =0    Use defalt initail values
C                 =1    Read intial values
C        AR(I):   AR coefficients (I=1,M3)
C     Parameters:
C        NMAX:    Adjustable dimension of Y
C        MJ:      Adjustable dimension of XF, VF, etc.
C     @TEST.FILTER2O    NOV.29,1990, SEP.02,1992
C
cc      PARAMETER( NMAX=1000,MJ=40 )
cc      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  AA(40), AR(40), PAR(40), CMA(40)
cc      DIMENSION  TAIC(0:20,0:20), TFF(0:20,0:20)
cc      DIMENSION  SAA(40,0:20,0:20)
cc      COMMON  /C92825/  OUTMIN, OUTMAX
cc      COMMON  /C92826/  Y(NMAX)
cc      COMMON  /C92907/  ALIMIT
cc      COMMON  /C92908/  M, L, N
cc      COMMON  /C92909/  FLK, SIG2
cc      COMMON  / CCC /  ISW, IPR, ISMT, IDIF
c
      INTEGER N, MMAX, LMAX, MLMAX, IER(3)
      DOUBLE PRECISION Y0(N), TSIG2(0:MMAX, 0:LMAX), TFF(0:MMAX,0:LMAX),
     1                 TAIC(0:MMAX, 0:LMAX), AR(MMAX, 0:MMAX, 0:LMAX),
     2                 CMA(LMAX, 0:MMAX, 0:LMAX)
c local
      INTEGER I, J, M, L, IPRAM, IOPT, IDIF, NSUM
      DOUBLE PRECISION Y(N), PI, SUM, YMEAN, YVAR, ALIMIT, OUTMIN,
     1                 OUTMAX, PAR(MLMAX), AA(MMAX+LMAX), SIG2, FLK,
     2                 AIC, SAA(MMAX+LMAX, 0:MMAX, 0:LMAX)
c
      EXTERNAL  FFARMA
      PI = 3.1415926530D0
cc      IPR = 0
c
      TSIG2(0:MMAX, 0:LMAX) = 0.0D0
      TFF(0:MMAX,0:LMAX) = 0.0D0
      TAIC(0:MMAX, 0:LMAX) = 0.0D0
      AR(1:MMAX, 0:MMAX, 0:LMAX) = 0.0D0
      CMA(1:LMAX, 0:MMAX, 0:LMAX) = 0.0D0
      PAR(1:MLMAX) = 0.0D0
      AA(1:(MMAX+LMAX)) = 0.0D0
      SAA(1:(MMAX+LMAX), 0:MMAX, 0:LMAX) = 0.0D0
      IER(1:3) = 0
C
C  ...  Read Time Series  ...
C
cc      CALL  READTS( 1,Y,N )
      DO 5 I = 1,N
      Y(I) = Y0(I)
    5 CONTINUE
C
C  ...  Subtrac Mean Value  ...
C
cc      SUM = 0.0D0
cc      DO 10 I=1,N
cccc   10 SUM = SUM + Y(I)
cc      SUM = SUM + Y(I)
cc   10 CONTINUE
cc      YMEAN = SUM/N
      CALL  MEAN( Y,N,-1.0D30,1.0D30,NSUM,YMEAN )
      DO 20 I=1,N
cc   20 Y(I) = Y(I) - YMEAN
      Y(I) = Y(I) - YMEAN
   20 CONTINUE

cc      DO 25 I=0,20
cc      DO 25 J=0,20
cc   25 TFF(I,J) = -1.0D10
      DO 26 I=0,MMAX
      DO 25 J=0,LMAX
      TFF(I,J) = -1.0D10
   25 CONTINUE
   26 CONTINUE
C
C  ...  Read Model Orders  ...
C
C      READ( 5,* )  M, L, IPARAM
cc      open( 3,file='temp.dat' )
cc      MMAX = 10
ccc     LMAX = 10
cc      DO 100  M=0,MMAX
      DO 101  M=0,MMAX
      DO 100  L=0,LMAX
      
cc      write(3,*) M, L
      IF( M.EQ.0 .and. L.EQ.0 )  THEN
cc         SUM = 0.0D0
cc         DO 110 I=1,N
cccc  110    SUM = SUM + Y(I)**2
cc         SUM = SUM + Y(I)**2
cc  110    CONTINUE
         CALL  MEAN2( Y,N,-1.0D30,1.0D30,SUM )
         YVAR = SUM/N
         TSIG2(0,0) = YVAR
         TFF(0,0) = -0.5D0*N*(DLOG(2*PI*YVAR) + 1)
         TAIC(0,0) = -2*TFF(0,0) + 2
         GO TO 100
      END IF
C
C  ...  Set Defalt Parameters  ...
C
cc      IPR   = 0
      IPRAM = 0
      IOPT  = 1
      IDIF  = 2
      ALIMIT = 0.95D0
      
      IF( L.EQ.0 .OR. M.EQ.0 ) THEN 
      
cc      CALL  SPARA1( M,L,AR,CMA,OUTMIN,OUTMAX,IOPT )
      CALL  SPARA1( M,L,MLMAX,AR(1,M,L),CMA(1,M,L),OUTMIN,OUTMAX,IOPT )
C      IF( IPARAM.EQ.1 )  THEN
C         READ( 5,* )  (AR(I),I=1,M)
C         READ( 5,* )  (CMA(I),I=1,L)
C      END IF
C
C      WRITE(6,*) M, L
C      WRITE(6,*) (AR(I),I=1,M)
C      WRITE(6,*) (CMA(I),I=1,L)
C
cc      CALL  PARCOR( AR,M,PAR )
      CALL  PARCOR( AR(1,M,L),M,PAR )
C      write(6,*) (PAR(I),I=1,M)
      DO 30 I=1,M
cc   30 AA(I) = DLOG( (ALIMIT+PAR(I))/(ALIMIT-PAR(I)) )
      AA(I) = DLOG( (ALIMIT+PAR(I))/(ALIMIT-PAR(I)) )
   30 CONTINUE
cc      CALL  PARCOR( CMA,L,PAR )
      CALL  PARCOR( CMA(1,M,L),L,PAR )
C      write(6,*) (PAR(I),I=1,L)
      DO 40 I=1,L
cc   40 AA(M+I) = DLOG( (ALIMIT+PAR(I))/(ALIMIT-PAR(I)) )
      AA(M+I) = DLOG( (ALIMIT+PAR(I))/(ALIMIT-PAR(I)) )
   40 CONTINUE
   
      ELSE
         AA(M+L) = 0.0000D0
C         IF( M.EQ.0.AND.L.EQ.1 )  AA(1) = 0.18D0
C         IF( M.GT.0 )  THEN
         IF( TFF(M-1,L).GT.TFF(M,L-1) )  THEN
            DO 70 I=1,M-1
cc   70       AA(I) = SAA(I,M-1,L)
            AA(I) = SAA(I,M-1,L) 
   70       CONTINUE
           AA(M) = 0.0000D0
            DO 80 I=1,L
cc   80       AA(M+I) = SAA(M-1+I,M-1,L)         
            AA(M+I) = SAA(M-1+I,M-1,L)         
   80       CONTINUE         
         END IF
C         END IF
      END IF
      
cc      WRITE(6,*) M, L
cc      WRITE(6,666) (AA(I),I=1,M+L)
cc  666 FORMAT( 10F10.5 )
  
C
C  ...  Maximum Likelihood Method  ...
C
cc      write(3,*)  (AA(I),I=1,M+L)
cc      IF( IOPT.EQ.1 )  CALL  DAVIDN( FFARMA,AA,M+L,2 )
      IF( IOPT.EQ.1 )  CALL  DAVIDN( FFARMA,AA,M+L,2,
     *    Y,N,M,L,MLMAX,OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER(1) )
      if( ier(1).ne.0 ) then
         ier(2) = m
         ier(3) = l
         return
      end if
C
C
      DO 50 I=1,M
cc   50 PAR(I)  = ALIMIT*(DEXP(AA(I))-1.0D0)/(DEXP(AA(I))+1.0D0)
      PAR(I)  = ALIMIT*(DEXP(AA(I))-1.0D0)/(DEXP(AA(I))+1.0D0)
   50 CONTINUE
cc      CALL  ARCOEF( PAR,M,AR )
      CALL  ARCOEF( PAR,M,AR(1,M,L) )
      DO 60 I=1,L
cc   60 PAR(I)  = ALIMIT*(DEXP(AA(M+I))-1.0D0)/(DEXP(AA(M+I))+1.0D0)
      PAR(I)  = ALIMIT*(DEXP(AA(M+I))-1.0D0)/(DEXP(AA(M+I))+1.0D0)
   60 CONTINUE
cc      CALL  ARCOEF( PAR,L,CMA )
      CALL  ARCOEF( PAR,L,CMA(1,M,L) )
C
      AIC = -2*FLK + 2*(M+L+1)
C
C  ...  Print out the Maximum Likelihood Estimates  ...
C
cc      WRITE(6,666) (AA(I),I=1,M+L)

cc      CALL  PRARMA( M,L,AR,CMA,SIG2,FLK,AIC )
      TFF(M,L) = FLK
      TAIC(M,L) = AIC
      TSIG2(M,L) = SIG2
      DO 90 I=1,M+L
cc   90 SAA(I,M,L) = AA(I)
      SAA(I,M,L) = AA(I)
   90 CONTINUE
cc      WRITE(3,33)  FLK, AIC
cc   33 FORMAT( 5f12.4 )
cc      WRITE(3,*)  (AA(I),I=1,M+L)
  100 CONTINUE
  101 CONTINUE
C
cc      open( 2,file='temp_table.dat' )
cc      WRITE(6,*)  'log-Likelihood(M,L)'
cc      WRITE(2,*)  'log-Likelihood(M,L)'
cc      DO 200 I=0,MMAX
cc      WRITE(2,600)  (TFF(I,J),J=0,LMAX)
cc  200 WRITE(6,600)  (TFF(I,J),J=0,LMAX)
cc      WRITE(2,*)  'AIC(M,L)'
cc      WRITE(6,*)  'AIC(M,L)'
cc      DO 210 I=0,MMAX
cc      WRITE(2,600)  (TAIC(I,J),J=0,LMAX)
cc  210 WRITE(6,600)  (TAIC(I,J),J=0,LMAX)
cc      close( 2 )
cc      close( 3 )
cc  600 FORMAT( 11F10.4 )
cc      STOP
      RETURN
      E N D
      SUBROUTINE  MEAN2( Y,N,OUTMIN,OUTMAX,YSUM )
C
C     Inputs:
C        Y(I):    time series
C        N:       data length
C        OUTMIN:  bound for outliers in low side
C        OUTMAX:  bound for outliers in high side
C     Outputs:
C        YSUM:    sum of non-outlier observations
C
      INTEGER N
      DOUBLE PRECISION Y(N), OUTMIN, OUTMAX, YSUM
c local
      INTEGER I
C
      YSUM = 0.0D0
      DO 10 I=1,N
      IF( Y(I) .GT.OUTMIN .AND. Y(I).LT.OUTMAX ) THEN
         YSUM = YSUM + Y(I)**2
      END IF
   10 CONTINUE
C
      RETURN
      E N D
