C     PROGRAM  8.1  LSAR1
      SUBROUTINE LSAR1( Y,N,LAG,NS0,NB,NF0,NNS,NN0,NN1,
c     * IIF,MS,SDS,AICS,MP,SDP,AICP,AS,MFS,SIG2S,NNF )
     * IIF,MS,SDS,AICS,MP,SDP,AICP,AS,MFS,SIG2S,NNF,IER )
C
      INCLUDE 'TSSS.h'
C
C  ...  Decomposition of time interval to stationary subintevals  ...
C
C     Inputs:
C        LAG:     Highest order od AR model
C        NS:      Basic local span
C     The following inputs are required in the subroutine READTS.
C        TITLE:   Caption of the data set
C        N:       Data length
C        Y(I):    Time series, (i=1,...,N)
C     Parameters:
C        IDEV:    Input device for time series
C        NMAX:    Adjustable dimension of Y (NMAX.GE.N)
C        KMAX,MJ2:  Adjustable dimensions
C        NF:      Number of frequencies for computing spectrum
C
cc      PARAMETER( NMAX=3000,KMAX=20,MJ2=KMAX+1,MJ1=100,NF=100,IDEV=1)
cx      PARAMETER( MJ1=100 )
cxx      IMPLICIT   REAL*8(A-H,O-Z)
cc      DIMENSION  Y(NMAX), X(MJ1,MJ2), D(MJ1), U(MJ2,MJ2)
cc      DIMENSION  AS(KMAX), A(KMAX)
cx      DIMENSION  Y(N), X(MJ1,LAG+1), U(2*(LAG+1),LAG+1)
cxx      DIMENSION  Y(N), X(3*(LAG+1),LAG+1), U(LAG+1,LAG+1)
cxx      DIMENSION  AS(LAG,NB), A(LAG)
c
cxx      DIMENSION  NNS(NB), NN0(NB), NN1(NB), IIF(NB)
cxx      DIMENSION  MS(NB), SDS(NB), AICS(NB), MP(NB), SDP(NB), AICP(NB)
cxx      DIMENSION  MFS(NB), SIG2S(NB), NNF(NB)
C
      INTEGER N, LAG, NS0, NB, NF0, NNS(NB), NN0(NB), NN1(NB), IIF(NB),
     1        MS(NB), MP(NB), MFS(NB), NNF(NB), IER
      DOUBLE PRECISION Y(N), SDS(NB), AICS(NB), SDP(NB), AICP(NB),
     1                 AS(LAG,NB), SIG2S(NB)
c local
      INTEGER I, II, IF, L, LK, MF, MJ1, NF, NS, NBLOCK
      DOUBLE PRECISION X(3*(LAG+1),LAG+1), U(LAG+1,LAG+1), A(LAG),
     1                 SIG2, AICF
c
      EXTERNAL   SETXAR
C
C     ISW = 0
cc      READ( 5,* )  LAG, NS
      NF = NF0
      NS = NS0
      MJ1 = 3*(LAG+1)
C
C  ...  Read time series  ...
C
cc      CALL  READTS( IDEV,Y,N )
C
      IF = 0
      IIF(1) = 0
      AICF = 0
      IER = 0
C
      NBLOCK = N/NS
      DO 100 II=1,NBLOCK
C
      L = NS*(II-1)
cc      WRITE(6,600) L
      IF( II.EQ.NBLOCK )  NS = N - NS*(II-1) - LAG
      LK  = L + LAG
c
      NNS(II) = NS
      NN0(II) = LK + 1
      NN1(II) = LK + NS
C
C  ...  Locally stationary time series  ...
C
cc      CALL  LOCAL( SETXAR,Y,X,U,D,LAG,L,NS,LAG,IF,MJ1,MJ2,A,MF,SIG2)
cx      CALL  LOCAL( SETXAR,Y,X,U,LAG,NF,L,NS,LAG,IF,MJ1,
cx     *  A,MF,SIG2,MS(II),SDS(II),AICS(II),MP(II),SDP(II),AICP(II),AICF)
      CALL  LOCAL( SETXAR,Y,N,X,U,LAG,NF,L,NS,LAG,IF,MJ1,A,MF,SIG2,
     *  MS(II),SDS(II),AICS(II),MP(II),SDP(II),AICP(II),AICF, IER)
c 2023/01/11
      IF (IER .NE. 0) RETURN
      IIF(II) = IF
C
      NNF(II) = NF
C
cc      IF( IF .EQ. 2 )     LK0 = LK + 1
      IF( IF.EQ.2 .AND. II.GT.1 )  THEN
C         CALL  ARMASP( AS,MFS,B,0,SIG2S,NF,SP )
C         CALL  PTLSP( Y,N,SP,NF,LK0S,LK2S )
      END IF
cc      MFS = MF
cc      SIG2S = SIG2
      MFS(II) = MF
      SIG2S(II) = SIG2
cc      LK0S = LK0
cc      LK2S = LK + NS
      DO 20 I=1,MF
cc   20 AS(I) = A(I)
cxx   20 AS(I,II) = A(I)
      AS(I,II) = A(I)
   20 CONTINUE
C
  100 CONTINUE
C      CALL  ARMASP( AS,MFS,B,0,SIG2S,NF,SP )
C      CALL  PTLSP( Y,N,SP,NF,LK0S,LK2S )
C      CALL  PLOTE
C      call  plot( 0.0,0.0,999 )
C
cc      STOP
      RETURN
cxx  600 FORMAT( 1H ,'L =',I5 )
      E N D
cc      SUBROUTINE  LOCAL( SETX,Z,X,U,D,LAG,N0,NS,K,IF,MJ1,MJ2,
cc     *                   A,MF,SDF )
cx      SUBROUTINE  LOCAL( SETX,Z,X,U,LAG,NF,N0,NS,K,IF,MJ1,
cx     *   A,MF,SDF,MS,SDS,AICS,MP,SDP,AICP,AICF)
      SUBROUTINE  LOCAL( SETX,Z,N,X,U,LAG,NF,N0,NS,K,IF,MJ1,
     *   A,MF,SDF,MS,SDS,AICS,MP,SDP,AICP,AICF,IER)
C
C  ...  Locally stationary AR model  ...
C
C     Inputs:
C        SETX:    Name of the subroutine for making X(I,J)
C        Z(I):    Data vector
C        D,U:     Working area
C        LAG:     Highest order of the model
C        N0:      Time point of the previous set ofobservations
C        NS:      Number of new observations
C        MJ1,MJ2: Adjustable dimension
C     Output:
C        A(I):    AR coefficients of the current model
C        MF:      Order of the current model
C        SDF:     Innovation variance of the current model
C
cxx      IMPLICIT  REAL*8 (A-H,O-Z)
cx      DIMENSION  Z(1)
cc      DIMENSION  X(MJ1,1), U(MJ2,1), D(1), A(1)
cx      DIMENSION  X(MJ1,K+1), U(2*(LAG+1),LAG+1), A(LAG)
cc      DIMENSION  B(20), AA(20,20), AIC(0:20), SIG2(0:20)
cxx      DIMENSION  Z(N)
cxx      DIMENSION  X(MJ1,K+1), U(LAG+1,LAG+1), A(LAG)
cxx      DIMENSION  B(LAG), AA(LAG,LAG), AIC(0:LAG), SIG2(0:LAG)
C
      INTEGER N, LAG, NF, N0, NS, K, IF, MJ1, MF, MS, MP, IER
      DOUBLE PRECISION Z(N), X(MJ1,K+1), U(LAG+1,LAG+1), A(LAG), SDF,
     1                 SDS, AICS, SDP, AICP, AICF
c local
      INTEGER I, K1, K2, L, MJ2, NN0, NN1, NP
      DOUBLE PRECISION B(LAG), AA(LAG,LAG), AIC(0:LAG), SIG2(0:LAG),
     1                 AIC0
C
      EXTERNAL   SETX
C
      K1 = K + 1
      K2 = K1*2
C
cx      MJ2 = K2
      MJ2 = K1
C
      NN0 = N0 + LAG + 1
      NN1 = N0 + LAG + NS
cc      WRITE(6,600)  NN0, NN1
C
cc      CALL  REDUCT( SETX,Z,D,NS,N0,K,MJ1,X )
cc      CALL  REGRES( X,K,NS,MJ1,20,AA,SIG2,AIC,MS )
c 2023/01/11
      L = MIN0(NS, MJ1)
      if (L .GE. K+1) THEN

      CALL  REDUCT( SETX,Z,NS,N0,K,MJ1,X )
      CALL  REGRES( X,K,NS,MJ1,AA,SIG2,AIC,MS )
C
      SDS = SIG2(MS)
      DO 10 I=1,MS
cxx   10 B(I) = AA(I,MS)
      B(I) = AA(I,MS)
   10 CONTINUE
      IF( IF.EQ.0 )  THEN
      CALL  COPY( X,K1,0,0,MJ1,MJ2,U )
      AICS = AIC(MS)
      AIC0 = AIC(MS)
cc      WRITE(6,610)  NS, MS, SDS, AICS
c-----
      AICP = 0.0d0
      SDP = 0.0d0
c-----
      ELSE
C
      AICS = AIC(MS) + AICF
      AIC0 = AIC(MS)
cc      WRITE(6,620)  NF, NS, MS, SDS, AICS
      CALL  COPY( X,K1,0,K2,MJ1,MJ1,X )
      CALL  COPY( U,K1,0,K1,MJ2,MJ1,X )
cc      CALL  HUSHLD( X,D,MJ1,K2,K1 )
      CALL  HUSHLD( X,MJ1,K2,K1 )
      NP = NF + NS
cc      CALL  REGRES( X,K,NP,MJ1,20,AA,SIG2,AIC,MP )
      CALL  REGRES( X,K,NP,MJ1,AA,SIG2,AIC,MP )
C
      AICP = AIC(MP)
      SDP  = SIG2(MP)
      DO 20 I=1,MP
cxx   20 A(I) = AA(I,MP)
      A(I) = AA(I,MP)
   20 CONTINUE
cc      WRITE(6,630)  NP, MP, SDP, AICP
      IF( AICS.GE.AICP )  GO TO 40
cc      WRITE(6,640)
      CALL  COPY( X,K1,K2,0,MJ1,MJ2,U )
      END IF
C
      IF = 2
      NF = NS
      MF = MS
      AICF = AIC0
      DO 30  I=1,MF
cxx   30 A(I) = B(I)
      A(I) = B(I)
   30 CONTINUE
      SDF = SDS
C
      GO TO 50
   40 IF = 1
      CALL  COPY( X,K1,0,0,MJ1,MJ2,U )
cc      WRITE(6,650)
      SDF = SDP
      MF = MP
      AICF = AICP
      NF = NF + NS
   50 CONTINUE
C
      ELSE
         IER = -1
      END IF
      RETURN
cxx  600 FORMAT( 1H0,'<<< NEW DATA (N =',I4,' ---',I4,')  >>>' )
cxx  610 FORMAT( 1H ,'  INITIAL MODEL: NS =',I4,/,10X,'MS =',I2,3X,
cxx     1  'SDS =',D13.6,3X,'AICS =',F12.3 )
cxx  620 FORMAT( 1H ,'  SWITCHED MODEL: (NF =',I4,', NS =',I4,1H),
cxx     2  /,10X,'MS =',I2,3X,'SDS =',D13.6,3X,'AICS =',F12.3 )
cxx  630 FORMAT( 1H ,'  POOLED MODEL:   (NP =',I4,1H),
cxx     3  /,10X,'MP =',I2,3X,'SDP =',D13.6,3X,'AICP =',F12.3 )
cxx  640 FORMAT( 1H ,30X,'***  SWITCHED MODEL ACCEPTED  ***' )
cxx  650 FORMAT( 1H ,30X,'***   POOLED MODEL ACCEPTED   ***' )
      E N D
      SUBROUTINE  COPY( X,K,II,JJ,MJ1,MJ2,Y )
C
C  ...  Make a copy of X on Y
C
cxx      IMPLICIT  REAL*8( A-H,O-Z )
cc      DIMENSION  X(MJ1,1) , Y(MJ2,1)
cxx      DIMENSION  X(MJ1,K) , Y(MJ2,K)
C
      INTEGER K, II, JJ, MJ1, MJ2 
      DOUBLE PRECISION X(MJ1,K) , Y(MJ2,K)
c local
      INTEGER I, I1, I2, J
C
cxx      DO 10  I=1,K
      DO 11  I=1,K
      I1 = I + II
      I2 = I + JJ
      DO 10  J=1,K
cxx   10 Y(I2,J) = X(I1,J)
      Y(I2,J) = X(I1,J)
   10 CONTINUE
   11 CONTINUE
C
      RETURN
      E N D

