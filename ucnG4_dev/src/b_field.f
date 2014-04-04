!------------------------------------------------------------------
      SUBROUTINE B_FIELD(BX,BY,BZ,x,y,z)
!------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,J-M,O-Z), INTEGER*4(I,N)
      DOUBLE PRECISION MAGZ(250,50)
      DOUBLE PRECISION ZA(250),PA(50),MAGP(250,50)
      DOUBLE PRECISION ZIN(4),ZIN1(4),ZIN2(4),ZIN12(4),PIN(4),
     1 PIN1(4),PIN2(4),PIN12(4)
      COMMON/FIELD/ZA,PA,MAGZ,MAGP
!-----------------------------------------------------------------
      P=DSQRT(X**2+Y**2)
      psi=DASIN(Y/P)
c      ABZ=DABS(Z)
      
c      NZ = 1 + INT(ABZ)
c      NP = 1 + INT(P)

c      CALL GETARRAYS(ZIN,ZIN1,ZIN2,ZIN12,MAGZ,ZA,PA,NZ,NP)
c      CALL GETARRAYS(PIN,PIN1,PIN2,PIN12,MAGP,ZA,PA,NZ,NP)
c
c      CALL BCUINT(ZIN,ZIN1,ZIN2,ZIN12,ZA(NZ),ZA(NZ+1),PA(NP),PA(NP+1),
c     1      ABZ,P,BZ,BZG1,BZG2)
c
c      CALL BCUINT(PIN,PIN1,PIN2,PIN12,ZA(NZ),ZA(NZ+1),PA(NP),PA(NP+1),
c     1      ABZ,P,BP,BPG1,BPG2)

      CALL CALCULATE_BZ(BZ,X,Y,Z)
      CALL CALCULATE_BP(BP,X,Y,Z)

      IF(X.EQ.0.0.AND.Y.EQ.0.0)THEN
         BX=0.0
         BY=0.0
      ELSE
         BX=BP*DCOS(psi)
         BY=BP*DSIN(psi)
         IF(X.LT.0.0)then
            BX=-1.0*DABS(BX)
         endif
         IF(Y.LT.0.0)then 
            BY=-1.0*DABS(BY)
         endif
      ENDIF
      
      IF(Z.lt.0.0) Then
        BX=-1.0*BX
        BY=-1.0*BY
      ENDIF
c
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE GETARRAYS(Y,Y1,Y2,Y12,XX,Z,P,iZ,iP)
      DOUBLE PRECISION Y(4),Y1(4),Y2(4),Y12(4),XX(250,250),Z(250),
     1 P(250)
      INTEGER iZ,iP
      DOUBLE PRECISION XXDIF1,XXDIF2,XXDIF3,XXDIF4,XXPDIF1,XXPDIF2,
     1 XXPDIF3,XXPDIF4,ZAA,PAA
c----------------------------------------------------------------------c  
C    GETARRAYS DETERMINES  THE MAGNETIC FIELD AND THE GRADIENT OF THE  C
C    FIELD ON THE FOUR GRID POINTS SURROUNDING THE ELECTRON POSTION.   C
C    THESE POINTS ARE THE INPUT TO THE BICUBIC INTERPOLATION ROUTINE.  C
C     
C     Y(i) = B_j(i) ; WHERE I IS THE IS GRID POINT, STARTING AT THE    C
C                     THE MIN(Z),MIN(P) POINT AND MOVING AROUND        C
C                     COUNTER CLOCKWISE, AND J IS THE Z OR P FIELD     C
C                     COMPONENT.                                       C
C     
C    Y1(i) = dB_j(i)/dZ ; Y1 IS THE FIELD GRADIENT IN THE Z DIRECTION. C
C                        
C    Y2(i) = dB_j(i)/dP ; Y2 IS THE FIELD GRADIENT IN THE P DIRECTION. C
C
C    Y12(i) = d^2 B_j(i) / dZdP ; Y12 IS THE DOULBE DIREVATIVE.        C
C----------------------------------------------------------------------C

      Y(1) = XX(iZ,iP)
      Y(2) = XX(iZ+1,iP)
      Y(3) = XX(iZ+1,iP+1)
      Y(4) = XX(iZ,iP+1)

      IF(iZ.EQ.1)THEN
         XXDIF1 = XX(2,iP)
         XXDIF2 = XX(2,iP+1)
         IF(iP.GT.1)THEN
            XXDIF3 = XX(2,iP-1)
         ELSE
            XXDIF3 = XX(2,2)
         ENDIF
         XXDIF4 = XX(2,iP+2) 
         ZAA = -Z(2)
      ELSEIF(iZ.GT.1)THEN
         XXDIF1 = XX(IZ-1,IP)
         XXDIF2 = XX(IZ-1,IP+1)
         IF(iP.GT.1)THEN
            XXDIF3 = XX(iZ-1,iP-1)
         ELSE
            XXDIF3 = XX(IZ-1,2)
         ENDIF
         XXDIF4 = XX(IZ-1,IP+2)
         ZAA = Z(IZ-1)
      ENDIF

      IF(IP.EQ.1)THEN
         XXPDIF1 = XX(IZ,2)
         XXPDIF2 = XX(IZ+1,2)
         IF(IZ.GT.1)THEN
           XXPDIF3 = XX(IZ-1,2)
         ELSE 
           XXPDIF3 = XX(2,2)
         ENDIF
         XXPDIF4 = XX(IZ+2,2)
         PAA = -P(2)
      ELSEIF(IP.GT.1)THEN
         XXPDIF1 = XX(IZ,IP-1)
         XXPDIF2 = XX(IZ+1,IP-1)

         IF(IZ.GT.1)THEN
             XXPDIF3 = XX(IZ-1,IP-1)
         ELSE
             XXPDIF3 = XX(2,IP-1)
         ENDIF

         XXPDIF4 = XX(IZ+2,IP-1)
         PAA = P(IP-1)
      ENDIF

      DIFS = 2.
 
      Y1(1) = (XX(IZ+1,IP)   - XXDIF1)     /DIFS !(Z(IZ+1)-ZAA)
      Y1(2) = (XX(IZ+2,IP)   - XX(IZ,IP))  /DIFS !(Z(IZ+2)-Z(IZ))
      Y1(3) = (XX(IZ+2,IP+1) - XX(IZ,IP+1))/DIFS !(Z(IZ+2)-Z(IZ))
      Y1(4) = (XX(IZ+1,IP+1) - XXDIF2)     /DIFS !(Z(IZ+1)-ZAA)

      Y2(1) = (XX(IZ,IP+1)   - XXPDIF1)    /DIFS !(P(IP+1)-PAA)
      Y2(2) = (XX(IZ+1,IP+1) - XXPDIF2)    /DIFS !(P(IP+1)-PAA)
      Y2(3) = (XX(IZ+1,IP+2) - XX(IZ+1,IP))/DIFS !(P(IP+2)-P(IP))
      Y2(4) = (XX(IZ,IP+2)   - XX(IZ,IP))  /DIFS !(P(IP+2)-P(IP))
c
c     double derivative for point nz,np
c
      Y12(1) = (XX(IZ+1,IP+1)-XXPDIF2-XXDIF2+XXDIF3)/(DIFS**2)!((Z(IZ+1)-ZAA)*
C     1              (P(IP+1)-PAA))
c
c     double derivative for point nz+1,np
c
      Y12(2) = (XX(IZ+2,IP+1) - XXPDIF4 - XX(IZ,IP+1) +
     1             XXPDIF1)/(DIFS**2)!((Z(IZ+2)-Z(IZ))*
C     1              (P(IP+1)-PAA))
c
c     double derivative for point nz+1,np+1
c
      Y12(3) = (XX(IZ+2,IP+2)-XX(IZ+2,IP)-XX(IZ,IP+2)+
     1              XX(IZ,IP))/(DIFS**2)!((Z(IZ+2)-Z(IZ))*
C     1              (P(IP+2)-P(IP)))
c
c     double derivative for point nz,np+1
c
      Y12(4) = (XX(IZ+1,IP+2)-XX(IZ+1,IP)-XXDIF4+
     1              XXDIF1)/(DIFS**2)!((Z(IZ+1)-ZAA)*
C     1              (P(IP+2)-P(IP)))

      RETURN
      END

c----------------------------------------------------------------------c
      SUBROUTINE BCUINT(Y,Y1,Y2,Y12,X1L,X1U,X2L,X2U,X1,X2,ANSY,ANSY1,
     *    ANSY2)
      DOUBLE PRECISION ANSY,ANSY1,ANSY2,X1,X1L,X1U,X2,X2L,X2U,Y(4),Y1(4)
     * ,Y12(4), Y2(4)
      INTEGER I
      DOUBLE PRECISION T,U,C(4,4)
      CALL BCUCOF(Y,Y1,Y2,Y12,X1U-X1L,X2U-X2L,C)
      IF(X1U.EQ.X1L.OR.X2U.EQ.X2L)THEN
        PRINT*,X1U,X1L,X2U,X2L,x1,x2
c        PAUSE 'BAD INPUT IN BCUINT'
      ENDIF
      T=(X1-X1L)/(X1U-X1L)
      U=(X2-X2L)/(X2U-X2L)
      ANSY=0.
      ANSY1=0.
      ANSY2=0.
      DO  I=4,1,-1
           ANSY=T*ANSY+((C(I,4)*U+C(I,3))*U+C(I,2))*U+C(I,1)
           ANSY2=T*ANSY2+(3.*C(I,4)*U+2.*C(I,3))*U+C(I,2)
           ANSY1=T*ANSY1+(3.*C(4,I)*T+2.*C(3,I))*T+C(2,I)
      ENDDO
      ANSY1=ANSY1/(X1U-X1L)
      ANSY2=ANSY2/(X2U-X2L)
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE BCUCOF(Y,Y1,Y2,Y12,D1,D2,C)
C----------------------------------------------------------------------C
      DOUBLE PRECISION D1,D2,C(4,4),Y(4),Y1(4),Y12(4),Y2(4)
      INTEGER I,J,K,L
      DOUBLE PRECISION D1D2,XX,CL(16),WT(16,16),X(16)
      SAVE WT
      DATA WT/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4
     *  ,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4
     *  ,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2
     *  ,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2
     *  ,0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2
     *  ,10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2
     *  ,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1
     *  ,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/
      D1D2=D1*D2

      DO  I=1,4
           X(I)=Y(I)
           X(I+4)=Y1(I)*D1
           X(I+8)=Y2(I)*D2
           X(I+12)=Y12(I)*D1D2
      ENDDO 
      DO  I=1,16
           XX = 0.
           DO  K=1,16
                XX=XX+WT(I,K)*X(K)
           ENDDO 
           CL(I)=XX
      ENDDO 
      L = 0 
      DO  I=1,4
           DO J=1,4
                L=L+1
                C(I,J)=CL(L)
           ENDDO 
      ENDDO 
      RETURN
      END
C-----------------------------------------------------------------------C
      SUBROUTINE CALCULATE_BZ(BZ,X,Y,Z)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z),INTEGER*4(I-N)
      PARAMETER(A = 0.796,B = 0.204,C=192.0,D=0.08491)

      P    = DSQRT(X*X + Y*Y)
      ARGZ = D*(ABS(Z) - C)

      CH   = COSH(ARGZ)*COSH(ARGZ)
      TH   = TANH(ARGZ)
      
      D0BZ = A - B*TH

      D2BZ = 2.*A*D**2*TH/CH

      D4BZ = 8.*A*D**4*(TH*TH*TH/(CH) - 2.*TH/(CH*CH))

      BZ   = D0BZ - 0.25*D2BZ*P**2 + D4BZ*P**4/32. 

      RETURN
      END

C----------------------------------------------------------------------C
      SUBROUTINE CALCULATE_BP(BP,X,Y,Z)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z),INTEGER*4(I-N)
      PARAMETER(A = 0.796,B = 0.204,C=192.0,D=0.08491)

      P    = DSQRT(X*X + Y*Y)
      ARGZ = D*(ABS(Z) - C)

      CH   = COSH(ARGZ)*COSH(ARGZ)
      TH   = TANH(ARGZ)*TANH(ARGZ)

      DBZ  = -B*D/CH

      D3BZ = 2*B*D*D*D/CH*(1./CH- 2.*TH)
     
      D5BZ = 8*B*D*D*D*D*D/CH*
     1        (6*TH/CH - (TH - 2./CH) * (2*TH + 1./CH))

      BP   = -DBZ*P/2. + D3BZ*PP**3/16. + D5BZ*P**5/C
    
      RETURN 
      END
