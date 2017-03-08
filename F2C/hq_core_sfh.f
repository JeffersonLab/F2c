      SUBROUTINE SETIPDF(IPDF1)
Cf2py intent(in) ipdf
      integer IPDF1
      COMMON /PDFI/ IPDF
      IPDF=IPDF1
      END

*DECK SFH
      SUBROUTINE SFH(FH, X1, QQ1, ISF1)
Cf2py intent(in) x1
Cf2py intent(in) qq1
Cf2py intent(in) isf1
Cf2py intent(out) fh
C
C     PACKAGE HQ
C     HEAVY-QUARK PRODUCTION IN LO APPROXIMATION
C     CORE ROUTINE
C
C     CALCULATE NUCLEON HEAVY-QUARK STRUCTURE FUNCTION 
C     IN LO APPROX
C
C     INPUT:
C     X     X VARIABLE FOR NUCLEON
C     QQ    Q2 (GEV**2)
C     ISF   SWITCH STRUCTURE FUNCTION
C           ISF = 1   F2H
C                 2   FLH
C
C     OUTPUT:
C     FH   NUCLEON HEAVY QUARK STRUCTURE FUNCTION
C
C     EXTERNAL PARAMETERS:
C     UH    HEAVY QUARK MASS (GEV)
C     EH    HEAVY QUARK CHARGE (IN UNITS OF E)
C
C     USES SLATEC LIBRARY ROUTINE DQAGS
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
      integer ISF1
C
C     ...EXTERNAL PARAMETERS SFH
C
      COMMON /SFHX/ UH, EH
C
      COMMON /SFHR/ X, QQ
      COMMON /SFHI/ ISF
C
C     WORKSPACE AND TOLERANCES FOR SLATEC ROUTINE DQAGS
C
      PARAMETER (LIMIT = 100, LENW  = 400)
      !PARAMETER (LIMIT = 100, LENW  = 4*LIMIT)
      DIMENSION IWORK(LIMIT), WORK(LENW)
      PARAMETER (EPSABS = 1.D-6, EPSREL = 1.D-6)
C
      EXTERNAL SFH1
C
      X  = X1
      QQ = QQ1
C
      ISF = ISF1
C
      A = 1 + 4*UH**2/QQ
      YMIN = A*X
C
C     ...CHECK WHETHER  X  WITHIN LIMITS
C
      IF (YMIN.LT.1.D0) THEN
C
C     ...COMPUTE INTEGRAL OVER  Y
C
         CALL DQAGS(SFH1, YMIN, 1.D0, EPSABS, EPSREL, RES, 
     *        ABSERR, NEVAL, IER, LIMIT, LENW, LAST, IWORK, WORK)
C
         FH = RES
C
      ELSE
C
         FH = 0.D0
C
      ENDIF
C
      END
C
C---------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION SFH1(Y)
C
C     NUCLEON HEAVY-QUARK STRUCTURE FUNCTION IN LO APPROX
C     INTEGRAND FOR SFH
C
C     ISF = 1   F2H
C           2   FLH
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
      PARAMETER (PI = 0.314159265358979D+01)
C
      COMMON /SFHX/ UH, EH
      COMMON /SFHR/ X, QQ
      COMMON /SFHI/ ISF
C
      XI   = QQ/UH**2
      Z    = X/Y
C
      BETA = SQRT(-4./XI*Z/(1 - Z) + 1.D0)
      ARG  = (1 + BETA)/(1 - BETA)
C
      IF      (ISF.EQ.1) THEN
C
         A0 =   -Z**2*(1 -   Z) + Z/2 
         A1 =  2*Z**2*(1 - 3*Z)
         A2 = -4*Z**3
C
         B0 =  4*Z**2*(1 - Z) - Z/2
         B1 = -2*Z**2*(1 - Z)
C
      ELSE IF (ISF.EQ.2) THEN
C
         A0 = 0
         A1 = -4*Z**3
         A2 = 0
C
         B0 =  2*Z**2*(1 - Z)
         B1 = 0
C
      ENDIF
C
      CG = ((A0 + A1/XI + A2/XI**2)*LOG(ARG)
     *   +  (B0 + B1/XI)*BETA) *4*PI/XI
C
C     ...QCD COUPLING  (SCALE 4*UH**2)
C
      UM2 = 4*UH**2
      CALL RC(ALPHA, UM2)
C
C     ...GLUON DENSITY (SCALE 4*UH**2)
C
      CALL PDF(Y, UM2, UVAL, DVAL, USEA, DSEA, SSEA, GLU)
C
      SFH1 = ALPHA/UH**2 *QQ/4./PI**2 *EH**2 *CG *GLU/Y
C
      END
C
C---------------------------------------------------------------------
C
      SUBROUTINE RC(ALPHA, UM2)
C
C     QCD RUNNING COUPLING ALPHA-S
C
C     INPUT:
C     UM2     SCALE MU**2
C
C     OUTPUT:
C     ALPHA   RUNNING COUPLING ALPHA-S
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
C
      PARAMETER (PI = 0.31415 92653 58979 D+01)
C
      COMMON /ALQCD/ ALQCD
      COMMON /NFL/   NFL
C
C     LAMBDA PARAMETERS FROM GRV98
C
C      PARAMETER (ALQCD = 0.2041D0, NFL = 3)
C      PARAMETER (ALQCD = 0.1750D0, NFL = 4
C
C      PARAMETER (ALQCD = 0.192D0, NFL = 4)
C      PARAMETER (ALQCD = 0.232D0, NFL = 3)
C
      BETA0 = 11. - (2.D0/3.D0)*NFL
C
      ALPHA = 4*PI/BETA0/LOG(UM2/ALQCD**2)
C
      END
C
C---------------------------------------------------------------------
C
*DECK RCX1
      BLOCK DATA RCX1
C
C     SET GLOBAL PARAMETERS FOR RC
C
C     ALQCD  LAMBDA-QDC
C     NFL    NUMBER OF ACTIVE FLAVORS
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
C
      COMMON /ALQCD/ ALQCD
      COMMON /NFL/   NFL
C
      DATA ALQCD /0.2041D0/
      DATA NFL   /3/
C
      END
C
C---------------------------------------------------------------------
C
*DECK SFHX1
      BLOCK DATA SFHX1
C
C     SET GLOBAL PARAMETERS FOR SF
C
C     UH  HEAVY QUARK MASS (GEV)
C     EF  HEAVY QUARK CHARGE (IN UNITS OF E)
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
C
      COMMON /SFHX/ UH, EH
C
      DATA UH    /1.5D0/
      DATA EH    /.6666666D0/
C
      END
C
C---------------------------------------------------------------------
C
