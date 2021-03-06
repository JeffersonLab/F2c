*DECK PDF
      SUBROUTINE PDF(X, QQ, UV, DV, US, DS, SS, GL)
Cf2py intent(in) x
Cf2py intent(in) qq
Cf2py intent(out) uv
Cf2py intent(out) dv
Cf2py intent(out) us
Cf2py intent(out) ds
Cf2py intent(out) ss
Cf2py intent(out) gl
C
C     GENERIC PARTON DENSITY
C     CALLS SPECIFIC PARAMETRIZTION
C
C     INPUT:
C     X     X VARIABLE FOR NUCLEON
C     QQ    Q2 (GEV**2)
C
C     OUTPUT:
C     UV   U QUARK, VALENCE
C     DV   D QUARK, VALENCE
C     US   U QUARK, SEA
C     DS   D QUARK, SEA
C     SS   S QUARK, SEA (SS = S = SBAR)
C     GL   GLUON
C
C     ALL DENSITIES ARE RETURNED AS  X*F(X)
C
C     HERE: NUCLEAR VARIATION
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
C
      COMMON /PDFI/ IPDF
C
C     ...GRV98 PDF PARAMETRIZATION, LO
C
      ISET = 1
      CALL GRV98PA(ISET, X, QQ, UV, DV, US, DS, SS, GL)
C
      IF (IPDF.NE.0) THEN
C
C     ...NULCEAR RATIO EPS09
C
         IORD = 1
         IA = 56
         Q = SQRT(QQ)
         CALL EPS09(IORD, IPDF, IA, X, Q,
     *          RUV, RDV, RU, RD, RS, RC, RB, RG)
C
C     ...HERE: GLUON ONLY
C
         GL = GL*RG
C
      ENDIF
C
      END
C
C---------------------------------------------------------------------
C
*DECK PDFX
      BLOCK DATA PDFX
C
C     EXTERNAL PARAMETERS FOR PDF
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
      CHARACTER*70 PDFC
C
      COMMON /PDFC/ PDFC
C
      DATA PDFC /'PDFS GRV98LO / EPS09LO'/
C
      END
C
C---------------------------------------------------------------------
C
*DECK GRV98
      BLOCK DATA GRV98
C
C     INITIALIZATION OF GRV98 PDF PARAMETRIZATION
C
C     NOTE: IF MULTIPLE PDF SETS ARE USED IN THE SAME RUN,
C     THE INITIALIZATION MUST BE RESET BY THE CALLING PROGRAM.
C
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
C
      COMMON /INTINIP/ IINIP
      COMMON /INTINIF/ IINIF
      DATA IINIP /0/
      DATA IINIF /0/
C
      END
C
C---------------------------------------------------------------------
C
