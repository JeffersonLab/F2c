*DECK DQAGS
      SUBROUTINE DQAGS(F,A,B,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,IER,
     1   LIMIT,LENW,LAST,IWORK,WORK)
C***BEGIN PROLOGUE  DQAGS
C***DATE WRITTEN   800101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  H2A1A1
C***KEYWORDS  (END-POINT) SINGULARITIES,AUTOMATIC INTEGRATOR,
C             EXTRAPOLATION,GENERAL-PURPOSE,GLOBALLY ADAPTIVE
C***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C           DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C***PURPOSE  The routine calculates an approximation result to a given
C            Definite integral  I = Integral of F over (A,B),
C            Hopefully satisfying following claim for accuracy
C            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
C***DESCRIPTION
C
C        Computation of a definite integral
C        Standard fortran subroutine
C        Double precision version
C
C
C        PARAMETERS
C         ON ENTRY
C            F      - Double precision
C                     Function subprogram defining the integrand
C                     Function F(X). The actual name for F needs to be
C                     Declared E X T E R N A L in the driver program.
C
C            A      - Double precision
C                     Lower limit of integration
C
C            B      - Double precision
C                     Upper limit of integration
C
C            EPSABS - Double precision
C                     Absolute accuracy requested
C            EPSREL - Double precision
C                     Relative accuracy requested
C                     If  EPSABS.LE.0
C                     And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                     The routine will end with IER = 6.
C
C         ON RETURN
C            RESULT - Double precision
C                     Approximation to the integral
C
C            ABSERR - Double precision
C                     Estimate of the modulus of the absolute error,
C                     which should equal or exceed ABS(I-RESULT)
C
C            NEVAL  - Integer
C                     Number of integrand evaluations
C
C            IER    - Integer
C                     IER = 0 Normal and reliable termination of the
C                             routine. It is assumed that the requested
C                             accuracy has been achieved.
C                     IER.GT.0 Abnormal termination of the routine
C                             The estimates for integral and error are
C                             less reliable. It is assumed that the
C                             requested accuracy has not been achieved.
C            ERROR MESSAGES
C                     IER = 1 Maximum number of subdivisions allowed
C                             has been achieved. One can allow more sub-
C                             divisions by increasing the value of LIMIT
C                             (and taking the according dimension
C                             adjustments into account. However, if
C                             this yields no improvement it is advised
C                             to analyze the integrand in order to
C                             determine the integration difficulties. If
C                             the position of a local difficulty can be
C                             determined (E.G. SINGULARITY,
C                             DISCONTINUITY WITHIN THE INTERVAL) one
C                             will probably gain from splitting up the
C                             interval at this point and calling the
C                             integrator on the subranges. If possible,
C                             an appropriate special-purpose integrator
C                             should be used, which is designed for
C                             handling the type of difficulty involved.
C                         = 2 The occurrence of roundoff error is detec-
C                             ted, which prevents the requested
C                             tolerance from being achieved.
C                             The error may be under-estimated.
C                         = 3 Extremely bad integrand behaviour
C                             occurs at some points of the integration
C                             interval.
C                         = 4 The algorithm does not converge.
C                             Roundoff error is detected in the
C                             Extrapolation table. It is presumed that
C                             the requested tolerance cannot be
C                             achieved, and that the returned result is
C                             the best which can be obtained.
C                         = 5 The integral is probably divergent, or
C                             slowly convergent. It must be noted that
C                             divergence can occur with any other value
C                             of IER.
C                         = 6 The input is invalid, because
C                             (EPSABS.LE.0 AND
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28)
C                             OR LIMIT.LT.1 OR LENW.LT.LIMIT*4.
C                             RESULT, ABSERR, NEVAL, LAST are set to
C                             zero.Except when LIMIT or LENW is invalid,
C                             IWORK(1), WORK(LIMIT*2+1) and
C                             WORK(LIMIT*3+1) are set to zero, WORK(1)
C                             is set to A and WORK(LIMIT+1) TO B.
C
C         DIMENSIONING PARAMETERS
C            LIMIT - Integer
C                    DIMENSIONING PARAMETER FOR IWORK
C                    LIMIT determines the maximum number of subintervals
C                    in the partition of the given integration interval
C                    (A,B), LIMIT.GE.1.
C                    IF LIMIT.LT.1, the routine will end with IER = 6.
C
C            LENW  - Integer
C                    DIMENSIONING PARAMETER FOR WORK
C                    LENW must be at least LIMIT*4.
C                    If LENW.LT.LIMIT*4, the routine will end
C                    with IER = 6.
C
C            LAST  - Integer
C                    On return, LAST equals the number of subintervals
C                    produced in the subdivision process, detemines the
C                    number of significant elements actually in the WORK
C                    Arrays.
C
C         WORK ARRAYS
C            IWORK - Integer
C                    Vector of dimension at least LIMIT, the first K
C                    elements of which contain pointers
C                    to the error estimates over the subintervals
C                    such that WORK(LIMIT*3+IWORK(1)),... ,
C                    WORK(LIMIT*3+IWORK(K)) form a decreasing
C                    sequence, with K = LAST IF LAST.LE.(LIMIT/2+2),
C                    and K = LIMIT+1-LAST otherwise
C
C            WORK  - Double precision
C                    Vector of dimension at least LENW
C                    on return
C                    WORK(1), ..., WORK(LAST) contain the left
C                     end-points of the subintervals in the
C                     partition of (A,B),
C                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain
C                     the right end-points,
C                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
C                     the integral approximations over the subintervals,
C                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST)
C                     contain the error estimates.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DQAGSE,XERROR
C***END PROLOGUE  DQAGS
C
C
      DOUBLE PRECISION A,ABSERR,B,EPSABS,EPSREL,F,RESULT,WORK
      INTEGER IER,IWORK,LAST,LENW,LIMIT,LVL,L1,L2,L3,NEVAL
C
      DIMENSION IWORK(LIMIT),WORK(LENW)
C
!      EXTERNAL F
C
C         CHECK VALIDITY OF LIMIT AND LENW.
C
C***FIRST EXECUTABLE STATEMENT  DQAGS
      IER = 6
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      IF(LIMIT.LT.1.OR.LENW.LT.LIMIT*4) GO TO 10
C
C         PREPARE CALL FOR DQAGSE.
C
      L1 = LIMIT+1
      L2 = LIMIT+L1
      L3 = LIMIT+L2
C
      CALL DQAGSE(F,A,B,EPSABS,EPSREL,LIMIT,RESULT,ABSERR,NEVAL,
     1  IER,WORK(1),WORK(L1),WORK(L2),WORK(L3),IWORK,LAST)
C
C         CALL ERROR HANDLER IF NECESSARY.
C
      LVL = 0
10    IF(IER.EQ.6) LVL = 1
      IF(IER.NE.0) CALL XERROR( 'ABNORMAL RETURN FROM DQAGS',26,IER,LVL)
      RETURN
      END
*DECK DQAGSE
      SUBROUTINE DQAGSE(F,A,B,EPSABS,EPSREL,LIMIT,RESULT,ABSERR,NEVAL,
     1   IER,ALIST,BLIST,RLIST,ELIST,IORD,LAST)
C***BEGIN PROLOGUE  DQAGSE
C***DATE WRITTEN   800101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  H2A1A1
C***KEYWORDS  (END POINT) SINGULARITIES,AUTOMATIC INTEGRATOR,
C             EXTRAPOLATION,GENERAL-PURPOSE,GLOBALLY ADAPTIVE
C***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C           DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C***PURPOSE  The routine calculates an approximation result to a given
C            definite integral I = Integral of F over (A,B),
C            hopefully satisfying following claim for accuracy
C            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
C***DESCRIPTION
C
C        Computation of a definite integral
C        Standard fortran subroutine
C        Double precision version
C
C        PARAMETERS
C         ON ENTRY
C            F      - Double precision
C                     Function subprogram defining the integrand
C                     function F(X). The actual name for F needs to be
C                     declared E X T E R N A L in the driver program.
C
C            A      - Double precision
C                     Lower limit of integration
C
C            B      - Double precision
C                     Upper limit of integration
C
C            EPSABS - Double precision
C                     Absolute accuracy requested
C            EPSREL - Double precision
C                     Relative accuracy requested
C                     If  EPSABS.LE.0
C                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                     the routine will end with IER = 6.
C
C            LIMIT  - Integer
C                     Gives an upperbound on the number of subintervals
C                     in the partition of (A,B)
C
C         ON RETURN
C            RESULT - Double precision
C                     Approximation to the integral
C
C            ABSERR - Double precision
C                     Estimate of the modulus of the absolute error,
C                     which should equal or exceed ABS(I-RESULT)
C
C            NEVAL  - Integer
C                     Number of integrand evaluations
C
C            IER    - Integer
C                     IER = 0 Normal and reliable termination of the
C                             routine. It is assumed that the requested
C                             accuracy has been achieved.
C                     IER.GT.0 Abnormal termination of the routine
C                             the estimates for integral and error are
C                             less reliable. It is assumed that the
C                             requested accuracy has not been achieved.
C            ERROR MESSAGES
C                         = 1 Maximum number of subdivisions allowed
C                             has been achieved. One can allow more sub-
C                             divisions by increasing the value of LIMIT
C                             (and taking the according dimension
C                             adjustments into account). However, if
C                             this yields no improvement it is advised
C                             to analyze the integrand in order to
C                             determine the integration difficulties. If
C                             the position of a local difficulty can be
C                             determined (e.g. singularity,
C                             discontinuity within the interval) one
C                             will probably gain from splitting up the
C                             interval at this point and calling the
C                             integrator on the subranges. If possible,
C                             an appropriate special-purpose integrator
C                             should be used, which is designed for
C                             handling the type of difficulty involved.
C                         = 2 The occurrence of roundoff error is detec-
C                             ted, which prevents the requested
C                             tolerance from being achieved.
C                             The error may be under-estimated.
C                         = 3 Extremely bad integrand behaviour
C                             occurs at some points of the integration
C                             interval.
C                         = 4 The algorithm does not converge.
C                             Roundoff error is detected in the
C                             extrapolation table.
C                             It is presumed that the requested
C                             tolerance cannot be achieved, and that the
C                             returned result is the best which can be
C                             obtained.
C                         = 5 The integral is probably divergent, or
C                             slowly convergent. It must be noted that
C                             divergence can occur with any other value
C                             of IER.
C                         = 6 The input is invalid, because
C                             EPSABS.LE.0 and
C                             EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28).
C                             RESULT, ABSERR, NEVAL, LAST, RLIST(1),
C                             IORD(1) and ELIST(1) are set to zero.
C                             ALIST(1) and BLIST(1) are set to A and B
C                             respectively.
C
C            ALIST  - Double precision
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the left end points
C                     of the subintervals in the partition of the
C                     given integration range (A,B)
C
C            BLIST  - Double precision
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the right end points
C                     of the subintervals in the partition of the given
C                     integration range (A,B)
C
C            RLIST  - Double precision
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the integral
C                     approximations on the subintervals
C
C            ELIST  - Double precision
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the moduli of the
C                     absolute error estimates on the subintervals
C
C            IORD   - Integer
C                     Vector of dimension at least LIMIT, the first K
C                     elements of which are pointers to the
C                     error estimates over the subintervals,
C                     such that ELIST(IORD(1)), ..., ELIST(IORD(K))
C                     form a decreasing sequence, with K = LAST
C                     If LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST
C                     otherwise
C
C            LAST   - Integer
C                     Number of subintervals actually produced in the
C                     subdivision process
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,DQELG,DQK21,DQPSRT
C***END PROLOGUE  DQAGSE
C
      DOUBLE PRECISION A,ABSEPS,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,
     1  A2,B,BLIST,B1,B2,CORREC,DABS,DEFABS,DEFAB1,DEFAB2,D1MACH,DMAX1,
     2  DRES,ELIST,EPMACH,EPSABS,EPSREL,ERLARG,ERLAST,ERRBND,ERRMAX,
     3  ERROR1,ERROR2,ERRO12,ERRSUM,ERTEST,F,OFLOW,RESABS,RESEPS,RESULT,
     4  RES3LA,RLIST,RLIST2,SMALL,UFLOW
      INTEGER ID,IER,IERRO,IORD,IROFF1,IROFF2,IROFF3,JUPBND,K,KSGN,
     1  KTMIN,LAST,LIMIT,MAXERR,NEVAL,NRES,NRMAX,NUMRL2
      LOGICAL EXTRAP,NOEXT
C
      DIMENSION ALIST(LIMIT),BLIST(LIMIT),ELIST(LIMIT),IORD(LIMIT),
     1 RES3LA(3),RLIST(LIMIT),RLIST2(52)
C
!      EXTERNAL F
C
C            THE DIMENSION OF RLIST2 IS DETERMINED BY THE VALUE OF
C            LIMEXP IN SUBROUTINE DQELG (RLIST2 SHOULD BE OF DIMENSION
C            (LIMEXP+2) AT LEAST).
C
C            LIST OF MAJOR VARIABLES
C            -----------------------
C
C           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
C                       (ALIST(I),BLIST(I))
C           RLIST2    - ARRAY OF DIMENSION AT LEAST LIMEXP+2 CONTAINING
C                       THE PART OF THE EPSILON TABLE WHICH IS STILL
C                       NEEDED FOR FURTHER COMPUTATIONS
C           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
C           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST ERROR
C                       ESTIMATE
C           ERRMAX    - ELIST(MAXERR)
C           ERLAST    - ERROR ON THE INTERVAL CURRENTLY SUBDIVIDED
C                       (BEFORE THAT SUBDIVISION HAS TAKEN PLACE)
C           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
C           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
C           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
C                       ABS(RESULT))
C           *****1    - VARIABLE FOR THE LEFT INTERVAL
C           *****2    - VARIABLE FOR THE RIGHT INTERVAL
C           LAST      - INDEX FOR SUBDIVISION
C           NRES      - NUMBER OF CALLS TO THE EXTRAPOLATION ROUTINE
C           NUMRL2    - NUMBER OF ELEMENTS CURRENTLY IN RLIST2. IF AN
C                       APPROPRIATE APPROXIMATION TO THE COMPOUNDED
C                       INTEGRAL HAS BEEN OBTAINED IT IS PUT IN
C                       RLIST2(NUMRL2) AFTER NUMRL2 HAS BEEN INCREASED
C                       BY ONE.
C           SMALL     - LENGTH OF THE SMALLEST INTERVAL CONSIDERED UP
C                       TO NOW, MULTIPLIED BY 1.5
C           ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER
C                       THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW
C           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE IS
C                       ATTEMPTING TO PERFORM EXTRAPOLATION I.E. BEFORE
C                       SUBDIVIDING THE SMALLEST INTERVAL WE TRY TO
C                       DECREASE THE VALUE OF ERLARG.
C           NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION
C                       IS NO LONGER ALLOWED (TRUE VALUE)
C
C            MACHINE DEPENDENT CONSTANTS
C            ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  DQAGSE
      EPMACH = D1MACH(4)
C
C            TEST ON VALIDITY OF PARAMETERS
C            ------------------------------
      IER = 0
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      ALIST(1) = A
      BLIST(1) = B
      RLIST(1) = 0.0D+00
      ELIST(1) = 0.0D+00
      IF(EPSABS.LE.0.0D+00.AND.EPSREL.LT.DMAX1(0.5D+02*EPMACH,0.5D-28))
     1   IER = 6
      IF(IER.EQ.6) GO TO 999
C
C           FIRST APPROXIMATION TO THE INTEGRAL
C           -----------------------------------
C
      UFLOW = D1MACH(1)
      OFLOW = D1MACH(2)
      IERRO = 0
      CALL DQK21(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
C
C           TEST ON ACCURACY.
C
      DRES = DABS(RESULT)
      ERRBND = DMAX1(EPSABS,EPSREL*DRES)
      LAST = 1
      RLIST(1) = RESULT
      ELIST(1) = ABSERR
      IORD(1) = 1
      IF(ABSERR.LE.1.0D+02*EPMACH*DEFABS.AND.ABSERR.GT.ERRBND) IER = 2
      IF(LIMIT.EQ.1) IER = 1
      IF(IER.NE.0.OR.(ABSERR.LE.ERRBND.AND.ABSERR.NE.RESABS).OR.
     1  ABSERR.EQ.0.0D+00) GO TO 140
C
C           INITIALIZATION
C           --------------
C
      RLIST2(1) = RESULT
      ERRMAX = ABSERR
      MAXERR = 1
      AREA = RESULT
      ERRSUM = ABSERR
      ABSERR = OFLOW
      NRMAX = 1
      NRES = 0
      NUMRL2 = 2
      KTMIN = 0
      EXTRAP = .FALSE.
      NOEXT = .FALSE.
      IROFF1 = 0
      IROFF2 = 0
      IROFF3 = 0
      KSGN = -1
      IF(DRES.GE.(0.1D+01-0.5D+02*EPMACH)*DEFABS) KSGN = 1
C
C           MAIN DO-LOOP
C           ------------
C
      DO 90 LAST = 2,LIMIT
C
C           BISECT THE SUBINTERVAL WITH THE NRMAX-TH LARGEST ERROR
C           ESTIMATE.
C
        A1 = ALIST(MAXERR)
        B1 = 0.5D+00*(ALIST(MAXERR)+BLIST(MAXERR))
        A2 = B1
        B2 = BLIST(MAXERR)
        ERLAST = ERRMAX
        CALL DQK21(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        CALL DQK21(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
C
C           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
C           AND ERROR AND TEST FOR ACCURACY.
C
        AREA12 = AREA1+AREA2
        ERRO12 = ERROR1+ERROR2
        ERRSUM = ERRSUM+ERRO12-ERRMAX
        AREA = AREA+AREA12-RLIST(MAXERR)
        IF(DEFAB1.EQ.ERROR1.OR.DEFAB2.EQ.ERROR2) GO TO 15
        IF(DABS(RLIST(MAXERR)-AREA12).GT.0.1D-04*DABS(AREA12)
     1  .OR.ERRO12.LT.0.99D+00*ERRMAX) GO TO 10
        IF(EXTRAP) IROFF2 = IROFF2+1
        IF(.NOT.EXTRAP) IROFF1 = IROFF1+1
   10   IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF3 = IROFF3+1
   15   RLIST(MAXERR) = AREA1
        RLIST(LAST) = AREA2
        ERRBND = DMAX1(EPSABS,EPSREL*DABS(AREA))
C
C           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
C
        IF(IROFF1+IROFF2.GE.10.OR.IROFF3.GE.20) IER = 2
        IF(IROFF2.GE.5) IERRO = 3
C
C           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF SUBINTERVALS
C           EQUALS LIMIT.
C
        IF(LAST.EQ.LIMIT) IER = 1
C
C           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
C           AT A POINT OF THE INTEGRATION RANGE.
C
        IF(DMAX1(DABS(A1),DABS(B2)).LE.(0.1D+01+0.1D+03*EPMACH)*
     1  (DABS(A2)+0.1D+04*UFLOW)) IER = 4
C
C           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
C
        IF(ERROR2.GT.ERROR1) GO TO 20
        ALIST(LAST) = A2
        BLIST(MAXERR) = B1
        BLIST(LAST) = B2
        ELIST(MAXERR) = ERROR1
        ELIST(LAST) = ERROR2
        GO TO 30
   20   ALIST(MAXERR) = A2
        ALIST(LAST) = A1
        BLIST(LAST) = B1
        RLIST(MAXERR) = AREA2
        RLIST(LAST) = AREA1
        ELIST(MAXERR) = ERROR2
        ELIST(LAST) = ERROR1
C
C           CALL SUBROUTINE DQPSRT TO MAINTAIN THE DESCENDING ORDERING
C           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
C           WITH NRMAX-TH LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
C
   30   CALL DQPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
C ***JUMP OUT OF DO-LOOP
        IF(ERRSUM.LE.ERRBND) GO TO 115
C ***JUMP OUT OF DO-LOOP
        IF(IER.NE.0) GO TO 100
        IF(LAST.EQ.2) GO TO 80
        IF(NOEXT) GO TO 90
        ERLARG = ERLARG-ERLAST
        IF(DABS(B1-A1).GT.SMALL) ERLARG = ERLARG+ERRO12
        IF(EXTRAP) GO TO 40
C
C           TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE
C           SMALLEST INTERVAL.
C
        IF(DABS(BLIST(MAXERR)-ALIST(MAXERR)).GT.SMALL) GO TO 90
        EXTRAP = .TRUE.
        NRMAX = 2
   40   IF(IERRO.EQ.3.OR.ERLARG.LE.ERTEST) GO TO 60
C
C           THE SMALLEST INTERVAL HAS THE LARGEST ERROR.
C           BEFORE BISECTING DECREASE THE SUM OF THE ERRORS OVER THE
C           LARGER INTERVALS (ERLARG) AND PERFORM EXTRAPOLATION.
C
        ID = NRMAX
        JUPBND = LAST
        IF(LAST.GT.(2+LIMIT/2)) JUPBND = LIMIT+3-LAST
        DO 50 K = ID,JUPBND
          MAXERR = IORD(NRMAX)
          ERRMAX = ELIST(MAXERR)
C ***JUMP OUT OF DO-LOOP
          IF(DABS(BLIST(MAXERR)-ALIST(MAXERR)).GT.SMALL) GO TO 90
          NRMAX = NRMAX+1
   50   CONTINUE
C
C           PERFORM EXTRAPOLATION.
C
   60   NUMRL2 = NUMRL2+1
        RLIST2(NUMRL2) = AREA
        CALL DQELG(NUMRL2,RLIST2,RESEPS,ABSEPS,RES3LA,NRES)
        KTMIN = KTMIN+1
        IF(KTMIN.GT.5.AND.ABSERR.LT.0.1D-02*ERRSUM) IER = 5
        IF(ABSEPS.GE.ABSERR) GO TO 70
        KTMIN = 0
        ABSERR = ABSEPS
        RESULT = RESEPS
        CORREC = ERLARG
        ERTEST = DMAX1(EPSABS,EPSREL*DABS(RESEPS))
C ***JUMP OUT OF DO-LOOP
        IF(ABSERR.LE.ERTEST) GO TO 100
C
C           PREPARE BISECTION OF THE SMALLEST INTERVAL.
C
   70   IF(NUMRL2.EQ.1) NOEXT = .TRUE.
        IF(IER.EQ.5) GO TO 100
        MAXERR = IORD(1)
        ERRMAX = ELIST(MAXERR)
        NRMAX = 1
        EXTRAP = .FALSE.
        SMALL = SMALL*0.5D+00
        ERLARG = ERRSUM
        GO TO 90
   80   SMALL = DABS(B-A)*0.375D+00
        ERLARG = ERRSUM
        ERTEST = ERRBND
        RLIST2(2) = AREA
   90 CONTINUE
C
C           SET FINAL RESULT AND ERROR ESTIMATE.
C           ------------------------------------
C
  100 IF(ABSERR.EQ.OFLOW) GO TO 115
      IF(IER+IERRO.EQ.0) GO TO 110
      IF(IERRO.EQ.3) ABSERR = ABSERR+CORREC
      IF(IER.EQ.0) IER = 3
      IF(RESULT.NE.0.0D+00.AND.AREA.NE.0.0D+00) GO TO 105
      IF(ABSERR.GT.ERRSUM) GO TO 115
      IF(AREA.EQ.0.0D+00) GO TO 130
      GO TO 110
  105 IF(ABSERR/DABS(RESULT).GT.ERRSUM/DABS(AREA)) GO TO 115
C
C           TEST ON DIVERGENCE.
C
  110 IF(KSGN.EQ.(-1).AND.DMAX1(DABS(RESULT),DABS(AREA)).LE.
     1 DEFABS*0.1D-01) GO TO 130
      IF(0.1D-01.GT.(RESULT/AREA).OR.(RESULT/AREA).GT.0.1D+03
     1 .OR.ERRSUM.GT.DABS(AREA)) IER = 6
      GO TO 130
C
C           COMPUTE GLOBAL INTEGRAL SUM.
C
  115 RESULT = 0.0D+00
      DO 120 K = 1,LAST
         RESULT = RESULT+RLIST(K)
  120 CONTINUE
      ABSERR = ERRSUM
  130 IF(IER.GT.2) IER = IER-1
  140 NEVAL = 42*LAST-21
  999 RETURN
      END
*DECK DQK21
      SUBROUTINE DQK21(F,A,B,RESULT,ABSERR,RESABS,RESASC)
C***BEGIN PROLOGUE  DQK21
C***DATE WRITTEN   800101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  H2A1A2
C***KEYWORDS  21-POINT GAUSS-KRONROD RULES
C***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C           DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C***PURPOSE  To compute I = Integral of F over (A,B), with error
C                           estimate
C                       J = Integral of ABS(F) over (A,B)
C***DESCRIPTION
C
C           Integration rules
C           Standard fortran subroutine
C           Double precision version
C
C           PARAMETERS
C            ON ENTRY
C              F      - Double precision
C                       Function subprogram defining the integrand
C                       FUNCTION F(X). The actual name for F needs to be
C                       Declared E X T E R N A L in the driver program.
C
C              A      - Double precision
C                       Lower limit of integration
C
C              B      - Double precision
C                       Upper limit of integration
C
C            ON RETURN
C              RESULT - Double precision
C                       Approximation to the integral I
C                       RESULT is computed by applying the 21-POINT
C                       KRONROD RULE (RESK) obtained by optimal addition
C                       of abscissae to the 10-POINT GAUSS RULE (RESG).
C
C              ABSERR - Double precision
C                       Estimate of the modulus of the absolute error,
C                       which should not exceed ABS(I-RESULT)
C
C              RESABS - Double precision
C                       Approximation to the integral J
C
C              RESASC - Double precision
C                       Approximation to the integral of ABS(F-I/(B-A))
C                       over (A,B)
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH
C***END PROLOGUE  DQK21
C
      DOUBLE PRECISION A,ABSC,ABSERR,B,CENTR,DABS,DHLGTH,DMAX1,DMIN1,
     1  D1MACH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,RESABS,RESASC,
     2  RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK
      INTEGER J,JTW,JTWM1
!      EXTERNAL F
C
      DIMENSION FV1(10),FV2(10),WG(5),WGK(11),XGK(11)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 21-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 10-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 10-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 21-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 10-POINT GAUSS RULE
C
C
C GAUSS QUADRATURE WEIGHTS AND KRONRON QUADRATURE ABSCISSAE AND WEIGHTS
C AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
C BELL LABS, NOV. 1981.
C
      DATA WG  (  1) / 0.0666713443 0868813759 3568809893 332 D0 /
      DATA WG  (  2) / 0.1494513491 5058059314 5776339657 697 D0 /
      DATA WG  (  3) / 0.2190863625 1598204399 5534934228 163 D0 /
      DATA WG  (  4) / 0.2692667193 0999635509 1226921569 469 D0 /
      DATA WG  (  5) / 0.2955242247 1475287017 3892994651 338 D0 /
C
      DATA XGK (  1) / 0.9956571630 2580808073 5527280689 003 D0 /
      DATA XGK (  2) / 0.9739065285 1717172007 7964012084 452 D0 /
      DATA XGK (  3) / 0.9301574913 5570822600 1207180059 508 D0 /
      DATA XGK (  4) / 0.8650633666 8898451073 2096688423 493 D0 /
      DATA XGK (  5) / 0.7808177265 8641689706 3717578345 042 D0 /
      DATA XGK (  6) / 0.6794095682 9902440623 4327365114 874 D0 /
      DATA XGK (  7) / 0.5627571346 6860468333 9000099272 694 D0 /
      DATA XGK (  8) / 0.4333953941 2924719079 9265943165 784 D0 /
      DATA XGK (  9) / 0.2943928627 0146019813 1126603103 866 D0 /
      DATA XGK ( 10) / 0.1488743389 8163121088 4826001129 720 D0 /
      DATA XGK ( 11) / 0.0000000000 0000000000 0000000000 000 D0 /
C
      DATA WGK (  1) / 0.0116946388 6737187427 8064396062 192 D0 /
      DATA WGK (  2) / 0.0325581623 0796472747 8818972459 390 D0 /
      DATA WGK (  3) / 0.0547558965 7435199603 1381300244 580 D0 /
      DATA WGK (  4) / 0.0750396748 1091995276 7043140916 190 D0 /
      DATA WGK (  5) / 0.0931254545 8369760553 5065465083 366 D0 /
      DATA WGK (  6) / 0.1093871588 0229764189 9210590325 805 D0 /
      DATA WGK (  7) / 0.1234919762 6206585107 7958109831 074 D0 /
      DATA WGK (  8) / 0.1347092173 1147332592 8054001771 707 D0 /
      DATA WGK (  9) / 0.1427759385 7706008079 7094273138 717 D0 /
      DATA WGK ( 10) / 0.1477391049 0133849137 4841515972 068 D0 /
      DATA WGK ( 11) / 0.1494455540 0291690566 4936468389 821 D0 /
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 10-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 21-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  DQK21
      EPMACH = D1MACH(4)
      UFLOW = D1MACH(1)
C
      CENTR = 0.5D+00*(A+B)
      HLGTH = 0.5D+00*(B-A)
      DHLGTH = DABS(HLGTH)
C
C           COMPUTE THE 21-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      RESG = 0.0D+00
      FC = F(CENTR)
      RESK = WGK(11)*FC
      RESABS = DABS(RESK)
      DO 10 J=1,5
        JTW = 2*J
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(DABS(FVAL1)+DABS(FVAL2))
   10 CONTINUE
      DO 15 J = 1,5
        JTWM1 = 2*J-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(DABS(FVAL1)+DABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*0.5D+00
      RESASC = WGK(11)*DABS(FC-RESKH)
      DO 20 J=1,10
        RESASC = RESASC+WGK(J)*(DABS(FV1(J)-RESKH)+DABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = DABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00)
     1  ABSERR = RESASC*DMIN1(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
      IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = DMAX1
     1  ((EPMACH*0.5D+02)*RESABS,ABSERR)
      RETURN
      END
*DECK DQELG
      SUBROUTINE DQELG(N,EPSTAB,RESULT,ABSERR,RES3LA,NRES)
C***BEGIN PROLOGUE  DQELG
C***REFER TO  DQAGIE,DQAGOE,DQAGPE,DQAGSE
C***ROUTINES CALLED  D1MACH
C***REVISION DATE  830518   (YYMMDD)
C***KEYWORDS  CONVERGENCE ACCELERATION,EPSILON ALGORITHM,EXTRAPOLATION
C***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C           DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C***PURPOSE  The routine determines the limit of a given sequence of
C            approximations, by means of the Epsilon algorithm of
C            P.Wynn. An estimate of the absolute error is also given.
C            The condensed Epsilon table is computed. Only those
C            elements needed for the computation of the next diagonal
C            are preserved.
C***DESCRIPTION
C
C           Epsilon algorithm
C           Standard fortran subroutine
C           Double precision version
C
C           PARAMETERS
C              N      - Integer
C                       EPSTAB(N) contains the new element in the
C                       first column of the epsilon table.
C
C              EPSTAB - Double precision
C                       Vector of dimension 52 containing the elements
C                       of the two lower diagonals of the triangular
C                       epsilon table. The elements are numbered
C                       starting at the right-hand corner of the
C                       triangle.
C
C              RESULT - Double precision
C                       Resulting approximation to the integral
C
C              ABSERR - Double precision
C                       Estimate of the absolute error computed from
C                       RESULT and the 3 previous results
C
C              RES3LA - Double precision
C                       Vector of dimension 3 containing the last 3
C                       results
C
C              NRES   - Integer
C                       Number of calls to the routine
C                       (should be zero at first call)
C***END PROLOGUE  DQELG
C
      DOUBLE PRECISION ABSERR,DABS,DELTA1,DELTA2,DELTA3,DMAX1,D1MACH,
     1  EPMACH,EPSINF,EPSTAB,ERROR,ERR1,ERR2,ERR3,E0,E1,E1ABS,E2,E3,
     2  OFLOW,RES,RESULT,RES3LA,SS,TOL1,TOL2,TOL3
      INTEGER I,IB,IB2,IE,INDX,K1,K2,K3,LIMEXP,N,NEWELM,NRES,NUM
      DIMENSION EPSTAB(52),RES3LA(3)
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           E0     - THE 4 ELEMENTS ON WHICH THE COMPUTATION OF A NEW
C           E1       ELEMENT IN THE EPSILON TABLE IS BASED
C           E2
C           E3                 E0
C                        E3    E1    NEW
C                              E2
C           NEWELM - NUMBER OF ELEMENTS TO BE COMPUTED IN THE NEW
C                    DIAGONAL
C           ERROR  - ERROR = ABS(E1-E0)+ABS(E2-E1)+ABS(NEW-E2)
C           RESULT - THE ELEMENT IN THE NEW DIAGONAL WITH LEAST VALUE
C                    OF ERROR
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
C           LIMEXP IS THE MAXIMUM NUMBER OF ELEMENTS THE EPSILON
C           TABLE CAN CONTAIN. IF THIS NUMBER IS REACHED, THE UPPER
C           DIAGONAL OF THE EPSILON TABLE IS DELETED.
C
C***FIRST EXECUTABLE STATEMENT  DQELG
      EPMACH = D1MACH(4)
      OFLOW = D1MACH(2)
      NRES = NRES+1
      ABSERR = OFLOW
      RESULT = EPSTAB(N)
      IF(N.LT.3) GO TO 100
      LIMEXP = 50
      EPSTAB(N+2) = EPSTAB(N)
      NEWELM = (N-1)/2
      EPSTAB(N) = OFLOW
      NUM = N
      K1 = N
      DO 40 I = 1,NEWELM
        K2 = K1-1
        K3 = K1-2
        RES = EPSTAB(K1+2)
        E0 = EPSTAB(K3)
        E1 = EPSTAB(K2)
        E2 = RES
        E1ABS = DABS(E1)
        DELTA2 = E2-E1
        ERR2 = DABS(DELTA2)
        TOL2 = DMAX1(DABS(E2),E1ABS)*EPMACH
        DELTA3 = E1-E0
        ERR3 = DABS(DELTA3)
        TOL3 = DMAX1(E1ABS,DABS(E0))*EPMACH
        IF(ERR2.GT.TOL2.OR.ERR3.GT.TOL3) GO TO 10
C
C           IF E0, E1 AND E2 ARE EQUAL TO WITHIN MACHINE
C           ACCURACY, CONVERGENCE IS ASSUMED.
C           RESULT = E2
C           ABSERR = ABS(E1-E0)+ABS(E2-E1)
C
        RESULT = RES
        ABSERR = ERR2+ERR3
C ***JUMP OUT OF DO-LOOP
        GO TO 100
   10   E3 = EPSTAB(K1)
        EPSTAB(K1) = E1
        DELTA1 = E1-E3
        ERR1 = DABS(DELTA1)
        TOL1 = DMAX1(E1ABS,DABS(E3))*EPMACH
C
C           IF TWO ELEMENTS ARE VERY CLOSE TO EACH OTHER, OMIT
C           A PART OF THE TABLE BY ADJUSTING THE VALUE OF N
C
        IF(ERR1.LE.TOL1.OR.ERR2.LE.TOL2.OR.ERR3.LE.TOL3) GO TO 20
        SS = 0.1D+01/DELTA1+0.1D+01/DELTA2-0.1D+01/DELTA3
        EPSINF = DABS(SS*E1)
C
C           TEST TO DETECT IRREGULAR BEHAVIOUR IN THE TABLE, AND
C           EVENTUALLY OMIT A PART OF THE TABLE ADJUSTING THE VALUE
C           OF N.
C
        IF(EPSINF.GT.0.1D-03) GO TO 30
   20   N = I+I-1
C ***JUMP OUT OF DO-LOOP
        GO TO 50
C
C           COMPUTE A NEW ELEMENT AND EVENTUALLY ADJUST
C           THE VALUE OF RESULT.
C
   30   RES = E1+0.1D+01/SS
        EPSTAB(K1) = RES
        K1 = K1-2
        ERROR = ERR2+DABS(RES-E2)+ERR3
        IF(ERROR.GT.ABSERR) GO TO 40
        ABSERR = ERROR
        RESULT = RES
   40 CONTINUE
C
C           SHIFT THE TABLE.
C
   50 IF(N.EQ.LIMEXP) N = 2*(LIMEXP/2)-1
      IB = 1
      IF((NUM/2)*2.EQ.NUM) IB = 2
      IE = NEWELM+1
      DO 60 I=1,IE
        IB2 = IB+2
        EPSTAB(IB) = EPSTAB(IB2)
        IB = IB2
   60 CONTINUE
      IF(NUM.EQ.N) GO TO 80
      INDX = NUM-N+1
      DO 70 I = 1,N
        EPSTAB(I)= EPSTAB(INDX)
        INDX = INDX+1
   70 CONTINUE
   80 IF(NRES.GE.4) GO TO 90
      RES3LA(NRES) = RESULT
      ABSERR = OFLOW
      GO TO 100
C
C           COMPUTE ERROR ESTIMATE
C
   90 ABSERR = DABS(RESULT-RES3LA(3))+DABS(RESULT-RES3LA(2))
     1  +DABS(RESULT-RES3LA(1))
      RES3LA(1) = RES3LA(2)
      RES3LA(2) = RES3LA(3)
      RES3LA(3) = RESULT
  100 ABSERR = DMAX1(ABSERR,0.5D+01*EPMACH*DABS(RESULT))
      RETURN
      END
*DECK DQPSRT
      SUBROUTINE DQPSRT(LIMIT,LAST,MAXERR,ERMAX,ELIST,IORD,NRMAX)
C***BEGIN PROLOGUE  DQPSRT
C***REFER TO  DQAGE,DQAGIE,DQAGPE,DQAWSE
C***ROUTINES CALLED  (NONE)
C***REVISION DATE  810101   (YYMMDD)
C***KEYWORDS  SEQUENTIAL SORTING
C***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C           DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C***PURPOSE  This routine maintains the descending ordering in the
C            list of the local error estimated resulting from the
C            interval subdivision process. At each call two error
C            estimates are inserted using the sequential search
C            method, top-down for the largest error estimate and
C            bottom-up for the smallest error estimate.
C***DESCRIPTION
C
C           Ordering routine
C           Standard fortran subroutine
C           Double precision version
C
C           PARAMETERS (MEANING AT OUTPUT)
C              LIMIT  - Integer
C                       Maximum number of error estimates the list
C                       can contain
C
C              LAST   - Integer
C                       Number of error estimates currently in the list
C
C              MAXERR - Integer
C                       Maxerr points to the NRMAX-th largest error
C                       estimate currently in the list
C
C              ERMAX  - Double precision
C                       NRMAX-th largest error estimate
C                       ERMAX = ELIST(MAXERR)
C
C              ELIST  - Double precision
C                       Vector of dimension LAST containing
C                       the error estimates
C
C              IORD   - Integer
C                       Vector of dimension LAST, the first K elements
C                       of which contain pointers to the error
C                       estimates, such that
C                       ELIST(IORD(1)),...,  ELIST(IORD(K))
C                       form a decreasing sequence, with
C                       K = LAST if LAST.LE.(LIMIT/2+2), and
C                       K = LIMIT+1-LAST otherwise
C
C              NRMAX  - Integer
C                       MAXERR = IORD(NRMAX)
C***END PROLOGUE  DQPSRT
C
      DOUBLE PRECISION ELIST,ERMAX,ERRMAX,ERRMIN
      INTEGER I,IBEG,IDO,IORD,ISUCC,J,JBND,JUPBN,K,LAST,LIMIT,MAXERR,
     1  NRMAX
      DIMENSION ELIST(LAST),IORD(LAST)
C
C           CHECK WHETHER THE LIST CONTAINS MORE THAN
C           TWO ERROR ESTIMATES.
C
C***FIRST EXECUTABLE STATEMENT  DQPSRT
      IF(LAST.GT.2) GO TO 10
      IORD(1) = 1
      IORD(2) = 2
      GO TO 90
C
C           THIS PART OF THE ROUTINE IS ONLY EXECUTED IF, DUE TO A
C           DIFFICULT INTEGRAND, SUBDIVISION INCREASED THE ERROR
C           ESTIMATE. IN THE NORMAL CASE THE INSERT PROCEDURE SHOULD
C           START AFTER THE NRMAX-TH LARGEST ERROR ESTIMATE.
C
   10 ERRMAX = ELIST(MAXERR)
      IF(NRMAX.EQ.1) GO TO 30
      IDO = NRMAX-1
      DO 20 I = 1,IDO
        ISUCC = IORD(NRMAX-1)
C ***JUMP OUT OF DO-LOOP
        IF(ERRMAX.LE.ELIST(ISUCC)) GO TO 30
        IORD(NRMAX) = ISUCC
        NRMAX = NRMAX-1
   20    CONTINUE
C
C           COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO BE MAINTAINED
C           IN DESCENDING ORDER. THIS NUMBER DEPENDS ON THE NUMBER OF
C           SUBDIVISIONS STILL ALLOWED.
C
   30 JUPBN = LAST
      IF(LAST.GT.(LIMIT/2+2)) JUPBN = LIMIT+3-LAST
      ERRMIN = ELIST(LAST)
C
C           INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN,
C           STARTING COMPARISON FROM THE ELEMENT ELIST(IORD(NRMAX+1)).
C
      JBND = JUPBN-1
      IBEG = NRMAX+1
      IF(IBEG.GT.JBND) GO TO 50
      DO 40 I=IBEG,JBND
        ISUCC = IORD(I)
C ***JUMP OUT OF DO-LOOP
        IF(ERRMAX.GE.ELIST(ISUCC)) GO TO 60
        IORD(I-1) = ISUCC
   40 CONTINUE
   50 IORD(JBND) = MAXERR
      IORD(JUPBN) = LAST
      GO TO 90
C
C           INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP.
C
   60 IORD(I-1) = MAXERR
      K = JBND
      DO 70 J=I,JBND
        ISUCC = IORD(K)
C ***JUMP OUT OF DO-LOOP
        IF(ERRMIN.LT.ELIST(ISUCC)) GO TO 80
        IORD(K+1) = ISUCC
        K = K-1
   70 CONTINUE
      IORD(I) = LAST
      GO TO 90
   80 IORD(K+1) = LAST
C
C           SET MAXERR AND ERMAX.
C
   90 MAXERR = IORD(NRMAX)
      ERMAX = ELIST(MAXERR)
      RETURN
      END







