C     FILE SLMACH.FOR
C
C     SLATEC LIBRARY MACHINE DEPENDENT ROUTINES
C
C     CONTENTS:
C     I1MACH - INTEGER MACHINE CONSTANTS
C     D1MACH - DOUBLE PRECISION MACHINE CONSTANTS
C     XERROR - ERROR HANDLING
C
C     NOTE: NOT ALL MACHINE CONSTANTS HAVE BEEN SET
C
C---------------------------------------------------------------------
C
*DECK I1MACH
      INTEGER FUNCTION I1MACH(I)
C***BEGIN PROLOGUE  I1MACH
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  840405   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns integer machine dependent constants
C***DESCRIPTION
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C   These machine constant routines must be activated for
C   a particular environment.
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     I1MACH can be used to obtain machine-dependent parameters
C     for the local machine environment.  It is a function
C     subroutine with one (input) argument, and can be called
C     as follows, for example
C
C          K = I1MACH(I)
C
C     where I=1,...,16.  The (output) value of K above is
C     determined by the (input) value of I.  The results for
C     various values of I are discussed below.
C
C  I/O unit numbers.
C    I1MACH( 1) = the standard input unit.
C    I1MACH( 2) = the standard output unit.
C    I1MACH( 3) = the standard punch unit.
C    I1MACH( 4) = the standard error message unit.
C
C  Words.
C    I1MACH( 5) = the number of bits per integer storage unit.
C    I1MACH( 6) = the number of characters per integer storage unit.
C
C  Integers.
C    assume integers are represented in the S-digit, base-A form
C
C               sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               where 0 .LE. X(I) .LT. A for I=0,...,S-1.
C    I1MACH( 7) = A, the base.
C    I1MACH( 8) = S, the number of base-A digits.
C    I1MACH( 9) = A**S - 1, the largest magnitude.
C
C  Floating-Point Numbers.
C    Assume floating-point numbers are represented in the T-digit,
C    base-B form
C               sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               where 0 .LE. X(I) .LT. B for I=1,...,T,
C               0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
C    I1MACH(10) = B, the base.
C
C  Single-Precision
C    I1MACH(11) = T, the number of base-B digits.
C    I1MACH(12) = EMIN, the smallest exponent E.
C    I1MACH(13) = EMAX, the largest exponent E.
C
C  Double-Precision
C    I1MACH(14) = T, the number of base-B digits.
C    I1MACH(15) = EMIN, the smallest exponent E.
C    I1MACH(16) = EMAX, the largest exponent E.
C
C  To alter this function for a particular environment,
C  the desired set of DATA statements should be activated by
C  removing the C from column 1.  Also, the values of
C  I1MACH(1) - I1MACH(4) should be checked for consistency
C  with the local operating system.
C***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
C                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  I1MACH
C
      INTEGER IMACH(16),OUTPUT
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR THE VAX 11/780
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    5 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) /2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -127 /
C     DATA IMACH(13)/  127 /
C     DATA IMACH(14)/   56 /
C     DATA IMACH(15)/ -127 /
C     DATA IMACH(16)/  127 /
C
C     MACHINE CONSTANTS FOR THE IBM AIX XL FORTRAN COMPILER/6000
C     (IMPLICIT INTEGER, REAL AND DOUBLE PRECISION TYPE)
C
      DATA IMACH( 1) /     5/
      DATA IMACH( 2) /     6/
      DATA IMACH( 3) /     5/
      DATA IMACH( 4) /     0/
      DATA IMACH( 5) /    32/
      DATA IMACH( 6) /     4/
      DATA IMACH( 7) /     2/
      DATA IMACH( 8) /    31/
      DATA IMACH( 9) /2147483647/
      DATA IMACH(10) /     2/
      DATA IMACH(11) /    24/
      DATA IMACH(12) /  -127/
      DATA IMACH(13) /   127/
      DATA IMACH(14) /    52/
      DATA IMACH(15) / -1023/
      DATA IMACH(16) /  1023/
C
C
C***FIRST EXECUTABLE STATEMENT  I1MACH
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10
C
      I1MACH=IMACH(I)
      RETURN
C
   10 CONTINUE
      WRITE(OUTPUT,9000)
9000  FORMAT('1ERROR    1 IN I1MACH - I OUT OF BOUNDS ')
C
C     CALL FDUMP
C
C
      STOP
      END
C
C---------------------------------------------------------------------
C
*DECK D1MACH
      DOUBLE PRECISION FUNCTION D1MACH(I)
C
C     ADAPTED
C
C***BEGIN PROLOGUE  D1MACH
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  831014   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns double precision machine dependent constants
C***DESCRIPTION
C
C     D1MACH can be used to obtain machine-dependent parameters
C     for the local machine environment.  It is a function
C     subprogram with one (input) argument, and can be called
C     as follows, for example
C
C          D = D1MACH(I)
C
C     where I=1,...,5.  The (output) value of D above is
C     determined by the (input) value of I.  The results for
C     various values of I are discussed below.
C
C  Double-precision machine constants
C  D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C  D1MACH( 3) = B**(-T), the smallest relative spacing.
C  D1MACH( 4) = B**(1-T), the largest relative spacing.
C  D1MACH( 5) = LOG10(B)
C***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
C                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  D1MACH
C
      DOUBLE PRECISION DMACH(5)
C
C     MACHINE CONSTANTS FOR PC AND PTDEC1 (ROUTINES RX, MH)
C
C      DATA DMACH(1) /1.D-37/
C      DATA DMACH(2) /1.D+37/
C      DATA DMACH(3) /1.D-16/
C      DATA DMACH(4) /1.D-15/
C      DATA DMACH(5) /1.D0/
C
C     MACHINE CONSTANTS FOR THE IBM AIX XL FORTRAN COMPILER/6000
C     (IMPLICIT DOUBLE PRECISION)
C     OVERFLOW AND UNDERFLOW BOUNDS CONSERVATIVE 
C
      DATA DMACH(1) /1.D-307/
      DATA DMACH(2) /1.D+307/
      DATA DMACH(3) /2.22D-16/
      DATA DMACH(4) /4.44D-16/
      DATA DMACH(5) /0.30103D0/
C
C***FIRST EXECUTABLE STATEMENT  D1MACH
      IF (I .LT. 1  .OR.  I .GT. 5)
     1   CALL XERROR( 'D1MACH -- I OUT OF BOUNDS',25,1,2)
C
      D1MACH = DMACH(I)
      RETURN
C
      END
C
C--------------------------------------------------------------------
C
*DECK XERROR
      SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
C
C     ADAPTED FROM ORIGINAL SLATEC ROUTINE XERROR
C     PRINTS ERROR MESSAGE
C
      CHARACTER MESSG(NMESSG)
C***FIRST EXECUTABLE STATEMENT  XERROR
      PRINT*
      PRINT*, '*** SLATEC XERROR'
      PRINT*, 'ERROR NR ', NERR, '  LEVEL ', LEVEL
      PRINT*, MESSG
      PRINT*
      RETURN
      END
C
C--------------------------------------------------------------------
C
