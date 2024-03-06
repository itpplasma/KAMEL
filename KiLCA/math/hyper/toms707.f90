MODULE Confluent_Hypergeometric

!      ALGORITHM 707, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 3, SEPTEMBER, 1992, PP. 345-349.
!     ****************************************************************
!     *                                                              *
!     *      SOLUTION TO THE CONFLUENT HYPERGEOMETRIC FUNCTION       *
!     *                                                              *
!     *                           by                                 *
!     *                                                              *
!     *                      MARK NARDIN,                            *
!     *                                                              *
!     *              W. F. PERGER and ATUL BHALLA                    *
!     *                                                              *
!     *                                                              *
!     *  Michigan Technological University, Copyright 1989           *
!     *                                                              *
!     *                                                              *
!     *  Description : A numerical evaluator for the confluent       *
!     *    hypergeometric function for complex arguments with large  *
!     *    magnitudes using a direct summation of the Kummer series. *
!     *    The method used allows an accuracy of up to thirteen      *
!     *    decimal places through the use of large real arrays       *
!     *    and a single final division.  LNCHF is a variable which   *
!     *    selects how the result should be represented.  A '0' will *
!     *    return the value in standard exponential form.  A '1'     *
!     *    will return the LOG of the result.  IP is an integer      *
!     *    variable that specifies how many array positions are      *
!     *    desired (usually 10 is sufficient).  Setting IP=0 causes  *
!     *    the program to estimate the number of array positions.    *
!     *                                                              *
!     *    The confluent hypergeometric function is the solution to  *
!     *    the differential equation:                                *
!     *                                                              *
!     *             zf"(z) + (a-z)f'(z) - bf(z) = 0                  *
!     *                                                              *
!     *  Subprograms called: BITS, CHGF                              *
!     *                                                              *
!     ****************************************************************


! Code converted using TO_F90 by Alan Miller
! Date: 2002-03-12  Time: 12:42:32

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

CONTAINS


FUNCTION conhyp(a, b, z, lnchf, ip) RESULT(fn_val)

COMPLEX (dp), INTENT(IN)  :: a
COMPLEX (dp), INTENT(IN)  :: b
COMPLEX (dp), INTENT(IN)  :: z
INTEGER, INTENT(IN)       :: lnchf
INTEGER, INTENT(IN)       :: ip
COMPLEX (dp)              :: fn_val

INTEGER    :: i
REAL (dp)  :: nterm, fx, term1, MAX, term2, ang

IF (ABS(z) /= 0.0D0) THEN
  ang = ATAN2(AIMAG(z), REAL(z, KIND=dp))
ELSE
  ang = 1.0D0
END IF
IF (ABS(ang) < 3.14159265358979_dp*0.5) THEN
  ang = 1.0D0
ELSE
  ang = SIN(ABS(ang) - 3.14159265358979_dp*0.5D0) + 1.0D0
END IF
MAX = 0
nterm = 0
fx = 0
term1 = 0

10 nterm = nterm + 1
term2 = ABS((a+nterm-1)*z/((b+nterm-1)*nterm))
IF (term2 /= 0.0D0) THEN
  IF (term2 < 1.0D0) THEN
    IF ((DBLE(a)+nterm-1) > 1.0D0) THEN
      IF ((DBLE(b)+nterm-1) > 1.0D0) THEN
        IF (term2-term1 < 0.0D0) THEN
          GO TO 20
        END IF
      END IF
    END IF
  END IF
  fx = fx + LOG(term2)
  IF (fx > MAX) MAX = fx
  term1 = term2
  GO TO 10
END IF

20 MAX = MAX * 2 / (digits(0.0_dp)*6.93147181D-1)
i = INT(MAX*ang) + 7
IF (i < 5) i = 5
IF (ip > i) i = ip
fn_val = chgf(a, b, z, i, lnchf)
RETURN
END FUNCTION conhyp



!     ****************************************************************
!     *                                                              *
!     *                   FUNCTION CHGF                              *
!     *                                                              *
!     *                                                              *
!     *  Description : Function that sums the Kummer series and      *
!     *    returns the solution of the confluent hypergeometric      *
!     *    function.                                                 *
!     *                                                              *
!     *  Subprograms called: ARMULT, ARYDIV, BITS, CMPADD, CMPMUL    *
!     *                                                              *
!     ****************************************************************

FUNCTION chgf(a, b, z, l, lnchf) RESULT(fn_val)

COMPLEX (dp), INTENT(IN)  :: a
COMPLEX (dp), INTENT(IN)  :: b
COMPLEX (dp), INTENT(IN)  :: z
INTEGER, INTENT(IN)       :: l
INTEGER, INTENT(IN)       :: lnchf
COMPLEX (dp)              :: fn_val

INTEGER, PARAMETER  :: length = 777
INTEGER       :: i, bit
COMPLEX (dp)  :: final
REAL (dp)     :: ar, ai, cr, ci, xr, xi, cnt, sigfig, mx1, mx2, rmax
REAL (dp)     :: ar2, ai2, cr2, ci2, xr2, xi2
REAL (dp)     :: sumr(-1:length), sumi(-1:length), numr(-1:length)
REAL (dp)     :: numi(-1:length), denomr(-1:length), denomi(-1:length)
REAL (dp)     :: qr1(-1:length), qr2(-1:length), qi1(-1:length), qi2(-1:length)

bit = digits(0.0_dp)
rmax = 2.0D0 ** (bit/2)
sigfig = 2.0D0 ** (bit/4)
ar2 = DBLE(a) * sigfig
ar = AINT(ar2)
ar2 = ANINT((ar2-ar)*rmax)
ai2 = AIMAG(a) * sigfig
ai = AINT(ai2)
ai2 = ANINT((ai2-ai)*rmax)
cr2 = DBLE(b) * sigfig
cr = AINT(cr2)
cr2 = ANINT((cr2-cr)*rmax)
ci2 = AIMAG(b) * sigfig
ci = AINT(ci2)
ci2 = ANINT((ci2-ci)*rmax)
xr2 = DBLE(z) * sigfig
xr = AINT(xr2)
xr2 = ANINT((xr2-xr)*rmax)
xi2 = AIMAG(z) * sigfig
xi = AINT(xi2)
xi2 = ANINT((xi2-xi)*rmax)
sumr(-1) = 1.0D0
sumi(-1) = 1.0D0
numr(-1) = 1.0D0
numi(-1) = 1.0D0
denomr(-1) = 1.0D0
denomi(-1) = 1.0D0
DO  i = 0, l + 1
  sumr(i) = 0.0D0
  sumi(i) = 0.0D0
  numr(i) = 0.0D0
  numi(i) = 0.0D0
  denomr(i) = 0.0D0
  denomi(i) = 0.0D0
END DO
sumr(1) = 1.0D0
numr(1) = 1.0D0
denomr(1) = 1.0D0
cnt = sigfig
20 IF (sumr(1) < 0.5) THEN
  mx1 = sumi(l+1)
ELSE IF (sumi(1) < 0.5) THEN
  mx1 = sumr(l+1)
ELSE
  mx1 = MAX(sumr(l+1),sumi(l+1))
END IF
IF (numr(1) < 0.5) THEN
  mx2 = numi(l+1)
ELSE IF (numi(1) < 0.5) THEN
  mx2 = numr(l+1)
ELSE
  mx2 = MAX(numr(l+1),numi(l+1))
END IF
IF (mx1-mx2 > 2.0) THEN
  IF (cr > 0.0D0) THEN
    IF (ABS(CMPLX(ar, ai, KIND=dp)*CMPLX(xr, xi, KIND=dp) /  &
            (CMPLX(cr, ci, KIND=dp)*cnt)) <= 1.0D0) GO TO 30
  END IF
END IF
CALL cmpmul(sumr, sumi, cr, ci, qr1, qi1, l, rmax)
CALL cmpmul(sumr, sumi, cr2, ci2, qr2, qi2, l, rmax)
qr2(l+1) = qr2(l+1) - 1
qi2(l+1) = qi2(l+1) - 1
CALL cmpadd(qr1, qi1, qr2, qi2, sumr, sumi, l, rmax)

CALL armult(sumr, cnt, sumr, l, rmax)
CALL armult(sumi, cnt, sumi, l, rmax)
CALL cmpmul(denomr, denomi, cr, ci, qr1, qi1, l, rmax)
CALL cmpmul(denomr, denomi, cr2, ci2, qr2, qi2, l, rmax)
qr2(l+1) = qr2(l+1) - 1
qi2(l+1) = qi2(l+1) - 1
CALL cmpadd(qr1, qi1, qr2, qi2, denomr, denomi, l, rmax)

CALL armult(denomr, cnt, denomr, l, rmax)
CALL armult(denomi, cnt, denomi, l, rmax)
CALL cmpmul(numr, numi, ar, ai, qr1, qi1, l, rmax)
CALL cmpmul(numr, numi, ar2, ai2, qr2, qi2, l, rmax)
qr2(l+1) = qr2(l+1) - 1
qi2(l+1) = qi2(l+1) - 1
CALL cmpadd(qr1, qi1, qr2, qi2, numr, numi, l, rmax)

CALL cmpmul(numr, numi, xr, xi, qr1, qi1, l, rmax)
CALL cmpmul(numr, numi, xr2, xi2, qr2, qi2, l, rmax)
qr2(l+1) = qr2(l+1) - 1
qi2(l+1) = qi2(l+1) - 1
CALL cmpadd(qr1, qi1, qr2, qi2, numr, numi, l, rmax)

CALL cmpadd(sumr, sumi, numr, numi, sumr, sumi, l, rmax)
cnt = cnt + sigfig
ar = ar + sigfig
cr = cr + sigfig
GO TO 20
30 CALL arydiv(sumr, sumi, denomr, denomi, final, l, lnchf, rmax, bit)
fn_val = final
RETURN
END FUNCTION chgf


!     ****************************************************************
!     *                                                              *
!     *                 SUBROUTINE ARADD                             *
!     *                                                              *
!     *                                                              *
!     *  Description : Accepts two arrays of numbers and returns     *
!     *    the sum of the array.  Each array is holding the value    *
!     *    of one number in the series.  The parameter L is the      *
!     *    size of the array representing the number and RMAX is     *
!     *    the actual number of digits needed to give the numbers    *
!     *    the desired accuracy.                                     *
!     *                                                              *
!     *  Subprograms called: none                                    *
!     *                                                              *
!     ****************************************************************

SUBROUTINE aradd(a, b, c, l, rmax)

REAL (dp), INTENT(IN)   :: a(-1:)
REAL (dp), INTENT(IN)   :: b(-1:)
REAL (dp), INTENT(OUT)  :: c(-1:)
INTEGER, INTENT(IN)     :: l
REAL (dp), INTENT(IN)   :: rmax

REAL (dp)  :: z(-1:777)
INTEGER    :: ediff, i, j

DO  i = 0, l + 1
  z(i) = 0.0D0
END DO
ediff = ANINT(a(l+1)-b(l+1))
IF (ABS(a(1)) >= 0.5 .AND. ediff > -l) THEN
  IF (ABS(b(1)) < 0.5 .OR. ediff >= l) GO TO 30
  GO TO 50
END IF
DO  i = -1, l + 1
  c(i) = b(i)
END DO
GO TO 300
30 DO  i = -1, l + 1
  c(i) = a(i)
END DO
GO TO 300
50 z(-1) = a(-1)
IF (ABS(a(-1)-b(-1)) >= 0.5) THEN
  IF (ediff > 0) THEN
    z(l+1) = a(l+1)
    GO TO 160
  END IF
  IF (ediff < 0) THEN
    z(l+1) = b(l+1)
    z(-1) = b(-1)
    GO TO 200
  END IF
  DO  i = 1, l
    IF (a(i) > b(i)) THEN
      z(l+1) = a(l+1)
      GO TO 160
    END IF
    IF (a(i) < b(i)) THEN
      z(l+1) = b(l+1)
      z(-1) = b(-1)
      GO TO 200
    END IF
  END DO
  GO TO 280
END IF

IF (ediff <= 0) THEN
  IF (ediff < 0) GO TO 120
  z(l+1) = a(l+1)
  DO  i = l, 1, -1
    z(i) = a(i) + b(i) + z(i)
    IF (z(i) >= rmax) THEN
      z(i) = z(i) - rmax
      z(i-1) = 1.0D0
    END IF
  END DO
  IF (z(0) > 0.5) THEN
    DO  i = l, 1, -1
      z(i) = z(i-1)
    END DO
    z(l+1) = z(l+1) + 1.0D0
    z(0) = 0.0D0
  END IF
  GO TO 280
END IF
z(l+1) = a(l+1)
DO  i = l, 1 + ediff, -1
  z(i) = a(i) + b(i-ediff) + z(i)
  IF (z(i) >= rmax) THEN
    z(i) = z(i) - rmax
    z(i-1) = 1.0D0
  END IF
END DO
DO  i = ediff, 1, -1
  z(i) = a(i) + z(i)
  IF (z(i) >= rmax) THEN
    z(i) = z(i) - rmax
    z(i-1) = 1.0D0
  END IF
END DO
IF (z(0) > 0.5) THEN
  DO  i = l, 1, -1
    z(i) = z(i-1)
  END DO
  z(l+1) = z(l+1) + 1
  z(0) = 0.0D0
END IF
GO TO 280

120 z(l+1) = b(l+1)
DO  i = l, 1 - ediff, -1
  z(i) = a(i+ediff) + b(i) + z(i)
  IF (z(i) >= rmax) THEN
    z(i) = z(i) - rmax
    z(i-1) = 1.0D0
  END IF
END DO
DO  i = 0 - ediff, 1, -1
  z(i) = b(i) + z(i)
  IF (z(i) >= rmax) THEN
    z(i) = z(i) - rmax
    z(i-1) = 1.0D0
  END IF
END DO
IF (z(0) > 0.5) THEN
  DO  i = l, 1, -1
    z(i) = z(i-1)
  END DO
  z(l+1) = z(l+1) + 1.0D0
  z(0) = 0.0D0
END IF
GO TO 280

160 IF (ediff <= 0) THEN
  DO  i = l, 1, -1
    z(i) = a(i) - b(i) + z(i)
    IF (z(i) < 0.0D0) THEN
      z(i) = z(i) + rmax
      z(i-1) = -1.0D0
    END IF
  END DO
  GO TO 240
END IF
DO  i = l, 1 + ediff, -1
  z(i) = a(i) - b(i-ediff) + z(i)
  IF (z(i) < 0.0D0) THEN
    z(i) = z(i) + rmax
    z(i-1) = -1.0D0
  END IF
END DO
DO  i = ediff, 1, -1
  z(i) = a(i) + z(i)
  IF (z(i) < 0.0D0) THEN
    z(i) = z(i) + rmax
    z(i-1) = -1.0D0
  END IF
END DO
GO TO 240

200 IF (ediff >= 0) THEN
  DO  i = l, 1, -1
    z(i) = b(i) - a(i) + z(i)
    IF (z(i) < 0.0D0) THEN
      z(i) = z(i) + rmax
      z(i-1) = -1.0D0
    END IF
  END DO
ELSE
  DO  i = l, 1 - ediff, -1
    z(i) = b(i) - a(i+ediff) + z(i)
    IF (z(i) < 0.0D0) THEN
      z(i) = z(i) + rmax
      z(i-1) = -1.0D0
    END IF
  END DO
  DO  i = 0 - ediff, 1, -1
    z(i) = b(i) + z(i)
    IF (z(i) < 0.0D0) THEN
      z(i) = z(i) + rmax
      z(i-1) = -1.0D0
    END IF
  END DO
END IF

240 IF (z(1) <= 0.5) THEN
  i = 1
  250   i = i + 1
  IF (z(i) < 0.5 .AND. i < l+1) GO TO 250
  IF (i == l+1) THEN
    z(-1) = 1.0D0
    z(l+1) = 0.0D0
    GO TO 280
  END IF
  DO  j = 1, l + 1 - i
    z(j) = z(j+i-1)
  END DO
  DO  j = l + 2 - i, l
    z(j) = 0.0D0
  END DO
  z(l+1) = z(l+1) - i + 1
END IF
280 DO  i = -1, l + 1
  c(i) = z(i)
END DO
300 IF (c(1) < 0.5) THEN
  c(-1) = 1.0D0
  c(l+1) = 0.0D0
END IF
RETURN
END SUBROUTINE aradd


!     ****************************************************************
!     *                                                              *
!     *                 SUBROUTINE ARSUB                             *
!     *                                                              *
!     *                                                              *
!     *  Description : Accepts two arrays and subtracts each element *
!     *    in the second array from the element in the first array   *
!     *    and returns the solution.  The parameters L and RMAX are  *
!     *    the size of the array and the number of digits needed for *
!     *    the accuracy, respectively.                               *
!     *                                                              *
!     *  Subprograms called: ARADD                                   *
!     *                                                              *
!     ****************************************************************

SUBROUTINE arsub(a,b,c,l,rmax)

REAL (dp), INTENT(IN OUT)  :: a(-1:)
REAL (dp), INTENT(IN)      :: b(-1:)
REAL (dp), INTENT(IN OUT)  :: c(-1:)
INTEGER, INTENT(IN)        :: l
REAL (dp), INTENT(IN OUT)  :: rmax

INTEGER :: i
REAL (dp) :: b2(-1:777)

DO  i = -1, l + 1
  b2(i) = b(i)
END DO
b2(-1) = (-1.0D0) * b2(-1)
CALL aradd(a, b2, c, l, rmax)
RETURN
END SUBROUTINE arsub


!     ****************************************************************
!     *                                                              *
!     *                 SUBROUTINE ARMULT                            *
!     *                                                              *
!     *                                                              *
!     *  Description : Accepts two arrays and returns the product.   *
!     *    L and RMAX are the size of the arrays and the number of   *
!     *    digits needed to represent the numbers with the required  *
!     *    accuracy.                                                 *
!     *                                                              *
!     *  Subprograms called: none                                    *
!     *                                                              *
!     ****************************************************************

SUBROUTINE armult(a, b, c, l, rmax)

REAL (dp), INTENT(IN)      :: a(-1:)
REAL (dp), INTENT(IN OUT)  :: b
REAL (dp), INTENT(OUT)     :: c(-1:)
INTEGER, INTENT(IN)        :: l
REAL (dp), INTENT(IN)      :: rmax

REAL (dp)  :: z(-1:777), b2, carry
INTEGER    :: i

z(-1) = SIGN(1.0D0,b) * a(-1)
b2 = ABS(b)
z(l+1) = a(l+1)
DO  i = 0, l
  z(i) = 0.0D0
END DO
IF (b2 <= 1.0D-10 .OR. a(1) <= 1.0D-10) THEN
  z(-1) = 1.0D0
  z(l+1) = 0.0D0
  GO TO 40
END IF
DO  i = l, 1, -1
  z(i) = a(i) * b2 + z(i)
  IF (z(i) >= rmax) THEN
    carry = AINT(z(i)/rmax)
    z(i) = z(i) - carry * rmax
    z(i-1) = carry
  END IF
END DO
IF (z(0) >= 0.5) THEN
  DO  i = l, 1, -1
    z(i) = z(i-1)
  END DO
  z(l+1) = z(l+1) + 1.0D0
  z(0) = 0.0D0
END IF


40 DO  i = -1, l + 1
  c(i) = z(i)
END DO
IF (c(1) < 0.5) THEN
  c(-1) = 1.0D0
  c(l+1) = 0.0D0
END IF
RETURN
END SUBROUTINE armult


!     ****************************************************************
!     *                                                              *
!     *                 SUBROUTINE CMPADD                            *
!     *                                                              *
!     *                                                              *
!     *  Description : Takes two arrays representing one real and    *
!     *    one imaginary part, and adds two arrays representing      *
!     *    another complex number and returns two array holding the  *
!     *    complex sum.                                              *
!     *              (CR,CI) = (AR+BR, AI+BI)                        *
!     *                                                              *
!     *  Subprograms called: ARADD                                   *
!     *                                                              *
!     ****************************************************************

SUBROUTINE cmpadd(ar,ai,br,bi,cr,ci,l,rmax)

REAL (dp), INTENT(IN OUT)  :: ar(-1:)
REAL (dp), INTENT(IN OUT)  :: ai(-1:)
REAL (dp), INTENT(IN OUT)  :: br(-1:)
REAL (dp), INTENT(IN OUT)  :: bi(-1:)
REAL (dp), INTENT(IN OUT)  :: cr(-1:)
REAL (dp), INTENT(IN OUT)  :: ci(-1:)
INTEGER, INTENT(IN)        :: l
REAL (dp), INTENT(IN OUT)  :: rmax

CALL aradd(ar,br,cr,l,rmax)
CALL aradd(ai,bi,ci,l,rmax)
RETURN
END SUBROUTINE cmpadd


!     ****************************************************************
!     *                                                              *
!     *                 SUBROUTINE CMPSUB                            *
!     *                                                              *
!     *                                                              *
!     *  Description : Takes two arrays representing one real and    *
!     *    one imaginary part, and subtracts two arrays representing *
!     *    another complex number and returns two array holding the  *
!     *    complex sum.                                              *
!     *              (CR,CI) = (AR+BR, AI+BI)                        *
!     *                                                              *
!     *  Subprograms called: ARADD                                   *
!     *                                                              *
!     ****************************************************************

SUBROUTINE cmpsub(ar,ai,br,bi,cr,ci,l,rmax)

REAL (dp), INTENT(IN OUT)  :: ar(-1:)
REAL (dp), INTENT(IN OUT)  :: ai(-1:)
REAL (dp), INTENT(IN OUT)  :: br(-1:)
REAL (dp), INTENT(IN OUT)  :: bi(-1:)
REAL (dp), INTENT(IN OUT)  :: cr(-1:)
REAL (dp), INTENT(IN OUT)  :: ci(-1:)
INTEGER, INTENT(IN OUT)           :: l
REAL (dp), INTENT(IN OUT)  :: rmax

CALL arsub(ar,br,cr,l,rmax)
CALL arsub(ai,bi,ci,l,rmax)
RETURN
END SUBROUTINE cmpsub


!     ****************************************************************
!     *                                                              *
!     *                 SUBROUTINE CMPMUL                            *
!     *                                                              *
!     *                                                              *
!     *  Description : Takes two arrays representing one real and    *
!     *    one imaginary part, and multiplies it with two arrays     *
!     *    representing another complex number and returns the       *
!     *    complex product.                                          *
!     *                                                              *
!     *  Subprograms called: ARMULT, ARSUB, ARADD                    *
!     *                                                              *
!     ****************************************************************

SUBROUTINE cmpmul(ar,ai,br,bi,cr,ci,l,rmax)

REAL (dp), INTENT(IN OUT)  :: ar(-1:)
REAL (dp), INTENT(IN OUT)  :: ai(-1:)
REAL (dp), INTENT(IN OUT)  :: br
REAL (dp), INTENT(IN OUT)  :: bi
REAL (dp), INTENT(IN OUT)  :: cr(-1:)
REAL (dp), INTENT(IN OUT)  :: ci(-1:)
INTEGER, INTENT(IN)        :: l
REAL (dp), INTENT(IN OUT)  :: rmax

REAL (dp)  :: d1(-1:777), d2(-1:777)

CALL armult(ar, br, d1, l, rmax)
CALL armult(ai, bi, d2, l, rmax)
CALL arsub(d1, d2, cr, l, rmax)
CALL armult(ar, bi, d1, l, rmax)
CALL armult(ai, br, d2, l, rmax)
CALL aradd(d1, d2, ci, l, rmax)
RETURN
END SUBROUTINE cmpmul


!     ****************************************************************
!     *                                                              *
!     *                 SUBROUTINE ARYDIV                            *
!     *                                                              *
!     *                                                              *
!     *  Description : Returns the REAL (dp) complex number          *
!     *    resulting from the division of four arrays, representing  *
!     *    two complex numbers.  The number returned will be in one  *
!     *    two different forms.  Either standard scientific or as    *
!     *    the log of the number.                                    *
!     *                                                              *
!     *  Subprograms called: CONV21, CONV12, EADD, ECPDIV, EMULT     *
!     *                                                              *
!     ****************************************************************

SUBROUTINE arydiv(ar, ai, br, bi, c, l, lnchf, rmax, bit)

REAL (dp), INTENT(IN)      :: ar(-1:)
REAL (dp), INTENT(IN)      :: ai(-1:)
REAL (dp), INTENT(IN)      :: br(-1:)
REAL (dp), INTENT(IN)      :: bi(-1:)
COMPLEX (dp), INTENT(OUT)  :: c
INTEGER, INTENT(IN)        :: l
INTEGER, INTENT(IN)        :: lnchf
REAL (dp), INTENT(IN)      :: rmax
INTEGER, INTENT(IN)        :: bit

INTEGER    :: rexp, ir10, ii10
REAL (dp)  :: phi, n1, n2, n3, e1, e2, e3, rr10, ri10, x
REAL (dp)  :: x1, x2, dum1, dum2
REAL (dp)  :: ae(2,2), be(2,2), ce(2,2)

rexp = bit / 2
x = rexp * (ar(l+1)-2)
rr10 = x * LOG10(2.0D0) / LOG10(10.0D0)
ir10 = INT(rr10)
rr10 = rr10 - ir10
x = rexp * (ai(l+1)-2)
ri10 = x * LOG10(2.0D0) / LOG10(10.0D0)
ii10 = ri10
ri10 = ri10 - ii10
dum1 = SIGN(ar(1)*rmax*rmax + ar(2)*rmax + ar(3), ar(-1))
dum2 = SIGN(ai(1)*rmax*rmax + ai(2)*rmax + ai(3), ai(-1))
dum1 = dum1 * 10 ** rr10
dum2 = dum2 * 10 ** ri10
CALL conv12(CMPLX(dum1, dum2, KIND=dp), ae)
ae(1,2) = ae(1,2) + ir10
ae(2,2) = ae(2,2) + ii10
x = rexp * (br(l+1)-2)
rr10 = x * LOG10(2.0D0) / LOG10(10.0D0)
ir10 = INT(rr10)
rr10 = rr10 - ir10
x = rexp * (bi(l+1)-2)
ri10 = x * LOG10(2.0D0) / LOG10(10.0D0)
ii10 = INT(ri10)
ri10 = ri10 - ii10
dum1 = SIGN(br(1)*rmax*rmax + br(2)*rmax + br(3), br(-1))
dum2 = SIGN(bi(1)*rmax*rmax + bi(2)*rmax + bi(3), bi(-1))
dum1 = dum1 * 10 ** rr10
dum2 = dum2 * 10 ** ri10
CALL conv12(CMPLX(dum1, dum2, KIND=dp), be)
be(1,2) = be(1,2) + ir10
be(2,2) = be(2,2) + ii10
CALL ecpdiv(ae, be, ce)
IF (lnchf == 0) THEN
  CALL conv21(ce, c)
ELSE
  CALL emult(ce(1,1), ce(1,2), ce(1,1), ce(1,2), n1, e1)
  CALL emult(ce(2,1), ce(2,2), ce(2,1), ce(2,2), n2, e2)
  CALL eadd(n1, e1, n2, e2, n3, e3)
  n1 = ce(1,1)
  e1 = ce(1,2) - ce(2,2)
  x2 = ce(2,1)
  IF (e1 > 74.0D0) THEN
    x1 = 1.0D75
  ELSE IF (e1 < -74.0D0) THEN
    x1 = 0
  ELSE
    x1 = n1 * (10**e1)
  END IF
  phi = ATAN2(x2,x1)
  c = CMPLX(0.50D0*(LOG(n3) + e3*LOG(10.0D0)), phi, KIND=dp)
END IF
RETURN
END SUBROUTINE arydiv


!     ****************************************************************
!     *                                                              *
!     *                 SUBROUTINE EMULT                             *
!     *                                                              *
!     *                                                              *
!     *  Description : Takes one base and exponent and multiplies it *
!     *    by another numbers base and exponent to give the product  *
!     *    in the form of base and exponent.                         *
!     *                                                              *
!     *  Subprograms called: none                                    *
!     *                                                              *
!     ****************************************************************

SUBROUTINE emult(n1, e1, n2, e2, nf, ef)

REAL (dp), INTENT(IN)   :: n1
REAL (dp), INTENT(IN)   :: e1
REAL (dp), INTENT(IN)   :: n2
REAL (dp), INTENT(IN)   :: e2
REAL (dp), INTENT(OUT)  :: nf
REAL (dp), INTENT(OUT)  :: ef

nf = n1 * n2
ef = e1 + e2
IF (ABS(nf) >= 10.0D0) THEN
  nf = nf / 10.0D0
  ef = ef + 1.0D0
END IF
RETURN
END SUBROUTINE emult


!     ****************************************************************
!     *                                                              *
!     *                 SUBROUTINE EDIV                              *
!     *                                                              *
!     *                                                              *
!     *  Description : returns the solution in the form of base and  *
!     *    exponent of the division of two exponential numbers.      *
!     *                                                              *
!     *  Subprograms called: none                                    *
!     *                                                              *
!     ****************************************************************

SUBROUTINE ediv(n1, e1, n2, e2, nf, ef)

REAL (dp), INTENT(IN)   :: n1
REAL (dp), INTENT(IN)   :: e1
REAL (dp), INTENT(IN)   :: n2
REAL (dp), INTENT(IN)   :: e2
REAL (dp), INTENT(OUT)  :: nf
REAL (dp), INTENT(OUT)  :: ef

nf = n1 / n2
ef = e1 - e2
IF (ABS(nf) < 1.0D0 .AND. nf /= 0.0D0) THEN
  nf = nf * 10.0D0
  ef = ef - 1.0D0
END IF
RETURN
END SUBROUTINE ediv


!     ****************************************************************
!     *                                                              *
!     *                 SUBROUTINE EADD                              *
!     *                                                              *
!     *                                                              *
!     *  Description : Returns the sum of two numbers in the form    *
!     *    of a base and an exponent.                                *
!     *                                                              *
!     *  Subprograms called: none                                    *
!     *                                                              *
!     ****************************************************************

SUBROUTINE eadd(n1, e1, n2, e2, nf, ef)

REAL (dp), INTENT(IN)   :: n1
REAL (dp), INTENT(IN)   :: e1
REAL (dp), INTENT(IN)   :: n2
REAL (dp), INTENT(IN)   :: e2
REAL (dp), INTENT(OUT)  :: nf
REAL (dp), INTENT(OUT)  :: ef

REAL (dp) :: ediff

ediff = e1 - e2
IF (ediff > 36.0D0) THEN
  nf = n1
  ef = e1
ELSE IF (ediff < -36.0D0) THEN
  nf = n2
  ef = e2
ELSE
  nf = n1 * (10.0D0**ediff) + n2
  ef = e2
  10   IF (ABS(nf) >= 10.0D0) THEN
    nf = nf / 10.0D0
    ef = ef + 1.0D0
    GO TO 10
  END IF
  20 IF (ABS(nf) >= 1.0D0 .OR. nf == 0.0D0) GO TO 30
  nf = nf * 10.0D0
  ef = ef - 1.0D0
  GO TO 20
END IF
30 RETURN
END SUBROUTINE eadd


!     ****************************************************************
!     *                                                              *
!     *                 SUBROUTINE ESUB                              *
!     *                                                              *
!     *                                                              *
!     *  Description : Returns the solution to the subtraction of    *
!     *    two numbers in the form of base and exponent.             *
!     *                                                              *
!     *  Subprograms called: EADD                                    *
!     *                                                              *
!     ****************************************************************

SUBROUTINE esub(n1, e1, n2, e2, nf, ef)

REAL (dp), INTENT(IN OUT)  :: n1
REAL (dp), INTENT(IN OUT)  :: e1
REAL (dp), INTENT(IN OUT)  :: n2
REAL (dp), INTENT(IN OUT)  :: e2
REAL (dp), INTENT(IN OUT)  :: nf
REAL (dp), INTENT(IN OUT)  :: ef

CALL eadd(n1, e1, DBLE(-n2), e2, nf, ef)
RETURN
END SUBROUTINE esub


!     ****************************************************************
!     *                                                              *
!     *                 SUBROUTINE CONV12                            *
!     *                                                              *
!     *                                                              *
!     *  Description : Converts a number from complex notation to a  *
!     *    form of a 2x2 real array.                                 *
!     *                                                              *
!     *  Subprograms called: none                                    *
!     *                                                              *
!     ****************************************************************

SUBROUTINE conv12(cn, cae)

COMPLEX (dp), INTENT(IN)  :: cn
REAL (dp), INTENT(OUT)    :: cae(2,2)

cae(1,1) = DBLE(cn)
cae(1,2) = 0.0D0
10 IF (ABS(cae(1,1)) >= 10.0D0) THEN
  cae(1,1) = cae(1,1) / 10.0D0
  cae(1,2) = cae(1,2) + 1.0D0
  GO TO 10
END IF
20 IF (.NOT.(ABS(cae(1,1)) >= 1.0D0 .OR. cae(1,1) == 0.0D0)) THEN
  cae(1,1) = cae(1,1) * 10.0D0
  cae(1,2) = cae(1,2) - 1.0D0
  GO TO 20
END IF
cae(2,1) = AIMAG(cn)
cae(2,2) = 0.0D0
30 IF (ABS(cae(2,1)) >= 10.0D0) THEN
  cae(2,1) = cae(2,1) / 10.0D0
  cae(2,2) = cae(2,2) + 1.0D0
  GO TO 30
END IF
40 IF (.NOT.(ABS(cae(2,1)) >= 1.0D0 .OR. cae(2,1) == 0.0D0)) THEN
  cae(2,1) = cae(2,1) * 10.0D0
  cae(2,2) = cae(2,2) - 1.0D0
  GO TO 40
END IF
RETURN
END SUBROUTINE conv12


!     ****************************************************************
!     *                                                              *
!     *                 SUBROUTINE CONV21                            *
!     *                                                              *
!     *                                                              *
!     *  Description : Converts a number represented in a 2x2 real   *
!     *    array to the form of a complex number.                    *
!     *                                                              *
!     *  Subprograms called: none                                    *
!     *                                                              *
!     ****************************************************************

SUBROUTINE conv21(cae, cn)

REAL (dp), INTENT(IN)      :: cae(2,2)
COMPLEX (dp), INTENT(OUT)  :: cn

IF (cae(1,2) > 75 .OR. cae(2,2) > 75) THEN
  cn = CMPLX(1.0D75, 1.0D75, KIND=dp)
ELSE IF (cae(2,2) < -75) THEN
  cn = CMPLX(cae(1,1)*(10**cae(1,2)), 0D0, KIND=dp)
ELSE
  cn = CMPLX(cae(1,1)*(10**cae(1,2)), cae(2,1)*(10**cae(2,2)), KIND=dp)
END IF
RETURN
END SUBROUTINE conv21


!     ****************************************************************
!     *                                                              *
!     *                 SUBROUTINE ECPMUL                            *
!     *                                                              *
!     *                                                              *
!     *  Description : Multiplies two numbers which are each         *
!     *    represented in the form of a two by two array and returns *
!     *    the solution in the same form.                            *
!     *                                                              *
!     *  Subprograms called: EMULT, ESUB, EADD                       *
!     *                                                              *
!     ****************************************************************

SUBROUTINE ecpmul(a,b,c)

REAL (dp), INTENT(IN)   :: a(2,2)
REAL (dp), INTENT(IN)   :: b(2,2)
REAL (dp), INTENT(OUT)  :: c(2,2)

REAL (dp) :: n1, e1, n2, e2, c2(2,2)

CALL emult(a(1,1), a(1,2), b(1,1), b(1,2), n1, e1)
CALL emult(a(2,1), a(2,2), b(2,1), b(2,2), n2, e2)
CALL esub(n1, e1, n2, e2, c2(1,1), c2(1,2))
CALL emult(a(1,1), a(1,2), b(2,1), b(2,2), n1, e1)
CALL emult(a(2,1), a(2,2), b(1,1), b(1,2), n2, e2)
CALL eadd(n1, e1, n2, e2, c(2,1), c(2,2))
c(1,1) = c2(1,1)
c(1,2) = c2(1,2)
RETURN
END SUBROUTINE ecpmul


!     ****************************************************************
!     *                                                              *
!     *                 SUBROUTINE ECPDIV                            *
!     *                                                              *
!     *                                                              *
!     *  Description : Divides two numbers and returns the solution. *
!     *    All numbers are represented by a 2x2 array.               *
!     *                                                              *
!     *  Subprograms called: EADD, ECPMUL, EDIV, EMULT               *
!     *                                                              *
!     ****************************************************************

SUBROUTINE ecpdiv(a,b,c)

REAL (dp), INTENT(IN)   :: a(2,2)
REAL (dp), INTENT(IN)   :: b(2,2)
REAL (dp), INTENT(OUT)  :: c(2,2)

REAL (dp)  :: n1, e1, n2, e2, b2(2,2), n3, e3, c2(2,2)

b2(1,1) = b(1,1)
b2(1,2) = b(1,2)
b2(2,1) = -1.0D0 * b(2,1)
b2(2,2) = b(2,2)
CALL ecpmul(a, b2, c2)
CALL emult(b(1,1), b(1,2), b(1,1), b(1,2), n1, e1)
CALL emult(b(2,1), b(2,2), b(2,1), b(2,2), n2, e2)
CALL eadd(n1, e1, n2, e2, n3, e3)
CALL ediv(c2(1,1), c2(1,2), n3, e3, c(1,1), c(1,2))
CALL ediv(c2(2,1), c2(2,2), n3, e3, c(2,1), c(2,2))
RETURN
END SUBROUTINE ecpdiv

END MODULE Confluent_Hypergeometric




PROGRAM sample
USE Confluent_Hypergeometric
IMPLICIT NONE
COMPLEX (dp)  :: a, b, z, chf
INTEGER       :: lnchf, ip
REAL (dp)     :: ar, ai

OPEN (5, FILE='SAMPLE.DAT', STATUS='OLD')
READ (5,*) a
WRITE (6,*) ' A= ', a
READ (5,*) b
WRITE (6,*) ' B= ', b
READ (5,*) z
WRITE (6,*) ' Z= ', z
READ (5,*) lnchf
WRITE (6,*) ' LNCHF= ', lnchf
READ (5,*) ip
WRITE (6,*) ' IP= ', ip
chf = conhyp(a, b, z, lnchf, ip)
WRITE (6,5000) chf
ar = 2.31145634403113D-12
ai = -1.96169649634905D-11
WRITE (6,5100) ar, ai
STOP

5000 FORMAT (' RESULT FROM CONHYP=', 2g25.12)
5100 FORMAT (' EXPECTED RESULT=   ', 2g25.12)
END PROGRAM sample



! Contents of the file: SAMPLE.DAT

(-15.0D0,55.0D0)
(20.0D0,25.0D0)
(-100.0D0,200.0D0)
0
10
