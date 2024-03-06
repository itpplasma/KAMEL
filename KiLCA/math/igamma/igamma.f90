!-----------------------------------------------------------------------------!

subroutine igamma(s, z, abs_err, rel_err, scal, G, est_err, Nterms, status)

! Evaluates G(s, z) by the continued fractions method
! abs_err, rel_err - desired tolerances
! scal - scaling: if scal == 0 return G, else return exp(z)*G
! G = G(s, z) the function value
! est_err - estimated relative error
! Nterms - number of terms used to reach the requested accuracy
! if status == 0 then Ok
! WARNING: the continued fraction diverges at real negative values of z

implicit none;

complex(8), intent(in)   :: s, z;
real(8),    intent(in)   :: abs_err, rel_err;
integer,    intent(in)   :: scal;
complex(8), intent(out)  :: G;
real(8),    intent(out)  :: est_err;
integer,    intent(out)  :: Nterms;
integer,    intent(out)  :: status;

integer,       parameter :: prec = 8;
integer,       parameter :: Nmin = 2, Nmax = 4194304;
complex(prec), parameter :: NaN = cmplx(1.0d100, 1.0d100, prec);
complex(prec), parameter :: c0 = cmplx(0.0d0, 0.0d0, prec), c1 = cmplx(1.0d0, 0.0d0, prec);

complex(prec), dimension(Nmax) :: a, b;

complex(prec) :: term, F, Fp;

integer :: k, j, M, N;

logical :: goon;

M = 0;
N = 4;

F  = NaN;
Fp = NaN;

goon = .true.;

do while (N <= Nmax .and. goon)

    term = c0; ! tail value

    do k = N,Nmin,-1

        if (k > M) then

            j = k / 2;

            a(k) = cmplx(j, 0, prec);

            if (j * 2 .eq. k) a(k) = a(k) - s;

            b(k) = c1;

        end if

        term = a(k) / z / (b(k) + term);

        ! print *, k, j, a(k), b(k), term

    end do

    Fp = F;

    F = c1 / z / (c1 + term); ! top level term

    goon = ( abs(F - Fp) > abs_err + rel_err * abs(F) );

    if (N > M) M = N;

    N = N * 2;

end do

est_err = abs(c1 - Fp/F);

if (goon) then
    status = 1;
else
    status = 0;
end if

Nterms = N / 2;

! return with or without exponent:
if (scal == 0) then

    G = z**s * exp(-z) * F;

else

    G = z**s * F;

end if

end subroutine

!-----------------------------------------------------------------------------!
