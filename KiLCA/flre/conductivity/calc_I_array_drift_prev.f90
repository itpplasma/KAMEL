!*******************************************************************!

module spec_data

use constants, only: dpc, dp;

implicit none;

! for subintegral function:
integer :: npara, nperp;
complex(dpc) :: bx1, bx2, bx4, gamma, gr;
complex(dpc), allocatable, dimension(:) :: arrfac1, arrfac2, arrfac3, arrfacl;

! for ode solver:
integer :: MAXNEQ, NEQ;
complex(dpc), allocatable, dimension(:) :: III;

integer :: LRW, LIW;
integer, dimension(15) :: INFO;
real(dp), allocatable, dimension(:) :: RWORK;
integer, allocatable, dimension(:) :: IWORK;

! for contour:
real(dp) :: phi;
complex(dpc) :: dtaudt;
real(dp), parameter :: argmin = -40.0d0;

! for high moments:
integer                                   :: mpara, mperp;
real(dp), allocatable, dimension(:,:,:,:) :: MC;
real(dp), allocatable, dimension(:)       :: dblfac;
real(dp), allocatable, dimension(:,:)     :: binocoef;
real(dp), allocatable, dimension(:)       :: pow2;
complex(dpc), allocatable, dimension(:)   :: powmibx1;
complex(dpc), allocatable, dimension(:)   :: gamarr;

integer :: feval;

end module

!*******************************************************************!

subroutine init_spec_data(npar, nper)

use spec_data, only: mpara, mperp;
use spec_data, only: MC, dblfac, binocoef, pow2, powmibx1, gamarr;
use spec_data, only: MAXNEQ, III, LRW, LIW, INFO, RWORK, IWORK;
use spec_data, only: arrfac1, arrfac2, arrfac3, arrfacl;

implicit none;

integer :: npar, nper;

integer :: i, j;

integer :: m, n, j1, j2, k1, k2, l1, l2;

mpara = npar;
mperp = nper;

allocate(MC(0:mpara, 0:mpara, 0:mpara, 0:mpara));

allocate(dblfac(-1:2*mpara), binocoef(0:mpara, 0:mpara), pow2(0:2*mpara), powmibx1(0:2*mpara), gamarr(0:mperp));

allocate(arrfac1(0:mpara), arrfac2(0:mpara), arrfac3(0:2*mpara), arrfacl(1:mperp));

! double factorial:
dblfac(-1) = 1.0d0;
dblfac(0)  = 1.0d0;

do i = 1,2*mpara

    dblfac(i) = dblfac(i-2) * i;

end do

! binomial coefficients C^i_j:
binocoef = 0.0d0;
do i = 0,mpara

    binocoef(i,0) = 1.0d0;

    do j = 1,i

        binocoef(i,j) = binocoef(i,j-1) * dble(i-j+1) / dble(j);

    end do

end do

! powers of two:
pow2(0) = 1.0d0;

do i = 1,2*mpara

    pow2(i) = pow2(i-1) * 2.0d0;

end do

! M coefficients:
do m = 0,mpara
    do n = 0,mpara
        do j1 = 0,mpara
            do j2 = 0,mpara

                MC(m, n, j1, j2) = 0.0d0;

                do k1 = 0,mpara
                    do k2 = 0,mpara
                        do l1 = 0,k1
                            do l2 = 0,k2

                                if(2 * j1 == l1 + l2 .and. 2 * j2 == k1 + k2 - l1 - l2) then

                                MC(m, n, j1, j2) = MC(m, n, j1, j2) + binocoef(m,k1) * binocoef(n,k2) * &
                                                   binocoef(k1, l1) * binocoef(k2, l2) * (-1.0d0)**(k1-l1);

                                end if

                            end do
                        end do
                    end do
                end do

                MC(m, n, j1, j2) = MC(m, n, j1, j2) * dblfac(2*j1 - 1) * dblfac(2*j2 - 1) / pow2(j1 + j2);

            end do
        end do
    end do
end do

powmibx1 = 1.0d100;

do m = 0,mperp

    gamarr(m) = gamma(1.0d0 + dble(m)/2.0d0);

end do

! for ode solver:
MAXNEQ = (mperp) * (mpara + 1) * (mpara + 2); ! max number of real equations (l starts from 1)

allocate(III(0:MAXNEQ/2-1));

INFO(1) = 0;
INFO(2) = 0;
INFO(3) = 0;
INFO(4) = 0;

LRW = 130 + 21 * MAXNEQ;
LIW = 51;

allocate(RWORK(LRW), IWORK(LIW));

end subroutine

!*******************************************************************!

subroutine close_spec_data()

use spec_data, only: MC, dblfac, binocoef, pow2, powmibx1, gamarr;
use spec_data, only: arrfac1, arrfac2, arrfac3, arrfacl;
use spec_data, only: III, RWORK, IWORK;

implicit none;

deallocate(MC, dblfac, binocoef, pow2, powmibx1, gamarr);
deallocate(arrfac1, arrfac2, arrfac3, arrfacl);
deallocate(III, RWORK, IWORK);

end subroutine

!*******************************************************************!

subroutine argument(t, a)

use constants, only: dp, dpc;
use spec_data, only: III;

implicit none;

real(dp), intent(in)  :: t;
complex(dpc), intent(out) :: a;

integer :: ind;

call ode_func(t, III, III, 0.0d0, 0);

ind = maxloc(abs(III),1) - 1; ! array index starts from 0!

a = log(III(ind));

end subroutine

!*******************************************************************!

subroutine func(t, f)

use constants, only: dp;
use spec_data, only: argmin;
use spec_data, only: III;

implicit none;

real(dp), intent(in)  :: t;
real(dp), intent(out) :: f;

call ode_func(t, III, III, 0.0d0, 0);

f = log(maxval(abs(III))) - argmin;

end subroutine

!*******************************************************************!

subroutine jacob(t, j)

use constants, only: dp;
use spec_data, only: III;

implicit none;

real(dp), intent(in)  :: t;
real(dp), intent(out) :: j;

real(dp) :: dt = 1.0d-6;
real(dp) :: fm, fp;

call ode_func(t - dt, III, III, 0.0d0, 0); fm = log(maxval(abs(III)));
call ode_func(t + dt, III, III, 0.0d0, 0); fp = log(maxval(abs(III)));

j = (fp - fm) / 2.0d0 / dt;

end subroutine

!*******************************************************************!

subroutine find_opt_contour(tmax, status)

use constants, only: dp, im, pi;
use spec_data, only: gr, phi, dtaudt, argmin;
use spec_data, only: III;

implicit none;

external func, jacob;

real(dp) :: tmax;
integer  :: status;

integer, parameter :: Nc = 8, deg = 7, nmax = 16;

! for the Newton:
real(dp)                       :: zmin, zmax;
integer, parameter             :: n_start = 1;
real(dp), dimension(1:n_start) :: z_start;
real(dp)                       :: abs_eps_z = 1.0d-6, rel_eps_z = 1.0d-3;
real(dp)                       :: abs_eps_f = 1.0d-1, rel_eps_f = 0.0d-3;
integer                        :: ier, max_iter_num = 32;
real(dp)                       :: rest;

real(dp), dimension(-deg:nmax) :: taut, argt;

real(dp) :: der, tmax_est, arg0;

real(dp) :: phi0, phie, phim, phi1, phi2, phit, dphi;
real(dp) :: gr0, gre, grm, gro;

integer :: k, n, dflag;

status = 0;

phi0 =   atan(real(gr) / imag(gr));
phie = - atan(imag(gr) / real(gr));

der = - real(gr) * sin(phi0) - imag(gr) * cos(phi0);

gr0 = real(gr*(cos(phi0) + sin(phi0)*im));
gre = real(gr*(cos(phie) + sin(phie)*im));

if (gre < 0.0d0) then

    phim = phie;

else if (der < 0.0d0) then

    phim =   pi / 2.0d0;

else

    phim = - pi / 2.0d0;

end if

grm = real(gr*(cos(phim) + sin(phim)*im));

phi1 = phim;
phi2 = phi0;

dphi = (phi2 - phi1) / Nc;

phi = pi;

do k = 0,Nc-1

    phit = phi1 + dphi * k;

    dtaudt = cos(phit) + sin(phit)*im;

    call ode_func(0.0d0, III, III, 0.0d0, 0);

    arg0 = log(maxval(abs(III)));

    dflag = -1;

    n = -(deg + 1);

    tmax_est = 2.0d0**n * (argmin / real(gr * dtaudt));

    do while (n < nmax)

        n = n + 1;

        tmax_est = 2.0d0 * tmax_est;

        call ode_func(tmax_est, III, III, 0.0d0, 0);

        taut(n) = tmax_est;
        argt(n) = log(maxval(abs(III)));
        !Lahey:
        !argt(n) = maxval(abs(III));
        !if(argt(n) > 0.0d0) then
        !    argt(n) = log(argt(n));
        !else
        !    argt(n) = 0.0d0;
        !end if

        !print *, 'info: find_opt_contour: ', 'phit =', phit, 'grt =', real(gr*dtaudt), 'taut =', taut(n), 'argt =', argt(n);

        if (argt(n) > 0.0d0 .and. argt(n) > min(1.2d0 * arg0, arg0 + 1.0d0)) then

            dflag = 1;
            exit;

        end if

        if (argt(n) < argmin .and. n >= deg) exit;

    end do

    if (dflag < 0) then

        phi = phit;
        exit;

    end if

end do

if (phi == pi) then

    status = 1;
    print *, 'error: find_opt_contour: unable to match all requirements for the contour!';
    goto 100;

end if

dtaudt = cos(phi) + sin(phi)*im;

gro = real(gr * dtaudt);

! find tmax:
z_start(1) = 0.0d0;

do k = n-1,-deg,-1

    if ( (argt(k+1) - argmin) < 0.0d0 .and. (argt(k) - argmin) > 0.0d0 ) then

        z_start(1) = (taut(k) + taut(k+1)) / 2.0d0;
        zmin = taut(k);
        zmax = taut(k+1);
        exit;

    end if

end do

if (z_start(1) == 0.0d0) then

    status = 1;
    print *, 'error: find_opt_contour: unable to locate taumax!';
    goto 100;

end if

call solve_system_newton(func, jacob, zmin, zmax, n_start, z_start, max_iter_num, &
                         abs_eps_z, rel_eps_z, abs_eps_f, rel_eps_f, tmax, rest, ier);

if(ier /= 0) then

    call solve_system_bisection(func, zmin, zmax, max_iter_num, abs_eps_z, rel_eps_z, abs_eps_f, rel_eps_f, tmax, rest, ier);

end if

if(ier /= 0) then

    call func(tmax, rest);

    status = 1;

    print *, 'error: find_opt_contour: failed to find taumax!';

    goto 100;

end if

return;

100 continue;

print *, 'find_opt_contour summary:';
print *, 'phi0 =', phi0, 'gr0 =', gr0;
print *, 'phie =', phie, 'gre =', gre;
print *, 'phim =', phim, 'grm =', grm;
print *, 'phio =', phi,  'gro =', gro;
print *, 'taut =', taut;
print *, 'argt =', argt;
print *, 'z_start =', z_start;
print *, 'tmax =', tmax, 'rest =', rest;

end subroutine

!*******************************************************************!

subroutine ode_func(t, u, uprime, rpar, ipar)

use constants, only: dpc, dp, im;
use spec_data, only: npara, nperp;
use spec_data, only: bx1, bx4, gamma, gr;
use spec_data, only: MC, powmibx1;
use spec_data, only: arrfac1, arrfac2, arrfac3, arrfacl;
use spec_data, only: dtaudt;
use spec_data, only: NEQ;
use spec_data, only: feval;

implicit none;

real(dp) :: t;
complex(dpc), dimension(0:NEQ/2-1) :: u, uprime;
real(dp) :: rpar;
integer :: ipar;

complex(dpc), parameter :: one = cmplx(1.0d0, 0.0d0, dpc);

complex(dpc) :: tau, expon, onepluexp, oneminexp, denplu, denmin, fac1, fac2, fac3, facl, comm, summa;

integer :: j1, j2, jmax, ind;

integer :: mm, nn, ll;

feval = feval + 1;

tau = dtaudt * t;

expon = exp(-tau);

onepluexp = one + expon;
oneminexp = one - expon;

denplu = one - gamma * onepluexp;
denmin = one - gamma * oneminexp;

fac1 = onepluexp / denplu;
fac2 = oneminexp / denmin;
fac3 = oneminexp / denplu;
facl = one + bx4 * tau * im;

! precompute arrays:
do mm = 0,npara

    arrfac1(mm) = fac1**mm;
    arrfac2(mm) = fac2**mm;

end do

do mm = 0,2*npara

    arrfac3(mm) = fac3**mm;

end do

do ll = 1,nperp

    arrfacl(ll) = facl**(dble(ll+1)/2.0d0);

end do

comm = dtaudt / sqrt(denmin * denplu) * exp( gr * tau + fac3 * (one - 2.0d0 * gamma) * bx1 * bx1 );

do mm = 0,npara

    do nn = 0,mm

        jmax = (mm + nn) / 2;

        summa = cmplx(0.0d0, 0.0d0, dpc);

        do j1 = 0,jmax

            do j2 = 0,jmax

                ind = mm + nn - 2 * (j1 + j2);

                if (ind < 0) cycle;

                summa = summa + MC(mm,nn,j1,j2) * powmibx1(ind) * arrfac1(j1) * arrfac2(j2) * arrfac3(ind);

            end do

        end do

        do ll = 1,nperp

            ind = (mm - nn) + (nn * (2 * npara + 3 - nn)) / 2 + (ll-1) * (((npara + 1) * (npara + 2)) / 2);

            uprime(ind) = comm * summa / arrfacl(ll);

        end do

    end do

end do

end subroutine

!*******************************************************************!

subroutine calc_imnl_quad(npar, nper, x1, x2, x3, x4, Imnl, status)

use constants, only: dpc, dp, im, sqrt2p;
use spec_data, only: bx1, bx2, bx4, gamma, gr, powmibx1;
use spec_data, only: mpara, mperp, npara, nperp;
use spec_data, only: III, NEQ, LRW, LIW, INFO, RWORK, IWORK;
use spec_data, only: feval;

implicit none;

external ode_func;

integer :: npar, nper; ! maximum parallel and perpendicular (nper = lperp = 2*mperp+1) indices
complex(dpc) :: x1, x2, x3, x4; ! arguments (without scaling)
complex(dpc), dimension(0:npar,0:npar,1:nper) :: Imnl; ! the array of I-function values
integer :: status; ! evaluation status - 0 means OK

complex(dpc) :: alpha, ep2a, factor, arg;

! for debugging:
!integer :: i, j, k;
!real(dp) :: tt;
!complex(dpc) :: zz, ff;

! for ode solver:
real(dp) :: T, TOUT;
real(dp) :: RTOL, ATOL;
real(dp) :: RPAR;
integer  :: IPAR;
integer  :: IDID;

integer :: k, mm, nn, ll, ind;

status = 0;

feval = 0;

Imnl = cmplx(0.0d0, 0.0d0, dp);

if(npar > mpara .or. nper > mperp) then

    print *, 'error: calc_imnl_quad: the values of (npar, nper) do not match to (mpara, mperp):', npar, nper, mpara, mperp;
    status = 1;
    return;

end if

! set for subintegral function:
npara = npar;
nperp = nper;

alpha = sqrt(0.25d0 - x3*im) - 0.5d0;

ep2a = 1.0d0 + 2.0d0*alpha;

gamma = alpha / ep2a;

bx1 = x1 / ep2a**1.5d0;

bx2 = x2 / ep2a + gamma * im;

bx4 = x4 / ep2a;

gr = - bx1 * bx1 + im * bx2;

! set array of powers:
do k = 0,2*mpara

    powmibx1(k) = (-im * bx1)**k;

end do

! solve all quadratures by ode integration:
NEQ = (nperp) * (npara + 1) * (npara + 2); ! number of real equations (l starts from 1)

call find_opt_contour(TOUT, status);

if (status /= 0) then

    print *, 'error: calc_imnl_quad: failed to find the integration contour!';
    return;

end if

T = 0.0d0;

INFO(1) = 0;

RTOL = 1.0d-10;
ATOL = 1.0d-10;

IDID = 0;

do while(IDID < 2)

    call ddeabm(ODE_FUNC, NEQ, T, III, TOUT, INFO, RTOL, ATOL, IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR);

    !call dderkf(ODE_FUNC, NEQ, T, III, TOUT, INFO, RTOL, ATOL, IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR);

    INFO(1) = 1;

    if (IDID < -1) then
        print *, 'warning: calc_imnl_quad: solver has troubles to integrate the equations:', TOUT, T, IDID;
    end if

end do

do mm = 0,npara

    do nn = 0,mm

        factor = cmplx(1.0d0, 0.0d0, dpc) / ep2a**(0.5d0*(3 + mm + nn));

        do ll = 1,nperp

            ind = (mm - nn) + (nn * (2 * npara + 3 - nn)) / 2 + (ll-1) * (((npara + 1) * (npara + 2)) / 2);

            Imnl(mm,nn,ll) = factor * III(ind);

            Imnl(nn,mm,ll) = Imnl(mm,nn,ll);

        end do

    end do

end do

if (feval > 50000) then

    print *, 'warning: calc_imnl_quad: large feval =', feval;
    print *, 'taumax =', TOUT;
    print *, 'x =', x1, x2, x3, x4;

end if

!debugging:
! TOUT = 1.0d3;
! open(100);
! do k=1,10000
!     t = (k-1)*TOUT/9999;
!     call argument(t, arg);
!     call ode_func(t, III, III, 0.0d0, 0);
!     write(100,*) t, real(arg), imag(arg), real(III), imag(III);
! end do
! close(100);

end subroutine

!*******************************************************************!

subroutine calc_imnkl_quad(econs, npar, nper, x1, x2, x3, x4, Imnkl, status)

use constants, only: dpc, dp;
use spec_data, only: gamarr;

implicit none;

logical :: econs; ! flag specifies whether to use the energy conservation formulas
integer :: npar, nper; ! maximum parallel and perpendicular indices
complex(dpc) :: x1, x2, x3, x4; ! arguments (without scaling)
complex(dpc), dimension(0:npar,0:npar,0:nper,0:nper) :: Imnkl; ! the array of I-function values
integer :: status; ! evaluation status - 0 means OK

complex(dpc), allocatable, dimension(:,:,:) :: Imnl;

integer :: mm, nn, kk, ll, lper;

complex(dpc) :: den;

lper = 2 * nper + 1;

allocate(Imnl(0:npar,0:npar,1:lper));

call calc_imnl_quad(npar, lper, x1, x2, x3, x4, Imnl, status);

if (status /= 0) then

    print *, 'error: calc_imnkl_quad: failed to evaluate Imnkl functions:', econs, npar, nper, x1, x2, x3, x4;
    deallocate(Imnl);
    return;

end if

if (econs) then

    den = cmplx(1.0d0, 0.0d0, dpc) + ( 2.0d0 * Imnl(2,0,1) - Imnl(2,2,1) - Imnl(0,0,1) );

    do mm = 0,npar
        do nn = 0,npar
            do ll = 0,nper
                do kk = 0,nper

                    Imnkl(mm,nn,kk,ll) = Imnl(mm,nn,kk+ll+1) + (gamarr(kk)*gamarr(ll)/gamarr(kk+ll)) * &

                    ((Imnl(mm,2,kk+1) - Imnl(mm,0,kk+1)) * (Imnl(nn,2,ll+1) - Imnl(nn,0,ll+1)) / den);

                end do
            end do
        end do
    end do

else

    do mm = 0,npar
        do nn = 0,npar
            do ll = 0,nper
                do kk = 0,nper

                    Imnkl(mm,nn,kk,ll) = Imnl(mm,nn,kk+ll+1);

                end do
            end do
        end do
    end do

end if

deallocate(Imnl);

end subroutine

!*******************************************************************!
