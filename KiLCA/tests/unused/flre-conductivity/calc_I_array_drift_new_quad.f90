!*******************************************************************!

module spec_data

use constants, only: dpc, dp;

implicit none;

! for subintegral function:
integer :: mm, nn, ll;

complex(dpc) :: bx1, bx2, bx4, gamma, gr;

! for quadrature:
integer, parameter :: limit = 1024, lenw = 4*limit;

integer, dimension(limit) :: iwork;
real(dp), dimension(lenw) :: work;

integer :: part;

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

end module

!*******************************************************************!

subroutine init_spec_data(npara, nperp)

use spec_data, only: mpara, mperp;
use spec_data, only: MC, dblfac, binocoef, pow2, powmibx1;

implicit none;

integer :: npara, nperp;

integer :: i, j;

integer :: m, n, j1, j2, k1, k2, l1, l2;

mpara = npara;
mperp = nperp;

allocate(MC(0:mpara, 0:mpara, 0:mpara, 0:mpara));

allocate(dblfac(-1:2*mpara), binocoef(0:mpara, 0:mpara), pow2(0:mpara), powmibx1(0:2*mpara));

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

do i = 1,mpara

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

            end do
        end do
    end do
end do

powmibx1 = 1.0d100;

end subroutine

!*******************************************************************!

subroutine close_spec_data()

use spec_data, only: MC, dblfac, binocoef, pow2, powmibx1;

implicit none;

deallocate(MC, dblfac, binocoef, pow2, powmibx1);

end subroutine

!*******************************************************************!

subroutine exparg(t, arg)

use constants, only: dpc, dp, im;
use spec_data, only: bx1, bx2, bx4, gamma, gr;
use spec_data, only: MC, dblfac, pow2, powmibx1;
use spec_data, only: mm, nn, ll;
use spec_data, only: dtaudt;

implicit none;

real(dp) :: t;
complex(dpc) :: arg;

complex(dpc), parameter :: one = cmplx(1.0d0, 0.0d0, dpc);

complex(dpc) :: tau, expon, onepluexp, oneminexp, denplu, denmin, fac1, fac2, fac3, facm, ans, summa;

integer :: j1, j2, jmax, ind;

tau = dtaudt * t;

expon = exp(-tau);

onepluexp = one + expon;
oneminexp = one - expon;

denplu = one - gamma * onepluexp;
denmin = one - gamma * oneminexp;

fac1 = onepluexp / denplu;
fac2 = oneminexp / denmin;
fac3 = oneminexp / denplu;

ans = - (ll+1) * log(one + bx4 * tau * im) - 0.5d0 * log(denmin * denplu) + &

        gr * tau + fac3 * (one - 2.0d0 * gamma) * bx1 * bx1;

jmax = (mm + nn) / 2;

summa = cmplx(0.0d0, 0.0d0, dpc);

do j1 = 0,jmax

    do j2 = 0,jmax

        ind = mm + nn - 2 * (j1 + j2);

        if (ind < 0) cycle;

        facm = MC(mm,nn,j1,j2) * dblfac(2*j1 - 1) * dblfac(2*j2 - 1) / pow2(j1 + j2) * powmibx1(ind);

        summa = summa + facm * fac1**j1 * fac2**j2 * fac3**ind;

    end do

end do

arg = ans + log(summa);

end subroutine

!*******************************************************************!

function subint(t)

use constants, only: dpc, dp, im;
use spec_data, only: bx1, bx2, bx4, gamma, gr;
use spec_data, only: MC, dblfac, pow2, powmibx1;
use spec_data, only: mm, nn, ll;
use spec_data, only: dtaudt;
use spec_data, only: part;

implicit none;

real(dp) :: t;
real(dp) :: subint;

complex(dpc), parameter :: one = cmplx(1.0d0, 0.0d0, dpc);

complex(dpc) :: tau, expon, onepluexp, oneminexp, denplu, denmin, fac1, fac2, fac3, facm, ans, summa;

integer :: j1, j2, jmax, ind;

tau = dtaudt * t;

expon = exp(-tau);

onepluexp = one + expon;
oneminexp = one - expon;

denplu = one - gamma * onepluexp;
denmin = one - gamma * oneminexp;

fac1 = onepluexp / denplu;
fac2 = oneminexp / denmin;
fac3 = oneminexp / denplu;

ans = one / (one + bx4 * tau * im)**(ll+1) / sqrt(denmin * denplu) * &
      exp( gr * tau + fac3 * (one - 2.0d0 * gamma) * bx1 * bx1 );

jmax = (mm + nn) / 2;

summa = cmplx(0.0d0, 0.0d0, dpc);

do j1 = 0,jmax

    do j2 = 0,jmax

        ind = mm + nn - 2 * (j1 + j2);

        if (ind < 0) cycle;

        facm = MC(mm,nn,j1,j2) * dblfac(2*j1 - 1) * dblfac(2*j2 - 1) / pow2(j1 + j2) * powmibx1(ind);

        summa = summa + facm * fac1**j1 * fac2**j2 * fac3**ind;

    end do

end do

ans = ans * summa;

if (part == 0) then
    subint = real(dtaudt * ans);
else
    subint = imag(dtaudt * ans);
end if

end function

!*******************************************************************!

subroutine func(t, f)

use constants, only: dpc;
use spec_data, only: argmin;

implicit none;

integer, parameter :: prcsn = 8

real(prcsn), intent(in)  :: t;
real(prcsn), intent(out) :: f;

complex(dpc) :: arg;

call exparg(t, arg);

f = real(arg) - argmin;

end subroutine

!*******************************************************************!

subroutine jacob(t, j)

use constants, only: dpc;

implicit none;

integer, parameter :: prcsn = 8

real(prcsn), intent(in)  :: t;
real(prcsn), intent(out) :: j;

real(prcsn) :: dt = 1.0d-6;
complex(dpc) :: fm, fp;

call exparg(t - dt, fm);
call exparg(t + dt, fp);

j = real(fp - fm) / 2.0d0 / dt;

end subroutine

!*******************************************************************!

subroutine find_opt_contour()

use constants, only: dpc, dp, im, pi;
use spec_data, only: bx1, bx2, gr, phi, dtaudt;

implicit none;

real(dp) ::  dargdt, grt, grmin;

real(dp) :: phit, dphi;

integer :: Nc = 17, k;

dphi = pi / (Nc-1);

phit = - atan(imag(gr) / real(gr));

grmin = 0.0d0;
phi   = 0.0d0;

do k = 1,Nc+1

    dtaudt = cos(phit) + sin(phit)*im;

    grt = real(gr * dtaudt);

    if(grt <= 0.0d0) then

        call jacob(0.0d0, dargdt);

        if(dargdt <= 0.0d0 .and. grt < grmin) then

            grmin = grt;
            phi   = phit;

        end if

    end if

!    print *, 'contour:', phit, grt, dargdt;

    phit = - pi / 2.0d0 + (k-1)*dphi;

end do

dtaudt = cos(phi) + sin(phi)*im;

if(grmin == 0.0d0) print *, 'warning: find_opt_contour: failed to find a good contour:', gr, phi;

end subroutine

!*******************************************************************!

subroutine calc_imnl_quad(npara, nperp, x1, x2, x3, x4, Imnl, status)

use constants, only: dpc, dp, im, sqrt2p;
use spec_data, only: bx1, bx2, bx4, gamma, gr, part;
use spec_data, only: mpara, powmibx1;
use spec_data, only: mm, nn, ll;
use spec_data, only: limit, lenw, iwork, work;
use spec_data, only: phi, dtaudt, argmin;

implicit none;

external subint, func, jacob;

integer :: npara, nperp;
complex(dpc) :: x1, x2, x3, x4;
complex(dpc), dimension(0:npara,0:npara,0:nperp) :: Imnl;
integer :: status;

! for quadrature:
real(dp) :: epsabs = 1.0d-10, epsrel = 1.0d-10;
real(dp) :: abserr;
real(dp) :: a, b;
real(dp) :: int_re, int_im;
integer  :: key = 20, neval, ier, last;

! for the Newton:
real(dp)                     :: zmin = 0.0d0, zmax = 1.0d6;
integer, parameter           :: n_start = 5;
real(dp), dimension(n_start) :: z_start;
real(dp)                     :: abs_eps_z = 1.0d-3, rel_eps_z = 1.0d-3;
real(dp)                     :: abs_eps_f = 1.0d-3, rel_eps_f = 1.0d-3;
integer                      :: max_iter_num = 32;
real(dp)                     :: rest;

complex(dpc) :: alpha, ep2a, factor;

! for debugging:
integer :: i, j, k;
real(dp) :: t;
complex(dpc) :: z, arg;

status = 0;

if(npara > mpara) then

    print *, 'warning: the values of npara and mpara do not match:', npara, mpara;
    status = 1;
    return;

end if

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

a = 0.0d0;

do mm = 0,npara
    do nn = 0,mm

!print *, '(m,n) =', mm, nn;

    factor = - sqrt2p / ep2a**(0.5d0*(3 + mm + nn));

    ll = 0;
    call find_opt_contour();

!print *, 'phi =', phi;

    ! for b search:
    z_start(1) = 1.0d2  * argmin / real(gr * dtaudt);
    z_start(2) = 1.0d-1 * z_start(1);
    z_start(3) = 1.0d-1 * z_start(2);
    z_start(4) = 1.0d-1 * z_start(3);
    z_start(5) = 1.0d-1 * z_start(4);

    do ll = 0,nperp

!print *, 'll =', ll;

        call solve_system_newton(func, jacob, zmin, zmax, n_start, z_start, max_iter_num, &
                                 abs_eps_z, rel_eps_z, abs_eps_f, rel_eps_f, b, rest, ier);

        if(ier /= 0) then
            status = 1;
            print *, 'warning: solve_system_newton failed to find taumax:', phi, b, rest;
        end if

        z_start(1) = b;

!print *, 'b =', b;

        part = 0;

        call dqag(subint, a, b, epsabs, epsrel, key, int_re, abserr, neval, ier, limit, lenw, last, iwork, work);

        if(ier /= 0) then
            status = 1
            print *, 'warning: dqag failed to evaluate the integral:', int_re, abserr, neval, ier;
        end if

        part = 1;

        call dqag(subint, a, b, epsabs, epsrel, key, int_im, abserr, neval, ier, limit, lenw, last, iwork, work);

        if(ier /= 0) then
            status = 1
            print *, 'warning: dqag failed to evaluate the integral:', int_im, abserr, neval, ier;
        end if

        Imnl(mm,nn,ll) = factor * (int_re + int_im * im);
        Imnl(nn,mm,ll) = Imnl(mm,nn,ll);

!         open(100);
!         do i=1,10000
!             t = (i-1)*b/9999;
!             call exparg(t, arg);
!             write(100,*) t, real(arg), imag(arg);
!         end do
!         close(100);
! 
!         open(200);
!         do i=1,100
!             do j=1,100
! 
!                 z = -0.0d0 + (i-1)*20.0d0/99 + (-10.0d0 + (j-1)*20.0d0/99)*im;
!                 call subintcmplx(z, arg);
!                 write(200,*) real(z), imag(z), real(arg), imag(arg);
! 
!             end do
!         end do
!         close(200);

    end do

end do
end do

end subroutine

!*******************************************************************!

subroutine subintcmplx(tau, val)

use constants, only: dpc, dp, im;
use spec_data, only: bx1, bx2, bx4, gamma, gr;
use spec_data, only: MC, dblfac, pow2, powmibx1;
use spec_data, only: mm, nn, ll;

implicit none;

complex(dpc) :: tau, val;

complex(dpc), parameter :: one = cmplx(1.0d0, 0.0d0, dpc);

complex(dpc) :: expon, onepluexp, oneminexp, denplu, denmin, fac1, fac2, fac3, facm, ans, summa;

integer :: j1, j2, jmax, ind;

expon = exp(-tau);

onepluexp = one + expon;
oneminexp = one - expon;

denplu = one - gamma * onepluexp;
denmin = one - gamma * oneminexp;

fac1 = onepluexp / denplu;
fac2 = oneminexp / denmin;
fac3 = oneminexp / denplu;

ans = one / (one + bx4 * tau * im)**(ll+1) / sqrt(denmin * denplu) * &
      exp( gr * tau + fac3 * (one - 2.0d0 * gamma) * bx1 * bx1 );

jmax = (mm + nn) / 2;

summa = cmplx(0.0d0, 0.0d0, dpc);

do j1 = 0,jmax

    do j2 = 0,jmax

        ind = mm + nn - 2 * (j1 + j2);

        if (ind < 0) cycle;

        facm = MC(mm,nn,j1,j2) * dblfac(2*j1 - 1) * dblfac(2*j2 - 1) / pow2(j1 + j2) * powmibx1(ind);

        summa = summa + facm * fac1**j1 * fac2**j2 * fac3**ind;

    end do

end do

val = ans * summa;

end subroutine

!*******************************************************************!
