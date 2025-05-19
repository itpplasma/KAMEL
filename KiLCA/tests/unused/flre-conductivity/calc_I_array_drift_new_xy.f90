!*******************************************************************!

module quad_data

use constants, only: dpc, dp;

implicit none;

! for subintegral function:
integer :: ll;

complex(dpc) :: bx1, bx2, bx4, gamma;

real(dp) :: x, y;

integer :: part;

! for quadrature:
integer, parameter :: limit = 1024, lenw = 4*limit;

integer, dimension(limit) :: iwork;
real(dp), dimension(lenw) :: work;

! for contour:
real(dp) :: phi;
complex(dpc) :: dtaudt;

! for mesh:
integer, parameter :: N = 3;
real(dp), parameter :: delta = 1.0d-6;
real(dp), dimension(-N:N) :: node = (/-3.0d0*delta, -2.0d0*delta, -1.0d0*delta, 0.0d0, 1.0d0*delta, 2.0d0*delta, 3.0d0*delta/);
complex(dpc), dimension(-N:N,-N:N) :: Qij;

end module

!*******************************************************************!

subroutine exparg(t, arg)

use constants, only: dpc, dp, im;
use quad_data, only: bx1, bx2, bx4, gamma, ll, x, y, part;
use quad_data, only: dtaudt;

implicit none;

real(dp) :: t;
real(dp) :: arg;

complex(dpc) :: tau, gme, ex, gex, gex2, den, r2, ans;

tau = dtaudt * t;

gme = gamma - 1.0d0;

ex = exp(-tau);

gex = gamma * ex;

gex2 = gex * gex;

den = gme * gme - gex2;

r2 = x*x + y*y;

ans = - (ll+1) * log(1.0d0 + bx4*tau*im) - 0.5d0 * log(den) + &

        (-bx1*bx1 + bx2*im) * tau + &

        (1.0d0 - ex) / (gme + gex) * bx1 * ((2.0d0*gamma - 1.0d0)*bx1 + (x + y)*im) + &

        (r2 * (gex*ex - gme) + 2.0d0*x*y*ex) / 2.0d0 / den;

if (part == 0) then
    arg = real(ans);
else
    arg = imag(ans);
end if

end subroutine

!*******************************************************************!

function subint(t)

use constants, only: dpc, dp, sqrt2p, im;
use quad_data, only: bx1, bx2, bx4, gamma, ll, x, y, part;
use quad_data, only: dtaudt;

implicit none;

real(dp) :: t;
real(dp) :: subint;

complex(dpc) :: tau, gme, ex, gex, gex2, den, r2, ans;

tau = dtaudt * t;

gme = gamma - 1.0d0;

ex = exp(-tau);

gex = gamma * ex;

gex2 = gex * gex;

den = gme * gme - gex2;

r2 = x*x + y*y;

ans = - sqrt2p / (1.0d0 + bx4*tau*im)**(ll+1) / sqrt(den) * &

        exp( (-bx1*bx1 + bx2*im) * tau + &

             (1.0d0 - ex) / (gme + gex) * bx1 * ((2.0d0*gamma - 1.0d0)*bx1 + (x + y)*im) + &

             (r2 * (gex*ex - gme) + 2.0d0*x*y*ex) / 2.0d0 / den);

if (part == 0) then
    subint = real(dtaudt * ans);
else
    subint = imag(dtaudt * ans);
end if

end function

!*******************************************************************!

subroutine func(t, f)

use quad_data, only: part;

implicit none;

integer, parameter :: prcsn = 8

real(prcsn), intent(in)  :: t;
real(prcsn), intent(out) :: f;

part = 0;

call exparg(t, f);

f = f + 40.0d0;

end subroutine

!*******************************************************************!

subroutine jacob(t, j)

use quad_data, only: part;

implicit none;

integer, parameter :: prcsn = 8

real(prcsn), intent(in)  :: t;
real(prcsn), intent(out) :: j;

real(prcsn) :: fm, fp, dt = 1.0d-6;

part = 0;

call exparg(t - dt, fm);
call exparg(t + dt, fp);

j = (fp - fm) / 2.0d0 / dt;

end subroutine

!*******************************************************************!

subroutine find_opt_contour()

use constants, only: dpc, dp, im, pi;
use quad_data, only: bx1, bx2, phi, dtaudt;

implicit none;

complex(dpc) :: arg;

real(dp) ::  dargdt, gr, grmin;

real(dp) :: phit, dphi;

integer :: Nc = 17, k;

arg = - bx1 * bx1 + im * bx2;

dphi = pi / (Nc-1);

phit = - atan(imag(arg) / real(arg));

grmin = 0.0d0;
phi   = 0.0d0;

do k = 1,Nc+1

    dtaudt = cos(phit) + sin(phit)*im;

    gr = real(arg * dtaudt);

    if(gr <= 0.0d0) then

        call jacob(0.0d0, dargdt);

        if(dargdt <= 0.0d0 .and. gr < grmin) then

            grmin = gr;
            phi   = phit;

        end if

    end if

    phit = - pi / 2.0d0 + (k-1)*dphi;

end do

dtaudt = cos(phi) + sin(phi)*im;

if(grmin == 0.0d0) print *, 'warning: find_opt_contour: failed to find a good contour:', arg, phi;

end subroutine

!*******************************************************************!

subroutine calc_imnl_quad(m_max, n_max, l_max, x1, x2, x3, x4, Imnl, status)

use constants, only: dpc, dp, im;
use quad_data, only: bx1, bx2, bx4, gamma, ll, x, y, part;
use quad_data, only: limit, lenw, iwork, work;
use quad_data, only: phi, dtaudt;

implicit none;

external subint, func, jacob;

integer :: m_max, n_max, l_max;
complex(dpc) :: x1, x2, x3, x4;
complex(dpc), dimension(0:m_max,0:n_max,0:l_max) :: Imnl;
integer :: status;

complex(dpc) :: alpha, ep2a, factor, arg;

! for quadrature:
real(dp), parameter :: epsabs = 1.0d-12, epsrel = 1.0d-12;

real(dp) :: abserr;

real(dp) :: a, b;

real(dp) :: int_re, int_im;

integer :: key = 20, neval, ier, last;

! for the Newton:
real(dp)                     :: zmin = 0.0d0, zmax = 1.0d8;
integer, parameter           :: n_start = 3;
real(dp), dimension(n_start) :: z_start;
real(dp)                     :: abs_eps_z = 1.0d-3, rel_eps_z = 1.0d-3;
real(dp)                     :: abs_eps_f = 1.0d-3, rel_eps_f = 1.0d-3;
integer                      :: max_iter_num = 32;
real(dp)                     :: rest;

integer :: m, n;

status = 0;

alpha = sqrt(0.25d0 - x3*im) - 0.5d0;

ep2a = 1.0d0 + 2.0d0*alpha;

gamma = alpha / ep2a;

bx1 = x1 / ep2a**1.5d0;

bx2 = x2 / ep2a + gamma * im;

bx4 = x4 / ep2a;

arg = - bx1 * bx1 + im * bx2;

call find_opt_contour();

a = 0.0d0;

! for b:
z_start(1) = - 40.0d0 / real(arg * dtaudt);
z_start(2) = z_start(1);
z_start(3) = 1.0d-6;

do ll = 0,l_max

    call solve_system_newton(func, jacob, zmin, zmax, n_start, z_start, max_iter_num, &
                             abs_eps_z, rel_eps_z, abs_eps_f, rel_eps_f, b, rest, ier);

    if(ier /= 0) then
        status = 1;
        print *, 'warning: solve_system_newton failed to find taumax:', phi, b, rest;
    end if

    z_start(1) = b;

    ! mesh:
    do j = -N,N

        y = node(j);

        do i = j,N

            x = node(i);

            print *, x, y; ! check real(arg)

            part = 0;

            call dqag(subint, a, b, epsabs, epsrel, key, int_re, abserr, neval, ier, limit, lenw, last, iwork, work);

            if(ier /= 0) then
                status = 1
                print *, 'warning: dqag failed to evaluate the integral:', ier, neval;
            end if

            part = 1;

            call dqag(subint, a, b, epsabs, epsrel, key, int_im, abserr, neval, ier, limit, lenw, last, iwork, work);

            if(ier /= 0) then
                status = 1
                print *, 'warning: dqag failed to evaluate the integral:', ier, neval;
            end if

            Qij(i,j) = int_re + int_im * im;
            Qij(j,i) = Qij(i,j);

        end do

    end do

    call eval_neville_polynom_ (int * dim, const double * xg, const double * yg, int * deg,
                            double * x, int * Dmin, int * Dmax, int * ind, double * R);

    do m = 0,mmax
        do n = 0,nmax

            factor = 1.0d0 / ep2a**(0.5d0*(3 + m + n));

            call eval_func_grid(node, Qij, 0.0d0, 0.0d0, m, n, res);

            Imnl(m,n,ll) = factor * res;
            Imnl(n,m,ll) = Imnl(m,n,ll);

        end do
    end do

end do

end subroutine

!*******************************************************************!
