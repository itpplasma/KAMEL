!<The definitions of some constants used throughout the library.

module constants

implicit none

integer, parameter :: pp = 8 !integers size for C pointers, 32 bit - 4, 64 bit - 8 and so on.

integer, parameter :: sp = 4
integer, parameter :: dp = 8
integer, parameter :: qp = 16 !quad precision
integer, parameter :: vp = dp !

integer, parameter :: spc = 4
integer, parameter :: dpc = 8
integer, parameter :: qpc = 16 !quad precision
integer, parameter :: vpc = dpc !
integer, parameter :: lgt = kind (.true.)

real(dp), parameter :: pi    = 3.141592653589793238462643383279502884197_dp
real(dp), parameter :: pio2  = 1.57079632679489661923132169163975144209858_dp
real(dp), parameter :: twopi = 6.283185307179586476925286766559005768394_dp
real(dp), parameter :: sqrt2 = 1.41421356237309504880168872420969807856967_dp
real(dp), parameter :: euler = 0.5772156649015328606065120900824024310422_dp
real(dp), parameter :: boltz = 1.60216428d-12 ! erg/ev

real(dp), parameter :: srpi = sqrt(pi)
real(dp), parameter :: sr2 = sqrt(2.0d0)
real(dp), parameter :: sqrt2p = sqrt(2.0d0*pi)
real(dp), parameter :: isqrt2pi = 1.0d0/sqrt2p
real(dp), parameter :: tppoh = (2.0d0*pi)**(1.5d0)

complex(dpc), parameter :: im = (0.0d0, 1.0d0)

real(dp), parameter :: c  = 29979245800.0
real(dp), parameter :: mp = 1.67262158d-24
real(dp), parameter :: me = mp/1.8361526675d3
real(dp), parameter :: e  = 4.8032d-10

end module constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_pointer_precision (ppp)

use constants, only: pp;

integer, intent(out) :: ppp;

ppp = pp;

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!<The module is designed to store some FLRE zone properties used in Fortran part of the code.

module flre_sett

use constants, only: dp, dpc

implicit none;

!Conductivity settings:
integer :: flre_order = 5;  !order of FLR expansion
integer :: Nmax = 10;       !highest cyclotron harmonic
integer :: Nbmax = 8;       !max number of terms in Bessel expansions

end module

!--------------------------------------------------------------------

!<The module containes basic arrays for the evaluation of sigma (k) matrices.

module conduct_arrays

use constants, only: dpc, dp;

implicit none;

real(dp), allocatable, dimension(:) :: factor;

complex(dpc), allocatable, dimension(:,:,:,:) :: D; !m; n; b; l;

real(dp), allocatable, dimension(:) :: x;

end module

!------------------------------------------------------------------------------

subroutine allocate_and_set_conductivity_arrays ()

use flre_sett, only: Nmax, Nbmax, flre_order;
use conduct_arrays, only: D;
use conduct_arrays, only: factor;
use conduct_arrays, only: x
use constants, only: dp;

implicit none;


integer :: lc, b, Nf;

real(dp) :: NaN = 1.0d100;

!factorial:
Nf = max(flre_order, 2+Nmax+2*Nbmax);
allocate (factor(0:Nf));

factor(0) = 1;
do lc = 1,Nf
    factor(lc) = lc * factor(lc-1);
end do

allocate (D(0:1, 0:flre_order+1, -Nmax:Nmax, 0:Nbmax)); !m;n;lc;b;

D = NaN;

! allocate (Jlb(0:Nmax,0:Nbmax));
!
! ! only l>=0 is considered here since: Jlb(-l,b) = (-)^l Jlb(l,b)
! do lc = 0,Nmax
!
!     do b = 0,Nbmax
!
!         !print *, lc, b, factor(b), factor(lc+b);
!         !print *, (-1.0d0)**b
!         !print *,2.0d0**(lc+2*b)
!         Jlb(lc,b) = (-1.0d0)**b / factor(b) / factor(lc+b) / 2.0d0**(lc+2*b);
!         !print *, Jlb(lc,b);
!     end do
!
! end do

allocate (x(0:Nmax+2*Nbmax));

end subroutine

!------------------------------------------------------------------------------
