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

real(dp), parameter :: c  = 29979245800.0   !speed of light in vacuum
real(dp), parameter :: mp = 1.67262158d-24  !proton mass
real(dp), parameter :: me = mp/1.8361526675d3 !electron mass
real(dp), parameter :: e  = 4.8032d-10        !elementary charge

integer :: PInf, NaN

data Pinf /B'01111111100000000000000000000000'/      ! +Infinity
data NaN  /B'01111111100000100000000000000000'/      ! NaN

end module constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_pointer_precision (ppp)

use constants, only: pp;

integer, intent(out) :: ppp;

ppp = pp;

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
