subroutine current_density (j_surf)

! This routine calculates the surface current density produced by the
! antenna for a specific mode. The mode numbers are imported from the
! module mode_data.
! j_surf ... vector of length 2 containing jtheta and jz (jphi)

use constants, only: dp, dpc, pi, im
use antenna_data, only: I0, ra, flag_debug
use background, only: rtor
use mode_data, only: m, n

implicit none;

complex(dpc), dimension(2), intent(out) :: j_surf

real(dp) :: F_s, arg, thc
complex(dpc) :: F_n

integer :: ml, nl

integer :: iftest = 1; !test spectrum: all harmonics have equal ampls (used for DIIID)

ml = m
nl = n

!if q>0 then for resonance n,m must have different signs
!Here we make thc<0 to fit the antenna to q>0 case

thc = -pi/3.0
arg = ml-pi*nl/thc

if (arg == 0) then
    F_s = 1.0d0/(2.0d0*pi*pi)*abs(thc)
else
    F_s = 1.0d0/(2.0d0*pi*pi)*sin(arg*abs(thc))/arg
end if

if (mod(nl-1,4)==0 .and. iftest == 0) then
    !configuration dependent:
    !(3,1) mode:
    if (flag_debug) then
        print *,'information: antenna is in the (3,1) configuration of DED!'
    end if
    F_n = 4.0*exp(-3.0*IM*pi*nl/16.0)*sin(pi*n/4.0)/sin(pi*n/16.0)
    j_surf(1) = (-1)**(ml+nl)*I0/rtor*F_s*F_n
    j_surf(2) = (-1)**(ml+nl)*I0/ra*pi/abs(thc)*F_s*F_n
else if (mod(nl-4,16)==0 .and. iftest == 0) then
    !configuration dependent:
    !(12,4) mode:
    if (flag_debug) then
        print *,'information: antenna is in the (12,4) configuration of DED!'
    end if
    F_n = 16.0
    j_surf(1) = (-1)**(ml+nl)*I0/rtor*F_s*F_n
    j_surf(2) = (-1)**(ml+nl)*I0/ra*pi/abs(thc)*F_s*F_n
else
    !use test spectrum:
    if (flag_debug) then
        print *,'information: antenna is in the test configuration!'
    end if
    F_s = 1.0d0/(2.0d0*pi*pi)*abs(thc)
    F_n = 1.0
    j_surf(1) = I0/rtor*F_s*F_n
    j_surf(2) = - (dble(ml)/dble(nl))*(rtor/ra)*j_surf(1) !to ensure zero div for any harm.
end if

if (flag_debug) then
    print *, 'antenna: n =', nl, 'm =', ml
    print *, 'antenna: j_th =', j_surf(1)
    print *, 'antenna: j_z =', j_surf(2)
    print *, 'antenna: check divergence: div =', ml*j_surf(1)/ra + nl*j_surf(2)/rtor
    print *
end if

end subroutine
