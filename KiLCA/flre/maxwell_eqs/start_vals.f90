!<A collection of subroutines used to calculate starting values for Maxwell system (ODE) solution.

!------------------------------------------------------------------------------

subroutine solution_in_center (r, kper, pvec, zvec)

use constants, only: dp, dpc, im
use flre_sett, only: Nwaves
use mode_data, only: m
use maxwell_equations, only: iErsp_state, iErsp_sys, dim_Ersp_state, num_vars, num_eqs

implicit none

interface
    recursive double complex function besselj (nu, zarg, n) result (res)
        double complex, intent(in) :: zarg
        integer, intent(in) :: nu, n
    end function
end interface

real(dp),                        intent(in)  :: r
complex(dpc),                    intent(in)  :: kper !perpendicular wave number
complex(dpc), dimension(Nwaves), intent(in)  :: pvec !Cartesian polarisation vector
complex(dpc), dimension(Nwaves), intent(out) :: zvec !starting vector

complex(dpc), allocatable, dimension(:) :: v_sys
complex(dpc), allocatable, dimension(:) :: rhs

complex(dpc) :: vp, vm

integer :: k, i, nu

nu = m

allocate (v_sys(num_vars), rhs(num_eqs))

!transforms cartesian polvec to cylinder coordinates and compute the starting values

!bessel derivs of any order should be used or transformation of derivs:
!near r=0: s~~theta,p~~z

v_sys = cmplx(0.0d0, 0.0d0, 8);
rhs   = cmplx(0.0d0, 0.0d0, 8);

call state2sys (r, 'w', pvec, v_sys, rhs) !WKB case is here!

vp = v_sys(iErsp_sys(1)) + im*v_sys(iErsp_sys(2))
vm = v_sys(iErsp_sys(1)) - im*v_sys(iErsp_sys(2))

!Er derivatives:
k = 1
do i=0,dim_Ersp_state(k)-1
    zvec(iErsp_state(k)+i) = -0.5d0*im*(kper**i)* &
    (vm*besselj(nu-1,kper*r,i) - vp*besselj(nu+1,kper*r,i))
end do

!Es derivatives:
k = 2
do i=0,dim_Ersp_state(k)-1
    zvec(iErsp_state(k)+i) = 0.5d0*(kper**i)* &
    (vm*besselj(nu-1,kper*r,i) + vp*besselj(nu+1,kper*r,i))
end do

!Ep derivatives:
k = 3
do i=0,dim_Ersp_state(k)-1
    zvec(iErsp_state(k)+i) = pvec(iErsp_state(k))*(kper**i)*besselj (nu,kper*r,i)
end do

deallocate (v_sys, rhs);

end subroutine

!------------------------------------------------------------------------------

subroutine calc_start_values_center_with_correct_asymptotic (rstart, zstart)

!This subroutine provides the routines for calculating the starting values for
!integration.
!rstart - radial value where to calculate start values; determine
!whether value lies inside the antenna or outside
!zstart - matrix containing start vectors
!dimension 1: number of fundamental solution to be integrated
!dimension 2: number of components per solution

use constants, only: dp, dpc, c, im
use flre_sett, only: flre_order, Nwaves
use maxwell_equations, only: iErsp_state

implicit none

real(dp), intent(in) :: rstart
complex(dpc), dimension(Nwaves,Nwaves/2), intent(out) :: zstart

integer :: i

complex(dpc), dimension(Nwaves) :: kperp
complex(dpc), dimension(Nwaves,Nwaves) :: polv

call calc_dispersion (rstart, 'w', 0, kperp, polv) !check if k's are sorted!

! do i=1,Nwaves/2
!     call solution_in_center (rstart, kperp(2*i-1), polv(:,2*i-1), zstart(:,i))
! end do

do i=1,Nwaves/2
    call solution_in_center (rstart, kperp(i), polv(:,i), zstart(:,i))
end do

! print *, 'calc_start_values: starting values at r =', rstart
! do i=1,Nwaves/2
!     print *, 'z(',i,') =', zstart(:,i)
! end do

end subroutine

!------------------------------------------------------------------------------

subroutine calc_start_values_ideal_wall_low_derivs (rstart, zstart)

!This subroutine provides the routines for calculating the starting values for
!integration.
!rstart - radial value where to calculate start values; determine
!whether value lies inside the antenna or outside
!zstart - matrix containing start vectors
!dimension 1: number of fundamental solution to be integrated
!dimension 2: number of components per solution

use constants, only: dp, dpc, c, im
use flre_sett, only: flre_order, Nwaves
use maxwell_equations, only: iErsp_state

implicit none

real(dp), intent(in) :: rstart
complex(dpc), dimension(Nwaves,Nwaves/2), intent(out) :: zstart

integer, dimension(3) :: imin, imax
integer :: i, j, ind

imin(1) = iErsp_state(1)
imax(1) = iErsp_state(1)+flre_order-2

imin(2) = iErsp_state(2)+1
imax(2) = iErsp_state(2)+flre_order

imin(3) = iErsp_state(3)+1
imax(3) = iErsp_state(3)+flre_order

!print *, 'calc_start_values: imin=', imin
!print *, 'calc_start_values: imax=', imax

ind = 0
do i=1,3
    do j=imin(i),imax(i)
        ind = ind+1
        zstart(:,ind) = cmplx(0.d0,0.d0,dp)
        zstart(j,ind) = cmplx(1.d0,0.d0,dp)
    end do
end do

if (ind /= Nwaves/2) then
    print *, 'error: calc_start_values: calc_start: dim zstart vector =', ind
end if

!print *, 'calc_start_values: starting values at r =', rstart
!do i=1,Nwaves/2
!    print *, 'z(',i,') =', zstart(:,i)
!end do

end subroutine

!------------------------------------------------------------------------------

subroutine calc_start_values_anywhere_low_derivs (rstart, zstart)

!This subroutine provides the routines for calculating the starting values for
!integration.
!rstart - radial value where to calculate start values; determine
!whether value lies inside the antenna or outside
!zstart - matrix containing start vectors
!dimension 1: number of fundamental solution to be integrated
!dimension 2: number of components per solution

use constants, only: dp, dpc, c, im
use flre_sett, only: flre_order, Nwaves
use maxwell_equations, only: iErsp_state

implicit none

real(dp), intent(in) :: rstart
complex(dpc), dimension(Nwaves,Nwaves/2), intent(out) :: zstart

integer, dimension(3) :: imin, imax
integer :: i, j, ind

imin(1) = iErsp_state(1)
imax(1) = iErsp_state(1)+flre_order-2

imin(2) = iErsp_state(2)+1
imax(2) = iErsp_state(2)+flre_order

imin(3) = iErsp_state(3)+1
imax(3) = iErsp_state(3)+flre_order

!print *, 'calc_start_values: imin=', imin
!print *, 'calc_start_values: imax=', imax

ind = 0

do i=1,3
    do j=imin(i),imax(i)
        ind = ind+1
        zstart(:,ind) = cmplx(0.d0,0.d0,dp)
        zstart(j,ind) = cmplx(1.d0,0.d0,dp)
    end do
end do

if (ind /= Nwaves/2) then
    print *, 'error: calc_start_values: calc_start: dim zstart vector =', ind
end if

!print *, 'calc_start_values: starting values at r =', rstart
!do i=1,Nwaves/2
!    print *, 'z(',i,') =', zstart(:,i)
!end do

end subroutine

!------------------------------------------------------------------------------

!implement starting vectors with supressed fake modes...
!problem: starting vectors must be the same in all moving frames, so we probably need
!conductivity and rhs matrix in lab frame to compute starting vectors without fake modes
!the same is true for starting  vectors with correct asymptotics!!!
!there must be a flag to evaluate conductivity in a lab frame...
