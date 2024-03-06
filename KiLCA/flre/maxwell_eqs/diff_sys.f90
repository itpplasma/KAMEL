!<Subroutines used to obtain ODE system matrix A in the direct form: u' = A*u.

!------------------------------------------------------------------------------
!Subs below do not depend (hopefully :)) on the number of equations:
!------------------------------------------------------------------------------

subroutine state2sys (r, flagback, v_state, v_sys, rhs)

!flagback: flag for background
!r: radial value; calculate background and derivatives at this value
!v_state: vector of state;
!v_sys: state at radial value r

use flre_sett, only: Nwaves, flag_debug
use maxwell_equations

implicit none

real(dp), intent(in) :: r
character(1), intent(in) :: flagback
complex(dpc), dimension(Nwaves), intent(in) :: v_state
complex(dpc), dimension(num_vars), intent(out) :: v_sys
complex(dpc), dimension(num_eqs), intent(in) :: rhs

complex(dpc), dimension(num_eqs) :: dev
integer :: k, i, dimA(1)

real(dp) :: err

dimA(1) = num_vars*num_eqs

!compute matrix for current r-value:
call eval_maxwell_system_matrix (r, flagback)

v_sys = cmplx(0.d0,0.d0,dpc)

!filling system and state vector (partly):
do k = 1,3
	do i=0,dim_Ersp_state(k)-1
		v_sys(iErsp_sys(k)+i) = v_state(iErsp_state(k)+i)
	end do
end do

!only last derivatives of E are left unknown: solve the system:
call solve_sub_system (num_eqs, num_vars, A, num_eqs, elist, ulist, v_sys, rhs)

!check the solution:
if (flag_debug == 1) then
    dev = matmul(A, v_sys)-rhs
    err = maxval(abs(dev))/maxval(abs(v_sys))/maxval(reshape(abs(A),dimA))

    if (err > 1.0d-8) then
        print *, 'state2sys: warning: error estimation =', err, 'at r =', r
    end if
end if

end subroutine state2sys

!------------------------------------------------------------------------------

subroutine calc_diff_sys_matrix_slow (r, flagback, Dmat)

!backcase: case for background
!r: radial value; calculate background and derivatives at this value
!D: matrix

use constants, only: dp, dpc
use flre_sett, only: Nwaves

implicit none

real(dp), intent(in) :: r
character(1), intent(in) :: flagback
complex(dpc), dimension(Nwaves, Nwaves), intent(out) :: Dmat

complex(dpc), dimension(Nwaves) :: x

integer :: i

x = cmplx(0.0d0,0.0d0,dpc)

do i=1,Nwaves
    x(i) = cmplx(1.0d0,0.0d0,dpc)
    call calc_dstate_vector (r, flagback, x, Dmat(:,i))
    x(i) = cmplx(0.0d0,0.0d0,dpc)
end do

end subroutine

!------------------------------------------------------------------------------

subroutine calc_diff_sys_matrix (r, flagback, Dmat)

!r: radial value; calculate background and derivatives at this value
!flagback: case for background
!D: ODE matrix

use constants, only: dp, dpc
use flre_sett, only: Nwaves, flag_debug
use maxwell_equations, only: iErsp_state, num_vars, num_eqs, iErsp_sys, dim_Ersp_state
use maxwell_equations, only: A, elist, ulist

implicit none;

real(dp), intent(in) :: r
character(1), intent(in) :: flagback
complex(dpc), dimension(Nwaves, Nwaves), intent(out) :: Dmat

complex(dpc), dimension(Nwaves) :: x
complex(dpc), dimension(num_eqs) :: rhs
complex(dpc), dimension(num_vars) :: v_sys

integer :: i, k, j

complex(dpc), dimension(num_eqs) :: dev
integer :: dimA(1)
real(dp) :: err

dimA(1) = num_vars*num_eqs

!compute matrix for current r-value:
call eval_maxwell_system_matrix (r, flagback)

x = cmplx(0.0d0, 0.0d0, dpc)
rhs = cmplx(0.0d0, 0.0d0, dpc)
v_sys = cmplx(0.d0, 0.d0, dpc)

do i=1,Nwaves
    x(i) = cmplx(1.0d0, 0.0d0, dpc)

    !call calc_dstate_vector (r, flagback, x, Dmat(:,i))
    do k = 1,3
        do j=0,dim_Ersp_state(k)-2
            Dmat(iErsp_state(k)+j, i) = x(iErsp_state(k)+j+1) !derivatives
            v_sys(iErsp_sys(k)+j) = x(iErsp_state(k)+j)
        end do
        j = dim_Ersp_state(k)-1
        if (dim_Ersp_state(k) > 0) then
            v_sys(iErsp_sys(k)+j) = x(iErsp_state(k)+j)
        end if
    end do

    !only last derivatives of E are left unknown: solve the system: might be optimized:
    call solve_sub_system (num_eqs, num_vars, A, num_eqs, elist, ulist, v_sys, rhs)

    !check the solution:
    if (flag_debug == 1) then
        dev = matmul(A, v_sys)-rhs
        err = maxval(abs(dev))/maxval(abs(v_sys))/maxval(reshape(abs(A),dimA))

        if (err > 1.0d-15) then
            print *, 'state2sys: warning: error estimation: ', err
        end if
    end if

    !set last derivatives:
    do k = 1,3
        j = dim_Ersp_state(k)-1
        if (dim_Ersp_state(k) > 0) then
            Dmat(iErsp_state(k)+j, i) = v_sys(iErsp_sys(k)+j+1)
        end if
    end do

    x(i) = cmplx(0.0d0, 0.0d0, dpc)
end do

!debug: add perturbations
! do i=1,Nwaves
!     do j=1,Nwaves
!         Dmat(i,j) = Dmat(i,j) + ((-1.0)**i)*((-1.0)**j)*(1.0d-5)*Dmat(i,j);
!     end do;
! end do;

end subroutine

!------------------------------------------------------------------------------

subroutine calc_dstate_vector (r, flagback, v_state, dv_state)

!r: radial value;
!v_state: vector of state; use to calculate derivatives
!output: dv_state: derivatives of v_state at radial value r

!v_state=Er,...,Er^(2N-3),Es,...,Es^(2N-1),Ep,...,Ep^(2N-1): 6N-2 modes
!v_state=Er,...,Er^(2N-3),Es,...,Es^(2N-3),Ep,...,Ep^(2N-1): 6N-4 modes

use constants, only: dp, dpc
use flre_sett, only: Nwaves
use maxwell_equations, only: iErsp_state, num_vars, num_eqs, iErsp_sys, dim_Ersp_state

implicit none

real(dp), intent(in) :: r
character(1), intent(in) :: flagback
complex(dpc), dimension(Nwaves), intent(in) :: v_state
complex(dpc), dimension(Nwaves), intent(out) :: dv_state

complex(dpc), dimension(num_vars) :: v_sys
complex(dpc), dimension(num_eqs) :: rhs

integer :: k, i

do k = 1,3
	do i=0,dim_Ersp_state(k)-2
		dv_state(iErsp_state(k)+i) = v_state(iErsp_state(k)+i+1)
	end do
end do

rhs = cmplx(0.0d0, 0.0d0, dpc)

!only last derivatives of E are left unknown: solve the system:
call state2sys (r, flagback, v_state, v_sys, rhs)

!set last derivatives:
do k = 1,3
	i = dim_Ersp_state(k)-1
	if (dim_Ersp_state(k) > 0) then
		dv_state(iErsp_state(k)+i) = v_sys(iErsp_sys(k)+i+1)
	end if
end do

end subroutine

!------------------------------------------------------------------------------

subroutine eval_maxwell_system_coeffs (r, flagback)

use constants, only: dp, dpc, pi, im, c
use core, only: bp_ptr
use maxwell_equations
use conduct_parameters

implicit none

real(dp), intent(in) :: r
character(1), intent(in) :: flagback

integer :: i, j, k
complex(dpc) :: fac

real(dp) :: dhs, ddhs, dddhs
real(dp), dimension(2) :: dddh;

call eval_and_set_background_parameters_spec_independent (r, flagback);
call eval_and_set_wave_parameters (r, flagback);

cio = c/(im*omega_)

Ns = c/omega_*ks_
Np = c/omega_*kp_

dNs = c/omega_*dks_
dNp = c/omega_*dkp_

call eval_hthz (r, 3, 3, bp_ptr, dddh); !3th derivative dddhth dddhz

dhs = hz_*dht_ - ht_*dhz_
ddhs = hz_*ddht_- ht_*ddhz_
dddhs = dhz_*ddht_ + hz_*dddh(1)- dht_*ddhz_- ht_*dddh(2)

N1  = cio*(dhs-ht_*hz_/r_)
N2  = cio*ht_*ht_/r_
N3  = cio*(dhs+ht_*hz_/r_)
N4  = cio*hz_*hz_/r_

dN1 = cio*(ddhs-dht_*hz_/r_-ht_*dhz_/r_+ht_*hz_/r_/r_)
dN2 = cio*(2.0d0*ht_*dht_/r_-ht_*ht_/r_/r_)
dN3 = cio*(ddhs+dht_*hz_/r_+ht_*dhz_/r_-ht_*hz_/r_/r_)
dN4 = cio*(2.0d0*dhz_*hz_/r_-hz_*hz_/r_/r_)

ddNs = c/omega_*ddks_

ddN3 = cio*(dddhs+ddht_*hz_/r_+2.0d0*dht_*dhz_/r_-2.0d0*dht_*hz_/r_/r_+ht_*ddhz_/r_ - 2.0d0*ht_*dhz_/r_/r_+2.0d0*ht_*hz_/r_/r_/r_);

ddN4 = cio*(2.0d0*dhz_*dhz_/r_-4.0d0*hz_*dhz_/r_/r_+2.0d0*hz_*ddhz_/r_+2.0d0*hz_*hz_/r_/r_/r_);

!conductivity:
fac = 4.0d0*pi/im/omega_

!index ordering is changed to fit to the splines ordering:
call ctensor (0, 0, r, flagback, cti) !conductivity for ions, type=1 (ct)
call ctensor (1, 0, r, flagback, cte) !conductivity for electrons, type=1 (ct)

do i=1,3
    do j=1,3
        do k=0,der_order(i,j)
            epst(0,k,i,j) = unit3(k,i,j) - fac*(cti(j,i,k,0) + cte(j,i,k,0))
            epst(1,k,i,j) = - fac*(cti(j,i,k,1) + cte(j,i,k,1))
        end do
    end do
end do

end subroutine

!------------------------------------------------------------------------------
