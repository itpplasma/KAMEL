!<Subroutines used to build Maxwell system matrix in the original form.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_and_set_maxwell_system_parameters_module ()

!system of 9 equations is for any N>0 and inhomogenious background and
!for N=1 and homogenious background

!system vector for 9 eqs: N=1 and inhomogenious case for N>1
!Er, dEr,...,dEr^(2*flre_order-1)
!Es, dEs,...,dEs^(2*flre_order)
!Ep, dEp,...,dEp^(2*flre_order)
!Br, Bs, dBs, Bp, dBp

!state vector for diff_sys:
!v_state=Er,...,Er^(2N-3),Es,...,Es^(2N-1),Ep,...,Ep^(2N-1): 6N-2 modes

use flre_sett, only: flre_order, Nwaves
use maxwell_equations

implicit none

integer :: k,i,l,n

num_vars = 6*flre_order+7
num_eqs = 9

!allocate (v_sys(1:num_vars));

allocate (A(1:num_eqs,1:num_vars));
A = cmplx(0.0d0,0.0d0, dpc)

allocate (D(1:Nwaves,1:Nwaves));

allocate (cti(1:3,1:3,0:2*flre_order,0:1));
allocate (cte(1:3,1:3,0:2*flre_order,0:1));

allocate (epst(0:1,0:2*flre_order,1:3,1:3));
epst = cmplx(0.0d0,0.0d0, dpc)

allocate (unit3(0:2*flre_order,1:3,1:3));
unit3 = (0.d0,0.d0)
unit3(0,1,1) = (1.d0,0.d0)
unit3(0,2,2) = (1.d0,0.d0)
unit3(0,3,3) = (1.d0,0.d0)

!maximal order of derivs
max_der_order_Ersp(1) = 2*flre_order-1
max_der_order_Ersp(2) = 2*flre_order
max_der_order_Ersp(3) = 2*flre_order

!order of derivs of E in a current (rsp)
der_order(1,1) = 2*flre_order-2
der_order(1,2) = 2*flre_order-1
der_order(1,3) = 2*flre_order-1

der_order(2,1) = 2*flre_order-1
der_order(2,2) = 2*flre_order
der_order(2,3) = 2*flre_order

der_order(3,1) = 2*flre_order-1
der_order(3,2) = 2*flre_order
der_order(3,3) = 2*flre_order

!system vector:
!Er, dEr,...,dEr^(2*flre_order-1)
!Es, dEs,...,dEs^(2*flre_order)
!Ep, dEp,...,dEp^(2*flre_order)
!Br, Bs, dBs, Bp, dBp

!indices:

iEr  = 1
iEs  = iEr+max_der_order_Ersp(1)+1
iEp  = iEs+max_der_order_Ersp(2)+1
iBr  = iEp+max_der_order_Ersp(3)+1
idBr = -1 !not needed for 9 eqs
iBs  = iBr+1
idBs = iBs+1
iBp  = idBs+1
idBp = iBp+1
iddBp= -1 !not needed for 9 eqs

dim_Ersp_sys(1) = iEs-iEr
dim_Ersp_sys(2) = iEp-iEs
dim_Ersp_sys(3) = iBr-iEp

dim_Brsp_sys(1) = iBs-iBr
dim_Brsp_sys(2) = iBp-iBs
dim_Brsp_sys(3) = num_vars-iBp+1

!put in arrays:
iErsp_sys(1) = iEr
iErsp_sys(2) = iEs
iErsp_sys(3) = iEp
iBrsp_sys(1) = iBr
iBrsp_sys(2) = iBs
iBrsp_sys(3) = iBp

allocate (names_sys(1:num_vars), names_state(Nwaves));

do k = 1,3
    do i=0,dim_Ersp_sys(k)-1
        write (names_sys(iErsp_sys(k)+i),'(a,a,i1)') 'E', names_comp(k), i
    end do
end do

names_sys(iBrsp_sys(1))  = 'Br0'
names_sys(iBrsp_sys(2))  = 'Bs0'
names_sys(iBrsp_sys(2)+1)= 'Bs1'
names_sys(iBrsp_sys(3))  = 'Bp0'
names_sys(iBrsp_sys(3)+1)= 'Bp1'

! print *, 'system vector:'
! do k = 1,num_vars
!     print *, 'sys(',k,') = ', names_sys(k)
! end do

!state vector's stuff for diff_sys:
!v_state=Er,...,Er^(2N-3),Es,...,Es^(2N-1),Ep,...,Ep^(2N-1): 6N-2 modes

dim_Ersp_state(1) = max_der_order_Ersp(1)-1
dim_Ersp_state(2) = max_der_order_Ersp(2)
dim_Ersp_state(3) = max_der_order_Ersp(3)

!indices of quants in a state vector (are different from sys vector!):

iErsp_state(1) = 1 !doesn't matter if dim(1)=0
iErsp_state(2) = iErsp_state(1)+dim_Ersp_state(1)
iErsp_state(3) = iErsp_state(2)+dim_Ersp_state(2)

if(iErsp_state(3) + dim_Ersp_state(3)-1 /= Nwaves) then
    print *, 'set_parameters: error: wrong state vector!..'
end if

do k = 1,3
    do i=0,dim_Ersp_state(k)-1
        write (names_state(iErsp_state(k)+i),'(a,a,i1)') 'E', names_comp(k), i
    end do
end do

! print *, 'state vector:'
! do i=1,Nwaves
!     print *, 'state(',i,') = ', names_state(i)
! end do

!all 9 eqs are solved: optimization is possible!
allocate (elist(num_eqs), ulist(num_eqs))

elist = (/1,2,3,4,5,6,7,8,9/)

ulist = (/iErsp_sys(2)-2, iErsp_sys(2)-1, iErsp_sys(3)-1, &
          iBrsp_sys(1)-1, iBrsp_sys(1), iBrsp_sys(2), &
          iBrsp_sys(2)+1, iBrsp_sys(3), iBrsp_sys(3)+1/)

!indices of state fields in system fields arrray:
allocate (sys_ind(1:Nwaves));

n = 1;
do k = 1,3
    do i=0, dim_Ersp_state(k)-1

        l = iErsp_state(k) + i;
        sys_ind(l) = iErsp_sys(k) + i;

        if (l /= n) then
            print *, 'error: calc_and_set_maxwell_system_parameters_module: indexing error: n, l =', n, l;
        end if

        n = n + 1;

    end do
end do

if (l /= Nwaves) then
    print *, 'error: calc_and_set_maxwell_system_parameters_module: indexing error: l =', l;
    stop;
end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine eval_maxwell_system_matrix (r, flagback)

use constants, only: dp
use maxwell_equations

implicit none

real(dp), intent(in) :: r
character(*), intent(in) :: flagback

integer :: i, k

call eval_maxwell_system_coeffs (r, flagback)

!1: Br
A(1,iEs) = -Np
A(1,iEp) = Ns
A(1,iBr) = -1.d0

!2: dEp
A(2,iEp+1) = -cio
A(2,iEr) = Np
A(2,iEs) = N1
A(2,iEp) = -N2
A(2,iBs) = -1.d0

!3: dEs
A(3,iEs+1) = cio
A(3,iEr) = -Ns
A(3,iEp) = N3
A(3,iEs) = N4
A(3,iBp) = -1.d0

!4:
A(4,iBs) = -Np
A(4,iBp) = Ns
do k=1,3 !r,s,p
	do i=0,der_order(1,k)
        	A(4,iErsp_sys(k)+i) = epst(0,i,1,k)
	end do
end do

!5: dBp
A(5,idBp) = -cio
A(5,iBr) = Np
A(5,iBs) = N1
A(5,iBp) = -N2

do k=1,3 !r,s,p
	do i=0,der_order(2,k)
        	A(5,iErsp_sys(k)+i) = epst(0,i,2,k)
	end do
end do

!6: dBs
A(6,idBs) = cio
A(6,iBr) = -Ns
A(6,iBp) = N3
A(6,iBs) = N4

do k=1,3 !r,s,p
	do i=0,der_order(3,k)
        	A(6,iErsp_sys(k)+i) = epst(0,i,3,k)
	end do
end do

!DIFF, ddEp
A(7,iEp+2) = -cio
A(7,iEr)   = dNp
A(7,iEr+1) = Np
A(7,iEs)   = dN1
A(7,iEs+1) = N1
A(7,iEp)   = -dN2
A(7,iEp+1) = -N2
A(7,idBs)  = -1.d0

!DIFF, ddEs
A(8,iEs+2) = cio
A(8,iEr)   = -dNs
A(8,iEr+1) = -Ns
A(8,iEp)   = dN3
A(8,iEp+1) = N3
A(8,iEs)   = dN4
A(8,iEs+1) = N4
A(8,idBp)  = -1.d0

!DIFF, ddEr
A(9,:) = (0.0d0,0.0d0)
A(9,idBs) = -Np
A(9,iBs)  = -dNp
A(9,idBp) = Ns
A(9,iBp)  = dNs

do k=1,3 !r,s,p
	do i=0,der_order(1,k)
        	A(9,iErsp_sys(k)+i+1) = epst(0,i,1,k)
		A(9,iErsp_sys(k)+i) = A(9,iErsp_sys(k)+i) + epst(1,i,1,k)
	end do
end do

!print *, 'r:', r
!print *, 'A1:', A(1,:)
!print *, 'A2:', A(2,:)
!print *, 'A3:', A(3,:)
!print *, 'A4:', A(4,:)
!print *, 'A5:', A(5,:)
!print *, 'A6:', A(6,:)
!print *, 'A7:', A(7,:)
!print *, 'A8:', A(8,:)
!print *, 'A9:', A(9,:)

end subroutine

!------------------------------------------------------------------------------

subroutine calc_diff_sys_rhs_vector (r, flagback, f)

!r: radial value; calculate background and derivatives at this value
!flagback: case for background
!du: rhs vector for the inhomogenious system u' = Dmat*u + f

use constants, only: dp, dpc, c, im, pi, isqrt2pi
use mode_data, only: omov
use flre_sett, only: Nwaves
use antenna_data, only: ra, wa
use maxwell_equations, only: iErsp_state, num_vars, num_eqs, iErsp_sys, dim_Ersp_state

implicit none;

real(dp), intent(in) :: r
character(*), intent(in) :: flagback
complex(dpc), dimension(Nwaves), intent(out) :: f

complex(dpc), dimension(num_eqs) :: rhs
complex(dpc), dimension(num_vars) :: v_sys
complex(dpc), dimension(Nwaves) :: v_state

real(dp), parameter :: fpc = 4.0d0*pi/c

complex(dpc), dimension(2) :: Ja_rsp

complex(dpc) :: cio;

integer :: i, k

cio = c/im/omov;

!(r,s,p) antenna current density Fourier amplutudes:
call calc_current_density_r_s_p (r, Ja_rsp);

rhs = cmplx(0.0d0, 0.0d0, dpc);
rhs(5) = cio*fpc*Ja_rsp(1);
rhs(6) = cio*fpc*Ja_rsp(2);

v_state = cmplx(0.0d0, 0.0d0, dpc);

call state2sys (r, flagback, v_state, v_sys, rhs);

!set state vector derivatives:
do k=1,3
    do i=0,dim_Ersp_state(k)-1
        f(iErsp_state(k)+i) = v_sys(iErsp_sys(k)+i+1)
    end do
end do

end subroutine

!------------------------------------------------------------------------------
