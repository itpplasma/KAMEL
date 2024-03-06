!<Everything needed to build the system of Maxwell equations.

module maxwell_equations

use constants, only: dp, dpc

implicit none

integer :: num_vars, num_eqs

integer, dimension(3) :: max_der_order_Ersp !maximal order of derivs

integer, dimension(3,3) :: der_order !order of derivs of E in a current (rsp)

character (len=3), allocatable, dimension(:) :: names_sys
character (len=3), allocatable, dimension(:) :: names_state

character (len=1), dimension(3) :: names_comp = (/'r','s','p'/)

!indices:
integer :: iEr, iEs, iEp, iBr, idBr, iBs, idBs, iBp, idBp, iddBp

!number of an appropriate components in a state vector:
integer, dimension(3) :: dim_Ersp_state, iErsp_state

!put in arrays:
integer, dimension(3) :: dim_Ersp_sys, iErsp_sys, dim_Brsp_sys, iBrsp_sys

!complex(dpc), allocatable, dimension(:) :: v_sys !system vector

complex(dpc), allocatable, dimension(:,:) :: A

complex(dpc), allocatable, dimension(:,:) :: D !complex matrix: u' = D(r)*u + f(r)

complex(dpc), allocatable, dimension(:,:,:,:) :: cti, cte
complex(dpc), allocatable, dimension(:,:,:,:) :: epst

complex(dpc), allocatable, dimension(:,:,:) :: unit3

!for the system of Maxwell equations:
complex(dpc) :: Ns, Np, dNs, dNp, ddNs
complex(dpc) :: cio, N1, N2, N3, N4, dN1, dN2, dN3, dN4, ddN3, ddN4

!optimization is possible!
integer, allocatable, dimension(:) :: elist, ulist

!indices of a state vector components in a system vector:
integer, allocatable, dimension(:) :: sys_ind

end module

!------------------------------------------------------------------------------

subroutine copy_module_data_to_maxwell_eqs_data_struct_f (num_vars_p, num_eqs_p, &
dim_Ersp_sys_p, iErsp_sys_p, dim_Brsp_sys_p, iBrsp_sys_p, der_order_p)

use maxwell_equations, only: num_vars, num_eqs
use maxwell_equations, only: dim_Ersp_sys, iErsp_sys
use maxwell_equations, only: dim_Brsp_sys, iBrsp_sys
use maxwell_equations, only: der_order

implicit none;

integer, intent(out) :: num_vars_p, num_eqs_p
integer, dimension(3), intent(out) :: dim_Ersp_sys_p, iErsp_sys_p
integer, dimension(3), intent(out) :: dim_Brsp_sys_p, iBrsp_sys_p
integer, dimension(3,3), intent(out) :: der_order_p

num_vars_p = num_vars
num_eqs_p = num_eqs

dim_Ersp_sys_p = dim_Ersp_sys
iErsp_sys_p = iErsp_sys

dim_Brsp_sys_p = dim_Brsp_sys
iBrsp_sys_p = iBrsp_sys

der_order_p = transpose(der_order) !C and Fortran have different array's elements ordering

end subroutine

!--------------------------------------------------------------------

subroutine get_Ersp_state_indices_and_dims_f (dim_Ersp_state_p, iErsp_state_p)

use maxwell_equations, only: dim_Ersp_state, iErsp_state

implicit none;

integer, dimension(3), intent(out) :: dim_Ersp_state_p, iErsp_state_p

dim_Ersp_state_p = dim_Ersp_state
iErsp_state_p = iErsp_state

end subroutine

!--------------------------------------------------------------------

subroutine get_sys_ind_array_f (sys_ind_p)

use flre_sett, only: Nwaves
use maxwell_equations, only: sys_ind

integer, dimension(Nwaves), intent(out) :: sys_ind_p

sys_ind_p = sys_ind

end subroutine

!--------------------------------------------------------------------

subroutine clean_maxwell_system_parameters_module ()

use maxwell_equations;

implicit none;

num_vars = 0;
num_eqs = 0;

max_der_order_Ersp = 0;

der_order = 0;

deallocate (names_sys, names_state);

dim_Ersp_state = 0;
iErsp_state = 0;

dim_Ersp_sys = 0;
iErsp_sys = 0;

dim_Brsp_sys = 0;
iBrsp_sys = 0;

deallocate (A, D, cti, cte, epst, unit3);

deallocate (elist, ulist);

deallocate (sys_ind);

end subroutine

!--------------------------------------------------------------------
