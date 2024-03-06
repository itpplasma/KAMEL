
subroutine solve_sub_system (dim_eqs, dim_vars, A, dim_unk, elist, ulist, x, rhs)

!routine for generating & solving subsystem of the system Ax = b
!dim_eqs, dim_vars, A: system matrix and its dimensions
!elist: index list of desired equations to solve
!ulist: index list of desired unknowns to solve
!x:  full state vector;
!on entry: filled with knowns (at right positions), rest has to be zero
!on exit: positions corresponding to contents of ulist filled with computed solution
!rhs: right hand side vector

use constants, only: dpc

implicit none

integer, intent(in) :: dim_eqs, dim_vars, dim_unk
complex(dpc), dimension(dim_eqs, dim_vars), intent(in) :: A
integer, dimension(dim_unk), intent(in) :: elist, ulist
complex(dpc), dimension(dim_vars), intent(inout) :: x
complex(dpc), dimension(dim_eqs), intent(in) :: rhs

integer :: i, j, info, dimbb = 1

complex(dpc), dimension(dim_unk, dim_unk) :: AA
complex(dpc), dimension(dim_unk) :: bb

integer, dimension(dim_unk) :: IPIV

!make sure that the unknowns are zero
do i=1,dim_unk
    x(ulist(i)) = cmplx(0.0d0, 0.0d0, dpc)
end do

!if only one equation:
if (dim_unk==1) then
    x(ulist(1)) = (rhs(elist(1))-sum(A(elist(1),:)*x))/A(elist(1),ulist(1))
    return
end if

!linear system matrix:
do i=1,dim_unk

    do j=1,dim_unk
        AA(i,j) = A(elist(i),ulist(j))
    end do

    !right hand side
    bb(i) = rhs(elist(i))-sum(A(elist(i),:)*x)
end do

!solve general complex equation:
call zgesv (dim_unk, dimbb, AA, dim_unk, IPIV, bb, dim_unk, info)

if (info /= 0) then
    print *, 'error: solve_sub_system: failed to solve the system: ierr=', info
end if

!assign result back to a state vector
do j=1,dim_unk
    x(ulist(j)) = bb(j)
end do

end subroutine

!------------------------------------------------------------------------------
