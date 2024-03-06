!------------------------------------------------------------------------------

subroutine eigsys (dim, mat, evl, evc)

!Finds eigenvalues and eigenvectors of mat

implicit none;

integer, parameter :: dp = 8, dpc = 8

integer, intent(in) :: dim
complex(dpc), dimension(dim,dim), intent(in)  :: mat
complex(dpc), dimension(dim),   intent(out) :: evl
complex(dpc), dimension(dim,dim), intent(out) :: evc

character :: jobvl = 'N', jobvr = 'V'
integer :: info, lwork

complex(dpc), dimension(2*dim) :: tork !trial work

real(dp), allocatable :: rwork(:)
complex(dpc), allocatable :: A(:,:), work(:)

allocate (A(dim,dim), rwork(2*dim))

A = mat

lwork = -1; !query for optimal lwork

call zgeevt (jobvl, jobvr, dim, A, dim, evl, evc, dim, evc, dim, tork, lwork, rwork, info)

lwork = tork(1) !optimal lwork

allocate (work(lwork))

call zgeevt (jobvl, jobvr, dim, A, dim, evl, evc, dim, evc, dim, work, lwork, rwork, info)

if (info /= 0) then
    print *, 'warning : eigsys: failed to find solution: info=', info
end if

deallocate (A, work, rwork)

end subroutine

!------------------------------------------------------------------------------
