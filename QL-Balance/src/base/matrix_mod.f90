
module matrix_mod
    integer :: isw_rhs
    integer :: nz,nsize
    integer,          dimension(:),   allocatable :: irow,icol
    double precision, dimension(:),   allocatable :: amat,rhsvec
end module matrix_mod
