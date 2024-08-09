
module recstep_mod
    integer :: nstack
    double precision :: tol
    double precision, dimension(:),   allocatable :: tim_stack
    double precision, dimension(:),   allocatable :: timstep_arr
    double precision, dimension(:,:), allocatable :: y_stack
end module recstep_mod