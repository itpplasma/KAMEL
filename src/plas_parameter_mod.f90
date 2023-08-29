module plas_parameter

    implicit none
    integer :: iprof_length
    double precision, allocatable :: r_prof(:)
    double precision, allocatable :: n_prof(:)
    double precision, allocatable :: Te_prof(:)
    double precision, allocatable :: Ti_prof(:)
    double precision, allocatable :: Er_prof(:)
    double precision, allocatable :: q_prof(:)

    integer, dimension(:), allocatable :: Zi ! ion charge number
    integer, dimension(:), allocatable :: Ai ! ion mass number
    

end module plas_parameter