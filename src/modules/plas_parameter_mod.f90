module plas_parameter

    implicit none
    integer :: iprof_length
    double precision, allocatable :: r_prof(:)
    double precision, allocatable :: n_prof(:)
    double precision, allocatable :: Te_prof(:) 
    double precision, allocatable :: Er_prof(:)
    double precision, allocatable :: q_prof(:)

    double precision, allocatable :: ni_prof(:, :)
    double precision, allocatable :: Ti_prof(:, :)
    integer, dimension(:), allocatable :: Zi ! ion charge number
    integer, dimension(:), allocatable :: Ai ! ion mass number
    
    double precision, allocatable :: dndr_prof(:)
    double precision, allocatable :: dTedr_prof(:)
    double precision, allocatable :: dTidr_prof(:,:)
    double precision, allocatable :: dqdr_prof(:)
    double precision, allocatable :: dnidr_prof(:, :)

end module plas_parameter