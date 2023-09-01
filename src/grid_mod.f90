module grid
    implicit none

    integer :: k_space_dim
    logical :: reduce_r
    integer :: reduced_r_dim
    double precision, dimension(:), allocatable :: kr  ! radial wavenumber
    double precision, dimension(:), allocatable :: krp ! radial wavenumber prime

end module