module grid
    implicit none

    integer :: k_space_dim
    logical :: reduce_r
    integer :: reduced_r_dim
    integer :: spline_base
    integer :: grid_spacing
    double precision, dimension(:), allocatable :: kr  ! radial wavenumber
    double precision, dimension(:), allocatable :: krp ! radial wavenumber prime

    double complex, dimension(:,:), allocatable :: varphi_lkr
    
end module