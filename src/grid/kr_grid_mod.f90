module kr_grid

    integer :: k_space_dim ! dimension of kr space grid
    double precision, dimension(:), allocatable :: kr  ! radial wavenumber
    double precision, dimension(:), allocatable :: krp ! radial wavenumber prime
    double precision :: kr_grid_ampl_res
    double precision :: kr_grid_width_res
    double precision :: kr_res = 0.0d0

    integer :: closest_kr_ind_lower
    integer :: closest_kr_ind_upper

end module