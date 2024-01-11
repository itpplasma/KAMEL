module grid
    implicit none

    integer :: k_space_dim ! dimension of kr space grid
    integer :: l_space_dim ! dimension of spline grid
    integer :: r_space_dim ! dimension of r grid
    logical :: reduce_r
    integer :: reduced_r_dim
    integer :: spline_base
    integer :: grid_spacing
    integer :: num_gengrid_points

    double precision, dimension(:), allocatable :: kr  ! radial wavenumber
    double precision, dimension(:), allocatable :: krp ! radial wavenumber prime

    double precision, dimension(:), allocatable :: xl  ! xl grid (real space)

    double complex, dimension(:,:), allocatable :: varphi_lkr

    integer :: npoib, npoic, npoi_der, nbaleqs, neqset, iboutype
    integer :: mwind
    double precision :: rmin,rmax
    double precision :: gg_factor = 50
    double precision :: gg_width = 2 
    double precision :: gg_r_res = 95.34
    double precision :: rb_cut_in, re_cut_in, rb_cut_out, re_cut_out
    integer,          dimension(:),   allocatable :: ipbeg, ipend
    double precision, dimension(:),   allocatable :: rb, rc, Sb, Sc
    double precision, dimension(:),   allocatable :: y, dery, dery_equisource
    double precision, dimension(:),   allocatable :: r_resonant
    double precision, dimension(:,:), allocatable :: deriv_coef, reint_coef
    
end module