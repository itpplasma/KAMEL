module setup_m

    use KIM_kinds_m, only: dp

    implicit none

    real(dp) :: btor   ! toroidal magnetic field
    real(dp) :: R0     ! major radius
    integer :: set_profiles_constant
    integer          :: m_mode ! poloidal mode number
    integer          :: n_mode ! toroidal mode number
    real(dp) :: omega  ! perturbation frequency
    real(dp) :: cut_off_fac ! factor for cut off distance in creating the sparse matrix
    real(dp) :: kr_cut_off_fac !
    integer          :: type_br_field ! type of the radial magnetic field, 1 for constant, 2 for point charge case (no multiplication with kernel), 3 for linear increase
    logical          :: collisions_off
    integer :: bc_type ! boundary condition. 0: None; 1: Neuman left, Dirichlet right
    integer :: spline_base
    integer :: mphi_max ! maximum number of cyclotron harmonics to include in kernel calculations
    real(dp) :: Br_boundary_re = 1.0d0  ! real part of Br at right boundary
    real(dp) :: Br_boundary_im = 0.0d0  ! imaginary part of Br at right boundary

    ! fourier_periodic run type: resonant radius, resonant-layer and
    ! transition half-widths in cm, mode count, one-period quadrature
    ! points, and the periodized-background tabulation size
    real(dp) :: fp_r_res = -1.0d0
    real(dp) :: fp_dr_layer = -1.0d0
    real(dp) :: fp_dr_transition = -1.0d0
    integer :: fp_n_modes = 128
    integer :: fp_n_quad = 512
    integer :: fp_grid_points = 2048

end module
