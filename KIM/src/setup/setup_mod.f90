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

end module
