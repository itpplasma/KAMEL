module setup

    use KIM_kinds, only: dp

    implicit none

    real(dp) :: btor   ! toroidal magnetic field
    real(dp) :: R0     ! major radius
    real(dp) :: r_plas
    integer :: set_profiles_constant
    integer          :: m_mode ! poloidal mode number
    integer          :: n_mode ! toroidal mode number
    real(dp) :: omega  ! perturbation frequency
    real(dp) :: cut_off_fac ! factor for cut off distance in creating the sparse matrix
    real(dp) :: kr_cut_off_fac ! 
    integer          :: type_br_field ! type of the radial magnetic field, 1 for constant, 2 for point charge case (no multiplication with kernel), 3 for linear increase
    logical          :: collisions_off
    real(dp) :: eps_reg

end module