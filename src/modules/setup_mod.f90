module setup

    implicit none

    double precision :: btor   ! toroidal magnetic field
    double precision :: R0     ! major radius
    integer          :: m_mode ! poloidal mode number
    integer          :: n_mode ! toroidal mode number
    double precision :: omega  ! perturbation frequency
    double precision :: cut_off_fac ! factor for cut off distance in creating the sparse matrix
    double precision :: kr_cut_off_fac ! 
    integer          :: type_br_field ! type of the radial magnetic field, 1 for constant, 2 for point charge case (no multiplication with kernel), 3 for linear increase
    logical          :: collisions_off

end module