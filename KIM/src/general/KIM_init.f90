subroutine kim_init

    use grid, only: reduce_r
    use plasma_parameter, only: read_profiles
    use equilibrium, only: calculate_equil

    implicit none

    call read_config
    call read_profiles(reduce_r)
    ! calculate equilibrium B field and J
    call calculate_equil(.true.)
    ! calculate quantities used for the kernels, e.g. A1, A2, dndr, omega_c,...
    call allocate_backs
    call calculate_backs(.true.)

end subroutine