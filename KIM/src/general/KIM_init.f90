subroutine kim_init

    use grid, only: reduce_r
    use species, only: read_profiles, allocate_plasma, init_plasma, plasma
    use equilibrium, only: calculate_equil

    implicit none

    call read_config
    call allocate_plasma
    call init_plasma(plasma)
    call read_profiles(reduce_r)
    ! calculate equilibrium B field and J
    !call calculate_equil(.false.)
    ! calculate quantities used for the kernels, e.g. A1, A2, dndr, omega_c,...
    !call allocate_backs
    !call calculate_backs(.false.)

end subroutine