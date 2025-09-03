subroutine kim_init
    use species_m, only: read_profiles, allocate_plasma, init_plasma, plasma

    implicit none

    call read_config
    call allocate_plasma
    call init_plasma(plasma)
    call read_profiles()

end subroutine
