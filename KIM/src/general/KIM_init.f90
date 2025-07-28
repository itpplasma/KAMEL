subroutine kim_init

    use grid, only: reduce_r
    use species, only: read_profiles, allocate_plasma, init_plasma, plasma

    implicit none

    call read_config
    call allocate_plasma
    call init_plasma(plasma)
    call read_profiles(reduce_r)

end subroutine