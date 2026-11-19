program test_kernel_taper_config
    ! Reads the namelist given as argument; kim_read_config must stop on a
    ! kernel_taper_skip_threshold outside (0, 1) and accept the default.
    use config_m, only: nml_config_path
    use grid_m, only: kernel_taper_skip_threshold

    implicit none

    character(len=512) :: config_path

    if (command_argument_count() /= 1) error stop 'usage: test_kernel_taper_config <config-path>'
    call get_command_argument(1, config_path)
    nml_config_path = trim(config_path)
    call kim_read_config
    print *, 'kernel_taper_skip_threshold accepted: ', kernel_taper_skip_threshold
end program test_kernel_taper_config
