program test_periodic_config_read
    use KIM_kinds_m, only: dp
    use config_m, only: nml_config_path, periodic_Bparallel_ratio, &
        periodic_dr_asis_scale, periodic_dr_tr_scale, periodic_kmax_scale, &
        periodic_n_rg

    implicit none

    complex(dp), parameter :: expected_ratio = cmplx(0.25_dp, -0.5_dp, dp)
    character(len=512) :: config_path, second_config_path
    integer :: argument_count, environment_status
    logical :: defaults_reset_case, should_reject

    argument_count = command_argument_count()
    defaults_reset_case = argument_count == 0
    select case (argument_count)
    case (0)
        call get_environment_variable('KIM_TEST_CONFIG_PATH', config_path, &
            status=environment_status)
        if (environment_status /= 0 .or. len_trim(config_path) == 0) then
            error stop 'KIM_TEST_CONFIG_PATH is required for the reset test'
        end if
        call get_environment_variable('KIM_TEST_SECOND_CONFIG_PATH', &
            second_config_path, status=environment_status)
        if (environment_status /= 0 .or. len_trim(second_config_path) == 0) then
            error stop 'KIM_TEST_SECOND_CONFIG_PATH is required for the reset test'
        end if
    case (1)
        call get_command_argument(1, config_path)
    case default
        error stop 'usage: test_periodic_config_read <config-path>'
    end select

    nml_config_path = trim(config_path)
    should_reject = index(config_path, 'malformed') > 0 .or. &
        index(config_path, 'unterminated') > 0

    call kim_read_config
    if (.not. should_reject .and. periodic_Bparallel_ratio /= expected_ratio) then
        error stop 'valid KIM_PERIODIC ratio was not read exactly'
    end if

    if (defaults_reset_case) then
        if (periodic_dr_asis_scale /= 7.0_dp .or. &
                periodic_dr_tr_scale /= 13.0_dp .or. &
                periodic_kmax_scale /= 9.0_dp .or. periodic_n_rg /= 17) then
            error stop 'custom KIM_PERIODIC values were not read exactly'
        end if

        nml_config_path = trim(second_config_path)
        call kim_read_config
        if (periodic_dr_asis_scale /= 5.0_dp .or. &
                periodic_dr_tr_scale /= 10.0_dp .or. &
                periodic_kmax_scale /= 5.0_dp .or. periodic_n_rg /= 96 .or. &
                periodic_Bparallel_ratio /= (0.0_dp, 0.0_dp)) then
            error stop 'omitted KIM_PERIODIC group did not restore defaults'
        end if
        print *, 'periodic config default-reset assertions PASSED'
    end if
end program test_periodic_config_read
