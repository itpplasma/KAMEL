program test_periodic_checkpoint
    use QLBalance_kinds, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use KAMEL_hdf5_tools
    use periodic_amplitude_state_m, only: periodic_amplitude_state_t
    use periodic_checkpoint_m, only: periodic_restart_t, write_periodic_checkpoint, &
                                     read_periodic_checkpoint
    implicit none
    type(periodic_amplitude_state_t) :: state, restored
    type(periodic_restart_t) :: source, loaded
    integer(HID_T) :: file_id
    integer :: i
    logical :: found
    character(32) :: argument
    character(128) :: filename
    complex(dp) :: amps(2), currents(2), residuals(2)

    call get_command_argument(1, argument)
    filename = 'periodic_checkpoint_'//trim(argument)//'.h5'
    amps = [(1.0_dp, 2.0_dp), (3.0_dp, -4.0_dp)]
    currents = [(5.0_dp, 6.0_dp), (7.0_dp, 8.0_dp)]
    residuals = [(0.1_dp, 0.2_dp), (0.3_dp, 0.4_dp)]
    call state%initialize(amps, currents, residuals, [0, 0], 0.25_dp, 0.8_dp)
    source%active = .true.
    source%accepted_step = 7
    source%time = 0.12345678901234567_dp
    source%next_dt = 1.234567890123456e-5_dp
    source%tolerance = 0.031234567890123_dp
    source%rc = [1.5_dp, 2.5_dp]
    source%rb = [1.0_dp, 2.0_dp, 3.0_dp]
    source%params = reshape([(1.123456789012345_dp * real(i, dp), i=1, 8)], [4, 2])
    source%q = [2.0_dp, 3.0_dp, 4.0_dp]
    source%er = [-0.11_dp, -0.22_dp, -0.33_dp]
    source%vth = [1.0_dp, 2.0_dp, 3.0_dp]
    source%dery_equisource = [(1.23456789012345_dp * real(i, dp), i=1, 8)]
    source%accepted_floor = 1.0e-29_dp
    source%accepted_max_scale = 1.0e11_dp
    source%timstep_arr = [(source%next_dt, i=1, 8)]
    source%tim_stack = [(0.01_dp * real(i, dp), i=1, 8)]
    source%y_stack = reshape([(real(i, dp), i=1, 16)], [8, 2])
    source%br_abs = [(real(i, dp), i=1, 7)]
    source%br_abs_time = source%br_abs * 0.01_dp
    source%br_abs_antenna_factor = source%br_abs * 0.1_dp
    source%dqle22_res_time = source%br_abs * 0.2_dp
    source%dae22_res_time = source%br_abs * 0.3_dp
    source%bif_criterion = source%br_abs * 0.4_dp
    source%Ipar_time = cmplx(source%br_abs, -source%br_abs, dp)
    source%br_formfactor = source%Ipar_time * 0.2_dp
    source%br_vac_res = source%Ipar_time * 0.3_dp
    source%nstack = 2
    source%timscal_dql = 0.123456789012345_dp
    source%timscal_dqli = 0.098765432109876_dp
    source%rate_dql = 98765.432109876_dp
    source%ramp_mode = 6
    source%ramp_up_down = 1
    source%iexit = 2
    source%br_beta = 0.4_dp
    source%br_predicted = 0.7_dp
    source%t_hysteresis_turn = 0.08_dp
    source%antenna_factor = 0.81_dp
    source%target_current = 0.5_dp
    source%current_floor = 1.0e-30_dp
    source%max_scale = 1.0e12_dp
    source%relaxation = 0.6_dp
    call h5_init()
    call h5_create(trim(filename), file_id)
    call write_periodic_checkpoint(file_id, '/multi_mode/KinProfiles/1007/', &
                                   [-6, -7], [2, 2], state, source)
    source%params = source%params * 1.123456789012345_dp
    source%y_stack = reshape([(real(i, dp), i=1, 24)], [8, 3])
    source%nstack = 3
    call write_periodic_checkpoint(file_id, '/multi_mode/KinProfiles/1007', &
                                   [-6, -7], [2, 2], state, source)
    if (trim(argument) == 'version') then
        call h5_delete(file_id, '/multi_mode/KinProfiles/1007/PeriodicCheckpoint/version')
        call h5_add(file_id, '/multi_mode/KinProfiles/1007/PeriodicCheckpoint/version', 999)
    end if
    if (trim(argument) == 'phase') then
        call h5_delete(file_id, '/multi_mode/KinProfiles/1007/PeriodicCheckpoint/phase')
        call h5_add_string(file_id, &
                           '/multi_mode/KinProfiles/1007/PeriodicCheckpoint/phase', 'wrong')
    end if
    if (trim(argument) == 'status') then
        call h5_delete(file_id, '/multi_mode/KinProfiles/1007/PeriodicCheckpoint/status')
        call h5_add(file_id, '/multi_mode/KinProfiles/1007/PeriodicCheckpoint/status', &
                    [0, 2], [1], [2])
    end if
    if (trim(argument) == 'duplicate') then
        call h5_delete(file_id, '/multi_mode/KinProfiles/1007/PeriodicCheckpoint/m')
        call h5_add(file_id, '/multi_mode/KinProfiles/1007/PeriodicCheckpoint/m', &
                    [-6, -6], [1], [2])
    end if
    if (trim(argument) == 'nonfinite') then
        source%params(1, 1) = ieee_value(1.0_dp, ieee_quiet_nan)
        call h5_delete(file_id, '/multi_mode/KinProfiles/1007/PeriodicCheckpoint/params')
        call h5_add(file_id, '/multi_mode/KinProfiles/1007/PeriodicCheckpoint/params', &
                    source%params, [1, 1], shape(source%params))
    end if
    if (trim(argument) == 'partial') &
        call h5_delete(file_id, '/multi_mode/KinProfiles/1007/PeriodicCheckpoint/params')
    call h5_close(file_id)
    call h5_open(trim(filename), file_id)
    if (trim(argument) == 'modes') then
        call read_periodic_checkpoint(file_id, '/multi_mode/KinProfiles/1007/', &
                                      [-6, -8], [2, 2], restored, loaded, found)
    else
        call read_periodic_checkpoint(file_id, '/multi_mode/KinProfiles/1007/', &
                                      [-7, -6], [2, 2], restored, loaded, found)
    end if
    if (len_trim(argument) > 0) stop 0
    call require(found .and. loaded%active .and. restored%initialized, &
                 'checkpoint not restored')
    call require(all(restored%accepted == amps(2:1:-1)), 'signed mode remapping failed')
    call require(all(restored%accepted_current_unit == currents(2:1:-1)), &
                 'current remapping failed')
    call require(all(restored%accepted_residual == residuals(2:1:-1)), 'residual remapping failed')
    call require(all(restored%trial == restored%accepted), 'restart has an uncommitted trial')
    call require(all(loaded%params == source%params), 'double precision profiles changed')
    call require(all(loaded%rc == source%rc) .and. all(loaded%rb == source%rb), 'grid changed')
    call require(all(loaded%q == source%q) .and. all(loaded%er == source%er), &
                 'background changed')
    call require(loaded%time == source%time .and. loaded%next_dt == source%next_dt &
                 .and. loaded%tolerance == source%tolerance, 'adaptive clock changed')
    call require(all(loaded%timstep_arr == source%timstep_arr) &
                 .and. all(loaded%tim_stack == source%tim_stack) &
                 .and. all(loaded%y_stack == source%y_stack), 'adaptive history changed')
    call require(all(loaded%br_abs == source%br_abs) &
                 .and. all(loaded%Ipar_time == source%Ipar_time), 'stopping history changed')
    call require(loaded%timscal_dql == source%timscal_dql &
                 .and. loaded%timscal_dqli == source%timscal_dqli &
                 .and. loaded%rate_dql == source%rate_dql, 'QL adaptive time scales changed')
    call require(all(loaded%dery_equisource == source%dery_equisource), &
                 'fixed equilibrium source changed')
    call require(loaded%accepted_floor == source%accepted_floor &
            .and. loaded%accepted_max_scale == source%accepted_max_scale, 'accepted guards changed')
    call require(loaded%ramp_mode == 6 .and. loaded%ramp_up_down == 1 &
                 .and. loaded%br_beta == source%br_beta, 'ramp state changed')
    call require(loaded%target_current == 0.5_dp .and. loaded%relaxation == 0.6_dp, &
                 'next-step controls lost to accepted metadata')
    call require(restored%accepted_target_current == 0.25_dp &
                 .and. restored%accepted_relaxation == 0.8_dp, 'accepted metadata changed')
    call read_periodic_checkpoint(file_id, '/old_without_amplitudes/', [-6], [2], &
                                  restored, loaded, found)
    call require(.not. found .and. .not. loaded%active .and. .not. restored%initialized, &
                 'old checkpoint retained stale pending state')
    call h5_close(file_id)
    call h5_deinit()
    print *, 'periodic checkpoint roundtrip passed'
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) error stop message
    end subroutine require
end program test_periodic_checkpoint
