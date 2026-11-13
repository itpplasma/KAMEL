!> The NTV torque is scaled with the antenna factor the QL coefficients use in
!> the same step, repeated scaling does not compound, runs without NTV are
!> untouched, and a non-positive reference amplitude is rejected.
program test_ntv_antenna_scaling
    use QLBalance_kinds, only: dp
    use grid_mod, only: torque_density_ntv, torque_density_ntv_unscaled, &
                        ntv_antenna_factor_exponent
    use wave_code_data, only: antenna_factor
    use transp_coeffs_mod, only: scale_ntv_torque_by_antenna_factor
    use control_mod, only: is_time_evolution_run
    implicit none

    real(dp), parameter :: unscaled(3) = [2.0_dp, -1.5_dp, 0.25_dp]
    character(len=32) :: mode

    call get_command_argument(1, mode)
    allocate (torque_density_ntv(3))

    if (trim(mode) == 'zero_reference') then
        allocate (torque_density_ntv_unscaled(3))
        torque_density_ntv_unscaled = unscaled
        antenna_factor = 1.0_dp
        call scale_ntv_torque_by_antenna_factor(0.0_dp)
        stop 'zero reference was accepted'
    end if

    ! Without NTV the torque stays as it is.
    torque_density_ntv = 7.0_dp
    antenna_factor = 0.3_dp
    call scale_ntv_torque_by_antenna_factor(1.0_dp)
    call check(all(torque_density_ntv == 7.0_dp), 'runs without NTV are untouched')

    allocate (torque_density_ntv_unscaled(3))
    torque_density_ntv_unscaled = unscaled
    ntv_antenna_factor_exponent = 2.0_dp

    antenna_factor = 0.5_dp
    call scale_ntv_torque_by_antenna_factor(2.0_dp)
    call check(all(abs(torque_density_ntv - unscaled*0.0625_dp) < 1.0e-15_dp), &
               '(antenna_factor/reference)**exponent scaling')
    call scale_ntv_torque_by_antenna_factor(2.0_dp)
    call check(all(abs(torque_density_ntv - unscaled*0.0625_dp) < 1.0e-15_dp), &
               'repeated scaling does not compound')

    ! The step's final antenna factor (e.g. from I_par_toroidal) replaces the
    ! value seen when NEO-RT ran.
    antenna_factor = 1.0_dp
    call scale_ntv_torque_by_antenna_factor(2.0_dp)
    call check(all(abs(torque_density_ntv - unscaled*0.25_dp) < 1.0e-15_dp), &
               'scaling follows the updated antenna factor')

    call check(is_time_evolution_run('TimeEvolution'), 'TimeEvolution is a time evolution')
    call check(is_time_evolution_run('TimeEvolutionNTV '), 'TimeEvolutionNTV is a time evolution')
    call check(.not. is_time_evolution_run('ParameterScan'), 'ParameterScan is not')

    print *, 'test_ntv_antenna_scaling passed'

contains

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            print *, 'FAILED: ', message
            error stop 1
        end if
    end subroutine check

end program test_ntv_antenna_scaling
