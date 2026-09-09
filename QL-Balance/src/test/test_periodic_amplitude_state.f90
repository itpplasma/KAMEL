program test_periodic_amplitude_state
    use QLBalance_kinds, only: dp
    use periodic_amplitude_state_m, only: periodic_amplitude_state_t, periodic_normalization_version
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    implicit none
    type(periodic_amplitude_state_t) :: state
    complex(dp), parameter :: initial(2) = [(1.0_dp, 0.2_dp), (2.0_dp, -0.3_dp)]
    complex(dp), parameter :: next(2) = [(3.0_dp, 1.0_dp), (4.0_dp, -1.0_dp)]
    complex(dp), parameter :: unit(2) = [(5.0_dp, 0.0_dp), (6.0_dp, 0.1_dp)]
    complex(dp), parameter :: residual(2) = [(0.1_dp, 0.2_dp), (0.3_dp, 0.4_dp)]
    complex(dp) :: values(2), units(2), residuals(2)
    real(dp) :: bad(2), target, relaxation
    integer :: stat, item, variant, failures
    character(32) :: argument

    failures = 0
    bad = [ieee_value(0.0_dp, ieee_quiet_nan), ieee_value(0.0_dp, ieee_positive_inf)]
    call get_command_argument(1, argument)
    if (len_trim(argument) > 0) then
        call invalid_operation(trim(argument))
        print *, 'Invalid state operation was accepted'
        stop 0
    end if

    call state%begin_trial(initial, unit, residual, [0, 0], 9.0_dp, 0.75_dp)
    call check(.not. state%initialized, 'fresh trial does not create accepted state')
    call check(state%trial_ready, 'fresh trial is staged')
    call check(.not. allocated(state%accepted), 'fresh accepted storage remains absent')
    call state%accept(stat)
    call check(stat == 0 .and. state%initialized, 'fresh valid trial can be accepted')
    call check(.not. state%trial_ready, 'accept consumes the pending trial')
    call check(all(state%accepted == initial), 'accepted values match staged values')

    call state%begin_trial(next, 2.0_dp * unit, 2.0_dp * residual, [0, 0], 10.0_dp, 0.5_dp)
    call state%reject()
    call check(all(state%trial == initial), 'reject restores amplitude')
    call check(all(state%trial_current_unit == unit), 'reject restores current')
    call check(all(state%trial_residual == residual), 'reject restores residual')
    call check(all(state%trial_status == 0), 'reject restores statuses')
    call check(state%target_current == 9.0_dp .and. state%relaxation == 0.75_dp, &
               'reject restores accepted controls')
    call check(.not. state%trial_ready, 'reject consumes the pending trial')

    do item = 1, 3
        call state%begin_trial(next, 2.0_dp * unit, 2.0_dp * residual, [0, item], 10.0_dp, 0.5_dp)
        call state%accept(stat)
        call check(stat /= 0, 'normalization guard prevents commit')
        call check_accepted()
    end do
    ! Missing optional metadata comes from accepted state, never a rejected
    ! candidate, even when the caller replaces a failed trial without reject().
    call state%begin_trial(next)
    call check(all(state%trial_current_unit == unit) .and. all(state%trial_status == 0), &
               'replacement trial inherits accepted metadata')
    call state%reject()

    do variant = 1, size(bad)
        do item = 1, 8
            values = next
            units = unit
            residuals = residual
            target = 10.0_dp
            relaxation = 0.5_dp
            select case (item)
            case (1)
                values(1) = cmplx(bad(variant), 0.0_dp, dp)
            case (2)
                values(1) = cmplx(0.0_dp, bad(variant), dp)
            case (3)
                units(1) = cmplx(bad(variant), 0.0_dp, dp)
            case (4)
                units(1) = cmplx(0.0_dp, bad(variant), dp)
            case (5)
                residuals(1) = cmplx(bad(variant), 0.0_dp, dp)
            case (6)
                residuals(1) = cmplx(0.0_dp, bad(variant), dp)
            case (7)
                target = bad(variant)
            case (8)
                relaxation = bad(variant)
            end select
            call state%begin_trial(values, units, residuals, [0, 0], target, relaxation)
            call state%accept(stat)
            call check(stat /= 0, 'every nonfinite trial component prevents commit')
            call check_accepted()
        end do
    end do
    do item = 0, 2
        call state%begin_trial(next, relaxation=real(item, dp))
        if (item == 1) cycle
        call state%accept(stat)
        call check(stat /= 0, 'out-of-range relaxation prevents commit')
        call check_accepted()
    end do
    call state%begin_trial(next, status=[-2, 0])
    call state%accept(stat)
    call check(stat /= 0, 'unknown status prevents commit')
    call check_accepted()
    call state%begin_trial(next, status=[-1, -1], target_current=-1.0_dp)
    call state%accept(stat)
    call check(stat == 0 .and. all(state%accepted == next), 'finite manual state is admissible')

    call state%reset()
    call check(.not. state%initialized .and. .not. state%trial_ready, 'reset clears lifecycle')
    call check(.not. allocated(state%accepted) .and. .not. allocated(state%trial) .and. &
               .not. allocated(state%accepted_current_unit) .and. &
               .not. allocated(state%trial_current_unit) .and. &
        .not. allocated(state%accepted_residual) .and. .not. allocated(state%trial_residual) .and. &
               .not. allocated(state%accepted_status) .and. .not. allocated(state%trial_status), &
               'reset releases every snapshot array')
    call check(state%target_current == 0.0_dp .and. state%relaxation == 1.0_dp, &
               'reset clears controls')
    call check(periodic_normalization_version == 2, 'residual convention version is two')
    call state%accept(stat)
    call check(stat /= 0 .and. .not. state%initialized, 'no pending trial cannot commit')
    call state%begin_trial(next, status=[3, 0])
    call state%accept(stat)
 call check(stat /= 0 .and. .not. state%initialized, 'failed first trial creates no accepted state')
    call state%initialize(initial, unit, residual, [0, 0], 9.0_dp, 0.75_dp)
    call check_accepted()
    call state%initialize(next(:1))
    call check(size(state%accepted) == 1, 'explicit initialization may start a new mode set')
    if (failures /= 0) error stop 'periodic amplitude state regression'
    print *, 'periodic amplitude state tests passed'

contains
    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            failures = failures + 1
            print *, 'FAIL: ', label
        end if
    end subroutine check

    subroutine check_accepted()
        call check(state%initialized .and. all(state%accepted == initial), &
                   'failed trial leaves accepted amplitude unchanged')
        call check(all(state%accepted_current_unit == unit), 'accepted current unchanged')
        call check(all(state%accepted_residual == residual), 'accepted residual unchanged')
        call check(all(state%accepted_status == 0), 'accepted statuses unchanged')
        call check(state%accepted_target_current == 9.0_dp .and. &
                   state%accepted_relaxation == 0.75_dp, 'accepted controls unchanged')
    end subroutine check_accepted

    subroutine invalid_operation(which)
        character(*), intent(in) :: which
        call state%initialize(initial)
        select case (which)
        case ('mode_count')
            call state%begin_trial(next(:1))
        case ('metadata_shape')
            call state%begin_trial(next, current_unit=unit(:1))
        case ('guard_no_stat')
            call state%begin_trial(next, status=[0, 3])
            call state%accept()
        case ('invalid_initialize')
            call state%initialize(next, status=[0, 2])
        case ('reject_fresh')
            call state%reset()
            call state%reject()
        case default
            error stop 'unknown amplitude-state rejection test'
        end select
    end subroutine invalid_operation
end program test_periodic_amplitude_state
