program test_signed_resonances
    use KIM_kinds_m, only: dp
    use kim_resonances_m, only: r_res
    use config_m, only: type_of_run
    use setup_m, only: m_mode, n_mode, type_br_field
    use species_m, only: plasma
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    external :: kim_prepare_resonances
    integer :: failures
    character(32) :: argument

    failures = 0
    type_of_run = 'electrostatic_periodic'
    type_br_field = 1
    m_mode = -6
    n_mode = 2
    call get_command_argument(1, argument)
    call set_profile([2.0_dp, 2.5_dp, 3.5_dp, 4.0_dp])
    if (len_trim(argument) > 0) then
        select case (trim(argument))
        case ('nan_q')
            plasma%q(2) = ieee_value(0.0_dp, ieee_quiet_nan)
        case ('nan_r')
            plasma%r_grid(1) = ieee_value(0.0_dp, ieee_quiet_nan)
        case ('unordered_r')
            plasma%r_grid(2) = plasma%r_grid(1)
        case ('shape')
            plasma%grid_size = 3
        case ('missing_q')
            deallocate (plasma%q)
        case default
            error stop 'unknown signed-resonance rejection test'
        end select
        call kim_prepare_resonances()
        print *, 'Invalid resonance geometry was accepted'
        stop 0
    end if

    call expect_radius(25.0_dp, 'signed increasing positive q')
    m_mode = 6
    call expect_radius(0.0_dp, 'wrong sign has no resonance and clears previous radius')
    n_mode = -2
    call expect_radius(25.0_dp, 'negative n retains signed rational condition')
    n_mode = 2
    call set_profile([-2.0_dp, -2.5_dp, -3.5_dp, -4.0_dp])
    call expect_radius(25.0_dp, 'negative q and positive m locate physical resonance')
    m_mode = -6
    call expect_radius(0.0_dp, 'negative q rejects wrong signed mode')
    call set_profile([4.0_dp, 3.5_dp, 2.5_dp, 2.0_dp])
    call expect_radius(25.0_dp, 'descending q interpolates radius instead of returning qres')
    call set_profile([3.0_dp, 3.5_dp, 4.0_dp, 4.5_dp])
    call expect_radius(10.0_dp, 'exact first endpoint')
    call set_profile([1.0_dp, 2.0_dp, 2.5_dp, 3.0_dp])
    call expect_radius(40.0_dp, 'exact last endpoint')
    call set_profile([2.0_dp, 3.0_dp, 3.0_dp, 4.0_dp])
    call expect_radius(20.0_dp, 'plateau selects its innermost exact sample')
    call set_profile([4.0_dp, 2.0_dp, 4.0_dp, 2.0_dp])
    call expect_radius(15.0_dp, 'multiple crossings select innermost radius')
    type_br_field = 2
    call expect_radius(15.0_dp, 'periodic point-charge selector cannot fake resonance')
    type_of_run = 'electromagnetic'
    call expect_radius(20.0_dp, 'nonperiodic point-charge override remains available')
    type_of_run = 'electrostatic_periodic'
    type_br_field = 1
    m_mode = -12
    call expect_radius(0.0_dp, 'out-of-range rational value clears stale radius')
    m_mode = -6
    n_mode = 0
    call expect_radius(0.0_dp, 'zero toroidal mode has no rational resonance')
    if (failures /= 0) error stop 'signed resonance regression failed'
    print *, 'signed resonance regression passed'

contains

    subroutine set_profile(q)
        real(dp), intent(in) :: q(:)
        integer :: i
        if (allocated(plasma%q)) deallocate (plasma%q)
        if (allocated(plasma%r_grid)) deallocate (plasma%r_grid)
        plasma%grid_size = size(q)
        plasma%q = q
        plasma%r_grid = [(10.0_dp * real(i, dp), i=1, size(q))]
    end subroutine set_profile

    subroutine expect_radius(expected, label)
        real(dp), intent(in) :: expected
        character(*), intent(in) :: label
        r_res = 777.0_dp
        call kim_prepare_resonances()
        if (.not. ieee_is_finite(r_res) .or. abs(r_res - expected) > 1.0e-12_dp) then
            failures = failures + 1
            print *, 'FAIL: ', label, ' got ', r_res, ' expected ', expected
        end if
    end subroutine expect_radius
end program test_signed_resonances
