program test_profile_input
    !> Unit tests for coordinate-type detection in profile_input_m.
    !> Exercises the production classify_coordinate_type directly (no inline copy),
    !> so the COORD_THRESHOLD decision is guarded by this test.
    use KIM_kinds_m, only: dp
    use profile_input_m, only: classify_coordinate_type
    implicit none

    logical :: all_passed
    all_passed = .true.

    call check('r_eff well above threshold', classify_coordinate_type(45.0_dp),   'r_eff',     all_passed)
    call check('r_eff just above threshold', classify_coordinate_type(2.0001_dp), 'r_eff',     all_passed)
    call check('sqrt_psiN below threshold',  classify_coordinate_type(0.95_dp),   'sqrt_psiN', all_passed)
    call check('sqrt_psiN at threshold',     classify_coordinate_type(2.0_dp),    'sqrt_psiN', all_passed)

    if (all_passed) then
        print *, 'All tests PASSED'
        stop 0
    else
        print *, 'Some tests FAILED'
        stop 1
    end if

contains

    subroutine check(name, got, expected, passed)
        character(*), intent(in) :: name, got, expected
        logical, intent(inout) :: passed

        if (trim(got) == trim(expected)) then
            print *, 'PASS: ', name
        else
            print *, 'FAIL: ', name, ' expected=', trim(expected), ' got=', trim(got)
            passed = .false.
        end if
    end subroutine check

end program test_profile_input
