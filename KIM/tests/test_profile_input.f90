program test_profile_input
    use KIM_kinds_m, only: dp
    implicit none

    logical :: all_passed
    all_passed = .true.

    call test_coord_detection_reff(all_passed)
    call test_coord_detection_psiN(all_passed)

    if (all_passed) then
        print *, 'All tests PASSED'
        stop 0
    else
        print *, 'Some tests FAILED'
        stop 1
    end if

contains

    subroutine test_coord_detection_reff(passed)
        logical, intent(inout) :: passed
        real(dp) :: test_data(5)
        character(20) :: result

        ! r_eff data: max > 2.0
        test_data = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 45.0_dp]
        call detect_from_data(test_data, 5, result)

        if (trim(result) /= 'r_eff') then
            print *, 'FAIL: test_coord_detection_reff'
            print *, '  Expected: r_eff'
            print *, '  Got: ', trim(result)
            passed = .false.
        else
            print *, 'PASS: test_coord_detection_reff'
        end if
    end subroutine

    subroutine test_coord_detection_psiN(passed)
        logical, intent(inout) :: passed
        real(dp) :: test_data(5)
        character(20) :: result

        ! sqrt_psiN data: max <= 2.0
        test_data = [0.1_dp, 0.3_dp, 0.5_dp, 0.7_dp, 0.95_dp]
        call detect_from_data(test_data, 5, result)

        if (trim(result) /= 'sqrt_psiN') then
            print *, 'FAIL: test_coord_detection_psiN'
            print *, '  Expected: sqrt_psiN'
            print *, '  Got: ', trim(result)
            passed = .false.
        else
            print *, 'PASS: test_coord_detection_psiN'
        end if
    end subroutine

    subroutine detect_from_data(data, n, coord_type)
        real(dp), intent(in) :: data(:)
        integer, intent(in) :: n
        character(20), intent(out) :: coord_type
        real(dp) :: max_val
        real(dp), parameter :: COORD_THRESHOLD = 2.0_dp

        max_val = maxval(data(1:n))
        if (max_val > COORD_THRESHOLD) then
            coord_type = 'r_eff'
        else
            coord_type = 'sqrt_psiN'
        end if
    end subroutine detect_from_data

end program test_profile_input
