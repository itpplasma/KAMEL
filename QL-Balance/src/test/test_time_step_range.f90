program test_time_step_range
    use time_evolution, only: Nstorage, periodic_completed_steps, time_step_range

    implicit none

    integer :: first_step, last_step, num_failed

    num_failed = 0

    ! A fresh run covers 1..Nstorage.
    periodic_completed_steps = 0
    Nstorage = 5
    call time_step_range(first_step, last_step)
    call assert_true(first_step == 1 .and. last_step == 5, "fresh run range")

    ! A restart from a periodic checkpoint continues after the completed steps.
    periodic_completed_steps = 3
    call time_step_range(first_step, last_step)
    call assert_true(first_step == 4 .and. last_step == 8, "restart range")

    if (num_failed > 0) then
        print '(A,I0)', "FAILED assertions: ", num_failed
        stop 1
    end if
    print *, "ALL TESTS PASSED"

contains

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label

        if (.not. condition) then
            print '(A,A)', "FAIL: ", label
            num_failed = num_failed + 1
        end if
    end subroutine assert_true

end program test_time_step_range
