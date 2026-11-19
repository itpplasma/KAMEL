!> The Ipar history of TimeEvolution must start finite (zero) and take the
!> integrated parallel current that get_dql selected for the antenna factor.
program test_ipar_history
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use QLBalance_kinds, only: dp
    use grid_mod, only: Ipar
    use time_evolution, only: alloc_Br_Dqle_for_timeevol, record_parallel_current_history, &
                              Ipar_time, Nstorage, time_ind
    implicit none

    integer :: i
    complex(dp) :: expected(3)

    Nstorage = 3
    call alloc_Br_Dqle_for_timeevol()

    call check(size(Ipar_time) == 3, "history has one entry per step")
    do i = 1, 3
        call check(ieee_is_finite(real(Ipar_time(i))) .and. ieee_is_finite(aimag(Ipar_time(i))), &
                   "fresh history is finite")
        call check(abs(Ipar_time(i)) == 0.0_dp, "fresh history is zero")
    end do

    expected = [cmplx(1.5e3_dp, -2.0e2_dp, dp), cmplx(-4.0e-7_dp, 3.0e-7_dp, dp), &
                cmplx(0.0_dp, 9.0e5_dp, dp)]
    do time_ind = 1, 3
        Ipar = expected(time_ind)
        call record_parallel_current_history()
    end do

    do i = 1, 3
        call check(Ipar_time(i) == expected(i), "history equals the selected Ipar")
    end do

    print *, "test_ipar_history passed"

contains

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            print *, "FAILED: ", message
            error stop 1
        end if
    end subroutine check

end program test_ipar_history
