program test_findIndex
    !> Unit tests for find_index_m%findClosestIndex: returns the index of the
    !> array element nearest a target value (used for grid/profile lookups).
    use KIM_kinds_m, only: dp
    use find_index_m, only: findClosestIndex
    implicit none

    logical :: all_passed
    real(dp) :: arr(4)
    integer :: idx

    all_passed = .true.
    arr = [1.0_dp, 5.0_dp, 10.0_dp, 20.0_dp]

    call findClosestIndex(arr, 6.0_dp, idx)
    call report('closest to 6.0 is 5.0 (idx 2)', idx == 2, all_passed)

    call findClosestIndex(arr, -3.0_dp, idx)
    call report('below range -> first element (idx 1)', idx == 1, all_passed)

    call findClosestIndex(arr, 100.0_dp, idx)
    call report('above range -> last element (idx 4)', idx == 4, all_passed)

    call findClosestIndex(arr, 10.0_dp, idx)
    call report('exact match (idx 3)', idx == 3, all_passed)

    if (all_passed) then
        print *, 'All findIndex tests PASSED'
        stop 0
    else
        print *, 'Some findIndex tests FAILED'
        stop 1
    end if

contains

    subroutine report(name, ok, passed)
        character(*), intent(in) :: name
        logical, intent(in) :: ok
        logical, intent(inout) :: passed
        if (ok) then
            print *, 'PASS: ', name
        else
            print *, 'FAIL: ', name
            passed = .false.
        end if
    end subroutine report

end program test_findIndex
