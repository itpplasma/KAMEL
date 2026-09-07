program test_loading_bar
    use iso_fortran_env, only: dp => real64
    use loading_bar_m, only: drawLoadingBar, updateLoadingBar, updateLoadingBarWithETA

    implicit none

    integer(kind=8) :: start_count, count_rate

    call drawLoadingBar(50.0_dp, 10)
    write (*, '(A)') ' BAR_TEST_PASS'

    call system_clock(start_count, count_rate)
    call updateLoadingBar(50, 100, 'Test')
    write (*, '(A)') ''
    call updateLoadingBarWithETA(50, 100, start_count, count_rate, 'Test')
    write (*, '(A)') ''
    call updateLoadingBarWithETA(50, 100, start_count, count_rate)
    write (*, '(A)') ''
    write (*, '(A)') ' UPDATE_TEST_PASS'

end program test_loading_bar
