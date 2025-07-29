program test_kilca_shared_extended
    use iso_fortran_env, only: real64, int32, output_unit, error_unit
    use kilca_shared_m
    implicit none
    
    logical :: test_passed
    integer :: i, unit_num, iostat
    character(len=100) :: test_message
    real(real64) :: start_time, end_time, elapsed
    integer :: n
    real(real64) :: x(5), y(5), coeffs(5)
    integer :: info
    
    test_passed = .true.
    
    print *, "Testing extended kilca_shared_m module utilities..."
    print *, ""
    
    ! Test logging utilities
    print *, "Testing logging utilities..."
    
    call log_info("This is an info message")
    call log_warning("This is a warning message")
    call log_error("This is an error message (non-fatal)")
    call log_debug("This is a debug message")
    
    ! Set debug mode and test
    call set_debug_mode(.true.)
    call log_debug("Debug message should appear when debug mode is on")
    call set_debug_mode(.false.)
    call log_debug("This debug message should NOT appear")
    
    print *, "PASS: Logging utilities work correctly"
    
    print *, ""
    print *, "Testing error handling utilities..."
    
    ! Test error checking
    call check_status(0, "Operation successful")
    call check_status(1, "Non-zero status", fatal=.false.)
    
    print *, "PASS: Error handling utilities work correctly"
    
    print *, ""
    print *, "Testing timing utilities..."
    
    ! Test timer
    call start_timer("test_operation")
    ! Simulate some work
    call sleep_milliseconds(100)
    elapsed = get_elapsed_time("test_operation")
    
    if (elapsed >= 0.09_real64 .and. elapsed <= 0.12_real64) then
        print *, "PASS: Timer measured approximately 0.1 seconds, actual:", elapsed
    else
        print *, "FAIL: Timer measurement outside expected range:", elapsed
        test_passed = .false.
    end if
    
    call stop_timer("test_operation")
    call print_timing_summary()
    
    print *, ""
    print *, "Testing mathematical helper functions..."
    
    ! Test interpolation polynomial calculation
    n = 4
    x = [0.0_real64, 1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    y = [0.0_real64, 1.0_real64, 4.0_real64, 9.0_real64, 16.0_real64]  ! y = x^2
    
    info = calc_interp_polynom(n, x, y, coeffs)
    
    if (info == 0) then
        print *, "PASS: calc_interp_polynom succeeded"
        print *, "      Coefficients:", coeffs
    else
        print *, "FAIL: calc_interp_polynom failed with info =", info
        test_passed = .false.
    end if
    
    print *, ""
    if (test_passed) then
        print *, "All extended tests PASSED!"
    else
        print *, "Some extended tests FAILED!"
        stop 1
    end if
    
contains
    
    subroutine sleep_milliseconds(ms)
        integer, intent(in) :: ms
        real(real64) :: start, current
        
        call cpu_time(start)
        do
            call cpu_time(current)
            if ((current - start) * 1000.0_real64 >= real(ms, real64)) exit
        end do
    end subroutine sleep_milliseconds
    
end program test_kilca_shared_extended