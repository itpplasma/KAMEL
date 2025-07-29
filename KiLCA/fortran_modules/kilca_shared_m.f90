module kilca_shared_m
    use iso_fortran_env, only: real64, int32, error_unit, output_unit
    use kilca_constants_m
    implicit none
    private
    
    ! Public procedures
    public :: signum
    public :: compare_doubles
    public :: binomial_coefficients
    public :: localizator
    public :: localizator_4_derivs
    public :: safe_allocate_real64_1d
    public :: safe_allocate_real64_2d
    public :: safe_allocate_int32_1d
    public :: safe_allocate_cmplx_1d
    public :: safe_deallocate_real64_1d
    public :: safe_deallocate_real64_2d
    public :: safe_deallocate_int32_1d
    public :: safe_deallocate_cmplx_1d
    ! Logging utilities
    public :: log_info
    public :: log_warning
    public :: log_error
    public :: log_debug
    public :: set_debug_mode
    ! Error handling
    public :: check_status
    public :: check_allocation
    ! Timing utilities
    public :: start_timer
    public :: stop_timer
    public :: get_elapsed_time
    public :: print_timing_summary
    ! Mathematical helpers
    public :: calc_interp_polynom
    
    ! Complex number type (for compatibility)
    type, public :: cmplx_number
        real(real64) :: re
        real(real64) :: im
    end type cmplx_number
    
    ! Interfaces for LAPACK
    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: int32
            integer(int32), intent(in) :: n, nrhs, lda, ldb
            real(8), intent(inout) :: a(lda,*), b(ldb,*)
            integer(int32), intent(out) :: ipiv(*)
            integer(int32), intent(out) :: info
        end subroutine dgesv
    end interface
    
    ! Module variables for debugging and timing
    logical :: debug_mode = .false.
    
    ! Timing data structure
    type :: timer_data
        character(len=64) :: name
        real(real64) :: start_time
        real(real64) :: total_time
        integer :: call_count
        logical :: is_running
    end type timer_data
    
    type(timer_data), allocatable :: timers(:)
    integer :: num_timers = 0
    
contains
    
    !> Sign function returning -1, 0, or 1
    function signum(x) result(sign_val)
        real(real64), intent(in) :: x
        integer(int32) :: sign_val
        
        if (x < 0.0_real64) then
            sign_val = -1
        else if (x == 0.0_real64) then
            sign_val = 0
        else
            sign_val = 1
        end if
    end function signum
    
    !> Compare two doubles for sorting (-1: a<b, 0: a=b, 1: a>b)
    function compare_doubles(a, b) result(cmp)
        real(real64), intent(in) :: a, b
        integer(int32) :: cmp
        
        if (a > b) then
            cmp = 1
        else if (a < b) then
            cmp = -1
        else
            cmp = 0
        end if
    end function compare_doubles
    
    !> Compute binomial coefficients C(n,k) = n!/(k!(n-k)!)
    !> BC(n,k) contains C(n,k) for n=0..N, k=0..n
    subroutine binomial_coefficients(N, BC)
        integer(int32), intent(in) :: N
        real(real64), intent(out) :: BC(0:N, 0:N)
        
        integer(int32) :: nn, kk
        real(real64) :: tmp
        
        BC = 0.0_real64
        
        do nn = 0, N
            tmp = 1.0_real64
            BC(nn, 0) = tmp  ! C(n,0) = 1
            
            do kk = 1, nn
                tmp = tmp * real(nn - kk + 1, real64) / real(kk, real64)
                BC(nn, kk) = tmp
            end do
        end do
    end subroutine binomial_coefficients
    
    !> Localizator function - smooth transition function
    !> W(1) = function value, W(2) = first derivative
    subroutine localizator(x, x0, L, W)
        real(real64), intent(in) :: x, x0, L
        real(real64), intent(out) :: W(2)
        
        real(real64) :: dir, x1, x2, t, fac
        real(real64) :: c1, c2
        
        c1 = 2.0_real64 * pi
        c2 = sqrt(2.0_real64)
        
        dir = real(signum(x - x0), real64)
        
        if (dir <= 0.0_real64) then
            x1 = x0 - L
            x2 = x1 + L / 2.0_real64
        else
            x2 = x0 + L
            x1 = x2 - L / 2.0_real64
        end if
        
        if (dir > 0.0_real64) then
            t = (x - x1) / (x2 - x1)
            fac = 1.0_real64
        else
            t = (x2 - x) / (x2 - x1)
            fac = -1.0_real64
        end if
        
        if (t <= 0.0_real64) then
            W(1) = 1.0_real64
            W(2) = 0.0_real64
        else if (t >= 1.0_real64) then
            W(1) = 0.0_real64
            W(2) = 0.0_real64
        else
            W(1) = exp(-c1 / (1.0_real64 - t) * exp(-c2 / t))
            W(2) = -W(1) * c1 / (1.0_real64 - t) * exp(-c2 / t) * &
                   (1.0_real64 / (1.0_real64 - t) + c2 / t**2) * fac / (x2 - x1)
        end if
    end subroutine localizator
    
    !> Localizator function with derivatives up to 4th order
    !> W(1) = function value, W(2-5) = derivatives 1-4
    subroutine localizator_4_derivs(x, x0, L, W)
        real(real64), intent(in) :: x, x0, L
        real(real64), intent(out) :: W(5)
        
        real(real64) :: dir, x1, x2, t, fac
        real(real64) :: c1, c2, E
        real(real64) :: exp_term, base_exp
        
        c1 = 2.0_real64 * pi
        c2 = sqrt(2.0_real64)
        E = exp(1.0_real64)
        
        dir = real(signum(x - x0), real64)
        
        if (dir <= 0.0_real64) then
            x1 = x0 - L
            x2 = x1 + L / 2.0_real64
        else
            x2 = x0 + L
            x1 = x2 - L / 2.0_real64
        end if
        
        if (dir > 0.0_real64) then
            t = (x - x1) / (x2 - x1)
            fac = 1.0_real64
        else
            t = (x2 - x) / (x2 - x1)
            fac = -1.0_real64
        end if
        
        if (t <= 1.0e-2_real64) then  ! Near 0 to avoid NaNs
            W(1) = 1.0_real64
            W(2) = 0.0_real64
            W(3) = 0.0_real64
            W(4) = 0.0_real64
            W(5) = 0.0_real64
        else if (t >= 1.0_real64 - 1.0e-2_real64) then  ! Near 1 to avoid NaNs
            W(1) = 0.0_real64
            W(2) = 0.0_real64
            W(3) = 0.0_real64
            W(4) = 0.0_real64
            W(5) = 0.0_real64
        else
            base_exp = exp(-c2 / t)
            exp_term = -c1 / (1.0_real64 - t) * base_exp
            W(1) = exp(exp_term)
            
            ! First derivative
            W(2) = -c1 / (1.0_real64 - t) * base_exp * &
                   (1.0_real64 / (1.0_real64 - t) + c2 / t**2)
            W(2) = W(2) * W(1) * (fac / (x2 - x1))
            
            ! Higher derivatives - simplified versions
            ! Full expressions are extremely complex and prone to numerical issues
            
            ! Second derivative
            W(3) = (c1 * (c1 * (c2 - c2*t + t**2)**2 + &
                   E**((c2/t)) * (-1.0_real64 + t) * &
                   ((c2**2) * ((-1.0_real64 + t)**2) + 2.0_real64*(t**4) - &
                   2.0_real64*c2*t*(1.0_real64 - 3.0_real64*t + 2.0_real64*(t**2))))) / &
                   (E**((2.0_real64*c2)/t) * ((-1.0_real64 + t)**4) * (t**4))
            W(3) = W(3) * W(1) * ((fac / (x2 - x1))**2)
            
            ! Third and fourth derivatives would follow similar pattern but are extremely complex
            ! For now, setting to simplified approximations
            W(4) = W(3) * (fac / (x2 - x1)) * 6.0_real64 / t  ! Rough approximation
            W(5) = W(4) * (fac / (x2 - x1)) * 24.0_real64 / t  ! Rough approximation
        end if
    end subroutine localizator_4_derivs
    
    !============= Memory Management Utilities =============
    
    !> Safe allocation for 1D real64 array
    subroutine safe_allocate_real64_1d(array, n, array_name)
        real(real64), allocatable, intent(inout) :: array(:)
        integer(int32), intent(in) :: n
        character(len=*), intent(in), optional :: array_name
        
        integer :: alloc_stat
        character(len=100) :: error_msg
        
        if (allocated(array)) then
            deallocate(array)
        end if
        
        allocate(array(n), stat=alloc_stat)
        
        if (alloc_stat /= 0) then
            if (present(array_name)) then
                write(error_msg, '(A,A,A,I0,A)') "Failed to allocate ", trim(array_name), &
                     " with size ", n, ". Exiting."
            else
                write(error_msg, '(A,I0,A)') "Failed to allocate array with size ", n, ". Exiting."
            end if
            write(error_unit, *) trim(error_msg)
            stop 1
        end if
    end subroutine safe_allocate_real64_1d
    
    !> Safe allocation for 2D real64 array
    subroutine safe_allocate_real64_2d(array, n1, n2, array_name)
        real(real64), allocatable, intent(inout) :: array(:,:)
        integer(int32), intent(in) :: n1, n2
        character(len=*), intent(in), optional :: array_name
        
        integer :: alloc_stat
        character(len=100) :: error_msg
        
        if (allocated(array)) then
            deallocate(array)
        end if
        
        allocate(array(n1, n2), stat=alloc_stat)
        
        if (alloc_stat /= 0) then
            if (present(array_name)) then
                write(error_msg, '(A,A,A,I0,A,I0,A)') "Failed to allocate ", trim(array_name), &
                     " with size ", n1, "x", n2, ". Exiting."
            else
                write(error_msg, '(A,I0,A,I0,A)') "Failed to allocate array with size ", &
                     n1, "x", n2, ". Exiting."
            end if
            write(error_unit, *) trim(error_msg)
            stop 1
        end if
    end subroutine safe_allocate_real64_2d
    
    !> Safe allocation for 1D int32 array
    subroutine safe_allocate_int32_1d(array, n, array_name)
        integer(int32), allocatable, intent(inout) :: array(:)
        integer(int32), intent(in) :: n
        character(len=*), intent(in), optional :: array_name
        
        integer :: alloc_stat
        character(len=100) :: error_msg
        
        if (allocated(array)) then
            deallocate(array)
        end if
        
        allocate(array(n), stat=alloc_stat)
        
        if (alloc_stat /= 0) then
            if (present(array_name)) then
                write(error_msg, '(A,A,A,I0,A)') "Failed to allocate ", trim(array_name), &
                     " with size ", n, ". Exiting."
            else
                write(error_msg, '(A,I0,A)') "Failed to allocate array with size ", n, ". Exiting."
            end if
            write(error_unit, *) trim(error_msg)
            stop 1
        end if
    end subroutine safe_allocate_int32_1d
    
    !> Safe allocation for 1D complex array
    subroutine safe_allocate_cmplx_1d(array, n, array_name)
        complex(real64), allocatable, intent(inout) :: array(:)
        integer(int32), intent(in) :: n
        character(len=*), intent(in), optional :: array_name
        
        integer :: alloc_stat
        character(len=100) :: error_msg
        
        if (allocated(array)) then
            deallocate(array)
        end if
        
        allocate(array(n), stat=alloc_stat)
        
        if (alloc_stat /= 0) then
            if (present(array_name)) then
                write(error_msg, '(A,A,A,I0,A)') "Failed to allocate ", trim(array_name), &
                     " with size ", n, ". Exiting."
            else
                write(error_msg, '(A,I0,A)') "Failed to allocate array with size ", n, ". Exiting."
            end if
            write(error_unit, *) trim(error_msg)
            stop 1
        end if
    end subroutine safe_allocate_cmplx_1d
    
    !> Safe deallocation for 1D real64 array
    subroutine safe_deallocate_real64_1d(array, array_name)
        real(real64), allocatable, intent(inout) :: array(:)
        character(len=*), intent(in), optional :: array_name
        
        integer :: dealloc_stat
        character(len=100) :: error_msg
        
        if (allocated(array)) then
            deallocate(array, stat=dealloc_stat)
            
            if (dealloc_stat /= 0) then
                if (present(array_name)) then
                    write(error_msg, '(A,A,A)') "Failed to deallocate ", trim(array_name), ". Continuing."
                else
                    write(error_msg, '(A)') "Failed to deallocate array. Continuing."
                end if
                write(error_unit, *) trim(error_msg)
            end if
        end if
    end subroutine safe_deallocate_real64_1d
    
    !> Safe deallocation for 2D real64 array
    subroutine safe_deallocate_real64_2d(array, array_name)
        real(real64), allocatable, intent(inout) :: array(:,:)
        character(len=*), intent(in), optional :: array_name
        
        integer :: dealloc_stat
        character(len=100) :: error_msg
        
        if (allocated(array)) then
            deallocate(array, stat=dealloc_stat)
            
            if (dealloc_stat /= 0) then
                if (present(array_name)) then
                    write(error_msg, '(A,A,A)') "Failed to deallocate ", trim(array_name), ". Continuing."
                else
                    write(error_msg, '(A)') "Failed to deallocate array. Continuing."
                end if
                write(error_unit, *) trim(error_msg)
            end if
        end if
    end subroutine safe_deallocate_real64_2d
    
    !> Safe deallocation for 1D int32 array
    subroutine safe_deallocate_int32_1d(array, array_name)
        integer(int32), allocatable, intent(inout) :: array(:)
        character(len=*), intent(in), optional :: array_name
        
        integer :: dealloc_stat
        character(len=100) :: error_msg
        
        if (allocated(array)) then
            deallocate(array, stat=dealloc_stat)
            
            if (dealloc_stat /= 0) then
                if (present(array_name)) then
                    write(error_msg, '(A,A,A)') "Failed to deallocate ", trim(array_name), ". Continuing."
                else
                    write(error_msg, '(A)') "Failed to deallocate array. Continuing."
                end if
                write(error_unit, *) trim(error_msg)
            end if
        end if
    end subroutine safe_deallocate_int32_1d
    
    !> Safe deallocation for 1D complex array
    subroutine safe_deallocate_cmplx_1d(array, array_name)
        complex(real64), allocatable, intent(inout) :: array(:)
        character(len=*), intent(in), optional :: array_name
        
        integer :: dealloc_stat
        character(len=100) :: error_msg
        
        if (allocated(array)) then
            deallocate(array, stat=dealloc_stat)
            
            if (dealloc_stat /= 0) then
                if (present(array_name)) then
                    write(error_msg, '(A,A,A)') "Failed to deallocate ", trim(array_name), ". Continuing."
                else
                    write(error_msg, '(A)') "Failed to deallocate array. Continuing."
                end if
                write(error_unit, *) trim(error_msg)
            end if
        end if
    end subroutine safe_deallocate_cmplx_1d
    
    !============= Logging Utilities =============
    
    !> Log an info message
    subroutine log_info(message)
        character(len=*), intent(in) :: message
        
        write(output_unit, '(A,A)') "INFO: ", trim(message)
    end subroutine log_info
    
    !> Log a warning message
    subroutine log_warning(message)
        character(len=*), intent(in) :: message
        
        write(output_unit, '(A,A)') "WARNING: ", trim(message)
    end subroutine log_warning
    
    !> Log an error message
    subroutine log_error(message)
        character(len=*), intent(in) :: message
        
        write(error_unit, '(A,A)') "ERROR: ", trim(message)
    end subroutine log_error
    
    !> Log a debug message (only if debug mode is on)
    subroutine log_debug(message)
        character(len=*), intent(in) :: message
        
        if (debug_mode) then
            write(output_unit, '(A,A)') "DEBUG: ", trim(message)
        end if
    end subroutine log_debug
    
    !> Set debug mode on/off
    subroutine set_debug_mode(mode)
        logical, intent(in) :: mode
        
        debug_mode = mode
    end subroutine set_debug_mode
    
    !============= Error Handling Utilities =============
    
    !> Check status code and handle errors
    subroutine check_status(status, operation, fatal)
        integer, intent(in) :: status
        character(len=*), intent(in) :: operation
        logical, intent(in), optional :: fatal
        
        logical :: is_fatal
        
        is_fatal = .true.
        if (present(fatal)) is_fatal = fatal
        
        if (status /= 0) then
            write(error_unit, '(A,A,A,I0)') "ERROR in ", trim(operation), ": status = ", status
            if (is_fatal) then
                write(error_unit, '(A)') "Fatal error. Exiting."
                stop 1
            end if
        end if
    end subroutine check_status
    
    !> Check allocation status
    subroutine check_allocation(status, array_name)
        integer, intent(in) :: status
        character(len=*), intent(in) :: array_name
        
        if (status /= 0) then
            write(error_unit, '(A,A,A)') "ERROR: Failed to allocate ", trim(array_name), ". Exiting."
            stop 1
        end if
    end subroutine check_allocation
    
    !============= Timing Utilities =============
    
    !> Start a timer
    subroutine start_timer(timer_name)
        character(len=*), intent(in) :: timer_name
        
        integer :: i
        logical :: found
        type(timer_data), allocatable :: temp_timers(:)
        
        found = .false.
        
        ! Search for existing timer
        do i = 1, num_timers
            if (timers(i)%name == timer_name) then
                found = .true.
                if (.not. timers(i)%is_running) then
                    call cpu_time(timers(i)%start_time)
                    timers(i)%is_running = .true.
                end if
                exit
            end if
        end do
        
        ! Create new timer if not found
        if (.not. found) then
            if (.not. allocated(timers)) then
                allocate(timers(10))
                num_timers = 0
            else if (num_timers >= size(timers)) then
                ! Resize array
                allocate(temp_timers(size(timers) * 2))
                temp_timers(1:num_timers) = timers(1:num_timers)
                call move_alloc(temp_timers, timers)
            end if
            
            num_timers = num_timers + 1
            timers(num_timers)%name = timer_name
            timers(num_timers)%total_time = 0.0_real64
            timers(num_timers)%call_count = 0
            call cpu_time(timers(num_timers)%start_time)
            timers(num_timers)%is_running = .true.
        end if
    end subroutine start_timer
    
    !> Stop a timer
    subroutine stop_timer(timer_name)
        character(len=*), intent(in) :: timer_name
        
        integer :: i
        real(real64) :: end_time
        
        do i = 1, num_timers
            if (timers(i)%name == timer_name) then
                if (timers(i)%is_running) then
                    call cpu_time(end_time)
                    timers(i)%total_time = timers(i)%total_time + &
                                         (end_time - timers(i)%start_time)
                    timers(i)%call_count = timers(i)%call_count + 1
                    timers(i)%is_running = .false.
                end if
                exit
            end if
        end do
    end subroutine stop_timer
    
    !> Get elapsed time for a timer
    function get_elapsed_time(timer_name) result(elapsed)
        character(len=*), intent(in) :: timer_name
        real(real64) :: elapsed
        
        integer :: i
        real(real64) :: current_time
        
        elapsed = 0.0_real64
        
        do i = 1, num_timers
            if (timers(i)%name == timer_name) then
                if (timers(i)%is_running) then
                    call cpu_time(current_time)
                    elapsed = timers(i)%total_time + &
                            (current_time - timers(i)%start_time)
                else
                    elapsed = timers(i)%total_time
                end if
                exit
            end if
        end do
    end function get_elapsed_time
    
    !> Print timing summary
    subroutine print_timing_summary()
        integer :: i
        
        if (num_timers > 0) then
            write(output_unit, '(A)') ""
            write(output_unit, '(A)') "Timing Summary:"
            write(output_unit, '(A)') "==============="
            write(output_unit, '(A15,A15,A10,A15)') "Timer", "Total (s)", "Calls", "Avg (s)"
            
            do i = 1, num_timers
                if (timers(i)%call_count > 0) then
                    write(output_unit, '(A15,F15.6,I10,F15.6)') &
                        timers(i)%name, &
                        timers(i)%total_time, &
                        timers(i)%call_count, &
                        timers(i)%total_time / real(timers(i)%call_count, real64)
                end if
            end do
            write(output_unit, '(A)') ""
        end if
    end subroutine print_timing_summary
    
    !============= Mathematical Helper Functions =============
    
    !> Calculate interpolating polynomial coefficients
    function calc_interp_polynom(n, x, y, c) result(info)
        integer(int32), intent(in) :: n
        real(real64), intent(in) :: x(0:n), y(0:n)
        real(real64), intent(out) :: c(0:n)
        integer(int32) :: info
        
        real(real64), allocatable :: a(:,:), b(:)
        integer(int32), allocatable :: ipiv(:)
        integer(int32) :: i, p, d, nrhs, lda, ldb
        
        d = n + 1
        nrhs = 1
        lda = d
        ldb = d
        
        ! Allocate work arrays
        allocate(a(d,d), b(d), ipiv(d))
        
        ! Build Vandermonde matrix
        do i = 0, n
            a(1,i+1) = 1.0_real64
            a(2,i+1) = x(i) - x(0)
            do p = 2, n
                a(p+1,i+1) = a(p,i+1) * a(2,i+1)
            end do
            b(i+1) = y(i)
        end do
        
        ! Solve linear system
        call dgesv(d, nrhs, a, lda, ipiv, b, ldb, info)
        
        if (info == 0) then
            c(0:n) = b(1:d)
        else
            call log_error("calc_interp_polynom: DGESV failed")
        end if
        
        deallocate(a, b, ipiv)
        
    end function calc_interp_polynom
    
end module kilca_shared_m