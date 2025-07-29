program test_performance
    use iso_fortran_env, only: real64, int32, int64
    use kilca_complex_m
    use kilca_stability_m
    use kilca_spline_m
    implicit none
    
    ! Test parameters
    integer, parameter :: WARMUP_RUNS = 10
    integer, parameter :: TEST_RUNS = 100
    integer, parameter :: SMALL_SIZE = 100
    integer, parameter :: MEDIUM_SIZE = 1000
    integer, parameter :: LARGE_SIZE = 10000
    
    ! Timing variables
    real :: start_time, end_time
    real(real64) :: elapsed_time, avg_time
    
    ! Test data
    complex(real64), allocatable :: z(:), z_result(:)
    complex(real64), allocatable :: A(:,:), b(:)
    real(real64), allocatable :: x(:), y(:), eval_points(:), result(:)
    type(spline_data_t) :: spline
    
    integer :: i, j, n, ierr
    real(real64) :: dummy
    
    print *, "===================================================="
    print *, "KiLCA Performance Testing - Custom Routines Only"
    print *, "===================================================="
    print *, ""
    
    ! Test 1: Complex number operations performance
    call test_complex_performance()
    
    ! Test 2: Stability analysis performance
    call test_stability_performance()
    
    ! Test 3: Spline evaluation performance
    call test_spline_performance()
    
    print *, ""
    print *, "===================================================="
    print *, "Performance testing complete"
    print *, "===================================================="
    
contains

    !---------------------------------------------------------------------------
    ! Test complex number operations performance
    !---------------------------------------------------------------------------
    subroutine test_complex_performance()
        integer :: run
        complex(real64) :: z1, z2
        real(real64) :: r
        
        print *, "1. Complex Number Operations Performance"
        print *, "----------------------------------------"
        
        ! Test 1.1: Complex logarithm of absolute value (custom implementation)
        n = LARGE_SIZE
        allocate(z(n))
        
        ! Initialize test data
        do i = 1, n
            z(i) = cmplx(real(i, real64) * 1.0e-10_real64, real(i, real64) * 1.0e-10_real64, real64)
        end do
        
        ! Warmup
        do run = 1, WARMUP_RUNS
            do i = 1, n
                r = cmplx_logabs(z(i))
            end do
        end do
        
        ! Timed runs
        call cpu_time(start_time)
        do run = 1, TEST_RUNS
            do i = 1, n
                r = cmplx_logabs(z(i))
            end do
        end do
        call cpu_time(end_time)
        
        avg_time = (end_time - start_time) / real(TEST_RUNS * n, real64)
        print '(A,ES12.3,A)', "   cmplx_logabs:           ", avg_time * 1.0e9_real64, " ns/operation"
        
        ! Test 1.2: Safe complex power (custom implementation)
        z1 = cmplx(0.0_real64, 0.0_real64, real64)
        z2 = cmplx(2.0_real64, 3.0_real64, real64)
        
        call cpu_time(start_time)
        do run = 1, TEST_RUNS * 1000
            z1 = cmplx_pow_safe(z1, z2)
        end do
        call cpu_time(end_time)
        
        avg_time = (end_time - start_time) / real(TEST_RUNS * 1000, real64)
        print '(A,ES12.3,A)', "   cmplx_pow_safe:         ", avg_time * 1.0e9_real64, " ns/operation"
        
        ! Test 1.3: Fast operations
        allocate(z_result(n))
        
        ! Fast magnitude squared
        call cpu_time(start_time)
        do run = 1, TEST_RUNS
            do i = 1, n
                dummy = cmplx_abs2_fast(z(i))
            end do
        end do
        call cpu_time(end_time)
        
        avg_time = (end_time - start_time) / real(TEST_RUNS * n, real64)
        print '(A,ES12.3,A)', "   cmplx_abs2_fast:        ", avg_time * 1.0e9_real64, " ns/operation"
        
        deallocate(z, z_result)
        
    end subroutine test_complex_performance
    
    !---------------------------------------------------------------------------
    ! Test stability analysis performance
    !---------------------------------------------------------------------------
    subroutine test_stability_performance()
        integer :: run, size
        real(real64) :: cond_num, rcond
        
        print *, ""
        print *, "2. Stability Analysis Performance"
        print *, "----------------------------------------"
        
        do size = 1, 3
            select case(size)
            case(1)
                n = 10
                print *, "   Small matrices (10x10):"
            case(2)
                n = 100
                print *, "   Medium matrices (100x100):"
            case(3)
                n = 500
                print *, "   Large matrices (500x500):"
            end select
            
            allocate(A(n,n))
            
            ! Initialize well-conditioned matrix
            A = cmplx(0.0_real64, 0.0_real64, real64)
            do i = 1, n
                A(i,i) = cmplx(real(i, real64), 0.0_real64, real64)
                if (i > 1) A(i,i-1) = cmplx(0.5_real64, 0.0_real64, real64)
                if (i < n) A(i,i+1) = cmplx(0.5_real64, 0.0_real64, real64)
            end do
            
            ! Test condition number estimation
            call cpu_time(start_time)
            do run = 1, TEST_RUNS
                cond_num = estimate_condition_number(A, ierr)
            end do
            call cpu_time(end_time)
            
            avg_time = (end_time - start_time) / real(TEST_RUNS, real64)
            print '(A,ES12.3,A)', "     Condition number:     ", avg_time * 1.0e3_real64, " ms/operation"
            
            ! Test reciprocal condition number
            call cpu_time(start_time)
            do run = 1, TEST_RUNS
                rcond = estimate_rcond(A, ierr)
            end do
            call cpu_time(end_time)
            
            avg_time = (end_time - start_time) / real(TEST_RUNS, real64)
            print '(A,ES12.3,A)', "     Reciprocal cond:      ", avg_time * 1.0e3_real64, " ms/operation"
            
            deallocate(A)
        end do
        
    end subroutine test_stability_performance
    
    !---------------------------------------------------------------------------
    ! Test spline evaluation performance
    !---------------------------------------------------------------------------
    subroutine test_spline_performance()
        integer :: run, npts
        
        print *, ""
        print *, "3. Spline Evaluation Performance"
        print *, "----------------------------------------"
        
        ! Create spline data
        n = 100  ! Number of spline nodes
        allocate(x(n), y(n))
        
        do i = 1, n
            x(i) = real(i-1, real64) / real(n-1, real64)
            y(i) = sin(10.0_real64 * x(i))
        end do
        
        ! Create spline
        call spline_create(spline, 3, SPLINE_NATURAL, n, x, ierr)
        call spline_calc_coefficients(spline, y, ierr)
        
        ! Test different evaluation sizes
        do i = 1, 3
            select case(i)
            case(1)
                npts = 100
                print *, "   Small evaluation (100 points):"
            case(2)
                npts = 1000
                print *, "   Medium evaluation (1000 points):"
            case(3)
                npts = 10000
                print *, "   Large evaluation (10000 points):"
            end select
            
            allocate(eval_points(npts), result(npts))
            
            ! Create evaluation points
            do j = 1, npts
                eval_points(j) = real(j-1, real64) / real(npts-1, real64)
            end do
            
            ! Warmup
            do run = 1, WARMUP_RUNS
                call spline_eval(spline, eval_points, result, ierr)
            end do
            
            ! Timed runs
            call cpu_time(start_time)
            do run = 1, TEST_RUNS
                call spline_eval(spline, eval_points, result, ierr)
            end do
            call cpu_time(end_time)
            
            avg_time = (end_time - start_time) / real(TEST_RUNS * npts, real64)
            print '(A,ES12.3,A)', "     Spline eval:          ", avg_time * 1.0e9_real64, " ns/point"
            
            ! Test derivative evaluation
            call cpu_time(start_time)
            do run = 1, TEST_RUNS
                call spline_eval_deriv(spline, eval_points, 1, result, ierr)
            end do
            call cpu_time(end_time)
            
            avg_time = (end_time - start_time) / real(TEST_RUNS * npts, real64)
            print '(A,ES12.3,A)', "     Spline deriv:         ", avg_time * 1.0e9_real64, " ns/point"
            
            deallocate(eval_points, result)
        end do
        
        call spline_destroy(spline)
        deallocate(x, y)
        
    end subroutine test_spline_performance

end program test_performance