program test_kilca_complex_performance
    use iso_fortran_env, only: real64, int32
    use kilca_complex_m
    use kilca_constants_m, only: pi
    implicit none
    
    logical :: test_passed
    complex(real64) :: z1, z2, z3, result
    complex(real64), allocatable :: array1(:), array2(:), result_array(:)
    real(real64) :: start_time, end_time, elapsed_time
    integer :: i, n_ops, n_array
    real(real64) :: tol
    
    test_passed = .true.
    tol = epsilon(1.0_real64) * 100.0_real64
    n_ops = 100000  ! Number of operations for timing
    n_array = 10000  ! Size of arrays for bulk operations
    
    print *, "Testing complex number performance optimizations..."
    print *, ""
    
    ! Test 1: Fast magnitude squared calculation
    print *, "Testing fast magnitude squared calculation..."
    
    z1 = cmplx(3.0_real64, 4.0_real64, real64)
    
    ! Test cmplx_abs2_fast vs standard abs(z)**2
    call cpu_time(start_time)
    do i = 1, n_ops
        result = cmplx_abs2_fast(z1)
    end do
    call cpu_time(end_time)
    elapsed_time = end_time - start_time
    
    if (abs(cmplx_abs2_fast(z1) - 25.0_real64) > tol) then
        print *, "FAIL: cmplx_abs2_fast incorrect result"
        test_passed = .false.
    else
        print *, "PASS: cmplx_abs2_fast correct result"
        print *, "      Time for ", n_ops, " operations:", elapsed_time, " seconds"
    end if
    
    ! Test 2: Fast division by real number
    print *, ""
    print *, "Testing fast division by real number..."
    
    z1 = cmplx(6.0_real64, 8.0_real64, real64)
    z2 = cmplx_div_real_fast(z1, 2.0_real64)
    z3 = cmplx(3.0_real64, 4.0_real64, real64)
    
    if (abs(z2 - z3) > tol) then
        print *, "FAIL: cmplx_div_real_fast incorrect result"
        test_passed = .false.
    else
        print *, "PASS: cmplx_div_real_fast correct result"
    end if
    
    ! Test 3: Fast multiplication by imaginary unit
    print *, ""
    print *, "Testing fast multiplication by imaginary unit..."
    
    z1 = cmplx(3.0_real64, 4.0_real64, real64)
    z2 = cmplx_mult_i_fast(z1)
    z3 = cmplx(-4.0_real64, 3.0_real64, real64)  ! i * (3+4i) = -4+3i
    
    if (abs(z2 - z3) > tol) then
        print *, "FAIL: cmplx_mult_i_fast incorrect result"
        test_passed = .false.
    else
        print *, "PASS: cmplx_mult_i_fast correct result"
    end if
    
    ! Test 4: Fast conjugate multiplication for |z|^2
    print *, ""
    print *, "Testing fast conjugate multiplication..."
    
    z1 = cmplx(3.0_real64, 4.0_real64, real64)
    result = cmplx_mult_conj_fast(z1, z1)
    
    if (abs(real(result, real64) - 25.0_real64) > tol .or. abs(aimag(result)) > tol) then
        print *, "FAIL: cmplx_mult_conj_fast incorrect result"
        test_passed = .false.
    else
        print *, "PASS: cmplx_mult_conj_fast correct result"
    end if
    
    ! Test 5: Bulk array operations
    print *, ""
    print *, "Testing bulk array operations..."
    
    allocate(array1(n_array), array2(n_array), result_array(n_array))
    
    ! Initialize arrays
    do i = 1, n_array
        array1(i) = cmplx(real(i, real64), real(i+1, real64), real64)
        array2(i) = cmplx(real(i+2, real64), real(i+3, real64), real64)
    end do
    
    ! Test bulk addition
    call cpu_time(start_time)
    call cmplx_array_add_fast(array1, array2, result_array)
    call cpu_time(end_time)
    elapsed_time = end_time - start_time
    
    ! Verify first element
    z3 = array1(1) + array2(1)
    if (abs(result_array(1) - z3) > tol) then
        print *, "FAIL: cmplx_array_add_fast incorrect result"
        test_passed = .false.
    else
        print *, "PASS: cmplx_array_add_fast correct result"
        print *, "      Time for array size ", n_array, ":", elapsed_time, " seconds"
    end if
    
    ! Test 6: Optimized power functions for common exponents
    print *, ""
    print *, "Testing optimized power functions..."
    
    z1 = cmplx(2.0_real64, 1.0_real64, real64)
    
    ! Test z^2 optimization
    z2 = cmplx_square_fast(z1)
    z3 = z1 * z1
    if (abs(z2 - z3) > tol) then
        print *, "FAIL: cmplx_square_fast incorrect result"
        test_passed = .false.
    else
        print *, "PASS: cmplx_square_fast correct result"
    end if
    
    ! Test z^3 optimization
    z2 = cmplx_cube_fast(z1)
    z3 = z1 * z1 * z1
    if (abs(z2 - z3) > tol) then
        print *, "FAIL: cmplx_cube_fast incorrect result"
        test_passed = .false.
    else
        print *, "PASS: cmplx_cube_fast correct result"
    end if
    
    ! Test 7: Cached trigonometric values for common angles
    print *, ""
    print *, "Testing cached trigonometric values..."
    
    ! Test exp(i*pi/2) = i
    z2 = cmplx_exp_i_pi_half()
    z3 = cmplx_I
    if (abs(z2 - z3) > tol) then
        print *, "FAIL: cmplx_exp_i_pi_half incorrect result"
        test_passed = .false.
    else
        print *, "PASS: cmplx_exp_i_pi_half correct result"
    end if
    
    ! Test exp(i*pi) = -1
    z2 = cmplx_exp_i_pi()
    z3 = -cmplx_E
    if (abs(z2 - z3) > tol) then
        print *, "FAIL: cmplx_exp_i_pi incorrect result"
        test_passed = .false.
    else
        print *, "PASS: cmplx_exp_i_pi correct result"
    end if
    
    ! Test 8: Memory-efficient operations
    print *, ""
    print *, "Testing memory-efficient operations..."
    
    z1 = cmplx(3.0_real64, 4.0_real64, real64)
    z2 = cmplx(1.0_real64, 2.0_real64, real64)
    
    ! Test in-place operations
    call cmplx_add_inplace(z1, z2)  ! z1 = z1 + z2
    z3 = cmplx(4.0_real64, 6.0_real64, real64)
    if (abs(z1 - z3) > tol) then
        print *, "FAIL: cmplx_add_inplace incorrect result"
        test_passed = .false.
    else
        print *, "PASS: cmplx_add_inplace correct result"
    end if
    
    deallocate(array1, array2, result_array)
    
    print *, ""
    if (test_passed) then
        print *, "All complex number performance optimization tests PASSED!"
    else
        print *, "Some complex number performance optimization tests FAILED!"
        stop 1
    end if
    
end program test_kilca_complex_performance