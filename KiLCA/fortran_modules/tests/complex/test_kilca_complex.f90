program test_kilca_complex
    use iso_fortran_env, only: real64, int32
    use kilca_complex_m
    use kilca_constants_m, only: pi
    implicit none
    
    logical :: test_passed
    complex(real64) :: z1, z2, z3, result
    real(real64) :: r1, r2, tol
    complex(real64), allocatable :: z_array(:), z_matrix(:,:)
    integer :: i, j, n
    
    test_passed = .true.
    tol = epsilon(1.0_real64) * 10.0_real64
    
    print *, "Testing kilca_complex_m module..."
    print *, ""
    
    ! Test complex constants
    print *, "Testing complex constants..."
    
    if (abs(cmplx_O) > tol) then
        print *, "FAIL: cmplx_O should be (0,0)"
        test_passed = .false.
    else
        print *, "PASS: cmplx_O = ", cmplx_O
    end if
    
    if (abs(cmplx_E - cmplx(1.0_real64, 0.0_real64, real64)) > tol) then
        print *, "FAIL: cmplx_E should be (1,0)"
        test_passed = .false.
    else
        print *, "PASS: cmplx_E = ", cmplx_E
    end if
    
    if (abs(cmplx_I - cmplx(0.0_real64, 1.0_real64, real64)) > tol) then
        print *, "FAIL: cmplx_I should be (0,1)"
        test_passed = .false.
    else
        print *, "PASS: cmplx_I = ", cmplx_I
    end if
    
    print *, ""
    print *, "Testing complex creation functions..."
    
    ! Test complex creation
    z1 = cmplx_make(3.0_real64, 4.0_real64)
    if (abs(z1 - cmplx(3.0_real64, 4.0_real64, real64)) > tol) then
        print *, "FAIL: cmplx_make(3,4)"
        test_passed = .false.
    else
        print *, "PASS: cmplx_make(3,4) = ", z1
    end if
    
    ! Test polar form
    r1 = 2.0_real64
    r2 = pi / 4.0_real64
    z1 = cmplx_polar(r1, r2)
    z2 = cmplx(r1 * cos(r2), r1 * sin(r2), real64)
    if (abs(z1 - z2) > tol) then
        print *, "FAIL: cmplx_polar"
        test_passed = .false.
    else
        print *, "PASS: cmplx_polar(2, pi/4) = ", z1
    end if
    
    print *, ""
    print *, "Testing complex utility functions..."
    
    ! Test real and imaginary parts
    z1 = cmplx(3.0_real64, 4.0_real64, real64)
    r1 = cmplx_real(z1)
    r2 = cmplx_imag(z1)
    if (abs(r1 - 3.0_real64) > tol .or. abs(r2 - 4.0_real64) > tol) then
        print *, "FAIL: cmplx_real/cmplx_imag"
        test_passed = .false.
    else
        print *, "PASS: cmplx_real = ", r1, ", cmplx_imag = ", r2
    end if
    
    ! Test absolute value
    r1 = cmplx_abs(z1)
    if (abs(r1 - 5.0_real64) > tol) then
        print *, "FAIL: cmplx_abs of (3,4) should be 5"
        test_passed = .false.
    else
        print *, "PASS: cmplx_abs((3,4)) = ", r1
    end if
    
    ! Test argument
    z1 = cmplx(1.0_real64, 1.0_real64, real64)
    r1 = cmplx_arg(z1)
    if (abs(r1 - pi/4.0_real64) > tol) then
        print *, "FAIL: cmplx_arg of (1,1) should be pi/4"
        test_passed = .false.
    else
        print *, "PASS: cmplx_arg((1,1)) = ", r1
    end if
    
    ! Test conjugate
    z1 = cmplx(3.0_real64, 4.0_real64, real64)
    z2 = cmplx_conj(z1)
    if (abs(z2 - cmplx(3.0_real64, -4.0_real64, real64)) > tol) then
        print *, "FAIL: cmplx_conj"
        test_passed = .false.
    else
        print *, "PASS: cmplx_conj((3,4)) = ", z2
    end if
    
    print *, ""
    print *, "Testing complex arithmetic operations..."
    
    ! Test addition
    z1 = cmplx(1.0_real64, 2.0_real64, real64)
    z2 = cmplx(3.0_real64, 4.0_real64, real64)
    z3 = cmplx_add(z1, z2)
    if (abs(z3 - cmplx(4.0_real64, 6.0_real64, real64)) > tol) then
        print *, "FAIL: cmplx_add"
        test_passed = .false.
    else
        print *, "PASS: cmplx_add((1,2), (3,4)) = ", z3
    end if
    
    ! Test multiplication
    z3 = cmplx_mult(z1, z2)
    result = cmplx(-5.0_real64, 10.0_real64, real64)
    if (abs(z3 - result) > tol) then
        print *, "FAIL: cmplx_mult"
        test_passed = .false.
    else
        print *, "PASS: cmplx_mult((1,2), (3,4)) = ", z3
    end if
    
    ! Test division
    z3 = cmplx_div(z1, z2)
    result = cmplx(0.44_real64, 0.08_real64, real64)
    if (abs(z3 - result) > tol * 100.0_real64) then
        print *, "FAIL: cmplx_div"
        test_passed = .false.
    else
        print *, "PASS: cmplx_div((1,2), (3,4)) = ", z3
    end if
    
    print *, ""
    print *, "Testing complex transcendental functions..."
    
    ! Test exponential
    z1 = cmplx(0.0_real64, pi, real64)
    z2 = cmplx_exp(z1)
    if (abs(z2 - cmplx(-1.0_real64, 0.0_real64, real64)) > tol * 100.0_real64) then
        print *, "FAIL: cmplx_exp(i*pi) should be -1"
        test_passed = .false.
    else
        print *, "PASS: cmplx_exp(i*pi) = ", z2
    end if
    
    ! Test logarithm
    z1 = cmplx(-1.0_real64, 0.0_real64, real64)
    z2 = cmplx_log(z1)
    if (abs(cmplx_imag(z2) - pi) > tol * 100.0_real64) then
        print *, "FAIL: cmplx_log(-1) should have imaginary part pi"
        test_passed = .false.
    else
        print *, "PASS: cmplx_log(-1) = ", z2
    end if
    
    ! Test power
    z1 = cmplx_I
    z2 = cmplx_pow_real(z1, 2.0_real64)
    if (abs(z2 - cmplx(-1.0_real64, 0.0_real64, real64)) > tol * 10.0_real64) then
        print *, "FAIL: i^2 should be -1"
        test_passed = .false.
    else
        print *, "PASS: i^2 = ", z2
    end if
    
    print *, ""
    print *, "Testing complex array operations..."
    
    ! Test array allocation and operations
    n = 5
    call cmplx_allocate_1d(z_array, n)
    
    if (.not. allocated(z_array)) then
        print *, "FAIL: cmplx_allocate_1d failed"
        test_passed = .false.
    else if (size(z_array) /= n) then
        print *, "FAIL: cmplx_allocate_1d wrong size"
        test_passed = .false.
    else
        print *, "PASS: cmplx_allocate_1d"
        
        ! Fill array
        do i = 1, n
            z_array(i) = cmplx(real(i, real64), real(i, real64), real64)
        end do
        
        ! Test array sum
        z1 = cmplx_array_sum(z_array)
        result = cmplx(15.0_real64, 15.0_real64, real64)
        if (abs(z1 - result) > tol) then
            print *, "FAIL: cmplx_array_sum"
            test_passed = .false.
        else
            print *, "PASS: cmplx_array_sum = ", z1
        end if
    end if
    
    call cmplx_deallocate_1d(z_array)
    
    print *, ""
    print *, "Testing complex matrix operations..."
    
    ! Test matrix allocation
    call cmplx_allocate_2d(z_matrix, 3, 3)
    
    if (.not. allocated(z_matrix)) then
        print *, "FAIL: cmplx_allocate_2d failed"
        test_passed = .false.
    else
        print *, "PASS: cmplx_allocate_2d"
        
        ! Fill matrix with identity
        z_matrix = cmplx_O
        do i = 1, 3
            z_matrix(i,i) = cmplx_E
        end do
        
        ! Test trace
        z1 = cmplx_matrix_trace(z_matrix)
        if (abs(z1 - cmplx(3.0_real64, 0.0_real64, real64)) > tol) then
            print *, "FAIL: cmplx_matrix_trace of identity"
            test_passed = .false.
        else
            print *, "PASS: cmplx_matrix_trace = ", z1
        end if
    end if
    
    call cmplx_deallocate_2d(z_matrix)
    
    print *, ""
    if (test_passed) then
        print *, "All tests PASSED!"
    else
        print *, "Some tests FAILED!"
        stop 1
    end if
    
end program test_kilca_complex