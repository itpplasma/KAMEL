program test_kilca_shared
    use iso_fortran_env, only: real64, int32
    use kilca_shared_m
    implicit none
    
    logical :: test_passed
    integer :: i, j, n
    real(real64) :: x, result
    real(real64), allocatable :: BC(:,:), BC_ref(:,:)
    real(real64) :: W(5), x0, L
    integer :: compare_result
    real(real64) :: a(5), b(5)
    
    test_passed = .true.
    
    print *, "Testing kilca_shared_m module..."
    print *, ""
    
    ! Test signum function
    print *, "Testing signum function..."
    
    if (signum(-5.0_real64) /= -1) then
        print *, "FAIL: signum(-5.0) should be -1"
        test_passed = .false.
    else
        print *, "PASS: signum(-5.0) = ", signum(-5.0_real64)
    end if
    
    if (signum(0.0_real64) /= 0) then
        print *, "FAIL: signum(0.0) should be 0"
        test_passed = .false.
    else
        print *, "PASS: signum(0.0) = ", signum(0.0_real64)
    end if
    
    if (signum(3.14_real64) /= 1) then
        print *, "FAIL: signum(3.14) should be 1"
        test_passed = .false.
    else
        print *, "PASS: signum(3.14) = ", signum(3.14_real64)
    end if
    
    print *, ""
    print *, "Testing compare_doubles function..."
    
    ! Test compare_doubles
    a(1) = 3.14_real64; b(1) = 2.71_real64
    compare_result = compare_doubles(a(1), b(1))
    if (compare_result /= 1) then
        print *, "FAIL: compare_doubles(3.14, 2.71) should return 1"
        test_passed = .false.
    else
        print *, "PASS: compare_doubles(3.14, 2.71) = ", compare_result
    end if
    
    a(1) = 1.0_real64; b(1) = 1.0_real64
    compare_result = compare_doubles(a(1), b(1))
    if (compare_result /= 0) then
        print *, "FAIL: compare_doubles(1.0, 1.0) should return 0"
        test_passed = .false.
    else
        print *, "PASS: compare_doubles(1.0, 1.0) = ", compare_result
    end if
    
    a(1) = -1.0_real64; b(1) = 0.0_real64
    compare_result = compare_doubles(a(1), b(1))
    if (compare_result /= -1) then
        print *, "FAIL: compare_doubles(-1.0, 0.0) should return -1"
        test_passed = .false.
    else
        print *, "PASS: compare_doubles(-1.0, 0.0) = ", compare_result
    end if
    
    print *, ""
    print *, "Testing binomial coefficients..."
    
    ! Test binomial coefficients for n=4
    n = 4
    allocate(BC(0:n, 0:n))
    allocate(BC_ref(0:n, 0:n))
    
    ! Reference values: Pascal's triangle
    BC_ref = 0.0_real64
    BC_ref(0,0) = 1.0_real64
    BC_ref(1,0) = 1.0_real64; BC_ref(1,1) = 1.0_real64
    BC_ref(2,0) = 1.0_real64; BC_ref(2,1) = 2.0_real64; BC_ref(2,2) = 1.0_real64
    BC_ref(3,0) = 1.0_real64; BC_ref(3,1) = 3.0_real64; BC_ref(3,2) = 3.0_real64; BC_ref(3,3) = 1.0_real64
    BC_ref(4,0) = 1.0_real64; BC_ref(4,1) = 4.0_real64; BC_ref(4,2) = 6.0_real64; BC_ref(4,3) = 4.0_real64; BC_ref(4,4) = 1.0_real64
    
    call binomial_coefficients(n, BC)
    
    ! Check some values
    if (abs(BC(2,0) - 1.0_real64) > epsilon(1.0_real64)) then
        print *, "FAIL: C(2,0) should be 1"
        test_passed = .false.
    else
        print *, "PASS: C(2,0) = ", BC(2,0)
    end if
    
    if (abs(BC(4,2) - 6.0_real64) > epsilon(1.0_real64)) then
        print *, "FAIL: C(4,2) should be 6"
        test_passed = .false.
    else
        print *, "PASS: C(4,2) = ", BC(4,2)
    end if
    
    if (abs(BC(3,2) - 3.0_real64) > epsilon(1.0_real64)) then
        print *, "FAIL: C(3,2) should be 3"
        test_passed = .false.
    else
        print *, "PASS: C(3,2) = ", BC(3,2)
    end if
    
    deallocate(BC, BC_ref)
    
    print *, ""
    print *, "Testing localizator function..."
    
    ! Test localizator at different points
    x0 = 0.0_real64
    L = 1.0_real64
    
    ! Test at x = x0 (should be close to 1)
    x = x0
    call localizator(x, x0, L, W)
    if (abs(W(1) - 1.0_real64) > 0.1_real64) then
        print *, "FAIL: localizator at x=x0 should be close to 1"
        test_passed = .false.
    else
        print *, "PASS: localizator(x0, x0, L) W(0) = ", W(1)
    end if
    
    ! Test far from x0 (should be close to 0)
    x = x0 + 2.0_real64 * L
    call localizator(x, x0, L, W)
    if (abs(W(1)) > 0.1_real64) then
        print *, "FAIL: localizator far from x0 should be close to 0"
        test_passed = .false.
    else
        print *, "PASS: localizator(x0+2L, x0, L) W(0) = ", W(1)
    end if
    
    print *, ""
    print *, "Testing localizator_4_derivs function..."
    
    ! Test localizator_4_derivs at x = x0
    x = x0
    call localizator_4_derivs(x, x0, L, W)
    if (abs(W(1) - 1.0_real64) > 0.1_real64) then
        print *, "FAIL: localizator_4_derivs at x=x0 should be close to 1"
        test_passed = .false.
    else
        print *, "PASS: localizator_4_derivs(x0, x0, L) W(0) = ", W(1)
        print *, "      First derivative W(1) = ", W(2)
    end if
    
    print *, ""
    print *, "Testing memory management utilities..."
    
    ! Test safe allocation/deallocation
    block
        real(real64), allocatable :: test_array_1d(:)
        real(real64), allocatable :: test_array_2d(:,:)
        integer(int32), allocatable :: test_int_array(:)
        complex(real64), allocatable :: test_cmplx_array(:)
        
        ! Test 1D real allocation
        call safe_allocate_real64_1d(test_array_1d, 100, "test_array_1d")
        if (.not. allocated(test_array_1d)) then
            print *, "FAIL: safe_allocate_real64_1d failed"
            test_passed = .false.
        else if (size(test_array_1d) /= 100) then
            print *, "FAIL: safe_allocate_real64_1d wrong size"
            test_passed = .false.
        else
            print *, "PASS: safe_allocate_real64_1d"
            test_array_1d = 3.14_real64
        end if
        
        ! Test 2D real allocation
        call safe_allocate_real64_2d(test_array_2d, 10, 20, "test_array_2d")
        if (.not. allocated(test_array_2d)) then
            print *, "FAIL: safe_allocate_real64_2d failed"
            test_passed = .false.
        else if (size(test_array_2d, 1) /= 10 .or. size(test_array_2d, 2) /= 20) then
            print *, "FAIL: safe_allocate_real64_2d wrong size"
            test_passed = .false.
        else
            print *, "PASS: safe_allocate_real64_2d"
        end if
        
        ! Test integer allocation
        call safe_allocate_int32_1d(test_int_array, 50, "test_int_array")
        if (.not. allocated(test_int_array)) then
            print *, "FAIL: safe_allocate_int32_1d failed"
            test_passed = .false.
        else
            print *, "PASS: safe_allocate_int32_1d"
        end if
        
        ! Test complex allocation
        call safe_allocate_cmplx_1d(test_cmplx_array, 25, "test_cmplx_array")
        if (.not. allocated(test_cmplx_array)) then
            print *, "FAIL: safe_allocate_cmplx_1d failed"
            test_passed = .false.
        else
            print *, "PASS: safe_allocate_cmplx_1d"
        end if
        
        ! Test deallocations
        call safe_deallocate_real64_1d(test_array_1d, "test_array_1d")
        if (allocated(test_array_1d)) then
            print *, "FAIL: safe_deallocate_real64_1d failed"
            test_passed = .false.
        else
            print *, "PASS: safe_deallocate_real64_1d"
        end if
        
        call safe_deallocate_real64_2d(test_array_2d, "test_array_2d")
        call safe_deallocate_int32_1d(test_int_array, "test_int_array")
        call safe_deallocate_cmplx_1d(test_cmplx_array, "test_cmplx_array")
    end block
    
    print *, ""
    if (test_passed) then
        print *, "All tests PASSED!"
    else
        print *, "Some tests FAILED!"
        stop 1
    end if
    
end program test_kilca_shared