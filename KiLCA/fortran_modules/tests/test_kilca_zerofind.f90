program test_kilca_zerofind
    use iso_fortran_env
    use kilca_types_m
    implicit none
    
    integer :: total_tests, passed_tests, failed_tests
    logical :: test_passed
    character(len=100) :: test_name
    
    ! Initialize test counters
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    
    print *, ""
    print *, "========================================================"
    print *, "KiLCA Zero-Finding Library - Unit Tests"
    print *, "========================================================"
    print *, ""
    
    ! Run all test suites
    call test_basic_root_finding()
    call test_complex_root_finding()
    call test_multiple_roots()
    call test_winding_number()
    call test_rectangular_region()
    call test_newton_iterations()
    call test_argument_principle()
    call test_recursive_subdivision()
    call test_edge_cases()
    call test_convergence_criteria()
    call test_performance_settings()
    call test_determinant_zeros()
    
    ! Print summary
    print *, ""
    print *, "========================================================"
    print *, "TEST SUMMARY"
    print *, "========================================================"
    print '(A,I5)', " Total tests run:     ", total_tests
    print '(A,I5)', " Tests passed:        ", passed_tests  
    print '(A,I5)', " Tests failed:        ", failed_tests
    print '(A,F8.2,A)', " Success rate:        ", &
            100.0_dp * real(passed_tests, dp) / real(total_tests, dp), " %"
    
    if (failed_tests > 0) then
        print *, ""
        print *, " *** SOME TESTS FAILED! ***"
        stop 1
    else
        print *, ""
        print *, " *** ALL TESTS PASSED! ***"
    end if
    
contains

    !---------------------------------------------------------------------------
    ! Helper subroutines
    !---------------------------------------------------------------------------
    subroutine start_test(name)
        character(len=*), intent(in) :: name
        test_name = name
        write(*, '(A,A,A)', advance='no') "Testing ", trim(name), " ... "
    end subroutine start_test
    
    subroutine end_test(passed)
        logical, intent(in) :: passed
        total_tests = total_tests + 1
        if (passed) then
            passed_tests = passed_tests + 1
            print *, " PASSED"
        else
            failed_tests = failed_tests + 1
            print *, " FAILED"
        end if
    end subroutine end_test
    
    !---------------------------------------------------------------------------
    ! Test 1: Basic root finding for real polynomial
    !---------------------------------------------------------------------------
    subroutine test_basic_root_finding()
        use kilca_zerofind_m
        complex(dp) :: roots(10)
        integer :: n_roots, ierr
        real(dp) :: xmin, xmax, ymin, ymax
        type(zerofind_settings_t) :: settings
        
        call start_test("Basic root finding")
        test_passed = .true.
        
        ! Test finding roots of (z-1)(z-2)(z-3) = z^3 - 6z^2 + 11z - 6
        ! in region [-1, 4] x [-1, 1]
        xmin = -1.0_dp
        xmax = 4.0_dp
        ymin = -1.0_dp
        ymax = 1.0_dp
        
        ! Create settings
        call zerofind_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Set basic parameters
        settings%use_winding = .false.
        settings%max_partition_level = 64
        settings%eps_abs = 1.0e-12_dp
        settings%eps_rel = 1.0e-12_dp
        
        ! Find zeros
        n_roots = 0
        call zerofind_find_roots(polynomial_func, xmin, xmax, ymin, ymax, &
                                settings, roots, n_roots, ierr)
        
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (n_roots == 3)
        
        ! Check if roots are approximately 1, 2, 3
        if (n_roots == 3) then
            test_passed = test_passed .and. &
                         (any(abs(roots(1:3) - cmplx(1.0_dp, 0.0_dp, dp)) < 1.0e-10_dp))
            test_passed = test_passed .and. &
                         (any(abs(roots(1:3) - cmplx(2.0_dp, 0.0_dp, dp)) < 1.0e-10_dp))
            test_passed = test_passed .and. &
                         (any(abs(roots(1:3) - cmplx(3.0_dp, 0.0_dp, dp)) < 1.0e-10_dp))
        end if
        
        ! Clean up
        call zerofind_settings_destroy(settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_basic_root_finding
    
    ! Test function: polynomial (z-1)(z-2)(z-3)
    function polynomial_func(z) result(f)
        complex(dp), intent(in) :: z
        complex(dp) :: f
        f = (z - cmplx(1.0_dp, 0.0_dp, dp)) * &
            (z - cmplx(2.0_dp, 0.0_dp, dp)) * &
            (z - cmplx(3.0_dp, 0.0_dp, dp))
    end function polynomial_func
    
    !---------------------------------------------------------------------------
    ! Test 2: Complex root finding
    !---------------------------------------------------------------------------
    subroutine test_complex_root_finding()
        use kilca_zerofind_m
        complex(dp) :: roots(10)
        integer :: n_roots, ierr
        real(dp) :: xmin, xmax, ymin, ymax
        type(zerofind_settings_t) :: settings
        
        call start_test("Complex root finding")
        test_passed = .true.
        
        ! Test finding roots of (z - (1+i))(z - (1-i)) = z^2 - 2z + 2
        xmin = 0.0_dp
        xmax = 2.0_dp
        ymin = -2.0_dp
        ymax = 2.0_dp
        
        call zerofind_settings_create(settings, ierr)
        settings%use_winding = .true.
        settings%max_partition_level = 64
        
        n_roots = 0
        call zerofind_find_roots(complex_poly_func, xmin, xmax, ymin, ymax, &
                                settings, roots, n_roots, ierr)
        
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (n_roots == 2)
        
        ! Check roots are 1+i and 1-i
        if (n_roots == 2) then
            test_passed = test_passed .and. &
                         (any(abs(roots(1:2) - cmplx(1.0_dp, 1.0_dp, dp)) < 1.0e-10_dp))
            test_passed = test_passed .and. &
                         (any(abs(roots(1:2) - cmplx(1.0_dp, -1.0_dp, dp)) < 1.0e-10_dp))
        end if
        
        call zerofind_settings_destroy(settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_complex_root_finding
    
    ! Test function: (z - (1+i))(z - (1-i))
    function complex_poly_func(z) result(f)
        complex(dp), intent(in) :: z
        complex(dp) :: f
        f = (z - cmplx(1.0_dp, 1.0_dp, dp)) * &
            (z - cmplx(1.0_dp, -1.0_dp, dp))
    end function complex_poly_func
    
    !---------------------------------------------------------------------------
    ! Test 3: Multiple roots handling
    !---------------------------------------------------------------------------
    subroutine test_multiple_roots()
        use kilca_zerofind_m
        complex(dp) :: roots(10)
        integer :: n_roots, ierr
        type(zerofind_settings_t) :: settings
        
        call start_test("Multiple roots handling")
        test_passed = .true.
        
        ! Test finding roots of (z-1)^3 = z^3 - 3z^2 + 3z - 1
        call zerofind_settings_create(settings, ierr)
        settings%handle_multiplicities = .true.
        
        n_roots = 0
        call zerofind_find_roots(multiple_root_func, 0.0_dp, 2.0_dp, &
                                -1.0_dp, 1.0_dp, settings, roots, n_roots, ierr)
        
        test_passed = test_passed .and. (ierr == 0)
        ! Should find one root with multiplicity 3
        test_passed = test_passed .and. (n_roots >= 1)
        
        if (n_roots >= 1) then
            test_passed = test_passed .and. &
                         (abs(roots(1) - cmplx(1.0_dp, 0.0_dp, dp)) < 1.0e-10_dp)
        end if
        
        call zerofind_settings_destroy(settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_multiple_roots
    
    ! Test function: (z-1)^3
    function multiple_root_func(z) result(f)
        complex(dp), intent(in) :: z
        complex(dp) :: f
        complex(dp) :: z1
        z1 = z - cmplx(1.0_dp, 0.0_dp, dp)
        f = z1 * z1 * z1
    end function multiple_root_func
    
    !---------------------------------------------------------------------------
    ! Test 4: Winding number calculation
    !---------------------------------------------------------------------------
    subroutine test_winding_number()
        use kilca_zerofind_m
        integer :: winding_num, ierr
        real(dp) :: xmin, xmax, ymin, ymax
        type(zerofind_settings_t) :: settings
        
        call start_test("Winding number calculation")
        test_passed = .true.
        
        ! Test winding number around a simple pole
        xmin = -2.0_dp
        xmax = 2.0_dp
        ymin = -2.0_dp
        ymax = 2.0_dp
        
        call zerofind_settings_create(settings, ierr)
        
        ! Calculate winding number for function with zero at origin
        call zerofind_calc_winding_number(simple_zero_func, xmin, xmax, ymin, ymax, &
                                         settings, winding_num, ierr)
        
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (winding_num == 1)
        
        call zerofind_settings_destroy(settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_winding_number
    
    ! Test function: z (simple zero at origin)
    function simple_zero_func(z) result(f)
        complex(dp), intent(in) :: z
        complex(dp) :: f
        f = z
    end function simple_zero_func
    
    !---------------------------------------------------------------------------
    ! Test 5: Rectangular region handling
    !---------------------------------------------------------------------------
    subroutine test_rectangular_region()
        use kilca_zerofind_m
        type(rectangle_t) :: rect
        integer :: ierr
        logical :: contains_point
        complex(dp) :: test_point
        
        call start_test("Rectangular region handling")
        test_passed = .true.
        
        ! Create rectangle
        call rectangle_create(rect, -1.0_dp, 1.0_dp, -2.0_dp, 2.0_dp, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test containment
        test_point = cmplx(0.0_dp, 0.0_dp, dp)
        call rectangle_contains(rect, test_point, contains_point)
        test_passed = test_passed .and. contains_point
        
        test_point = cmplx(2.0_dp, 0.0_dp, dp)
        call rectangle_contains(rect, test_point, contains_point)
        test_passed = test_passed .and. (.not. contains_point)
        
        ! Test subdivision
        call rectangle_subdivide(rect, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call rectangle_destroy(rect, ierr)
        
        call end_test(test_passed)
    end subroutine test_rectangular_region
    
    !---------------------------------------------------------------------------
    ! Test 6: Newton iterations
    !---------------------------------------------------------------------------
    subroutine test_newton_iterations()
        use kilca_zerofind_m
        complex(dp) :: z0, z_root
        integer :: n_iter, ierr
        type(newton_settings_t) :: newton_settings
        
        call start_test("Newton iterations")
        test_passed = .true.
        
        ! Create Newton settings
        call newton_settings_create(newton_settings, ierr)
        newton_settings%max_iter = 50
        newton_settings%eps_abs = 1.0e-12_dp
        newton_settings%eps_rel = 1.0e-12_dp
        
        ! Start near a root of z^2 - 1
        z0 = cmplx(0.8_dp, 0.0_dp, dp)
        
        call newton_solve(quadratic_func, quadratic_deriv, z0, newton_settings, &
                         z_root, n_iter, ierr)
        
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (n_iter < 20)
        test_passed = test_passed .and. &
                     (abs(z_root - cmplx(1.0_dp, 0.0_dp, dp)) < 1.0e-10_dp)
        
        call newton_settings_destroy(newton_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_newton_iterations
    
    ! Test function: z^2 - 1
    function quadratic_func(z) result(f)
        complex(dp), intent(in) :: z
        complex(dp) :: f
        f = z*z - cmplx(1.0_dp, 0.0_dp, dp)
    end function quadratic_func
    
    ! Derivative: 2z
    function quadratic_deriv(z) result(df)
        complex(dp), intent(in) :: z
        complex(dp) :: df
        df = 2.0_dp * z
    end function quadratic_deriv
    
    !---------------------------------------------------------------------------
    ! Test 7: Argument principle
    !---------------------------------------------------------------------------
    subroutine test_argument_principle()
        use kilca_zerofind_m
        real(dp) :: arg_change
        integer :: n_zeros, ierr
        type(zerofind_settings_t) :: settings
        
        call start_test("Argument principle")
        test_passed = .true.
        
        call zerofind_settings_create(settings, ierr)
        
        ! Calculate argument change for rational function
        call zerofind_calc_arg_change(rational_func, -2.0_dp, 2.0_dp, &
                                     -2.0_dp, 2.0_dp, settings, &
                                     arg_change, ierr)
        
        test_passed = test_passed .and. (ierr == 0)
        
        ! Number of zeros from argument principle
        n_zeros = nint(arg_change / (2.0_dp * PI))
        test_passed = test_passed .and. (n_zeros == 2)
        
        call zerofind_settings_destroy(settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_argument_principle
    
    ! Test function: (z-0.5)(z+0.5)/(z-3)
    function rational_func(z) result(f)
        complex(dp), intent(in) :: z
        complex(dp) :: f
        complex(dp) :: num, den
        num = (z - cmplx(0.5_dp, 0.0_dp, dp)) * &
              (z + cmplx(0.5_dp, 0.0_dp, dp))
        den = z - cmplx(3.0_dp, 0.0_dp, dp)
        if (abs(den) > 1.0e-10_dp) then
            f = num / den
        else
            f = cmplx(1.0e10_dp, 0.0_dp, dp)  ! Large value near pole
        end if
    end function rational_func
    
    !---------------------------------------------------------------------------
    ! Test 8: Recursive subdivision
    !---------------------------------------------------------------------------
    subroutine test_recursive_subdivision()
        use kilca_zerofind_m
        type(rectangle_tree_t) :: tree
        integer :: ierr, n_leaves
        
        call start_test("Recursive subdivision")
        test_passed = .true.
        
        ! Create subdivision tree
        call rect_tree_create(tree, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Subdivide based on function properties
        call rect_tree_adaptive_subdivide(tree, polynomial_func, 4, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Count leaves
        call rect_tree_count_leaves(tree, n_leaves)
        test_passed = test_passed .and. (n_leaves > 1)
        
        call rect_tree_destroy(tree, ierr)
        
        call end_test(test_passed)
    end subroutine test_recursive_subdivision
    
    !---------------------------------------------------------------------------
    ! Test 9: Edge cases
    !---------------------------------------------------------------------------
    subroutine test_edge_cases()
        use kilca_zerofind_m
        complex(dp) :: roots(10)
        integer :: n_roots, ierr
        type(zerofind_settings_t) :: settings
        
        call start_test("Edge cases handling")
        test_passed = .true.
        
        call zerofind_settings_create(settings, ierr)
        
        ! Test with no zeros in region
        n_roots = 0
        call zerofind_find_roots(no_zeros_func, -1.0_dp, 1.0_dp, &
                                -1.0_dp, 1.0_dp, settings, roots, n_roots, ierr)
        
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (n_roots == 0)
        
        ! Test with zero on boundary
        n_roots = 0
        call zerofind_find_roots(boundary_zero_func, -1.0_dp, 1.0_dp, &
                                -1.0_dp, 1.0_dp, settings, roots, n_roots, ierr)
        
        test_passed = test_passed .and. (ierr == 0)
        
        call zerofind_settings_destroy(settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_edge_cases
    
    ! Test function: z^2 + 1 (no real zeros)
    function no_zeros_func(z) result(f)
        complex(dp), intent(in) :: z
        complex(dp) :: f
        f = z*z + cmplx(1.0_dp, 0.0_dp, dp)
    end function no_zeros_func
    
    ! Test function: z - 1 (zero on boundary at z=1)
    function boundary_zero_func(z) result(f)
        complex(dp), intent(in) :: z
        complex(dp) :: f
        f = z - cmplx(1.0_dp, 0.0_dp, dp)
    end function boundary_zero_func
    
    !---------------------------------------------------------------------------
    ! Test 10: Convergence criteria
    !---------------------------------------------------------------------------
    subroutine test_convergence_criteria()
        use kilca_zerofind_m
        type(zerofind_settings_t) :: settings
        logical :: converged
        complex(dp) :: z1, z2
        integer :: ierr
        
        call start_test("Convergence criteria")
        test_passed = .true.
        
        call zerofind_settings_create(settings, ierr)
        settings%eps_abs = 1.0e-10_dp
        settings%eps_rel = 1.0e-10_dp
        
        ! Test absolute convergence
        z1 = cmplx(1.0_dp, 0.0_dp, dp)
        z2 = cmplx(1.0_dp + 1.0e-11_dp, 0.0_dp, dp)
        call zerofind_check_convergence(z1, z2, settings, converged)
        test_passed = test_passed .and. converged
        
        ! Test relative convergence
        z1 = cmplx(1000.0_dp, 0.0_dp, dp)
        z2 = cmplx(1000.0_dp * (1.0_dp + 1.0e-11_dp), 0.0_dp, dp)
        call zerofind_check_convergence(z1, z2, settings, converged)
        test_passed = test_passed .and. converged
        
        call zerofind_settings_destroy(settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_convergence_criteria
    
    !---------------------------------------------------------------------------
    ! Test 11: Performance settings
    !---------------------------------------------------------------------------
    subroutine test_performance_settings()
        use kilca_zerofind_m
        type(zerofind_settings_t) :: settings
        integer :: ierr
        
        call start_test("Performance settings")
        test_passed = .true.
        
        call zerofind_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test setting various performance parameters
        settings%min_recursion_level = 4
        settings%max_partition_level = 128
        settings%interpolation_error = 1.0e-3_dp
        settings%jump_error = 1.0e-3_dp
        settings%n_split_x = 4
        settings%n_split_y = 4
        
        ! Validate settings
        call zerofind_settings_validate(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call zerofind_settings_destroy(settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_performance_settings
    
    !---------------------------------------------------------------------------
    ! Test 12: Determinant zeros (plasma physics application)
    !---------------------------------------------------------------------------
    subroutine test_determinant_zeros()
        use kilca_zerofind_m
        complex(dp) :: roots(20)
        integer :: n_roots, ierr
        type(zerofind_settings_t) :: settings
        
        call start_test("Determinant zeros finding")
        test_passed = .true.
        
        ! Simulate plasma dispersion determinant
        call zerofind_settings_create(settings, ierr)
        settings%use_winding = .true.
        settings%target_n_zeros = 5
        settings%eps_residual = 1.0e-10_dp
        
        ! Search in typical frequency range
        n_roots = 0
        call zerofind_find_roots(plasma_det_func, 0.1_dp, 10.0_dp, &
                                -1.0_dp, 1.0_dp, settings, roots, n_roots, ierr)
        
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (n_roots > 0)
        
        ! Verify residuals
        if (n_roots > 0) then
            test_passed = test_passed .and. &
                         (abs(plasma_det_func(roots(1))) < settings%eps_residual)
        end if
        
        call zerofind_settings_destroy(settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_determinant_zeros
    
    ! Test function: Simplified plasma dispersion determinant
    function plasma_det_func(z) result(f)
        complex(dp), intent(in) :: z
        complex(dp) :: f
        ! Simplified model: zeros at z = 1, 2, 3, 4, 5
        f = (z - cmplx(1.0_dp, 0.0_dp, dp)) * &
            (z - cmplx(2.0_dp, 0.0_dp, dp)) * &
            (z - cmplx(3.0_dp, 0.0_dp, dp)) * &
            (z - cmplx(4.0_dp, 0.0_dp, dp)) * &
            (z - cmplx(5.0_dp, 0.0_dp, dp))
    end function plasma_det_func
    
end program test_kilca_zerofind