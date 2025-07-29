program test_kilca_solver
    use iso_fortran_env, only: real64, int32, error_unit
    use kilca_solver_m
    use kilca_types_m
    implicit none
    
    ! Test counters
    integer :: total_tests, passed_tests, failed_tests
    logical :: test_passed
    
    ! Initialize test counters
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    
    print *, "========================================================"
    print *, "KiLCA Solver Framework - Unit Tests"
    print *, "========================================================"
    print *, ""
    
    ! Test 1: Solver settings creation and initialization
    call test_solver_settings()
    
    ! Test 2: Simple ODE integration (test problem)
    call test_simple_ode_integration()
    
    ! Test 3: Basis vector orthonormalization
    call test_basis_orthonormalization()
    
    ! Test 4: Basis vector superposition
    call test_basis_superposition()
    
    ! Test 5: Multi-vector ODE integration
    call test_multi_vector_integration()
    
    ! Test 6: Error handling and validation
    call test_solver_error_handling()
    
    ! Final summary
    print *, ""
    print *, "========================================================"
    print *, "TEST SUMMARY"
    print *, "========================================================"
    print *, "Total tests run:    ", total_tests
    print *, "Tests passed:       ", passed_tests
    print *, "Tests failed:       ", failed_tests
    print *, "Success rate:       ", real(passed_tests)/real(total_tests)*100.0, "%"
    print *, ""
    
    if (failed_tests == 0) then
        print *, "*** ALL TESTS PASSED! ***"
    else
        print *, "*** SOME TESTS FAILED! ***"
        stop 1
    end if
    
contains

    !---------------------------------------------------------------------------
    ! Test 1: Solver settings creation and initialization
    !---------------------------------------------------------------------------
    subroutine test_solver_settings()
        type(solver_settings_t) :: settings
        integer :: ierr
        
        call start_test("Solver settings creation")
        
        ! Test default initialization
        call solver_settings_create(settings, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            test_passed = (settings%Nort > 0)
            test_passed = test_passed .and. (settings%eps_rel > 0.0_real64)
            test_passed = test_passed .and. (settings%eps_abs > 0.0_real64)
            test_passed = test_passed .and. (settings%norm_fac > 1.0_real64)
            test_passed = test_passed .and. (settings%debug >= 0)
        end if
        
        ! Test custom settings
        call solver_settings_set(settings, 10, 1.0e-8_real64, 1.0e-12_real64, 100.0_real64, 1, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (settings%Nort == 10)
        test_passed = test_passed .and. (abs(settings%eps_rel - 1.0e-8_real64) < 1.0e-15_real64)
        
        call end_test(test_passed)
    end subroutine test_solver_settings
    
    !---------------------------------------------------------------------------
    ! Test 2: Simple ODE integration (test problem: y' = y, y(0) = 1)
    !---------------------------------------------------------------------------
    subroutine test_simple_ode_integration()
        type(solver_settings_t) :: settings
        real(real64), parameter :: r0 = 0.0_real64, r1 = 1.0_real64
        real(real64), parameter :: expected_result = exp(1.0_real64)  ! e^1
        real(real64), allocatable :: rvec(:), y(:)
        real(real64) :: error
        integer :: dim, ierr, i
        
        call start_test("Simple ODE integration")
        
        ! Setup test problem
        call solver_settings_create(settings, ierr)
        call solver_settings_set(settings, 5, 1.0e-10_real64, 1.0e-14_real64, 100.0_real64, 0, ierr)
        
        dim = 11
        allocate(rvec(dim), y(dim))
        
        ! Create grid from 0 to 1
        do i = 1, dim
            rvec(i) = real(i-1, real64) / real(dim-1, real64)
        end do
        
        ! Initial condition: y(0) = 1
        y(1) = 1.0_real64
        
        ! Integrate simple exponential ODE
        call solver_integrate_ode(simple_exponential_rhs, 1, 1, dim, rvec, y, settings, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            error = abs(y(dim) - expected_result)
            test_passed = (error < 1.0e-6_real64)  ! Reasonable tolerance for simple ODE
        end if
        
        if (.not. test_passed) then
            print *, "  Simple ODE integration failed, ierr =", ierr, ", error =", error
            print *, "  Expected:", expected_result, ", Got:", y(dim)
        end if
        
        deallocate(rvec, y)
        call end_test(test_passed)
    end subroutine test_simple_ode_integration
    
    !---------------------------------------------------------------------------
    ! Test 3: Basis vector orthonormalization
    !---------------------------------------------------------------------------
    subroutine test_basis_orthonormalization()
        real(real64), allocatable :: vectors(:, :), result(:, :)
        real(real64) :: dot_product, error
        integer :: n, nvec, i, j, ierr
        
        call start_test("Basis vector orthonormalization")
        
        ! Create test vectors that are not orthogonal
        n = 4      ! Vector dimension
        nvec = 3   ! Number of vectors
        allocate(vectors(n, nvec), result(n, nvec))
        
        ! Set up linearly independent but non-orthogonal vectors
        vectors(:, 1) = [1.0_real64, 1.0_real64, 0.0_real64, 0.0_real64]
        vectors(:, 2) = [1.0_real64, 0.0_real64, 1.0_real64, 0.0_real64]
        vectors(:, 3) = [1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64]
        
        ! Orthonormalize
        call solver_orthonormalize_vectors(n, nvec, vectors, result, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            ! Check orthogonality: <vi, vj> = δij
            do i = 1, nvec
                do j = 1, nvec
                    dot_product = sum(result(:, i) * result(:, j))
                    if (i == j) then
                        error = abs(dot_product - 1.0_real64)  ! Should be 1 (normalized)
                    else
                        error = abs(dot_product)  ! Should be 0 (orthogonal)
                    end if
                    test_passed = test_passed .and. (error < 1.0e-12_real64)
                end do
            end do
        end if
        
        deallocate(vectors, result)
        call end_test(test_passed)
    end subroutine test_basis_orthonormalization
    
    !---------------------------------------------------------------------------
    ! Test 4: Basis vector superposition
    !---------------------------------------------------------------------------
    subroutine test_basis_superposition()
        real(real64), allocatable :: basis(:, :), coeffs(:), result(:), expected(:)
        real(real64) :: error
        integer :: n, nvec, ierr
        
        call start_test("Basis vector superposition")
        
        ! Create test basis vectors
        n = 3      ! Vector dimension
        nvec = 2   ! Number of basis vectors
        allocate(basis(n, nvec), coeffs(nvec), result(n), expected(n))
        
        ! Set up simple basis vectors
        basis(:, 1) = [1.0_real64, 0.0_real64, 0.0_real64]
        basis(:, 2) = [0.0_real64, 1.0_real64, 0.0_real64]
        
        ! Coefficients for linear combination
        coeffs = [2.0_real64, 3.0_real64]
        
        ! Expected result: 2*e1 + 3*e2
        expected = [2.0_real64, 3.0_real64, 0.0_real64]
        
        ! Compute superposition
        call solver_superpose_vectors(n, nvec, basis, coeffs, result, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            error = maxval(abs(result - expected))
            test_passed = (error < 1.0e-14_real64)
        end if
        
        deallocate(basis, coeffs, result, expected)
        call end_test(test_passed)
    end subroutine test_basis_superposition
    
    !---------------------------------------------------------------------------
    ! Test 5: Multi-vector ODE integration (coupled system)
    !---------------------------------------------------------------------------
    subroutine test_multi_vector_integration()
        type(solver_settings_t) :: settings
        real(real64), allocatable :: rvec(:), y(:)
        real(real64) :: error
        integer :: dim, Nfs, Nw, ierr, i
        
        call start_test("Multi-vector ODE integration")
        
        ! Setup coupled system: y1' = y2, y2' = -y1 (harmonic oscillator)
        call solver_settings_create(settings, ierr)
        
        Nfs = 2  ! Two fundamental solutions
        Nw = 1   ! One wave
        dim = 21
        
        allocate(rvec(dim), y(2*Nfs*Nw*dim))  ! 2*Nfs*Nw = 4 equations per grid point
        
        ! Create grid from 0 to π
        do i = 1, dim
            rvec(i) = real(i-1, real64) * 3.14159265359_real64 / real(dim-1, real64)
        end do
        
        ! Initial conditions: y1(0)=1, y2(0)=0 for first solution
        !                     y1(0)=0, y2(0)=1 for second solution
        y = 0.0_real64
        y(1) = 1.0_real64  ! First solution: [1, 0]
        y(3) = 1.0_real64  ! Second solution: [0, 1]
        
        ! Integrate coupled system
        call solver_integrate_basis_vectors(harmonic_oscillator_rhs, Nfs, Nw, dim, rvec, y, settings, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            ! At r = π, harmonic oscillator should give: cos(π) = -1, sin(π) = 0
            error = abs(y(4*(dim-1)+1) + 1.0_real64)  ! y1(π) should be -1
            error = max(error, abs(y(4*(dim-1)+2)))    ! y2(π) should be 0
            test_passed = (error < 1.0e-3_real64)  ! Reasonable tolerance for numerical integration
        end if
        
        deallocate(rvec, y)
        call end_test(test_passed)
    end subroutine test_multi_vector_integration
    
    !---------------------------------------------------------------------------
    ! Test 6: Error handling and validation
    !---------------------------------------------------------------------------
    subroutine test_solver_error_handling()
        type(solver_settings_t) :: settings
        real(real64), allocatable :: rvec(:), y(:)
        integer :: ierr
        
        call start_test("Solver error handling")
        
        call solver_settings_create(settings, ierr)
        
        ! Test with invalid dimensions
        allocate(rvec(0), y(0))
        call solver_integrate_ode(simple_exponential_rhs, 1, 1, 0, rvec, y, settings, ierr)
        test_passed = (ierr /= 0)  ! Should fail
        deallocate(rvec, y)
        
        ! Test with invalid grid (non-monotonic)
        allocate(rvec(3), y(3))
        rvec = [0.0_real64, 2.0_real64, 1.0_real64]  ! Non-monotonic
        call solver_integrate_ode(simple_exponential_rhs, 1, 1, 3, rvec, y, settings, ierr)
        test_passed = test_passed .and. (ierr /= 0)  ! Should fail
        deallocate(rvec, y)
        
        call end_test(test_passed)
    end subroutine test_solver_error_handling
    
    !---------------------------------------------------------------------------
    ! Test RHS functions
    !---------------------------------------------------------------------------
    subroutine simple_exponential_rhs(r, y, ydot)
        real(real64), intent(in) :: r
        real(real64), intent(in) :: y(*)
        real(real64), intent(out) :: ydot(*)
        
        ! y' = y (exponential growth)
        ydot(1) = y(1)
    end subroutine simple_exponential_rhs
    
    subroutine harmonic_oscillator_rhs(r, y, ydot)
        real(real64), intent(in) :: r
        real(real64), intent(in) :: y(*)
        real(real64), intent(out) :: ydot(*)
        
        ! Harmonic oscillator: y1' = y2, y2' = -y1
        ! For multiple fundamental solutions, equations are interleaved
        integer :: i, base
        
        ! For 2 fundamental solutions, we have 4 equations total
        do i = 1, 2  ! Hard-coded for simplicity since we can't use size(y) with assumed-size
            base = 2*(i-1)
            ydot(base+1) = y(base+2)   ! y1' = y2
            ydot(base+2) = -y(base+1)  ! y2' = -y1
        end do
    end subroutine harmonic_oscillator_rhs
    
    !---------------------------------------------------------------------------
    ! Test utilities  
    !---------------------------------------------------------------------------
    subroutine start_test(name)
        character(len=*), intent(in) :: name
        total_tests = total_tests + 1
        write(*,'(A,A)', advance='no') "Testing ", name
        write(*,'(A)', advance='no') " ... "
    end subroutine start_test
    
    subroutine end_test(passed)
        logical, intent(in) :: passed
        if (passed) then
            print *, "PASSED"
            passed_tests = passed_tests + 1
        else
            print *, "FAILED"
            failed_tests = failed_tests + 1
        end if
    end subroutine end_test

end program test_kilca_solver