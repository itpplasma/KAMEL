program test_kilca_eigensolver
    use iso_fortran_env, only: real64, int32
    use kilca_eigensolver_m
    implicit none
    
    ! Test counters
    integer :: total_tests, passed_tests, failed_tests
    logical :: test_passed
    
    ! Initialize test counters
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    
    print *, "===================================================="
    print *, "KiLCA Eigensolver Module - Unit Tests"
    print *, "===================================================="
    print *, ""
    
    ! Test 1: eigsys_zgeev replacement
    call test_eigsys_zgeev_replacement()
    
    ! Test 2: Hermitian eigenvalue problem
    call test_hermitian_eigenproblem()
    
    ! Test 3: General real eigenvalue problem
    call test_general_real_eigenproblem()
    
    ! Test 4: Symmetric real eigenvalue problem
    call test_symmetric_real_eigenproblem()
    
    ! Test 5: Mode eigenvalue problem with sorting
    call test_mode_eigenvalue_problem()
    
    ! Test 6: Growth rate sorting
    call test_growth_rate_sorting()
    
    ! Final summary
    print *, ""
    print *, "===================================================="
    print *, "TEST SUMMARY"
    print *, "===================================================="
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
    ! Test 1: eigsys_zgeev as replacement for original eigsys
    !---------------------------------------------------------------------------
    subroutine test_eigsys_zgeev_replacement()
        complex(real64) :: matrix(3,3), eigenvals(3), eigenvecs(3,3)
        integer :: ierr
        real(real64) :: tol
        
        call start_test("eigsys_zgeev replacement functionality")
        
        tol = 1.0e-10_real64
        
        ! Create test matrix with known eigenvalues
        matrix = cmplx(0.0_real64, 0.0_real64, real64)
        matrix(1,1) = cmplx(2.0_real64, 0.0_real64, real64)
        matrix(2,2) = cmplx(3.0_real64, 0.0_real64, real64)
        matrix(3,3) = cmplx(1.0_real64, 0.0_real64, real64)
        ! Add small perturbations
        matrix(1,2) = cmplx(0.1_real64, 0.05_real64, real64)
        matrix(2,3) = cmplx(0.05_real64, -0.1_real64, real64)
        
        call eigsys_zgeev(3, matrix, eigenvals, eigenvecs, ierr)
        
        test_passed = (ierr == 0)
        
        if (.not. test_passed) then
            print *, "  eigsys_zgeev failed with ierr =", ierr
        else
            print *, "  Eigenvalues:", eigenvals
        end if
        
        call end_test(test_passed)
    end subroutine test_eigsys_zgeev_replacement
    
    !---------------------------------------------------------------------------
    ! Test 2: Hermitian eigenvalue problem
    !---------------------------------------------------------------------------
    subroutine test_hermitian_eigenproblem()
        complex(real64) :: H_matrix(3,3), eigenvecs(3,3)
        real(real64) :: eigenvals(3)
        integer :: ierr, i
        
        call start_test("Hermitian eigenvalue problem")
        
        ! Create Hermitian matrix
        H_matrix(1,1) = cmplx(1.0_real64, 0.0_real64, real64)
        H_matrix(2,2) = cmplx(2.0_real64, 0.0_real64, real64)
        H_matrix(3,3) = cmplx(3.0_real64, 0.0_real64, real64)
        H_matrix(1,2) = cmplx(0.5_real64, 0.1_real64, real64)
        H_matrix(2,1) = cmplx(0.5_real64, -0.1_real64, real64)  ! Hermitian
        H_matrix(1,3) = cmplx(0.2_real64, 0.3_real64, real64)
        H_matrix(3,1) = cmplx(0.2_real64, -0.3_real64, real64)  ! Hermitian
        H_matrix(2,3) = cmplx(0.0_real64, 0.2_real64, real64)
        H_matrix(3,2) = cmplx(0.0_real64, -0.2_real64, real64)  ! Hermitian
        
        call solve_hermitian_eigenproblem(3, H_matrix, eigenvals, eigenvecs, ierr)
        
        ! Check that eigenvalues are ordered
        test_passed = (ierr == 0)
        if (test_passed) then
            do i = 1, 2
                test_passed = test_passed .and. (eigenvals(i) <= eigenvals(i+1))
            end do
        end if
        
        if (.not. test_passed) then
            print *, "  Hermitian solver failed with ierr =", ierr
        else
            print *, "  Real eigenvalues (ordered):", eigenvals
        end if
        
        call end_test(test_passed)
    end subroutine test_hermitian_eigenproblem
    
    !---------------------------------------------------------------------------
    ! Test 3: General real eigenvalue problem
    !---------------------------------------------------------------------------
    subroutine test_general_real_eigenproblem()
        real(real64) :: A_matrix(3,3), eigenvals_real(3), eigenvals_imag(3)
        real(real64) :: eigenvecs(3,3)
        integer :: ierr
        
        call start_test("General real eigenvalue problem")
        
        ! Create test matrix
        A_matrix = 0.0_real64
        A_matrix(1,1) = 1.0_real64
        A_matrix(2,2) = 2.0_real64
        A_matrix(3,3) = 3.0_real64
        A_matrix(1,2) = 0.5_real64
        A_matrix(2,3) = 0.3_real64
        
        call solve_general_real_eigenproblem(3, A_matrix, eigenvals_real, eigenvals_imag, &
                                            eigenvecs, ierr)
        
        test_passed = (ierr == 0)
        
        if (.not. test_passed) then
            print *, "  Real general solver failed with ierr =", ierr
        else
            print *, "  Real parts:", eigenvals_real
            print *, "  Imaginary parts:", eigenvals_imag
        end if
        
        call end_test(test_passed)
    end subroutine test_general_real_eigenproblem
    
    !---------------------------------------------------------------------------
    ! Test 4: Symmetric real eigenvalue problem
    !---------------------------------------------------------------------------
    subroutine test_symmetric_real_eigenproblem()
        real(real64) :: S_matrix(3,3), eigenvals(3), eigenvecs(3,3)
        integer :: ierr, i
        
        call start_test("Symmetric real eigenvalue problem")
        
        ! Create symmetric matrix
        S_matrix(1,1) = 2.0_real64
        S_matrix(2,2) = 3.0_real64
        S_matrix(3,3) = 1.0_real64
        S_matrix(1,2) = 0.5_real64
        S_matrix(2,1) = 0.5_real64  ! Symmetric
        S_matrix(1,3) = 0.3_real64
        S_matrix(3,1) = 0.3_real64  ! Symmetric
        S_matrix(2,3) = 0.4_real64
        S_matrix(3,2) = 0.4_real64  ! Symmetric
        
        call solve_symmetric_real_eigenproblem(3, S_matrix, eigenvals, eigenvecs, ierr)
        
        ! Check that eigenvalues are ordered
        test_passed = (ierr == 0)
        if (test_passed) then
            do i = 1, 2
                test_passed = test_passed .and. (eigenvals(i) <= eigenvals(i+1))
            end do
        end if
        
        if (.not. test_passed) then
            print *, "  Symmetric solver failed with ierr =", ierr
        else
            print *, "  Eigenvalues (ordered):", eigenvals
        end if
        
        call end_test(test_passed)
    end subroutine test_symmetric_real_eigenproblem
    
    !---------------------------------------------------------------------------
    ! Test 5: Mode eigenvalue problem with sorting
    !---------------------------------------------------------------------------
    subroutine test_mode_eigenvalue_problem()
        complex(real64) :: mode_matrix(3,3)
        complex(real64), allocatable :: frequencies(:), eigenvectors(:,:)
        real(real64), allocatable :: growth_rates(:)
        integer :: ierr
        
        call start_test("Mode eigenvalue problem with growth rate sorting")
        
        ! Create mode matrix with different growth rates
        mode_matrix = cmplx(0.0_real64, 0.0_real64, real64)
        mode_matrix(1,1) = cmplx(1.0_real64, 0.5_real64, real64)  ! Growth rate = 0.5
        mode_matrix(2,2) = cmplx(2.0_real64, -0.3_real64, real64) ! Growth rate = -0.3 (stable)
        mode_matrix(3,3) = cmplx(3.0_real64, 1.2_real64, real64)  ! Growth rate = 1.2 (most unstable)
        
        call solve_mode_eigenvalue_problem(mode_matrix, frequencies, eigenvectors, &
                                          growth_rates, ierr)
        
        ! Check that growth rates are sorted in descending order
        test_passed = (ierr == 0) .and. allocated(frequencies) .and. &
                      allocated(eigenvectors) .and. allocated(growth_rates)
        
        if (test_passed) then
            ! Most unstable mode should be first
            test_passed = (growth_rates(1) >= growth_rates(2)) .and. &
                         (growth_rates(2) >= growth_rates(3))
        end if
        
        if (.not. test_passed) then
            print *, "  Mode solver failed with ierr =", ierr
        else
            print *, "  Growth rates (sorted):", growth_rates
            print *, "  Frequencies:", frequencies
        end if
        
        if (allocated(frequencies)) deallocate(frequencies)
        if (allocated(eigenvectors)) deallocate(eigenvectors)
        if (allocated(growth_rates)) deallocate(growth_rates)
        
        call end_test(test_passed)
    end subroutine test_mode_eigenvalue_problem
    
    !---------------------------------------------------------------------------
    ! Test 6: Growth rate sorting functionality
    !---------------------------------------------------------------------------
    subroutine test_growth_rate_sorting()
        complex(real64) :: frequencies(4), eigenvectors(4,4)
        real(real64) :: growth_rates(4)
        integer :: i
        
        call start_test("Growth rate sorting functionality")
        
        ! Create test data with unsorted growth rates
        frequencies(1) = cmplx(1.0_real64, 0.2_real64, real64)   ! Growth rate = 0.2
        frequencies(2) = cmplx(2.0_real64, -0.5_real64, real64)  ! Growth rate = -0.5
        frequencies(3) = cmplx(3.0_real64, 1.5_real64, real64)   ! Growth rate = 1.5 (highest)
        frequencies(4) = cmplx(4.0_real64, 0.8_real64, real64)   ! Growth rate = 0.8
        
        growth_rates = aimag(frequencies)
        
        ! Initialize eigenvectors with identity-like pattern for testing
        eigenvectors = cmplx(0.0_real64, 0.0_real64, real64)
        do i = 1, 4
            eigenvectors(i,i) = cmplx(real(i, real64), 0.0_real64, real64)
        end do
        
        call sort_eigenvalues_by_growth_rate(frequencies, eigenvectors, growth_rates)
        
        ! Check sorting (descending order: 1.5, 0.8, 0.2, -0.5)
        test_passed = (abs(growth_rates(1) - 1.5_real64) < 1e-10_real64) .and. &
                      (abs(growth_rates(2) - 0.8_real64) < 1e-10_real64) .and. &
                      (abs(growth_rates(3) - 0.2_real64) < 1e-10_real64) .and. &
                      (abs(growth_rates(4) - (-0.5_real64)) < 1e-10_real64)
        
        if (.not. test_passed) then
            print *, "  Sorting failed. Growth rates:", growth_rates
        else
            print *, "  Growth rates correctly sorted:", growth_rates
        end if
        
        call end_test(test_passed)
    end subroutine test_growth_rate_sorting
    
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

end program test_kilca_eigensolver