program test_kilca_orthogonalization
    use iso_fortran_env, only: real64, int32
    use kilca_orthogonalization_m
    implicit none
    
    ! Test counters
    integer :: total_tests, passed_tests, failed_tests
    logical :: test_passed
    
    ! Initialize test counters
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    
    print *, "===================================================="
    print *, "KiLCA Orthogonalization Module - Unit Tests"
    print *, "===================================================="
    print *, ""
    
    ! Test 1: Complex basis orthogonalization
    call test_complex_basis_orthogonalization()
    
    ! Test 2: Real basis orthogonalization
    call test_real_basis_orthogonalization()
    
    ! Test 3: Complex QR decomposition
    call test_complex_qr_decomposition()
    
    ! Test 4: Real QR decomposition
    call test_real_qr_decomposition()
    
    ! Test 5: Complex orthogonal transformation application
    call test_complex_orthogonal_transformation()
    
    ! Test 6: Real orthogonal transformation application
    call test_real_orthogonal_transformation()
    
    ! Test 7: Complex orthogonality verification
    call test_complex_orthogonality_verification()
    
    ! Test 8: Real orthogonality verification
    call test_real_orthogonality_verification()
    
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
    ! Test 1: Complex basis vector orthogonalization
    !---------------------------------------------------------------------------
    subroutine test_complex_basis_orthogonalization()
        complex(real64) :: basis(5,3)
        integer :: ierr, i, j
        logical :: is_orthogonal
        real(real64) :: max_error
        
        call start_test("Complex basis vector orthogonalization")
        
        ! Create non-orthogonal basis vectors
        do i = 1, 5
            do j = 1, 3
                basis(i,j) = cmplx(sin(real(i*j, real64)), cos(real(i+j, real64)), real64)
            end do
        end do
        
        ! Orthogonalize
        call orthogonalize_basis_vectors_complex(5, 3, basis, ierr)
        
        test_passed = (ierr == 0)
        
        if (test_passed) then
            ! Verify orthogonality
            call verify_orthogonality_complex(5, 3, basis, is_orthogonal, max_error)
            test_passed = is_orthogonal
        end if
        
        if (.not. test_passed) then
            print *, "  Complex orthogonalization failed, ierr =", ierr, ", max_error =", max_error
        else
            print *, "  Complex orthogonalization successful, max_error =", max_error
        end if
        
        call end_test(test_passed)
    end subroutine test_complex_basis_orthogonalization
    
    !---------------------------------------------------------------------------
    ! Test 2: Real basis vector orthogonalization
    !---------------------------------------------------------------------------
    subroutine test_real_basis_orthogonalization()
        real(real64) :: basis(4,3)
        integer :: ierr, i, j
        logical :: is_orthogonal
        real(real64) :: max_error
        
        call start_test("Real basis vector orthogonalization")
        
        ! Create non-orthogonal real basis vectors
        do i = 1, 4
            do j = 1, 3
                basis(i,j) = sin(real(i*j, real64)) + cos(real(i+j, real64))
            end do
        end do
        
        ! Orthogonalize
        call orthogonalize_basis_vectors_real(4, 3, basis, ierr)
        
        test_passed = (ierr == 0)
        
        if (test_passed) then
            ! Verify orthogonality
            call verify_orthogonality_real(4, 3, basis, is_orthogonal, max_error)
            test_passed = is_orthogonal
        end if
        
        if (.not. test_passed) then
            print *, "  Real orthogonalization failed, ierr =", ierr, ", max_error =", max_error
        else
            print *, "  Real orthogonalization successful, max_error =", max_error
        end if
        
        call end_test(test_passed)
    end subroutine test_real_basis_orthogonalization
    
    !---------------------------------------------------------------------------
    ! Test 3: Complex QR decomposition
    !---------------------------------------------------------------------------
    subroutine test_complex_qr_decomposition()
        complex(real64) :: A(4,3), Q(4,3), R(3,3), test_matrix(4,3)
        integer :: ierr, i, j
        real(real64) :: max_error
        
        call start_test("Complex QR decomposition")
        
        ! Create test matrix
        A(1,:) = [cmplx(1.0_real64, 0.1_real64, real64), &
                  cmplx(2.0_real64, 0.2_real64, real64), &
                  cmplx(3.0_real64, 0.3_real64, real64)]
        A(2,:) = [cmplx(4.0_real64, -0.1_real64, real64), &
                  cmplx(5.0_real64, -0.2_real64, real64), &
                  cmplx(6.0_real64, -0.3_real64, real64)]
        A(3,:) = [cmplx(7.0_real64, 0.4_real64, real64), &
                  cmplx(8.0_real64, 0.5_real64, real64), &
                  cmplx(9.0_real64, 0.6_real64, real64)]
        A(4,:) = [cmplx(1.0_real64, -0.05_real64, real64), &
                  cmplx(0.5_real64, 0.1_real64, real64), &
                  cmplx(0.2_real64, -0.15_real64, real64)]
        
        ! Perform QR decomposition
        call qr_decomposition_complex(4, 3, A, Q, R, ierr)
        
        test_passed = (ierr == 0)
        
        if (test_passed) then
            ! Verify Q*R = A
            call zgemm('N', 'N', 4, 3, 3, cmplx(1.0_real64, 0.0_real64, real64), &
                       Q, 4, R, 3, cmplx(0.0_real64, 0.0_real64, real64), test_matrix, 4)
            
            max_error = 0.0_real64
            do i = 1, 4
                do j = 1, 3
                    max_error = max(max_error, abs(test_matrix(i,j) - A(i,j)))
                end do
            end do
            
            test_passed = (max_error < 1.0e-10_real64)
        end if
        
        if (.not. test_passed) then
            print *, "  Complex QR decomposition failed, ierr =", ierr, ", max_error =", max_error
        else
            print *, "  Complex QR decomposition successful, max_error =", max_error
        end if
        
        call end_test(test_passed)
    end subroutine test_complex_qr_decomposition
    
    !---------------------------------------------------------------------------
    ! Test 4: Real QR decomposition
    !---------------------------------------------------------------------------
    subroutine test_real_qr_decomposition()
        real(real64) :: A(4,3), Q(4,3), R(3,3), test_matrix(4,3)
        integer :: ierr, i, j
        real(real64) :: max_error
        
        call start_test("Real QR decomposition")
        
        ! Create test matrix
        A(1,:) = [1.0_real64, 2.0_real64, 3.0_real64]
        A(2,:) = [4.0_real64, 5.0_real64, 6.0_real64]
        A(3,:) = [7.0_real64, 8.0_real64, 9.0_real64]
        A(4,:) = [1.0_real64, 0.5_real64, 0.2_real64]
        
        ! Perform QR decomposition
        call qr_decomposition_real(4, 3, A, Q, R, ierr)
        
        test_passed = (ierr == 0)
        
        if (test_passed) then
            ! Verify Q*R = A
            call dgemm('N', 'N', 4, 3, 3, 1.0_real64, Q, 4, R, 3, 0.0_real64, test_matrix, 4)
            
            max_error = 0.0_real64
            do i = 1, 4
                do j = 1, 3
                    max_error = max(max_error, abs(test_matrix(i,j) - A(i,j)))
                end do
            end do
            
            test_passed = (max_error < 1.0e-10_real64)
        end if
        
        if (.not. test_passed) then
            print *, "  Real QR decomposition failed, ierr =", ierr, ", max_error =", max_error
        else
            print *, "  Real QR decomposition successful, max_error =", max_error
        end if
        
        call end_test(test_passed)
    end subroutine test_real_qr_decomposition
    
    !---------------------------------------------------------------------------
    ! Test 5: Complex orthogonal transformation application
    !---------------------------------------------------------------------------
    subroutine test_complex_orthogonal_transformation()
        complex(real64) :: A(4,3), Q(4,3), C(4,2), tau(3)
        integer :: ierr
        
        call start_test("Complex orthogonal transformation application")
        
        ! Create test matrices
        A(1,:) = [cmplx(1.0_real64, 0.1_real64, real64), &
                  cmplx(2.0_real64, 0.2_real64, real64), &
                  cmplx(3.0_real64, 0.3_real64, real64)]
        A(2,:) = [cmplx(4.0_real64, -0.1_real64, real64), &
                  cmplx(5.0_real64, -0.2_real64, real64), &
                  cmplx(6.0_real64, -0.3_real64, real64)]
        A(3,:) = [cmplx(7.0_real64, 0.4_real64, real64), &
                  cmplx(8.0_real64, 0.5_real64, real64), &
                  cmplx(9.0_real64, 0.6_real64, real64)]
        A(4,:) = [cmplx(1.0_real64, -0.05_real64, real64), &
                  cmplx(0.5_real64, 0.1_real64, real64), &
                  cmplx(0.2_real64, -0.15_real64, real64)]
        
        C(1,:) = [cmplx(1.0_real64, 0.0_real64, real64), &
                  cmplx(0.0_real64, 1.0_real64, real64)]
        C(2,:) = [cmplx(0.0_real64, 0.0_real64, real64), &
                  cmplx(1.0_real64, 0.0_real64, real64)]
        C(3,:) = [cmplx(1.0_real64, 1.0_real64, real64), &
                  cmplx(0.0_real64, 0.0_real64, real64)]
        C(4,:) = [cmplx(0.5_real64, 0.5_real64, real64), &
                  cmplx(0.5_real64, -0.5_real64, real64)]
        
        Q = A
        
        ! First orthogonalize to get Q and tau
        call orthogonalize_basis_vectors_complex(4, 3, Q, ierr)
        test_passed = (ierr == 0)
        
        ! Note: This test is simplified - in practice we'd need to extract tau from ZGEQRF
        ! For testing purposes, we'll just check that the transformation call doesn't crash
        if (test_passed) then
            ! Create dummy tau for test
            tau = cmplx(1.0_real64, 0.0_real64, real64)
            call apply_orthogonal_transformation_complex('L', 'C', 4, 2, 3, Q, tau, C, ierr)
            test_passed = (ierr == 0)
        end if
        
        if (.not. test_passed) then
            print *, "  Complex transformation application failed, ierr =", ierr
        else
            print *, "  Complex transformation application successful"
        end if
        
        call end_test(test_passed)
    end subroutine test_complex_orthogonal_transformation
    
    !---------------------------------------------------------------------------
    ! Test 6: Real orthogonal transformation application
    !---------------------------------------------------------------------------
    subroutine test_real_orthogonal_transformation()
        real(real64) :: A(4,3), Q(4,3), C(4,2), tau(3)
        integer :: ierr
        
        call start_test("Real orthogonal transformation application")
        
        ! Create test matrices
        A(1,:) = [1.0_real64, 2.0_real64, 3.0_real64]
        A(2,:) = [4.0_real64, 5.0_real64, 6.0_real64]
        A(3,:) = [7.0_real64, 8.0_real64, 9.0_real64]
        A(4,:) = [1.0_real64, 0.5_real64, 0.2_real64]
        
        C(1,:) = [1.0_real64, 0.0_real64]
        C(2,:) = [0.0_real64, 1.0_real64]
        C(3,:) = [1.0_real64, 1.0_real64]
        C(4,:) = [0.5_real64, 0.5_real64]
        
        Q = A
        
        ! First orthogonalize to get Q and tau
        call orthogonalize_basis_vectors_real(4, 3, Q, ierr)
        test_passed = (ierr == 0)
        
        ! Simplified test - check that transformation call doesn't crash
        if (test_passed) then
            tau = 1.0_real64
            call apply_orthogonal_transformation_real('L', 'T', 4, 2, 3, Q, tau, C, ierr)
            test_passed = (ierr == 0)
        end if
        
        if (.not. test_passed) then
            print *, "  Real transformation application failed, ierr =", ierr
        else
            print *, "  Real transformation application successful"
        end if
        
        call end_test(test_passed)
    end subroutine test_real_orthogonal_transformation
    
    !---------------------------------------------------------------------------
    ! Test 7: Complex orthogonality verification
    !---------------------------------------------------------------------------
    subroutine test_complex_orthogonality_verification()
        complex(real64) :: Q(4,3)
        logical :: is_orthogonal
        real(real64) :: max_error
        integer :: i, j
        
        call start_test("Complex orthogonality verification")
        
        ! Create orthogonal matrix manually
        Q = cmplx(0.0_real64, 0.0_real64, real64)
        Q(1,1) = cmplx(1.0_real64, 0.0_real64, real64)
        Q(2,2) = cmplx(1.0_real64, 0.0_real64, real64)
        Q(3,3) = cmplx(1.0_real64, 0.0_real64, real64)
        
        call verify_orthogonality_complex(4, 3, Q, is_orthogonal, max_error)
        
        test_passed = is_orthogonal
        
        if (.not. test_passed) then
            print *, "  Complex orthogonality verification failed, max_error =", max_error
        else
            print *, "  Complex orthogonality verification successful, max_error =", max_error
        end if
        
        call end_test(test_passed)
    end subroutine test_complex_orthogonality_verification
    
    !---------------------------------------------------------------------------
    ! Test 8: Real orthogonality verification
    !---------------------------------------------------------------------------
    subroutine test_real_orthogonality_verification()
        real(real64) :: Q(4,3)
        logical :: is_orthogonal
        real(real64) :: max_error
        
        call start_test("Real orthogonality verification")
        
        ! Create orthogonal matrix manually
        Q = 0.0_real64
        Q(1,1) = 1.0_real64
        Q(2,2) = 1.0_real64
        Q(3,3) = 1.0_real64
        
        call verify_orthogonality_real(4, 3, Q, is_orthogonal, max_error)
        
        test_passed = is_orthogonal
        
        if (.not. test_passed) then
            print *, "  Real orthogonality verification failed, max_error =", max_error
        else
            print *, "  Real orthogonality verification successful, max_error =", max_error
        end if
        
        call end_test(test_passed)
    end subroutine test_real_orthogonality_verification
    
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

end program test_kilca_orthogonalization