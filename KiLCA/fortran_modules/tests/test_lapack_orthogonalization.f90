program test_lapack_orthogonalization
    use iso_fortran_env, only: real64, int32, error_unit
    implicit none
    
    ! Test counters
    integer :: total_tests, passed_tests, failed_tests
    logical :: test_passed
    
    ! Initialize test counters
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    
    print *, "===================================================="
    print *, "LAPACK Orthogonalization Tests"
    print *, "===================================================="
    print *, ""
    
    ! Test 1: ZGEQRF + ZUNGQR for complex QR decomposition
    call test_zgeqrf_zungqr()
    
    ! Test 2: DGEQRF + DORGQR for real QR decomposition
    call test_dgeqrf_dorgqr()
    
    ! Test 3: ZUNMQR for applying Q to another matrix
    call test_zunmqr_application()
    
    ! Test 4: Orthogonality verification
    call test_orthogonality_check()
    
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
    ! Test 1: Complex QR decomposition using ZGEQRF + ZUNGQR
    !---------------------------------------------------------------------------
    subroutine test_zgeqrf_zungqr()
        complex(real64), allocatable :: A(:,:), Q(:,:), R(:,:), tau(:), work(:)
        integer :: m, n, k, lwork, info
        real(real64) :: tol
        
        call start_test("ZGEQRF + ZUNGQR complex QR decomposition")
        
        m = 4
        n = 3
        k = min(m, n)
        tol = 1.0e-12_real64
        
        ! Allocate arrays
        allocate(A(m,n), Q(m,n), R(n,n), tau(k))
        
        ! Create test matrix
        A(1,1) = cmplx(1.0_real64, 0.5_real64, real64)
        A(2,1) = cmplx(2.0_real64, -0.3_real64, real64)
        A(3,1) = cmplx(1.5_real64, 0.8_real64, real64)
        A(4,1) = cmplx(0.5_real64, 0.2_real64, real64)
        
        A(1,2) = cmplx(0.8_real64, -0.2_real64, real64)
        A(2,2) = cmplx(1.2_real64, 0.6_real64, real64)
        A(3,2) = cmplx(0.9_real64, -0.4_real64, real64)
        A(4,2) = cmplx(1.1_real64, 0.3_real64, real64)
        
        A(1,3) = cmplx(0.3_real64, 0.7_real64, real64)
        A(2,3) = cmplx(0.6_real64, -0.1_real64, real64)
        A(3,3) = cmplx(0.4_real64, 0.5_real64, real64)
        A(4,3) = cmplx(0.7_real64, -0.3_real64, real64)
        
        ! Store original A for verification
        Q = A
        
        ! Query optimal workspace for ZGEQRF
        lwork = -1
        allocate(work(1))
        call zgeqrf(m, n, Q, m, tau, work, lwork, info)
        lwork = int(real(work(1)))
        deallocate(work)
        allocate(work(lwork))
        
        ! Perform QR factorization
        call zgeqrf(m, n, Q, m, tau, work, lwork, info)
        
        test_passed = (info == 0)
        
        if (test_passed) then
            ! Generate Q matrix using ZUNGQR
            call zungqr(m, n, k, Q, m, tau, work, lwork, info)
            test_passed = (info == 0)
        end if
        
        if (.not. test_passed) then
            print *, "  QR decomposition failed with info =", info
        else
            print *, "  QR decomposition successful"
        end if
        
        deallocate(A, Q, R, tau, work)
        call end_test(test_passed)
    end subroutine test_zgeqrf_zungqr
    
    !---------------------------------------------------------------------------
    ! Test 2: Real QR decomposition using DGEQRF + DORGQR
    !---------------------------------------------------------------------------
    subroutine test_dgeqrf_dorgqr()
        real(real64), allocatable :: A(:,:), Q(:,:), tau(:), work(:)
        integer :: m, n, k, lwork, info
        
        call start_test("DGEQRF + DORGQR real QR decomposition")
        
        m = 4
        n = 3
        k = min(m, n)
        
        ! Allocate arrays
        allocate(A(m,n), Q(m,n), tau(k))
        
        ! Create test matrix
        A(1,:) = [1.0_real64, 2.0_real64, 3.0_real64]
        A(2,:) = [4.0_real64, 5.0_real64, 6.0_real64]
        A(3,:) = [7.0_real64, 8.0_real64, 9.0_real64]
        A(4,:) = [1.0_real64, 0.5_real64, 0.2_real64]
        
        ! Copy to Q
        Q = A
        
        ! Query optimal workspace
        lwork = -1
        allocate(work(1))
        call dgeqrf(m, n, Q, m, tau, work, lwork, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        
        ! Perform QR factorization
        call dgeqrf(m, n, Q, m, tau, work, lwork, info)
        
        test_passed = (info == 0)
        
        if (test_passed) then
            ! Generate Q matrix
            call dorgqr(m, n, k, Q, m, tau, work, lwork, info)
            test_passed = (info == 0)
        end if
        
        if (.not. test_passed) then
            print *, "  Real QR decomposition failed with info =", info
        else
            print *, "  Real QR decomposition successful"
        end if
        
        deallocate(A, Q, tau, work)
        call end_test(test_passed)
    end subroutine test_dgeqrf_dorgqr
    
    !---------------------------------------------------------------------------
    ! Test 3: Apply Q matrix to another matrix using ZUNMQR
    !---------------------------------------------------------------------------
    subroutine test_zunmqr_application()
        complex(real64), allocatable :: A(:,:), Q(:,:), C(:,:), tau(:), work(:)
        integer :: m, n, k, nrhs, lwork, info
        
        call start_test("ZUNMQR application of Q to matrix")
        
        m = 4
        n = 3
        k = min(m, n)
        nrhs = 2
        
        ! Allocate arrays
        allocate(A(m,n), Q(m,n), C(m,nrhs), tau(k))
        
        ! Create test matrices
        A(1,:) = [cmplx(1.0_real64, 0.1_real64, real64), &
                  cmplx(2.0_real64, 0.2_real64, real64), &
                  cmplx(3.0_real64, 0.3_real64, real64)]
        A(2,:) = [cmplx(4.0_real64, 0.4_real64, real64), &
                  cmplx(5.0_real64, 0.5_real64, real64), &
                  cmplx(6.0_real64, 0.6_real64, real64)]
        A(3,:) = [cmplx(7.0_real64, 0.7_real64, real64), &
                  cmplx(8.0_real64, 0.8_real64, real64), &
                  cmplx(9.0_real64, 0.9_real64, real64)]
        A(4,:) = [cmplx(1.0_real64, 0.01_real64, real64), &
                  cmplx(0.5_real64, 0.05_real64, real64), &
                  cmplx(0.2_real64, 0.02_real64, real64)]
        
        C(1,:) = [cmplx(1.0_real64, 0.0_real64, real64), &
                  cmplx(0.0_real64, 1.0_real64, real64)]
        C(2,:) = [cmplx(0.0_real64, 0.0_real64, real64), &
                  cmplx(1.0_real64, 0.0_real64, real64)]
        C(3,:) = [cmplx(1.0_real64, 1.0_real64, real64), &
                  cmplx(0.0_real64, 0.0_real64, real64)]
        C(4,:) = [cmplx(0.5_real64, 0.5_real64, real64), &
                  cmplx(0.5_real64, -0.5_real64, real64)]
        
        Q = A
        
        ! Query workspace and perform QR factorization
        lwork = -1
        allocate(work(1))
        call zgeqrf(m, n, Q, m, tau, work, lwork, info)
        lwork = int(real(work(1)))
        deallocate(work)
        allocate(work(lwork))
        
        call zgeqrf(m, n, Q, m, tau, work, lwork, info)
        test_passed = (info == 0)
        
        if (test_passed) then
            ! Apply Q^H to C using ZUNMQR
            call zunmqr('L', 'C', m, nrhs, k, Q, m, tau, C, m, work, lwork, info)
            test_passed = (info == 0)
        end if
        
        if (.not. test_passed) then
            print *, "  ZUNMQR application failed with info =", info
        else
            print *, "  ZUNMQR application successful"
        end if
        
        deallocate(A, Q, C, tau, work)
        call end_test(test_passed)
    end subroutine test_zunmqr_application
    
    !---------------------------------------------------------------------------
    ! Test 4: Verify orthogonality of Q matrix
    !---------------------------------------------------------------------------
    subroutine test_orthogonality_check()
        complex(real64), allocatable :: A(:,:), Q(:,:), QtQ(:,:), tau(:), work(:)
        integer :: m, n, k, lwork, info, i, j
        real(real64) :: tol, max_error
        
        call start_test("Orthogonality verification of Q matrix")
        
        m = 5
        n = 3
        k = min(m, n)
        tol = 1.0e-10_real64
        
        ! Allocate arrays
        allocate(A(m,n), Q(m,n), QtQ(n,n), tau(k))
        
        ! Create well-conditioned test matrix
        do i = 1, m
            do j = 1, n
                A(i,j) = cmplx(sin(real(i*j, real64)), cos(real(i+j, real64)), real64)
            end do
        end do
        
        Q = A
        
        ! QR factorization
        lwork = -1
        allocate(work(1))
        call zgeqrf(m, n, Q, m, tau, work, lwork, info)
        lwork = int(real(work(1)))
        deallocate(work)
        allocate(work(lwork))
        
        call zgeqrf(m, n, Q, m, tau, work, lwork, info)
        test_passed = (info == 0)
        
        if (test_passed) then
            ! Generate orthogonal Q
            call zungqr(m, n, k, Q, m, tau, work, lwork, info)
            test_passed = (info == 0)
        end if
        
        if (test_passed) then
            ! Compute Q^H * Q (should be identity)
            call zgemm('C', 'N', n, n, m, cmplx(1.0_real64, 0.0_real64, real64), &
                       Q, m, Q, m, cmplx(0.0_real64, 0.0_real64, real64), QtQ, n)
            
            ! Check if Q^H * Q is identity
            max_error = 0.0_real64
            do i = 1, n
                do j = 1, n
                    if (i == j) then
                        max_error = max(max_error, abs(QtQ(i,j) - cmplx(1.0_real64, 0.0_real64, real64)))
                    else
                        max_error = max(max_error, abs(QtQ(i,j)))
                    end if
                end do
            end do
            
            test_passed = (max_error < tol)
        end if
        
        if (.not. test_passed) then
            print *, "  Orthogonality check failed, max error =", max_error
        else
            print *, "  Orthogonality verified, max error =", max_error
        end if
        
        deallocate(A, Q, QtQ, tau, work)
        call end_test(test_passed)
    end subroutine test_orthogonality_check
    
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

end program test_lapack_orthogonalization