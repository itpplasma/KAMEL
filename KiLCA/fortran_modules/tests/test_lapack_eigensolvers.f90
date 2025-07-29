program test_lapack_eigensolvers
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
    print *, "LAPACK Eigensolver Usage Tests"
    print *, "===================================================="
    print *, ""
    
    ! Test 1: ZGEEV for general complex eigenvalue problems
    call test_zgeev_general()
    
    ! Test 2: ZHEEV for Hermitian eigenvalue problems
    call test_zheev_hermitian()
    
    ! Test 3: DGEEV for general real eigenvalue problems
    call test_dgeev_general()
    
    ! Test 4: DSYEV for symmetric real eigenvalue problems
    call test_dsyev_symmetric()
    
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
    ! Test 1: ZGEEV for general complex matrices
    !---------------------------------------------------------------------------
    subroutine test_zgeev_general()
        complex(real64), allocatable :: A(:,:), VL(:,:), VR(:,:), W(:)
        complex(real64), allocatable :: work(:)
        real(real64), allocatable :: rwork(:)
        integer :: n, lwork, info
        complex(real64) :: expected_eval
        real(real64) :: tol
        
        call start_test("ZGEEV for general complex eigenvalue problem")
        
        n = 3
        tol = 1.0e-10_real64
        
        ! Allocate arrays
        allocate(A(n,n), VL(n,n), VR(n,n), W(n))
        allocate(rwork(2*n))
        
        ! Create test matrix with known eigenvalues
        ! Matrix with eigenvalues 1, 2, 3
        A = cmplx(0.0_real64, 0.0_real64, real64)
        A(1,1) = cmplx(1.0_real64, 0.0_real64, real64)
        A(2,2) = cmplx(2.0_real64, 0.0_real64, real64)
        A(3,3) = cmplx(3.0_real64, 0.0_real64, real64)
        ! Add some off-diagonal elements
        A(1,2) = cmplx(0.5_real64, 0.1_real64, real64)
        A(2,1) = cmplx(0.5_real64, -0.1_real64, real64)
        
        ! Query optimal workspace
        lwork = -1
        allocate(work(1))
        call zgeev('N', 'V', n, A, n, W, VL, n, VR, n, work, lwork, rwork, info)
        lwork = int(real(work(1)))
        deallocate(work)
        allocate(work(lwork))
        
        ! Compute eigenvalues and right eigenvectors
        call zgeev('N', 'V', n, A, n, W, VL, n, VR, n, work, lwork, rwork, info)
        
        test_passed = (info == 0)
        
        if (.not. test_passed) then
            print *, "  ZGEEV failed with info =", info
        else
            print *, "  Eigenvalues found:", W
        end if
        
        deallocate(A, VL, VR, W, work, rwork)
        call end_test(test_passed)
    end subroutine test_zgeev_general
    
    !---------------------------------------------------------------------------
    ! Test 2: ZHEEV for Hermitian matrices
    !---------------------------------------------------------------------------
    subroutine test_zheev_hermitian()
        complex(real64), allocatable :: A(:,:), work(:)
        real(real64), allocatable :: W(:), rwork(:)
        integer :: n, lwork, info
        real(real64) :: tol
        
        call start_test("ZHEEV for Hermitian eigenvalue problem")
        
        n = 3
        tol = 1.0e-10_real64
        
        ! Allocate arrays
        allocate(A(n,n), W(n))
        allocate(rwork(3*n-2))
        
        ! Create Hermitian test matrix
        A(1,1) = cmplx(1.0_real64, 0.0_real64, real64)
        A(2,2) = cmplx(2.0_real64, 0.0_real64, real64)
        A(3,3) = cmplx(3.0_real64, 0.0_real64, real64)
        A(1,2) = cmplx(0.5_real64, 0.1_real64, real64)
        A(2,1) = cmplx(0.5_real64, -0.1_real64, real64)  ! Hermitian
        A(1,3) = cmplx(0.0_real64, 0.0_real64, real64)
        A(3,1) = cmplx(0.0_real64, 0.0_real64, real64)
        A(2,3) = cmplx(0.3_real64, 0.2_real64, real64)
        A(3,2) = cmplx(0.3_real64, -0.2_real64, real64)  ! Hermitian
        
        ! Query optimal workspace
        lwork = -1
        allocate(work(1))
        call zheev('V', 'U', n, A, n, W, work, lwork, rwork, info)
        lwork = int(real(work(1)))
        deallocate(work)
        allocate(work(lwork))
        
        ! Compute eigenvalues and eigenvectors
        call zheev('V', 'U', n, A, n, W, work, lwork, rwork, info)
        
        test_passed = (info == 0) .and. (W(1) <= W(2)) .and. (W(2) <= W(3))
        
        if (.not. test_passed) then
            print *, "  ZHEEV failed with info =", info
        else
            print *, "  Eigenvalues (should be real and ordered):", W
        end if
        
        deallocate(A, W, work, rwork)
        call end_test(test_passed)
    end subroutine test_zheev_hermitian
    
    !---------------------------------------------------------------------------
    ! Test 3: DGEEV for general real matrices
    !---------------------------------------------------------------------------
    subroutine test_dgeev_general()
        real(real64), allocatable :: A(:,:), VL(:,:), VR(:,:)
        real(real64), allocatable :: WR(:), WI(:), work(:)
        integer :: n, lwork, info
        
        call start_test("DGEEV for general real eigenvalue problem")
        
        n = 3
        
        ! Allocate arrays
        allocate(A(n,n), VL(n,n), VR(n,n))
        allocate(WR(n), WI(n))
        
        ! Create test matrix
        A = 0.0_real64
        A(1,1) = 1.0_real64
        A(2,2) = 2.0_real64
        A(3,3) = 3.0_real64
        A(1,2) = 0.5_real64
        A(2,1) = 0.5_real64
        
        ! Query optimal workspace
        lwork = -1
        allocate(work(1))
        call dgeev('N', 'V', n, A, n, WR, WI, VL, n, VR, n, work, lwork, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        
        ! Compute eigenvalues and right eigenvectors
        call dgeev('N', 'V', n, A, n, WR, WI, VL, n, VR, n, work, lwork, info)
        
        test_passed = (info == 0)
        
        if (.not. test_passed) then
            print *, "  DGEEV failed with info =", info
        else
            print *, "  Real parts:", WR
            print *, "  Imaginary parts:", WI
        end if
        
        deallocate(A, VL, VR, WR, WI, work)
        call end_test(test_passed)
    end subroutine test_dgeev_general
    
    !---------------------------------------------------------------------------
    ! Test 4: DSYEV for symmetric real matrices
    !---------------------------------------------------------------------------
    subroutine test_dsyev_symmetric()
        real(real64), allocatable :: A(:,:), W(:), work(:)
        integer :: n, lwork, info
        
        call start_test("DSYEV for symmetric eigenvalue problem")
        
        n = 3
        
        ! Allocate arrays
        allocate(A(n,n), W(n))
        
        ! Create symmetric test matrix
        A(1,1) = 1.0_real64
        A(2,2) = 2.0_real64
        A(3,3) = 3.0_real64
        A(1,2) = 0.5_real64
        A(2,1) = 0.5_real64  ! Symmetric
        A(1,3) = 0.3_real64
        A(3,1) = 0.3_real64  ! Symmetric
        A(2,3) = 0.4_real64
        A(3,2) = 0.4_real64  ! Symmetric
        
        ! Query optimal workspace
        lwork = -1
        allocate(work(1))
        call dsyev('V', 'U', n, A, n, W, work, lwork, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        
        ! Compute eigenvalues and eigenvectors
        call dsyev('V', 'U', n, A, n, W, work, lwork, info)
        
        test_passed = (info == 0) .and. (W(1) <= W(2)) .and. (W(2) <= W(3))
        
        if (.not. test_passed) then
            print *, "  DSYEV failed with info =", info
        else
            print *, "  Eigenvalues (should be ordered):", W
        end if
        
        deallocate(A, W, work)
        call end_test(test_passed)
    end subroutine test_dsyev_symmetric
    
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

end program test_lapack_eigensolvers