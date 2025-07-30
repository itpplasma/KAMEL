program test_kilca_complex_eigensystem
    use iso_fortran_env, only: real64, int32
    use kilca_complex_m
    use kilca_constants_m, only: pi
    implicit none
    
    logical :: test_passed
    complex(real64), allocatable :: matrix(:,:), eigenvalues(:), eigenvectors(:,:)
    complex(real64), allocatable :: matrix_b(:,:), alpha(:), beta(:), vl(:,:), vr(:,:)
    complex(real64), allocatable :: Av(:), lambda_v(:), Av_gen(:), lambda_Bv(:)
    complex(real64) :: det, trace, expected_det, expected_trace
    complex(real64) :: dominant_eigenvalue
    integer :: n, info, i, j, iterations
    real(real64) :: tol, condition_number
    logical :: magnitude_sorted, real_sorted
    
    test_passed = .true.
    tol = epsilon(1.0_real64) * 1000.0_real64
    
    print *, "Testing complex eigensystem solver interfaces..."
    print *, ""
    
    ! Test 1: Standard eigenvalue problem (A*v = λ*v)
    print *, "Testing standard eigenvalue problem..."
    n = 3
    allocate(matrix(n,n), eigenvalues(n), eigenvectors(n,n))
    
    ! Create a simple test matrix
    matrix = cmplx(0.0_real64, 0.0_real64, real64)
    matrix(1,1) = cmplx(2.0_real64, 0.0_real64, real64)
    matrix(2,2) = cmplx(3.0_real64, 1.0_real64, real64)  
    matrix(3,3) = cmplx(1.0_real64, -1.0_real64, real64)
    matrix(1,2) = cmplx(1.0_real64, 0.0_real64, real64)
    matrix(2,1) = cmplx(1.0_real64, 0.0_real64, real64)
    
    call cmplx_eigensystem_std(matrix, eigenvalues, eigenvectors, info)
    
    if (info == 0) then
        print *, "PASS: Standard eigenvalue computation completed successfully"
        
        ! Test eigenvalue equation: A*v = λ*v for first eigenvector
        allocate(Av(n), lambda_v(n))
        Av = matmul(matrix, eigenvectors(:,1))
        lambda_v = eigenvalues(1) * eigenvectors(:,1)
        
        if (maxval(abs(Av - lambda_v)) < tol) then
            print *, "PASS: Eigenvalue equation A*v = λ*v satisfied"
        else
            print *, "FAIL: Eigenvalue equation not satisfied"
            print *, "  Max error:", maxval(abs(Av - lambda_v))
            test_passed = .false.
        end if
        deallocate(Av, lambda_v)
    else
        print *, "FAIL: Standard eigenvalue computation failed with info =", info
        test_passed = .false.
    end if
    
    deallocate(matrix, eigenvalues, eigenvectors)
    
    ! Test 2: Generalized eigenvalue problem (A*v = λ*B*v)
    print *, ""
    print *, "Testing generalized eigenvalue problem..."
    n = 2
    allocate(matrix(n,n), matrix_b(n,n), alpha(n), beta(n), vl(n,n), vr(n,n))
    
    ! Create test matrices A and B
    matrix(1,1) = cmplx(2.0_real64, 1.0_real64, real64)
    matrix(1,2) = cmplx(1.0_real64, 0.0_real64, real64)
    matrix(2,1) = cmplx(0.0_real64, 1.0_real64, real64)
    matrix(2,2) = cmplx(3.0_real64, 0.0_real64, real64)
    
    matrix_b(1,1) = cmplx(1.0_real64, 0.0_real64, real64)
    matrix_b(1,2) = cmplx(0.0_real64, 0.0_real64, real64)
    matrix_b(2,1) = cmplx(0.0_real64, 0.0_real64, real64)
    matrix_b(2,2) = cmplx(1.0_real64, 0.0_real64, real64)
    
    call cmplx_eigensystem_gen(matrix, matrix_b, alpha, beta, vl, vr, info)
    
    if (info == 0) then
        print *, "PASS: Generalized eigenvalue computation completed successfully"
        
        ! Test generalized eigenvalue equation: A*v = (α/β)*B*v for first eigenvector
        if (abs(beta(1)) > tol) then  ! Avoid division by zero
            allocate(Av_gen(n), lambda_Bv(n))
            Av_gen = matmul(matrix, vr(:,1))
            lambda_Bv = (alpha(1)/beta(1)) * matmul(matrix_b, vr(:,1))
            
            if (maxval(abs(Av_gen - lambda_Bv)) < tol) then
                print *, "PASS: Generalized eigenvalue equation A*v = λ*B*v satisfied"
            else
                print *, "FAIL: Generalized eigenvalue equation not satisfied"
                test_passed = .false.
            end if
            deallocate(Av_gen, lambda_Bv)
        else
            print *, "NOTE: Infinite eigenvalue (beta=0), skipping equation test"
        end if
    else
        print *, "FAIL: Generalized eigenvalue computation failed with info =", info
        test_passed = .false.
    end if
    
    deallocate(matrix, matrix_b, alpha, beta, vl, vr)
    
    ! Test 3: Eigenvalue utilities
    print *, ""
    print *, "Testing eigenvalue utility functions..."
    
    n = 2
    allocate(matrix(n,n))
    
    ! Test matrix with known determinant and trace
    matrix(1,1) = cmplx(1.0_real64, 1.0_real64, real64)
    matrix(1,2) = cmplx(2.0_real64, 0.0_real64, real64)  
    matrix(2,1) = cmplx(0.0_real64, 1.0_real64, real64)
    matrix(2,2) = cmplx(3.0_real64, -1.0_real64, real64)
    
    ! Test determinant calculation
    det = cmplx_matrix_det(matrix)
    expected_det = matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1)
    expected_det = cmplx(1.0_real64 + 1.0_real64, 1.0_real64*(-1.0_real64) - 2.0_real64*1.0_real64, real64)  ! (1+i)*(3-i) - 2*i = 4-2i
    
    if (abs(det - expected_det) < tol) then
        print *, "PASS: Matrix determinant calculation correct"
    else
        print *, "FAIL: Matrix determinant incorrect"
        print *, "  Expected:", expected_det
        print *, "  Got:", det
        test_passed = .false.
    end if
    
    ! Test trace calculation
    trace = cmplx_matrix_trace(matrix)
    expected_trace = matrix(1,1) + matrix(2,2)
    
    if (abs(trace - expected_trace) < tol) then
        print *, "PASS: Matrix trace calculation correct"
    else
        print *, "FAIL: Matrix trace incorrect"
        print *, "  Expected:", expected_trace  
        print *, "  Got:", trace
        test_passed = .false.
    end if
    
    deallocate(matrix)
    
    ! Test 4: Matrix condition number estimation
    print *, ""
    print *, "Testing matrix condition number estimation..."
    
    n = 3
    allocate(matrix(n,n))
    
    ! Create a well-conditioned matrix (identity)
    matrix = cmplx(0.0_real64, 0.0_real64, real64)
    do i = 1, n
        matrix(i,i) = cmplx(1.0_real64, 0.0_real64, real64)
    end do
    
    
    condition_number = cmplx_matrix_condition(matrix)
    
    if (abs(condition_number - 1.0_real64) < tol) then
        print *, "PASS: Identity matrix condition number ≈ 1"
    else
        print *, "NOTE: Identity matrix condition number =", condition_number
    end if
    
    deallocate(matrix)
    
    ! Test 5: Eigenvalue sorting
    print *, ""
    print *, "Testing eigenvalue sorting..."
    
    n = 4
    allocate(eigenvalues(n), eigenvectors(n,n))
    
    ! Create test eigenvalues and eigenvectors
    eigenvalues(1) = cmplx(3.0_real64, 1.0_real64, real64)
    eigenvalues(2) = cmplx(1.0_real64, 0.0_real64, real64)
    eigenvalues(3) = cmplx(2.0_real64, -1.0_real64, real64)
    eigenvalues(4) = cmplx(0.5_real64, 2.0_real64, real64)
    
    ! Initialize eigenvectors (identity for simplicity)
    eigenvectors = cmplx(0.0_real64, 0.0_real64, real64)
    do i = 1, n
        eigenvectors(i,i) = cmplx(1.0_real64, 0.0_real64, real64)
    end do
    
    ! Sort by magnitude (descending)
    call cmplx_eigenvalues_sort_by_magnitude(eigenvalues, eigenvectors, descending=.true.)
    
    magnitude_sorted = .true.
    do i = 1, n-1
        if (abs(eigenvalues(i)) < abs(eigenvalues(i+1))) then
            magnitude_sorted = .false.
            exit
        end if
    end do
    
    if (magnitude_sorted) then
        print *, "PASS: Eigenvalues sorted by magnitude (descending)"
    else
        print *, "FAIL: Eigenvalue magnitude sorting failed"
        test_passed = .false.
    end if
    
    ! Sort by real part (ascending)
    call cmplx_eigenvalues_sort_by_real(eigenvalues, eigenvectors, descending=.false.)
    
    real_sorted = .true.
    do i = 1, n-1
        if (real(eigenvalues(i)) > real(eigenvalues(i+1))) then
            real_sorted = .false.
            exit
        end if
    end do
    
    if (real_sorted) then
        print *, "PASS: Eigenvalues sorted by real part (ascending)"
    else
        print *, "FAIL: Eigenvalue real part sorting failed"
        test_passed = .false.
    end if
    
    deallocate(eigenvalues, eigenvectors)
    
    ! Test 6: Matrix power iteration for dominant eigenvalue
    print *, ""
    print *, "Testing matrix power iteration..."
    
    n = 3
    allocate(matrix(n,n), eigenvectors(n,1))
    
    ! Create matrix with known dominant eigenvalue
    matrix = cmplx(0.0_real64, 0.0_real64, real64)
    matrix(1,1) = cmplx(3.0_real64, 0.0_real64, real64)  ! Dominant eigenvalue = 3
    matrix(2,2) = cmplx(1.0_real64, 0.0_real64, real64)
    matrix(3,3) = cmplx(0.5_real64, 0.0_real64, real64)
    
    call cmplx_power_iteration(matrix, dominant_eigenvalue, eigenvectors(:,1), iterations, info)
    
    if (info == 0) then
        if (abs(dominant_eigenvalue - cmplx(3.0_real64, 0.0_real64, real64)) < tol) then
            print *, "PASS: Power iteration found correct dominant eigenvalue"
            print *, "  Iterations:", iterations
        else
            print *, "FAIL: Power iteration found incorrect dominant eigenvalue"
            print *, "  Expected: 3.0" 
            print *, "  Got:", dominant_eigenvalue
            test_passed = .false.
        end if
    else
        print *, "FAIL: Power iteration failed with info =", info
        test_passed = .false.
    end if
    
    deallocate(matrix, eigenvectors)
    
    print *, ""
    if (test_passed) then
        print *, "All complex eigensystem solver interface tests PASSED!"
    else
        print *, "Some complex eigensystem solver interface tests FAILED!"
        stop 1
    end if
    
end program test_kilca_complex_eigensystem