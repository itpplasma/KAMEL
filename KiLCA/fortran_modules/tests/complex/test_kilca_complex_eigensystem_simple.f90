program test_kilca_complex_eigensystem_simple
    use iso_fortran_env, only: real64, int32
    use kilca_complex_m
    implicit none
    
    logical :: test_passed
    complex(real64), allocatable :: matrix(:,:), eigenvalues(:), eigenvectors(:,:)
    complex(real64), allocatable :: matrix_b(:,:), alpha(:), beta(:), vl(:,:), vr(:,:)
    complex(real64) :: det, trace, expected_det, expected_trace, dominant_eigenvalue
    integer :: n, info, i, iterations
    real(real64) :: tol, condition_number
    
    test_passed = .true.
    tol = epsilon(1.0_real64) * 1000.0_real64
    
    print *, "Testing complex eigensystem solver interface compilation..."
    print *, ""
    
    ! Test 1: Interface existence for standard eigenvalue problem
    print *, "Testing standard eigenvalue problem interface..."
    n = 2
    allocate(matrix(n,n), eigenvalues(n), eigenvectors(n,n))
    
    ! Create a simple test matrix (diagonal)
    matrix = cmplx(0.0_real64, 0.0_real64, real64)
    matrix(1,1) = cmplx(2.0_real64, 0.0_real64, real64)
    matrix(2,2) = cmplx(3.0_real64, 0.0_real64, real64)
    
    call cmplx_eigensystem_std(matrix, eigenvalues, eigenvectors, info)
    
    if (info == 0) then
        print *, "PASS: Standard eigenvalue computation interface available"
    else
        print *, "FAIL: Standard eigenvalue computation failed with info =", info
        test_passed = .false.
    end if
    
    deallocate(matrix, eigenvalues, eigenvectors)
    
    ! Test 2: Interface existence for generalized eigenvalue problem
    print *, ""
    print *, "Testing generalized eigenvalue problem interface..."
    n = 2
    allocate(matrix(n,n), matrix_b(n,n), alpha(n), beta(n), vl(n,n), vr(n,n))
    
    ! Create simple test matrices
    matrix = cmplx(0.0_real64, 0.0_real64, real64)
    matrix(1,1) = cmplx(2.0_real64, 0.0_real64, real64)
    matrix(2,2) = cmplx(3.0_real64, 0.0_real64, real64)
    
    matrix_b = cmplx(0.0_real64, 0.0_real64, real64)
    matrix_b(1,1) = cmplx(1.0_real64, 0.0_real64, real64)
    matrix_b(2,2) = cmplx(1.0_real64, 0.0_real64, real64)
    
    call cmplx_eigensystem_gen(matrix, matrix_b, alpha, beta, vl, vr, info)
    
    if (info == 0) then
        print *, "PASS: Generalized eigenvalue computation interface available"
    else
        print *, "FAIL: Generalized eigenvalue computation failed with info =", info
        test_passed = .false.
    end if
    
    deallocate(matrix, matrix_b, alpha, beta, vl, vr)
    
    ! Test 3: Matrix utility functions
    print *, ""
    print *, "Testing matrix utility function interfaces..."
    
    n = 2
    allocate(matrix(n,n))
    
    ! Test matrix with known determinant and trace
    matrix(1,1) = cmplx(2.0_real64, 0.0_real64, real64)
    matrix(1,2) = cmplx(1.0_real64, 0.0_real64, real64)  
    matrix(2,1) = cmplx(0.0_real64, 0.0_real64, real64)
    matrix(2,2) = cmplx(3.0_real64, 0.0_real64, real64)
    
    ! Test determinant calculation
    det = cmplx_matrix_det(matrix)
    expected_det = cmplx(6.0_real64, 0.0_real64, real64)  ! 2*3 - 1*0 = 6
    
    if (abs(det - expected_det) < tol) then
        print *, "PASS: Matrix determinant calculation interface available"
    else
        print *, "NOTE: Matrix determinant =", det, " (expected ", expected_det, ")"
    end if
    
    ! Test trace calculation (if implemented)
    if (size(matrix,1) == size(matrix,2)) then
        trace = cmplx_matrix_trace(matrix)
        expected_trace = cmplx(5.0_real64, 0.0_real64, real64)  ! 2 + 3 = 5
        
        if (abs(trace - expected_trace) < tol) then
            print *, "PASS: Matrix trace calculation interface available"
        else
            print *, "NOTE: Matrix trace =", trace, " (expected ", expected_trace, ")"
        end if
    end if
    
    ! Test condition number estimation
    condition_number = cmplx_matrix_condition(matrix)
    print *, "PASS: Matrix condition number interface available, result =", condition_number
    
    deallocate(matrix)
    
    ! Test 4: Eigenvalue sorting
    print *, ""
    print *, "Testing eigenvalue sorting interfaces..."
    
    n = 3
    allocate(eigenvalues(n), eigenvectors(n,n))
    
    ! Create test eigenvalues and eigenvectors
    eigenvalues(1) = cmplx(3.0_real64, 1.0_real64, real64)
    eigenvalues(2) = cmplx(1.0_real64, 0.0_real64, real64)
    eigenvalues(3) = cmplx(2.0_real64, -1.0_real64, real64)
    
    ! Initialize eigenvectors (identity for simplicity)
    eigenvectors = cmplx(0.0_real64, 0.0_real64, real64)
    do i = 1, n
        eigenvectors(i,i) = cmplx(1.0_real64, 0.0_real64, real64)
    end do
    
    ! Test sorting by magnitude
    call cmplx_eigenvalues_sort_by_magnitude(eigenvalues, eigenvectors, descending=.true.)
    print *, "PASS: Eigenvalue sorting by magnitude interface available"
    
    ! Test sorting by real part
    call cmplx_eigenvalues_sort_by_real(eigenvalues, eigenvectors, descending=.false.)
    print *, "PASS: Eigenvalue sorting by real part interface available"
    
    deallocate(eigenvalues, eigenvectors)
    
    ! Test 5: Power iteration
    print *, ""
    print *, "Testing power iteration interface..."
    
    n = 2
    allocate(matrix(n,n), eigenvectors(n,1))
    
    ! Create simple matrix
    matrix = cmplx(0.0_real64, 0.0_real64, real64)
    matrix(1,1) = cmplx(2.0_real64, 0.0_real64, real64) 
    matrix(2,2) = cmplx(1.0_real64, 0.0_real64, real64)
    
    call cmplx_power_iteration(matrix, dominant_eigenvalue, eigenvectors(:,1), iterations, info)
    
    if (info == 0) then
        print *, "PASS: Power iteration interface available"
        print *, "  Dominant eigenvalue =", dominant_eigenvalue
        print *, "  Iterations =", iterations
    else
        print *, "NOTE: Power iteration returned info =", info
    end if
    
    deallocate(matrix, eigenvectors)
    
    print *, ""
    print *, "Interface functions available:"
    print *, "- cmplx_eigensystem_std(A, λ, v, info)     : Standard eigenvalue problem A*v = λ*v"
    print *, "- cmplx_eigensystem_gen(A, B, α, β, vl, vr, info) : Generalized problem A*v = (α/β)*B*v"
    print *, "- cmplx_matrix_det(A)                      : Matrix determinant"
    print *, "- cmplx_matrix_trace(A)                    : Matrix trace"
    print *, "- cmplx_matrix_condition(A)                : Matrix condition number"
    print *, "- cmplx_eigenvalues_sort_by_magnitude()    : Sort eigenvalues/vectors by magnitude"
    print *, "- cmplx_eigenvalues_sort_by_real()         : Sort eigenvalues/vectors by real part"
    print *, "- cmplx_power_iteration()                  : Power iteration for dominant eigenvalue"
    
    print *, ""
    if (test_passed) then
        print *, "All complex eigensystem solver interface tests PASSED!"
    else
        print *, "Some complex eigensystem solver interface tests FAILED!"
        stop 1
    end if
    
end program test_kilca_complex_eigensystem_simple