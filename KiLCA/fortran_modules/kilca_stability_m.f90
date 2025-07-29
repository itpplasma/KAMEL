module kilca_stability_m
    use iso_fortran_env, only: real64, int32
    use kilca_types_m
    implicit none
    
    private
    
    ! LAPACK interface declarations
    interface
        real(real64) function zlange(norm, m, n, a, lda, work)
            import :: real64
            character :: norm
            integer :: m, n, lda
            complex(real64) :: a(lda,*)
            real(real64) :: work(*)
        end function zlange
        
        real(real64) function dlange(norm, m, n, a, lda, work)
            import :: real64
            character :: norm
            integer :: m, n, lda
            real(real64) :: a(lda,*)
            real(real64) :: work(*)
        end function dlange
    end interface
    
    ! Public procedures
    public :: estimate_condition_number
    public :: estimate_rcond
    public :: check_matrix_stability
    public :: is_matrix_singular
    public :: compute_numerical_rank
    public :: compute_matrix_norm
    public :: compute_det_magnitude
    public :: check_system_stability
    public :: analyze_perturbation_sensitivity
    
    ! Interface for norm computation
    interface compute_matrix_norm
        module procedure compute_matrix_norm_complex
        module procedure compute_matrix_norm_real
    end interface compute_matrix_norm
    
contains

    !---------------------------------------------------------------------------
    ! Estimate condition number using SVD
    !---------------------------------------------------------------------------
    function estimate_condition_number(A, info) result(cond_num)
        complex(real64), intent(in) :: A(:,:)
        integer, intent(out) :: info
        real(real64) :: cond_num
        
        complex(real64), allocatable :: A_copy(:,:), U(:,:), VT(:,:), work(:)
        real(real64), allocatable :: S(:), rwork(:)
        integer :: m, n, lwork, lrwork
        
        info = 0
        m = size(A, 1)
        n = size(A, 2)
        
        ! Check for square matrix
        if (m /= n) then
            info = -1
            cond_num = -1.0_real64
            return
        end if
        
        ! Allocate workspace
        allocate(A_copy(m, n))
        allocate(S(min(m, n)))
        allocate(U(m, m))
        allocate(VT(n, n))
        
        ! Copy matrix (SVD destroys input)
        A_copy = A
        
        ! Query optimal workspace
        lwork = -1
        allocate(work(1))
        allocate(rwork(5*min(m,n)))
        call zgesvd('N', 'N', m, n, A_copy, m, S, U, m, VT, n, work, lwork, rwork, info)
        
        if (info /= 0) then
            cond_num = -1.0_real64
            return
        end if
        
        lwork = int(real(work(1)))
        deallocate(work)
        allocate(work(lwork))
        
        ! Compute SVD
        call zgesvd('N', 'N', m, n, A_copy, m, S, U, m, VT, n, work, lwork, rwork, info)
        
        if (info /= 0) then
            cond_num = -1.0_real64
            return
        end if
        
        ! Compute condition number
        if (S(n) > epsilon(1.0_real64)) then
            cond_num = S(1) / S(n)
        else
            cond_num = huge(1.0_real64)
        end if
        
        deallocate(A_copy, S, U, VT, work, rwork)
        
    end function estimate_condition_number
    
    !---------------------------------------------------------------------------
    ! Estimate reciprocal condition number using LAPACK's zgecon
    !---------------------------------------------------------------------------
    function estimate_rcond(A, info) result(rcond)
        complex(real64), intent(in) :: A(:,:)
        integer, intent(out) :: info
        real(real64) :: rcond
        
        complex(real64), allocatable :: A_copy(:,:), work(:)
        real(real64), allocatable :: rwork(:)
        integer, allocatable :: ipiv(:)
        integer :: n
        real(real64) :: anorm
        
        info = 0
        n = size(A, 1)
        
        ! Check for square matrix
        if (size(A, 2) /= n) then
            info = -1
            rcond = -1.0_real64
            return
        end if
        
        ! Allocate workspace
        allocate(A_copy(n, n))
        allocate(ipiv(n))
        allocate(work(2*n))
        allocate(rwork(2*n))
        
        ! Copy matrix and compute norm
        A_copy = A
        anorm = compute_matrix_norm_complex(A, '1')
        
        ! LU factorization
        call zgetrf(n, n, A_copy, n, ipiv, info)
        if (info /= 0) then
            rcond = 0.0_real64
            return
        end if
        
        ! Estimate reciprocal condition number
        call zgecon('1', n, A_copy, n, anorm, rcond, work, rwork, info)
        
        deallocate(A_copy, ipiv, work, rwork)
        
    end function estimate_rcond
    
    !---------------------------------------------------------------------------
    ! Check if matrix is numerically stable for computations
    !---------------------------------------------------------------------------
    function check_matrix_stability(A, tolerance) result(is_stable)
        complex(real64), intent(in) :: A(:,:)
        real(real64), intent(in) :: tolerance
        logical :: is_stable
        
        real(real64) :: rcond
        integer :: info
        
        rcond = estimate_rcond(A, info)
        
        if (info /= 0) then
            is_stable = .false.
        else
            is_stable = rcond > tolerance
        end if
        
    end function check_matrix_stability
    
    !---------------------------------------------------------------------------
    ! Check if matrix is singular
    !---------------------------------------------------------------------------
    function is_matrix_singular(A, tolerance) result(is_singular)
        complex(real64), intent(in) :: A(:,:)
        real(real64), intent(in) :: tolerance
        logical :: is_singular
        
        real(real64) :: rcond
        integer :: info
        
        rcond = estimate_rcond(A, info)
        
        if (info /= 0) then
            is_singular = .true.
        else
            is_singular = rcond < tolerance
        end if
        
    end function is_matrix_singular
    
    !---------------------------------------------------------------------------
    ! Compute numerical rank of matrix
    !---------------------------------------------------------------------------
    function compute_numerical_rank(A, tolerance, info) result(rank)
        complex(real64), intent(in) :: A(:,:)
        real(real64), intent(in) :: tolerance
        integer, intent(out) :: info
        integer :: rank
        
        complex(real64), allocatable :: A_copy(:,:), U(:,:), VT(:,:), work(:)
        real(real64), allocatable :: S(:), rwork(:)
        integer :: m, n, lwork, i
        
        info = 0
        m = size(A, 1)
        n = size(A, 2)
        
        ! Allocate workspace
        allocate(A_copy(m, n))
        allocate(S(min(m, n)))
        allocate(U(m, m))
        allocate(VT(n, n))
        
        ! Copy matrix
        A_copy = A
        
        ! Query optimal workspace
        lwork = -1
        allocate(work(1))
        allocate(rwork(5*min(m,n)))
        call zgesvd('N', 'N', m, n, A_copy, m, S, U, m, VT, n, work, lwork, rwork, info)
        
        if (info /= 0) then
            rank = -1
            return
        end if
        
        lwork = int(real(work(1)))
        deallocate(work)
        allocate(work(lwork))
        
        ! Compute SVD
        call zgesvd('N', 'N', m, n, A_copy, m, S, U, m, VT, n, work, lwork, rwork, info)
        
        if (info /= 0) then
            rank = -1
            return
        end if
        
        ! Count singular values above tolerance
        ! Note: SVD returns singular values in descending order
        rank = 0
        do i = 1, min(m, n)
            if (S(i) > tolerance) then
                rank = rank + 1
            end if
        end do
        
        deallocate(A_copy, S, U, VT, work, rwork)
        
    end function compute_numerical_rank
    
    !---------------------------------------------------------------------------
    ! Compute matrix norm (complex version)
    !---------------------------------------------------------------------------
    function compute_matrix_norm_complex(A, norm_type) result(norm_val)
        complex(real64), intent(in) :: A(:,:)
        character(len=1), intent(in) :: norm_type
        real(real64) :: norm_val
        
        real(real64), allocatable :: work(:)
        integer :: m, n
        
        m = size(A, 1)
        n = size(A, 2)
        
        if (norm_type == 'F' .or. norm_type == 'f') then
            ! Frobenius norm - simple implementation
            norm_val = sqrt(sum(abs(A)**2))
        else
            ! Use LAPACK for other norms
            allocate(work(max(m, n)))
            norm_val = zlange(norm_type, m, n, A, m, work)
            deallocate(work)
        end if
        
    end function compute_matrix_norm_complex
    
    !---------------------------------------------------------------------------
    ! Compute matrix norm (real version)
    !---------------------------------------------------------------------------
    function compute_matrix_norm_real(A, norm_type) result(norm_val)
        real(real64), intent(in) :: A(:,:)
        character(len=1), intent(in) :: norm_type
        real(real64) :: norm_val
        
        real(real64), allocatable :: work(:)
        integer :: m, n
        
        m = size(A, 1)
        n = size(A, 2)
        
        if (norm_type == 'F' .or. norm_type == 'f') then
            ! Frobenius norm
            norm_val = sqrt(sum(A**2))
        else
            ! Use LAPACK for other norms
            allocate(work(max(m, n)))
            norm_val = dlange(norm_type, m, n, A, m, work)
            deallocate(work)
        end if
        
    end function compute_matrix_norm_real
    
    !---------------------------------------------------------------------------
    ! Compute determinant magnitude for stability analysis
    !---------------------------------------------------------------------------
    function compute_det_magnitude(A, info) result(det_mag)
        complex(real64), intent(in) :: A(:,:)
        integer, intent(out) :: info
        real(real64) :: det_mag
        
        complex(real64), allocatable :: A_copy(:,:)
        integer, allocatable :: ipiv(:)
        complex(real64) :: det
        integer :: i, n
        
        info = 0
        n = size(A, 1)
        
        ! Check for square matrix
        if (size(A, 2) /= n) then
            info = -1
            det_mag = -1.0_real64
            return
        end if
        
        ! Allocate workspace
        allocate(A_copy(n, n))
        allocate(ipiv(n))
        
        ! Copy matrix
        A_copy = A
        
        ! LU factorization
        call zgetrf(n, n, A_copy, n, ipiv, info)
        
        if (info /= 0) then
            det_mag = 0.0_real64
            return
        end if
        
        ! Compute determinant from diagonal of U
        det = cmplx(1.0_real64, 0.0_real64, real64)
        do i = 1, n
            det = det * A_copy(i, i)
        end do
        
        ! Account for row swaps
        do i = 1, n
            if (ipiv(i) /= i) det = -det
        end do
        
        det_mag = abs(det)
        
        deallocate(A_copy, ipiv)
        
    end function compute_det_magnitude
    
    !---------------------------------------------------------------------------
    ! Check stability of linear system A*x = b
    !---------------------------------------------------------------------------
    function check_system_stability(A, b, tolerance) result(is_stable)
        complex(real64), intent(in) :: A(:,:), b(:)
        real(real64), intent(in) :: tolerance
        logical :: is_stable
        
        real(real64) :: cond_num, norm_A, norm_b
        integer :: info
        
        ! Check dimensions
        if (size(A, 1) /= size(b) .or. size(A, 1) /= size(A, 2)) then
            is_stable = .false.
            return
        end if
        
        ! Check matrix condition
        cond_num = estimate_condition_number(A, info)
        
        if (info /= 0 .or. cond_num < 0.0_real64) then
            is_stable = .false.
            return
        end if
        
        ! Check relative conditioning
        norm_A = compute_matrix_norm_complex(A, 'F')
        norm_b = sqrt(sum(abs(b)**2))
        
        ! System is stable if condition number is reasonable and norms are finite
        is_stable = (cond_num < 1.0_real64/tolerance) .and. &
                   (norm_A > epsilon(1.0_real64)) .and. &
                   (norm_b > epsilon(1.0_real64)) .and. &
                   (norm_A < huge(1.0_real64)/1000.0_real64) .and. &
                   (norm_b < huge(1.0_real64)/1000.0_real64)
        
    end function check_system_stability
    
    !---------------------------------------------------------------------------
    ! Analyze sensitivity to perturbations
    !---------------------------------------------------------------------------
    function analyze_perturbation_sensitivity(A, A_perturbed, info) result(sensitivity)
        complex(real64), intent(in) :: A(:,:), A_perturbed(:,:)
        integer, intent(out) :: info
        real(real64) :: sensitivity
        
        real(real64) :: norm_A, norm_delta, cond_A
        integer :: m, n
        
        info = 0
        m = size(A, 1)
        n = size(A, 2)
        
        ! Check dimensions
        if (m /= n .or. size(A_perturbed, 1) /= m .or. size(A_perturbed, 2) /= n) then
            info = -1
            sensitivity = -1.0_real64
            return
        end if
        
        ! Compute norms
        norm_A = compute_matrix_norm_complex(A, 'F')
        norm_delta = compute_matrix_norm_complex(A_perturbed - A, 'F')
        
        ! Compute condition number
        cond_A = estimate_condition_number(A, info)
        
        if (info /= 0) then
            sensitivity = -1.0_real64
            return
        end if
        
        ! Relative perturbation sensitivity
        if (norm_A > epsilon(1.0_real64)) then
            sensitivity = cond_A * (norm_delta / norm_A)
        else
            sensitivity = huge(1.0_real64)
        end if
        
    end function analyze_perturbation_sensitivity

end module kilca_stability_m