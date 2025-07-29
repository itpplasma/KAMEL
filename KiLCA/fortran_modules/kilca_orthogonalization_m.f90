module kilca_orthogonalization_m
    use iso_fortran_env, only: real64, int32, error_unit
    use kilca_types_m
    implicit none
    
    private
    
    ! Public procedures
    public :: orthogonalize_basis_vectors_complex
    public :: orthogonalize_basis_vectors_real
    public :: qr_decomposition_complex
    public :: qr_decomposition_real
    public :: apply_orthogonal_transformation_complex
    public :: apply_orthogonal_transformation_real
    public :: verify_orthogonality_complex
    public :: verify_orthogonality_real
    
contains

    !---------------------------------------------------------------------------
    ! High-level complex basis orthogonalization (matching C++ solver functionality)
    !---------------------------------------------------------------------------
    subroutine orthogonalize_basis_vectors_complex(m, n, basis_vectors, ierr)
        integer, intent(in) :: m, n
        complex(real64), intent(inout) :: basis_vectors(m, n)
        integer, intent(out) :: ierr
        
        complex(real64), allocatable :: tau(:), work(:)
        integer :: k, lwork
        
        ierr = 0
        k = min(m, n)
        
        if (k <= 0) then
            ierr = -1
            return
        end if
        
        allocate(tau(k))
        
        ! Workspace query
        lwork = -1
        allocate(work(1))
        call zgeqrf(m, n, basis_vectors, m, tau, work, lwork, ierr)
        
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'orthogonalize_basis_vectors_complex: workspace query failed, info=', ierr
            deallocate(tau, work)
            return
        end if
        
        lwork = int(real(work(1)))
        deallocate(work)
        allocate(work(lwork))
        
        ! QR factorization
        call zgeqrf(m, n, basis_vectors, m, tau, work, lwork, ierr)
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'orthogonalize_basis_vectors_complex: zgeqrf failed, info=', ierr
            deallocate(tau, work)
            return
        end if
        
        ! Generate orthogonal Q matrix
        call zungqr(m, n, k, basis_vectors, m, tau, work, lwork, ierr)
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'orthogonalize_basis_vectors_complex: zungqr failed, info=', ierr
        end if
        
        deallocate(tau, work)
        
    end subroutine orthogonalize_basis_vectors_complex
    
    !---------------------------------------------------------------------------
    ! High-level real basis orthogonalization
    !---------------------------------------------------------------------------
    subroutine orthogonalize_basis_vectors_real(m, n, basis_vectors, ierr)
        integer, intent(in) :: m, n
        real(real64), intent(inout) :: basis_vectors(m, n)
        integer, intent(out) :: ierr
        
        real(real64), allocatable :: tau(:), work(:)
        integer :: k, lwork
        
        ierr = 0
        k = min(m, n)
        
        if (k <= 0) then
            ierr = -1
            return
        end if
        
        allocate(tau(k))
        
        ! Workspace query
        lwork = -1
        allocate(work(1))
        call dgeqrf(m, n, basis_vectors, m, tau, work, lwork, ierr)
        
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'orthogonalize_basis_vectors_real: workspace query failed, info=', ierr
            deallocate(tau, work)
            return
        end if
        
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        
        ! QR factorization
        call dgeqrf(m, n, basis_vectors, m, tau, work, lwork, ierr)
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'orthogonalize_basis_vectors_real: dgeqrf failed, info=', ierr
            deallocate(tau, work)
            return
        end if
        
        ! Generate orthogonal Q matrix
        call dorgqr(m, n, k, basis_vectors, m, tau, work, lwork, ierr)
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'orthogonalize_basis_vectors_real: dorgqr failed, info=', ierr
        end if
        
        deallocate(tau, work)
        
    end subroutine orthogonalize_basis_vectors_real
    
    !---------------------------------------------------------------------------
    ! Complex QR decomposition with separate Q and R matrices
    !---------------------------------------------------------------------------
    subroutine qr_decomposition_complex(m, n, A, Q, R, ierr)
        integer, intent(in) :: m, n
        complex(real64), intent(in) :: A(m, n)
        complex(real64), intent(out) :: Q(m, n), R(n, n)
        integer, intent(out) :: ierr
        
        complex(real64), allocatable :: tau(:), work(:)
        integer :: k, lwork, i, j
        
        ierr = 0
        k = min(m, n)
        
        allocate(tau(k))
        
        ! Copy A to Q
        Q = A
        
        ! Workspace query
        lwork = -1
        allocate(work(1))
        call zgeqrf(m, n, Q, m, tau, work, lwork, ierr)
        
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'qr_decomposition_complex: workspace query failed, info=', ierr
            deallocate(tau, work)
            return
        end if
        
        lwork = int(real(work(1)))
        deallocate(work)
        allocate(work(lwork))
        
        ! QR factorization
        call zgeqrf(m, n, Q, m, tau, work, lwork, ierr)
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'qr_decomposition_complex: zgeqrf failed, info=', ierr
            deallocate(tau, work)
            return
        end if
        
        ! Extract R matrix (upper triangular part)
        R = cmplx(0.0_real64, 0.0_real64, real64)
        do j = 1, n
            do i = 1, min(j, m)
                R(i, j) = Q(i, j)
            end do
        end do
        
        ! Generate orthogonal Q matrix
        call zungqr(m, n, k, Q, m, tau, work, lwork, ierr)
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'qr_decomposition_complex: zungqr failed, info=', ierr
        end if
        
        deallocate(tau, work)
        
    end subroutine qr_decomposition_complex
    
    !---------------------------------------------------------------------------
    ! Real QR decomposition with separate Q and R matrices
    !---------------------------------------------------------------------------
    subroutine qr_decomposition_real(m, n, A, Q, R, ierr)
        integer, intent(in) :: m, n
        real(real64), intent(in) :: A(m, n)
        real(real64), intent(out) :: Q(m, n), R(n, n)
        integer, intent(out) :: ierr
        
        real(real64), allocatable :: tau(:), work(:)
        integer :: k, lwork, i, j
        
        ierr = 0
        k = min(m, n)
        
        allocate(tau(k))
        
        ! Copy A to Q
        Q = A
        
        ! Workspace query
        lwork = -1
        allocate(work(1))
        call dgeqrf(m, n, Q, m, tau, work, lwork, ierr)
        
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'qr_decomposition_real: workspace query failed, info=', ierr
            deallocate(tau, work)
            return
        end if
        
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        
        ! QR factorization
        call dgeqrf(m, n, Q, m, tau, work, lwork, ierr)
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'qr_decomposition_real: dgeqrf failed, info=', ierr
            deallocate(tau, work)
            return
        end if
        
        ! Extract R matrix
        R = 0.0_real64
        do j = 1, n
            do i = 1, min(j, m)
                R(i, j) = Q(i, j)
            end do
        end do
        
        ! Generate orthogonal Q matrix
        call dorgqr(m, n, k, Q, m, tau, work, lwork, ierr)
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'qr_decomposition_real: dorgqr failed, info=', ierr
        end if
        
        deallocate(tau, work)
        
    end subroutine qr_decomposition_real
    
    !---------------------------------------------------------------------------
    ! Apply complex orthogonal transformation Q or Q^H to matrix C
    !---------------------------------------------------------------------------
    subroutine apply_orthogonal_transformation_complex(side, trans, m, n, k, Q, tau, C, ierr)
        character, intent(in) :: side, trans
        integer, intent(in) :: m, n, k
        complex(real64), intent(in) :: Q(m, k), tau(k)
        complex(real64), intent(inout) :: C(m, n)
        integer, intent(out) :: ierr
        
        complex(real64), allocatable :: work(:)
        integer :: lwork, ldc
        
        ierr = 0
        ldc = m
        
        ! Workspace query
        lwork = -1
        allocate(work(1))
        call zunmqr(side, trans, m, n, k, Q, m, tau, C, ldc, work, lwork, ierr)
        
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'apply_orthogonal_transformation_complex: workspace query failed, info=', ierr
            deallocate(work)
            return
        end if
        
        lwork = int(real(work(1)))
        deallocate(work)
        allocate(work(lwork))
        
        ! Apply transformation
        call zunmqr(side, trans, m, n, k, Q, m, tau, C, ldc, work, lwork, ierr)
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'apply_orthogonal_transformation_complex: zunmqr failed, info=', ierr
        end if
        
        deallocate(work)
        
    end subroutine apply_orthogonal_transformation_complex
    
    !---------------------------------------------------------------------------
    ! Apply real orthogonal transformation Q or Q^T to matrix C
    !---------------------------------------------------------------------------
    subroutine apply_orthogonal_transformation_real(side, trans, m, n, k, Q, tau, C, ierr)
        character, intent(in) :: side, trans
        integer, intent(in) :: m, n, k
        real(real64), intent(in) :: Q(m, k), tau(k)
        real(real64), intent(inout) :: C(m, n)
        integer, intent(out) :: ierr
        
        real(real64), allocatable :: work(:)
        integer :: lwork, ldc
        
        ierr = 0
        ldc = m
        
        ! Workspace query
        lwork = -1
        allocate(work(1))
        call dormqr(side, trans, m, n, k, Q, m, tau, C, ldc, work, lwork, ierr)
        
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'apply_orthogonal_transformation_real: workspace query failed, info=', ierr
            deallocate(work)
            return
        end if
        
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        
        ! Apply transformation
        call dormqr(side, trans, m, n, k, Q, m, tau, C, ldc, work, lwork, ierr)
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'apply_orthogonal_transformation_real: dormqr failed, info=', ierr
        end if
        
        deallocate(work)
        
    end subroutine apply_orthogonal_transformation_real
    
    !---------------------------------------------------------------------------
    ! Verify orthogonality of complex Q matrix
    !---------------------------------------------------------------------------
    subroutine verify_orthogonality_complex(m, n, Q, is_orthogonal, max_error)
        integer, intent(in) :: m, n
        complex(real64), intent(in) :: Q(m, n)
        logical, intent(out) :: is_orthogonal
        real(real64), intent(out) :: max_error
        
        complex(real64), allocatable :: QtQ(:,:)
        real(real64), parameter :: tol = 1.0e-12_real64
        integer :: i, j
        
        allocate(QtQ(n, n))
        
        ! Compute Q^H * Q
        call zgemm('C', 'N', n, n, m, cmplx(1.0_real64, 0.0_real64, real64), &
                   Q, m, Q, m, cmplx(0.0_real64, 0.0_real64, real64), QtQ, n)
        
        ! Check deviation from identity
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
        
        is_orthogonal = (max_error < tol)
        deallocate(QtQ)
        
    end subroutine verify_orthogonality_complex
    
    !---------------------------------------------------------------------------
    ! Verify orthogonality of real Q matrix
    !---------------------------------------------------------------------------
    subroutine verify_orthogonality_real(m, n, Q, is_orthogonal, max_error)
        integer, intent(in) :: m, n
        real(real64), intent(in) :: Q(m, n)
        logical, intent(out) :: is_orthogonal
        real(real64), intent(out) :: max_error
        
        real(real64), allocatable :: QtQ(:,:)
        real(real64), parameter :: tol = 1.0e-12_real64
        integer :: i, j
        
        allocate(QtQ(n, n))
        
        ! Compute Q^T * Q
        call dgemm('T', 'N', n, n, m, 1.0_real64, Q, m, Q, m, 0.0_real64, QtQ, n)
        
        ! Check deviation from identity
        max_error = 0.0_real64
        do i = 1, n
            do j = 1, n
                if (i == j) then
                    max_error = max(max_error, abs(QtQ(i,j) - 1.0_real64))
                else
                    max_error = max(max_error, abs(QtQ(i,j)))
                end if
            end do
        end do
        
        is_orthogonal = (max_error < tol)
        deallocate(QtQ)
        
    end subroutine verify_orthogonality_real

end module kilca_orthogonalization_m