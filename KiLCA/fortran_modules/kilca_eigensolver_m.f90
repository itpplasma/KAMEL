module kilca_eigensolver_m
    use iso_fortran_env, only: real64, int32, error_unit
    use kilca_types_m
    implicit none
    
    private
    
    ! Public procedures
    public :: eigsys_zgeev
    public :: solve_hermitian_eigenproblem
    public :: solve_general_real_eigenproblem
    public :: solve_symmetric_real_eigenproblem
    public :: solve_mode_eigenvalue_problem
    public :: sort_eigenvalues_by_growth_rate
    
contains

    !---------------------------------------------------------------------------
    ! Direct replacement for eigsys.f90 using ZGEEV
    !---------------------------------------------------------------------------
    subroutine eigsys_zgeev(dim, mat, evl, evc, ierr)
        integer, intent(in) :: dim
        complex(real64), intent(in) :: mat(dim,dim)
        complex(real64), intent(out) :: evl(dim)
        complex(real64), intent(out) :: evc(dim,dim)
        integer, intent(out) :: ierr
        
        ! Local variables
        complex(real64), allocatable :: A(:,:), work(:), VL(:,:)
        real(real64), allocatable :: rwork(:)
        integer :: lwork
        character :: jobvl, jobvr
        
        ierr = 0
        jobvl = 'N'  ! Don't compute left eigenvectors
        jobvr = 'V'  ! Compute right eigenvectors
        
        ! Allocate arrays
        allocate(A(dim,dim), VL(dim,dim))
        allocate(rwork(2*dim))
        
        ! Make copy of input matrix (ZGEEV destroys it)
        A = mat
        
        ! Query optimal workspace
        lwork = -1
        allocate(work(1))
        call zgeev(jobvl, jobvr, dim, A, dim, evl, VL, dim, evc, dim, &
                   work, lwork, rwork, ierr)
        
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'eigsys_zgeev: workspace query failed, info=', ierr
            deallocate(A, VL, work, rwork)
            return
        end if
        
        ! Allocate optimal workspace
        lwork = int(real(work(1)))
        deallocate(work)
        allocate(work(lwork))
        
        ! Restore matrix (was modified during workspace query)
        A = mat
        
        ! Compute eigenvalues and right eigenvectors
        call zgeev(jobvl, jobvr, dim, A, dim, evl, VL, dim, evc, dim, &
                   work, lwork, rwork, ierr)
        
        if (ierr /= 0) then
            if (ierr < 0) then
                write(error_unit,'(A,I0)') 'eigsys_zgeev: illegal argument ', -ierr
            else
                write(error_unit,'(A,I0)') 'eigsys_zgeev: failed to converge, ', ierr, &
                    ' eigenvalues did not converge'
            end if
        end if
        
        deallocate(A, VL, work, rwork)
        
    end subroutine eigsys_zgeev
    
    !---------------------------------------------------------------------------
    ! Solve Hermitian eigenvalue problem using ZHEEV
    !---------------------------------------------------------------------------
    subroutine solve_hermitian_eigenproblem(n, H_matrix, eigenvals, eigenvecs, ierr)
        integer, intent(in) :: n
        complex(real64), intent(in) :: H_matrix(n,n)
        real(real64), intent(out) :: eigenvals(n)
        complex(real64), intent(out) :: eigenvecs(n,n)
        integer, intent(out) :: ierr
        
        ! Local variables
        complex(real64), allocatable :: A(:,:), work(:)
        real(real64), allocatable :: rwork(:)
        integer :: lwork
        character :: jobz, uplo
        
        ierr = 0
        jobz = 'V'  ! Compute eigenvalues and eigenvectors
        uplo = 'U'  ! Upper triangle of A is stored
        
        ! Allocate arrays
        allocate(A(n,n))
        allocate(rwork(3*n-2))
        
        ! Make copy of input matrix
        A = H_matrix
        
        ! Query optimal workspace
        lwork = -1
        allocate(work(1))
        call zheev(jobz, uplo, n, A, n, eigenvals, work, lwork, rwork, ierr)
        
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'solve_hermitian_eigenproblem: workspace query failed, info=', ierr
            deallocate(A, work, rwork)
            return
        end if
        
        ! Allocate optimal workspace
        lwork = int(real(work(1)))
        deallocate(work)
        allocate(work(lwork))
        
        ! Restore matrix
        A = H_matrix
        
        ! Compute eigenvalues and eigenvectors
        call zheev(jobz, uplo, n, A, n, eigenvals, work, lwork, rwork, ierr)
        
        if (ierr == 0) then
            ! Copy eigenvectors (stored in A after ZHEEV)
            eigenvecs = A
        else
            if (ierr < 0) then
                write(error_unit,'(A,I0)') 'solve_hermitian_eigenproblem: illegal argument ', -ierr
            else
                write(error_unit,'(A,I0)') 'solve_hermitian_eigenproblem: failed to converge'
            end if
        end if
        
        deallocate(A, work, rwork)
        
    end subroutine solve_hermitian_eigenproblem
    
    !---------------------------------------------------------------------------
    ! Solve general real eigenvalue problem using DGEEV
    !---------------------------------------------------------------------------
    subroutine solve_general_real_eigenproblem(n, A_matrix, eigenvals_real, eigenvals_imag, &
                                              eigenvecs, ierr)
        integer, intent(in) :: n
        real(real64), intent(in) :: A_matrix(n,n)
        real(real64), intent(out) :: eigenvals_real(n), eigenvals_imag(n)
        real(real64), intent(out) :: eigenvecs(n,n)
        integer, intent(out) :: ierr
        
        ! Local variables
        real(real64), allocatable :: A(:,:), work(:), VL(:,:)
        integer :: lwork
        character :: jobvl, jobvr
        
        ierr = 0
        jobvl = 'N'  ! Don't compute left eigenvectors
        jobvr = 'V'  ! Compute right eigenvectors
        
        ! Allocate arrays
        allocate(A(n,n), VL(n,n))
        
        ! Make copy of input matrix
        A = A_matrix
        
        ! Query optimal workspace
        lwork = -1
        allocate(work(1))
        call dgeev(jobvl, jobvr, n, A, n, eigenvals_real, eigenvals_imag, &
                   VL, n, eigenvecs, n, work, lwork, ierr)
        
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'solve_general_real_eigenproblem: workspace query failed, info=', ierr
            deallocate(A, VL, work)
            return
        end if
        
        ! Allocate optimal workspace
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        
        ! Restore matrix
        A = A_matrix
        
        ! Compute eigenvalues and eigenvectors
        call dgeev(jobvl, jobvr, n, A, n, eigenvals_real, eigenvals_imag, &
                   VL, n, eigenvecs, n, work, lwork, ierr)
        
        if (ierr /= 0) then
            if (ierr < 0) then
                write(error_unit,'(A,I0)') 'solve_general_real_eigenproblem: illegal argument ', -ierr
            else
                write(error_unit,'(A,I0)') 'solve_general_real_eigenproblem: failed to converge'
            end if
        end if
        
        deallocate(A, VL, work)
        
    end subroutine solve_general_real_eigenproblem
    
    !---------------------------------------------------------------------------
    ! Solve symmetric real eigenvalue problem using DSYEV
    !---------------------------------------------------------------------------
    subroutine solve_symmetric_real_eigenproblem(n, S_matrix, eigenvals, eigenvecs, ierr)
        integer, intent(in) :: n
        real(real64), intent(in) :: S_matrix(n,n)
        real(real64), intent(out) :: eigenvals(n)
        real(real64), intent(out) :: eigenvecs(n,n)
        integer, intent(out) :: ierr
        
        ! Local variables
        real(real64), allocatable :: A(:,:), work(:)
        integer :: lwork
        character :: jobz, uplo
        
        ierr = 0
        jobz = 'V'  ! Compute eigenvalues and eigenvectors
        uplo = 'U'  ! Upper triangle of A is stored
        
        ! Allocate arrays
        allocate(A(n,n))
        
        ! Make copy of input matrix
        A = S_matrix
        
        ! Query optimal workspace
        lwork = -1
        allocate(work(1))
        call dsyev(jobz, uplo, n, A, n, eigenvals, work, lwork, ierr)
        
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'solve_symmetric_real_eigenproblem: workspace query failed, info=', ierr
            deallocate(A, work)
            return
        end if
        
        ! Allocate optimal workspace
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        
        ! Restore matrix
        A = S_matrix
        
        ! Compute eigenvalues and eigenvectors
        call dsyev(jobz, uplo, n, A, n, eigenvals, work, lwork, ierr)
        
        if (ierr == 0) then
            ! Copy eigenvectors (stored in A after DSYEV)
            eigenvecs = A
        else
            if (ierr < 0) then
                write(error_unit,'(A,I0)') 'solve_symmetric_real_eigenproblem: illegal argument ', -ierr
            else
                write(error_unit,'(A,I0)') 'solve_symmetric_real_eigenproblem: failed to converge'
            end if
        end if
        
        deallocate(A, work)
        
    end subroutine solve_symmetric_real_eigenproblem
    
    !---------------------------------------------------------------------------
    ! Solve mode eigenvalue problem for plasma physics applications
    !---------------------------------------------------------------------------
    subroutine solve_mode_eigenvalue_problem(mode_matrix, frequencies, eigenvectors, &
                                            growth_rates, ierr)
        complex(real64), intent(in) :: mode_matrix(:,:)
        complex(real64), allocatable, intent(out) :: frequencies(:)
        complex(real64), allocatable, intent(out) :: eigenvectors(:,:)
        real(real64), allocatable, intent(out) :: growth_rates(:)
        integer, intent(out) :: ierr
        
        integer :: n
        
        n = size(mode_matrix, 1)
        
        ! Allocate output arrays
        allocate(frequencies(n), eigenvectors(n,n), growth_rates(n))
        
        ! Use ZGEEV for general complex eigenvalue problem
        call eigsys_zgeev(n, mode_matrix, frequencies, eigenvectors, ierr)
        
        if (ierr /= 0) return
        
        ! Extract growth rates (imaginary part of frequency)
        growth_rates = aimag(frequencies)
        
        ! Sort by growth rate (most unstable first)
        call sort_eigenvalues_by_growth_rate(frequencies, eigenvectors, growth_rates)
        
    end subroutine solve_mode_eigenvalue_problem
    
    !---------------------------------------------------------------------------
    ! Sort eigenvalues by growth rate (imaginary part) in descending order
    !---------------------------------------------------------------------------
    subroutine sort_eigenvalues_by_growth_rate(frequencies, eigenvectors, growth_rates)
        complex(real64), intent(inout) :: frequencies(:)
        complex(real64), intent(inout) :: eigenvectors(:,:)
        real(real64), intent(inout) :: growth_rates(:)
        
        integer :: n, i, j, max_idx
        real(real64) :: temp_growth
        complex(real64) :: temp_freq, temp_vec(size(eigenvectors,1))
        
        n = size(frequencies)
        
        ! Simple selection sort by growth rate (descending)
        do i = 1, n-1
            max_idx = i
            do j = i+1, n
                if (growth_rates(j) > growth_rates(max_idx)) then
                    max_idx = j
                end if
            end do
            
            if (max_idx /= i) then
                ! Swap growth rates
                temp_growth = growth_rates(i)
                growth_rates(i) = growth_rates(max_idx)
                growth_rates(max_idx) = temp_growth
                
                ! Swap frequencies
                temp_freq = frequencies(i)
                frequencies(i) = frequencies(max_idx)
                frequencies(max_idx) = temp_freq
                
                ! Swap eigenvectors
                temp_vec = eigenvectors(:,i)
                eigenvectors(:,i) = eigenvectors(:,max_idx)
                eigenvectors(:,max_idx) = temp_vec
            end if
        end do
        
    end subroutine sort_eigenvalues_by_growth_rate

end module kilca_eigensolver_m