module kilca_eigtransform_m
    !---------------------------------------------------------------------------
    ! KiLCA Eigenvalue Transformation Module
    !
    ! This module provides eigenvalue transformation functionality for the KiLCA
    ! solver system, translating from the C++ eigtransform.cpp implementation.
    !
    ! Key functions:
    ! 1. coeff_start_vals - Compute coefficient starting values
    ! 2. coeffs_to_solution - Transform coefficients to physical solutions
    ! 3. eig_mat_eval - Evaluate eigenvalue matrices at radial positions
    !
    ! Author: Claude (Anthropic)
    ! Date: 2024
    !
    ! Translation from: KiLCA/solver/eigtransform.cpp
    !---------------------------------------------------------------------------
    
    use iso_fortran_env, only: real64, int32, error_unit
    use kilca_types_m
    implicit none
    
    private
    
    ! Public procedures
    public :: coeff_start_vals
    public :: coeffs_to_solution  
    public :: eig_mat_eval
    
    ! Internal parameters
    integer, parameter :: DEBUG_PRINT = 0  ! Set to 1 for debug output like C++ version
    
contains

    !---------------------------------------------------------------------------
    ! Compute coefficient starting values
    ! 
    ! Translates: int coeff_start_vals (int Nw, int Nfs, double r, double *zstart)
    !
    ! This function solves the linear system: eigmat * zstart = original_zstart
    ! to determine the coefficient starting values at a given radial position.
    !---------------------------------------------------------------------------
    subroutine coeff_start_vals(Nw, Nfs, r, zstart, ierr)
        integer, intent(in) :: Nw                    ! Number of waves
        integer, intent(in) :: Nfs                   ! Number of fundamental solutions
        real(real64), intent(in) :: r               ! Radial position
        real(real64), intent(inout) :: zstart(*)    ! Coefficient array (2*Nw*Nfs)
        integer, intent(out) :: ierr                 ! Error code
        
        ! Local variables
        real(real64), allocatable :: eigmat(:)      ! Eigenvalue matrix (2*Nw*Nw)
        integer, allocatable :: ipiv(:)             ! Pivot indices for LAPACK
        integer :: k, j, info
        
        ierr = 0
        info = 0
        
        ! Validate inputs
        if (Nw <= 0 .or. Nfs <= 0) then
            ierr = -1
            return
        end if
        
        ! Allocate working arrays
        allocate(eigmat(2*Nw*Nw), ipiv(Nw), stat=ierr)
        if (ierr /= 0) then
            ierr = -2
            return
        end if
        
        ! Evaluate eigenvalue matrix at radial position r
        call eig_mat_eval(r, Nw, eigmat, ierr)
        if (ierr /= 0) then
            deallocate(eigmat, ipiv)
            return
        end if
        
        ! Debug output (like C++ version)
        if (DEBUG_PRINT == 1) then
            do k = 0, Nfs-1
                do j = 0, 2*Nw-1
                    write(*, '(A,I0,A,I0,A,ES25.15E3)') &
                        "coeffs: k=", k, " j=", j, ": zstart=", zstart(2*Nw*k+j+1)
                end do
            end do
            
            do k = 0, Nw-1
                do j = 0, Nw-1
                    write(*, '(A,I0,A,I0,A,ES25.15E3,A,ES25.15E3,A)') &
                        "eigmat: k=", k, " j=", j, &
                        ": val=(", eigmat(2*Nw*k+2*j+1), ", ", eigmat(2*Nw*k+2*j+2), ")"
                end do
            end do
        end if
        
        ! Solve linear system: eigmat * zstart = original_zstart
        ! Using LAPACK complex general solver (ZGESV)
        ! 
        ! Note: ZGESV expects complex matrices, but we store them as real arrays
        ! with interleaved real/imaginary parts
        call ZGESV(Nw, Nfs, eigmat, Nw, ipiv, zstart, Nw, info)
        
        if (info /= 0) then
            write(error_unit, '(A,I0)') &
                "ERROR: coeff_start_vals: ZGESV failed with info=", info
            ierr = info
        end if
        
        ! Debug output for results
        if (DEBUG_PRINT == 1 .and. info == 0) then
            do k = 0, Nfs-1
                do j = 0, 2*Nw-1
                    write(*, '(A,I0,A,I0,A,ES25.15E3)') &
                        "coeffs: k=", k, " j=", j, ": zstart=", zstart(2*Nw*k+j+1)
                end do
            end do
        end if
        
        ! Clean up
        deallocate(eigmat, ipiv)
        
    end subroutine coeff_start_vals
    
    !---------------------------------------------------------------------------
    ! Transform coefficients to physical solution
    !
    ! Translates: int coeffs2solution (int Nw, int Nfs, int dim, double *rgrid, double *sol)
    !
    ! Computes: sol = eigmat * coeff at each grid point
    ! Input: sol contains coefficients
    ! Output: sol contains physical solutions
    !---------------------------------------------------------------------------
    subroutine coeffs_to_solution(Nw, Nfs, dim, rgrid, sol, ierr)
        integer, intent(in) :: Nw                   ! Number of waves
        integer, intent(in) :: Nfs                  ! Number of fundamental solutions  
        integer, intent(in) :: dim                  ! Number of grid points
        real(real64), intent(in) :: rgrid(dim)     ! Radial grid
        real(real64), intent(inout) :: sol(*)      ! Solution array (2*Nw*Nfs*dim)
        integer, intent(out) :: ierr                ! Error code
        
        ! Local variables
        real(real64), allocatable :: eigmat(:)     ! Eigenvalue matrix (2*Nw*Nw)
        real(real64), allocatable :: tmp(:)        ! Temporary storage (Neq)
        real(real64), parameter :: alpha(2) = [1.0_real64, 0.0_real64]  ! Complex 1.0
        real(real64), parameter :: beta(2) = [0.0_real64, 0.0_real64]   ! Complex 0.0
        character, parameter :: trans = 'N'        ! No transpose
        integer :: info, i, k, Neq
        
        ierr = 0
        info = 0
        Neq = 2*Nw*Nfs
        
        ! Validate inputs  
        if (Nw <= 0 .or. Nfs <= 0 .or. dim <= 0) then
            ierr = -1
            return
        end if
        
        ! Allocate working arrays
        allocate(eigmat(2*Nw*Nw), tmp(Neq), stat=ierr)
        if (ierr /= 0) then
            ierr = -2
            return
        end if
        
        ! Transform coefficients to solutions at each grid point
        do i = 1, dim
            ! Evaluate eigenvalue matrix at current grid point
            call eig_mat_eval(rgrid(i), Nw, eigmat, ierr)
            if (ierr /= 0) then
                deallocate(eigmat, tmp)
                return
            end if
            
            ! Perform matrix multiplication: tmp = eigmat * sol
            ! Using LAPACK complex matrix-matrix multiply (ZGEMM)
            ! eigmat is Nw x Nw, sol segment is Nw x Nfs, result tmp is Nw x Nfs
            call ZGEMM(trans, trans, Nw, Nfs, Nw, &
                       alpha, eigmat, Nw, &
                       sol(Neq*(i-1)+1), Nw, &
                       beta, tmp, Nw)
            
            ! Copy result back to solution array
            do k = 1, Neq
                sol(Neq*(i-1) + k) = tmp(k)
            end do
        end do
        
        ! Clean up
        deallocate(eigmat, tmp)
        
    end subroutine coeffs_to_solution
    
    !---------------------------------------------------------------------------
    ! Evaluate eigenvalue matrix at radial position
    !
    ! This is a placeholder for the actual eigenvalue matrix evaluation.
    ! In the full KiLCA implementation, this would call physics-specific
    ! routines to compute the local eigenvalue matrix based on plasma
    ! parameters and magnetic field configuration.
    !
    ! For now, we implement a simple test matrix that produces reasonable
    ! behavior for testing purposes.
    !---------------------------------------------------------------------------
    subroutine eig_mat_eval(r, Nw, eigmat, ierr)
        real(real64), intent(in) :: r              ! Radial position
        integer, intent(in) :: Nw                  ! Number of waves
        real(real64), intent(out) :: eigmat(*)     ! Eigenvalue matrix (2*Nw*Nw)
        integer, intent(out) :: ierr               ! Error code
        
        ! Local variables
        integer :: i, j, idx
        real(real64) :: omega, damping
        
        ierr = 0
        
        ! Validate inputs
        if (Nw <= 0) then
            ierr = -1
            return
        end if
        
        ! Initialize matrix to zero
        eigmat(1:2*Nw*Nw) = 0.0_real64
        
        ! Create a simple test eigenvalue matrix
        ! This represents a basic plasma dispersion relation
        do i = 1, Nw
            do j = 1, Nw
                ! Calculate index for complex element (i,j)
                ! Real part: eigmat[2*Nw*(i-1) + 2*(j-1) + 1]
                ! Imag part: eigmat[2*Nw*(i-1) + 2*(j-1) + 2]
                idx = 2*Nw*(i-1) + 2*(j-1)
                
                if (i == j) then
                    ! Diagonal elements: frequency and damping dependent on radius
                    omega = real(i, real64) * (1.0_real64 + 0.1_real64 * r)
                    damping = 0.01_real64 * real(i, real64) * (1.0_real64 + r)
                    
                    eigmat(idx + 1) = omega      ! Real part (frequency)
                    eigmat(idx + 2) = -damping   ! Imaginary part (damping)
                else
                    ! Off-diagonal coupling terms
                    eigmat(idx + 1) = 0.1_real64 * sin(r * real(i+j, real64))
                    eigmat(idx + 2) = 0.05_real64 * cos(r * real(i+j, real64))
                end if
            end do
        end do
        
    end subroutine eig_mat_eval

end module kilca_eigtransform_m