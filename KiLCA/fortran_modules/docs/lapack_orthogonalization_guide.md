# LAPACK Orthogonalization Routines Guide for KiLCA

## Overview

This document provides guidance on using LAPACK orthogonalization routines in the Fortran translation of KiLCA. These routines are essential for basis vector orthogonalization in the solver framework, particularly for mode analysis and field decomposition.

## QR Decomposition Routines

### 1. Complex QR Decomposition

#### ZGEQRF - Complex QR Factorization
- **Purpose**: Computes QR factorization of complex matrix A = Q * R
- **Usage**: Basis vector orthogonalization in solver
- **Interface**:
```fortran
call zgeqrf(m, n, A, lda, tau, work, lwork, info)
```
- **Parameters**:
  - `m`: Number of rows of A
  - `n`: Number of columns of A  
  - `A`: Input matrix (overwritten with R and Householder vectors)
  - `tau`: Scalar factors of elementary reflectors
  - `work`: Workspace array
  - `lwork`: Size of workspace (-1 for query)

#### ZUNGQR - Generate Orthogonal Q Matrix
- **Purpose**: Generates the complex orthogonal matrix Q from ZGEQRF output
- **Usage**: Extract orthogonal basis vectors
- **Interface**:
```fortran
call zungqr(m, n, k, A, lda, tau, work, lwork, info)
```
- **Parameters**:
  - `k`: Number of elementary reflectors (min(m,n))
  - Other parameters same as ZGEQRF

#### ZUNMQR - Apply Q Matrix
- **Purpose**: Applies Q or Q^H to another matrix
- **Usage**: Transform vectors using orthogonal basis
- **Interface**:
```fortran
call zunmqr(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
```
- **Parameters**:
  - `side`: 'L' (left) or 'R' (right)
  - `trans`: 'N' (no transpose) or 'C' (conjugate transpose)

### 2. Real QR Decomposition

#### DGEQRF, DORGQR, DORMQR
- **Same functionality as complex versions but for real matrices**
- **Use for**: Real-valued problems or when precision allows

## Usage Patterns from C++ Code

### 1. Basis Vector Orthogonalization (solver.cpp)

The C++ code shows this pattern:
```cpp
// Step 1: QR factorization
zgeqrf_(&Nw, &Nfs, ydata, &Nw, taudata, WORK, &LWORK, &INFO);

// Step 2: Generate orthogonal Q matrix  
zungqr_(&Nw, &Nfs, &Nfs, ydata, &Nw, taudata, WORK, &LWORK, &INFO);
```

**Fortran equivalent**:
```fortran
! Step 1: QR factorization
call zgeqrf(Nw, Nfs, ydata, Nw, taudata, work, lwork, info)
if (info /= 0) then
    write(error_unit,'(A,I0)') 'zgeqrf failed with info=', info
    return
end if

! Step 2: Generate orthogonal basis
call zungqr(Nw, Nfs, Nfs, ydata, Nw, taudata, work, lwork, info)
if (info /= 0) then
    write(error_unit,'(A,I0)') 'zungqr failed with info=', info
    return
end if
```

### 2. Workspace Management

Always use workspace query for optimal performance:
```fortran
! Query optimal workspace
lwork = -1
allocate(work(1))
call zgeqrf(m, n, A, m, tau, work, lwork, info)
lwork = int(real(work(1)))
deallocate(work)
allocate(work(lwork))

! Perform factorization
call zgeqrf(m, n, A, m, tau, work, lwork, info)
```

## Fortran Implementation Module

### Orthogonalization Module Structure

```fortran
module kilca_orthogonalization_m
    use iso_fortran_env, only: real64, error_unit
    implicit none
    
    private
    public :: orthogonalize_basis_vectors
    public :: qr_decomposition_complex
    public :: apply_orthogonal_transformation
    
contains
    
    !---------------------------------------------------------------------------
    ! High-level basis orthogonalization (matching C++ solver functionality)
    !---------------------------------------------------------------------------
    subroutine orthogonalize_basis_vectors(m, n, basis_vectors, ierr)
        integer, intent(in) :: m, n
        complex(real64), intent(inout) :: basis_vectors(m, n)
        integer, intent(out) :: ierr
        
        complex(real64), allocatable :: tau(:), work(:)
        integer :: k, lwork
        
        k = min(m, n)
        allocate(tau(k))
        
        ! Workspace query
        lwork = -1
        allocate(work(1))
        call zgeqrf(m, n, basis_vectors, m, tau, work, lwork, ierr)
        
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'orthogonalize_basis_vectors: workspace query failed, info=', ierr
            deallocate(tau, work)
            return
        end if
        
        lwork = int(real(work(1)))
        deallocate(work)
        allocate(work(lwork))
        
        ! QR factorization
        call zgeqrf(m, n, basis_vectors, m, tau, work, lwork, ierr)
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'orthogonalize_basis_vectors: zgeqrf failed, info=', ierr
            deallocate(tau, work)
            return
        end if
        
        ! Generate orthogonal Q matrix
        call zungqr(m, n, k, basis_vectors, m, tau, work, lwork, ierr)
        if (ierr /= 0) then
            write(error_unit,'(A,I0)') 'orthogonalize_basis_vectors: zungqr failed, info=', ierr
        end if
        
        deallocate(tau, work)
        
    end subroutine orthogonalize_basis_vectors
    
end module kilca_orthogonalization_m
```

## Integration with Mode Solver

### Application in Mode Analysis

```fortran
! In mode solver context
subroutine solve_mode_system_with_orthogonal_basis(mode_data, ierr)
    type(mode_data_t), intent(inout) :: mode_data
    integer, intent(out) :: ierr
    
    ! Orthogonalize basis vectors for numerical stability
    do iz = 1, mode_data%Nzones
        if (allocated(mode_data%zones(iz)%basis)) then
            call orthogonalize_basis_vectors(mode_data%zones(iz)%dim, &
                                           mode_data%zones(iz)%Nwaves, &
                                           mode_data%zones(iz)%basis, ierr)
            if (ierr /= 0) return
        end if
    end do
    
    ! Continue with mode calculation...
    
end subroutine
```

## Performance Considerations

1. **Workspace Size**: Always query optimal workspace
   - LAPACK provides highly optimized block algorithms
   - Proper workspace size crucial for performance

2. **Matrix Layout**: Ensure column-major ordering
   - Fortran native layout matches LAPACK expectations
   - More efficient than row-major for these operations

3. **Reuse Factorizations**: When possible, reuse QR factorizations
   - Store tau factors for multiple applications
   - Use ZUNMQR to apply Q to multiple right-hand sides

## Error Handling

```fortran
if (info < 0) then
    write(error_unit,'(A,I0)') 'LAPACK orthogonalization: illegal argument ', -info
    ierr = info
    return
else if (info > 0) then
    write(error_unit,'(A)') 'LAPACK orthogonalization: numerical failure'
    ierr = info
    return
end if
```

## Verification of Orthogonality

```fortran
subroutine verify_orthogonality(m, n, Q, is_orthogonal, max_error)
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
    
end subroutine verify_orthogonality
```

## Summary

- Use ZGEQRF + ZUNGQR for complex basis orthogonalization
- Always query optimal workspace for best performance
- Handle errors appropriately with meaningful messages  
- Verify orthogonality when debugging numerical issues
- Direct LAPACK calls are more efficient than custom implementations