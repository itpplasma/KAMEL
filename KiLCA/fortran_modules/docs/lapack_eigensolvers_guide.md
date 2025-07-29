# LAPACK Eigensolver Usage Guide for KiLCA

## Overview

This document provides guidance on using LAPACK eigensolvers directly in the Fortran translation of KiLCA. LAPACK provides highly optimized routines for solving eigenvalue problems, which are essential for mode analysis and stability calculations in plasma physics.

## Eigensolver Selection Guide

### 1. Complex Matrices

#### ZGEEV - General Complex Eigenvalue Problems
- **Use for**: Non-Hermitian complex matrices
- **Computes**: All eigenvalues and optionally left/right eigenvectors
- **Interface**:
```fortran
call zgeev(jobvl, jobvr, n, A, lda, W, VL, ldvl, VR, ldvr, &
           work, lwork, rwork, info)
```
- **Key parameters**:
  - `jobvl`: 'N' (no left eigenvectors) or 'V' (compute left eigenvectors)
  - `jobvr`: 'N' (no right eigenvectors) or 'V' (compute right eigenvectors)
  - `A`: Input matrix (destroyed on output)
  - `W`: Complex eigenvalues
  - `VL`, `VR`: Left and right eigenvectors

#### ZHEEV - Hermitian Complex Eigenvalue Problems
- **Use for**: Hermitian matrices (A = A†)
- **Advantage**: Real eigenvalues, more efficient than ZGEEV
- **Interface**:
```fortran
call zheev(jobz, uplo, n, A, lda, W, work, lwork, rwork, info)
```
- **Key parameters**:
  - `jobz`: 'N' (eigenvalues only) or 'V' (eigenvalues and eigenvectors)
  - `uplo`: 'U' (upper triangle) or 'L' (lower triangle)
  - `W`: Real eigenvalues (always real for Hermitian matrices)

### 2. Real Matrices

#### DGEEV - General Real Eigenvalue Problems
- **Use for**: Non-symmetric real matrices
- **Computes**: Real and imaginary parts of eigenvalues
- **Interface**:
```fortran
call dgeev(jobvl, jobvr, n, A, lda, WR, WI, VL, ldvl, VR, ldvr, &
           work, lwork, info)
```
- **Key parameters**:
  - `WR`, `WI`: Real and imaginary parts of eigenvalues

#### DSYEV - Symmetric Real Eigenvalue Problems
- **Use for**: Symmetric matrices (A = A^T)
- **Advantage**: Real eigenvalues, more efficient
- **Interface**:
```fortran
call dsyev(jobz, uplo, n, A, lda, W, work, lwork, info)
```

## Workspace Query Pattern

All LAPACK eigensolvers support optimal workspace query:

```fortran
! Step 1: Query optimal workspace
lwork = -1
allocate(work(1))
call zgeev('N', 'V', n, A, n, W, VL, n, VR, n, work, lwork, rwork, info)
lwork = int(real(work(1)))

! Step 2: Allocate optimal workspace
deallocate(work)
allocate(work(lwork))

! Step 3: Perform actual computation
call zgeev('N', 'V', n, A, n, W, VL, n, VR, n, work, lwork, rwork, info)
```

## Replacement Strategy for C++ Code

### 1. Replace Custom Eigensolver (eigsys.f90)
The existing `eigsys` subroutine should be updated to:
```fortran
subroutine eigsys(dim, mat, evl, evc)
    integer, intent(in) :: dim
    complex(dp), intent(in) :: mat(dim,dim)
    complex(dp), intent(out) :: evl(dim)
    complex(dp), intent(out) :: evc(dim,dim)
    
    ! Local variables
    complex(dp), allocatable :: A(:,:), work(:)
    real(dp), allocatable :: rwork(:)
    integer :: lwork, info
    
    ! Make copy of input matrix (ZGEEV destroys it)
    allocate(A(dim,dim))
    A = mat
    
    ! Workspace query and allocation
    lwork = -1
    allocate(work(1), rwork(2*dim))
    call zgeev('N', 'V', dim, A, dim, evl, evc, dim, evc, dim, &
               work, lwork, rwork, info)
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))
    
    ! Compute eigenvalues and right eigenvectors
    call zgeev('N', 'V', dim, A, dim, evl, evc, dim, evc, dim, &
               work, lwork, rwork, info)
    
    if (info /= 0) then
        write(error_unit,*) 'eigsys: ZGEEV failed with info =', info
    end if
    
    deallocate(A, work, rwork)
end subroutine eigsys
```

### 2. Mode Analysis Applications

For plasma mode analysis, eigenvalue problems typically arise in:

1. **Dispersion relation solving**: Finding complex frequencies ω
2. **Stability analysis**: Determining growth rates (Im(ω))
3. **Mode structure**: Computing eigenvectors for field patterns

Example for dispersion matrix:
```fortran
module dispersion_solver
    use iso_fortran_env, only: real64
    implicit none
    
contains
    
    subroutine solve_dispersion(n, D_matrix, omega, eigenmodes, ierr)
        integer, intent(in) :: n
        complex(real64), intent(in) :: D_matrix(n,n)
        complex(real64), intent(out) :: omega(n)
        complex(real64), intent(out) :: eigenmodes(n,n)
        integer, intent(out) :: ierr
        
        ! ZGEEV for non-Hermitian dispersion matrix
        call zgeev('N', 'V', n, D_matrix, n, omega, ...)
        
        ! Sort by growth rate (imaginary part)
        call sort_by_growth_rate(omega, eigenmodes)
    end subroutine
    
end module dispersion_solver
```

## Performance Considerations

1. **Matrix Structure**: Use specialized routines when possible
   - Hermitian → ZHEEV (faster, real eigenvalues)
   - General → ZGEEV (handles all cases)

2. **Eigenvector Computation**: Only compute if needed
   - `jobvl = 'N', jobvr = 'N'` for eigenvalues only
   - Significant performance improvement for large matrices

3. **Workspace**: Always use optimal workspace query
   - Improves performance significantly
   - LAPACK's internal blocking optimized for cache

## Error Handling

```fortran
! Check info parameter after each call
if (info < 0) then
    write(error_unit,'(A,I0)') 'LAPACK error: illegal argument ', -info
else if (info > 0) then
    write(error_unit,'(A,I0)') 'LAPACK error: convergence failure at ', info
end if
```

## Common Pitfalls

1. **Matrix Destruction**: Most eigensolvers destroy input matrix
   - Always work with a copy if you need the original

2. **Array Ordering**: Fortran column-major vs C++ row-major
   - Ensure proper transposition when interfacing

3. **Complex Conjugate**: Hermitian routines expect A = A†
   - Verify matrix properties before choosing routine

## Integration with KiLCA

For the mode solver framework:
```fortran
module kilca_eigensolver_m
    use iso_fortran_env, only: real64
    use kilca_types_m
    implicit none
    
    private
    public :: solve_mode_eigenvalue_problem
    
contains
    
    subroutine solve_mode_eigenvalue_problem(mode_matrix, frequencies, &
                                            eigenvectors, growth_rates, ierr)
        complex(real64), intent(in) :: mode_matrix(:,:)
        complex(real64), allocatable, intent(out) :: frequencies(:)
        complex(real64), allocatable, intent(out) :: eigenvectors(:,:)
        real(real64), allocatable, intent(out) :: growth_rates(:)
        integer, intent(out) :: ierr
        
        integer :: n
        
        n = size(mode_matrix, 1)
        allocate(frequencies(n), eigenvectors(n,n), growth_rates(n))
        
        ! Use ZGEEV for general complex eigenvalue problem
        call compute_eigenvalues_zgeev(mode_matrix, frequencies, eigenvectors, ierr)
        
        ! Extract growth rates
        growth_rates = aimag(frequencies)
        
        ! Sort by growth rate (most unstable first)
        call sort_modes_by_growth_rate(frequencies, eigenvectors, growth_rates)
        
    end subroutine
    
end module kilca_eigensolver_m
```

## Summary

- Use LAPACK eigensolvers directly - no wrappers needed
- Choose the right routine for your matrix structure
- Always use workspace query for optimal performance
- Handle errors appropriately
- Make copies of matrices that need to be preserved