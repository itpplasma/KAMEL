# KiLCA Fortran Linear Algebra Guide

## Overview

This document describes the linear algebra operations used in KiLCA and their Fortran implementations, addressing Tasks 063A, 066A, 067A, 068A, and 074A from the translation backlog.

## Task 063A: Native Operations Replacing C++ Utilities

### Complex Number Operations
- **C++ `std::complex<double>`** → **Fortran `complex(real64)`**
  - All arithmetic operations (+, -, *, /) are native in Fortran
  - `std::abs()` → `abs()` 
  - `std::arg()` → `atan2(aimag(z), real(z))`
  - `std::conj()` → `conjg()`
  - `std::real()` → `real()`
  - `std::imag()` → `aimag()`

### Vector Operations
- **C++ custom vector ops** → **Fortran intrinsics**
  - Dot product: `dot_product(a, b)`
  - Vector norms: `norm2(v)` or `sqrt(sum(abs(v)**2))`
  - Element-wise operations: Native array syntax `a + b`, `a * b`

### Matrix Operations
- **C++ custom matrix ops** → **Fortran array intrinsics**
  - Matrix multiply: `matmul(A, B)`
  - Transpose: `transpose(A)`
  - Conjugate transpose: `conjg(transpose(A))`

## Task 066A: LAPACK Solver Usage

### Dense Linear Systems (Ax = b)

#### General Complex Systems
```fortran
! ZGESV: Solve general complex system
integer :: n, nrhs, lda, ldb, info
integer, allocatable :: ipiv(:)
complex(real64), allocatable :: A(:,:), b(:,:)

! Solve A*x = b, solution overwrites b
call zgesv(n, nrhs, A, lda, ipiv, b, ldb, info)
```

#### Hermitian Systems
```fortran
! ZHESV: Solve Hermitian complex system (more efficient)
character :: uplo = 'U'  ! Upper triangular stored
call zhesv(uplo, n, nrhs, A, lda, ipiv, b, ldb, work, lwork, info)
```

#### Positive Definite Systems
```fortran
! ZPOSV: Solve positive definite Hermitian system (most efficient)
call zposv(uplo, n, nrhs, A, lda, b, ldb, info)
```

### Error Handling
- `info = 0`: Success
- `info < 0`: Illegal value in argument -info
- `info > 0`: Singular matrix (for ZGESV, U(info,info) = 0)

## Task 067A: LAPACK Decomposition Usage

### LU Decomposition
```fortran
! ZGETRF: LU factorization with partial pivoting
call zgetrf(m, n, A, lda, ipiv, info)

! ZGETRS: Solve using LU factors
call zgetrs(trans, n, nrhs, A, lda, ipiv, b, ldb, info)
```

### Cholesky Decomposition
```fortran
! ZPOTRF: Cholesky factorization for positive definite
call zpotrf(uplo, n, A, lda, info)

! ZPOTRS: Solve using Cholesky factors
call zpotrs(uplo, n, nrhs, A, lda, b, ldb, info)
```

### QR Decomposition
```fortran
! ZGEQRF: QR factorization
call zgeqrf(m, n, A, lda, tau, work, lwork, info)

! ZUNMQR: Apply Q or Q^H
call zunmqr(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
```

### SVD Decomposition
```fortran
! ZGESVD: Singular value decomposition
call zgesvd(jobu, jobvt, m, n, A, lda, S, U, ldu, VT, ldvt, work, lwork, rwork, info)
```

## Task 068A: Sparse Matrix Library Choices

### Recommended Libraries

1. **SuiteSparse (UMFPACK)**
   - Best for general sparse systems
   - Direct solver with good performance
   - Already used in KiLCA C++ version

2. **MUMPS (Multifrontal Massively Parallel Solver)**
   - Excellent for parallel environments
   - Handles symmetric and unsymmetric systems
   - Native Fortran implementation

3. **PETSc**
   - Comprehensive suite of solvers
   - Good for iterative methods
   - Fortran interface available

### Interface Pattern
```fortran
! Example UMFPACK usage pattern
integer(c_long) :: symbolic, numeric
real(c_double) :: control(20), info_umf(90)

! Symbolic factorization
call umfpack_zi_symbolic(n, n, Ap, Ai, Ax, Az, symbolic, control, info_umf)

! Numeric factorization
call umfpack_zi_numeric(Ap, Ai, Ax, Az, symbolic, numeric, control, info_umf)

! Solve
call umfpack_zi_solve(UMFPACK_A, Ap, Ai, Ax, Az, x, xz, b, bz, numeric, control, info_umf)
```

## Task 069A: Matrix I/O Convenience Procedures

### HDF5 Format (if custom format needed)
```fortran
subroutine write_complex_matrix_hdf5(filename, matrix_name, A)
    character(len=*), intent(in) :: filename, matrix_name
    complex(real64), intent(in) :: A(:,:)
    ! Implementation using HDF5 library
end subroutine

subroutine read_complex_matrix_hdf5(filename, matrix_name, A)
    character(len=*), intent(in) :: filename, matrix_name
    complex(real64), allocatable, intent(out) :: A(:,:)
    ! Implementation using HDF5 library
end subroutine
```

### Text Format (for debugging)
```fortran
subroutine write_matrix_text(unit, A, format_spec)
    integer, intent(in) :: unit
    complex(real64), intent(in) :: A(:,:)
    character(len=*), intent(in), optional :: format_spec
    ! Write matrix in readable format
end subroutine
```

## Task 074A: Iterative Solver Library Usage

### GMRES (via PETSc or custom implementation)
```fortran
! Generalized Minimal Residual method
! Good for non-symmetric systems
call gmres(A, b, x, restart, max_iter, tol, info)
```

### BiCGSTAB (Biconjugate Gradient Stabilized)
```fortran
! Good for non-symmetric systems
! More stable than BiCG
call bicgstab(A, b, x, max_iter, tol, info)
```

### CG (Conjugate Gradient)
```fortran
! For symmetric positive definite systems only
! Most efficient when applicable
call cg(A, b, x, max_iter, tol, info)
```

### Preconditioning
- ILU (Incomplete LU): Good general purpose
- Jacobi: Simple but often effective
- SSOR: Symmetric Successive Over-Relaxation

## Performance Considerations

1. **Choose the right solver**:
   - Hermitian systems: Use ZHESV instead of ZGESV
   - Positive definite: Use ZPOSV (fastest)
   - Large sparse: Use iterative methods

2. **Workspace queries**:
   - Always query optimal workspace size first
   - Example:
   ```fortran
   lwork = -1
   call zgesv(n, nrhs, A, lda, ipiv, b, ldb, work, lwork, info)
   lwork = int(real(work(1)))
   allocate(work(lwork))
   ```

3. **Memory layout**:
   - Fortran is column-major (opposite of C++)
   - Ensure correct stride parameters (LDA, LDB)

## Migration Notes from C++

1. **LAPACK is native Fortran**: No wrappers needed
2. **Error handling**: Check INFO parameter after each call
3. **Complex numbers**: Use `complex(real64)` type
4. **Work arrays**: Many routines need workspace - query size first
5. **Pivoting arrays**: Integer arrays for row/column permutations

## Example: Complete Linear System Solution

```fortran
subroutine solve_complex_system(A, b, x, info)
    complex(real64), intent(in) :: A(:,:), b(:)
    complex(real64), intent(out) :: x(:)
    integer, intent(out) :: info
    
    integer :: n, lwork
    integer, allocatable :: ipiv(:)
    complex(real64), allocatable :: A_copy(:,:), work(:)
    
    n = size(A, 1)
    allocate(A_copy(n,n), ipiv(n))
    
    ! Copy input (LAPACK overwrites)
    A_copy = A
    x = b
    
    ! Query workspace
    allocate(work(1))
    lwork = -1
    call zgesv(n, 1, A_copy, n, ipiv, x, n, work, lwork, info)
    
    ! Solve with optimal workspace
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))
    
    call zgesv(n, 1, A_copy, n, ipiv, x, n, work, lwork, info)
    
    deallocate(A_copy, ipiv, work)
end subroutine
```