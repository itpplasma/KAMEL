# SLATEC Routines Usage Guide for KiLCA

## Overview

This document provides specific guidance on which SLATEC routines to use for various mathematical operations in the KiLCA Fortran translation. SLATEC (Sandia, Los Alamos, Air Force Weapons Laboratory Technical Exchange Committee) provides a comprehensive set of mathematical routines that can replace custom C++ implementations.

## 1. Ordinary Differential Equations (ODEs)

### DDEABM - Adams-Bashforth-Moulton Method
- **Use for**: Background equilibrium calculations
- **Replaces**: Custom ODE integration in `calc_back_slatec.cpp`
- **When to use**: Initial value problems with smooth solutions
- **Example**:
```fortran
call ddeabm(deriv_func, neq, t0, y0, tout, info, rtol, atol, &
            idid, rwork, lrw, iwork, liw, rpar, ipar)
```

### Alternative: DLSODE/DVODE
- **When to use**: If DDEABM proves insufficient for stiff problems
- **Note**: Part of ODEPACK, often bundled with SLATEC

## 2. Linear Algebra

### Dense Matrix Factorization and Solution

#### DGEFA/DGESL - General Matrix Factorization/Solution
- **Use for**: General dense linear systems Ax = b
- **Replaces**: Custom LU decomposition
- **Example**:
```fortran
! Factor matrix A
call dgefa(A, lda, n, ipvt, info)
! Solve Ax = b
call dgesl(A, lda, n, ipvt, b, job)
```

#### DPOFA/DPOSL - Positive Definite Matrix Factorization/Solution
- **Use for**: Symmetric positive definite systems
- **When to use**: Covariance matrices, metric tensors
- **More efficient than**: DGEFA for SPD matrices

#### DSIFA/DSISL - Symmetric Indefinite Factorization/Solution
- **Use for**: Symmetric indefinite systems
- **When to use**: Hamiltonian systems, saddle point problems

### Recommended Modern Alternatives (if available)
- **Use LAPACK instead**: DGETRF/DGETRS (replaces DGEFA/DGESL)
- **Use LAPACK instead**: DPOTRF/DPOTRS (replaces DPOFA/DPOSL)
- **Use LAPACK instead**: DSYTRF/DSYTRS (replaces DSIFA/DSISL)

## 3. Special Functions

### Bessel Functions

#### DBESJ - Bessel Function J (real argument)
```fortran
call dbesj(x, alpha, n, y, nz)
! x: argument
! alpha: order
! n: number of functions to compute
! y: array of results J_{alpha+k-1}(x), k=1,n
```

#### DBESI - Modified Bessel Function I (real argument)
```fortran
call dbesi(x, alpha, kode, n, y, nz)
! kode: 1 for I_n(x), 2 for exp(-|x|)*I_n(x)
```

#### DBESY - Bessel Function Y (real argument)
```fortran
call dbesy(x, fnu, n, y)
```

#### Complex Bessel Functions (AMOS library)
- **ZBESJ**: Complex argument Bessel J
- **ZBESI**: Complex argument modified Bessel I
- **ZBESY**: Complex argument Bessel Y
- **ZBESK**: Complex argument modified Bessel K

### Gamma and Related Functions

#### DGAMMA - Gamma Function
```fortran
result = dgamma(x)
! Returns Gamma(x)
```

#### DLGAMA - Log Gamma Function
```fortran
result = dlgama(x)
! Returns ln(|Gamma(x)|)
```

## 4. Interpolation

### DPCHIM/DPCHEV - Piecewise Cubic Hermite Interpolation
- **Use for**: Monotone interpolation of profiles
- **Preserves**: Monotonicity and shape
- **Example**:
```fortran
! Compute derivatives
call dpchim(n, x, f, d, incfd, ierr)
! Evaluate interpolant
call dpchev(n, x, f, d, incfd, xe, fe, de, inde, ierr)
```

### Alternative: Use Existing Spline Module
- If custom spline implementation provides needed features
- Consider FITPACK routines if available

## 5. Nonlinear Equations

### DNSQE - Nonlinear System Solver
- **Use for**: Finding zeros of nonlinear systems
- **Method**: Modified Powell hybrid method
- **Example**:
```fortran
call dnsqe(fcn, jac, iopt, n, x, fvec, tol, nprint, &
           info, wa, lwa)
```

## 6. Sorting and Searching

### DPSORT - Double Precision Sort
- **Use for**: Sorting arrays
- **Features**: Can sort multiple arrays by same permutation
```fortran
call dpsort(x, n, iperm, kflag, ier)
! kflag: 1 ascending, -1 descending
```

## 7. Quadrature (Numerical Integration)

### DQAG - Adaptive Integration
- **Use for**: General purpose integration
- **From**: QUADPACK (often included with SLATEC)

### DQAWO - Oscillatory Integrands
- **Use for**: Integrals with sin(ωx) or cos(ωx)
- **Important for**: Fourier transforms, wave problems

## 8. Special Considerations for Plasma Physics

### Plasma Dispersion Function
- **Not in SLATEC**: Implement using complex error function
- **Alternative**: Use relationship to Faddeeva function
```fortran
! Z(ζ) = i√π w(ζ)
! where w is the Faddeeva function
```

### Hypergeometric Functions
- **Check**: TOMS Algorithm 707 (already in Fortran at `math/hyper/toms707.f90`)
- **Use**: Existing implementation unless SLATEC version is needed

## 9. Recommended Usage Pattern

```fortran
module slatec_wrappers
    use iso_fortran_env
    implicit none
    
    ! Import only needed SLATEC routines
    interface
        subroutine ddeabm(...)
            ! Full interface
        end subroutine
        
        function dgamma(x)
            real(real64) :: dgamma
            real(real64), intent(in) :: x
        end function
    end interface
    
contains
    
    ! Wrapper with error handling
    subroutine solve_ode_system(...)
        ! Input validation
        ! Call SLATEC routine
        ! Error handling
        ! Post-processing
    end subroutine
    
end module
```

## 10. Build Configuration

```makefile
# Link SLATEC library
LIBS += -lslatec

# Or if using CMake
find_library(SLATEC_LIBRARY slatec)
target_link_libraries(kilca_fortran ${SLATEC_LIBRARY})
```

## Summary of Key Replacements

| C++ Implementation | SLATEC Routine | Purpose |
|-------------------|----------------|---------|
| Custom ODE solver | DDEABM | Equilibrium calculation |
| std::sort | DPSORT | Array sorting |
| Custom Bessel | DBESJ/I/Y, ZBESJ/I/Y | Special functions |
| Custom gamma | DGAMMA | Factorial-like calculations |
| Custom spline | DPCHIM/DPCHEV | Profile interpolation |
| Nonlinear solver | DNSQE | Root finding |

## Notes

1. **Prefer LAPACK over SLATEC** for linear algebra when possible
2. **Check existing Fortran code** before implementing new wrappers
3. **Validate numerical results** against C++ implementation
4. **Document any deviations** from original algorithms