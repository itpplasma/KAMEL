# SLATEC Error Handling Mapping Guide

## Overview

This document maps C++ error handling patterns for SLATEC routines to Fortran conventions. Since SLATEC is native Fortran, error handling is more direct in Fortran than through C++ wrappers.

## DDEABM Error Handling

### C++ Pattern
```cpp
// Initialize info array
for (i=0; i<15; i++) info[i] = 0;
info[3] = 1;  // Request intermediate output

// Call DDEABM
ddeabm_(..., &idid, ...);

// Check error
if (idid < 2) fprintf(stderr, "\nwarning: calculate_equilibrium: r=%le\tidid=%d\n", r1, idid);

// Continue with next step
info[0] = 1;  // Signal continuation
```

### Fortran Pattern
```fortran
module slatec_error_handling
    use iso_fortran_env, only: real64, error_unit
    implicit none
    
contains
    
    !---------------------------------------------------------------------------
    ! Handle DDEABM return codes
    !---------------------------------------------------------------------------
    subroutine handle_ddeabm_status(idid, r, ierr)
        integer, intent(in) :: idid
        real(real64), intent(in) :: r
        integer, intent(out) :: ierr
        
        ierr = 0
        
        select case(idid)
        case(1)
            ! Normal return - reached tout successfully
            return
            
        case(2)
            ! Intermediate output - normal for info(3)=1
            return
            
        case(3)
            ! Integration not completed due to more than MAXNUM steps
            write(error_unit,'(A,ES12.5,A,I0)') &
                "WARNING: DDEABM max steps at r=", r, ", idid=", idid
            ierr = 1
            
        case(-1)
            ! Invalid input parameters
            write(error_unit,'(A,I0)') &
                "ERROR: DDEABM invalid input, idid=", idid
            ierr = -1
            
        case(-2)
            ! Repeated error test failures
            write(error_unit,'(A,ES12.5,A,I0)') &
                "ERROR: DDEABM error test failures at r=", r, ", idid=", idid
            ierr = -2
            
        case(-3)
            ! Repeated convergence test failures
            write(error_unit,'(A,ES12.5,A,I0)') &
                "ERROR: DDEABM convergence failures at r=", r, ", idid=", idid
            ierr = -3
            
        case(-4)
            ! Singular matrix encountered
            write(error_unit,'(A,ES12.5,A,I0)') &
                "ERROR: DDEABM singular matrix at r=", r, ", idid=", idid
            ierr = -4
            
        case default
            ! Unknown error
            write(error_unit,'(A,I0)') &
                "ERROR: DDEABM unknown error, idid=", idid
            ierr = -99
        end select
        
    end subroutine handle_ddeabm_status
    
end module slatec_error_handling
```

## Usage in Background Equilibrium Calculation

### Fortran Implementation
```fortran
subroutine calculate_equilibrium_fortran(bg, ierr)
    use slatec_error_handling
    type(background_t), intent(inout) :: bg
    integer, intent(out) :: ierr
    
    ! DDEABM parameters
    real(real64) :: rtol = 1.0e-12_real64
    real(real64) :: atol = 1.0e-12_real64
    integer, parameter :: neq = 1, lrw = 151, liw = 51
    integer :: info(15), iwork(liw), idid
    real(real64) :: rwork(lrw)
    real(real64) :: r0, r1, uval
    integer :: i, local_ierr
    
    ! Initialize info array
    info = 0
    info(4) = 1  ! Fortran arrays are 1-indexed (C++ info[3] = Fortran info(4))
    rwork(1) = bg%x(bg%dimx)
    
    ! Initial conditions
    r0 = bg%x(1)
    call spline_eval(bg%spline, r0, bg%i_q, result)
    uval = bg%sd%background_settings%B0**2 * &
           (1.0_real64 + r0**2/bg%sd%background_settings%rtor**2/result(1)**2)
    
    ! Integration loop
    do i = 2, bg%dimx
        r1 = bg%x(i)
        
        ! Call DDEABM
        call ddeabm(deriv_eqn_equil, neq, r0, uval, r1, info, rtol, atol, &
                   idid, rwork, lrw, iwork, liw, bg%sd%background_settings%rtor, bg)
        
        ! Handle return status
        call handle_ddeabm_status(idid, r1, local_ierr)
        if (local_ierr < 0) then
            ierr = local_ierr
            return
        else if (local_ierr > 0 .and. idid < 2) then
            ! Warning but continue
            write(error_unit,'(A,ES12.5,A,I0)') &
                "WARNING: calculate_equilibrium: r=", r1, ", idid=", idid
        end if
        
        ! Store result
        bg%u(i) = uval
        
        ! Signal continuation for next step
        info(1) = 1
    end do
    
    ierr = 0
    
end subroutine calculate_equilibrium_fortran
```

## Error Code Mapping Summary

### DDEABM (idid parameter)
| C++ Check | Fortran idid | Meaning | Action |
|-----------|--------------|---------|---------|
| idid >= 2 | idid = 1,2 | Success | Continue |
| idid < 2 | idid = 3 | Max steps | Warning, may continue |
| idid < 2 | idid = -1 | Invalid input | Error, stop |
| idid < 2 | idid = -2 | Error test failures | Error, stop |
| idid < 2 | idid = -3 | Convergence failures | Error, stop |
| idid < 2 | idid = -4 | Singular matrix | Error, stop |

### General SLATEC Error Conventions
- Positive return codes: Success or informational
- Zero: Normal completion
- Negative return codes: Errors
- Use `error_unit` from `iso_fortran_env` for error messages

## Other SLATEC Routines Error Patterns

### DDASSL (Differential/Algebraic Solver)
```fortran
! idid return codes similar to DDEABM
! info array setup more complex
! Refer to SLATEC documentation for details
```

### DNSQE (Nonlinear Equations)
```fortran
! info return codes:
! info = 1: Solution found
! info = 0: Invalid input
! info = 2: Max iterations
! info = 3: Too small tolerance
! info = 4: Bad progress
```

## Best Practices

1. **Always check return codes** - Don't assume success
2. **Use structured error handling** - Create dedicated error handling modules
3. **Provide context** - Include variable values in error messages
4. **Allow graceful degradation** - Some warnings may allow continuation
5. **Document error codes** - Reference SLATEC documentation

## Integration with KiLCA Error System

```fortran
module kilca_error_codes
    implicit none
    
    ! SLATEC-related error codes
    integer, parameter :: KILCA_SUCCESS = 0
    integer, parameter :: KILCA_SLATEC_INPUT_ERROR = -100
    integer, parameter :: KILCA_SLATEC_CONVERGENCE_ERROR = -101
    integer, parameter :: KILCA_SLATEC_SINGULAR_ERROR = -102
    integer, parameter :: KILCA_SLATEC_UNKNOWN_ERROR = -199
    
end module kilca_error_codes
```

## References

- SLATEC Common Mathematical Library documentation
- Individual routine documentation (DDEABM, DDASSL, etc.)
- Fortran error handling best practices