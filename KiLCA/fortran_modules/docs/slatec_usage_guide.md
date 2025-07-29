# SLATEC Usage Patterns in KiLCA

## Overview

This document maps SLATEC (Sandia, Los Alamos, Air Force Weapons Laboratory Technical Exchange Committee) library usage from the C++ KiLCA code to guide the Fortran translation. SLATEC is a comprehensive mathematical software library written in Fortran 77, so these routines can be called directly from modern Fortran without wrappers.

## SLATEC Functions Used in KiLCA

### 1. ODE Solvers

#### DDEABM - Adams-Bashforth-Moulton ODE Solver
- **C++ Usage**: `background/calc_back_slatec.cpp` for equilibrium calculations
- **Purpose**: Solves initial value problems for systems of ODEs
- **Fortran Call**: Direct call to `ddeabm` from SLATEC library
- **Example**:
```fortran
! Fortran interface (SLATEC is Fortran native)
call ddeabm(deriv_func, neq, r0, uval, r1, info, rtol, atol, &
            idid, rwork, lrw, iwork, liw, rpar, ipar)
```

#### DDASSL - Differential/Algebraic System Solver  
- **C++ Declaration**: Found in `slatec.h`
- **Purpose**: Solves systems of differential/algebraic equations
- **Fortran Call**: Direct call available
- **Status**: Declaration only, usage not found in current codebase

### 2. Nonlinear Equation Solvers

#### DNSQE - Nonlinear System of Equations Solver
- **C++ Declaration**: Found in `slatec.h`
- **Purpose**: Solves systems of nonlinear equations
- **Fortran Call**: Direct call available
- **Status**: Declaration only, usage not found in current codebase

### 3. Sorting Utilities

#### DPSORT - Double Precision Sort
- **C++ Declaration**: Found in `slatec.h`
- **Purpose**: Sorts arrays in ascending/descending order
- **Fortran Call**: Direct call available
- **Status**: Declaration only, usage not found in current codebase

## Special Functions (Non-SLATEC)

The C++ code also uses special functions that may come from other sources:

### Bessel Functions
- **Location**: `math/bessel/` directory contains custom implementations
- **AMOS Library**: Complex Bessel functions (zbesi, zbesj, zbesy)
- **Recommendation**: Use SLATEC's Bessel functions or existing Fortran implementations

### Hypergeometric Functions  
- **Location**: `math/hyper/` directory with TOMS Algorithm 707
- **File**: `toms707.f90` - already in Fortran!
- **Recommendation**: Use existing Fortran implementation directly

### Gamma Functions
- **Usage**: Mathematical constants and special function evaluations
- **Recommendation**: Use SLATEC's dgamma function

## Migration Strategy

### Direct SLATEC Calls (No Translation Needed)

Since SLATEC is written in Fortran, the translation strategy is straightforward:

1. **Link SLATEC Library**: Add SLATEC to the Fortran build system
2. **Use Modules**: Create interface modules for type safety (optional)
3. **Direct Calls**: Call SLATEC routines directly without wrappers

### Example Interface Module

```fortran
module slatec_interfaces
    use iso_fortran_env, only: real64
    implicit none
    
    interface
        ! ODE solver
        subroutine ddeabm(df, neq, t, y, tout, info, rtol, atol, &
                          idid, rwork, lrw, iwork, liw, rpar, ipar)
            import :: real64
            external :: df
            integer, intent(in) :: neq, lrw, liw
            real(real64), intent(inout) :: t, y(*), tout
            integer, intent(inout) :: info(15), idid, iwork(liw)
            real(real64), intent(inout) :: rtol, atol, rwork(lrw)
            real(real64), intent(inout) :: rpar(*)
            integer, intent(inout) :: ipar(*)
        end subroutine ddeabm
    end interface
    
end module slatec_interfaces
```

## Error Handling

SLATEC uses specific error codes that need to be mapped:

### DDEABM Error Codes (idid parameter)
- `idid = 1`: Normal return, reached tout
- `idid = 2`: Output at intermediate point
- `idid = 3`: Integration not completed
- `idid < 0`: Error conditions

### Fortran Error Handling Pattern
```fortran
integer :: idid, ierr

call ddeabm(...)
if (idid < 2) then
    write(error_unit,'(A,I0)') "DDEABM warning/error: idid = ", idid
    ierr = idid
    return
end if
```

## Required SLATEC Modules

Based on the C++ code analysis, the minimal SLATEC requirements are:

1. **DDEABM** - For equilibrium ODE integration
2. **DPSORT** - If sorting is needed (optional)
3. **DDASSL** - If DAE solving is needed (optional)
4. **DNSQE** - If nonlinear system solving is needed (optional)

## Build System Integration

Add to CMakeLists.txt or Makefile:
```cmake
# Find or build SLATEC library
find_library(SLATEC_LIB slatec)
target_link_libraries(kilca_fortran ${SLATEC_LIB})
```

Or use package manager:
```bash
# Debian/Ubuntu
apt-get install libslatec-dev

# Or build from source
# SLATEC is available from netlib.org
```

## Summary

- SLATEC is already in Fortran - no translation needed
- Only DDEABM is actively used in the current C++ code
- Other SLATEC routines are declared but not used
- Special functions use custom implementations that can be replaced with SLATEC equivalents
- Direct Fortran calls are more efficient than C++ wrappers