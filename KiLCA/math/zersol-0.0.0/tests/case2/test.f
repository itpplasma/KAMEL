! Copyright (C) 2012 Ivan B. Ivanov. All rights reserved.

! Contact author: www.ivi.com or navi.adler@gmail.com.

! This file is part of ZeroSolver C++ library.

! ZeroSolver is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! ZeroSolver is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with ZeroSolver. If not, see <http://www.gnu.org/licenses/>.

! This program demonstrates how to use the ZerSol library
! to find all zeros of complex analytic function f(z) = exp(3*z) + 2*z*cos(z) - 1
! located in the rectangular region [-2,2] x [-2,3] of complex plane (z).
!
! Only the basic steps are shown with all default settings.
! See the following cases for more detailed and comprehensive examples.
!
! The test function is taken from the paper:
! P. Kravanja et al. / Computer Physics Communications 124 (2000) 212-232
!
! You can change the type of floating point numbers in zersolf.cpp in the line: #define _Real double
! Don't forget to rebuild the ZerSol Fortran bindings after that!
! Make sure that precision of floating point numbers in your Fortran code
! is consistent with that used in the ZerSol Fortran bindings.

! Please note, that pp parameter defined in the Fortran code must conform the architecture of your machine:
! integer, parameter :: pp = 8    ! pp = 4 for 32-bit machines and pp = 8 for 64-bit machines

! The complex analytic function f(z) and its derivative f'(z) should be defined as shown below:
!
! subroutine function_name (z, f, p)
! integer, parameter :: prec = 8   ! prec = 4 for single precision and prec = 8 for double precision
! complex(prec) :: z               ! z - value of independent complex variable z
! complex(prec) :: f               ! f - the function value at z
! * :: p                           ! p - any additional parameter: a number, an array, etc...

!********************************************************************!

! The definition of the function: f(z) = exp(3*z) + 2*z*cos(z) - 1
subroutine func(z, f, p)

implicit none

integer, parameter :: prec = 8   ! prec = 4 for single precision and prec = 8 for double precision
complex(prec) :: z, f            ! z - value of independent complex variable, f - the function value at z
integer :: p                     ! p - additional (dummy) parameter

f = exp(3.0*z) + 2.0 * z * cos(z) - 1.0

end subroutine

!********************************************************************!

! The definition of the function derivative: f'(z) = 3*exp(3*z) + 2*cos(z) - 2*z*sin(z)
subroutine dfunc(z, df, p)

implicit none

integer, parameter :: prec = 8   ! prec = 4 for single precision and prec = 8 for double precision
complex(prec) :: z, df           ! z - value of independent complex variable, df - the derivative value at z
integer :: p                     ! p - additional (dummy) parameter

df = 3.0 * exp(3.0*z) + 2.0 * cos(z) - 2.0 * z * sin(z)

end subroutine

!********************************************************************!

program test

implicit none

integer, parameter :: pp = 8    ! pp = 4 for 32 bit machines and pp = 8 for 64 bit machines
integer, parameter :: prec = 8  ! prec = 4 for single precision and prec = 8 for double precision

integer(pp) :: solver           ! holds pointer to the ZerSol object

integer, parameter :: max_n_zeros = 32         ! specify maximum number of wanted zeros
complex(prec), dimension(max_n_zeros) :: Z, V  ! arrays for zeros and values of the function

integer :: n_zeros                     ! number of zeros
integer :: status                      ! the search status
integer :: n                           ! a counter
integer :: data = 0                    ! additional (dummy) parameter for passing to f(z) and f'(z)
real(prec) :: xmin = -2.0, xmax = 2.0  ! x boundaries of the rectangular region
real(prec) :: ymin = -2.0, ymax = 3.0  ! y boundaries of the rectangular region
character(1) :: file = ''              ! name of a file to print the solver info: dimension = len(trim(file)) + 1

external func, dfunc                   ! external subroutines for f(z) and f'(z) evaluation

integer(pp) :: f, df                   ! addresses of f(z) and f'(z) subroutines

f = loc(func)    ! address of f(z) subroutine
df = loc(dfunc)  ! address of f'(z) subroutine, set df = 0 if f'(z) is unavailable (a numerical approximation will be used instead)

! pass the function and the region parameters:
call create_solver(solver, f, df, data, xmin, xmax, ymin, ymax)

! the solver sets Z, V, n_zeros and returns 0 if successful:
call find_zeros(solver, max_n_zeros, Z, V, n_zeros, status)

! the solver prints information about the search status on the screen since 'file' is an empty string
call print_status(solver, file, len(trim(file)))

! free the memory allocated by the solver:
call free_solver(solver)

! print the zeros:
write(*, *) 'The solver returned the following zeros of f(z) in the given region:'
write(*, '(A12,A42,A48)') '#', '(Re{zero}, Im{zero})', '(Re{f(zero)}, Im{f(zero)})'

do n = 1,min(n_zeros, max_n_zeros)
    write(*,*) n, Z(n), V(n)
end do

end program

!********************************************************************!
