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

!> \file
!! \brief The Fortran language bindings for the ZerSol C++ library functions.
!!
!! The following subroutines demonstrate the Fortran interface of the ZerSol library functions.
!! Make sure that precision of floating point numbers in your Fortran code
!!\code{.fortran}
!! integer, parameter :: prec = 8  ! prec = 4 for single precision and prec = 8 for double precision
!!\endcode
!! is consistent with that used in the ZerSol Fortran bindings.
!! You can change the type of floating point numbers used in the ZerSol interface
!! in the file zersolf.cpp in line:
!! \code{.c}
!! #define _Real double /* float, double (recommended) */
!! \endcode
!! Don't forget to rebuild the ZerSol Fortran bindings after that!
!!
!! Please note, that pp parameter defined in a user Fortran code (see test cases)
!! must conform the architecture of your machine:
!!\code{.fortran}
!! integer, parameter :: pp = 8    ! pp = 4 for 32-bit machines and pp = 8 for 64-bit machines
!!\endcode
!<

!> \brief The template subroutine for complex analytic function f(z) and its derivative f'(z).
!! \param z value of independent complex variable
!! \param f the function value at z
!! \param p any additional parameter: a number, an array, etc...
!<
subroutine template_subroutine(z, f, p)

complex(prec) :: z               ! z - value of independent complex variable z
complex(prec) :: f               ! f - the function value at z
anything :: p                    ! p - any additional parameter: a number, an array, etc...

end subroutine

!********************************************************************!

!> \brief Create a data structure needed for the subsequent calls of the solver functions.
!! \param solver      address of the solver structure
!! \param f           address of a Fortran subroutine that evaluates the complex function value
!! \param df          address of a Fortran subroutine that evaluates the complex function derivative
!! \param data        a structure (a number, an array, etc) that can be used to pass any data into the function and its derivative
!! \param xmin        minimum x-coordinate of the rectangular region
!! \param xmax        maximum x-coordinate of the rectangular region
!! \param ymin        minimum y-coordinate of the rectangular region
!! \param ymax        maximum y-coordinate of the rectangular region
!<
subroutine create_solver(solver, f, df, data, xmin, xmax, ymin, ymax)

integer(pp) :: solver
integer(pp) :: f, df
anything :: data
real(prec) :: xmin, xmax
real(prec) :: ymin, ymax

end subroutine

!********************************************************************!

!>
!!  \brief Run the zeros search procedure.
!!  \param solver  address of the solver structure
!!  \param N       length of the passed Z and V arrays
!!  \param Z       array for the function zeros
!!  \param V       array for the function values V = F(Z)
!!  \param n_zeros number of the zeros found
!!  \param status  search status
!<
subroutine find_zeros(solver, N, Z, V, n_zeros, status)

integer(pp) :: solver
integer :: N, n_zeros, status
complex(prec), dimension(N) :: Z, V

end subroutine

!********************************************************************!

!>
!!  \brief Print the status of the zeros search.
!!  \param solver    address of the solver structure
!!  \param file_name name of a file to store the search status
!!  \param length    length of the trimmed file_name character array
!<
subroutine print_status(solver, file_name, length)

integer(pp) :: solver
character(length+1) :: file_name
integer :: length

end subroutine

!********************************************************************!

!>
!!  \brief Print the dump of the solver data.
!!  \param solver    address of the solver structure
!!  \param file_name name of a file to store the solver dump
!!  \param length    length of the trimmed file_name character array
!<
subroutine print_dump(solver, file_name, length)

integer(pp) :: solver
character(length+1) :: file_name
integer :: length

end subroutine

!********************************************************************!

!>
!!  \brief Free the memory allocated by the solver.
!!  \param solver address of the solver structure
!<
subroutine free_solver(solver)

integer(pp) :: solver

end subroutine

!********************************************************************!

!>
!!  \brief Set a finite difference method for numerical evaluation of the function derivative.
!!  \param solver address of the solver structure
!!  \param method finite difference method: pass method = 0 for 3-point stencil approximation (second order) and method = 1 for 5-point stencil approximation (fourth order)
!<
subroutine set_nd_method(solver, method)

integer(pp) :: solver
integer :: method

end subroutine

!********************************************************************!

!>
!!  \brief Set a step size (dz) for numerical evaluation of the function derivative.
!!  \param solver address of the solver structure
!!  \param h step size
!<
subroutine set_nd_step(solver, h)

integer(pp) :: solver
real(prec) :: h

end subroutine

!********************************************************************!

!>
!!  \brief Set the minimum recursion level for the winding number evaluation procedure.
!!  \param solver address of the solver structure
!!  \param min_rec_lev minimum recursion level
!<
subroutine set_min_rec_lev(solver, min_rec_lev)

integer(pp) :: solver
integer :: min_rec_lev

end subroutine

!********************************************************************!

!>
!!  \brief Set the maximum recursion level for the winding number evaluation procedure.
!!  \param solver address of the solver structure
!!  \param max_rec_lev maximum recursion level
!<
subroutine set_max_rec_lev(solver, max_rec_lev)

integer(pp) :: solver
integer :: max_rec_lev

end subroutine

!********************************************************************!

!>
!!  \brief Set a tolerance of interpolation for the winding number evaluation procedure.
!!  \param solver address of the solver structure
!!  \param interp_err tolerance of interpolation
!<
subroutine set_interp_err(solver, interp_err)

integer(pp) :: solver
real(prec) :: interp_err

end subroutine

!********************************************************************!

!>
!!  \brief Set a tolerance of test for discontinuities for the winding number evaluation procedure.
!!  \param solver address of the solver structure
!!  \param jump_err tolerance of test for discontinuities
!<
subroutine set_jump_err(solver, jump_err)

integer(pp) :: solver
real(prec) :: jump_err

end subroutine

!********************************************************************!

!>
!!  \brief Set a desired number of the function zeros.
!!  \param solver address of the solver structure
!!  \param n_target target number of zeros
!<
subroutine set_n_target(solver, n_target)

integer(pp) :: solver
integer :: n_target

end subroutine

!********************************************************************!

!>
!!  \brief Set the array of a user supplied starting values for the Newton's iterations.
!!  \param solver address of the solver structure
!!  \param n_start number of starting values
!!  \param start array of starting values
!<
subroutine set_start_array(solver, n_start, start)

integer(pp) :: solver
integer :: n_start
complex(prec), dimension(n_start) :: start

end subroutine

!********************************************************************!

!>
!!  \brief Set the number of automatic starting points along Re{z}.
!!  \param solver address of the solver structure
!!  \param n_split_x number of starting points
!<
subroutine set_n_split_x(solver, n_split_x)

integer(pp) :: solver
integer :: n_split_x

end subroutine

!********************************************************************!

!>
!!  \brief Set the number of automatic starting points along Im{z}.
!!  \param solver address of the solver structure
!!  \param n_split_y number of starting points
!<
subroutine set_n_split_y(solver, n_split_y)

integer(pp) :: solver
integer :: n_split_y

end subroutine

!********************************************************************!

!>
!!  \brief Set the maximum number of the Newton's iterations for each starting point.
!!  \param solver address of the solver structure
!!  \param max_iter_num maximum number of iterations
!<
subroutine set_max_iter_num(solver, max_iter_num)

integer(pp) :: solver
integer :: max_iter_num

end subroutine

!********************************************************************!

!>
!!  \brief Set the multiplicity of a zero.
!!  \param solver address of the solver structure
!!  \param multiplicity expected multiplicity of a zero
!<
subroutine set_multiplicity(solver, multiplicity)

integer(pp) :: solver
integer :: multiplicity

end subroutine

!********************************************************************!

!>
!!  \brief Set the absolute and relative tolerances for the convergence condition on z.
!!  \param solver address of the solver structure
!!  \param abs_eps_z absolute tolerance
!!  \param rel_eps_z relative tolerance
!<
subroutine set_eps_for_arg(solver, abs_eps_z, rel_eps_z)

integer(pp) :: solver
real(prec) :: abs_eps_z, rel_eps_z

end subroutine

!********************************************************************!

!>
!!  \brief Set the absolute and relative tolerances for the convergence condition on f(z).
!!  \param solver address of the solver structure
!!  \param abs_eps_f absolute tolerance
!!  \param rel_eps_f relative tolerance
!<
subroutine set_eps_for_func(solver, abs_eps_f, rel_eps_f)

integer(pp) :: solver
real(prec) :: abs_eps_f, rel_eps_f

end subroutine

!********************************************************************!

!>
!!  \brief Set a flag that defines whether the solver will use the winding number evaluation or just Newton iterations.
!!  \param solver address of the solver structure
!!  \param use_winding flag: true - use the winding number evaluation, false - use the Newton's iterations only
!<
subroutine set_use_winding(solver, use_winding)

integer(pp) :: solver
integer :: use_winding

end subroutine

!********************************************************************!

!>
!!  \brief Set the maximum allowed level of the rectangle partition.
!!  \param solver address of the solver structure
!!  \param max_part_level partition level
!<
subroutine set_max_part_level(solver, max_part_level)

integer(pp) :: solver
integer :: max_part_level

end subroutine

!********************************************************************!

!>
!!  \brief Set the debug level (determines the amount of internal checks).
!!  \param solver address of the solver structure
!!  \param debug_level debug level
!<
subroutine set_debug_level(solver, debug_level)

integer(pp) :: solver
integer :: debug_level

end subroutine

!********************************************************************!

!>
!!  \brief Set the print level (determines the amount of details printed).
!!  \param solver address of the solver structure
!!  \param print_level print level
!<
subroutine set_print_level(solver, print_level)

integer(pp) :: solver
integer :: print_level

end subroutine

!********************************************************************!

!>
!!  \brief Get a finite difference method for numerical evaluation of the function derivative.
!!  \param solver address of the solver structure
!!  \param method finite difference method: method = 0 for 3-point stencil approximation (second order) and method = 1 for 5-point stencil approximation (fourth order)
!<
subroutine get_nd_method(solver, method)

integer(pp) :: solver
integer :: method

end subroutine

!********************************************************************!

!>
!!  \brief Get a step size (dz) for numerical evaluation of the function derivative.
!!  \param solver address of the solver structure
!!  \param h step size
!<
subroutine get_nd_step(solver, h)

integer(pp) :: solver
real(prec) :: h

end subroutine

!********************************************************************!

!>
!!  \brief Get the minimum recursion level for the winding number evaluation procedure.
!!  \param solver address of the solver structure
!!  \param min_rec_lev minimum recursion level
!<
subroutine get_min_rec_lev(solver, min_rec_lev)

integer(pp) :: solver
integer :: min_rec_lev

end subroutine

!********************************************************************!

!>
!!  \brief Get the maximum recursion level for the winding number evaluation procedure.
!!  \param solver address of the solver structure
!!  \param max_rec_lev maximum recursion level
!<
subroutine get_max_rec_lev(solver, max_rec_lev)

integer(pp) :: solver
integer :: max_rec_lev

end subroutine

!********************************************************************!

!>
!!  \brief Get a tolerance of interpolation for the winding number evaluation procedure.
!!  \param solver address of the solver structure
!!  \param interp_err tolerance of interpolation
!<
subroutine get_interp_err(solver, interp_err)

integer(pp) :: solver
real(prec) :: interp_err

end subroutine

!********************************************************************!

!>
!!  \brief Get a tolerance of test for discontinuities for the winding number evaluation procedure.
!!  \param solver address of the solver structure
!!  \param jump_err tolerance of test for discontinuities
!<
subroutine get_jump_err(solver, jump_err)

integer(pp) :: solver
real(prec) :: jump_err

end subroutine

!********************************************************************!

!>
!!  \brief Get a desired number of the function zeros.
!!  \param solver address of the solver structure
!!  \param n_target target number of zeros
!<
subroutine get_n_target(solver, n_target)

integer(pp) :: solver
integer :: n_target

end subroutine

!********************************************************************!

!>
!!  \brief Get the array of a user supplied starting values for the Newton's iterations.
!!  \param solver address of the solver structure
!!  \param n_start number of starting values
!!  \param start array of starting values to be filled in
!<
subroutine get_start_array(solver, n_start, start)

integer(pp) :: solver
integer :: n_start
complex(prec) :: start

end subroutine

!********************************************************************!

!>
!!  \brief Get the number of automatic starting points along Re{z}.
!!  \param solver address of the solver structure
!!  \param n_split_x number of starting points
!<
subroutine get_n_split_x(solver, n_split_x)

integer(pp) :: solver
integer :: n_split_x

end subroutine

!********************************************************************!

!>
!!  \brief Get the number of automatic starting points along Im{z}.
!!  \param solver address of the solver structure
!!  \param n_split_y number of starting points
!<
subroutine get_n_split_y(solver, n_split_y)

integer(pp) :: solver
integer :: n_split_y

end subroutine

!********************************************************************!

!>
!!  \brief Get the maximum number of the Newton's iterations for each starting point.
!!  \param solver address of the solver structure
!!  \param max_iter_num maximum number of iterations
!<
subroutine get_max_iter_num(solver, max_iter_num)

integer(pp) :: solver
integer :: max_iter_num

end subroutine

!********************************************************************!

!>
!!  \brief Get the multiplicity of a zero.
!!  \param solver address of the solver structure
!!  \param multiplicity expected multiplicity of a zero
!<
subroutine get_multiplicity(solver, multiplicity)

integer(pp) :: solver
integer :: multiplicity

end subroutine

!********************************************************************!

!>
!!  \brief Get the absolute and relative tolerances for the convergence condition on z.
!!  \param solver address of the solver structure
!!  \param abs_eps_z absolute tolerance
!!  \param rel_eps_z relative tolerance
!<
subroutine get_eps_for_arg(solver, abs_eps_z, rel_eps_z)

integer(pp) :: solver
real(prec) :: abs_eps_z, rel_eps_z

end subroutine

!********************************************************************!

!>
!!  \brief Get the absolute and relative tolerances for the convergence condition on f(z).
!!  \param solver address of the solver structure
!!  \param abs_eps_f absolute tolerance
!!  \param rel_eps_f relative tolerance
!<
subroutine get_eps_for_func(solver, abs_eps_f, rel_eps_f)

integer(pp) :: solver
real(prec) :: abs_eps_f, rel_eps_f

end subroutine

!********************************************************************!

!>
!!  \brief Get a flag that defines whether the solver will use the winding number evaluation or just Newton iterations.
!!  \param solver address of the solver structure
!!  \param use_winding flag: true - use the winding number evaluation, false - use the Newton's iterations only
!<
subroutine get_use_winding(solver, use_winding)

integer(pp) :: solver
integer :: use_winding

end subroutine

!********************************************************************!

!>
!!  \brief Get the maximum allowed level of the rectangle partition.
!!  \param solver address of the solver structure
!!  \param max_part_level maximum partition level
!<
subroutine get_max_part_level(solver, max_part_level)

integer(pp) :: solver
integer :: max_part_level

end subroutine

!********************************************************************!

!>
!!  \brief Get the debug level (determines the amount of internal checks).
!!  \param solver address of the solver structure
!!  \param debug_level debug level
!<
subroutine get_debug_level(solver, debug_level)

integer(pp) :: solver
integer :: debug_level

end subroutine

!********************************************************************!

!>
!!  \brief Get the print level (determines the amount of details printed).
!!  \param solver address of the solver structure
!!  \param print_level print level
!<
subroutine get_print_level(solver, print_level)

integer(pp) :: solver
integer :: print_level

end subroutine

!********************************************************************!
