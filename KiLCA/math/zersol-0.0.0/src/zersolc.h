/*
// Copyright (C) 2012 Ivan B. Ivanov. All rights reserved.

// Contact author: www.ivi.com or navi.adler@gmail.com.

// This file is part of ZeroSolver C++ library.

// ZeroSolver is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// ZeroSolver is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with ZeroSolver. If not, see <http://www.gnu.org/licenses/>.
*/

/*! \file
    \brief The C language bindings for the ZerSol C++ library functions.

It is recommended to use the ZerSol typedefs \_real\_ and \_complex\_ for real and complex numbers
to ensure the consistence of your data types with those used in ZerSol C bindings.

You can change the type of floating point numbers in zersolc.h in the line:
\code
#define _Real double
\endcode
Don't forget to rebuild the ZerSol C bindings after that!
*/

#ifndef ZER_SOL_C_H
#define ZER_SOL_C_H

#include <complex.h>

#ifdef __cplusplus
    extern "C"
    {
    #undef I
    #undef complex
#endif

/*!
*  \brief \_Real is a convenient alias for the real number type. Possible values: float, double (recommended).
*/
#define _Real double

/*!
*  \brief IU is a convenient alias for the imaginary unit.
*/
#define IU _Complex_I

/*!
*  \brief \_real\_ is typedef for real numbers for the ZerSol C interface.
*/
typedef _Real _real_;                 /* _real_: ZerSol type for real numbers */

/*!
*  \brief \_complex\_ is typedef for complex numbers for the ZerSol C interface.
*/
typedef _Real _Complex _complex_;     /* _complex_: ZerSol type for complex numbers */

/*!
*  \brief complex_function is typedef for complex functions for the ZerSol C interface.
*/
typedef _complex_ (*complex_function) (_complex_, void *);

/* main functions: */
/*!
*  \brief The template for a complex analytic function and its derivative.
*  \param z value of independent complex variable
*  \param p pointer to an object that can be used to pass any additional parameters into the function
*  \return  function value at z
*/
_complex_ template_function(_complex_ z, void * p);

/*!
*  \brief Create a data structure needed for the subsequent calls of the solver functions.
*  \param f           name of a C function that evaluates the complex function value
*  \param df          name of a C function that evaluates the complex function derivative
*  \param data        pointer to a structure that can be used to pass any data into the function and its derivative
*  \param description string that indicates the function purpose
*  \param xmin        minimum x-coordinate of the rectangular region
*  \param xmax        maximum x-coordinate of the rectangular region
*  \param ymin        minimum y-coordinate of the rectangular region
*  \param ymax        maximum y-coordinate of the rectangular region
*  \return            pointer to the solver structure
*/
void * create_solver(complex_function f,
                     complex_function df,
                     void * data,
                     char const * description,
                     _real_ xmin,
                     _real_ xmax,
                     _real_ ymin,
                     _real_ ymax
                    );

/*!
*  \brief Run the zeros search procedure.
*  \param solver  pointer to the solver structure
*  \param N       length of the passed Z and V arrays
*  \param Z       array for the function zeros
*  \param V       array for the function values V = F(Z)
*  \param n_zeros number of the zeros found
*  \return        search status
*/
int find_zeros(void * solver, int N, _complex_ * Z, _complex_ * V, int * n_zeros);

/*!
*  \brief Print the status of the zeros search.
*  \param solver    pointer to the solver structure
*  \param file_name name of a file to print in. The status is printed to stdout if file_name is empty string.
*/
void print_status(void * solver, const char * file_name);

/*!
*  \brief Print the dump of the solver data.
*  \param solver    pointer to the solver structure
*  \param file_name name of a file to print in. The dump is printed to stdout if file_name is empty string.
*/
void print_dump(void * solver, const char * file_name);

/*!
*  \brief Free the memory allocated by the solver.
*  \param solver pointer to the solver structure
*/
void free_solver(void * solver);

/* additional SET functions: */
/*!
*  \brief Set a finite difference method for numerical evaluation of the function derivative.
*  \param solver pointer to the solver structure
*  \param method finite difference method: pass method = 0 for 3-point stencil approximation (second order) and method = 1 for 5-point stencil approximation (fourth order)
*/
void set_nd_method(void * solver, int method);

/*!
*  \brief Set a step size (dz) for numerical evaluation of the function derivative.
*  \param solver pointer to the solver structure
*  \param h step size
*/
void set_nd_step(void * solver, _real_ h);

/*!
*  \brief Set the minimum recursion level for the winding number evaluation procedure.
*  \param solver pointer to the solver structure
*  \param min_rec_lev minimum recursion level
*/
void set_min_rec_lev(void * solver, int min_rec_lev);

/*!
*  \brief Set the maximum recursion level for the winding number evaluation procedure.
*  \param solver pointer to the solver structure
*  \param max_rec_lev maximum recursion level
*/
void set_max_rec_lev(void * solver, int max_rec_lev);

/*!
*  \brief Set a tolerance of interpolation for the winding number evaluation procedure.
*  \param solver pointer to the solver structure
*  \param interp_err tolerance of interpolation
*/
void set_interp_err(void * solver, _real_ interp_err);

/*!
*  \brief Set a tolerance of test for discontinuities for the winding number evaluation procedure.
*  \param solver pointer to the solver structure
*  \param jump_err tolerance of test for discontinuities
*/
void set_jump_err(void * solver, _real_ jump_err);

/*!
*  \brief Set a desired number of the function zeros.
*  \param solver pointer to the solver structure
*  \param n_target target number of zeros
*/
void set_n_target(void * solver, int n_target);

/*!
*  \brief Set the array of a user supplied starting values for the Newton's iterations.
*  \param solver pointer to the solver structure
*  \param n_start number of starting values
*  \param start   array of starting values
*/
void set_start_array(void * solver, int n_start, _complex_ * start);

/*!
*  \brief Set the number of automatic starting points along Re{z}.
*  \param solver pointer to the solver structure
*  \param n_split_x number of starting points
*/
void set_n_split_x(void * solver, int n_split_x);

/*!
*  \brief Set the number of automatic starting points along Im{z}.
*  \param solver pointer to the solver structure
*  \param n_split_y number of starting points
*/
void set_n_split_y(void * solver, int n_split_y);

/*!
*  \brief Set the maximum number of the Newton's iterations for each starting point.
*  \param solver pointer to the solver structure
*  \param max_iter_num maximum number of iterations
*/
void set_max_iter_num(void * solver, int max_iter_num);

/*!
*  \brief Set the multiplicity of a zero.
*  \param solver pointer to the solver structure
*  \param multiplicity expected multiplicity of a zero
*/
void set_multiplicity (void * solver, int multiplicity);

/*!
*  \brief Set the absolute and relative tolerances for the convergence condition on z.
*  \param solver pointer to the solver structure
*  \param abs_eps_z absolute tolerance
*  \param rel_eps_z relative tolerance
*/
void set_eps_for_arg(void * solver, _real_ abs_eps_z, _real_ rel_eps_z);

/*!
*  \brief Set the absolute and relative tolerances for the convergence condition on f(z).
*  \param solver pointer to the solver structure
*  \param abs_eps_f absolute tolerance
*  \param rel_eps_f relative tolerance
*/
void set_eps_for_func(void * solver, _real_ abs_eps_f, _real_ rel_eps_f);

/*!
*  \brief Set a flag that defines whether the solver will use the winding number evaluation or just Newton iterations.
*  \param solver pointer to the solver structure
*  \param use_winding flag: true - use the winding number evaluation, false - use the Newton's iterations only
*/
void set_use_winding(void * solver, int use_winding);

/*!
*  \brief Set the maximum allowed level of the rectangle partition.
*  \param solver pointer to the solver structure
*  \param max_part_level partition level
*/
void set_max_part_level(void * solver, int max_part_level);

/*!
*  \brief Set the debug level (determines the amount of internal checks).
*  \param solver pointer to the solver structure
*  \param debug_level debug level
*/
void set_debug_level(void * solver, int debug_level);

/*!
*  \brief Set the print level (determines the amount of details printed).
*  \param solver pointer to the solver structure
*  \param print_level print level
*/
void set_print_level(void * solver, int print_level);

/* additional GET functions: */
/*!
*  \brief Get a finite difference method for numerical evaluation of the function derivative.
*  \param solver pointer to the solver structure
*  \return finite difference method: method = 0 for 3-point stencil approximation (second order) and method = 1 for 5-point stencil approximation (fourth order)
*/
int get_nd_method(void * solver);

/*!
*  \brief Get a step size (dz) for numerical evaluation of the function derivative.
*  \param solver pointer to the solver structure
*  \return step size
*/
_real_ get_nd_step(void * solver);

/*!
*  \brief Get the minimum recursion level for the winding number evaluation procedure.
*  \param solver pointer to the solver structure
*  \return minimum recursion level
*/
int get_min_rec_lev(void * solver);

/*!
*  \brief Get the maximum recursion level for the winding number evaluation procedure.
*  \param solver pointer to the solver structure
*  \return maximum recursion level
*/
int get_max_rec_lev(void * solver);

/*!
*  \brief Get a tolerance of interpolation for the winding number evaluation procedure.
*  \param solver pointer to the solver structure
*  \return tolerance of interpolation
*/
_real_ get_interp_err(void * solver);

/*!
*  \brief Get a tolerance of test for discontinuities for the winding number evaluation procedure.
*  \param solver pointer to the solver structure
*  \return tolerance of test for discontinuities
*/
_real_ get_jump_err(void * solver);

/*!
*  \brief Get a desired number of the function zeros.
*  \param solver pointer to the solver structure
*  \return target number of zeros
*/
int get_n_target(void * solver);

/*!
*  \brief Get the array of a user supplied starting values for the Newton's iterations.
*  \param solver pointer to the solver structure
*  \param n_start number of starting values (pointer to the number)
*  \param start   array of starting values (pointer to the array)
*/
void get_start_array(void * solver, int * n_start, _complex_ ** start);

/*!
*  \brief Get the number of automatic starting points along Re{z}.
*  \param solver pointer to the solver structure
*  \return number of starting points
*/
int get_n_split_x(void * solver);

/*!
*  \brief Get the number of automatic starting points along Im{z}.
*  \param solver pointer to the solver structure
*  \return number of starting points
*/
int get_n_split_y(void * solver);

/*!
*  \brief Get the maximum number of the Newton's iterations for each starting point.
*  \param solver pointer to the solver structure
*  \return maximum number of iterations
*/
int get_max_iter_num(void * solver);

/*!
*  \brief Get the multiplicity of a zero.
*  \param solver pointer to the solver structure
*  \return expected multiplicity of a zero
*/
int get_multiplicity (void * solver);

/*!
*  \brief Get the absolute and relative tolerances for the convergence condition on z.
*  \param solver pointer to the solver structure
*  \param abs_eps_z absolute tolerance (pointer)
*  \param rel_eps_z relative tolerance (pointer)
*/
void get_eps_for_arg(void * solver, _real_ * abs_eps_z, _real_ * rel_eps_z);

/*!
*  \brief Get the absolute and relative tolerances for the convergence condition on f(z).
*  \param solver pointer to the solver structure
*  \param abs_eps_f absolute tolerance (pointer)
*  \param rel_eps_f relative tolerance (pointer)
*/
void get_eps_for_func(void * solver, _real_ * abs_eps_f, _real_ * rel_eps_f);

/*!
*  \brief Get a flag that defines whether the solver will use the winding number evaluation or just Newton iterations.
*  \param solver pointer to the solver structure
*  \return flag: true - use the winding number evaluation, false - use the Newton's iterations only
*/
int get_use_winding(void * solver);

/*!
*  \brief Get the maximum allowed level of the rectangle partition.
*  \param solver pointer to the solver structure
*  \return maximum partition level
*/
int get_max_part_level(void * solver);

/*!
*  \brief Get the debug level (determines the amount of internal checks).
*  \param solver pointer to the solver structure
*  \return debug level
*/
int get_debug_level(void * solver);

/*!
*  \brief Get the print level (determines the amount of details printed).
*  \param solver pointer to the solver structure
*  \return print level
*/
int get_print_level(void * solver);

#ifdef __cplusplus
    }
#endif

#endif
