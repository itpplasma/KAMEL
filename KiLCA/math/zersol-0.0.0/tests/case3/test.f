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
! to find all zeros of complex analytic function
! f(z) = z^50 + z^12 - 5 * sin(20*z) * cos(12*z) - 1
! located in the rectangular region [-20.3, 20.7] x [-20.3, 20.7] of complex plane (z).
!
! In this rather comprehensive and hard example the following steps are shown:
! 1. How to use an external source file to define the function and its derivative.
!
! 2. How to pass any additional parameters into the function and its derivative.
!
! 3. How to adjust all parameters of the solver.
! The minimum recursion level for the argument evaluation is set to 8
! since the function has strong small scale variations of the amplitude.
! This allows to find all 424 zeros contained in the given rectangular region.
!
! 4. How to dump the solver data to a file.
!
! 5. How to plot the final partition of the rectangular region together with the zeros found.
!
! The test function is taken from the paper:
! M. Dellnitz et al. / Journal of Computational and Applied Mathematics 138 (2002) 325-333

!Warning: it is recommended to use the ZerSol typedefs _real_ and _complex_ for real and complex numbers
!to ensure the consistence of your data types with those used in ZerSol C bindings.
!You can change the type of floating point numbers in zersolc.h in the line: #define _Real double
!Don't forget to rebuild the ZerSol C bindings after that!

!********************************************************************!

int main (int argc, char ** argv)
{
/*
The solver has default values for all settings.
Don't change them unless you want to improve the solution or fix a problem.
The commented lines show a typical use of set_*() functions for the solver customization.
*/

/*
When you define a rectangular region (box) in the complex plane
take into account that the recursive partitioning of the box should not cause crossings of
the function zeros by edges of the subrectangles
(rectangles are divided recursively into 2 equal parts);
the symmetric boxes are not generally recommended since the edges
of subrectangles will cross real and imaginary axes where zeros are frequently located
*/

struct FuncParams data = {50, 12, 5, 20, 12, 1}; /* structure that holds the function parameters */

void * solver = create_solver(f, df, (void *)&data, "test", -20.3, 20.7, -20.3, 20.7);  /* pass the function and the region parameters */

/* The winding number settings: */
set_min_rec_lev(solver, 8);             /* set the minimum recursion level for the function argument evaluation */
/* set_max_rec_lev(solver, 16);    */   /* set the maximum recursion level for the function argument evaluation */
/* set_interp_err(solver, 5.0e-2); */   /* set the tolerance of interpolation used in the function argument evaluation */
/* set_jump_err(solver, 5.0e-2);   */   /* set the tolerance of test for the function argument discontinuities */

/* Newton iterations settings: */
/* set_n_target(solver, 100); */                      /* set the target number of the function zeros */
/* int n_user = 4; */                                 /* number of a user supplied starting values */
/* _complex_ user[] = {1.0, IU, 2.0+2.0*IU, 10}; */   /* user supplied starting values */
/* set_start_array(solver, n_user, user); */          /* set user supplied starting array */
/* set_n_split_x(solver, 6); */                       /* set the number of automatic starting points along Re{z} */
/* set_n_split_y(solver, 6); */                       /* set the number of automatic starting points along Im{z} */
/* set_max_iter_num(solver, 24); */                   /* set the maximum number of Newton iterations for each starting point */
/* set_eps_for_arg(solver, 1.0e-12, 1.0e-12); */      /* set absolute and relative tolerances for convergence condition on z */
/* set_eps_for_func(solver, 1.0e-15, 1.0e-15); */     /* set absolute and relative tolerances for convergence condition on f(z) */

/* The region partition settings: */
/* set_use_winding(solver, 0); */      /* define whether the solver will use the winding number evaluation or just simple Newton's iterations */
/* set_max_part_level(solver, 64); */  /* set the maximum level of the rectangle partition */
/* set_debug_level(solver, 0); */      /* set the debug level (amount of internal checks) */
/* set_print_level(solver, 0); */      /* set the print level (amount of the details printed) */

int max_n_zeros = 1024, n_zeros = 0;              /* specify maximum and current number of wanted zeros */

_complex_ Z[max_n_zeros], V[max_n_zeros];         /* allocate arrays for zeros and values of the function */

int status = find_zeros(solver, max_n_zeros, Z, V, &n_zeros);  /* the solver sets Z, V, n_zeros and returns 0 if successful */

if (status)  /* do something: report an error, try another settings, etc... */
{
    print_dump(solver, "c.dump");                 /* dump the solver data to file */
}

print_status(solver, "");                         /* the solver prints information about the search status */

/*
to plot the final partition of the rectangular region together with the zeros found
run in the tests/ directory the following Matlab script: plot_dump('case3/c.dump', 0);
*/
print_dump(solver, "c.dump");                     /* dump the solver data to file */

free_solver(solver);                              /* free the memory allocated by the solver */

return 0;
}

/********************************************************************/
