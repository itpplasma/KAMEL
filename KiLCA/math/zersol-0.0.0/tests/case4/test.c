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

/*
This program demonstrates how to use the ZerSol library
to find all 27 zeros of complex analytic function
f(z) = sin( (z^2 + pi^2) / (z + pi (2 i - 3)) )
located in the rectangular region [-10.0, 10.0] x [-5.0, 10.0] of complex plane (z).

The test function is taken from the paper:
Mark J. Schaefer and Tilmann Bubeck. / Journal of Reliable Computing 1 (1995) 317-323
*/

#include <stdio.h>

#include "zersolc.h"  /* declarations of the solver C interface functions and macroses */

#include "fdf.h"      /* definitions of the analytic functions f(z) and f'(z) and their parameters */

/*
Warning: it is recommended to use the ZerSol typedefs _real_ and _complex_ for real and complex numbers
to ensure the consistence of your data types with those used in ZerSol C bindings.
You can change the type of floating point numbers in zersolc.h in the line: #define _Real double
Don't forget to rebuild the ZerSol C bindings after that!
*/

/********************************************************************/

int main (int argc, char ** argv)
{
/*
The solver has default values for all settings.
Don't change them unless you want to improve the solution or fix a problem.
*/

/*
When you define a rectangular region (box) in the complex plane
take into account that the recursive partitioning of the box should not cause crossings of
the function zeros by edges of the subrectangles
(rectangles are divided recursively into 2 equal parts);
the symmetric boxes are not generally recommended since the edges
of subrectangles will cross real and imaginary axes where zeros are frequently located
*/

void * solver = create_solver(f, df, 0, "f(z) = sin((z^2 + pi^2)/(z + pi(2i - 3)))", -10.1, 10.0, -5.0, 10.0);  /* pass the function and the region parameters */

/* set_min_rec_lev(solver, 6);                 */ /* set the minimum recursion level for the function argument evaluation */
/* set_eps_for_arg(solver, 1.0e-12, 1.0e-12);  */ /* set absolute and relative tolerances for convergence condition on z */
/* set_eps_for_func(solver, 1.0e-15, 1.0e-15); */ /* set absolute and relative tolerances for convergence condition on f(z) */

int max_n_zeros = 128, n_zeros = 0;               /* specify maximum and current number of wanted zeros */

_complex_ Z[max_n_zeros], V[max_n_zeros];         /* allocate arrays for zeros and values of the function */

int status = find_zeros(solver, max_n_zeros, Z, V, &n_zeros);  /* the solver sets Z, V, n_zeros and returns 0 if successful */

if (status)  /* do something: report an error, try another settings, etc... */
{
    print_dump(solver, "c.dump");                 /* dump the solver data to file */
}

print_status(solver, "");                         /* the solver prints information about the search status */

/*
to plot the final partition of the rectangular region together with the zeros found
run in the tests/ directory the following Matlab script: plot_dump('case4/c.dump', 0);
*/
print_dump(solver, "c.dump");                     /* dump the solver data to file */

free_solver(solver);                              /* free the memory allocated by the solver */

return 0;
}

/********************************************************************/
