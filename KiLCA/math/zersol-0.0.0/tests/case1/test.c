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
to find all zeros of complex analytic function f(z) = exp(3*z) + 2*z*cos(z) - 1
located in the rectangular region [-2,2] x [-2,3] of complex plane (z).

Only the basic steps are shown with all default settings.
See the following cases for more detailed and comprehensive examples.

The test function is taken from the paper:
P. Kravanja et al. / Computer Physics Communications 124 (2000) 212-232
*/

#include <stdio.h>

#include "zersolc.h"  /* declarations of the solver C interface functions and macroses */

/*
Warning: it is recommended to use the ZerSol typedefs _real_ and _complex_ for real and complex numbers
to ensure the consistence of your data types with those used in ZerSol C bindings.
You can change the type of floating point numbers in zersolc.h in the line: #define _Real double
Don't forget to rebuild the ZerSol C bindings after that!
*/

/*
The complex analytic function and its derivative should be defined as shown below:

_complex_ function_name (_complex_ z, void * p)
{
_complex_ z - value of independent complex variable z
void * p    - pointer to an object that can be used to pass any additional parameters into the function
return      - the function value as _complex_
}
*/

/********************************************************************/
/* The definition of the function: f(z) = exp(3*z) + 2*z*cos(z) - 1 */
_complex_ f (_complex_ z, void * p)
{
return cexp(3.0*z) + 2.0 * z * ccos(z) - 1.0;
}

/********************************************************************/
/* The definition of the function derivative: f'(z) = 3*exp(3*z) + 2*cos(z) - 2*z*sin(z) */
_complex_ df (_complex_ z, void * p)
{
return 3.0 * cexp(3.0*z) + 2.0 * ccos(z) - 2.0 * z * csin(z);
}

/********************************************************************/
/* The auxiliary function that is used to print the content of complex arrays Z and V (defined below) */

void print_zeros (int n_zeros, _complex_ * Z, _complex_ * V, FILE * out);

/********************************************************************/

int main (int argc, char ** argv)
{
void * solver = create_solver(f, df, 0, "", -2.0, 2.0, -2.0, 3.0);  /* pass the function and the region parameters */

int max_n_zeros = 32, n_zeros = 0;                    /* specify maximum and current number of wanted zeros */

_complex_ Z[max_n_zeros], V[max_n_zeros];             /* allocate arrays for zeros and values of the function */

find_zeros(solver, max_n_zeros, Z, V, &n_zeros);      /* the solver sets Z, V, n_zeros and returns 0 if successful */

print_status(solver, "");                             /* the solver prints information about the search status */

free_solver(solver);                                  /* free the memory allocated by the solver */

print_zeros(min(n_zeros, max_n_zeros), Z, V, stdout); /* print the content of arrays of zeros Z and function values V */

return 0;
}

/********************************************************************/

void print_zeros (int n_zeros,     /* the number of zeros */
                  _complex_ * Z,   /* array of zeros */
                  _complex_ * V,   /* array of function values at the zeros */
                  FILE * out)      /* output stream */
{
fprintf(out, "\nThe solver returned the following zeros of f(z) in the given region:\n");

fprintf(out, "%4s\t%50s\t%50s\n", "#", "(Re{zero}, Im{zero})", "(Re{f(zero)}, Im{f(zero)})");

int i;
for (i = 0; i < n_zeros; i+=1)
{
    fprintf(out, "%4d\t(%+22.16le, %+22.16le)\t(%+22.16le, %+22.16le)\n", i, creal(Z[i]), cimag(Z[i]), creal(V[i]), cimag(V[i]));
}
}

/********************************************************************/
