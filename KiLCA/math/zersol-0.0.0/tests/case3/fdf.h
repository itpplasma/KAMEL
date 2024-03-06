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
The definition of analytic function f(z) and its derivative
*/

#include <zersolc.h>

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

// Data structure that holds the function parameters:

struct FuncParams
{
    _real_ deg1;  /* degree 1    */
    _real_ deg2;  /* degree 2    */
    _real_ F;     /* factor      */
    _real_ O1;    /* frequency 1 */
    _real_ O2;    /* frequency 2 */
    _real_ C;     /* constant    */
};

/*
the typical function: deg1 = 50, deg2 = 12, F = 5, O1=20, O2=12, C = 1.
f(z) = z^50 + z^12 - 5 * sin(20*z) * cos(12*z) - 1
*/

/********************************************************************/
/* The definition of the function: f(z) */
_complex_ f (_complex_ z, void * p)
{
struct FuncParams * P = (struct FuncParams *)(p);

return cpow(z, P->deg1) + cpow(z, P->deg2) - P->F * csin(P->O1 * z) * ccos(P->O2 * z) - P->C;
}

/********************************************************************/
/* The definition of the function derivative: f'(z) */
_complex_ df (_complex_ z, void * p)
{
struct FuncParams * P = (struct FuncParams *)(p);

return P->deg1 * cpow(z, P->deg1-1) + P->deg2 * cpow(z, P->deg2-1) -
       P->F * (P->O1 * ccos(P->O1 * z) * ccos(P->O2 * z) - P->O2 * csin(P->O1 * z) * csin(P->O2 * z));
}

/********************************************************************/
