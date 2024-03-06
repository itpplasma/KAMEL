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

/* function: f(z) = sin( (z^2 + pi^2) / (z + pi (2 i - 3)) ) */

/********************************************************************/
/* The definition of the function: f(z) */
_complex_ f (_complex_ z, void * p)
{
const _real_    pi = 3.141592653589793238462643383279502884197;
const _complex_ C1 = pi*pi;
const _complex_ C2 = -3.0*pi + 2.0*pi*IU;

return csin( (z*z + C1) / (z + C2) );
}

/********************************************************************/
/* The definition of the function derivative: f'(z) */
_complex_ df (_complex_ z, void * p)
{
const _real_    pi = 3.141592653589793238462643383279502884197;
const _complex_ C1 = pi*pi;
const _complex_ C2 = -3.0*pi + 2.0*pi*IU;

return ccos( (z*z + C1) / (z + C2) ) * ( ( (2.0*z) * (z + C2) - (z*z + C1) )/(z + C2)/(z + C2) );
}

/********************************************************************/
