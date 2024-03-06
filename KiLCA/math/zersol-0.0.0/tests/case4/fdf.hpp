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

/*
The definition of analytic function f(z) and its derivative
*/

#include <complex>
#include <cmath>

/*
The complex analytic function and its derivative should be defined as a template function as shown below:

template <typename T> std::complex<T> function_name (const std::complex<T> & z, void * p)
{
const std::complex<T> & z - value of independent complex variable z
void * p                  - pointer to an object that can be used to pass any additional parameters into the function
return                    - the function value as std::complex<T> variable
}
*/

// function: f(z) = sin( (z^2 + pi^2) / (z + pi (2 i - 3)) )

/********************************************************************/

template <typename T> std::complex<T> f (const std::complex<T> & z, void * p)
{
const T pi = 3.141592653589793238462643383279502884197;
const std::complex<T> C1(pi*pi, 0);
const std::complex<T> C2(-3.0*pi, 2.0*pi);

return std::sin( (z*z + C1) / (z + C2) );
}

/********************************************************************/

template <typename T> std::complex<T> df (const std::complex<T> & z, void * p)
{
const T pi = 3.141592653589793238462643383279502884197;
const std::complex<T> C1(pi*pi, 0);
const std::complex<T> C2(-3.0*pi, 2.0*pi);

return std::cos( (z*z + C1) / (z + C2) ) * ( ( (2.0*z) * (z + C2) - (z*z + C1) )/(z + C2)/(z + C2) );
}

/********************************************************************/
