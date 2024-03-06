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

// Data structure that holds the function parameters:

template <typename T> struct FuncParams
{
    T deg1;  // degree 1
    T deg2;  // degree 2
    T F;     // factor
    T O1;    // frequency 1
    T O2;    // frequency 2
    T C;     // constant
};

// the typical function: deg1 = 50, deg2 = 12, F = 5, O1=20, O2=12, C = 1.
// f(z) = z^50 + z^12 - 5 * sin(20*z) * cos(12*z) - 1

/********************************************************************/

template <typename T> std::complex<T> f (const std::complex<T> & z, void * p)
{
FuncParams<T> * P = static_cast<FuncParams<T> *>(p);

return pow(z, P->deg1) + pow(z, P->deg2) - P->F * sin(P->O1 * z) * cos(P->O2 * z) - P->C;
}

/********************************************************************/

template <typename T> std::complex<T> df (const std::complex<T> & z, void * p)
{
FuncParams<T> * P = static_cast<FuncParams<T> *>(p);

return P->deg1 * pow(z, P->deg1-1) + P->deg2 * pow(z, P->deg2-1) -
       P->F * (P->O1 * cos(P->O1 * z) * cos(P->O2 * z) - P->O2 * sin(P->O1 * z) * sin(P->O2 * z));
}

/********************************************************************/
