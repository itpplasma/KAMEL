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
This program demonstrates how to use the ZerSol library
to find all zeros of complex analytic function f(z) = exp(3*z) + 2*z*cos(z) - 1
located in the rectangular region [-2,2] x [-2,3] of complex plane (z).

Only the basic steps are shown with all default settings.
See the following cases for more detailed and comprehensive examples.

The test function is taken from the paper:
P. Kravanja et al. / Computer Physics Communications 124 (2000) 212-232
*/

#include <complex>
#include <cmath>
#include <iostream>

#include "zerosolver.hpp"  // definitions of the solver functions

/*
The complex analytic function and its derivative should be defined as a template function as shown below:

template <typename T> std::complex<T> function_name (const std::complex<T> & z, void * p)
{
const std::complex<T> & z - value of independent complex variable z
void * p                  - pointer to an object that can be used to pass any additional parameters into the function
return                    - the function value as std::complex<T> variable
}
*/

/********************************************************************/
// The definition of the function: f(z) = exp(3*z) + 2*z*cos(z) - 1

template <typename T> std::complex<T> f (const std::complex<T> & z, void * p)
{
return std::exp(3.0*z) + 2.0 * z * std::cos(z) - 1.0;
}

/********************************************************************/
// The definition of the function derivative: f'(z) = 3*exp(3*z) + 2*cos(z) - 2*z*sin(z)

template <typename T> std::complex<T> df (const std::complex<T> & z, void * p)
{
return 3.0 * std::exp(3.0*z) + 2.0 * std::cos(z) - 2.0 * z * std::sin(z);
}

/********************************************************************/
// The auxiliary function that is used to print the content of complex arrays Z and V (defined below)

template <typename T> void print_zeros (const int & n_zeros,
                                        const std::complex<T> * Z,
                                        const std::complex<T> * V,
                                        std::ostream & output = std::cout);

/********************************************************************/

int main (int argc, char ** argv)
{
type::Function<double> F(f, df);            // define the complex function object (with double precision)

type::Box<double> B(-2.0, 2.0, -2.0, 3.0);  // define the rectangular region (box) in the complex plane

alg::zersol::Settings<double> S;            // define the default solver setings

alg::zersol::ZeroSolver<double> solver(F, B, S);  // define the solver with the specified F, B, S

int max_n_zeros = 32, n_zeros = 0;                // specify maximum and current number of wanted zeros

std::complex<double> Z[max_n_zeros], V[max_n_zeros];  // allocate arrays for zeros and values of the function

solver.FindZeros(max_n_zeros, Z, V, n_zeros);  // the solver sets Z, V, n_zeros and returns 0 if successful

solver.print_status();       // the solver prints information about the search status

print_zeros(std::min(n_zeros, max_n_zeros), Z, V);  // print the content of arrays of zeros Z and function values V

return 0;
}

/********************************************************************/

template <typename T> void print_zeros (const int & n_zeros,               // the number of zeros
                                        const std::complex<T> * Z,         // array of zeros
                                        const std::complex<T> * V,         // array of function values at the zeros
                                        std::ostream & output = std::cout) // output stream
{
output << "\nThe solver returned the following zeros of f(z) in the given region:\n";

output << std::setw(4) << "#" << std::setw(60) << "(Re{zero}, Im{zero})" << std::setw(60) << "(Re{f(zero)}, Im{f(zero)})" << "\n";

for (int i = 0; i < n_zeros; ++i)
{
    output << std::scientific << std::setprecision(16) << std::noshowpos
           << std::setw(4)  << i << std::showpos
           << std::setw(60) << Z[i]
           << std::setw(60) << V[i] << "\n";
}
}

/********************************************************************/
