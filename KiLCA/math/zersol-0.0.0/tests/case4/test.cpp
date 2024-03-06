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
to find all 27 zeros of complex analytic function
f(z) = sin( (z^2 + pi^2) / (z + pi (2 i - 3)) )
located in the rectangular region [-10.0, 10.0] x [-5.0, 10.0] of complex plane (z).

The test function is taken from the paper:
Mark J. Schaefer and Tilmann Bubeck. / Journal of Reliable Computing 1 (1995) 317-323
*/

#include <complex>
#include <cmath>
#include <iostream>

#include "zerosolver.hpp" // the solver functions

#include "fdf.hpp"        // definitions of the analytic functions f(z) and f'(z) and their parameters

/********************************************************************/

int main (int argc, char ** argv)
{
// The solver has default values for all settings.
// Don't change them unless you want to improve the solution or fix a problem.

typedef double P;                            // define double precision for the zeros search algorithm (recommended)

type::Function<P> F(f, df, 0, "f(z) = sin((z^2 + pi^2)/(z + pi(2i - 3)))");  // define the function object with analytic evaluation of the derivative
//type::Function<P> F(f, 0, 0, "f(z) = sin((z^2 + pi^2)/(z + pi(2i - 3)))"); // define the function object with numerical evaluation of the derivative
//F.set_nd_method(1);                                                        // set the the fourth order accurate method for numerical derivative

type::Box<P> B(-10.1, 10.0, -5.0, 10.0);     // define the rectangular region (box) in the complex plane;
                                             // the partitioning of the box should not cause crossings of
                                             // the function zeros by edges of the subrectangles
                                             // (rectangles are divided recursively into 2 equal parts);
                                             // the symmetric boxes are not generally recommended since the edges
                                             // of subrectangles will cross real and imaginary axes where zeros
                                             // are frequently located

alg::zersol::Settings<P> S;              // define the default solver setings
//S.set_min_rec_lev(6);                  // set the minimum recursion level for the function argument evaluation
//S.set_eps_for_arg(1.0e-12, 1.0e-12);   // set absolute and relative tolerances for convergence condition on z
//S.set_eps_for_func(1.0e-15, 1.0e-15);  // set absolute and relative tolerances for convergence condition on f(z)

alg::zersol::ZeroSolver<P> solver(F, B, S);      // define the solver with the specified F, B, S

int max_n_zeros = 128, n_zeros = 0;              // specify maximum and current number of wanted zeros

std::complex<P> Z[max_n_zeros], V[max_n_zeros];  // allocate arrays for zeros and values of the function

int status = solver.FindZeros(max_n_zeros, Z, V, n_zeros);  // the solver sets Z, V, n_zeros and returns 0 if successful

if (status)  // do something: throw an exception, try another settings, etc...
{
    //std::clog << solver;     // dump the solver data to log stream. The log stream can be redirected to a file if desired.
}

solver.print_status();       // the solver prints information about the search status

// dump the solver data to file:
std::ofstream file("cpp.dump", std::ofstream::out);
file << solver;
file.close();

// to plot the final partition of the rectangular region together with the zeros found
// run in the tests/ directory the following Matlab script: plot_dump('case4/cpp.dump', 0);

return 0;
}

/********************************************************************/
