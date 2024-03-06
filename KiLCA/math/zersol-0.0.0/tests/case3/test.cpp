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
to find all zeros of complex analytic function
f(z) = z^50 + z^12 - 5 * sin(20*z) * cos(12*z) - 1
located in the rectangular region [-20.3, 20.7] x [-20.3, 20.7] of complex plane (z).

In this rather comprehensive and hard example the following steps are shown:
1. How to use an external source file to define the function and its derivative.

2. How to pass any additional parameters into the function and its derivative.

3. How to adjust all parameters of the solver.
The minimum recursion level for the argument evaluation is set to 8
since the function has strong small scale variations of the amplitude.
This allows to find all 424 zeros contained in the given rectangular region.

4. How to dump the solver data to a file.

5. How to plot the final partition of the rectangular region together with the zeros found.

The test function is taken from the paper:
M. Dellnitz et al. / Journal of Computational and Applied Mathematics 138 (2002) 325-333
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
// The commented lines show a typical use of set_*() functions for the solver customization.

typedef double P;                            // define double precision for the zeros search algorithm (recommended)

FuncParams<P> data = {50, 12, 5, 20, 12, 1}; // structure that holds the function parameters

type::Function<P> F(f, df, &data, "test");   // define the function object with analytic evaluation of the derivative

//type::Function<P> F(f, 0, &data, "test");  // define the function object with numerical evaluation of the derivative
//F.set_nd_method(1);                        // set the the fourth order accurate method for numerical derivative
//F.set_nd_step(1.0e-6);                     // set the stepsize for evaluation of numerical derivative

type::Box<P> B(-20.3, 20.7, -20.3, 20.7);    // define the rectangular region (box) in the complex plane
                                             // the partitioning of the box should not cause crossings of
                                             // the function zeros by edges of the subrectangles
                                             // (rectangles are divided recursively into 2 equal parts)

alg::zersol::Settings<P> S;      // define the default solver setings

//the winding number settings:
S.set_min_rec_lev(8);            // set the minimum recursion level for the function argument evaluation
//S.set_max_rec_lev(16);         // set the maximum recursion level for the function argument evaluation
//S.set_interp_err(5.0e-2);      // set the tolerance of interpolation used in the function argument evaluation
//S.set_jump_err(5.0e-2);        // set the tolerance of test for the function argument discontinuities

//Newton iterations settings:
//S.set_n_target(100);           // set the target number of the function zeros
//int n_user = 4;                // number of a user supplied starting values
//std::complex<P> user[] = {std::complex<P>(1,0),
//                          std::complex<P>(0,1),
//                          std::complex<P>(2,2),
//                          std::complex<P>(10, 0)}; // user supplied starting values
//S.set_start_array(n_user, user);       // set user supplied starting array
//S.set_n_split_x(6);                    // set the number of automatic starting points along Re{z}
//S.set_n_split_y(6);                    // set the number of automatic starting points along Im{z}
//S.set_max_iter_num(24);                // set the maximum number of Newton iterations for each starting point
//S.set_eps_for_arg(1.0e-12, 1.0e-12);   // set absolute and relative tolerances for convergence condition on z
//S.set_eps_for_func(1.0e-15, 1.0e-15);  // set absolute and relative tolerances for convergence condition on f(z)

//the region partition settings:
//S.set_use_winding(false);   // define whether the solver will use the winding number evaluation or just simple Newton's iterations
//S.set_max_part_level(64);   // set the maximum level of the rectangle partition
//S.set_debug_level(0);       // set the debug level (amount of internal checks)
//S.set_print_level(0);       // set the print level (amount of the details printed)

alg::zersol::ZeroSolver<P> solver(F, B, S);      // define the solver with the specified F, B, S

int max_n_zeros = 1024, n_zeros = 0;             // specify maximum and current number of wanted zeros

std::complex<P> Z[max_n_zeros], V[max_n_zeros];  // allocate arrays for zeros and values of the function

int status = solver.FindZeros(max_n_zeros, Z, V, n_zeros);  // the solver sets Z, V, n_zeros and returns 0 if successful

if (status)  // do something: throw an exception, try another settings, etc...
{
    std::clog << solver;     // dump the solver data to log stream. The log stream can be redirected to a file if desired.
}

solver.print_status();       // the solver prints information about the search status

// dump the solver data to file:
std::ofstream file("cpp.dump", std::ofstream::out);
file << solver;
file.close();

// to plot the final partition of the rectangular region together with the zeros
// run in the tests/ directory the following Matlab script: plot_dump('case3/cpp.dump', 0);

return 0;
}

/********************************************************************/
