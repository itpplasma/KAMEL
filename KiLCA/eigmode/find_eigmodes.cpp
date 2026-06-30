/*! \file find_eigmodes.cpp
    \brief The implementation of the reliable zeros search algorithm.
*/

#include <complex>
#include <cmath>
#include <iostream>

#include "zerosolver.hpp"

#include "eigmode_sett.h"
#include "calc_eigmode.h"
#include "inout.h"
#include "mode.h"

/**********************************************************************************/

template <typename T> std::complex<T> determinant (const std::complex<T> & freq, void * params)
{
int ind = ((det_params *) params)->ind;
int m   = ((det_params *) params)->m;
int n   = ((det_params *) params)->n;

core_data *cd = ((det_params *) params)->cd;

{
    std::complex<T> olab_local = 2.0*pi*freq;
    cd->mda[ind] = mode_data_create_ (m, n, real(olab_local), imag(olab_local), (intptr_t)cd->sd, (intptr_t)cd->bp, cd->sd->path2project);
}

mode_data_calc_all_mode_data_ (cd->mda[ind], 0);

complex<double> det = complex<double>(wave_data_get_det_re_(mode_data_get_wd_(cd->mda[ind])), wave_data_get_det_im_(mode_data_get_wd_(cd->mda[ind])));

//if (DEBUG_FLAG)
//{
    FILE *out;
    if (!(out = fopen ("det.dat", "a")))
    {
        fprintf (stderr, "\nFailed to open file %s\a\n", "det.dat");
    }

    fprintf (out, "\n%.20le %.20le\t%.20le %.20le", real(freq), imag(freq), real(det), imag(det));

    fclose (out);
//}

//clean up:
mode_data_destroy_ (cd->mda[ind]);
cd->mda[ind] = NULL;
clear_all_data_in_mode_data_module_ (); //clean up fortran module data

return det;
}

/**********************************************************************************/

int find_eigmodes (int ind, int m, int n, core_data *cd)
{
using namespace std;
using namespace type;
using namespace alg;

typedef double P; // define double precision for the zeros search algorithm (recommended)

double es_delta = get_eigmode_delta_ ();

det_params p = {ind, m, n, es_delta, cd};

type::Function<P> F(determinant, 0, &p, "determinant"); // define the function object with numerical evaluation of the derivative
F.set_nd_method(1);                                     // set the fourth order accurate finite difference approximation
if (es_delta > 0.0) F.set_nd_step(es_delta);            // set the stepsize for evaluation of numerical derivative

type::Box<P> B(get_eigmode_rfmin_ (), get_eigmode_rfmax_ (), get_eigmode_ifmin_ (), get_eigmode_ifmax_ ());

alg::zersol::Settings<P> S;      // define the default solver setings

//the winding number settings:
S.set_min_rec_lev(8);            // set the minimum recursion level for the function argument evaluation
//S.set_max_rec_lev(16);         // set the maximum recursion level for the function argument evaluation
//S.set_interp_err(1.0e-2);      // set the tolerance of interpolation used in the function argument evaluation
//S.set_jump_err(1.0e-2);        // set the tolerance of test for the function argument discontinuities

//Newton iterations settings:
int es_Nguess = get_eigmode_Nguess_ ();
complex<double> *es_fstart = new complex<double>[es_Nguess];
for (int k = 0; k < es_Nguess; k++)
{
    es_fstart[k] = complex<double> (get_eigmode_fstart_re_ (k), get_eigmode_fstart_im_ (k));
}

S.set_n_target(get_eigmode_n_zeros_ ());            // set the target number of the function zeros
S.set_start_array(es_Nguess, es_fstart);            // set user supplied starting array
S.set_n_split_x(4);                           // set the number of automatic starting points along Re{z}
S.set_n_split_y(4);                           // set the number of automatic starting points along Im{z}
S.set_max_iter_num(24);                       // set the maximum number of Newton iterations for each starting point
S.set_eps_for_arg(get_eigmode_eps_abs_ (), get_eigmode_eps_rel_ ());  // set absolute and relative tolerances for convergence condition on z
S.set_eps_for_func(0.0, get_eigmode_eps_res_ ());     // set absolute and relative tolerances for convergence condition on f(z)

//the region partition settings:
S.set_use_winding(get_eigmode_use_winding_ ());  // define whether the solver will use the winding number evaluation or just simple Newton's iterations
//S.set_max_part_level(64);          // set the maximum level of the rectangle partition
S.set_debug_level(1);                // set the debug level (amount of internal checks)
S.set_print_level(1);                // set the print level (amount of the details printed)

alg::zersol::ZeroSolver<P> solver(F, B, S); // define the solver with the specified F, B, S

if (DEBUG_FLAG) std::cout << S;

int max_n_zeros = 128, n_zeros = 0;              // specify maximum and current number of wanted zeros

std::complex<P> Z[max_n_zeros], V[max_n_zeros];  // allocate arrays for zeros and values of the function

int status = solver.FindZeros(max_n_zeros, Z, V, n_zeros);  // the solver sets Z, V, n_zeros and returns 0 if successful

if (status)  // do something: throw an exception, try another settings, etc...
{
    std::clog << solver;     // dump the solver data to log stream. The log stream can be redirected to a file if desired.
}

solver.print_status();       // the solver prints information about the search status

//if (DEBUG_FLAG) // dump the solver data to file
//{
    std::ofstream file("search.dump", std::ofstream::out);
    file << solver;
    file.close();
//}

delete [] es_fstart;

//output file:
char *full_name = new char[1024];
char es_fname[1024];
get_eigmode_fname_ (es_fname);
sprintf (full_name, "%s%s", cd->sd->path2project, es_fname);

FILE *out;
if (!(out = fopen (full_name, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", full_name);
}

fprintf (out, "%%#\tRe(f)\t\t\t\tIm(f)\t\t\t\tRe(det)\t\t\t\tIm(det)");

for (int i = 0; i < n_zeros; ++i)
{
    fprintf (out, "\n%5u\t%.20le  %.20le\t%.20le  %.20le", i, real(Z[i]), imag(Z[i]), real(V[i]), imag(V[i]));
}

fclose (out);

delete [] full_name;

return 0;
}

/**********************************************************************************/
