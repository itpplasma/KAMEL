/*! \file calc_eigmode.h
    \brief The declarations of functions used to find eigenmode.
*/

#include <complex>

#include "core.h"
#include "mode.h"

/**********************************************************************************/

typedef complex<double> (*cmplx_func)(complex<double>, void *params);

/**********************************************************************************/

struct det_params
{
    int ind;
    int m;
    int n;
    double delta;
    intptr_t cd;
};

/**********************************************************************************/

// Residual on the fortnum_vector_fn ABI: x = (Re f, Im f), f = (Re det, Im det).
void eval_det (int n, const double *x, double *f, void *params);

extern "C" int find_det_zeros (int ind, int m, int n, intptr_t cd);

complex<double> calc_circle_integral(complex<double> center, double radius, cmplx_func inv_det, void *params);

complex<double> inv_det(complex<double> z, void *params);

complex<double> test_func(complex<double> z, void *params);

extern "C" int loop_over_frequences (int ind, int m, int n, intptr_t cd);

void calc_det (double *freq, double *absdet, void *params);

extern "C" int find_eigmodes (int ind, int m, int n, intptr_t cd);

/**********************************************************************************/
