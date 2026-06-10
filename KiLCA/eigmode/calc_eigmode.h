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
    core_data *cd;
};

/**********************************************************************************/

int eval_det (const double x[2], void *params, double f[2]);
int eval_jac (const double x[2], void *params, double J[4]);
int eval_det_jac (const double x[2], void *params, double f[2], double J[4]);

int find_det_zeros (int ind, int m, int n, core_data *cd);

void print_search_state (size_t iter, const double x[2], const double f[2]);

complex<double> calc_circle_integral(complex<double> center, double radius, cmplx_func inv_det, void *params);

complex<double> inv_det(complex<double> z, void *params);

complex<double> test_func(complex<double> z, void *params);

int loop_over_frequences (int ind, int m, int n, core_data *cd);

void calc_det (double *freq, double *absdet, void *params);

int find_eigmodes (int ind, int m, int n, core_data *cd);

/**********************************************************************************/
