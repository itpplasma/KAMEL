/*! \file calc_eigmode.h
    \brief The declarations of functions used to find eigenmode.
*/

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include <complex>
#include <gsl/gsl_integration.h>

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

struct quad_params
{
    cmplx_func func;
    void *params;
    complex<double> center;
    double radius;
    int part;
};

/**********************************************************************************/

int eval_det (const gsl_vector *x, void *params, gsl_vector *f);
int eval_jac (const gsl_vector *x, void *params, gsl_matrix *J);
int eval_det_jac (const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J);

int find_det_zeros (int ind, int m, int n, core_data *cd);

void print_f_search_state (size_t iter, gsl_multiroot_fsolver *s);
void print_fdf_search_state (size_t iter, gsl_multiroot_fdfsolver *s);

complex<double> calc_circle_integral(complex<double> center, double radius, cmplx_func inv_det, void *params);

complex<double> calc_circle_integral_gkq (complex<double> center, double radius, cmplx_func f, void *params);

complex<double> inv_det(complex<double> z, void *params);

complex<double> test_func(complex<double> z, void *params);

double func_on_circle (double phi, void *params);

int loop_over_frequences (int ind, int m, int n, core_data *cd);

void calc_det (double *freq, double *absdet, void *params);

int find_eigmodes (int ind, int m, int n, core_data *cd);

/**********************************************************************************/
