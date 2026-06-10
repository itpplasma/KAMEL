/*! \file deriv_central.cpp
    \brief The implementation of the functions declared in deriv_central.h.
*/

#include <math.h>
#include <float.h>

#include "deriv_central.h"

/*******************************************************************/

//5-point central difference with truncation and rounding error estimates
static void central_diff (deriv_function f, void *params, double x, double h,
                          double *result, double *abserr_trunc,
                          double *abserr_round)
{
double fm1 = f (x - h, params);
double fp1 = f (x + h, params);
double fmh = f (x - h/2.0, params);
double fph = f (x + h/2.0, params);

double r3 = 0.5*(fp1 - fm1);
double r5 = (4.0/3.0)*(fph - fmh) - (1.0/3.0)*r3;

//rounding error of the function values, amplified by the difference rule
double e3 = (fabs (fp1) + fabs (fm1))*DBL_EPSILON;
double e5 = 2.0*(fabs (fph) + fabs (fmh))*DBL_EPSILON + e3;

//rounding error due to the finite precision of x +/- h
double dy = fmax (fabs (r3/h), fabs (r5/h))*(fabs (x)/h)*DBL_EPSILON;

*result = r5/h;
*abserr_trunc = fabs ((r5 - r3)/h);
*abserr_round = fabs (e5/h) + dy;
}

/*******************************************************************/

int deriv_central (deriv_function f, void *params, double x, double h,
                   double *result, double *abserr)
{
double r0, trunc0, round0;

central_diff (f, params, x, h, &r0, &trunc0, &round0);

*result = r0;
*abserr = trunc0 + round0;

if (round0 < trunc0 && round0 > 0.0 && trunc0 > 0.0)
{
    //balance rounding against truncation error with an improved step size
    double h_opt = h*pow (round0/(2.0*trunc0), 1.0/3.0);

    double r_opt, trunc_opt, round_opt;

    central_diff (f, params, x, h_opt, &r_opt, &trunc_opt, &round_opt);

    double err_opt = trunc_opt + round_opt;

    if (err_opt < *abserr && fabs (r_opt - r0) < 4.0*(*abserr))
    {
        *result = r_opt;
        *abserr = err_opt;
    }
}

return 0;
}
