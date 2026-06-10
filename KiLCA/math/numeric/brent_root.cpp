/*! \file brent_root.cpp
    \brief The implementation of the functions declared in brent_root.h.
*/

#include <math.h>
#include <float.h>

#include "brent_root.h"

/*******************************************************************/

int brent_root (root_function f, void *params, double a, double b,
                double epsabs, double epsrel, int max_iter, double *root)
{
double fa = f (a, params);
double fb = f (b, params);

*root = 0.5*(a + b);

if ((fa < 0.0 && fb < 0.0) || (fa > 0.0 && fb > 0.0)) return ROOT_EINVAL;

//c is the previous best bracket endpoint opposite to b
double c = a, fc = fa;
double d = b - a; //current step
double e = b - a; //previous step

for (int iter = 0; iter < max_iter; iter++)
{
    if (fabs (fc) < fabs (fb))
    {
        //b must carry the smaller function value
        a = b; b = c; c = a;
        fa = fb; fb = fc; fc = fa;
    }

    double tol = 2.0*DBL_EPSILON*fabs (b);
    double m = 0.5*(c - b);

    if (fb == 0.0)
    {
        *root = b;
        return ROOT_SUCCESS;
    }

    //convergence test on the bracket, as in gsl_root_test_interval
    double lo = (b < c) ? b : c;
    double hi = (b < c) ? c : b;
    double abs_lo = fabs (lo), abs_hi = fabs (hi);
    double min_abs = (lo > 0.0 || hi < 0.0) ? fmin (abs_lo, abs_hi) : 0.0;

    if (hi - lo < epsabs + epsrel*min_abs)
    {
        *root = b;
        return ROOT_SUCCESS;
    }

    if (fabs (e) < tol || fabs (fa) <= fabs (fb))
    {
        d = m; //bisection
        e = m;
    }
    else
    {
        //attempt inverse quadratic interpolation or secant step
        double s = fb/fa;
        double p, q;

        if (a == c)
        {
            p = 2.0*m*s;
            q = 1.0 - s;
        }
        else
        {
            double r = fb/fc;
            double t = fa/fc;
            p = s*(2.0*m*t*(t - r) - (b - a)*(r - 1.0));
            q = (t - 1.0)*(r - 1.0)*(s - 1.0);
        }

        if (p > 0.0) q = -q;
        else p = -p;

        if (2.0*p < fmin (3.0*m*q - fabs (tol*q), fabs (e*q)))
        {
            e = d;
            d = p/q;
        }
        else
        {
            d = m; //interpolation rejected, fall back to bisection
            e = m;
        }
    }

    a = b;
    fa = fb;

    if (fabs (d) > tol) b += d;
    else b += (m > 0.0 ? tol : -tol);

    fb = f (b, params);

    //keep the bracket: c must stay on the other side of the root
    if ((fb < 0.0 && fc < 0.0) || (fb > 0.0 && fc > 0.0))
    {
        c = a;
        fc = fa;
        d = b - a;
        e = b - a;
    }

    *root = b;
}

return ROOT_EMAXITER;
}
