/* Library libcerf:
 *   Compute complex error functions, based on a new implementation of
 *   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
 *
 * File testtool.h
 *   Auxiliary functions and preprocessor macros for numeric tests.
 *
 * Copyright:
 *   (C) 2012 Massachusetts Institute of Technology
 *   (C) 2013 Forschungszentrum Jülich GmbH
 *
 * Licence:
 *   ../LICENSE
 *
 * Authors:
 *   Steven G. Johnson, Massachusetts Institute of Technology, 2012, core author
 *   Joachim Wuttke, Forschungszentrum Jülich, 2013, package maintainer
 *
 * Website:
 *   http://apps.jcns.fz-juelich.de/libcerf
 *
 * Revision history:
 *   ../CHANGELOG
 */

#ifdef __cplusplus
#include <cassert>
#include <cmath>
using std::isfinite;
#else
#include <assert.h>
#include <math.h>
#endif
#include "defs.h" // defines _cerf_cmplx, CMPLX, NaN, Inf
#include <float.h>
#include <stdio.h>


typedef struct {
    int failed;
    int total;
} result_t;

// Compute relative error |b-a|/|b+offs|, handling case of NaN and Inf.
// The tiny offset offs ensures a resonable return value for a=(almost underflowing), b=0.
static double relerr(double a, double b)
{
    if (!isfinite(a))
        return !isfinite(b) ? 0 : Inf;
    if (!isfinite(b)) {
        assert(isfinite(a)); // implied by the above
        return Inf;
    }
   return fabs((b - a)) / (fabs(b) + 1e-300);
}

// Test whether real numbers 'computed' and 'expected' agree within relative error bound 'limit'
void rtest(result_t* result, double limit, double computed, double expected, const char* name)
{
    ++result->total;
    const double re = relerr(computed, expected);
    if (re > limit) {
        printf("failure in subtest %i: %s\n", result->total, name);
        printf("- fct value %20.15g\n", computed);
        printf("- expected  %20.15g\n", expected);
        printf("=> error %6.2g above limit %6.2g\n", re, limit);
        ++result->failed;
    }
}

// Test whether complex numbers 'computed' and 'expected' agree within relative error bound 'limit'
void ztest(
    result_t* result, double limit, _cerf_cmplx computed, _cerf_cmplx expected, const char* name)
{
    ++result->total;
    const double re_r = relerr(creal(computed), creal(expected));
    const double re_i = relerr(cimag(computed), cimag(expected));
    if (re_r > limit || re_i > limit) {
        printf("failure in subtest %i: %s\n", result->total, name);
        printf("- fct value %20.15g%+20.15g\n", creal(computed), cimag(computed));
        printf("- expected  %20.15g%+20.15g\n", creal(expected), cimag(expected));
        printf("=> error %6.2g or %6.2g above limit %6.2g\n", re_r, re_i, limit);
        ++result->failed;
    }
}

// Wrap rtest; use preprocessor stringification to print the calling 'function_val' verbatim
#define RTEST(result, limit, function_val, expected_val)                                           \
    rtest(&result, limit, function_val, expected_val, #function_val);

// Wrap ztest; use preprocessor stringification to print the calling 'function_val' verbatim
#define ZTEST(result, limit, function_val, expected_val)                                           \
    ztest(&result, limit, function_val, expected_val, #function_val);
