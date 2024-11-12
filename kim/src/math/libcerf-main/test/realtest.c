/* Library libcerf:
 *   Compute complex error functions, based on a new implementation of
 *   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
 *
 * File realtest.c
 *   Compare complex and real function for various real arguments
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

#include "cerf.h"
#include "testtool.h"

// For testing the Dawson and error functions for the special case of a real argument

void xTest(
    result_t* result, const char* name, _cerf_cmplx (*F)(_cerf_cmplx), double (*FRE)(double),
    double xmin, double xmax)
{
    char info[30];
    const int n=10000;
    for (int i = 0; i < n; ++i) {
        double x = pow(10., -300. + i * 600. / (n - 1));
        if (x<xmin || x>xmax)
            continue;
        snprintf(info, 30, "%s(%g)", name, x);
        rtest(result, 1e-13, creal(F(C(x, 0))), FRE(x), info);
        rtest(result, 1e-10, creal(F(C(x, x * 1e-10))), FRE(x), info);
        rtest(result, 1e-6, creal(F(C(x, x * 1e-6))), FRE(x), info);
    }
}

// For testing the Dawson and error functions for the special case of an infinite argument

void iTest(result_t* result, const char* name, _cerf_cmplx (*F)(_cerf_cmplx), double (*FRE)(double))
{
    char info[30];
    snprintf(info, 30, "%s(Inf)", name);
    rtest(result, 1e-13, creal(F(C(Inf, 0.))), FRE(Inf), info );
    snprintf(info, 30, "%s(-Inf)", name);
    rtest(result, 1e-13, creal(F(C(-Inf, 0.))), FRE(-Inf), info);
    snprintf(info, 30, "%s(NaN)", name);
    rtest(result, 1e-13, creal(F(C(NaN, 0.))), FRE(NaN), info);
}

int main(void)
{
    result_t result = {0, 0};

    xTest(&result, "erf", cerf, erf, 1e-300, 1e300);
    iTest(&result, "erf", cerf, erf);

    xTest(&result, "erfi", cerfi, erfi, 1e-300, 1e300);
    iTest(&result, "erfi", cerfi, erfi);

    xTest(&result, "erfc", cerfc, erfc, 1e-300, 1e300);
    iTest(&result, "erfc", cerfc, erfc);

    xTest(&result, "erfcx", cerfcx, erfcx, 1e-300, 1e300);
    iTest(&result, "erfcx", cerfcx, erfcx);

    xTest(&result, "dawson", cdawson, dawson, 1e-300, 1e150);
    iTest(&result, "dawson", cdawson, dawson);

    printf("%i/%i tests failed\n", result.failed, result.total);
    return result.failed;
}
