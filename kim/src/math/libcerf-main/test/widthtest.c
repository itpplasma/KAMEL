/* Library libcerf:
 *   Compute complex error functions, based on a new implementation of
 *   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
 *
 * File widthtest.c:
 *   Test function voigt_hwhm
 *
 * Copyright:
 *   (C) 2021 Forschungszentrum Jülich GmbH
 *
 * Licence:
 *   ../LICENSE
 *
 * Authors:
 *   Joachim Wuttke, Forschungszentrum Jülich, 2021
 *
 * Website:
 *   http://apps.jcns.fz-juelich.de/libcerf
 *
 * Revision history:
 *   ../CHANGELOG
 *
 * More information:
 *   man 3 voigt_hwhm
 */

#include "cerf.h"
#include "testtool.h"

// excellent approximation [Olivero & Longbothum, 1977], used as starting value in voigt_hwhm
double hwhm0(double sigma, double gamma)
{
    return .5*(1.06868*gamma+sqrt(0.86743*gamma*gamma+4*2*log(2)*sigma*sigma));
}

void widtest(result_t* result, double limit, double sigma, double gamma)
{
    ++result->total;
    const double expected = sigma*hwhm0(1., gamma/sigma);
    const double computed = voigt_hwhm(sigma, gamma);
    const double re = relerr(computed, expected);
    if (re > limit) {
        printf("failure in subtest %i: sigma=%14.9e, gamma=%14.9e\n", result->total, sigma, gamma);
        printf("- fct value %20.15g\n", computed);
        printf("- expected  %20.15g\n", expected);
        printf("=> error %6.2g above limit %6.2g\n", re, limit);
        ++result->failed;
    }
}


int main()
{
    result_t result = {0, 0};
    const int N = 100;
    const int M = 10000;
    for (int i = 0; i<=N; ++i ) {
        const double sigma = pow(10., 180*(i-N/2)/(N/2));
        for (int j = 0; j<=M; ++j ) {
            const double gamma = sigma * pow(10., 17*(j-M/2)/(M/2));
            widtest(&result, 1e-2, sigma, gamma);
        }
    }

    printf("%i/%i tests failed\n", result.failed, result.total);
    return result.failed ? 1 : 0;
}
