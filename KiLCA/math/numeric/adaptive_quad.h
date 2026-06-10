/*! \file adaptive_quad.h
    \brief Adaptive Gauss-Kronrod quadrature with QUADPACK QAG/QAGI semantics.

    Thread-safe replacement for gsl_integration_qag/qagi/qagiu: all state
    lives on the stack of the caller, no workspace allocation is needed.
    Algorithm and error estimate follow Piessens, de Doncker-Kapenga,
    Ueberhuber, Kahaner, "QUADPACK" (Springer 1983), public domain.
*/

#ifndef ADAPTIVE_QUAD_INCLUDE
#define ADAPTIVE_QUAD_INCLUDE

#include <stddef.h>

//integrand with user parameters, same shape as gsl_function
typedef double (*quad_integrand) (double x, void *params);

//status codes follow the QUADPACK ier convention
enum
{
    QUAD_SUCCESS = 0,
    QUAD_EMAXITER = 1, //maximum number of subdivisions reached
    QUAD_EROUND = 2,   //roundoff error prevents requested tolerance
    QUAD_ESINGULAR = 3, //bad integrand behavior in the integration interval
    QUAD_ENOCONV = 4,  //roundoff error in the extrapolation table
    QUAD_EDIVERGE = 5, //integral probably divergent or slowly convergent
    QUAD_EINVAL = 6    //invalid tolerances or rule key
};

/*adaptive integration of f over [a, b]; key selects the Gauss-Kronrod
  rule (15, 21 or 31 points), limit the maximum number of subintervals*/
int quad_qag (quad_integrand f, void *params, double a, double b,
              double epsabs, double epsrel, size_t limit, int key,
              double *result, double *abserr);

/*adaptive integration of f over (-infinity, +infinity), mapped onto
  (0, 1) and extrapolated with the epsilon algorithm (QUADPACK dqagie,
  implemented in quadpack_qagi.cpp)*/
int quad_qagi (quad_integrand f, void *params,
               double epsabs, double epsrel, size_t limit,
               double *result, double *abserr);

/*adaptive integration of f over [a, +infinity), mapped onto (0, 1)*/
int quad_qagiu (quad_integrand f, void *params, double a,
                double epsabs, double epsrel, size_t limit,
                double *result, double *abserr);

#endif
