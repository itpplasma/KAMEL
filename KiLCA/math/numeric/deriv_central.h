/*! \file deriv_central.h
    \brief Adaptive central finite-difference derivative.

    Replacement for gsl_deriv_central: 5-point rule at x-h, x-h/2, x+h/2,
    x+h with the error estimated from the difference to the 3-point rule,
    followed by one step-size optimization that balances truncation
    against rounding error. Thread-safe, no global state.
*/

#ifndef DERIV_CENTRAL_INCLUDE
#define DERIV_CENTRAL_INCLUDE

typedef double (*deriv_function) (double x, void *params);

/*derivative of f at x with initial step size h; returns 0 on success*/
int deriv_central (deriv_function f, void *params, double x, double h,
                   double *result, double *abserr);

#endif
