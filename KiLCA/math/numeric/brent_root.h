/*! \file brent_root.h
    \brief Scalar root bracketing solver using the Brent-Dekker method.

    Replacement for the gsl_root_fsolver brent solver; algorithm after
    R. P. Brent, "Algorithms for Minimization without Derivatives"
    (Prentice-Hall 1973), ch. 4. Thread-safe, no global state.
*/

#ifndef BRENT_ROOT_INCLUDE
#define BRENT_ROOT_INCLUDE

typedef double (*root_function) (double x, void *params);

enum
{
    ROOT_SUCCESS = 0,
    ROOT_EMAXITER = 1, //no convergence within max_iter iterations
    ROOT_EINVAL = 2    //f(a) and f(b) do not bracket a root
};

/*find the root of f in the bracket [a, b]; iterates until the bracket
  width satisfies |b - a| < epsabs + epsrel*min(|a|, |b|)*/
int brent_root (root_function f, void *params, double a, double b,
                double epsabs, double epsrel, int max_iter, double *root);

#endif
