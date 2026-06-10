/*! \file ode_rk8pd.h
    \brief Adaptive embedded Runge-Kutta integrator of order 8(7).

    Thread-safe replacement for the gsl_odeiv rk8pd stepper with the
    standard y-based step-size control: all state lives in an
    ode_rk8pd_state owned by the caller, no global data.

    The coefficients are the RK8(7)13M pair of P. J. Prince and
    J. R. Dormand, "High order embedded Runge-Kutta formulae",
    J. Comp. Appl. Math. 7 (1981) 67-75.
*/

#ifndef ODE_RK8PD_INCLUDE
#define ODE_RK8PD_INCLUDE

#include <stddef.h>

enum
{
    ODE_SUCCESS = 0,
    ODE_FAILURE = -1 //step size underflow, no progress possible
};

//right-hand side dy/dt = f(t, y); returns ODE_SUCCESS on success
typedef int (*ode_rhs) (double t, const double y[], double dydt[], void *params);

//workspace is malloc-based so that no C++ runtime symbols are required
struct ode_rk8pd_state
{
    size_t dim;
    double eps_abs;
    double eps_rel;
    double *k;    //13 stage derivatives, dim each
    double *ytmp; //stage argument
    double *y0;   //saved state for step retries
    double *yerr; //embedded error estimate

    ode_rk8pd_state (size_t n, double epsabs, double epsrel);
    ~ode_rk8pd_state ();

private:
    ode_rk8pd_state (const ode_rk8pd_state &);
    ode_rk8pd_state &operator= (const ode_rk8pd_state &);
};

/*advance the solution from *t towards t1 by one adaptive step; on success
  *t and y are updated and *h holds the suggested next step size*/
int ode_rk8pd_evolve (ode_rk8pd_state &s, ode_rhs f, void *params,
                      double *t, double t1, double *h, double y[]);

/*apply a single fixed step of size h at t; y is updated in place and
  yerr receives the embedded error estimate*/
int ode_rk8pd_step_once (ode_rk8pd_state &s, ode_rhs f, void *params,
                         double t, double h, double y[], double yerr[]);

#endif
